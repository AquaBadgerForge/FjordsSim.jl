module Simulations

export SimulationConfig,
    simulation_architecture,
    build_simulation,
    run_simulation,
    coupled_hydrostatic_simulation

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using NumericalEarth

using ..Configs:
    AbstractSimulationConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    FjordConfig,
    bathymetry_path,
    forcing_path,
    river_forcing_path,
    results_path
using ..Utils: recursive_merge, progress, cell_advection_timescale_coupled_model
using ..Atmospheres: prescribed_atmosphere, prescribed_radiation
using ..Forcing: forcing_from_file, interpolation_architecture
using ..BoundaryConditions: top_bottom_boundary_conditions
using ..Grids: ImmersedBoundaryGrid

"""
    SimulationConfig(; results_root, buoyancy, closure, coriolis, ...)

Simulation configuration: how a setup whose data is already prepared is run.

**No field has a default.** Every knob is a scientific choice about a particular fjord, and a
default would let one setup silently inherit another's — so a setup that forgets one gets an
`UndefKeywordError` naming it rather than a plausible-looking run. The setup file is therefore the
complete statement of what the simulation is; nothing about it is hidden in this struct.

The split between the two halves of the field list is not stylistic. `buoyancy`, `closure`,
`tracer_advection`, `momentum_advection`, `tracers`, `initial_conditions`, `coriolis`, `sea_ice`
and `biogeochemistry` are grid-independent, so they are stored as the objects
`coupled_hydrostatic_simulation` consumes. The free surface, the boundary conditions, the forcing
and the atmosphere are *not* — they need the grid, the architecture, or another config — so what
is stored here is the knob and `build_simulation` does the construction. Nothing that could be
derived from another config is a field at all: which forcing file, which open boundary and which
atmosphere reader all come from `forcing_config` and `atmosphere_config`.

# Fields
- `results_root`: Directory the output resolves against, via `results_path`.
- `output_file`: The snapshot file, via `results_path`.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `simulation_architecture`. A `Symbol`
  rather than a live `CPU()`/`GPU()` so this config's field types stay concrete and a setup file
  loads on a machine with no GPU.
- `buoyancy`, `closure`, `tracer_advection`, `momentum_advection`, `tracers`,
  `initial_conditions`, `coriolis`, `sea_ice`, `biogeochemistry`: passed to
  `coupled_hydrostatic_simulation` unchanged.
- `free_surface_cfl`: `SplitExplicitFreeSurface` needs the grid, so only the CFL is stored.
- `bottom_drag_coefficient`: passed to `top_bottom_boundary_conditions`, which needs the grid.
- `stop_time`: seconds of simulated time to run for.
- `output_interval`, `overwrite_existing`: the snapshot writer's schedule and clobber policy.
- `progress_interval`: how often the `progress` callback reports.
- `time_step_cfl`, `max_time_step`, `max_time_step_change`: the time-step wizard's settings.

A relative `output_file` resolves against `results_root`; an absolute one relocates just that
file.

Every duration is in seconds, and must be written as a `Float64` — with `Oceananigans.Units`
(`365days`, `1hour`, `3minutes`), whose constants already are one. `Base.@kwdef` on a
*parametric* struct generates a constructor that does not convert, unlike the non-parametric
configs elsewhere in FjordSim, so `stop_time = 3600` is a `MethodError` where `stop_time = 1hour`
is fine.
"""
Base.@kwdef mutable struct SimulationConfig{B,C,TA,MA,TR,I,CO,SI,BG} <: AbstractSimulationConfig
    results_root::String
    output_file::String
    architecture::Symbol
    buoyancy::B
    closure::C
    tracer_advection::TA
    momentum_advection::MA
    tracers::TR
    initial_conditions::I
    coriolis::CO
    sea_ice::SI
    biogeochemistry::BG
    free_surface_cfl::Float64
    bottom_drag_coefficient::Float64
    stop_time::Float64
    output_interval::Float64
    progress_interval::Float64
    overwrite_existing::Bool
    time_step_cfl::Float64
    max_time_step::Float64
    max_time_step_change::Float64
end

"""
    simulation_architecture(config)

Resolve `config.architecture` to a live `CPU()` or `GPU()`.

Shares the `Val`-dispatched methods of `interpolation_architecture`, so `:auto` picks the GPU when
`CUDA.functional()` and `:gpu` errors rather than silently falling back — the same behavior, and
the same message, as the forcing interpolation.
"""
simulation_architecture(config::AbstractSimulationConfig) =
    interpolation_architecture(Val(config.architecture))

"""
    open_boundary_conditions(::Val{edge})

The open boundary condition for the forcing config's `relaxation_edge`, as a nested named tuple
`recursive_merge` can merge into `top_bottom_boundary_conditions`. The open component is the
velocity normal to that edge.

Derived rather than configured: the relaxation edge is by construction the edge the regional
domain is open on, so a setup that named it twice could only ever disagree with itself.
"""
open_boundary_conditions(::Val{:south}) = (v = (south = NormalFlowBoundaryCondition(nothing),),)
open_boundary_conditions(::Val{:north}) = (v = (north = NormalFlowBoundaryCondition(nothing),),)
open_boundary_conditions(::Val{:west}) = (u = (west = NormalFlowBoundaryCondition(nothing),),)
open_boundary_conditions(::Val{:east}) = (u = (east = NormalFlowBoundaryCondition(nothing),),)

open_boundary_conditions(::Val{edge}) where {edge} = throw(
    ArgumentError("relaxation_edge must be one of (:south, :north, :west, :east), got :$edge"),
)

"""
    simulation_forcing_path(config)

Which prepared forcing file the simulation reads: the rivers-augmented copy `add_rivers` wrote
when the setup names rivers, so they are never silently dropped, and the file `prepare_forcing`
wrote when it does not. Dispatched on `config.rivers` exactly like `add_rivers`.
"""
simulation_forcing_path(config::FjordConfig) =
    simulation_forcing_path(config.forcing_config, config.forcing_config.rivers)

simulation_forcing_path(config::AbstractForcingConfig, ::Nothing) = forcing_path(config)
simulation_forcing_path(::AbstractForcingConfig, rivers::AbstractRiverConfig) =
    river_forcing_path(rivers)

"""
    forcing_prerequisite(rivers)

The subcommand that writes the file `simulation_forcing_path` picked, for its error message.
"""
forcing_prerequisite(::Nothing) = "prepare_forcing"
forcing_prerequisite(::AbstractRiverConfig) = "add_rivers"

"""
    build_simulation(config::FjordConfig)

Assemble the coupled simulation a setup describes, with its output writer, progress callback and
time-step wizard attached, and return it without running it.

Split from `run_simulation` so a run can be inspected or stepped by hand:

```julia
simulation = build_simulation(oslofjorden())
run!(simulation)
```

Every input comes from the setup's own prepared files — `bathymetry_path`, the forcing file
`simulation_forcing_path` picks, and the atmosphere the `prescribed_atmosphere` and
`prescribed_radiation` hooks read. A missing prerequisite is reported naming the step that writes
it, rather than as a read error from deep inside NetCDF.

Returns `nothing` for a setup naming no simulation config.
"""
function build_simulation(config::FjordConfig)
    simulation_config = config.simulation_config
    isnothing(simulation_config) && return nothing

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    forcing_file = simulation_forcing_path(config)
    isfile(forcing_file) || error(
        "Prepared forcing $forcing_file does not exist. Run `julia --project -m FjordSim " *
        "$(forcing_prerequisite(config.forcing_config.rivers))` for this setup first.",
    )

    architecture = simulation_architecture(simulation_config)
    grid = ImmersedBoundaryGrid(bathymetry_file, architecture, config.grid_config.halo)

    forcing = forcing_from_file(;
        grid = grid,
        filepath = forcing_file,
        tracers = simulation_config.tracers,
    )
    boundary_conditions = map(
        x -> FieldBoundaryConditions(; x...),
        recursive_merge(
            top_bottom_boundary_conditions(;
                grid = grid,
                bottom_drag_coefficient = simulation_config.bottom_drag_coefficient,
            ),
            open_boundary_conditions(Val(config.forcing_config.relaxation_edge)),
        ),
    )

    simulation = coupled_hydrostatic_simulation(
        grid,
        simulation_config.buoyancy,
        simulation_config.closure,
        simulation_config.tracer_advection,
        simulation_config.momentum_advection,
        simulation_config.tracers,
        simulation_config.initial_conditions,
        SplitExplicitFreeSurface(grid, cfl = simulation_config.free_surface_cfl),
        simulation_config.coriolis,
        forcing,
        boundary_conditions,
        prescribed_atmosphere(config.atmosphere_config, architecture),
        prescribed_radiation(config.atmosphere_config, architecture),
        simulation_config.sea_ice,
        simulation_config.biogeochemistry;
        results_dir = simulation_config.results_root,
        stop_time = simulation_config.stop_time,
    )

    simulation.callbacks[:progress] =
        Callback(progress, TimeInterval(simulation_config.progress_interval))

    ocean_sim = simulation.model.ocean
    ocean_model = ocean_sim.model
    ocean_sim.output_writers[:ocean] = NetCDFWriter(
        ocean_model,
        (
            T = ocean_model.tracers.T,
            S = ocean_model.tracers.S,
            u = ocean_model.velocities.u,
            v = ocean_model.velocities.v,
        );
        filename = results_path(simulation_config),
        schedule = TimeInterval(simulation_config.output_interval),
        overwrite_existing = simulation_config.overwrite_existing,
    )

    conjure_time_step_wizard!(
        simulation;
        cfl = simulation_config.time_step_cfl,
        max_Δt = simulation_config.max_time_step,
        max_change = simulation_config.max_time_step_change,
        cell_advection_timescale = cell_advection_timescale_coupled_model,
    )

    return simulation
end

"""
    run_simulation(config::FjordConfig)

Build the simulation a setup describes and run it to `stop_time`.

This is the setup-level driver, the same shape as `prepare_forcing(config::FjordConfig)`, and the
last step of the pipeline: it needs `prepare_bathymetry`, `prepare_forcing`, `add_rivers` if the
setup names rivers, and `prepare_atmosphere` if it names an atmosphere.

Returns `nothing` for a setup naming no simulation config, which is how `FjordSim.CLI.main`
reports a step a setup opts out of.
"""
function run_simulation(config::FjordConfig)
    simulation = build_simulation(config)
    isnothing(simulation) && return nothing

    output_file = results_path(config.simulation_config)
    run!(simulation)

    @info "Simulated $(prettytime(simulation.model.clock.time)) of model time"
    @info "Snapshots saved to $output_file"

    return (; simulation, output_file)
end

"""
    coupled_hydrostatic_simulation(grid, buoyancy, closure, ...)

Assemble a `HydrostaticFreeSurfaceModel` inside a NumericalEarth `OceanSeaIceModel` and return the
coupled `Simulation`.

The low-level entry point, taking every component already built: `build_simulation` is the one
that takes a `FjordConfig` and derives them.
"""
function coupled_hydrostatic_simulation(
    grid,
    buoyancy,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    initial_conditions,
    free_surface,
    coriolis,
    forcing,
    boundary_conditions,
    atmosphere,
    downwelling_radiation,
    sea_ice,
    biogeochemistry;
    results_dir = joinpath(homedir(), "FjordSim_results"),
    stop_time = 365days,
)
    isdir(results_dir) || mkpath(results_dir)

    println("Start compiling HydrostaticFreeSurfaceModel")
    ocean_model = HydrostaticFreeSurfaceModel(
        grid;
        buoyancy,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        free_surface,
        coriolis,
        forcing,
        boundary_conditions,
        biogeochemistry,
    )
    println("Done compiling HydrostaticFreeSurfaceModel")
    set!(ocean_model; initial_conditions...)
    Δt = 1second
    ocean_sim = Simulation(ocean_model; Δt, stop_time)
    coupled_model = OceanSeaIceModel(ocean_sim, sea_ice; atmosphere, radiation = downwelling_radiation)
    println("Initialized coupled model")
    coupled_simulation = Simulation(coupled_model; Δt, stop_time)
    return coupled_simulation
end  # function coupled_hydrostatic_simulation

end  # module Simulations
