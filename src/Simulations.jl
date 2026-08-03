module Simulations

export SimulationConfig,
    FromForcing,
    FromResults,
    simulation_architecture,
    build_simulation,
    run_simulation,
    coupled_hydrostatic_simulation

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using Oceananigans.TimeSteppers: reset!, update_state!
using NumericalEarth
using Dates: DateTime, Second
using NCDatasets

using ..Configs:
    AbstractSimulationConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    FjordConfig,
    bathymetry_path,
    forcing_path,
    river_forcing_path,
    results_path,
    run_tag
using ..Utils: recursive_merge, progress, cell_advection_timescale_coupled_model
using ..Atmospheres: prescribed_atmosphere, prescribed_radiation, atmosphere_date_range
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
`tracer_advection`, `momentum_advection`, `tracers`, `coriolis`, `sea_ice` and `biogeochemistry`
are grid-independent, so they are stored as the objects `coupled_hydrostatic_simulation` consumes;
`initial_conditions` sits with them but is resolved on the way through, since reading a state
needs the grid. The free surface, the boundary conditions, the forcing
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
- `buoyancy`, `closure`, `tracer_advection`, `momentum_advection`, `tracers`, `coriolis`,
  `sea_ice`, `biogeochemistry`: passed to `coupled_hydrostatic_simulation` unchanged.
- `initial_conditions`: where the ocean state starts from — a `NamedTuple` of constants, fields or
  functions, a `FromForcing`, or a `FromResults`. `build_simulation` puts it through
  `resolve_initial_conditions`, so what reaches `coupled_hydrostatic_simulation` is always a
  `NamedTuple` `set!` accepts.
- `free_surface_cfl`: `SplitExplicitFreeSurface` needs the grid, so only the CFL is stored.
- `bottom_drag_coefficient`: passed to `top_bottom_boundary_conditions`, which needs the grid.
- `start_date`: the calendar instant model time zero stands for.
- `stop_time`: seconds of simulated time one pass through the run window lasts.
- `loops`: how many times to run that window, carrying the ocean state across each restart. `1`
  runs it once. See `run_simulation`.
- `output_interval`, `overwrite_existing`: the snapshot writer's schedule and clobber policy.
- `progress_interval`: how often the `progress` callback reports.
- `checkpoint_interval`: how often to write a JLD2 checkpoint; `0.0` attaches no `Checkpointer`.
- `pickup`: resume from the newest checkpoint under `results_root` instead of starting fresh.
- `time_step_cfl`, `max_time_step`, `max_time_step_change`: the time-step wizard's settings.

A relative `output_file` resolves against `results_root`; an absolute one relocates just that
file. Either way `results_path` inserts `run_tag` — the `start_date` — before the extension, and a
looped run appends the loop index on top of that, so runs and repetitions land in separate files.

`start_date` is a field rather than something derived from an input file because every prepared
file has its *own* first record, and each reader used to zero its own axis there: the Oslofjord
forcing starts at 12:00 and its atmosphere at 00:00, which silently ran the two twelve hours out
of phase. One stated instant is the only thing they can all agree on, and
`validate_time_coverage` checks each file actually spans `[start_date, start_date + stop_time]`
rather than letting `Cyclical()` wrap the shortfall.

Every duration is in seconds, and must be written as a `Float64` — with `Oceananigans.Units`
(`365days`, `1hour`, `3minutes`), whose constants already are one. `Base.@kwdef` on a
*parametric* struct generates a constructor that does not convert, unlike the non-parametric
configs elsewhere in FjordSim, so `stop_time = 3600` is a `MethodError` where `stop_time = 1hour`
is fine, and `checkpoint_interval = 0` is one where `0.0` is fine. `loops` is the exception, being
an `Int` already.
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
    start_date::DateTime
    stop_time::Float64
    loops::Int
    output_interval::Float64
    progress_interval::Float64
    overwrite_existing::Bool
    checkpoint_interval::Float64
    pickup::Bool
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
The global attribute a snapshot file records its `start_date` in.

The snapshot writer's own time axis is seconds from model zero, which says nothing about the
calendar instant zero stood for — so a later `FromResults(path, date)` would have no way to turn a
date into a record. Writing the instant into the file it describes keeps that knowledge with the
data instead of making it a second config field that could disagree with the first.
"""
const RESULTS_START_DATE_ATTRIBUTE = "start_date"

"""
    FromForcing(date = nothing)

Initial conditions read from the prepared forcing file the simulation is already reading, at
`date` — or at the simulation's `start_date` when `date` is `nothing`.

The file is on the model grid by construction (`prepare_forcing` regrids onto it), so this is a
plain read: no regridding, no inpainting, no `NumericalEarth` dataset wrapper.
"""
struct FromForcing{D}
    date::D
end

FromForcing() = FromForcing(nothing)

"""
    FromResults(path, date = nothing)

Initial conditions read from a previous run's snapshot file, at `date` — or from its last record
when `date` is `nothing`.

A relative `path` resolves against `results_root`, like `output_file` does. Naming a `date` requires
the file to carry the `$RESULTS_START_DATE_ATTRIBUTE` attribute `build_simulation` writes, since a
snapshot's time axis is seconds from its own model zero; files written before that attribute existed
can only be read by their last record.
"""
struct FromResults{D}
    path::String
    date::D
end

FromResults(path::String) = FromResults(path, nothing)

"""
    initial_conditions_date(initial_conditions, start_date)

Which instant to read, defaulting an unnamed date to the run's own `start_date`.
"""
initial_conditions_date(::FromForcing{Nothing}, start_date) = start_date
initial_conditions_date(initial_conditions::FromForcing, start_date) = initial_conditions.date

"""
    resolve_initial_conditions(initial_conditions, grid, forcing_file, config)

Turn a setup's `initial_conditions` into the `NamedTuple` `set!(model; ...)` consumes.

Dispatched on the kind of source, so `coupled_hydrostatic_simulation` never learns that there is
more than one: it still receives something splattable and still applies it with a single `set!`.
A `NamedTuple` — constants, functions or fields — passes straight through, which is what every
setup did before the other two existed.

What gets set is every tracer `config.tracers` names plus `u` and `v`, intersected with what the
source file actually carries — see `state_variables`. Nothing is enumerated here, so adding a
biogeochemical tracer to a setup is enough to have it read back.

The free surface `η` and any closure-owned tracer the config does not name (CATKE's `e`) keep their
defaults, so reading a state in is a *warm start*, not a restart: the barotropic mode and the
turbulence field re-adjust over the first hours. Use `pickup` for an exact continuation.
"""
resolve_initial_conditions(initial_conditions::NamedTuple, grid, forcing_file, config) =
    initial_conditions

function resolve_initial_conditions(initial_conditions::FromForcing, grid, forcing_file, config)
    date = initial_conditions_date(initial_conditions, config.start_date)
    @info "Initial conditions from forcing $forcing_file at $date"
    return forcing_state(forcing_file, grid, date, config.tracers)
end

function resolve_initial_conditions(initial_conditions::FromResults, grid, forcing_file, config)
    filepath = results_state_path(initial_conditions, config)
    isfile(filepath) || error("Results file $filepath does not exist.")
    @info "Initial conditions from results $filepath"
    return results_state(filepath, grid, initial_conditions.date, config.tracers)
end

"""
    results_state_path(initial_conditions, config)

Where `FromResults` reads from: its `path` as given when absolute, resolved against `results_root`
when relative — the same rule `output_file` follows, so a previous run's output can be named by its
filename alone.
"""
results_state_path(initial_conditions::FromResults, config) =
    isabspath(initial_conditions.path) ? initial_conditions.path :
    joinpath(config.results_root, initial_conditions.path)

"""
    forcing_state(filepath, grid, date, tracers)

The state variables a prepared forcing file holds at `date`.

The `_lambda` twins are relaxation rates rather than state, so only the bare names are read.
"""
function forcing_state(filepath, grid, date, tracers)
    return NCDataset(filepath) do ds
        validate_state_dimensions(filepath, ds, grid, ("Nx", "Ny", "Nz"))
        dates = ds["time"][:]
        index = findfirst(==(date), dates)
        isnothing(index) && error(
            "No forcing record at $date in $filepath, whose axis runs $(first(dates)) to " *
            "$(last(dates)). Name a date on that axis, or prepare forcing covering $date.",
        )
        return state_variables(ds, index, eltype(grid), tracers)
    end
end

"""
    results_state(filepath, grid, date, tracers)

The state variables a previous run's snapshot file holds at `date`, or at its last record when
`date` is `nothing`.
"""
function results_state(filepath, grid, date, tracers)
    return NCDataset(filepath) do ds
        validate_state_dimensions(filepath, ds, grid, ("λ_caa", "φ_aca", "z_aac"))
        return state_variables(ds, results_record_index(filepath, ds, date), eltype(grid), tracers)
    end
end

"""
    results_record_index(filepath, ds, date)

Which record of a snapshot file `date` names. Its time axis is seconds from the run's own model
zero, so the instant that zero stood for has to come from the file's own
`$RESULTS_START_DATE_ATTRIBUTE` attribute; `nothing` takes the last record and needs no attribute.
"""
results_record_index(filepath, ds, ::Nothing) = ds.dim["time"]

function results_record_index(filepath, ds, date::DateTime)
    haskey(ds.attrib, RESULTS_START_DATE_ATTRIBUTE) || error(
        "$filepath carries no `$RESULTS_START_DATE_ATTRIBUTE` attribute, so its time axis " *
        "(seconds from that run's model zero) cannot be turned into a date. Read its last record " *
        "with `FromResults(path)` instead, or re-run the simulation that wrote it.",
    )

    source_start = DateTime(ds.attrib[RESULTS_START_DATE_ATTRIBUTE])
    seconds = Second(date - source_start).value
    times = ds["time"][:]
    index = argmin(abs.(times .- seconds))

    abs(times[index] - seconds) <= 1 || error(
        "No record within 1 s of $date in $filepath, whose axis runs $source_start to " *
        "$(source_start + Second(round(Int, last(times)))).",
    )

    return index
end

"""
    state_variables(ds, index, FT, tracers)

The state variables `ds` carries at time `index`, as a `NamedTuple` of `FT` arrays.

Which variables those are comes from the simulation config's `tracers` plus the two horizontal
velocities — never a list written out here, so a setup that adds a biogeochemical tracer gets it read
back without this module being touched. The rule is
`(map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)`, the same one `forcing_from_file` uses to decide
which forcing terms to build, so the two cannot disagree about what the state is.

Intersecting with the file's own variables is what lets one reader serve a forcing file and a
snapshot file: they need not carry the same set, and a tracer the source lacks is simply left at its
default rather than being an error.
"""
function state_variables(ds, index, FT, tracers)
    names = (map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)
    return (; (Symbol(name) => finite_slab(ds[name][:, :, :, index], FT) for name in names)...)
end

"""
    finite_slab(data, FT)

`data` as `FT` with every missing or non-finite cell replaced by zero.

Both sources mark land that way, and in both it is exactly the cells this grid immerses — the
forcing file's mask comes from `peripheral_node` on this very grid, and a snapshot was written from
a field on it — so zeroing them fills only cells the model never reads.

`FT` is the grid's element type, not the file's: `set!` moves the array to the model's architecture
and `copyto!`s it into the field, and a `Union{Missing,Float32}` array cannot become a `CuArray` at
all while a mismatched element type need not convert on the device.
"""
finite_slab(data, FT) =
    map(value -> ismissing(value) || !isfinite(value) ? zero(FT) : convert(FT, value), data)

"""
    validate_state_dimensions(filepath, ds, grid, names)

Check the file's horizontal and vertical extents are the grid's, so a state file from another setup
is a clear error rather than a silent misread.
"""
function validate_state_dimensions(filepath, ds, grid, names)
    expected = size(grid)
    found = ntuple(index -> ds.dim[names[index]], 3)
    found == expected || throw(
        DimensionMismatch("$filepath is $found but the simulation grid is $expected"),
    )
    return nothing
end

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
    validate_time_coverage(label, range, start_date, stop_time, remedy)

Check a prepared file's date `range` contains the run's whole interval.

Both readers use `Cyclical()` time indexing, which does not fail outside the data it was given —
it wraps, so a run that outlasts its forcing quietly replays the beginning and a run that starts
before it quietly reads the end. This is what makes that unreachable instead of merely unlikely.
A `nothing` range is a source that cannot report its dates, and is skipped.
"""
validate_time_coverage(label, ::Nothing, start_date, stop_time, remedy) = nothing

function validate_time_coverage(label, range, start_date, stop_time, remedy)
    first_available, last_available = range
    end_date = start_date + Second(round(Int, stop_time))

    start_date >= first_available || error(
        "The run starts at $start_date but the $label only begins at $first_available. Move the " *
        "simulation config's `start_date` to $first_available or later, or $remedy.",
    )
    end_date <= last_available || error(
        "The run ends at $end_date but the $label stops at $last_available. Shorten `stop_time` " *
        "to $(Second(last_available - start_date).value) seconds or less, or $remedy.",
    )

    return nothing
end

"""
    forcing_date_range(filepath)

First and last date on a prepared forcing file's time axis.
"""
forcing_date_range(filepath) = NCDataset(filepath) do ds
    dates = ds["time"][:]
    (first(dates), last(dates))
end

"""
    loop_output_path(config, loop)

The snapshot file for one repetition: the plain run-tagged name for a setup that runs its window
once, and the loop-indexed one when it repeats. A single run therefore keeps the shorter name
rather than gaining a `_loop01` nobody asked for.
"""
loop_output_path(config::AbstractSimulationConfig, loop) =
    config.loops == 1 ? results_path(config) : results_path(config, loop)

"""
    checkpoint_prefix(config, loop)

Prefix for one repetition's checkpoints, under which `Checkpointer` writes
`<prefix>_iteration<N>.jld2`.

The loop index is in the *name* because it is not in the checkpoint: the state records the clock but
not which repetition produced it, so without this `pickup` could not tell a loop-3 checkpoint from a
loop-1 one and would replay the whole spin-up. `resume_loop` reads the index back out.
"""
checkpoint_prefix(config::AbstractSimulationConfig, loop) =
    string("checkpoint_", run_tag(config), "_loop", lpad(loop, 2, '0'))

"""
    checkpointed_loops(config)

Every loop index that has a checkpoint under `results_root` for this run tag.
"""
function checkpointed_loops(config::AbstractSimulationConfig)
    isdir(config.results_root) || return Int[]
    pattern = Regex(string("^checkpoint_", run_tag(config), "_loop(\\d+)_iteration\\d+\\.jld2\$"))
    found = filter(!isnothing, match.(pattern, readdir(config.results_root)))
    return [parse(Int, captured[1]) for captured in found]
end

"""
    resume_loop(config)

Which repetition to start at: the highest one that has a checkpoint when `pickup` is set, and `1`
otherwise.

`build_simulation` attaches its writers for this loop rather than unconditionally for the first, so
a resumed run writes into the file it was already writing.
"""
function resume_loop(config::AbstractSimulationConfig)
    config.pickup || return 1

    loops = checkpointed_loops(config)
    if isempty(loops)
        @warn "`pickup` is set but no checkpoint for this run was found in " *
              "$(config.results_root); starting from the beginning"
        return 1
    end

    return maximum(loops)
end

"""
    attach_writers!(simulation, config, loop; picking_up = false)

Point the snapshot writer, and the checkpointer if there is one, at `loop`'s own files.

Called once by `build_simulation` and again by `restart_loop!` for each subsequent repetition. A
*fresh* writer each time rather than a renamed one, because a `TimeInterval` accumulates its
actuation count: reused after a clock reset it would believe it was thousands of records ahead and
never fire again. The old one is closed first — it holds an open `NCDataset`, so replacing it
silently would leak a handle and an unflushed tail per loop.

The snapshot file records `start_date` as a global attribute so `FromResults(path, date)` can later
turn a date into a record — see `RESULTS_START_DATE_ATTRIBUTE`.

`picking_up` suppresses `overwrite_existing`, because a checkpoint restores this writer's actuation
count: clobbering the file it is meant to continue would leave a schedule thousands of records ahead
of an empty one.
"""
function attach_writers!(
    simulation,
    config::AbstractSimulationConfig,
    loop;
    picking_up = false,
)
    ocean_sim = simulation.model.ocean
    ocean_model = ocean_sim.model

    haskey(ocean_sim.output_writers, :ocean) &&
        close(pop!(ocean_sim.output_writers, :ocean))

    ocean_sim.output_writers[:ocean] = NetCDFWriter(
        ocean_model,
        (
            T = ocean_model.tracers.T,
            S = ocean_model.tracers.S,
            u = ocean_model.velocities.u,
            v = ocean_model.velocities.v,
        );
        filename = loop_output_path(config, loop),
        schedule = TimeInterval(config.output_interval),
        overwrite_existing = config.overwrite_existing && !picking_up,
        global_attributes = Dict(RESULTS_START_DATE_ATTRIBUTE => string(config.start_date)),
    )

    return attach_checkpointer!(simulation, config, loop)
end

"""
    attach_checkpointer!(simulation, config, loop)

Attach a `Checkpointer` for `loop`, or nothing at all when `checkpoint_interval` is zero.

It goes on the *coupled* simulation, not the ocean one, for two reasons that both fail otherwise:
`prognostic_state` of the `OceanSeaIceModel` is what a resumable state actually is, and
`run!(…; pickup)` looks for its checkpointer in `simulation.output_writers`.

`NumericalEarth` keeps this cheap: a coupled checkpoint holds the ocean's prognostic fields, CATKE's
diffusivities and every component's clock, but the prescribed atmosphere and radiation contribute
only their clock. So none of the things that cannot be serialized — `ForcingFromFile` and its
`FieldTimeSeries`, the `FreshwaterExchange` in the tracer top boundary conditions — is ever written.
"""
function attach_checkpointer!(simulation, config::AbstractSimulationConfig, loop)
    config.checkpoint_interval > 0 || return simulation

    simulation.output_writers[:checkpointer] = Checkpointer(
        simulation.model;
        schedule = TimeInterval(config.checkpoint_interval),
        dir = config.results_root,
        prefix = checkpoint_prefix(config, loop),
        overwrite_existing = config.overwrite_existing,
        cleanup = true,
    )

    return simulation
end

"""
    rewind_clock!(component)

Send one coupled-model component's clock back to zero, or do nothing for a component that has no
clock.

`NumericalEarth`'s own `reset_clock!(::EarthSystemModel)` cannot be used for this. Its per-component
fallback is `reset!(getproperty(component, :clock))` and `components` includes `sea_ice`, so a
`FreezingLimitedOceanTemperature` — which is a liquidus and nothing else — makes it throw
`has no field clock`. Dispatching on whether the component *has* a clock keeps this module free of
any named component, so a setup that later names a real sea-ice or land model works unchanged, and
`Val` resolves the question at compile time.
"""
rewind_clock!(::Nothing) = nothing
rewind_clock!(component::Simulation) = rewind_clock!(component.model)
rewind_clock!(component) = rewind_clock!(component, Val(hasfield(typeof(component), :clock)))
rewind_clock!(component, ::Val{true}) = reset!(component.clock)
rewind_clock!(component, ::Val{false}) = nothing

"""
    restart_loop!(simulation, config, loop)

Send `simulation` back to the start of its window for repetition `loop`, keeping the ocean state.

Every clock in the coupled model goes back to zero — the coupled model's own, and each component's,
found by sweeping its properties rather than naming them. Nothing else is reset: the ocean's
velocities, tracers, CATKE diffusivities and free surface all carry over, which is the whole point.
`stop_time` is untouched, and so is `simulation.Δt`, so the time-step wizard keeps the step it
converged on instead of re-ramping from one second.

Two things then have to be prodded by hand, and both are silent if they are not.

The prescribed atmosphere and radiation hold their `FieldTimeSeries` window, and rewinding a clock
does not refill it. `time_step!(::EarthSystemModel)` assembles surface fluxes in
`maybe_prepare_first_time_step!` *before* it steps the atmosphere, so without this the first step of
each loop would be forced by last December — which is exactly why `NumericalEarth`'s own
`reset_clock!` ends with an `update_state!` on the atmosphere.

And the ocean is a `Simulation` of its own, which `run!` never touches: `run!` clears `initialized`
on the coupled simulation only, so `time_step!(ocean_sim)` would skip `initialize!` and the fresh
snapshot writer would never have its schedule initialized or its `t = 0` record written.

`clock.last_Δt` comes back as `Inf`, which makes the first step of each loop a forward Euler step.
That is both harmless and right: `G⁻` refers to a step taken under forcing a year away.
"""
function restart_loop!(simulation, config::AbstractSimulationConfig, loop)
    model = simulation.model

    rewind_clock!(model)
    for name in propertynames(model)
        rewind_clock!(getproperty(model, name))
    end

    update_state!(model.atmosphere)
    update_state!(model.radiation)

    model.ocean.initialized = false

    return attach_writers!(simulation, config, loop)
end

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

The initial conditions go through `resolve_initial_conditions`, so a `FromForcing` or `FromResults`
is read here — where the grid and the forcing path are known — and what
`coupled_hydrostatic_simulation` receives is always a plain `NamedTuple`. Note that `pickup`
supersedes them: the checkpoint restores the state that `set!` had just written.

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

    simulation_config.loops >= 1 ||
        throw(ArgumentError("loops must be at least 1, got $(simulation_config.loops)"))
    simulation_config.checkpoint_interval >= 0 || throw(
        ArgumentError(
            "checkpoint_interval must not be negative, got " *
            "$(simulation_config.checkpoint_interval)",
        ),
    )

    start_date = simulation_config.start_date
    stop_time = simulation_config.stop_time
    # One window, not `loops` of them: every repetition replays the same interval, so what the
    # prepared files have to span does not grow with the loop count.
    validate_time_coverage(
        "forcing $forcing_file",
        forcing_date_range(forcing_file),
        start_date,
        stop_time,
        "prepare more years of forcing",
    )
    validate_time_coverage(
        "atmosphere",
        atmosphere_date_range(config.atmosphere_config),
        start_date,
        stop_time,
        "prepare more years of atmosphere",
    )

    architecture = simulation_architecture(simulation_config)
    grid = ImmersedBoundaryGrid(bathymetry_file, architecture, config.grid_config.halo)

    # Every time axis is zeroed at the same instant, so the components stay in phase: each
    # prepared file has its own first record, and left to itself each reader would zero there.
    forcing = forcing_from_file(;
        grid = grid,
        filepath = forcing_file,
        tracers = simulation_config.tracers,
        reference_date = start_date,
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

    initial_conditions = resolve_initial_conditions(
        simulation_config.initial_conditions,
        grid,
        forcing_file,
        simulation_config,
    )

    simulation = coupled_hydrostatic_simulation(
        grid,
        simulation_config.buoyancy,
        simulation_config.closure,
        simulation_config.tracer_advection,
        simulation_config.momentum_advection,
        simulation_config.tracers,
        initial_conditions,
        SplitExplicitFreeSurface(grid, cfl = simulation_config.free_surface_cfl),
        simulation_config.coriolis,
        forcing,
        boundary_conditions,
        prescribed_atmosphere(config.atmosphere_config, architecture; reference_date = start_date),
        prescribed_radiation(config.atmosphere_config, architecture; reference_date = start_date),
        simulation_config.sea_ice,
        simulation_config.biogeochemistry;
        results_dir = simulation_config.results_root,
        stop_time,
    )

    simulation.callbacks[:progress] =
        Callback(progress, TimeInterval(simulation_config.progress_interval))

    # For the loop `run_simulation` will start at, not unconditionally the first: a resumed run has
    # to keep writing into the file it was already writing, and must not truncate it.
    start_loop = resume_loop(simulation_config)
    attach_writers!(
        simulation,
        simulation_config,
        start_loop;
        picking_up = simulation_config.pickup,
    )

    @info "Run $(run_tag(simulation_config)): snapshots to " *
          "$(loop_output_path(simulation_config, start_loop))"

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

Build the simulation a setup describes and run its window `loops` times.

This is the setup-level driver, the same shape as `prepare_forcing(config::FjordConfig)`, and the
last step of the pipeline: it needs `prepare_bathymetry`, `prepare_forcing`, `add_rivers` if the
setup names rivers, and `prepare_atmosphere` if it names an atmosphere.

# Looping

Each repetition replays `[start_date, start_date + stop_time]` with the ocean state carried over, so
`loops > 1` is a spin-up: one forcing year run again and again until the deep basins stop drifting.
`restart_loop!` does the carrying over — the clock goes back to zero, the state does not — and each
repetition writes its own file, so the loops can be compared rather than overwriting one another.

The clock is reset rather than left to run on for `loops * stop_time` because that keeps the run
inside `[t¹, tᴺ]` of every prepared file. Both readers use `Cyclical()` time indexing, and a
monotonic clock would lean on its wrap-around, whose period is inferred from the file's own axis and
so stops matching the loop the moment a file is padded to a different window.

Returns `nothing` for a setup naming no simulation config, which is how `FjordSim.CLI.main`
reports a step a setup opts out of; otherwise `(; simulation, output_files)`.
"""
function run_simulation(config::FjordConfig)
    simulation = build_simulation(config)
    isnothing(simulation) && return nothing

    simulation_config = config.simulation_config
    # Asked again rather than threaded out of `build_simulation`, whose contract is the bare
    # `Simulation`: it is a directory scan, and both entry points have to agree on the answer
    # anyway. The one visible cost is that a `pickup` with no checkpoints warns twice.
    start_loop = resume_loop(simulation_config)
    output_files = String[]

    for loop = start_loop:simulation_config.loops
        loop > start_loop && restart_loop!(simulation, simulation_config, loop)

        output_file = loop_output_path(simulation_config, loop)
        push!(output_files, output_file)
        @info "Loop $loop of $(simulation_config.loops): snapshots to $output_file"

        # `pickup` applies to the loop we resume into and to no other: the loops after it start from
        # a reset clock, and picking up there would send them back to the checkpoint every time.
        run!(
            simulation;
            pickup = simulation_config.pickup && loop == start_loop,
            checkpoint_at_end = simulation_config.checkpoint_interval > 0,
        )

        @info "Loop $loop reached $(prettytime(simulation.model.clock.time)) of model time"
    end

    @info "Snapshots saved to $(join(output_files, ", "))"

    return (; simulation, output_files)
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
