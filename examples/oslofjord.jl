# Oslofjord variant using an implicit free surface and a true (radiating) open boundary condition
# on velocity, instead of the split-explicit free surface and closed-wall lateral velocity
# `oslofjorden()` (`src/Setups/oslofjorden.jl`) uses. See that setup for the rationale behind every
# field this one copies verbatim.
#
# This is the stripped-down, idealized form of that variant: forcing, bathymetry and the boundary
# momentum sponge are all reduced to the minimum needed to test one thing — whether the open
# boundary + `ImplicitFreeSurface` + `ZStarCoordinate` combination is numerically stable on its own,
# with every other confound removed. Concretely:
#
#   - Ocean forcing is disabled entirely (`forcing_config = nothing`, a first-class option
#     documented on `FjordConfig`) — no NorKyst data, no rivers. The open boundary instead relaxes
#     towards a quiescent exterior (zero velocity, the run's own initial T/S) via plain constants.
#   - Bathymetry is flattened to one constant deep depth everywhere at sea (`flatten_bathymetry!`),
#     keeping only the land/sea mask from the real processed file, with the open edge's boundary
#     band forced to water regardless of what the mask said there.
#   - The boundary momentum sponge tried in an earlier pass of this variant is removed: tested
#     empirically, it did not help stability (Oceananigans' own closest reference case,
#     `validation/open_boundaries/barotropic_soliton.jl`, does use a sponge for its own 150-day run,
#     but that is not carried over here).
#
# Everything new is defined below, entirely outside the package, following FjordSim's own "Adding a
# new source" extension model — subtype a supertype, overload its hook:
#
#   - `ImplicitFreeSurfaceConfig` wraps Oceananigans' `ImplicitFreeSurface`, the free-surface family
#     open-boundary support (github.com/CliMA/Oceananigans.jl/issues/5229) is being built against
#     first, since its predictor/corrector structure mirrors `NonhydrostaticModel`'s pressure
#     correction.
#   - `RadiatingLateralBoundary` is a genuinely open lateral boundary — a radiating
#     `NormalFlowBoundaryCondition`/`ValueBoundaryCondition` pair using Oceananigans'
#     `PerturbationAdvection` scheme, matching how Oceananigans' own validation examples build a
#     constant/analytic open boundary (`validation/open_boundaries/flow_over_hill.jl`'s
#     `NormalFlowBoundaryCondition(U; scheme)`, `tracer_outflow.jl`'s
#     `ValueBoundaryCondition(c̄; scheme)`) rather than FjordSim's closed-wall `OpenLateralBoundary`.
#   - `OpenBoundaryHydrostaticSimulation` replaces `CoupledHydrostaticSimulation`: it builds the
#     `HydrostaticFreeSurfaceModel` on a grid with a *mutable* vertical coordinate
#     (`vertical_coordinate = ZStarCoordinate()`) instead of FjordSim's usual static `z`.
#
# That last piece is the one this variant's first pass was missing, and why it diverged (a
# `DomainError` inside the `ImplicitFreeSurface` PCG solver) once actually run end-to-end rather than
# just built. Both of Oceananigans' own reference implementations of this exact combination —
# `HydrostaticFreeSurfaceModel` + `ImplicitFreeSurface` + a genuinely open lateral velocity —
# (`validation/open_boundaries/flow_over_hill.jl` and `validation/open_boundaries/barotropic_soliton.jl`,
# both merged against #5229 well before the installed Oceananigans v0.110.13) run on a
# `MutableVerticalDiscretization` grid, not a static one: on a static grid, the free-surface
# height η is solved but never fed back into cell thickness, so a persistent one-directional inflow
# at the open boundary has no way to reconcile into the water column and the solver eventually goes
# indefinite. `ZStarCoordinate` is what makes η change the column depth those cells actually
# integrate over — see `ImmersedBoundaries/mutable_immersed_grid.jl`'s `column_depth` methods.
#
# `target_transport` (also new in #5229's PRs) is deliberately not used here: it only matters for
# `NonhydrostaticModel`'s net-zero-flux solvability condition across *paired* open boundaries (see
# `flow_over_hill.jl`, which has west+east, and `targeted_flux_xy.jl`, whose real invariant is that
# at least one *untargeted* boundary must remain to absorb the net imbalance). This domain has
# exactly one open boundary and no other named lateral boundary at all, and
# `HydrostaticFreeSurfaceModels/boundary_targeted_transport.jl` states plainly that the free surface
# is free to absorb any net imbalance on its own — which is exactly what `ZStarCoordinate` makes
# possible.
#
# FjordSim's own grid loader (`ImmersedBoundaryGrid(filepath, architecture, halo)` in `src/Grids.jl`)
# is not a hook: `build_simulation` calls it directly, independent of the model config, so there is
# no existing extension point to change the grid's vertical coordinate for one variant. Rather than
# touch `src/Grids.jl`, `OpenBoundaryHydrostaticSimulation`'s `coupled_simulation` method (the one
# hook `AbstractCoupledSimulationConfig` already grants for grid-dependent model assembly) rebuilds an
# equivalent grid itself, with the example's own `EvenGrid` fields for geometry and the already-built
# static grid's land/sea mask carried over (then flattened) — see `implicit_free_surface_grid` below.
#
# This file is itself a config, not a runner — it reads the same prepared atmosphere
# `oslofjorden()`'s prepare steps write (same `data_root`), so those steps have to have run first:
#
#   julia --project -m FjordSim download_atmosphere --config oslofjorden
#   julia --project -m FjordSim prepare_atmosphere  --config oslofjorden
#
# Bathymetry is *not* shared with `oslofjorden()`, on purpose: `bathymetry_config` points at
# `bathymetry_boundary.nc`, a one-off copy of `oslofjorden()`'s processed bathymetry, prepared with
# its own call:
#
#   julia --project -m FjordSim prepare_bathymetry --config examples/oslofjord.jl
#
# No forcing preparation is needed at all — `forcing_config = nothing` below means neither
# `download_forcing`, `prepare_forcing` nor `add_rivers` has anything to do for this variant.
#
# Then run this variant by pointing `--config` at this file directly:
#
#   julia --project -m FjordSim run_simulation --config examples/oslofjord.jl
#
# To step through the assembly instead of running it, build the simulation without starting it:
#
#   using FjordSim
#   using Oceananigans: run!
#
#   config = fjord_config("examples/oslofjord.jl")
#   simulation = build_simulation(config)
#   run!(simulation)
#

using FjordSim
using Oceananigans
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units
using Dates: DateTime
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature, OceanSeaIceModel

# --- ImplicitFreeSurfaceConfig: a new AbstractFreeSurfaceConfig ------------------------------

# Holds only the knobs, not the object — `ImplicitFreeSurface(; solver_method, ...)` needs no grid,
# but `free_surface` still takes one for symmetry with every other `AbstractFreeSurfaceConfig`
# method (`SplitExplicitFreeSurfaceConfig`'s does need it). `reltol`/`abstol`/`maxiter` are stated
# explicitly rather than left at Oceananigans' PCG defaults — but *not* tightened to
# `validation/open_boundaries/barotropic_soliton.jl`'s `1e-10`, which was tried first and traced
# (by instrumenting the solver's own iteration count and residual right up to the failure) to be
# the actual cause of the `DomainError`: the CG solver was converging cleanly in a steady 9
# iterations every step, nowhere near `maxiter`, with no anomalous growth in `u`/`v`/`η` or column
# depth beforehand — the signature of a tolerance pushed past what this larger, immersed-boundary
# grid's conditioning can resolve in Float64, so the solver's internal inner products occasionally
# round off to a value just below zero where an exact non-negative one was expected, and `sqrt` of
# it throws. `1e-7` matches `PCGImplicitFreeSurfaceSolver`'s own computed default
# (`min(1e-7, 10 * sqrt(eps(Float64)))`) rather than guessing a number.
struct ImplicitFreeSurfaceConfig <: AbstractFreeSurfaceConfig
    solver_method::Symbol
    gravitational_acceleration::Float64
    reltol::Float64
    abstol::Float64
    maxiter::Int
end

ImplicitFreeSurfaceConfig(; solver_method, gravitational_acceleration, reltol, abstol, maxiter) =
    ImplicitFreeSurfaceConfig(
        solver_method,
        Float64(gravitational_acceleration),
        Float64(reltol),
        Float64(abstol),
        Int(maxiter),
    )

# Extends FjordSim's exported `free_surface` hook on a type this file owns — not type piracy, the
# same pattern `src/Forcing/norkyst.jl` etc. use for their own hooks from inside the package.
FjordSim.free_surface(config::ImplicitFreeSurfaceConfig, grid) = ImplicitFreeSurface(
    solver_method              = config.solver_method,
    gravitational_acceleration = config.gravitational_acceleration,
    reltol                     = config.reltol,
    abstol                     = config.abstol,
    maxiter                    = config.maxiter,
)

# --- RadiatingLateralBoundary: a new AbstractBoundaryConditionConfig -------------------------

# A quiescent, forcing-free open boundary: velocity relaxes towards `exterior_velocity` (a still
# exterior ocean, `0.0`) and each tracer towards its own constant in `exterior_tracers`, both via
# plain `PerturbationAdvection`-scheme conditions — no `discrete_form`, no `FieldTimeSeries`,
# matching how Oceananigans' own validation examples build a constant/analytic open boundary
# (`flow_over_hill.jl`'s `NormalFlowBoundaryCondition(U; scheme)`, `tracer_outflow.jl`'s
# `ValueBoundaryCondition(c̄; scheme)`). An earlier pass of this variant read real, time-varying
# NorKyst data at the exact boundary column instead, which needed `discrete_form` and a
# `FieldTimeSeries`; with forcing disabled (`forcing_config = nothing`) there is nothing
# time-varying left to read, so that whole apparatus is gone.
#
# Carries its own `edge`/`relaxation_timescale` rather than reading a forcing config's, since
# `forcing_config` is `nothing` here. No `target_transport`: that only matters for balancing
# *paired* open boundaries, and this domain has one.
struct RadiatingLateralBoundary <: AbstractBoundaryConditionConfig
    edge::Symbol
    relaxation_timescale::Float64
    exterior_velocity::Float64
    exterior_tracers::NamedTuple
end

RadiatingLateralBoundary(; edge, relaxation_timescale, exterior_velocity, exterior_tracers) =
    RadiatingLateralBoundary(edge, Float64(relaxation_timescale), Float64(exterior_velocity), exterior_tracers)

open_boundary_scheme(config::RadiatingLateralBoundary) = PerturbationAdvection(;
    inflow_timescale = config.relaxation_timescale,
    outflow_timescale = config.relaxation_timescale,
)

# The open-velocity counterpart of FjordSim's internal `open_boundary_conditions`, towards a
# constant exterior value rather than the closed wall (`NormalFlowBoundaryCondition(nothing)`) that
# builds.
open_velocity_boundary_conditions(::Val{:south}, config::RadiatingLateralBoundary) =
    (v = (south = NormalFlowBoundaryCondition(config.exterior_velocity; scheme = open_boundary_scheme(config)),),)
open_velocity_boundary_conditions(::Val{:north}, config::RadiatingLateralBoundary) =
    (v = (north = NormalFlowBoundaryCondition(config.exterior_velocity; scheme = open_boundary_scheme(config)),),)
open_velocity_boundary_conditions(::Val{:west}, config::RadiatingLateralBoundary) =
    (u = (west = NormalFlowBoundaryCondition(config.exterior_velocity; scheme = open_boundary_scheme(config)),),)
open_velocity_boundary_conditions(::Val{:east}, config::RadiatingLateralBoundary) =
    (u = (east = NormalFlowBoundaryCondition(config.exterior_velocity; scheme = open_boundary_scheme(config)),),)
open_velocity_boundary_conditions(::Val{edge}, ::RadiatingLateralBoundary) where {edge} =
    throw(ArgumentError("edge must be one of (:south, :north, :west, :east), got :$edge"))

# One tracer's constant-relaxing open boundary condition, the tracer counterpart of
# `open_velocity_boundary_conditions` above.
function tracer_boundary_condition(config::RadiatingLateralBoundary, name)
    haskey(config.exterior_tracers, name) || error(
        "edge names an open tracer boundary but exterior_tracers does not name `$name`; add it to " *
        "`RadiatingLateralBoundary`'s `exterior_tracers`.",
    )
    return ValueBoundaryCondition(getproperty(config.exterior_tracers, name); scheme = open_boundary_scheme(config))
end

open_tracer_boundary_conditions(::Val{edge}, config::RadiatingLateralBoundary, tracer_names) where {edge} =
    NamedTuple(
        name => (; Symbol(edge) => tracer_boundary_condition(config, name))
        for name in tracer_names
    )

# Extends FjordSim's `boundary_conditions` hook, which is deliberately not re-exported (it would
# shadow `Oceananigans.Fields.boundary_conditions`) — reached the way `CLAUDE.md`'s Import
# Conventions section prescribes for extending a function you don't own: qualify the module, never
# `import Mod: foo`. `field_boundary_conditions`'s internal call to `boundary_conditions(config, ...)`
# dispatches to this method exactly as it does to the built-in ones, since a generic function has one
# global method table regardless of which module adds a method to it. `forcing`/`forcing_config` are
# unused — this variant runs with `forcing_config = nothing`, the same pattern `TopBottomFluxes`
# already uses for the parameters it doesn't need.
FjordSim.BoundaryConditions.boundary_conditions(config::RadiatingLateralBoundary, grid, forcing, forcing_config, tracers) =
    recursive_merge(
        open_velocity_boundary_conditions(Val(config.edge), config),
        open_tracer_boundary_conditions(Val(config.edge), config, tracers),
    )

# --- OpenBoundaryHydrostaticSimulation: a new AbstractCoupledSimulationConfig ------------------

# Flattens every wet cell in the whole domain to one constant `depth`, leaving land
# (`bottom_height >= 0`) untouched, then unconditionally forces the `boundary_rows`-wide band along
# `edge` to that same depth too — including any land there, so the open edge always has a full
# water boundary to work with regardless of what the real bathymetry's mask said. Dispatch on
# `Val(edge)` mirrors the removed `boundary_sponge_mask`.
#
# Isolates the open-boundary/`ImplicitFreeSurface`/`ZStarCoordinate` mechanism from bathymetry-driven
# artifacts entirely, rather than only floor the boundary row as an earlier pass of this variant
# did: a uniform column depth everywhere at sea means every wet cell has the same inertia and the
# same number of active vertical layers, so an instability can no longer be attributed to how deep
# or shallow any particular cell happens to be.
function flatten_bathymetry!(bottom, ::Val{:south}, boundary_rows, depth)
    bottom .= ifelse.(bottom .< 0, -depth, bottom)
    bottom[:, 1:boundary_rows] .= -depth
    return bottom
end
function flatten_bathymetry!(bottom, ::Val{:north}, boundary_rows, depth)
    bottom .= ifelse.(bottom .< 0, -depth, bottom)
    bottom[:, end-boundary_rows+1:end] .= -depth
    return bottom
end
function flatten_bathymetry!(bottom, ::Val{:west}, boundary_rows, depth)
    bottom .= ifelse.(bottom .< 0, -depth, bottom)
    bottom[1:boundary_rows, :] .= -depth
    return bottom
end
function flatten_bathymetry!(bottom, ::Val{:east}, boundary_rows, depth)
    bottom .= ifelse.(bottom .< 0, -depth, bottom)
    bottom[end-boundary_rows+1:end, :] .= -depth
    return bottom
end

# `implicit_free_surface_grid` rebuilds an equivalent grid with a mutable vertical coordinate from
# the ordinary static one `build_simulation` already built. Geometry (`size`, `halo`, `longitude`,
# `latitude`, `z_faces`) comes from the example's own `grid_config`, not by re-deriving it from
# `static_grid` — `src/Grids.jl`'s own doc comment flags exactly this face-vs-center count as an
# easy off-by-one. Only the land/sea mask is carried over from `static_grid`'s real, processed
# bathymetry; every wet depth is then overwritten by `flatten_bathymetry!` at the grid's own deepest
# `z_faces` level, including a forced-water band along the open edge.
function implicit_free_surface_grid(arch, grid_config::EvenGrid, static_grid, relaxation_edge, boundary_band_cells)
    underlying_grid = LatitudeLongitudeGrid(
        arch;
        size      = grid_config.size,
        halo      = grid_config.halo,
        longitude = grid_config.longitude,
        latitude  = grid_config.latitude,
        z         = MutableVerticalDiscretization(grid_config.z_faces),
    )
    Nx, Ny, _ = size(static_grid)
    # `static_grid.immersed_boundary.bottom_height` is the raw halo-including `OffsetArray`
    # Oceananigans stores, not a `Field` — `interior` doesn't apply to it directly. Indexing with a
    # range materializes a fresh (non-offset) array, so mutating it below does not touch
    # `static_grid`.
    bottom = static_grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]
    depth = abs(first(grid_config.z_faces))
    flatten_bathymetry!(bottom, Val(relaxation_edge), boundary_band_cells, depth)
    bottom_height = Field{Center,Center,Nothing}(underlying_grid)
    set!(bottom_height, bottom)
    fill_halo_regions!(bottom_height)
    return ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = false)
end

# Same nine fields as `CoupledHydrostaticSimulation`, plus `grid_config`, `relaxation_edge` and
# `boundary_band_cells`: `coupled_simulation` receives only the already-built static `grid`, and
# needs these to rebuild the mutable-vertical-coordinate one described above, including the
# flattened bathymetry's forced-water band at the open edge.
struct OpenBoundaryHydrostaticSimulation{B,C,TA,MA,TR,CO,SI,BG,FS} <: AbstractCoupledSimulationConfig
    grid_config::EvenGrid
    relaxation_edge::Symbol
    boundary_band_cells::Int
    buoyancy::B
    closure::C
    tracer_advection::TA
    momentum_advection::MA
    tracers::TR
    coriolis::CO
    sea_ice::SI
    biogeochemistry::BG
    free_surface::FS
end

OpenBoundaryHydrostaticSimulation(;
    grid_config,
    relaxation_edge,
    boundary_band_cells,
    buoyancy,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    coriolis,
    sea_ice,
    biogeochemistry,
    free_surface,
) = OpenBoundaryHydrostaticSimulation(
    grid_config,
    relaxation_edge,
    boundary_band_cells,
    buoyancy,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    coriolis,
    sea_ice,
    biogeochemistry,
    free_surface,
)

FjordSim.model_tracers(model::OpenBoundaryHydrostaticSimulation) = model.tracers

# Extends FjordSim's exported `coupled_simulation` hook — the one place `AbstractCoupledSimulationConfig`
# is handed the grid, which is why the vertical-coordinate swap happens here rather than anywhere
# upstream of it. Otherwise identical to `CoupledHydrostaticSimulation`'s method in
# `src/Simulations.jl`.
function FjordSim.coupled_simulation(
    model::OpenBoundaryHydrostaticSimulation,
    grid;
    forcing,
    boundary_conditions,
    initial_conditions,
    atmosphere,
    radiation,
    stop_time,
    initial_time_step,
)
    arch = architecture(grid)
    zstar_grid =
        implicit_free_surface_grid(arch, model.grid_config, grid, model.relaxation_edge, model.boundary_band_cells)

    @info "Compiling HydrostaticFreeSurfaceModel"
    ocean_model = HydrostaticFreeSurfaceModel(
        zstar_grid;
        buoyancy = model.buoyancy,
        closure = model.closure,
        tracer_advection = model.tracer_advection,
        momentum_advection = model.momentum_advection,
        tracers = model.tracers,
        free_surface = FjordSim.free_surface(model.free_surface, zstar_grid),
        vertical_coordinate = ZStarCoordinate(),
        coriolis = model.coriolis,
        forcing = forcing,
        boundary_conditions = boundary_conditions,
        biogeochemistry = model.biogeochemistry,
    )
    @info "Compiled HydrostaticFreeSurfaceModel"

    set!(ocean_model; initial_conditions...)

    Δt = initial_time_step
    ocean_sim = Simulation(ocean_model; Δt, stop_time)
    coupled_model = OceanSeaIceModel(ocean_sim, model.sea_ice; atmosphere, radiation)
    @info "Initialized coupled model"

    return Simulation(coupled_model; Δt, stop_time)
end

# --- The config itself ------------------------------------------------------------------------

data_root = joinpath(homedir(), "FjordSim_data", "oslofjorden")
FT = Oceananigans.defaults.FloatType

# Named rather than inlined into `FjordConfig` below: `OpenBoundaryHydrostaticSimulation` needs the
# same geometry to rebuild its mutable-vertical-coordinate grid in `coupled_simulation`.
grid_config = EvenGrid(
    size      = (240, 520, 18),
    halo      = (7, 7, 7),
    longitude = (10.2, 11.02),
    latitude  = (59.0, 59.93),
    z_faces   = [
        -450.0, -400.0, -350.0, -300.0, -250.0, -200.0, -150.0, -100.0,
        -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0,
    ],
)

# Named rather than inlined: `RadiatingLateralBoundary`, `OpenBoundaryHydrostaticSimulation` and
# `implicit_free_surface_grid`'s bathymetry flatten all need to agree on which edge is open and how
# wide its forced-water band is.
relaxation_edge = :south
boundary_band_cells = 10
relaxation_timescale = 86400.0

# Both the run's initial state and the open boundary's quiescent exterior state, so the two agree
# — see `RadiatingLateralBoundary` below.
initial_conditions = (T = 5.0, S = 33.0)

FjordConfig(
    # Same grid and bathymetry as `oslofjorden()` — this variant reads the same prepared
    # bathymetry rather than downloading its own; only the land/sea mask is kept, since
    # `implicit_free_surface_grid` flattens every wet depth. Forcing is disabled entirely so the
    # open boundary relaxes towards constants instead of real NorKyst data — isolating the
    # open-boundary/`ImplicitFreeSurface`/`ZStarCoordinate` mechanism from both bathymetry- and
    # forcing-driven confounds. See the header comment for the full rationale.
    grid_config = grid_config,
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry_boundary.nc",
        plot_file             = "bathymetry_boundary.png",
        raw_resolution_factor = 2,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        minimum_depth         = 2.0,
        spike_ratio           = 0.5,
        max_slope_factor      = 0.5,
        geonorge_cache        = true,
        regrid_cache          = false,
    ),
    forcing_config = nothing,
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        years            = [2020],
    ),
    simulation_config = SimulationConfig(
        # Separate from `oslofjorden()`'s results root: checkpoints are shared per `results_root`
        # (one resumable run per directory), so this variant needs its own to avoid clobbering
        # `oslofjorden()`'s checkpoints and snapshots.
        results_root       = joinpath(homedir(), "FjordSim_results", "oslofjorden_implicit_open"),
        architecture       = :auto,
        model              = OpenBoundaryHydrostaticSimulation(
            grid_config         = grid_config,
            relaxation_edge     = relaxation_edge,
            boundary_band_cells = boundary_band_cells,
            buoyancy            = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
            closure             = (
                CATKEVerticalDiffusivity(minimum_tke = 7e-6),
                HorizontalScalarBiharmonicDiffusivity(ν = 1e5, κ = 1e4),
            ),
            tracer_advection    = (T = WENO(), S = WENO()),
            momentum_advection  = WENOVectorInvariant(FT),
            tracers             = (:T, :S),
            coriolis            = HydrostaticSphericalCoriolis(FT),
            sea_ice             = FreezingLimitedOceanTemperature(),
            biogeochemistry     = nothing,
            # Implicit, on a mutable-vertical-coordinate grid: no barotropic substepping, solved
            # through a preconditioned conjugate gradient, with η actually changing column depth
            # (`ZStarCoordinate`, wired up in `coupled_simulation` above) rather than only being
            # diagnosed — the ingredient the open boundary below needs to stay stable. See the
            # header comment for why.
            free_surface        = ImplicitFreeSurfaceConfig(
                solver_method              = :PreconditionedConjugateGradient,
                gravitational_acceleration = 9.80665,
                reltol                     = 1e-7,
                abstol                     = 1e-7,
                maxiter                    = 100,
            ),
        ),
        # `RadiatingLateralBoundary` in place of `OpenLateralBoundary`: same surface fluxes and
        # bottom drag, but the open southern edge now radiates velocity instead of holding it at a
        # closed wall — towards a quiescent exterior (`exterior_velocity = 0`) and the same T/S the
        # run starts at, rather than real forcing data.
        boundary_conditions = (
            TopBottomFluxes(bottom_drag_coefficient = 0.003),
            RadiatingLateralBoundary(
                edge                 = relaxation_edge,
                relaxation_timescale = relaxation_timescale,
                exterior_velocity    = 0.0,
                exterior_tracers     = initial_conditions,
            ),
        ),
        writers = (
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v),
                interval           = 1hour,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 30days, overwrite_existing = true, cleanup = true),
        ),
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.1,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        initial_conditions = initial_conditions,
        start_date         = DateTime(2020, 1, 1),
        stop_time          = 366days,
        loops              = 1,
        progress_interval  = 1hour,
        pickup             = false,
    ),
)
