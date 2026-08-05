# Oslofjord variant using an implicit free surface and a true (radiating) open boundary condition
# on velocity, instead of the split-explicit free surface and closed-wall lateral velocity
# `oslofjorden()` (`src/Setups/oslofjorden.jl`) uses. See that setup for the rationale behind every
# field this one copies verbatim.
#
# Everything new is defined below, entirely outside the package, following FjordSim's own "Adding a
# new source" extension model — subtype a supertype, overload its hook:
#
#   - `ImplicitFreeSurfaceConfig` wraps Oceananigans' `ImplicitFreeSurface`, the free-surface family
#     open-boundary support (github.com/CliMA/Oceananigans.jl/issues/5229) is being built against
#     first, since its predictor/corrector structure mirrors `NonhydrostaticModel`'s pressure
#     correction.
#   - `RadiatingLateralBoundary` replaces `OpenLateralBoundary`'s closed wall on the velocity normal
#     to the forcing config's `relaxation_edge` with a radiating `NormalFlowBoundaryCondition`
#     (`PerturbationAdvection`), so flow can genuinely cross the domain's open edge rather than only
#     tracers relaxing across it.
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
# `flow_over_hill.jl`, which has west+east). This domain has exactly one open boundary, and
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
# static grid's bottom height carried over directly — see `implicit_free_surface_grid` below.
#
# This file is itself a config, not a runner — it reads the same prepared atmosphere
# `oslofjorden()`'s prepare steps write (same `data_root`), so those steps have to have run first:
#
#   julia --project -m FjordSim download_forcing    --config oslofjorden   # raw NorKyst files only
#   julia --project -m FjordSim download_atmosphere --config oslofjorden
#   julia --project -m FjordSim prepare_atmosphere  --config oslofjorden
#
# Bathymetry and forcing are *not* shared with `oslofjorden()`, on purpose: `bathymetry_config`
# points at `bathymetry_boundary_fixed.nc`, a one-off copy of `oslofjorden()`'s processed bathymetry
# with the south open-boundary band floored to 15 m (see `floor_boundary_depth!` below, and its
# comment for why — a thin water column there was found to blow up numerically). `forcing_config`
# is regridded fresh against that fixed bathymetry rather than reusing `oslofjorden()`'s
# `forcing.nc`, because `prepare_forcing`'s land/water mask is a property of the bathymetry it is
# built against — reusing the original regrid would leave `NaN` at the vertical levels the floor
# newly made active, right at the open boundary. So, once `bathymetry_boundary_fixed.nc` exists,
# this variant needs its own forcing prepared before it can run:
#
#   julia --project -m FjordSim prepare_forcing --config examples/oslofjord.jl
#   julia --project -m FjordSim add_rivers      --config examples/oslofjord.jl
#
# Then run this variant by pointing `--config` at this file directly:
#
#   julia --project -m FjordSim run_simulation --config examples/oslofjord.jl
#
# To step through the assembly instead of running it, build the simulation without starting it:
#
#   config = fjord_config("examples/oslofjord.jl")
#   simulation = build_simulation(config)
#   run!(simulation)

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

# Discrete-form boundary functions reading a velocity's own exterior `FieldTimeSeries` at the exact
# domain edge — the same shape as FjordSim's internal `x_boundary_tracer_value`/
# `y_boundary_tracer_value`, but for a Face-located velocity rather than a Center-located tracer.
@inline function x_boundary_velocity_value(j, k, grid, clock, model_fields, parameters)
    return @inbounds parameters.fts[parameters.boundary_index, j, k, Time(clock.time)]
end

@inline function y_boundary_velocity_value(i, k, grid, clock, model_fields, parameters)
    return @inbounds parameters.fts[i, parameters.boundary_index, k, Time(clock.time)]
end

# One velocity component's true, radiating open boundary condition at the exact domain edge: a
# `NormalFlowBoundaryCondition` with a `PerturbationAdvection` scheme relaxing towards `name`'s
# (`:u` or `:v`) own exterior series — in place of `OpenLateralBoundary`'s closed wall
# (`NormalFlowBoundaryCondition(nothing)`). No `target_transport`: that only matters for balancing
# *paired* open boundaries, and this domain has one.
function velocity_open_boundary_condition(boundary_function, boundary_index, forcing, name, config::AbstractForcingConfig)
    haskey(forcing, name) || error(
        "relaxation_edge names an open velocity boundary but the forcing does not prepare `$name`; " *
        "add it to `config.parameters`.",
    )
    scheme = PerturbationAdvection(;
        inflow_timescale = config.relaxation_timescale,
        outflow_timescale = config.relaxation_timescale,
    )
    parameters = (fts = getproperty(forcing, name).fts_value, boundary_index)
    return NormalFlowBoundaryCondition(boundary_function; discrete_form = true, parameters, scheme)
end

# The open-velocity counterpart of FjordSim's internal `open_boundary_conditions`. The boundary
# index is `1` for `:west`/`:south`; for `:east`/`:north` it is `grid.Nx + 1`/`grid.Ny + 1`, one past
# a tracer's `grid.Nx`/`grid.Ny`, because the normal velocity is Face-located and its own east/north
# domain boundary is genuinely one index further out.
open_velocity_boundary_conditions(::Val{:south}, grid, forcing, config) =
    (v = (south = velocity_open_boundary_condition(y_boundary_velocity_value, 1, forcing, :v, config),),)
open_velocity_boundary_conditions(::Val{:north}, grid, forcing, config) =
    (v = (north = velocity_open_boundary_condition(y_boundary_velocity_value, grid.Ny + 1, forcing, :v, config),),)
open_velocity_boundary_conditions(::Val{:west}, grid, forcing, config) =
    (u = (west = velocity_open_boundary_condition(x_boundary_velocity_value, 1, forcing, :u, config),),)
open_velocity_boundary_conditions(::Val{:east}, grid, forcing, config) =
    (u = (east = velocity_open_boundary_condition(x_boundary_velocity_value, grid.Nx + 1, forcing, :u, config),),)
open_velocity_boundary_conditions(::Val{edge}, grid, forcing, config) where {edge} =
    throw(ArgumentError("relaxation_edge must be one of (:south, :north, :west, :east), got :$edge"))

# The domain's open edge, with a genuinely open velocity: a radiating `NormalFlowBoundaryCondition`
# on the velocity normal to the forcing config's `relaxation_edge`, plus the same relaxing
# `ValueBoundaryCondition` on every simulated tracer `OpenLateralBoundary` builds. Carries no fields,
# for the same reason `OpenLateralBoundary` does: everything it needs is already in the forcing
# config and its prepared file.
struct RadiatingLateralBoundary <: AbstractBoundaryConditionConfig end

# Extends FjordSim's `boundary_conditions` hook, which is deliberately not re-exported (it would
# shadow `Oceananigans.Fields.boundary_conditions`) — reached the way `CLAUDE.md`'s Import
# Conventions section prescribes for extending a function you don't own: qualify the module, never
# `import Mod: foo`. `field_boundary_conditions`'s internal call to `boundary_conditions(config, ...)`
# dispatches to this method exactly as it does to the built-in ones, since a generic function has one
# global method table regardless of which module adds a method to it.
FjordSim.BoundaryConditions.boundary_conditions(::RadiatingLateralBoundary, grid, forcing, forcing_config, tracers) =
    recursive_merge(
        open_velocity_boundary_conditions(Val(forcing_config.relaxation_edge), grid, forcing, forcing_config),
        lateral_tracer_open_boundary_conditions(grid, forcing, forcing_config, tracers),
    )

# --- OpenBoundaryHydrostaticSimulation: a new AbstractCoupledSimulationConfig ------------------

# Floors the water-column depth in a `boundary_rows`-wide band along `edge`, leaving land cells
# (`bottom_height >= 0`) untouched — dispatch on `Val(edge)` mirrors `boundary_sponge_mask` above.
#
# Confirmed empirically, by manual time-stepping with per-cell diagnostics: a shallow (~5 m,
# barely more than the bathymetry's own `minimum_depth = 2.0` floor) column sitting directly on the
# open boundary row drove salinity to unphysical values (up to +204 psu, against a forcing exterior
# value of ~27–31 psu at the same cell and time — so the forcing data itself was clean) within the
# first minute of model time, well before the `ImplicitFreeSurface` solver's `DomainError` this was
# chased from. A thin water column has too little inertia to absorb the open boundary condition's
# momentum and tracer perturbations — the standard reason production regional models (ROMS, NEMO)
# keep open boundaries in reasonably deep water rather than running them up to the coastline. `15`
# meters keeps several full vertical layers wet even at the shallowest boundary cell (this grid's
# `z_faces` put four ~1–2 m layers in the top 5 m alone), comfortably deeper than the ~5 m columns
# found unstable, while staying modest next to the fjord's real depths (up to 450 m) so it does not
# distort the interior. The band matches `relaxation_cells`, the same width the boundary sponge and
# the interior relaxation zone already use.
floor_boundary_depth!(bottom, ::Val{:south}, boundary_rows, minimum_depth) =
    bottom[:, 1:boundary_rows] .=
        ifelse.(bottom[:, 1:boundary_rows] .< 0, min.(bottom[:, 1:boundary_rows], -minimum_depth), bottom[:, 1:boundary_rows])
floor_boundary_depth!(bottom, ::Val{:north}, boundary_rows, minimum_depth) =
    bottom[:, end-boundary_rows+1:end] .= ifelse.(
        bottom[:, end-boundary_rows+1:end] .< 0,
        min.(bottom[:, end-boundary_rows+1:end], -minimum_depth),
        bottom[:, end-boundary_rows+1:end],
    )
floor_boundary_depth!(bottom, ::Val{:west}, boundary_rows, minimum_depth) =
    bottom[1:boundary_rows, :] .=
        ifelse.(bottom[1:boundary_rows, :] .< 0, min.(bottom[1:boundary_rows, :], -minimum_depth), bottom[1:boundary_rows, :])
floor_boundary_depth!(bottom, ::Val{:east}, boundary_rows, minimum_depth) =
    bottom[end-boundary_rows+1:end, :] .= ifelse.(
        bottom[end-boundary_rows+1:end, :] .< 0,
        min.(bottom[end-boundary_rows+1:end, :], -minimum_depth),
        bottom[end-boundary_rows+1:end, :],
    )

# `implicit_free_surface_grid` rebuilds an equivalent grid with a mutable vertical coordinate from
# the ordinary static one `build_simulation` already built. Geometry (`size`, `halo`, `longitude`,
# `latitude`, `z_faces`) comes from the example's own `grid_config`, not by re-deriving it from
# `static_grid` — `src/Grids.jl`'s own doc comment flags exactly this face-vs-center count as an
# easy off-by-one. Bottom height is carried over from `static_grid` — so the two grids agree on
# where the seafloor is without a second bathymetry NetCDF read — except along the open boundary
# row, deepened by `floor_boundary_depth!` for the reason above. That makes this variant's grid
# deliberately diverge from the shared `bathymetry.nc` right at the boundary rather than editing
# that file, which `oslofjorden()`'s own setup also reads from the same `data_root`.
function implicit_free_surface_grid(arch, grid_config::EvenGrid, static_grid, relaxation_edge, relaxation_cells)
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
    floor_boundary_depth!(bottom, Val(relaxation_edge), relaxation_cells, 15.0)
    bottom_height = Field{Center,Center,Nothing}(underlying_grid)
    set!(bottom_height, bottom)
    fill_halo_regions!(bottom_height)
    return ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)
end

# Same nine fields as `CoupledHydrostaticSimulation`, plus `grid_config` and the forcing config's
# own `relaxation_edge`/`relaxation_cells`: `coupled_simulation` receives only the already-built
# static `grid`, and needs these to rebuild the mutable-vertical-coordinate one described above,
# including the boundary depth floor.
struct OpenBoundaryHydrostaticSimulation{B,C,TA,MA,TR,CO,SI,BG,FS} <: AbstractCoupledSimulationConfig
    grid_config::EvenGrid
    relaxation_edge::Symbol
    relaxation_cells::Int
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
    relaxation_cells,
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
    relaxation_cells,
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
        implicit_free_surface_grid(arch, model.grid_config, grid, model.relaxation_edge, model.relaxation_cells)

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

# --- Boundary momentum sponge ------------------------------------------------------------------

# A spatially-varying horizontal viscosity ramped up towards the open edge, the one ingredient
# `validation/open_boundaries/barotropic_soliton.jl` has beyond `ZStarCoordinate` that the simpler
# `flow_over_hill.jl` doesn't need: that validation case runs 150 days against an analytical
# solution and needs its sponge to keep the radiating boundary from accumulating spurious momentum
# over a long integration, which is exactly this variant's own failure mode — stable for ~600
# one-second steps (~10 minutes) on the `ZStarCoordinate` grid alone, but not yet proven for the
# full 366-day run.
#
# `PiecewiseLinearMask{D}` ramps linearly from 1 at the domain edge to 0 at `width` away — here
# `relaxation_cells` cells in from the forcing config's own `relaxation_edge`, so the sponge and
# the open boundary condition cover the same physical band. Dispatch on `Val(edge)` mirrors
# `open_velocity_boundary_conditions` above.
boundary_sponge_mask(::Val{:south}, grid_config::EvenGrid, relaxation_cells) = PiecewiseLinearMask{:y}(
    center = grid_config.latitude[1],
    width  = relaxation_cells * (grid_config.latitude[2] - grid_config.latitude[1]) / grid_config.size[2],
)
boundary_sponge_mask(::Val{:north}, grid_config::EvenGrid, relaxation_cells) = PiecewiseLinearMask{:y}(
    center = grid_config.latitude[2],
    width  = relaxation_cells * (grid_config.latitude[2] - grid_config.latitude[1]) / grid_config.size[2],
)
boundary_sponge_mask(::Val{:west}, grid_config::EvenGrid, relaxation_cells) = PiecewiseLinearMask{:x}(
    center = grid_config.longitude[1],
    width  = relaxation_cells * (grid_config.longitude[2] - grid_config.longitude[1]) / grid_config.size[1],
)
boundary_sponge_mask(::Val{:east}, grid_config::EvenGrid, relaxation_cells) = PiecewiseLinearMask{:x}(
    center = grid_config.longitude[2],
    width  = relaxation_cells * (grid_config.longitude[2] - grid_config.longitude[1]) / grid_config.size[1],
)

# --- The config itself ------------------------------------------------------------------------

data_root = joinpath(homedir(), "FjordSim_data", "oslofjorden")
FT = Oceananigans.defaults.FloatType

# Named rather than inlined into `FjordConfig` below: `OpenBoundaryHydrostaticSimulation` needs the
# same geometry to rebuild its mutable-vertical-coordinate grid in `coupled_simulation`, and the
# boundary sponge needs it to place its mask.
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

# Named rather than inlined into `NorKystConfig` below: the boundary sponge needs the same values
# to place its mask on the same edge and band the open boundary condition uses.
relaxation_edge = :south
relaxation_cells = 10

# A `let` block, not a bare top-level `sponge_ν(λ, φ, z, t) = ...` reading a global `sponge_mask`:
# this closure runs inside a GPU-compiled kernel (`HorizontalScalarDiffusivity`'s `ν`), and a
# non-`const` top-level global has no fixed type, so a kernel referencing one directly fails GPU
# compilation with `InvalidIRError` rather than merely running slowly — confirmed empirically here.
# The `let` captures `mask` and `ν_peak` as genuinely typed closure fields instead. Peak viscosity
# is bounded well under the diffusive stability limit (ν·Δt/Δx² ≲ 0.5) at the time-stepping
# config's own `max_time_step` below (3 minutes) and this grid's ~190 m cell spacing:
# 50 × 180 / 190² ≈ 0.25.
sponge_ν = let mask = boundary_sponge_mask(Val(relaxation_edge), grid_config, relaxation_cells),
               ν_peak = 50.0
    (λ, φ, z, t) -> ν_peak * mask(λ, φ, z)
end

FjordConfig(
    # Same grid, bathymetry, forcing and atmosphere as `oslofjorden()` — this variant reads the same
    # prepared files rather than downloading or regridding its own.
    grid_config = grid_config,
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        # `bathymetry_boundary_fixed.nc` is a one-off copy of `oslofjorden()`'s own `bathymetry.nc`
        # with the south open-boundary band floored to 15 m (see `floor_boundary_depth!` above);
        # not `prepare_bathymetry`'s regular output, and not shared with that setup, so this variant
        # is the only thing that reads or writes it.
        output_file           = "bathymetry_boundary_fixed.nc",
        plot_file             = "bathymetry_boundary_fixed.png",
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
    forcing_config = NorKystConfig(
        data_root            = data_root,
        output_directory     = "norkyst",
        # Regridded fresh against the fixed bathymetry above rather than sharing `oslofjorden()`'s
        # `forcing.nc`: `prepare_forcing`'s land/water mask is built from the bathymetry at
        # regrid time (`water_mask`, via `peripheral_node`), so a copy regridded against the
        # original shallow depths would still carry `NaN` at the vertical levels the floor newly
        # made active — exactly the mismatch this variant's boundary instability was traced to.
        output_file          = "forcing_boundary_fixed.nc",
        plot_file            = "forcing_boundary_fixed.png",
        relaxation_edge      = relaxation_edge,
        relaxation_cells     = relaxation_cells,
        relaxation_timescale = 86400.0,
        architecture         = :auto,
        parameters           = ["temperature", "salinity", "u_eastward", "v_northward"],
        years                = [2020],
        rivers               = OF800RiversConfig(
            data_root   = data_root,
            output_file = "forcing_rivers_boundary_fixed.nc",
        ),
    ),
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
            grid_config        = grid_config,
            relaxation_edge    = relaxation_edge,
            relaxation_cells   = relaxation_cells,
            buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
            closure            = (
                CATKEVerticalDiffusivity(minimum_tke = 7e-6),
                HorizontalScalarBiharmonicDiffusivity(ν = 1e5, κ = 1e4),
                # Boundary momentum sponge: extra harmonic viscosity ramped up towards the open
                # edge, on top of the constant biharmonic viscosity above. See its definition
                # under "Boundary momentum sponge" for why.
                HorizontalScalarDiffusivity(ν = sponge_ν),
            ),
            tracer_advection   = (T = WENO(), S = WENO()),
            momentum_advection = WENOVectorInvariant(FT),
            tracers            = (:T, :S),
            coriolis           = HydrostaticSphericalCoriolis(FT),
            sea_ice            = FreezingLimitedOceanTemperature(),
            biogeochemistry    = nothing,
            # Implicit, on a mutable-vertical-coordinate grid: no barotropic substepping, solved
            # through a preconditioned conjugate gradient, with η actually changing column depth
            # (`ZStarCoordinate`, wired up in `coupled_simulation` above) rather than only being
            # diagnosed — the ingredient the open boundary below needs to stay stable. See the
            # header comment for why.
            free_surface       = ImplicitFreeSurfaceConfig(
                solver_method              = :PreconditionedConjugateGradient,
                gravitational_acceleration = 9.80665,
                reltol                     = 1e-7,
                abstol                     = 1e-7,
                maxiter                    = 100,
            ),
        ),
        # `RadiatingLateralBoundary` in place of `OpenLateralBoundary`: same surface fluxes and
        # bottom drag, but the open southern edge now radiates velocity instead of holding it at a
        # closed wall.
        boundary_conditions = (
            TopBottomFluxes(bottom_drag_coefficient = 0.003),
            RadiatingLateralBoundary(),
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
        initial_conditions = (T = 5.0, S = 33.0),
        start_date         = DateTime(2020, 1, 1),
        stop_time          = 366days,
        loops              = 1,
        progress_interval  = 1hour,
        pickup             = false,
    ),
)
