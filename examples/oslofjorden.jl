# Oslofjord variant testing whether a radiating open lateral boundary is stable with
# `ImplicitFreeSurface`. Everything else is stripped to isolate that: no ocean forcing, and the
# bathymetry flattened to one constant depth at sea, keeping only the land/sea mask.
#
# `ZStarCoordinate` is load-bearing, not decoration: on a static grid η is solved but never fed back
# into cell thickness, so a persistent inflow at the open boundary cannot reconcile into the water
# column and the PCG solver eventually throws a `DomainError`. Oceananigans states the coordinate in
# two places — `MutableVerticalDiscretization` on the grid's `z` and `vertical_coordinate` on the
# model — so this file does too.
#
# The atmosphere is shared with `oslofjorden()`; the bathymetry is not. Prepare both, then run:
#
#   julia --project -m FjordSim download_atmosphere --config oslofjorden
#   julia --project -m FjordSim prepare_atmosphere  --config oslofjorden
#   julia --project -m FjordSim prepare_bathymetry  --config examples/oslofjorden.jl
#   julia --project -m FjordSim run_simulation      --config examples/oslofjorden.jl
#

using FjordSim
using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units
using Dates: DateTime
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature

# --- ImplicitFreeSurfaceConfig: a new AbstractFreeSurfaceConfig ------------------------------

# `1e-7` matches `PCGImplicitFreeSurfaceSolver`'s own computed default. Tightening it to `1e-10` (as
# `validation/open_boundaries/barotropic_soliton.jl` does) throws a `DomainError`: the solver
# converges in ~9 iterations either way, but at that tolerance an inner product rounds just below
# zero on this grid and `sqrt` of it fails.
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

FjordSim.free_surface(config::ImplicitFreeSurfaceConfig, grid) = ImplicitFreeSurface(
    solver_method              = config.solver_method,
    gravitational_acceleration = config.gravitational_acceleration,
    reltol                     = config.reltol,
    abstol                     = config.abstol,
    maxiter                    = config.maxiter,
)

# --- RadiatingLateralBoundary: a new AbstractBoundaryConditionConfig -------------------------

# Replaces FjordSim's data-driven `OpenLateralBoundaryFromData` with a radiating one relaxing towards
# a quiescent exterior. Carries its own `edge` and `relaxation_timescale` rather than reading them off
# a boundary data config, since this setup names neither `forcing_config` nor `boundary_config` — its
# exterior state is a pair of constants, not a prepared file.
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

# The fourth argument is the setup's boundary data config, which this piece has no use for: it
# carries its own edge and its own constant exterior state, so it implements the five-argument form
# and never sees the prepared boundary series either.
FjordSim.boundary_condition_sides(config::RadiatingLateralBoundary, grid, forcing, boundary_config, tracers) =
    recursive_merge(
        open_velocity_boundary_conditions(Val(config.edge), config),
        open_tracer_boundary_conditions(Val(config.edge), config, tracers),
    )

# --- FlatMutableGrid: a new AbstractGridConfig -----------------------------------------------

# A uniform column depth at sea gives every wet cell the same inertia and the same number of active
# vertical layers, so an instability cannot be attributed to bathymetry. The band along `edge` is
# forced to water even where the mask says land, so the open edge always has a full water column.
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

struct FlatMutableGrid <: AbstractGridConfig
    geometry::EvenGrid
    open_edge::Symbol
    boundary_band_cells::Int
end

FlatMutableGrid(; geometry, open_edge, boundary_band_cells) =
    FlatMutableGrid(geometry, open_edge, Int(boundary_band_cells))

# The prepare pipelines regrid onto the static geometry; only the simulation grid needs a mutable
# vertical coordinate.
FjordSim.domain_grid(config::FlatMutableGrid, architecture) = domain_grid(config.geometry, architecture)

function FjordSim.simulation_grid(config::FlatMutableGrid, bathymetry_file, architecture)
    geometry = config.geometry
    underlying_grid = LatitudeLongitudeGrid(
        architecture;
        size      = geometry.size,
        halo      = geometry.halo,
        longitude = geometry.longitude,
        latitude  = geometry.latitude,
        z         = MutableVerticalDiscretization(geometry.z_faces),
    )

    static_grid = ImmersedBoundaryGrid(bathymetry_file, architecture, geometry.halo)
    Nx, Ny, _ = size(static_grid)
    # `bottom_height` is the raw halo-including `OffsetArray`, not a `Field`, so `interior` does not
    # apply. Range indexing materializes a fresh array, so the flatten below leaves `static_grid` be.
    bottom = static_grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]
    depth = abs(first(geometry.z_faces))
    flatten_bathymetry!(bottom, Val(config.open_edge), config.boundary_band_cells, depth)

    bottom_height = Field{Center,Center,Nothing}(underlying_grid)
    set!(bottom_height, bottom)
    fill_halo_regions!(bottom_height)
    return ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = false)
end

# --- The config itself ------------------------------------------------------------------------

data_root = joinpath(homedir(), "FjordSim_data", "oslofjorden")
FT = Oceananigans.defaults.FloatType

geometry = EvenGrid(
    size      = (240, 520, 18),
    halo      = (7, 7, 7),
    longitude = (10.2, 11.02),
    latitude  = (59.0, 59.93),
    z_faces   = [
        -450.0, -400.0, -350.0, -300.0, -250.0, -200.0, -150.0, -100.0,
        -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0,
    ],
)

# Shared so the boundary condition and the bathymetry flatten agree on which edge is open.
relaxation_edge = :south
boundary_band_cells = 10
relaxation_timescale = 86400.0

# Both the initial state and the open boundary's exterior state, so the two agree.
initial_conditions = (T = 5.0, S = 33.0)

FjordConfig(
    grid_config = FlatMutableGrid(
        geometry            = geometry,
        open_edge           = relaxation_edge,
        boundary_band_cells = boundary_band_cells,
    ),
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
        max_island_cells      = 6,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        minimum_cell_fraction = 0.2,
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
        # Its own root: checkpoints are shared per `results_root`, so sharing `oslofjorden()`'s
        # would clobber that run's.
        results_root       = joinpath(homedir(), "FjordSim_results", "oslofjorden_implicit_open"),
        architecture       = :auto,
        model              = CoupledHydrostaticSimulation(
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
            free_surface        = ImplicitFreeSurfaceConfig(
                solver_method              = :PreconditionedConjugateGradient,
                gravitational_acceleration = 9.80665,
                reltol                     = 1e-7,
                abstol                     = 1e-7,
                maxiter                    = 100,
            ),
            # The model half of the vertical-coordinate swap; the grid half is `FlatMutableGrid`'s
            # `MutableVerticalDiscretization`. The inner trailing comma is required — without it
            # this is a parenthesized assignment rather than a `NamedTuple`.
            extra_kwargs        = (ocean_model = (vertical_coordinate = ZStarCoordinate(),),),
        ),
        boundary_conditions = MergedBoundaryConditions(
            AirSeaFluxes(),
            QuadraticBottomDrag(coefficient = 0.003),
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
            CheckpointWriter(interval = 30days, cleanup = true),
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
        callbacks          = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
        pickup             = false,
    ),
)
