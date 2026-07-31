using Dates: DateTime
using Oceananigans
using Oceananigans.Units
using NumericalEarth
using SeawaterPolynomials.TEOS10
using FjordSim
using FjordSim.Datasets
using CUDA  # it should be here to make GPU() not throw an error

const FT = Oceananigans.defaults.FloatType

architecture = GPU()
grid = ImmersedBoundaryGrid(
    joinpath(homedir(), "FjordSim_data", "oslofjord", "bathymetry_105to232.nc"),
    architecture,
    (7, 7, 7),
)
buoyancy = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT))
closure = (
    CATKEVerticalDiffusivity(minimum_tke = 7e-6),
    Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν = 15, κ = 10),
)
tracer_advection = (T = WENO(), S = WENO())
momentum_advection = WENOVectorInvariant(FT)
tracers = (:T, :S)
# dataset = ResultsDataset(
#     "snapshots_ocean_2.nc",
#     joinpath(homedir(), "FjordSim_results", "oslofjord");
#     start_date_time = DateTime(2025, 1, 1),
# )
# initial_conditions = (
#     T = Metadatum(:temperature; dataset, date = last_date(dataset, :temperature)),
#     S = Metadatum(:salinity; dataset, date = last_date(dataset, :salinity)),
#     u = Metadatum(:u_velocity; dataset, date = last_date(dataset, :u_velocity)),
#     v = Metadatum(:v_velocity; dataset, date = last_date(dataset, :v_velocity)),
# )
initial_conditions = (T = 5.0, S = 33.0)
free_surface = SplitExplicitFreeSurface(grid, cfl = 0.7)
coriolis = HydrostaticSphericalCoriolis(FT)
forcing = forcing_from_file(;
    grid = grid,
    filepath = joinpath(homedir(), "FjordSim_data", "oslofjord", "forcing_105to232.nc"),
    tracers = tracers,
)
tbbc = top_bottom_boundary_conditions(; grid = grid, bottom_drag_coefficient = 0.003)
sobc = (v = (south = OpenBoundaryCondition(nothing),),)
boundary_conditions = map(x -> FieldBoundaryConditions(; x...), recursive_merge(tbbc, sobc))
# atmosphere = JRA55PrescribedAtmosphere(
#     architecture,
#     FT;
#     latitude = (58.98, 59.94),
#     longitude = (10.18, 11.03),
#     dir = joinpath(homedir(), "FjordSim_data", "JRA55"),
# )
atmosphere = NORA3PrescribedAtmosphere(architecture)
downwelling_radiation = NORA3PrescribedRadiation(architecture)
sea_ice = FreezingLimitedOceanTemperature()
biogeochemistry = nothing
results_dir = joinpath(homedir(), "FjordSim_results", "oslofjord")
stop_time = 365days

simulation = coupled_hydrostatic_simulation(
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
    results_dir,
    stop_time,
)

simulation.callbacks[:progress] = Callback(progress, TimeInterval(1hour))
ocean_sim = simulation.model.ocean
ocean_model = ocean_sim.model

prefix = joinpath(results_dir, "snapshots_ocean")
ocean_sim.output_writers[:ocean] = NetCDFWriter(
    ocean_model,
    (T = ocean_model.tracers.T, S = ocean_model.tracers.S, u = ocean_model.velocities.u, v = ocean_model.velocities.v);
    filename = "$prefix",
    schedule = TimeInterval(1hour),
    overwrite_existing = true,
)

conjure_time_step_wizard!(
    simulation;
    cfl = 0.1,
    max_Δt = 3minutes,
    max_change = 1.01,
    cell_advection_timescale = cell_advection_timescale_coupled_model,
)
run!(simulation)
