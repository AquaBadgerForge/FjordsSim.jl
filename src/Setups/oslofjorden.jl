"""
    oslofjorden()

The Oslofjord setup: Geonorge Sjøkart Dybdedata bathymetry, NorKyst-800m forcing with OF800
rivers, and NORA3 atmosphere, on a 240 x 520 x 18 grid covering 10.2-11.02°E, 59.0-59.93°N.

The river config gets the same `data_root` as the rest of the setup, so the river data downloads
alongside it rather than being shared from elsewhere.

It also names a `simulation_config`, so `run_simulation` works once the preparation steps have
run. `SimulationConfig` has no defaults, so every knob the simulation uses — buoyancy, closure,
advection, tracers, coriolis, sea ice, the run length and the wizard settings — is stated here
and nowhere else. Results go to `~/FjordSim_results/oslofjorden/`, separate from the input data
root.
"""
function oslofjorden()
    data_root = joinpath(homedir(), "FjordSim_data", "oslofjorden")
    FT = Oceananigans.defaults.FloatType

    return FjordConfig(
        grid_config = EvenGrid(
            size      = (240, 520, 18),
            halo      = (7, 7, 7),
            longitude = (10.2, 11.02),
            latitude  = (59.0, 59.93),
            z_faces   = [
                -450.0, -400.0, -350.0, -300.0, -250.0, -200.0, -150.0, -100.0,
                -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0,
            ],
        ),
        bathymetry_config = DybdedataConfig(
            data_root             = data_root,
            output_file           = "bathymetry.nc",
            plot_file             = "bathymetry.png",
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
            output_file          = "forcing.nc",
            plot_file            = "forcing.png",
            relaxation_edge      = :south,
            relaxation_cells     = 10,
            relaxation_timescale = 86400.0,
            architecture         = :auto,
            parameters           = ["temperature", "salinity", "u_eastward", "v_northward"],
            years                = [2020],
            rivers               = OF800RiversConfig(data_root = data_root),
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
            results_root            = joinpath(homedir(), "FjordSim_results", "oslofjorden"),
            output_file             = "snapshots_ocean.nc",
            architecture            = :auto,
            buoyancy                = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
            closure                 = (
                CATKEVerticalDiffusivity(minimum_tke = 7e-6),
                HorizontalScalarBiharmonicDiffusivity(ν = 15, κ = 10),
            ),
            tracer_advection        = (T = WENO(), S = WENO()),
            momentum_advection      = WENOVectorInvariant(FT),
            tracers                 = (:T, :S),
            initial_conditions      = (T = 5.0, S = 33.0),
            coriolis                = HydrostaticSphericalCoriolis(FT),
            sea_ice                 = FreezingLimitedOceanTemperature(),
            biogeochemistry         = nothing,
            free_surface_cfl        = 0.7,
            bottom_drag_coefficient = 0.003,
            stop_time               = 365days,
            output_interval         = 1hour,
            progress_interval       = 1hour,
            overwrite_existing      = true,
            time_step_cfl           = 0.1,
            max_time_step           = 3minutes,
            max_time_step_change    = 1.01,
        ),
    )
end
