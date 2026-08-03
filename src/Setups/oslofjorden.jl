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
root, named after `start_date` so runs do not overwrite each other.

`start_date` and `stop_time` also decide what the prepare steps write: they pad the forcing and
atmosphere time axes to span exactly this window, so changing either means re-running
`prepare_forcing`, `add_rivers` and `prepare_atmosphere`. The window here is the whole of 2020, which
needs that padding at both ends — files prepared against the old 12:00-anchored window will be
rejected by `validate_time_coverage` until those three steps have been re-run.
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
                HorizontalScalarBiharmonicDiffusivity(ν = 1e5, κ = 1e4),
            ),
            tracer_advection        = (T = WENO(), S = WENO()),
            momentum_advection      = WENOVectorInvariant(FT),
            tracers                 = (:T, :S),
            # FromForcing(): the NorKyst state at `start_date` rather than a uniform water 
            # column: every tracer named above plus u and v, whichever of them the forcing 
            # file carries. A literal NamedTuple (`(T = 5.0, S = 33.0)`) still works, and
            # `FromResults("snapshots_ocean_<tag>.nc")` continues from a previous run instead.
            initial_conditions      = (T = 5.0, S = 33.0),
            coriolis                = HydrostaticSphericalCoriolis(FT),
            sea_ice                 = FreezingLimitedOceanTemperature(),
            biogeochemistry         = nothing,
            free_surface_cfl        = 0.7,
            bottom_drag_coefficient = 0.003,
            # The whole calendar year 2020, which is a leap year — so 366 days from midnight on
            # 1 January lands exactly on midnight a year later. Neither prepared file has a record
            # at either end natively (NorKyst's are daily at 12:00, NORA3's hourly from 00:00 to
            # 23:00), so both prepare steps pad their axes to reach them: 12 hours at each end of
            # the forcing and one hour at the end of the atmosphere, each within the one-record-
            # spacing bound. That is what lets the window be a round year instead of being pinned
            # to whichever instant the forcing happened to start at.
            start_date              = DateTime(2020, 1, 1),
            stop_time               = 366days,
            # One pass. Raise it to spin the deep basins up on the same forcing year, carrying the
            # state over; each repetition writes its own `_loopNN` file.
            loops                   = 1,
            output_interval         = 1hour,
            progress_interval       = 1hour,
            overwrite_existing      = true,
            checkpoint_interval     = 30days,
            pickup                  = false,
            time_step_cfl           = 0.1,
            max_time_step           = 3minutes,
            max_time_step_change    = 1.01,
        ),
    )
end
