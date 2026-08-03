"""
    drammensfjorden()

The Drammensfjord setup: Geonorge Sjøkart Dybdedata bathymetry, NorKyst-800m forcing with OF800
rivers, and NORA3 atmosphere on a 150 x 200 x 11 grid covering 10.20-10.45°E, 59.58-59.75°N.

The river config gets the same `data_root` as the rest of the setup, as `oslofjorden()`'s does, so
the river data downloads alongside it rather than being shared from elsewhere. Only Drammenselva,
outlet 19 of the dataset's 25, is inside this grid — the rest are Oslofjord outlets far to the east
and south, and Lierelva at 59.7502°N is 0.000625° north of the last latitude node. `river_cells`
drops each of the 24 with an `@info`, which is not an error: `add_rivers` fails only when *no*
outlet lands. It is still a required step, since `simulation_forcing_path` picks the
rivers-augmented copy once a setup names rivers.

It also names an `atmosphere_config` and a `simulation_config` — the same buoyancy, closure,
advection, tracers, coriolis, sea ice and time-step wizard settings as `oslofjorden()`, so it
exercises the same physics. The two setups differ in scope, not in kind: this one runs a 30-day
window (`start_date`/`stop_time`) rather than a full year, which is what makes it a fast first run
to try the whole pipeline on, and its `initial_conditions` are `FromForcing()` rather than a literal
`NamedTuple`, reading the state at `start_date` from that rivers-augmented copy once `add_rivers`
has run. Results go to `~/FjordSim_results/drammensfjorden/`, separate from the input data root.
"""
function drammensfjorden()
    data_root = joinpath(homedir(), "FjordSim_data", "drammensfjorden")
    FT = Oceananigans.defaults.FloatType

    return FjordConfig(
        grid_config = EvenGrid(
            size      = (150, 200, 11),
            halo      = (7, 7, 7),
            longitude = (10.20, 10.45),
            latitude  = (59.58, 59.75),
            z_faces   = [-100.0, -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0],
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
            geonorge_cache        = false,
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
            results_root            = joinpath(homedir(), "FjordSim_results", "drammensfjorden"),
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
            # The NorKyst state at `start_date`: every tracer named above plus u and v, whichever
            # of them the prepared forcing file carries.
            initial_conditions      = FromForcing(),
            coriolis                = HydrostaticSphericalCoriolis(FT),
            sea_ice                 = FreezingLimitedOceanTemperature(),
            biogeochemistry         = nothing,
            free_surface_cfl        = 0.7,
            bottom_drag_coefficient = 0.003,
            # A 30-day window rather than oslofjorden's full year, so this setup's `run_simulation`
            # finishes quickly as a first end-to-end run. Both prepare steps still pad their axes
            # to reach it.
            start_date              = DateTime(2020, 1, 1),
            stop_time               = 30days,
            loops                   = 1,
            output_interval         = 1hour,
            progress_interval       = 1hour,
            overwrite_existing      = true,
            checkpoint_interval     = 0.0,
            pickup                  = false,
            time_step_cfl           = 0.1,
            max_time_step           = 3minutes,
            max_time_step_change    = 1.01,
        ),
    )
end
