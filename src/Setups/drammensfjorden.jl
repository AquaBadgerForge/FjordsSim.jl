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

It also names an `atmosphere_config` and a `simulation_config` — the same model, boundary conditions
and time-step wizard settings as `oslofjorden()`, so it exercises the same physics. The three
differences are all scope, not kind: it runs a 30-day window (`start_date`/`stop_time`) rather than a
full year, which is what makes it a fast first run to try the whole pipeline on; its
`initial_conditions` are `FromForcing()` rather than a literal `NamedTuple`, reading the state at
`start_date` from that rivers-augmented copy once `add_rivers` has run; and it names no
`CheckpointWriter`, so it writes snapshots and nothing else. Results go to
`~/FjordSim_results/drammensfjorden/`, separate from the input data root.
"""
function drammensfjorden()
    data_root = joinpath(homedir(), "FjordSim_data", "drammensfjorden")
    FT = Oceananigans.defaults.FloatType

    return FjordConfig(
        # Overloads LatitudeLongitudeGrid(architecture, config) — the one grid hook.
        grid_config = EvenGrid(
            size      = (150, 200, 11),
            halo      = (7, 7, 7),
            longitude = (10.20, 10.45),
            latitude  = (59.58, 59.75),
            z_faces   = [-100.0, -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0],
        ),
        # Overloads bathymetry_dataset (required), regrid_options and smoothing_options (both
        # optional) — the hooks prepare_bathymetry dispatches on.
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
        # Overloads forcing_time_steps, forcing_source_grid, forcing_variable_names and
        # download_forcing (required) — the hooks prepare_forcing and download_forcing dispatch
        # on. simulation_forcing is left at its default, since this is the FjordSim NetCDF
        # forcing contract build_simulation already reads.
        forcing_config = NorKystConfig(
            data_root        = data_root,
            output_directory = "norkyst",
            output_file      = "forcing.nc",
            plot_file        = "forcing.png",
            architecture     = :auto,
            parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
            years            = [2020],
            rivers           = OF800RiversConfig(data_root = data_root),
        ),
        # The hourly exterior state along the open southern edge, as in `oslofjorden()`, and the one
        # place that edge is stated. Its own `data_root`, so this setup downloads its own band rather
        # than sharing Oslofjord's: the band is derived from *this* grid's southern edge, which is
        # 20 km further north.
        boundary_config = NorKystBoundariesConfig(
            data_root        = data_root,
            output_directory = "norkyst_hourly",
            output_file      = "boundaries.nc",
            plot_file        = "boundaries.png",
            open_edges       = :south,
            margin           = 0.05,
            architecture     = :auto,
            parameters       = [
                "temperature", "salinity", "u_eastward", "v_northward", "zeta", "ubar", "vbar",
            ],
            years            = [2020],
        ),
        # Overloads atmosphere_time_steps, atmosphere_source_grid, atmosphere_variable_names,
        # download_atmosphere, prescribed_atmosphere and prescribed_radiation — the hooks
        # prepare_atmosphere, download_atmosphere and build_simulation dispatch on.
        atmosphere_config = NORA3Config(
            data_root        = data_root,
            output_directory = "nora3",
            output_file      = "atmosphere.nc",
            plot_file        = "atmosphere.png",
            resolution       = 0.02,
            padding          = 0.1,
            years            = [2020],
        ),
        # SimulationConfig itself has no hooks — everything below dispatches through one of its
        # four nested configs instead.
        simulation_config = SimulationConfig(
            results_root       = joinpath(homedir(), "FjordSim_results", "drammensfjorden"),
            architecture       = :auto,
            # Overloads coupled_simulation and model_tracers — the model hooks build_simulation
            # dispatches on.
            model              = CoupledHydrostaticSimulation(
                buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
                closure            = (
                    CATKEVerticalDiffusivity(minimum_tke = 7e-6),
                    HorizontalScalarBiharmonicDiffusivity(ν = 1e5, κ = 1e4),
                ),
                tracer_advection   = (T = WENO(), S = WENO()),
                momentum_advection = WENOVectorInvariant(FT),
                tracers            = (:T, :S),
                coriolis           = HydrostaticSphericalCoriolis(FT),
                sea_ice            = FreezingLimitedOceanTemperature(),
                biogeochemistry    = nothing,
                # Overloads free_surface(config, grid) — its own hook, called from inside
                # coupled_simulation once the grid exists.
                free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
                # Anything else the four constructors coupled_simulation calls accept, one slot
                # each: :ocean_model, :ocean_simulation, :coupled_model, :coupled_simulation.
                # Nothing extra here, but stated rather than defaulted like every other field.
                extra_kwargs       = (;),
            ),
            # MergedBoundaryConditions overloads field_boundary_conditions; each piece inside it
            # overloads boundary_condition_sides. Air-sea fluxes and quadratic bottom drag are
            # separate pieces, so either can be swapped alone, plus the genuinely open southern edge,
            # whose edge and exterior state both come from `boundary_config` above. The two timescales
            # are Marchesiello et al. (2001).
            boundary_conditions = MergedBoundaryConditions(
                AirSeaFluxes(),
                QuadraticBottomDrag(coefficient = 0.003),
                OpenLateralBoundaryFromData(
                    inflow_timescale  = 1day,
                    outflow_timescale = 360days,
                ),
            ),
            # Overloads attach_writer!. Snapshots only: a 30-day trial run has nothing worth
            # resuming, and naming no `CheckpointWriter` is how checkpointing is turned off. The
            # trailing comma is what makes this a one-element tuple rather than a bare writer.
            writers = (
                SnapshotWriter(
                    name               = :ocean,
                    output_file        = "snapshots_ocean.nc",
                    variables          = (:T, :S, :u, :v),
                    interval           = 1hour,
                    overwrite_existing = true,
                ),
            ),
            # Overloads attach_callback! — what the run reports while it runs. `report` is the
            # function itself, so a model whose tracers omit :T (which `progress` reads) names its
            # own here instead. An empty tuple runs silently.
            callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
            # Overloads attach_time_stepping! and initial_time_step.
            time_stepping = AdaptiveTimeStep(
                initial_time_step    = 1second,
                cfl                  = 0.1,
                max_time_step        = 3minutes,
                max_time_step_change = 1.01,
            ),
            # The NorKyst state at `start_date`: every tracer the model names plus u and v,
            # whichever of them the prepared forcing file carries.
            initial_conditions = FromForcing(),
            # A 30-day window rather than oslofjorden's full year, so this setup's `run_simulation`
            # finishes quickly as a first end-to-end run. Both prepare steps still pad their axes
            # to reach it.
            start_date         = DateTime(2020, 1, 1),
            stop_time          = 30days,
            loops              = 1,
            pickup             = false,
        ),
    )
end
