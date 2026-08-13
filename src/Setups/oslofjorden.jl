"""
    oslofjorden()

The Oslofjord setup: Geonorge Sjøkart Dybdedata bathymetry, NorKyst-800m forcing with OF800
rivers, and NORA3 atmosphere, on a 240 x 520 x 18 grid covering 10.2-11.02°E, 59.0-59.93°N.

The river config gets the same `data_root` as the rest of the setup, so the river data downloads
alongside it rather than being shared from elsewhere.

It also names a `simulation_config`, so `run_simulation` works once the preparation steps have
run. Nothing about the run has a default, so every knob is stated here and nowhere else — split
across four nested configs by what dispatches on it: `CoupledHydrostaticSimulation` is what the
model is, the `boundary_conditions` tuple is what bounds it, the `writers` tuple is what it writes,
and `AdaptiveTimeStep` is how its clock advances. Results go to `~/FjordSim_results/oslofjorden/`,
separate from the input data root, and each writer's file carries the launch tag so runs do not
overwrite each other.

`start_date` and `stop_time` also decide what the prepare steps write: they pad the forcing and
atmosphere time axes to span exactly this window, so changing either means re-running
`prepare_forcing`, `add_rivers` and `prepare_atmosphere`. The window here is the whole of 2020, which
needs that padding at both ends — files prepared against the old 12:00-anchored window will be
rejected by `validate_time_coverage` until those three steps have been re-run.
"""
function oslofjorden()
    data_root = joinpath(homedir(), "FjordSim_data", "oslofjorden")
    FT = Oceananigans.defaults.FloatType

    # FjordConfig overload entry points like run_simulation(), download_forcing(), etc.
    return FjordConfig(
        # Overloads LatitudeLongitudeGrid(architecture, config) — the one grid hook.
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
            minimum_depth         = 2.0,
            spike_ratio           = 0.5,
            max_slope_factor      = 0.5,
            geonorge_cache        = true,
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
            open_edge        = :south,
            architecture     = :auto,
            parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
            years            = [2020],
            rivers           = OF800RiversConfig(data_root = data_root),
            # The exterior state along the open southern edge, from the *hourly* NorKyst
            # collection: a Flather boundary compares the model's own η against the exterior one,
            # and the daily means above have the tide averaged out of them. Only the boundary row
            # is prepared at that cadence, so the file is a few hundred MB.
            boundaries       = NorKystBoundariesConfig(
                data_root        = data_root,
                output_directory = "norkyst_hourly",
                output_file      = "boundaries.nc",
                plot_file        = "boundaries.png",
                margin           = 0.05,
                architecture     = :auto,
                parameters       = [
                    "temperature", "salinity", "u_eastward", "v_northward", "zeta", "ubar", "vbar",
                ],
                years            = [2020],
            ),
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
            results_root       = joinpath(homedir(), "FjordSim_results", "oslofjorden"),
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
            # separate pieces, so either can be swapped alone, plus the genuinely open southern
            # edge — which edge comes from the forcing config, and the exterior state from the
            # boundary dataset hanging off it. Dropping it would close the domain; argument order is
            # merge precedence.
            #
            # The two timescales are Marchesiello et al. (2001): nudge hard towards the data on
            # inflow, let radiation do the work on outflow.
            boundary_conditions = MergedBoundaryConditions(
                AirSeaFluxes(),
                QuadraticBottomDrag(coefficient = 0.003),
                OpenLateralBoundaryFromForcing(
                    inflow_timescale  = 1day,
                    outflow_timescale = 360days,
                ),
            ),
            # Both overload attach_writer! — what the run writes. `variables` may name anything
            # `Oceananigans.fields` exposes on the ocean model — add `:w`, `:η` or a
            # biogeochemical tracer by naming it here. Dropping the `CheckpointWriter` turns
            # checkpointing off entirely.
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
            # FromForcing(): the NorKyst state at `start_date` rather than a uniform water
            # column: every tracer the model names plus u and v, whichever of them the forcing
            # file carries. A literal NamedTuple (`(T = 5.0, S = 33.0)`) still works, and
            # `FromResults("snapshots_ocean_<tag>.nc")` continues from a previous run instead.
            initial_conditions = (T = 5.0, S = 33.0),
            # The whole calendar year 2020, which is a leap year — so 366 days from midnight on
            # 1 January lands exactly on midnight a year later. Neither prepared file has a record
            # at either end natively (NorKyst's are daily at 12:00, NORA3's hourly from 00:00 to
            # 23:00), so both prepare steps pad their axes to reach them: 12 hours at each end of
            # the forcing and one hour at the end of the atmosphere, each within the one-record-
            # spacing bound. That is what lets the window be a round year instead of being pinned
            # to whichever instant the forcing happened to start at.
            start_date         = DateTime(2020, 1, 1),
            stop_time          = 366days,
            # One pass. Raise it to spin the deep basins up on the same forcing year, carrying the
            # state over; each repetition writes its own `_loopNN` file.
            loops              = 1,
            pickup             = false,
        ),
    )
end
