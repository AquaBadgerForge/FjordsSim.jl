"""
    oslofjorden()

The Oslofjord setup: Geonorge Sjøkart Dybdedata bathymetry, NorKyst-800m forcing with NVE rivers,
and NORA3 atmosphere, on a 240 x 520 x 24 grid covering 10.2-11.02°E, 59.0-59.93°N.

The rivers are **discovered**, not stated: `minimum_discharge = 0.5` has `add_rivers` read NVE's
ELVIS river network and REGINE catchments for this domain and place every mouth at or above that
size — 21 of them, carrying 99.6 % of the domain's 1179 m³/s — so no river coordinate appears in
this file. All 21 are forced whether or not they appear below: an outlet discovery finds already
gets a plume and a catchment-normal size on its own, with or without an override. The `outlets`
list is eight overrides, keyed by `vassdragsnr`, and their almost-universal job is attaching a real
HydAPI gauge — `discharge_station`, sometimes `temperature_station` — to a mouth discovery already
found, so it is forced by an observed series instead of by its catchment normal; changing
`plume_depth` or splitting a discharge with `discharge_fraction` is the exception three of the
eight need (Glomma's two mouths, Drammenselva), not what the list is for. The other thirteen
discovered mouths run freshwater-only, sized by their own catchment. Only the gauge half needs a
credential, an NVE HydAPI key in `NVE_API_KEY` (free at https://hydapi.nve.no/Users); the map
services are open, and nothing else here needs one either.

The river config gets the same `data_root` as the rest of the setup, so the cached NVE responses
land alongside it rather than being shared from elsewhere, and it writes `forcing_rivers_nve.nc`
rather than `forcing_rivers.nc` so the OF800 file this setup used to produce is left intact for
comparison. `simulation_forcing_path` returns whichever file the named river config points at, so
that is the one the run reads.

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

    # FjordConfig is what the driver-level generics (run_simulation(), download_forcing(), etc.)
    # dispatch on — one method per subcommand, shared by every setup. Each driver then calls the
    # finer-grained hooks the nested configs below overload.
    return FjordConfig(
        # EvenGrid overloads the two grid hooks — domain_grid(config, architecture), which every
        # prepare_* pipeline calls, and simulation_grid(config, bathymetry_file, architecture),
        # which build_simulation reads the model grid through — both dispatching on this config's
        # type. Implemented here by wrapping LatitudeLongitudeGrid / ImmersedBoundaryGrid.
        grid_config = EvenGrid(
            size      = (240, 520, 24),
            halo      = (7, 7, 7),
            longitude = (10.2, 11.02),
            latitude  = (59.0, 59.93),
            # A geometric stretch, ratio 1.25 capped at 33.5 m, from a 1 m surface cell to a
            # deepest face at -400 m — the deepest sounding is 395.1 m, so nothing is clipped and
            # nothing is wasted. The predecessor's 18 levels ran to -450 m, where the first layer
            # held no water at all, and reached the seabed in 25 m and 50 m steps: the median
            # column's floor sat in a 25 m layer and 61% of columns were floored in a layer thicker
            # than 20 m, so the deepest 25 m of a column was one cell with no vertical structure.
            # That is where the tracer extremes of the 2020 run sat — 96% of them in a bottom cell,
            # almost all at the two 25 m layers. Here the median floor layer is 9 m, 24% of columns
            # exceed 20 m, and no two adjacent layers differ by more than a factor 1.33 (it was
            # 2.50, at the 25 m to 10 m step).
            z_faces   = [
                -400.0, -366.5, -333.0, -299.5, -266.0, -232.5, -199.0, -165.5,
                -132.0, -105.0, -83.0, -66.0, -52.0, -41.0, -32.0, -25.0,
                -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0,
            ],
        ),
        # DybdedataConfig overloads bathymetry_dataset (required), plus regrid_options and
        # smoothing_options (both optional) — hooks prepare_bathymetry calls, dispatching on this
        # config's type.
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
            # Flood land within five rows of the open southern edge, so the boundary condition is
            # not reconciling a radiation scheme, a prescribed exterior state and a headland within
            # a couple of cells of each other. On this domain the band is nearly clear already and
            # it moves seven cells, none in the boundary row itself — cheap insurance rather than a
            # correction. The open edges come from `boundary_config` below.
            open_boundary_land_cells = 5,
            # Six and no more. Every remaining interior land component is 7 cells or larger with a
            # minimum width of 2 to 4 cells — compact skerries, not the one-cell ridges this stage
            # exists to remove. Raising it would delete real topography.
            max_island_cells      = 6,
            close_narrow_passages = true,
            spike_ratio           = 0.5,
            minimum_cell_fraction = 0.2,
            # 0.25, not the 0.5 this used to be. Measured on this bathymetry, tightening the
            # limit costs almost nothing: the Drøbak sill goes from 65.2 m to 64.6 m, the deepest
            # point is untouched at 395.1 m and water volume changes by 0.000%, while adjacent
            # pairs steeper than r = 0.3 go from 4620 to none. The worry that a tight limit
            # flattens a genuine sill does not survive contact with this domain.
            max_slope_factor      = 0.25,
            geonorge_cache        = true,
            regrid_cache          = false,
        ),
        # NorKystConfig overloads forcing_time_steps, forcing_source_grid and
        # forcing_variable_names — hooks prepare_forcing calls — plus download_forcing itself,
        # which the driver-level download_forcing(config::FjordConfig) dispatches straight to.
        # All dispatch on this config's type. simulation_forcing is left at its default, since
        # this is the FjordSim NetCDF forcing contract build_simulation already reads.
        forcing_config = NorKystConfig(
            data_root        = data_root,
            output_directory = "norkyst",
            output_file      = "forcing.nc",
            plot_file        = "forcing.png",
            architecture     = :auto,
            parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
            years            = [2020],
            # NVERiversConfig overloads river_locations, river_series and download_rivers (all
            # required, hooks add_rivers calls), plus river_minimum_levels, river_plume_depth and
            # river_lambdas (all optional). All dispatch on this config's type. It replaced
            # OF800RiversConfig, which read a fixed published artifact: this one queries NVE, so
            # the window follows `years` and every value has provenance. The HydAPI half needs an
            # API key — `NVE_API_KEY`, free at https://hydapi.nve.no/Users.
            #
            # `standalone = false`, so `add_rivers` patches a copy of the NorKyst forcing above
            # rather than writing a river-only file; `initial_conditions = FromForcing()` below
            # therefore still has a 3D ocean state to read.
            rivers           = NVERiversConfig(
                data_root  = data_root,
                # Not `forcing_rivers.nc`: the OF800 file of that name is left intact so the two
                # sources can be compared. `simulation_forcing_path` returns whichever the named
                # river config points at, so this is the file the run reads.
                output_file = "forcing_rivers_nve.nc",
                plot_file   = "forcing_rivers_nve.png",
                years       = [2020],
                # Every outlet below is *discovered* from NVE's ELVIS river network and sized by
                # its REGINE catchment, so this setup states no river coordinate at all. It used
                # to state nine, copied from the OF800 dataset's outlet table, and measured
                # against NVE's own mouths they were wrong in ways nothing checked: Mosseleva 7 km
                # and Gjersjøelva 6 km from their rivers, Drammenselva 602 m from the nearest
                # water in a column at the 2 m `minimum_depth` floor, Glomma on a shelf beside
                # both of its beds rather than in either, and Glomma's whole western arm missing.
                #
                # 0.5 m³/s gives 21 mouths carrying 99.6 % of the domain's 1179 m³/s. Lowering it
                # to 0.15 would add twelve more streams for the remaining 0.3 %.
                minimum_discharge = 0.5,
                # 5 m is a surface plume — 4 cells, 5.5 m, on this vertical grid. The two rivers
                # whose model cell is inside a river bed rather than an estuary say `Inf` below.
                default_plume_depth = 5.0,
                # Zero, where the OF800 config used five. The two are alternatives, not a pair:
                # `minimum_levels` fixes a shallow outlet by *relocating* it to a column deep
                # enough to resolve an estuarine exchange, `river_plume_depth` by *filling* the
                # column so there is no partial cell left to concentrate salt in. Every outlet
                # has a plume, so relocation would only move rivers away from their real mouths
                # for a problem the plume has already solved.
                minimum_levels = 0,
                # 600 s is the floor on the timescale `river_lambdas` derives from discharge, and
                # it is a stability bound rather than a preference: λ = Q̄/V is the physically
                # correct dilution rate, but at Drammensfjord-scale cells a large river's peak
                # discharge reaches λ·Δt > 1, and `ForcingFromFile` reads λ > 1 as an x-flux. On
                # this grid it binds for Glomma and Drammenselva and leaves the small streams at
                # their own rate — a spread one shared timescale cannot express at all.
                minimum_relaxation_timescale = 600.0,
                # Overrides, keyed by the `vassdragsnr` of each mouth's terminal ELVIS segment.
                # They never add a river — a mouth discovery finds is forced with or without one —
                # and an override matching no discovered mouth is an error. An override's usual job
                # is attaching a gauge (`discharge_station`, sometimes `temperature_station`) to a
                # mouth that would otherwise run on its catchment normal alone; only three of the
                # eight below also change `plume_depth` or `discharge_fraction`, which they need for
                # geometry reasons (a river bed rather than an estuary, a split mouth), not because
                # that is what an override is for. The eleven mouths named nowhere below —
                # Aulivassdraget, Årosvassdraget, Hølenelva, Alna, Istreelva, Askerelva,
                # Årungelva, Sageneelva, Borreelva, Selvikelva and Ljanselva — run freshwater-only
                # at 5 m with λ from their catchment normal, which is the whole point of
                # discovery: a river with no gauge still scales by its own size.
                outlets = [
                    # Glomma, Norway's largest catchment, reaching the sea through **two** mouths
                    # at Fredrikstad: Østerelva east of Kråkerøy and Vesterelva west of it. Both
                    # beds are resolved on this grid, cells (222, 99) and (198, 99). ELVIS files
                    # Vesterelva under the small Seutelva catchment, so discovery finds it as a
                    # 1.7 m³/s stream; these two overrides put Solbergfoss's observed series
                    # across both instead, and the 2/3–1/3 split is a stated assumption — NVE's
                    # own Nedre Glomma flood report describes the division and publishes no
                    # fraction, and none was found elsewhere.
                    #
                    # Discharge at Solbergfoss and temperature 40 km downstream at Sarpfossen,
                    # because no station carries both.
                    #
                    # `Inf`: both model cells are inside the river channel rather than an estuary,
                    # so the whole wet column is river and is relaxed to S = 0. It asks for the
                    # column that is there — 4.2 m east, 13.4 m west — not for a fixed depth.
                    NVERiver(
                        vassdragsnr = "002.A21", name = "Glomma (Osterelva)",
                        discharge_station = "2.605.0", temperature_station = "2.1087.0",
                        discharge_fraction = 2 // 3, plume_depth = Inf,
                    ),
                    NVERiver(
                        vassdragsnr = "002.2A", name = "Glomma (Vesterelva)",
                        discharge_station = "2.605.0", temperature_station = "2.1087.0",
                        discharge_fraction = 1 // 3, plume_depth = Inf,
                    ),
                    # Drammenselva, the second largest input. Mjøndalen bru carries discharge and
                    # temperature on one id at 4 m a.s.l., ~5 km above the fjord. `Inf` for the
                    # same reason as Glomma — the cell is river bed.
                    NVERiver(
                        vassdragsnr = "012.A2", name = "Drammenselva",
                        discharge_station = "12.534.0", temperature_station = "12.534.0",
                        plume_depth = Inf,
                    ),
                    # Mossevassdraget at Moss dam, 693 km². No water temperature series.
                    NVERiver(vassdragsnr = "003.A4", discharge_station = "3.23.0"),
                    # Nordmarkvassdraget is Akerselva's catchment. Temperature deliberately
                    # omitted: the 2020 series averages 21.2 °C and peaks at 31.1 °C, which is a
                    # dry or sun-exposed sensor, not an Oslo river. The row comes back as NaN,
                    # which `ForcingFromFile` reads as its land sentinel — so the river still
                    # freshens its cells and simply does not force their temperature.
                    NVERiver(
                        vassdragsnr = "006.A10", name = "Akerselva",
                        discharge_station = "6.38.0",
                    ),
                    # Lysakerelva. Temperature omitted for the same reason and worse: mean
                    # 19.6 °C, max 37.2 °C in 2020.
                    NVERiver(vassdragsnr = "007.A0", discharge_station = "7.29.0"),
                    # Sandvikselva, 226 km², and the only small stream here whose temperature
                    # survives inspection: mean 8.4 °C, range 1.0-20.9 °C in 2020.
                    NVERiver(
                        vassdragsnr = "008.A2",
                        discharge_station = "8.2.0", temperature_station = "8.2.0",
                    ),
                    # Lierelva at Oppsal, 223 km² of a 310 km² catchment. No temperature series.
                    NVERiver(vassdragsnr = "011.A0", discharge_station = "11.6.0"),
                ],
            ),
        ),
        # The exterior state along the open southern edge, from the *hourly* NorKyst collection: a
        # Flather boundary compares the model's own η against the exterior one, and the daily means
        # the forcing config reads have the tide averaged out of them. Only the boundary row is
        # prepared at that cadence, so the file is a few hundred MB.
        #
        # A `FjordConfig` field of its own, independent of the forcing above: this is a separate
        # file from a separate collection read by a separate pipeline, and it is what states the open
        # edge — for the boundary steps, for the forcing land mask, and for the boundary condition.
        # NorKystBoundariesConfig overloads boundary_time_steps, boundary_source_grid,
        # boundary_variable_names (all required) and download_boundaries — hooks
        # prepare_boundaries and download_boundaries call, dispatching on this config's type.
        boundary_config = NorKystBoundariesConfig(
            data_root        = data_root,
            output_directory = "norkyst_hourly",
            output_file      = "boundaries.nc",
            plot_file        = "boundaries.png",
            # One `Symbol` for one open edge; a collection for several — a region in the open
            # ocean writes `(:south, :north, :west, :east)` and gets one boundary file with four
            # sides in it. Oslofjord connects to the Skagerrak through its southern edge alone.
            open_edges       = :south,
            margin           = 0.05,
            architecture     = :auto,
            parameters       = [
                "temperature", "salinity", "u_eastward", "v_northward", "zeta", "ubar", "vbar",
            ],
            years            = [2020],
        ),
        # NORA3Config overloads atmosphere_time_steps, atmosphere_source_grid and
        # atmosphere_variable_names — hooks prepare_atmosphere calls — plus download_atmosphere
        # (called by the driver-level download_atmosphere) and prescribed_atmosphere /
        # prescribed_radiation (called directly by build_simulation, since the setup is
        # simulated). All dispatch on this config's type.
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
            # CoupledHydrostaticSimulation overloads coupled_simulation and model_tracers — the
            # model hooks build_simulation calls, dispatching on this config's type.
            model              = CoupledHydrostaticSimulation(
                buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
                # A BoundarySponge rather than a bare closure tuple: it adds a harmonic viscosity
                # and diffusivity ramping to zero over 16 cells (~3.2 km) inward from whichever
                # edges `boundary_config` opens, and passes `base` through untouched everywhere
                # else. Without it nothing damped what the open boundary radiates — the 2020 run
                # carried velocities of 0.4 m/s std and 2 m/s peak on the boundary row against
                # 0.13 and 0.43 in the interior, and grid-scale salinity roughness 70x the interior
                # value, which accumulated for 50 days in the near-boundary bottom cells and then
                # ran away. Viscous rather than a tracer relaxation band, so it cannot fight the
                # open boundary's own nudging towards the same data.
                #
                # 30 m² s⁻¹ is set by explicit-diffusion stability, not by taste: Δt ≤ Δx²/4ν is
                # 310 s on this 193 m cell, comfortably clear of `max_time_step` below. Raising it
                # much past 50 would cap the time step instead of the CFL doing it.
                closure            = BoundarySponge(
                    # `minimum_tke` is not a safety net on this grid — it *is* the vertical
                    # closure over most of the column. CATKE takes `w★ = sqrt(max(minimum_tke, e))`,
                    # and measured on the 2020 run at day 20.5 the prognostic `e` sits below this
                    # floor in 85 % of wet cells: 89 % below 6 m, 95-100 % below 9 m. So κ there is
                    # the documented background κ = Cʰⁱᶜ·e_min/N = 0.098·e_min/N with
                    # ν = 0.242·e_min/N, which the domain-median κᵤ/κᶜ of 2.38 against
                    # Cʰⁱᵤ/Cʰⁱᶜ = 2.47 confirms.
                    #
                    # 7e-6 is ROMS' `GLS_KMIN` default (7.6e-6), and it does *not* carry over
                    # unchanged: ROMS floors its length scale too (`GLS_PMIN`) where CATKE's is
                    # diagnostic and unfloored, so the same TKE floor buys a much larger κ. It is
                    # kept because it lands in the right place for *this* fjord, not because it came
                    # from ROMS. Median κᶜ at 94-118 m is 8-9e-5 m² s⁻¹ against basin-mean values
                    # measured here by density budget — 1.0-1.8e-4 in Bunnefjorden, 4.9-7.6e-4 in
                    # Vestfjorden, 2.5e-3 just inside the Drøbak sill (Gade 1970; Staalstrøm et al.
                    # 2012, Ocean Sci. 8, 525) — so the basin water is at the low end already and
                    # lowering the floor would slow the deep-water renewal this domain is judged on.
                    # The fjord relation is κ ∝ N^-1.5 (Stigebrandt & Aure 1989) where the floor
                    # gives κ ∝ N^-1, so one number cannot fit both the pycnocline and the basin;
                    # the basin is the one to fit.
                    #
                    # `Cᵇ` stays at Oceananigans' 0.28 rather than NumericalEarth's regional 0.01
                    # (`default_ocean_closure`), deliberately: it scales the near-bottom mixing
                    # length, and Stigebrandt attributed the high basin-mean diffusivity here to
                    # exactly that boundary mixing.
                    base = (
                        CATKEVerticalDiffusivity(minimum_tke = 7e-6),
                        # 2e4 m⁴ s⁻¹, down from 1e5. Biharmonic damping of the 2Δx mode goes as
                        # ν₄·16/Δx⁴, an e-folding of 14.5 min at 1e5 on this 193 m cell — the commit
                        # that raised it from 15 aimed at "~1.7 hours" and was out by 7x. 2e4 gives
                        # the 72 min it meant, keeps the scale selectivity that leaves an 8 km
                        # baroclinic eddy alone (56 h at 8Δx), and is still ~1300x what it replaced.
                        # For scale, Norkyst-800 — the ROMS system feeding this run — applies no
                        # explicit interior viscosity at all and a 10 m² s⁻¹ harmonic tracer
                        # diffusivity, where κ = 2e3 here is 0.2 m² s⁻¹ at 2Δx.
                        #
                        # It also matters for a stability limit nothing else watches: the explicit
                        # biharmonic needs Δt ≤ Δx⁴/32ν₄, which is 434 s at 1e5 and 2170 s here, and
                        # `AdaptiveTimeStep` measures only the advective CFL
                        # (`cell_advection_timescale_coupled_model`).
                        HorizontalScalarBiharmonicDiffusivity(ν = 2e4, κ = 2e3),
                    ),
                    width_cells = 16,
                    viscosity   = 30.0,
                    diffusivity = 15.0,
                ),
                # One scheme for every tracer, not a NamedTuple naming T and S. Oceananigans gives
                # any tracer a NamedTuple omits the `Centered()` default, and CATKE contributes an
                # `e` the setup never names — so `e` was being advected by an unbounded centered
                # scheme, which makes the vertical diffusivity noisy exactly where the bottom cells
                # were failing. A scalar covers whatever the closure adds, now and later.
                tracer_advection   = WENO(),
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
            # edge — which edge, and whose exterior state, come from `boundary_config` above.
            # Dropping it would close the domain; argument order is merge precedence.
            #
            # The two timescales are Marchesiello et al. (2001): nudge hard towards the data on
            # inflow, let radiation do the work on outflow.
            boundary_conditions = MergedBoundaryConditions(
                AirSeaFluxes(),
                QuadraticBottomDrag(coefficient = 0.003),
                OpenLateralBoundaryFromData(
                    # 3 hours, not the 1 day this was. At the time steps this run actually takes
                    # (~10 s), a one-day timescale relaxes by 1.2e-4 per step against an advective
                    # rate into the boundary cell of ~2.5e-3 s⁻¹ — 200x weaker, which is why the
                    # nudging never caught the drift it was there to catch. Both numbers are
                    # tuning knobs: too strong an inflow nudge over-constrains an open boundary
                    # and reflects.
                    inflow_timescale  = 3hours,
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
                    # `e` is CATKE's TKE tracer, in for diagnosis rather than for science: the
                    # 2026-08-25 run died at day 11.5 on a TEOS10 `sqrt(Sᴬ + 32)` DomainError —
                    # one cell below -32 psu — and the file carried nothing that could show what
                    # led there. `e` drives the vertical diffusivity, so it is the field most
                    # likely to hold the precursor.
                    #
                    # `η` belongs here too and cannot go here. Oceananigans' NetCDF writer cannot
                    # write a `(Center, Center, Nothing)` *user output* at all: the free surface
                    # asks for a singleton `z_aaf = [0.0]` while the grid's own vertical coordinate
                    # is already `z_aaf` with the 25 faces, and the writer raises rather than
                    # reconciling them. Measured: it fails with the 3D fields, in a file of its own,
                    # and with `include_grid_metrics = false`. `bottom_height` is the same location
                    # and *is* written, because grid metrics take a different path — so this is an
                    # upstream bug in the user-output path, not a configuration mistake.
                    variables          = (:T, :S, :u, :v, :e),
                    interval           = 3hour,
                    overwrite_existing = true,
                ),
                # `η` in JLD2 because it cannot go in the NetCDF file: Oceananigans' NetCDF writer
                # cannot emit a `(Center, Center, Nothing)` user output at all — the free surface
                # asks for a singleton `z_aaf = [0.0]` against the grid's own 25-face `z_aaf`, and
                # it fails beside the 3D fields, alone in its own file, and with grid metrics off.
                # It is what the Flather boundary and the barotropic solver actually exchange, so
                # it is worth a second file rather than going unwritten. See `FieldSnapshotWriter`.
                FieldSnapshotWriter(
                    name               = :surface,
                    output_file        = "snapshots_surface.jld2",
                    variables          = (:η,),
                    interval           = 3hour,
                    overwrite_existing = true,
                ),
                CheckpointWriter(interval = 12hours, cleanup = true),
            ),
            # Overloads attach_callback! — what the run reports while it runs. `report` is the
            # function itself, so a model whose tracers omit :T (which `progress` reads) names its
            # own here instead. An empty tuple runs silently.
            callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
            # Overloads attach_time_stepping! and initial_time_step.
            time_stepping = AdaptiveTimeStep(
                initial_time_step    = 1second,
                cfl                  = 0.3,
                max_time_step        = 3minutes,
                max_time_step_change = 1.01,
            ),
            # FromForcing(): the NorKyst state at `start_date` rather than a uniform water
            # column: every tracer the model names plus u and v, whichever of them the forcing
            # file carries. A literal NamedTuple (`(T = 5.0, S = 33.0)`) still works, and
            # `FromResults("snapshots_ocean_<tag>.nc")` continues from a previous run instead.
            initial_conditions = FromForcing(),  # (T = 5.0, S = 33.0),
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
