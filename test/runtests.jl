using FjordSim
using FjordSim.Bathymetry: write_bathymetry_file
using FjordSim.Configs: open_edges
using Dates: DateTime, Hour
using Test
using ArchGDAL
using NCDatasets
using Oceananigans
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: AbstractField
using NumericalEarth

# Shared stub configs, fixture factories and the run-tag helper. Included at top level because a
# @testset body is a `let` scope, where neither `struct` nor `const` may appear.
include("utilities.jl")

@testset "Public interface" begin
    @testset "exports" begin
        exported_symbols = [
            :ImmersedBoundaryGrid,
            :FjordConfig,
            :AbstractGridConfig,
            :AbstractBathymetryConfig,
            :AbstractForcingConfig,
            :AbstractRiverConfig,
            :AbstractSimulationConfig,
            :EvenGrid,
            :DybdedataConfig,
            :NorKystConfig,
            :OF800RiversConfig,
            :forcing_from_file,
            :simulation_forcing,
            :forcing_path,
            :forcing_directory,
            :river_forcing_path,
            :plot_path,
            :bathymetry_path,
            :prepare_bathymetry,
            :prepare_forcing,
            :download_forcing,
            :add_rivers,
            :download_rivers,
            :interpolation_architecture,
            :plot_bathymetry,
            :plot_forcing,
            :plot_boundaries,
            :download_boundaries,
            :prepare_boundaries,
            # extension hooks a new config subtype overloads
            :bathymetry_dataset,
            :regrid_options,
            :forcing_time_steps,
            :forcing_source_grid,
            :forcing_variable_names,
            :forcing_monthly_filename,
            :river_locations,
            :river_series,
            :river_search_radius,
            :boundary_time_steps,
            :boundary_source_grid,
            :boundary_variable_names,
            :boundary_date_range,
            :boundary_series,
            :boundary_data_path,
            :boundary_data_directory,
            :NorKystBoundariesConfig,
            :AbstractBoundaryDataConfig,
            :ProjectedSourceGrid,
            :RiverLocation,
            :geodatabase_path,
            :top_bottom_boundary_conditions,
            :air_sea_flux_boundary_conditions,
            :quadratic_bottom_drag_boundary_conditions,
            :boundary_condition_sides,
            :field_boundary_conditions,
            :AirSeaFluxes,
            :QuadraticBottomDrag,
            :OpenLateralBoundaryFromData,
            :MergedBoundaryConditions,
            :domain_grid,
            :simulation_grid,
            :forcing_date_range,
            :ProgressCallback,
            :attach_callback!,
            :coupled_simulation,
            :recursive_merge,
            :progress,
            :cell_advection_timescale_coupled_model,
            :NORA3PrescribedAtmosphere,
            :NORA3PrescribedRadiation,
            :MultiYearNORA3,
            :AbstractAtmosphereConfig,
            :prepare_atmosphere,
            :download_atmosphere,
            :plot_atmosphere,
            :atmosphere_path,
            :atmosphere_directory,
            :atmosphere_time_steps,
            :atmosphere_source_grid,
            :atmosphere_variable_names,
            :atmosphere_target_axes,
            :ProjectedAtmosphereGrid,
            :AtmosphereRecord,
            :NORA3Config,
            :fjord_config,
            :setup_names,
            :oslofjorden,
            :drammensfjorden,
            # simulation
            :SimulationConfig,
            :AbstractCoupledSimulationConfig,
            :AbstractFreeSurfaceConfig,
            :AbstractBoundaryConditionConfig,
            :AbstractWriterConfig,
            :AbstractTimeSteppingConfig,
            :CoupledHydrostaticSimulation,
            :SplitExplicitFreeSurfaceConfig,
            :free_surface,
            :SnapshotWriter,
            :CheckpointWriter,
            :AdaptiveTimeStep,
            :model_tracers,
            :attach_writer!,
            :attach_time_stepping!,
            :initial_time_step,
            :build_simulation,
            :run_simulation,
            :simulation_architecture,
            :results_path,
            :prescribed_atmosphere,
            :prescribed_radiation,
        ]

        for sym in exported_symbols
            @test isdefined(FjordSim, sym)
        end

        submodule_exports = [
            (:Utils, :compute_faces),
            (:Utils, :safe_execute),
            (:Utils, :extract_z_faces),
            (:Utils, :netcdf_to_jld2),
            (:Utils, :save_fts),
            (:Datasets, :ForcingDataset),
            (:Datasets, :ResultsDataset),
            (:Datasets, :last_date),
            (:Setups, :fjord_config),
            (:Setups, :setup_names),
            (:Simulations, :SimulationConfig),
            (:Simulations, :run_simulation),
            # Moved out of FjordSim.jl into Simulations; still re-exported, so the entry above in
            # `exported_symbols` and this one together pin both halves of that move.
            (:Simulations, :coupled_simulation),
            (:BoundaryConditions, :boundary_condition_sides),
            (:CLI, :parse_arguments),
        ]

        # `getfield` already errors on a missing module, so the loop covers every module it names.
        for (module_name, sym) in submodule_exports
            @test isdefined(getfield(FjordSim, module_name), sym)
        end

        @test isdefined(FjordSim, :Bathymetry)   # the one module no export above reaches through
        @test parentmodule(coupled_simulation) === FjordSim.Simulations
        # The per-piece boundary hook used to be called `boundary_conditions`, which collided with
        # `Oceananigans.Fields.boundary_conditions` and so was the one hook FjordSim could not
        # re-export. Renamed to `boundary_condition_sides`, it is exported like every other hook,
        # and the old name is gone rather than deprecated.
        @test :boundary_condition_sides ∈ names(FjordSim)
        @test :boundary_conditions ∉ names(FjordSim)
        @test !isdefined(FjordSim.BoundaryConditions, :boundary_conditions)
    end

    @testset "boundary condition signatures" begin
        mktempdir() do tmp
            # A legacy-format bathymetry file — positive depths, `lat` on the x dimension — since the
            # loader still accepts one, and nothing else here exercises that path.
            bathymetry_file = joinpath(tmp, "bathymetry.nc")
            ds = NCDataset(bathymetry_file, "c")
            defDim(ds, "x", 2)
            defDim(ds, "y", 2)
            defDim(ds, "zf", 3)
            z_faces = defVar(ds, "z_faces", Float64, ("zf",))
            h = defVar(ds, "h", Float64, ("x", "y"))
            lat = defVar(ds, "lat", Float64, ("x",))
            lon = defVar(ds, "lon", Float64, ("y",))
            z_faces[:] = [-20.0, -10.0, 0.0]
            h[:, :] = [10.0 10.0; 10.0 10.0]
            lat[:] = [59.0, 60.0]
            lon[:] = [10.0, 11.0]
            close(ds)

            boundary_conditions =
                @test_nowarn top_bottom_boundary_conditions(;
                    grid = ImmersedBoundaryGrid(bathymetry_file, CPU(), (1, 1, 1)),
                    bottom_drag_coefficient = 0.003,
                )

            # Velocities get a condition at both ends, tracers only at the top. Asserting the condition
            # is not `nothing` subsumes a `haskey`: reaching a field that is not there is an error.
            for name in (:u, :v), side in (:top, :bottom)
                @test !isnothing(getproperty(getproperty(boundary_conditions, name), side))
            end
            for name in (:T, :S)
                @test !isnothing(getproperty(boundary_conditions, name).top)
            end

            # Each velocity also gets a stress on the *immersed* seabed, which is the half that
            # actually acts on a fjord: a `bottom` condition applies at the underlying grid's floor,
            # and `u_quadratic_bottom_drag` reads `Φ.u[i, j, 1]` there. This fixture is itself a case
            # in point — its bottom sits at -10 m with faces at -20, -10 and 0, so level k = 1 is
            # immersed and the `bottom` condition alone would put drag nowhere at all.
            for name in (:u, :v)
                immersed = getproperty(boundary_conditions, name).immersed
                @test immersed isa ImmersedBoundaryCondition
                @test !isnothing(immersed.bottom)
            end

            # The tracer top conditions must carry a freshwater exchange NumericalEarth can read back
            # out of them: `net_fluxes` pulls the volume flux and its heat content straight from the
            # boundary condition, so a bare FluxBoundaryCondition raises a MethodError as soon as the
            # coupled model assembles fluxes — which nothing else here would catch.
            extract_freshwater_flux = NumericalEarth.Oceans.extract_freshwater_flux
            for name in (:T, :S)
                condition = getproperty(boundary_conditions, name).top.condition
                @test extract_freshwater_flux(condition) isa AbstractField
            end
            exchange = NumericalEarth.Oceans.freshwater_exchange(boundary_conditions.T.top.condition)
            @test exchange.content_flux isa AbstractField
        end
    end

    @testset "open lateral boundary" begin
        mktempdir() do tmp
            (; grid) = immersed_test_grid(joinpath(tmp, "bathymetry.nc"); size = (2, 3, 2))
            boundaries_config = test_boundaries_config(data_root = tmp)
            write_prepared_boundaries(
                FjordSim.boundary_data_path(boundaries_config);
                size = (2, 3, 2),
                edges = :south,
                value = (name, index) -> name == "T" ? 6.0f0 : name == "S" ? 30.0f0 : 0.5f0,
            )
            boundaries = FjordSim.boundary_series(boundaries_config, grid, DateTime(2020, 1, 1))

            # Every part of the state gets its own condition. The predecessor of this config put a
            # closed wall on the normal velocity and opened only the tracers, so the barotropic `V`
            # is the row that says the boundary is genuinely open.
            sides = boundary_condition_sides(
                OpenLateralBoundaryFromData(inflow_timescale = 1day, outflow_timescale = 360days),
                grid, NamedTuple(), boundaries_config, (:T, :S), boundaries,
            )
            @test Set(keys(sides)) == Set((:v, :u, :V, :T, :S))
            for field in (:v, :u, :V, :T, :S)
                @test keys(getproperty(sides, field)) == (:south,)
            end

            # The normal velocity radiates rather than being walled off, and the tracers and the
            # tangential velocity share its scheme.
            @test sides.v.south.classification isa
                  Oceananigans.BoundaryConditions.NormalFlow{<:Oceananigans.BoundaryConditions.NormalRadiation}
            @test sides.T.south.classification isa
                  Oceananigans.BoundaryConditions.Value{<:Oceananigans.BoundaryConditions.NormalRadiation}
            @test sides.u.south.classification isa
                  Oceananigans.BoundaryConditions.Value{<:Oceananigans.BoundaryConditions.NormalRadiation}
            @test sides.v.south.classification.scheme.inflow_timescale == 1day
            @test sides.v.south.classification.scheme.outflow_timescale == 360days

            # The barotropic transport is a Flather condition, and its condition is the series
            # itself for the velocities and tracers — a reduced `FieldTimeSeries` is a boundary
            # condition Oceananigans indexes by the two tangential indices on its own.
            @test sides.V.south.classification isa
                  Oceananigans.BoundaryConditions.NormalFlow{<:Oceananigans.BoundaryConditions.GravityWaveRadiation}
            @test sides.T.south.condition === boundaries.south.T
            @test sides.v.south.condition === boundaries.south.v

            # The stored transport function actually reads the prepared exterior state and turns it
            # into a transport with the model's own column depth — the check that closes the gap
            # between "a scheme is attached" and "it reads the right data".
            clock = (; time = 0.0)
            η = Field{Center,Center,Nothing}(grid)
            transport = sides.V.south.condition
            depth = Oceananigans.Grids.column_depthᶜᶠᵃ(1, 1, grid.Nz + 1, grid, η)
            exterior = transport.func(1, 1, grid, clock, (; η), transport.parameters)
            @test all(isapprox.(exterior, (0.5 * depth, 0.5)))

            # Each edge puts its conditions on the field normal to it, and nowhere else.
            for (edge, normal, tangential, barotropic) in (
                (:north, :v, :u, :V), (:west, :u, :v, :U), (:east, :u, :v, :U),
            )
                edge_config = test_boundaries_config(data_root = tmp, open_edges = edge)
                write_prepared_boundaries(
                    FjordSim.boundary_data_path(edge_config); size = (2, 3, 2), edges = edge,
                )
                edge_boundaries = FjordSim.boundary_series(edge_config, grid, DateTime(2020, 1, 1))
                edge_sides = boundary_condition_sides(
                    OpenLateralBoundaryFromData(inflow_timescale = 1day, outflow_timescale = 360days),
                    grid, NamedTuple(), edge_config, (:T,), edge_boundaries,
                )
                @test Set(keys(edge_sides)) == Set((normal, tangential, barotropic, :T))
                @test keys(getproperty(edge_sides, normal)) == (edge,)
                @test keys(getproperty(edge_sides, barotropic)) == (edge,)
            end

            # A bad edge is rejected where it is stated, at construction, rather than surviving as far
            # as a silently empty tuple of conditions.
            @test_throws ArgumentError test_boundaries_config(data_root = tmp, open_edges = :middle)
            # ...and the `Val` fallback behind that check is itself an `ArgumentError` rather than a
            # missing method, which is what makes an out-of-tree caller's mistake legible.
            @test_throws ArgumentError FjordSim.Forcing.boundary_along_axis(Val(:middle))
            @test_throws ArgumentError FjordSim.BoundaryConditions.open_normal_velocity_boundary_conditions(
                Val(:middle), boundaries.south, nothing,
            )

            # Several edges at once: a domain in the open ocean opens all four, and every group lands
            # on every side. `u` picks up a normal condition on the west and east and a tangential one
            # on the south and north; each tracer picks up all four.
            all_edges = (:south, :north, :west, :east)
            open_ocean = test_boundaries_config(data_root = tmp, open_edges = all_edges)
            write_prepared_boundaries(
                FjordSim.boundary_data_path(open_ocean); size = (2, 3, 2), edges = all_edges,
            )
            ocean_series = FjordSim.boundary_series(open_ocean, grid, DateTime(2020, 1, 1))
            @test Set(keys(ocean_series)) == Set(all_edges)

            ocean_sides = boundary_condition_sides(
                OpenLateralBoundaryFromData(inflow_timescale = 1day, outflow_timescale = 360days),
                grid, NamedTuple(), open_ocean, (:T, :S), ocean_series,
            )
            @test Set(keys(ocean_sides)) == Set((:u, :v, :U, :V, :T, :S))
            for field in (:T, :S, :u, :v)
                @test Set(keys(getproperty(ocean_sides, field))) == Set(all_edges)
            end
            # The barotropic transports carry only the sides they are normal to, which is what pairs
            # Oceananigans' Chapman condition onto each of the four.
            @test Set(keys(ocean_sides.V)) == Set((:south, :north))
            @test Set(keys(ocean_sides.U)) == Set((:west, :east))
        end
    end

    @testset "open lateral boundary without prepared data" begin
        mktempdir() do tmp
            (; grid) = immersed_test_grid(joinpath(tmp, "bathymetry.nc"); size = (2, 3, 2))
            open_boundary =
                OpenLateralBoundaryFromData(inflow_timescale = 1day, outflow_timescale = 360days)

            # The exterior state is data, so a setup naming no boundary dataset is a stated error
            # rather than a boundary that quietly relaxes towards nothing. Two independent statements
            # are missing there — the config that names the edge, and the series read from its file —
            # so each has its own error.
            @test_throws ErrorException boundary_condition_sides(
                open_boundary, grid, NamedTuple(), nothing, (:T, :S), nothing,
            )

            # ...and so is a tracer the prepared file does not carry.
            boundaries_config = test_boundaries_config(data_root = tmp)
            write_prepared_boundaries(
                FjordSim.boundary_data_path(boundaries_config); size = (2, 3, 2),
            )
            boundaries = FjordSim.boundary_series(boundaries_config, grid, DateTime(2020, 1, 1))
            @test_throws ErrorException boundary_condition_sides(
                open_boundary, grid, NamedTuple(), nothing, (:T, :not_prepared), boundaries,
            )
            @test_throws ErrorException boundary_condition_sides(
                open_boundary, grid, NamedTuple(), boundaries_config, (:T, :not_prepared), boundaries,
            )
        end
    end

    @testset "boundary condition configs" begin
        mktempdir() do tmp
            (; grid) = immersed_test_grid(joinpath(tmp, "bathymetry.nc"); size = (2, 3, 2))
            boundaries_config = test_boundaries_config(data_root = tmp)
            forcing_config = test_forcing_config(data_root = tmp)
            write_prepared_forcing(forcing_path(forcing_config); size = (2, 3, 2))
            write_prepared_boundaries(
                FjordSim.boundary_data_path(boundaries_config); size = (2, 3, 2),
            )
            forcing = forcing_from_file(forcing_config; grid, tracers = (:T, :S))
            boundaries = FjordSim.boundary_series(boundaries_config, grid, DateTime(2020, 1, 1))

            # The surface and the bottom are separate pieces now, so each must contribute its own
            # half and nothing else — that separability is the whole point of the split.
            surface = boundary_condition_sides(
                AirSeaFluxes(), grid, forcing, boundaries_config, (:T, :S),
            )
            @test Set(keys(surface)) == Set((:u, :v, :T, :S))
            @test keys(surface.u) == (:top,)
            @test keys(surface.T) == (:top,)

            drag = boundary_condition_sides(
                QuadraticBottomDrag(coefficient = 0.003), grid, forcing, boundaries_config, (:T, :S),
            )
            @test Set(keys(drag)) == Set((:u, :v))
            @test Set(keys(drag.u)) == Set((:bottom, :immersed))

            # Together they are still exactly what `top_bottom_boundary_conditions` returns, which
            # is what the "boundary condition signatures" testset asserts the shape of.
            both = FjordSim.recursive_merge(surface, drag)
            @test Set(keys(both)) == Set((:u, :v, :T, :S))
            @test Set(keys(both.u)) == Set((:top, :bottom, :immersed))
            @test keys(both.T) == (:top,)

            # The open edge merges four groups into one contribution — normal velocity, tangential
            # velocity, barotropic transport and every tracer. That merge is what left
            # `build_simulation`, and nothing but a full model build would otherwise exercise it.
            open = boundary_condition_sides(
                OpenLateralBoundaryFromData(inflow_timescale = 1day, outflow_timescale = 360days),
                grid, forcing, boundaries_config, (:T, :S), boundaries,
            )
            @test Set(keys(open)) == Set((:v, :u, :V, :T, :S))
            @test keys(open.v) == (:south,)
            @test keys(open.T) == (:south,)

            # `AirSeaFluxes` and `QuadraticBottomDrag` ignore the exterior state, so they implement
            # the five-argument hook and reach the six-argument one through its forwarding
            # fallback — which is what keeps an out-of-tree boundary config working unchanged.
            @test boundary_condition_sides(
                AirSeaFluxes(), grid, forcing, boundaries_config, (:T, :S), boundaries,
            ) |> keys |> Set == Set((:u, :v, :T, :S))

            # The set config materializes every piece into `FieldBoundaryConditions`, which is what
            # a `HydrostaticFreeSurfaceModel` consumes. `v` gets its top, bottom *and* south, so the
            # merge is a merge rather than a last-wins overwrite of whole variables.
            materialized = field_boundary_conditions(
                MergedBoundaryConditions(
                    AirSeaFluxes(),
                    QuadraticBottomDrag(coefficient = 0.003),
                    OpenLateralBoundaryFromData(
                        inflow_timescale = 1day, outflow_timescale = 360days,
                    ),
                ),
                grid, forcing, boundaries_config, (:T, :S), boundaries,
            )
            @test materialized.v isa FieldBoundaryConditions
            @test !isnothing(materialized.v.top)
            @test !isnothing(materialized.v.bottom)
            @test materialized.v.south isa
                  BoundaryCondition{<:Oceananigans.BoundaryConditions.NormalFlow}
            @test materialized.T.south.classification.scheme isa
                  Oceananigans.BoundaryConditions.NormalRadiation
            # The barotropic transport is materialized under its own key, which the model
            # regularizes through `assumed_field_location` and hands to the free surface.
            @test materialized.V isa FieldBoundaryConditions
            @test materialized.V.south.classification.scheme isa
                  Oceananigans.BoundaryConditions.GravityWaveRadiation

            # Naming no pieces is a valid statement — the model's own defaults everywhere — rather
            # than an error, so a setup can opt out of every piece.
            @test field_boundary_conditions(
                MergedBoundaryConditions(), grid, forcing, boundaries_config, (:T, :S),
            ) == (;)

            # The set is a struct rather than a bare `Tuple`, which is what makes an alternative
            # merge strategy expressible at all.
            @test MergedBoundaryConditions(AirSeaFluxes()) isa AbstractBoundaryConditionSetConfig
        end
    end

    @testset "MultiYearNORA3 fields" begin
        mktempdir() do tmp
            nora3_filename = "NORA3_test.nc"
            write_nora3_stub(joinpath(tmp, nora3_filename))

            # The four fields the read side reaches for, and their types. Asserting the type subsumes a
            # `hasfield`: a missing field is an error rather than a `false`.
            nora3 = MultiYearNORA3(nora3_filename, tmp)
            for (name, expected) in (
                :metadata_filename => String,
                :default_download_directory => String,
                :size => Tuple,
                :all_dates => Vector,
            )
                @test getproperty(nora3, name) isa expected
            end
        end
    end

    @testset "NORA3 dataset interface" begin
        # The interface must sit on NumericalEarth's own generics so its machinery dispatches into this
        # backend instead of falling back to the generic `Metadata` defaults. Defined under a bare
        # `using NumericalEarth`, each one becomes a module-local binding that shadows the upstream
        # function, and the failure is silent: every fallback is plausible.
        #
        # The time axis has to be CF-encoded here, unlike in the struct-field test above: a `Metadatum`
        # is a `Metadata` whose date is an `AnyDateTime`, so a raw Float64 axis would not produce one
        # and `size` would fall through to the vector method and pass for the wrong reason.
        DW = NumericalEarth.DataWrangling
        NORA3 = FjordSim.Atmospheres.NORA3

        mktempdir() do tmp
            filename = "NORA3_interface.nc"
            dates = [DateTime(2020, 1, 1), DateTime(2020, 1, 1, 1), DateTime(2020, 1, 1, 2)]
            ds = NCDataset(joinpath(tmp, filename), "c")
            defDim(ds, "lon", 3)
            defDim(ds, "lat", 4)
            defDim(ds, "time", length(dates))
            defVar(ds, "time", dates, ("time",))
            temperature = defVar(ds, "air_temperature_2m", Float32, ("lon", "lat", "time"))
            temperature[:, :, :] = fill(280.0f0, (3, 4, length(dates)))
            close(ds)

            nora3 = MultiYearNORA3(filename, tmp)
            vector_metadata = DW.Metadata(:temperature; dataset = nora3, dates, dir = tmp)
            scalar_metadata = DW.Metadata(:temperature; dataset = nora3, dates = first(dates), dir = tmp)

            @test !DW.is_three_dimensional(vector_metadata)     # upstream default is `true`
            @test DW.dataset_variable_name(vector_metadata) == "air_temperature_2m"
            @test DW.available_variables(nora3) == NORA3.NORA3_variable_names
            @test DW.all_dates(nora3) == dates
            @test DW.first_date(nora3) == first(dates)
            @test DW.last_date(nora3) == last(dates)
            @test DW.metadata_filename(nora3) == filename
            @test Oceananigans.location(vector_metadata) == (Center, Center, Center)

            # A Metadatum carries one scalar date, so it needs its own `size`: the vector method's
            # `length(metadata.dates)` cannot count a scalar.
            @test scalar_metadata isa NORA3.NORA3Metadatum
            @test !(vector_metadata isa NORA3.NORA3Metadatum)
            @test size(vector_metadata) == (3, 4, length(dates))
            @test size(scalar_metadata) == (3, 4, 1)
        end
    end

    @testset "save_fts" begin
        # `save_fts` takes the on-disk locations from the in-memory series it is handed; it used to
        # name three undefined variables, so every call was an UndefVarError.
        mktempdir() do tmp
            grid = RectilinearGrid(CPU(); size = (2, 2, 2), x = (0, 1), y = (0, 1), z = (0, 1))
            times = [0.0, 3600.0]
            jld2_filepath = joinpath(tmp, "series.jld2")

            for location_types in ((Center, Center, Center), (Face, Center, Nothing))
                fts = FieldTimeSeries{location_types...}(grid, times)
                rm(jld2_filepath; force = true)
                FjordSim.Utils.save_fts(;
                    jld2_filepath,
                    fts_name = "T",
                    fts,
                    grid,
                    times,
                    boundary_conditions = fts.boundary_conditions,
                )
                @test isfile(jld2_filepath)
                reloaded = FieldTimeSeries(jld2_filepath, "T")
                @test Oceananigans.location(reloaded) == location_types
                @test reloaded.times == times
            end
        end
    end
end

@testset "Configs" begin
    @testset "required keywords and path resolution" begin
        data_root = joinpath(tempdir(), "fjordsim_config_test")

        # A setup must name its own output files; only the Geonorge FileGDB name is defaulted.
        @test_throws UndefKeywordError DybdedataConfig(data_root = data_root)
        @test_throws UndefKeywordError DybdedataConfig(data_root = data_root, output_file = "b.nc")
        @test_throws UndefKeywordError DybdedataConfig(output_file = "b.nc", plot_file = "b.png")

        bathymetry_config = DybdedataConfig(
            data_root = data_root,
            output_file = "bathymetry_105to232.nc",
            plot_file = "bathymetry_105to232.png",
            padding_cells = 4,
        )
        # Chosen names resolve against data_root...
        @test bathymetry_path(bathymetry_config) == joinpath(data_root, "bathymetry_105to232.nc")
        @test plot_path(bathymetry_config) == joinpath(data_root, "bathymetry_105to232.png")
        @test geodatabase_path(bathymetry_config) ==
              joinpath(data_root, "Basisdata_0000_Norge_25833_Dybdedata_FGDB.gdb")
        @test isdir(bathymetry_config.raw_directory)  # scratch space created by __init__

        # ...while an absolute path overrides data_root, so one FileGDB copy can be shared.
        shared = test_bathymetry_config(
            data_root = data_root,
            geodatabase_file = "/shared/Dybdedata.gdb",
        )
        @test geodatabase_path(shared) == "/shared/Dybdedata.gdb"
        @test bathymetry_path(shared) == joinpath(data_root, "bathymetry.nc")  # others unaffected

        # A setup states its own forcing directory, variables and years; only the NorKyst
        # endpoints are defaulted.
        @test_throws UndefKeywordError NorKystConfig(data_root = data_root, years = [2020])
        @test_throws UndefKeywordError NorKystConfig(
            data_root = data_root,
            output_directory = "norkyst",
            years = [2020],
        )
        @test_throws UndefKeywordError NorKystConfig(
            data_root = data_root,
            parameters = ["temperature"],
            years = [2020],
        )

        forcing_config = NorKystConfig(
            data_root = data_root,
            output_directory = "norkyst",
            parameters = ["temperature", "salinity"],
            years = [2020],
        )
        @test forcing_directory(forcing_config) == joinpath(data_root, "norkyst")
        @test forcing_config.parameters == ["temperature", "salinity"]

        # The prepared forcing file and its plot are defaulted, since there is one prepared forcing
        # file per setup. The open edge is not a field here at all: a forcing dataset says nothing
        # about which edge of the domain is open.
        @test forcing_path(forcing_config) == joinpath(data_root, "forcing.nc")
        @test plot_path(forcing_config) == joinpath(data_root, "forcing.png")
        @test !hasfield(typeof(forcing_config), :open_edges)
        @test !hasfield(typeof(forcing_config), :boundaries)

        # Open-boundary data is opt-in and independent of forcing: a `FjordConfig` field of its own,
        # `nothing` unless named, and it is what states the open edge. Its own `boundary_data_*`
        # helpers resolve its file and download directory.
        boundaries_config = test_boundaries_config(; data_root)
        @test boundaries_config isa AbstractBoundaryDataConfig
        @test boundary_data_path(boundaries_config) == joinpath(data_root, "boundaries.nc")
        @test boundary_data_directory(boundaries_config) == joinpath(data_root, "norkyst_hourly")
        @test plot_path(boundaries_config) == joinpath(data_root, "boundaries.png")
        @test open_edges(boundaries_config) == [:south]
        @test open_edges(test_boundaries_config(; data_root, open_edges = :west)) == [:west]
        # No default, so every setup states which edge it opens rather than inheriting a side.
        @test_throws UndefKeywordError NorKystBoundariesConfig(
            data_root = data_root,
            output_directory = "norkyst_hourly",
            parameters = ["temperature"],
            years = [2020],
        )
        # A setup naming none has no open edge at all, which is how a closed domain is stated.
        @test isempty(open_edges(nothing))

        # Where the interpolation runs is a config field, not a command-line flag. `:auto` is the
        # default so one setup runs on a GPU machine and a laptop alike; `:cpu` is honoured
        # regardless of the hardware present, which is what makes this assertion machine-independent.
        @test forcing_config.architecture === :auto
        @test interpolation_architecture(test_forcing_config(architecture = :cpu)) == CPU()
        # An unknown selector is rejected up front rather than falling back to some default.
        @test_throws ArgumentError interpolation_architecture(
            test_forcing_config(architecture = :tpu),
        )

        # ...and an absolute path relocates just that file, as for the bathymetry config.
        relocated = test_forcing_config(
            data_root = data_root,
            output_file = "/shared/forcing.nc",
        )
        @test forcing_path(relocated) == "/shared/forcing.nc"
        @test plot_path(relocated) == joinpath(data_root, "forcing.png")
        @test forcing_directory(test_forcing_config(output_directory = "/nk")) == "/nk"
        @test forcing_monthly_filename(forcing_config, 2020, 3) == "NorKyst-800m_ZDEPTHS_avg_202003.nc"
        @test occursin("thredds.met.no", forcing_config.catalog_url)

        # A setup states its own atmosphere directory and years; only the NORA3 endpoint, the
        # prepared file names and the target resolution are defaulted.
        @test_throws UndefKeywordError NORA3Config(data_root = data_root, years = [2020])
        @test_throws UndefKeywordError NORA3Config(data_root = data_root, output_directory = "nora3")
        @test_throws UndefKeywordError NORA3Config(output_directory = "nora3", years = [2020])

        atmosphere_config = test_atmosphere_config(data_root = data_root)
        @test atmosphere_config isa AbstractAtmosphereConfig
        @test all(isconcretetype, fieldtypes(typeof(atmosphere_config)))
        @test atmosphere_path(atmosphere_config) == joinpath(data_root, "atmosphere.nc")
        @test atmosphere_directory(atmosphere_config) == joinpath(data_root, "nora3")
        @test plot_path(atmosphere_config) == joinpath(data_root, "atmosphere.png")
        @test atmosphere_config.resolution == 0.02
        @test atmosphere_config.padding == 0.1
        @test occursin("thredds.met.no", atmosphere_config.opendap_url)

        # ...and an absolute path relocates just that file, as for the other configs.
        relocated_atmosphere = NORA3Config(
            data_root = data_root,
            output_directory = "nora3",
            output_file = "/shared/atmosphere.nc",
            years = [2020],
        )
        @test atmosphere_path(relocated_atmosphere) == "/shared/atmosphere.nc"
        @test atmosphere_directory(relocated_atmosphere) == joinpath(data_root, "nora3")
        @test atmosphere_directory(
            test_atmosphere_config(data_root = data_root, output_directory = "/nora3"),
        ) == "/nora3"

        grid_config = EvenGrid(
            size = (2, 3, 2),
            halo = (1, 1, 1),
            longitude = (10.0, 12.0),
            latitude = (59.0, 62.0),
            z_faces = [-20.0, -10.0, 0.0],
        )

        boundary_config = test_boundaries_config(; data_root)
        config = FjordConfig(;
            grid_config, bathymetry_config, forcing_config, boundary_config, atmosphere_config,
        )
        grid = domain_grid(config.grid_config, CPU())

        # The open-boundary config is a `FjordConfig` field of its own, independent of the forcing
        # one: either is nameable without the other, which is what lets a setup run with no interior
        # forcing but a data-driven open boundary, or the reverse.
        @test config.boundary_config === boundary_config
        @test isnothing(FjordConfig(; grid_config, bathymetry_config).boundary_config)

        # `native_region!` derives the padded native region from the target grid: 4 cells of
        # 1 degree either side, refined by the default raw_resolution_factor of 4.
        FjordSim.Bathymetry.native_region!(config.bathymetry_config, grid)
        @test config.bathymetry_config.longitude == (6.0, 16.0)
        @test config.bathymetry_config.latitude == (55.0, 66.0)
        @test config.bathymetry_config.size == (40, 44, 1)
        @test dirname(config.bathymetry_config.raw_file) == config.bathymetry_config.raw_directory
        @test occursin("_40x44", config.bathymetry_config.raw_file)

        # Every config is editable in place.
        config.grid_config.size = (4, 6, 2)
        config.bathymetry_config.padding_cells = 0
        config.forcing_config.years = [2020, 2021]
        config.atmosphere_config.years = [2020, 2021]
        @test config.grid_config.size == (4, 6, 2)
        @test config.bathymetry_config.padding_cells == 0
        @test config.forcing_config.years == [2020, 2021]
        @test config.atmosphere_config.years == [2020, 2021]
    end

    @testset "extensibility" begin
        data_root = joinpath(tempdir(), "fjordsim_extensibility_test")

        # FjordConfig accepts the alternative subtypes without any change to its definition.
        config = FjordConfig(
            grid_config = SingleColumnGrid(120.0),
            bathymetry_config = MinimalBathymetry(data_root, "column.nc", "column.png"),
            forcing_config = ConstantForcing(data_root, "column_forcing.nc", "column_forcing.png", 8.0, nothing),
            atmosphere_config = MinimalAtmosphere(
                data_root,
                "column_atmosphere.nc",
                "column_atmosphere.png",
                "column_source",
                0.05,
                0.1,
            ),
            simulation_config = MinimalSimulation(data_root, DateTime(2021, 3, 4, 5), 7200.0),
        )
        rivers = MinimalRivers(data_root, "column_rivers.nc", 3600.0, 10)

        @test config.grid_config isa AbstractGridConfig
        @test config.bathymetry_config isa AbstractBathymetryConfig
        @test config.forcing_config isa AbstractForcingConfig
        @test config.atmosphere_config isa AbstractAtmosphereConfig
        @test config.simulation_config isa AbstractSimulationConfig
        @test rivers isa AbstractRiverConfig

        # Field types are still concrete, so the struct stays type-stable per instantiation.
        @test all(isconcretetype, fieldtypes(typeof(config)))
        @test fieldtype(typeof(config), :grid_config) === SingleColumnGrid
        @test fieldtype(typeof(config), :atmosphere_config) === MinimalAtmosphere
        @test fieldtype(typeof(config), :simulation_config) === MinimalSimulation
        # Unnamed, so its own type parameter resolves to `Nothing` — concrete, like the rest.
        @test fieldtype(typeof(config), :boundary_config) === Nothing

        # ...and the built-in types remain valid, i.e. FjordConfig is genuinely generic.
        @test FjordConfig(
            grid_config = EvenGrid(
                size = (2, 3, 2),
                halo = (1, 1, 1),
                longitude = (10.0, 12.0),
                latitude = (59.0, 62.0),
                z_faces = [-20.0, -10.0, 0.0],
            ),
            bathymetry_config = test_bathymetry_config(data_root = data_root),
            forcing_config = test_forcing_config(data_root = data_root),
        ) isa FjordConfig

        # Path resolution is inherited from the config supertypes — no new methods needed, and
        # `plot_path` serves bathymetry and forcing configs through separate methods.
        @test bathymetry_path(config.bathymetry_config) == joinpath(data_root, "column.nc")
        @test plot_path(config.bathymetry_config) == joinpath(data_root, "column.png")
        @test forcing_path(config.forcing_config) == joinpath(data_root, "column_forcing.nc")
        @test plot_path(config.forcing_config) == joinpath(data_root, "column_forcing.png")
        @test river_forcing_path(rivers) == joinpath(data_root, "column_rivers.nc")
        @test atmosphere_path(config.atmosphere_config) == joinpath(data_root, "column_atmosphere.nc")
        @test atmosphere_directory(config.atmosphere_config) == joinpath(data_root, "column_source")
        @test plot_path(config.atmosphere_config) == joinpath(data_root, "column_atmosphere.png")
        # The simulation config resolves against `results_root`, being the one config that writes, and
        # inherits the run tag, the loop-indexed variant and the coverage window for free. The
        # *filename* is the writer's, so this also proves a built-in writer type composes with a
        # foreign simulation config.
        @test run_tag(config.simulation_config) == FjordSim.Configs.LAUNCH_TAG[]
        column_writer = test_snapshot_writer(output_file = "column_snapshots.nc")
        with_run_tag() do
            @test results_path(column_writer, config.simulation_config) ==
                  joinpath(data_root, "column_snapshots_$PINNED_RUN_TAG.nc")
            @test results_path(column_writer, config.simulation_config, 12) ==
                  joinpath(data_root, "column_snapshots_$(PINNED_RUN_TAG)_loop12.nc")
        end
        # A writer that names no output file says so, rather than failing with `has no field
        # output_file` from inside `tagged_path`.
        @test_throws ArgumentError results_path(test_checkpoint_writer(), config.simulation_config)
        @test coverage_window(config.simulation_config) ==
              (DateTime(2021, 3, 4, 5), DateTime(2021, 3, 4, 7))

        # `river_search_radius` is the river pipeline's one optional hook.
        @test river_search_radius(rivers) == 10

        # A method overloaded on the new grid config is picked up by existing call sites.
        grid = domain_grid(config.grid_config, CPU())
        @test size(grid) == (1, 1, 2)
        @test collect(Oceananigans.Grids.znodes(grid, Face())) == [-120.0, -60.0, 0.0]
        # The grid contract is two hooks, and `SingleColumnGrid` implements only the first, so the
        # second is missing the same discoverable way every other unimplemented hook is. This is
        # what the grid used to have no way of expressing: `build_simulation` reached into
        # `config.grid_config.halo` and built an `ImmersedBoundaryGrid` itself.
        @test_throws MethodError simulation_grid(config.grid_config, "unused.nc", CPU())

        # `regrid_options` is the one optional hook, so a source that configures nothing inherits
        # an empty option set rather than having to implement it.
        @test regrid_options(config.bathymetry_config) === (;)

        # A forcing hook overloaded on the new config is picked up by the generic pipeline...
        @test forcing_variable_names(config.forcing_config) == Dict("temperature" => "T")

        # ...while a required hook that is not implemented fails as a missing method, naming what
        # the new subtype still owes, rather than silently doing nothing.
        @test_throws MethodError bathymetry_dataset(grid, config.bathymetry_config)
        @test_throws MethodError prepare_bathymetry(grid, config.bathymetry_config)
        @test_throws MethodError forcing_time_steps(config.forcing_config)
        @test_throws MethodError forcing_source_grid(config.forcing_config, "unused.nc")
        @test_throws MethodError river_locations(rivers)
        @test_throws MethodError river_series(rivers, [DateTime(2020, 1, 1)])
        @test_throws MethodError atmosphere_time_steps(config.atmosphere_config)
        @test_throws MethodError atmosphere_source_grid(config.atmosphere_config, "unused.nc")
        @test_throws MethodError atmosphere_variable_names(config.atmosphere_config)
        @test_throws MethodError prepare_atmosphere(grid, config.atmosphere_config)
        # The two read-side atmosphere hooks are missing the same way, so a source that prepares a
        # file but cannot be simulated says which methods it still owes.
        @test_throws MethodError prescribed_atmosphere(config.atmosphere_config, CPU())
        @test_throws MethodError prescribed_radiation(config.atmosphere_config, CPU())

        # The supertypes `SimulationConfig` nests own one hook contract each, and each is
        # missing the same discoverable way.
        @test MinimalModel() isa AbstractCoupledSimulationConfig
        @test MinimalFreeSurface() isa AbstractFreeSurfaceConfig
        @test MinimalBoundary() isa AbstractBoundaryConditionConfig
        @test MinimalBoundarySet() isa AbstractBoundaryConditionSetConfig
        @test MinimalWriter() isa AbstractWriterConfig
        @test MinimalCallback() isa AbstractCallbackConfig
        @test MinimalTimeStepping() isa AbstractTimeSteppingConfig

        @test_throws MethodError model_tracers(MinimalModel())
        @test_throws MethodError coupled_simulation(MinimalModel(), grid)
        @test_throws MethodError free_surface(MinimalFreeSurface(), grid)
        @test_throws MethodError boundary_condition_sides(
            MinimalBoundary(), grid, (;), config.boundary_config, (:T,),
        )
        # The set level is its own contract, so a set config that does not materialize its pieces
        # is missing a method rather than falling back to `MergedBoundaryConditions`' merge.
        @test_throws MethodError field_boundary_conditions(
            MinimalBoundarySet(), grid, (;), config.boundary_config, (:T,),
        )
        @test_throws MethodError attach_writer!(
            nothing, MinimalWriter(), config.simulation_config, 1,
        )
        @test_throws MethodError attach_callback!(
            nothing, MinimalCallback(), config.simulation_config,
        )
        @test_throws MethodError attach_time_stepping!(nothing, MinimalTimeStepping())
        @test_throws MethodError initial_time_step(MinimalTimeStepping())

        # Both writer traits default, so a new writer type only overloads the one that is not true
        # of it — a plain writer neither checkpoints nor names a reported file.
        @test !FjordSim.Simulations.checkpoints(MinimalWriter())
        @test FjordSim.Simulations.output_path_trait(MinimalWriter()) isa
              FjordSim.Simulations.NamesNoOutputFile

        # A forcing config carrying no rivers skips the step rather than needing a river dataset.
        @test isnothing(add_rivers(grid, config.forcing_config))

        # ...and a setup naming no atmosphere skips both atmosphere steps the same way, which is why
        # `atmosphere_config` defaults to `nothing`.
        @test isnothing(prepare_atmosphere(grid, nothing))
        @test isnothing(download_atmosphere(grid, nothing))
        @test isnothing(prescribed_atmosphere(nothing, CPU()))
        @test isnothing(prescribed_radiation(nothing, CPU()))
        bare = FjordConfig(
            grid_config = SingleColumnGrid(120.0),
            bathymetry_config = MinimalBathymetry(data_root, "column.nc", "column.png"),
            forcing_config = ConstantForcing(data_root, "f.nc", "f.png", 8.0, nothing),
        )
        @test isnothing(bare.atmosphere_config)
        @test isnothing(download_atmosphere(bare))
        # ...and the same for a setup that is only prepared and never run.
        @test isnothing(bare.simulation_config)
        @test isnothing(build_simulation(bare))
        @test isnothing(run_simulation(bare))

        # `atmosphere_target_axes` is generic over the grid config: it reads the domain through
        # `x_domain`/`y_domain`, so the stub grid works, and both axes must come out uniformly spaced
        # because `compute_faces` infers the spacing from the first difference.
        longitude, latitude = atmosphere_target_axes(grid, config.atmosphere_config)
        @test all(isapprox(longitude[2] - longitude[1]), diff(longitude))
        @test all(isapprox(latitude[2] - latitude[1]), diff(latitude))
        @test first(longitude) <= 10.0 && last(longitude) >= 11.0   # covers the stub grid's domain
        @test first(latitude) <= 59.0 && last(latitude) >= 60.0
    end
end

@testset "Setups" begin
    # Sorted, because Dict order is unspecified — asserted as an invariant rather than as a literal
    # list, so registering a fjord does not break it.
    @test issorted(setup_names())
    @test sort(collect(keys(FjordSim.Setups.SETUPS))) == setup_names()

    # A misspelled name must say which setups exist, not fall through to the file branch.
    @test_throws ArgumentError fjord_config("nordfjorden")
    message = try
        fjord_config("nordfjorden")
    catch exception
        exception.msg
    end
    @test all(occursin(name, message) for name in setup_names())

    # A path that looks like a config file but is not there fails as a file, not as a bad name.
    @test_throws ArgumentError fjord_config(joinpath(tempdir(), "no_such_setup.jl"))

    # An out-of-tree config file is loaded by path, so a fjord need not live in the package.
    mktempdir() do tmp
        external = joinpath(tmp, "hardangerfjorden.jl")
        write(
            external,
            """
            using FjordSim
            FjordConfig(
                grid_config = EvenGrid(
                    size = (8, 8, 2), halo = (3, 3, 3),
                    longitude = (5.5, 5.6), latitude = (60.0, 60.1),
                    z_faces = [-20.0, -10.0, 0.0],
                ),
                bathymetry_config = DybdedataConfig(
                    data_root = $(repr(tmp)), output_file = "b.nc", plot_file = "b.png",
                ),
                forcing_config = NorKystConfig(
                    data_root = $(repr(tmp)), output_directory = "norkyst",
                    output_file = "f.nc", plot_file = "f.png",
                    parameters = ["temperature"], years = [2020],
                ),
            )
            """,
        )

        external_config = fjord_config(external)
        @test external_config isa FjordConfig
        @test external_config.grid_config.size == (8, 8, 2)
        @test bathymetry_path(external_config.bathymetry_config) == joinpath(tmp, "b.nc")
        @test isnothing(external_config.atmosphere_config)  # defaults apply to a file config too

        # A file whose last expression is not a FjordConfig has to say so, rather than failing
        # later inside a pipeline with an unrelated error.
        not_a_config = joinpath(tmp, "not_a_config.jl")
        write(not_a_config, "42\n")
        @test_throws ArgumentError fjord_config(not_a_config)
    end

    for name in setup_names()
        config = fjord_config(name)
        data_root = joinpath(homedir(), "FjordSim_data", name)

        @test config isa FjordConfig
        # Every field is concrete or parameterized. `isconcretetype(typeof(config))` would be
        # vacuous — `typeof` never returns an abstract type.
        @test all(isconcretetype, fieldtypes(typeof(config)))
        @test config.bathymetry_config.data_root == data_root
        @test config.forcing_config.data_root == data_root
        @test startswith(forcing_directory(config.forcing_config), data_root)

        # Built inside the function, so the scratch path `__init__` fills in is already there. A
        # config built at precompile time would carry an empty `raw_directory` instead.
        @test isdir(config.bathymetry_config.raw_directory)

        grid_config = config.grid_config
        @test length(grid_config.z_faces) == grid_config.size[3] + 1
        @test issorted(grid_config.z_faces)  # bottom to top, as ImmersedBoundaryGrid expects
        @test last(grid_config.z_faces) == 0.0
        @test grid_config.longitude[1] < grid_config.longitude[2]
        @test grid_config.latitude[1] < grid_config.latitude[2]
    end

    # Each call returns a fresh config: the structs are mutable and `native_region!` edits the
    # bathymetry config, so a cached `const` instance would leak state between steps.
    first_call = oslofjorden()
    second_call = oslofjorden()
    @test first_call !== second_call
    @test first_call.bathymetry_config !== second_call.bathymetry_config
    first_call.grid_config.size = (1, 1, 1)
    @test second_call.grid_config.size != (1, 1, 1)

    # Conventions every registered setup follows, asserted over all of them. Which dataset a setup
    # picks, and how long it runs for, are its author's to change — so nothing here names a value,
    # only the conventions a setup getting them wrong would break.
    for name in setup_names()
        config = fjord_config(name)

        # A nested config downloads under the setup's own root rather than sharing from elsewhere.
        rivers = config.forcing_config.rivers
        if !isnothing(rivers)
            @test rivers.data_root == config.forcing_config.data_root
        end

        if !isnothing(config.atmosphere_config)
            @test startswith(
                atmosphere_path(config.atmosphere_config),
                config.atmosphere_config.data_root,
            )
        end

        # Results are rooted away from the input data, the simulation config being the only one that
        # writes rather than reads.
        if !isnothing(config.simulation_config)
            @test !startswith(
                config.simulation_config.results_root,
                config.bathymetry_config.data_root,
            )
        end
    end
end

@testset "CLI" begin
    @testset "argument parsing" begin
        parse_arguments = FjordSim.CLI.parse_arguments

        # `--config VALUE` and `--config=VALUE` are the same, in either order around the subcommand.
        @test parse_arguments(["prepare_forcing", "--config", "oslofjorden"]) ==
              (subcommand = "prepare_forcing", config = "oslofjorden", help = false)
        @test parse_arguments(["prepare_forcing", "--config=oslofjorden"]) ==
              parse_arguments(["prepare_forcing", "--config", "oslofjorden"])
        @test parse_arguments(["--config", "oslofjorden", "prepare_forcing"]) ==
              parse_arguments(["prepare_forcing", "--config", "oslofjorden"])

        @test_throws ArgumentError parse_arguments(String[])                              # no subcommand
        @test_throws ArgumentError parse_arguments(["prepare_forcing"])                   # no --config
        @test_throws ArgumentError parse_arguments(["prepare_forcing", "--config"])       # no value
        @test_throws ArgumentError parse_arguments(["prepare_forcing", "--config="])      # empty value
        @test_throws ArgumentError parse_arguments(["prepare_forcing", "--gpu", "--config", "oslofjorden"])
        @test_throws ArgumentError parse_arguments(["prepare_forcing", "add_rivers", "--config", "oslofjorden"])
        @test_throws ArgumentError FjordSim.CLI.subcommand_driver("frobnicate")

        # `--help` returns a sentinel instead of exiting, which is what makes it testable at all.
        @test parse_arguments(["-h"]).help
        @test parse_arguments(["--help"]).help
        @test parse_arguments(["prepare_forcing", "--config", "oslofjorden", "--help"]).help
    end

    @testset "subcommand table" begin
        subcommands = FjordSim.CLI.subcommand_names()
        @test subcommands == [
            "prepare_bathymetry",
            "download_forcing",
            "prepare_forcing",
            "add_rivers",
            "download_boundaries",
            "prepare_boundaries",
            "download_atmosphere",
            "prepare_atmosphere",
            "run_simulation",
        ]
        for name in subcommands
            @test occursin(name, FjordSim.CLI.USAGE)  # the usage text cannot drift from the table
        end
        for name in setup_names()
            @test occursin(name, FjordSim.CLI.USAGE)
        end

        # Every subcommand resolves to a driver that takes a whole setup. This is the assertion that
        # catches a driver method that was never defined, without running any of them.
        config = drammensfjorden()
        for name in subcommands
            @test hasmethod(FjordSim.CLI.subcommand_driver(name), Tuple{typeof(config)})
        end
    end

    @testset "main exit codes" begin
        @test isdefined(FjordSim, :main)
        @test hasmethod(FjordSim.main, Tuple{Vector{String}})
        # Exporting `main` would make Julia run the CLI after the body of any script that says
        # `using FjordSim` — including this test file.
        @test !(:main in names(FjordSim))
        @test FjordSim.main(["--help"]) === 0            # must be an exit code, not a driver's result
        @test FjordSim.main(["frobnicate", "--config", "oslofjorden"]) == 2
        @test FjordSim.main(["prepare_forcing", "--config", "nordfjorden"]) == 2
        @test FjordSim.main(String[]) == 2
    end

    @testset "opt-out no-ops" begin
        # A step the setup opts out of is a no-op through the FjordConfig arity too, matching
        # `add_rivers(grid, config, ::Nothing)` and `prepare_atmosphere(grid, ::Nothing)`. Rooted in a
        # temporary directory so nothing can touch the real ~/FjordSim_data.
        mktempdir() do tmp
            # `drammensfjorden()` itself now names rivers, an atmosphere and a simulation config, so a
            # config naming none of them is assembled from its grid and bathymetry plus a rivers-free
            # forcing config — `atmosphere_config` and `simulation_config` default to `nothing` when
            # omitted, and `rivers` is a type parameter of NorKystConfig, so it has to be left unnamed at
            # construction rather than unset afterwards.
            base = drammensfjorden()
            config = FjordConfig(
                grid_config = base.grid_config,
                bathymetry_config = base.bathymetry_config,
                forcing_config = test_forcing_config(data_root = tmp),
            )
            config.bathymetry_config.data_root = tmp

            @test isnothing(add_rivers(config))
            @test isnothing(prepare_atmosphere(config))
            @test isnothing(download_atmosphere(config))
            # Both return on the `nothing` simulation config, before any prerequisite is checked, so
            # neither touches the empty temporary root.
            @test isnothing(build_simulation(config))
            @test isnothing(run_simulation(config))
            @test FjordSim.main(["add_rivers", "--config", joinpath(tmp, "unused.jl")]) == 2
        end
    end

    @testset "tee_output" begin
        # The tee returns the closure's value — that is how `main` gets its exit code back out — and
        # captures both streams into one file in the order they were written.
        mktempdir() do tmp
            log_file = joinpath(tmp, "tee.log")
            @test FjordSim.CLI.tee_output(log_file) do
                println("out-line")
                println(stderr, "err-line")
                42
            end == 42

            logged = read(log_file, String)
            @test occursin("out-line", logged)
            @test occursin("err-line", logged)
        end
    end

    @testset "the run log" begin
        # Only `run_simulation` is logged, so a failing prepare step is exit code 1 and leaves no file
        # behind at all. An out-of-tree config rooted in `tmp` rather than the registered name, so the
        # assertion does not depend on whether ~/FjordSim_data/drammensfjorden happens to exist.
        mktempdir() do tmp
            config_file = joinpath(tmp, "config.jl")
            write(
                config_file,
                """
                using FjordSim
                config = drammensfjorden()
                config.bathymetry_config.data_root = raw"$tmp"
                config.forcing_config.data_root = raw"$tmp"
                config
                """,
            )

            cd(tmp) do
                @test FjordSim.main(["prepare_forcing", "--config", config_file]) == 1
                @test isempty(filter(name -> startswith(name, "fjordsim"), readdir(tmp)))
            end
        end

        # The log goes beside the output it describes, tagged with the launch instant so runs stay
        # distinguishable, and takes the simulation config because that is where both come from.
        oslo = oslofjorden()
        @test with_run_tag(() -> FjordSim.CLI.log_path(oslo.simulation_config)) ==
              joinpath(oslo.simulation_config.results_root, "fjordsim_$PINNED_RUN_TAG.log")

        # End to end: a failing `run_simulation` writes its log into the setup's results root and not
        # into the working directory, creating the directory if it does not exist yet. `build_simulation`
        # stops at the bathymetry prerequisite, which is as far as a run gets without data or a GPU.
        mktempdir() do tmp
            results = joinpath(tmp, "results")
            config_file = joinpath(tmp, "config.jl")
            write(
                config_file,
                """
                using FjordSim
                config = oslofjorden()
                config.bathymetry_config.data_root = raw"$tmp"
                config.forcing_config.data_root = raw"$tmp"
                config.simulation_config.results_root = raw"$results"
                config
                """,
            )

            cd(tmp) do
                @test !isdir(results)
                @test FjordSim.main(["run_simulation", "--config", config_file]) == 1
                @test isempty(filter(name -> startswith(name, "fjordsim"), readdir(tmp)))

                log_file = joinpath(results, "fjordsim_$(FjordSim.Configs.LAUNCH_TAG[]).log")
                logged = read(log_file, String)
                @test occursin("Processed bathymetry", logged)
                @test occursin("Stacktrace", logged)
            end
        end

        # A setup naming no simulation config has nothing to log and no results root to log into, so
        # `run_simulation` stays a no-op that writes nothing. `drammensfjorden()` itself now names one,
        # so this config is assembled from its grid/bathymetry/forcing instead, leaving
        # `simulation_config` at its default of `nothing`.
        mktempdir() do tmp
            config_file = joinpath(tmp, "config.jl")
            write(
                config_file,
                """
                using FjordSim
                base = drammensfjorden()
                config = FjordConfig(
                    grid_config = base.grid_config,
                    bathymetry_config = base.bathymetry_config,
                    forcing_config = base.forcing_config,
                )
                config.bathymetry_config.data_root = raw"$tmp"
                config.forcing_config.data_root = raw"$tmp"
                config
                """,
            )

            cd(tmp) do
                @test FjordSim.main(["run_simulation", "--config", config_file]) == 0
                @test isempty(filter(name -> startswith(name, "fjordsim"), readdir(tmp)))
            end
        end
    end

    @testset "abbreviated stacktraces" begin
        # Stacktrace frames are reported with their type parameters abbreviated, the way the REPL does
        # it: a bare `showerror` spells out every one of them, which is what makes a trace through the
        # coupled model unreadable in the terminal and in the log alike.
        # Deep enough that the frame is past `STACKTRACE_WIDTH` spelled out: `type_depth_limit` only
        # abbreviates what does not fit.
        inner = Nested{Nested{Int,Float64},Nested{Float32,Bool}}
        exception, backtrace = try
            nested_failure(Nested{Nested{inner,inner},Nested{inner,inner}}())
        catch caught
            caught, catch_backtrace()
        end

        compact = sprint(io -> FjordSim.CLI.show_compact_error(io, exception, backtrace))
        verbose = sprint(io -> showerror(io, exception, backtrace))

        @test occursin("nested boom", compact)
        @test occursin("{…}", compact)
        @test occursin("abbreviated", compact)
        @test !occursin("{…}", verbose)
        @test length(compact) < length(verbose)
    end
end

@testset "Bathymetry" begin
    @testset "writer" begin
        mktempdir() do tmp
            architecture = CPU()
            z_faces = [-20.0, -10.0, 0.0]
            grid = LatitudeLongitudeGrid(
                architecture;
                size = (2, 3, 2),
                halo = (1, 1, 1),
                longitude = (10.0, 12.0),
                latitude = (59.0, 62.0),
                z = z_faces,
            )

            bottom_height = Field{Center, Center, Nothing}(grid)
            set!(bottom_height, [-15.0 -16.0 -17.0; -18.0 -19.0 -20.0])

            bathymetry_file = joinpath(tmp, "bathymetry_written.nc")
            @test_nowarn write_bathymetry_file(bathymetry_file, grid, bottom_height)

            ds = NCDataset(bathymetry_file)
            @test ds["lon"][:] == [10.5, 11.5]
            @test ds["lat"][:] == [59.5, 60.5, 61.5]
            @test ds["h"][:, :] == Float32[-15.0 -16.0 -17.0; -18.0 -19.0 -20.0]
            # z_faces must survive the write unchanged: `vertical_faces` previously offset its
            # slice of the OffsetArray `grid.z.cᵃᵃᶠ` by the halo, dropping the deepest face and
            # appending one above the surface.
            @test ds["z_faces"][:] == z_faces
            close(ds)

            # Full round-trip: grid -> file -> grid must preserve the vertical extent.
            reloaded = ImmersedBoundaryGrid(bathymetry_file, architecture, (1, 1, 1))
            @test collect(Oceananigans.Grids.znodes(reloaded.underlying_grid, Face())) == z_faces

            @test FjordSim.Bathymetry.vertical_faces(grid) == z_faces
            # Independent of the halo size, which is what the old indexing got wrong.
            for halo_size in (1, 2, 3)
                deep_faces = [-450.0, -200.0, -50.0, 0.0]
                deep_grid = LatitudeLongitudeGrid(
                    architecture;
                    size = (4, 4, 3),
                    halo = (halo_size, halo_size, halo_size),
                    longitude = (10.0, 11.0),
                    latitude = (59.0, 60.0),
                    z = deep_faces,
                )
                @test FjordSim.Bathymetry.vertical_faces(deep_grid) == deep_faces
            end

            @test FjordSim.Bathymetry.contour_point_indices(10, 3) == [0, 3, 6, 9]
            @test FjordSim.Bathymetry.contour_point_indices(5, 3) == [0, 3, 4]
        end
    end

    @testset "gap filling" begin
        fill_diagonal_pairs = FjordSim.Bathymetry.fill_diagonal_pairs
        fill_secondary_diagonal_pairs = FjordSim.Bathymetry.fill_secondary_diagonal_pairs
        remove_isolated_sea_cells = FjordSim.Bathymetry.remove_isolated_sea_cells
        fill_isolated_land_cells = FjordSim.Bathymetry.fill_isolated_land_cells

        # top-left/bottom-right sea, top-right/bottom-left land -> fill top-right
        h = [-5.0 1.0; 1.0 -3.0]
        filled = fill_diagonal_pairs(h)
        @test filled[1, 2] == -4.0
        @test filled[1, 1] == -5.0
        @test filled[2, 1] == 1.0
        @test filled[2, 2] == -3.0
        @test fill_secondary_diagonal_pairs(h) == h  # opposite diagonal pattern doesn't match

        # top-right/bottom-left sea, top-left/bottom-right land -> fill top-left
        h2 = [1.0 -5.0; -3.0 1.0]
        filled2 = fill_secondary_diagonal_pairs(h2)
        @test filled2[1, 1] == -4.0
        @test filled2[1, 2] == -5.0
        @test filled2[2, 1] == -3.0
        @test filled2[2, 2] == 1.0
        @test fill_diagonal_pairs(h2) == h2  # opposite diagonal pattern doesn't match

        # interior sea cell surrounded by land on all 4 sides -> becomes land (0)
        h3 = [0.0 0.0 0.0; 0.0 -5.0 0.0; 0.0 0.0 0.0]
        replaced = remove_isolated_sea_cells(h3)
        @test replaced[2, 2] == 0.0
        @test replaced[1, :] == h3[1, :]  # boundary row untouched

        # a boundary sea cell is left alone even if it would otherwise qualify
        h4 = [0.0 -5.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        @test remove_isolated_sea_cells(h4) == h4

        # interior land cell surrounded by sea on all 4 sides -> filled with neighbor mean
        h5 = [-1.0 -2.0 -1.0; -6.0 0.0 -8.0; -1.0 -4.0 -1.0]
        filled5 = fill_isolated_land_cells(h5)
        @test filled5[2, 2] == -5.0  # mean(-6, -2, -8, -4)
        @test filled5[1, :] == h5[1, :]  # boundary row untouched

        # a boundary land cell is left alone even if it would otherwise qualify
        h6 = [0.0 -1.0 0.0; -1.0 -1.0 -1.0; 0.0 -1.0 0.0]
        @test fill_isolated_land_cells(h6) == h6

        remove_narrow_passages = FjordSim.Bathymetry.remove_narrow_passages
        wet_component = FjordSim.Bathymetry.wet_component

        # The flood fill honours `blocked`, which is what lets one candidate be tested without
        # mutating the field first.
        channel = [-1.0 -1.0 -1.0 -1.0 -1.0]
        @test count(wet_component(channel, (1, 1), (0, 0))) == 5
        @test count(wet_component(channel, (1, 1), (1, 3))) == 2

        # A basin split by a peninsula along j = 4, cut through at (4, 4) by a one-cell canal — the
        # geometry of the real Oslofjord defect. The canal closes a loop, since the water either side
        # of it is joined around both ends of the peninsula, so it goes.
        #
        # Two details of the fixture are load-bearing. The water is at least two cells thick
        # everywhere else, because in a one-cell-wide ring *every* cell is a one-cell-wide passage.
        # And the peninsula is two cells thick, because `fill_isolated_land_cells` would flood a
        # one-cell-thick one before this stage ever saw it — which is plausibly where the real canals
        # came from.
        loop = [
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -2.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
        ]
        opened = remove_narrow_passages(loop)
        @test opened[4, 4] == 0.0
        @test count(<(0), opened) == count(<(0), loop) - 1  # only the canal closed
        @test remove_narrow_passages(opened) == opened  # nothing left to close

        # The same one-cell canal when it is the *only* link between the two sides: kept, because
        # closing it would delete that water from the domain. This is what the connectivity test buys.
        pocket = [
            0.0  0.0  0.0  0.0  0.0  0.0 0.0
            0.0 -9.0 -9.0  0.0 -9.0 -9.0 0.0
            0.0 -9.0 -9.0 -2.0 -9.0 -9.0 0.0
            0.0 -9.0 -9.0  0.0 -9.0 -9.0 0.0
            0.0  0.0  0.0  0.0  0.0  0.0 0.0
        ]
        @test remove_narrow_passages(pocket) == pocket

        # Two canals through one isthmus, each redundant only while the other is open. Testing
        # candidates one at a time against the partially closed field closes exactly one; testing
        # them all at once against the input would close both and sever the two sides.
        pair = [
            0.0  0.0  0.0  0.0 0.0
            0.0 -9.0  0.0 -9.0 0.0
            0.0 -9.0 -2.0 -9.0 0.0
            0.0 -9.0  0.0 -9.0 0.0
            0.0 -9.0 -2.0 -9.0 0.0
            0.0 -9.0  0.0 -9.0 0.0
            0.0  0.0  0.0  0.0 0.0
        ]
        unlooped = remove_narrow_passages(pair)
        @test count(<(0), unlooped) == count(<(0), pair) - 1
        @test count(wet_component(unlooped, (2, 2), (0, 0))) == count(<(0), unlooped)  # still one basin

        # A two-cell-wide canal through the same isthmus is not narrow, so it is left alone: the
        # stage targets width, and two cells is the narrowest a resolved flow can occupy.
        wide = [
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -2.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -2.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
        ]
        @test remove_narrow_passages(wide) == wide

        fill_small_islands = FjordSim.Bathymetry.fill_small_islands

        # The dual of the canal: a three-cell, one-cell-wide *island* in mid-channel, the geometry
        # that carried the Oslofjord domain velocity maximum. Flow splits around it and closes a loop
        # just as a one-cell passage does, and one cell resolves nothing about that flow.
        # `fill_isolated_land_cells` cannot take it — every cell of it has a land neighbour, so none
        # has the four wet ones that stage needs.
        island = [
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0  0.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
        ]
        @test fill_isolated_land_cells(island) == island
        @test fill_small_islands(island; max_cells = 3) == fill(-9.0, 7, 7)

        # The threshold is a hard one, and `0` — the default — disables the stage outright.
        @test fill_small_islands(island; max_cells = 2) == island
        @test fill_small_islands(island; max_cells = 0) == island

        # The patch takes one depth, the mean over the sea cells around it. For a single cell that is
        # exactly what `fill_isolated_land_cells` computes, which is the sense in which this stage
        # generalizes it from one cell to `max_cells`.
        uneven = [
            -2.0 -4.0 -2.0
            -6.0  0.0 -8.0
            -2.0 -4.0 -2.0
        ]
        @test fill_small_islands(uneven; max_cells = 1)[2, 2] == -5.5  # mean(-4, -6, -8, -4)
        @test fill_small_islands(uneven; max_cells = 1) == fill_isolated_land_cells(uneven)

        # Size is measured per component: the three-cell island goes and the six-cell reef beside it
        # stays. And the result needs no cleanup pass — every neighbour of a flooded patch is sea, so
        # nothing new is isolated either way.
        reef = [
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
            -9.0  0.0 -9.0 -9.0  0.0  0.0 -9.0
            -9.0  0.0 -9.0 -9.0  0.0  0.0 -9.0
            -9.0  0.0 -9.0 -9.0  0.0  0.0 -9.0
            -9.0 -9.0 -9.0 -9.0 -9.0 -9.0 -9.0
        ]
        patched = fill_small_islands(reef; max_cells = 3)
        @test all(<(0), patched[:, 2])                # the three-cell island is flooded
        @test patched[2:4, 5:6] == reef[2:4, 5:6]     # the six-cell reef is kept whole
        @test count(>=(0), patched) == 6
        @test fill_isolated_land_cells(remove_isolated_sea_cells(patched)) == patched

        # A patch of the very same size on the domain edge is kept whatever the threshold: land that
        # continues outside the domain can have only a few cells inside it, and flooding them would
        # open the domain into water that is not there.
        shore = [
             0.0  0.0 -9.0 -9.0 -9.0
             0.0 -9.0 -9.0 -9.0 -9.0
            -9.0 -9.0 -9.0 -9.0 -9.0
        ]
        @test fill_small_islands(shore; max_cells = 10) == shore

        snap_partial_bottom_cells = FjordSim.Bathymetry.snap_partial_bottom_cells

        # `PartialCellBottom` will not make a bottom cell thinner than `minimum_cell_fraction` of its
        # layer; it pushes the bottom *down* until the cell is that thick, so a sounding just below a
        # face gives a cell whose floor is somewhere the sounding never was, holding a fraction of the
        # water its neighbours hold. This stage moves the bottom up to the face instead, ending the
        # column a layer higher. These faces are the 25 m layer that produced the Oslofjord runaway.
        faces = [-100.0, -75.0, -50.0, -25.0, -10.0, 0.0]
        soundings = [-51.31 -54.49 -57.72 -44.05]
        snapped = snap_partial_bottom_cells(soundings, faces; minimum_cell_fraction = 0.2)
        @test snapped[1] == -50.0   # fraction 0.052, raised to the face
        @test snapped[2] == -50.0   # fraction 0.180, raised
        @test snapped[3] == -57.72  # fraction 0.309, left alone
        @test snapped[4] == -44.05  # fraction 0.762, left alone

        # Land is never touched, and `0` disables the stage.
        mixed = [-51.31 0.0 5.0]
        @test snap_partial_bottom_cells(mixed, faces; minimum_cell_fraction = 0.2)[2:3] == [0.0, 5.0]
        @test snap_partial_bottom_cells(soundings, faces; minimum_cell_fraction = 0.0) == soundings

        # Two guards, so the stage can never dry a cell out or undercut the depth floor. A sliver in
        # the topmost layer would snap to z = 0, and one in the layer above `minimum_depth` would snap
        # through it; both keep their sliver instead and are floored by `PartialCellBottom` as before.
        guarded = [-1.5 -26.0]
        @test snap_partial_bottom_cells(guarded, faces; minimum_cell_fraction = 0.2) == [-1.5 -25.0]
        @test snap_partial_bottom_cells(
            guarded, faces; minimum_cell_fraction = 0.2, minimum_depth = 30.0,
        ) == guarded

        fill_shallow_spikes = FjordSim.Bathymetry.fill_shallow_spikes
        limit_bottom_slope = FjordSim.Bathymetry.limit_bottom_slope

        # A cell far shallower than its neighbours is lifted to their median depth. This is the
        # artifact `minimum_depth` leaves behind: a sliver floored to a constant depth, still far
        # shallower than the water around it.
        spike = [-10.0 -10.0 -10.0; -10.0 -1.0 -10.0; -10.0 -10.0 -10.0]
        despiked = fill_shallow_spikes(spike; ratio = 0.5)
        @test despiked[2, 2] == -10.0
        @test despiked[1, :] == spike[1, :]  # everything else untouched

        # 60% of the neighbour median is above a 0.5 threshold, so it is left alone: the filter must
        # not shave genuine topography, only cells that disagree with their surroundings.
        mild = [-10.0 -10.0 -10.0; -10.0 -6.0 -10.0; -10.0 -10.0 -10.0]
        @test fill_shallow_spikes(mild; ratio = 0.5) == mild

        # A uniformly shallow region is not a spike — its neighbours are shallow too.
        shelf = [-2.0 -2.0 -2.0; -2.0 -1.5 -2.0; -2.0 -2.0 -2.0]
        @test fill_shallow_spikes(shelf; ratio = 0.5) == shelf

        # Land is never touched, and a cell with fewer than `min_neighbours` sea neighbours is skipped
        # so a single-cell inlet is not swallowed by the one neighbour it has.
        inlet = [0.0 0.0 0.0; 0.0 -1.0 -10.0; 0.0 0.0 0.0]
        @test fill_shallow_spikes(inlet; ratio = 0.5) == inlet
        @test fill_shallow_spikes(inlet; ratio = 0.5, min_neighbours = 1)[2, 2] == -10.0

        # Slope limiting moves each offending pair symmetrically, so the pair's depth sum — and hence
        # the domain's water volume — is unchanged, and the resulting r sits exactly at the limit.
        steep = reshape([-1.0, -9.0], 2, 1)
        limited = limit_bottom_slope(steep; max_slope_factor = 0.5)
        d1, d2 = -limited[1, 1], -limited[2, 1]
        @test d1 + d2 ≈ 10.0
        @test abs(d1 - d2) / (d1 + d2) ≈ 0.5
        @test d1 > 1.0 && d2 < 9.0  # shallow deepened, deep shoaled

        # Already within the limit: returned untouched.
        gentle = reshape([-9.0, -10.0], 2, 1)
        @test limit_bottom_slope(gentle; max_slope_factor = 0.5) == gentle

        # The floor wins over the slope limit — flattening a slope must not undo `minimum_depth`.
        @test all(<=(-2.0), limit_bottom_slope(steep; max_slope_factor = 0.5, minimum_depth = 2.0))

        # Land stays land, and the whole field converges below the limit within the pass cap.
        mixed = [-1.0 -40.0 0.0; -30.0 -2.0 -50.0; 0.0 -60.0 -1.0]
        smoothed = limit_bottom_slope(mixed; max_slope_factor = 0.4)
        @test smoothed[1, 3] == 0.0 && smoothed[3, 1] == 0.0
        for i in 1:3, j = 1:3
            smoothed[i, j] < 0 || continue
            for (di, dj) in ((1, 0), (0, 1))
                ii, jj = i + di, j + dj
                (ii <= 3 && jj <= 3) || continue
                smoothed[ii, jj] < 0 || continue
                here, there = -smoothed[i, j], -smoothed[ii, jj]
                @test abs(here - there) / (here + there) <= 0.4 + 1e-9
            end
        end

        # All three stages are opt-in: the default keeps the topological cleanup and nothing else,
        # which is what leaves every other source's bathymetry unchanged. Asserted on configs built
        # here rather than on the built-in setups', whose numbers are a per-fjord choice.
        smoothing_options = FjordSim.Bathymetry.smoothing_options
        @test smoothing_options(test_bathymetry_config()) == (
            max_island_cells = 0,
            close_narrow_passages = false,
            spike_ratio = 0.0,
            minimum_cell_fraction = 0.0,
            max_slope_factor = 0.0,
            minimum_depth = 0.0,
        )
        @test smoothing_options(test_bathymetry_config(
            max_island_cells = 6,
            close_narrow_passages = true,
            spike_ratio = 0.5,
            minimum_cell_fraction = 0.2,
            max_slope_factor = 0.4,
            minimum_depth = 2.0,
        )) == (
            max_island_cells = 6,
            close_narrow_passages = true,
            spike_ratio = 0.5,
            minimum_cell_fraction = 0.2,
            max_slope_factor = 0.4,
            minimum_depth = 2.0,
        )

        # smooth_bathymetry_gaps! round-trips a Field through the same pipeline
        architecture = CPU()
        grid = LatitudeLongitudeGrid(
            architecture;
            size = (3, 3, 1),
            halo = (1, 1, 1),
            longitude = (10.0, 11.0),
            latitude = (59.0, 60.0),
            z = [-10.0, 0.0],
        )
        raw = Float32[0.0 0.0 0.0; 0.0 -5.0 0.0; 0.0 0.0 0.0]
        bottom_height = Field{Center, Center, Nothing}(grid)
        set!(bottom_height, raw)

        FjordSim.Bathymetry.smooth_bathymetry_gaps!(bottom_height)

        expected =
            FjordSim.Bathymetry.fill_secondary_diagonal_pairs(FjordSim.Bathymetry.fill_diagonal_pairs(copy(raw)))
        for _ in 1:FjordSim.Bathymetry.BATHYMETRY_GAP_FILL_PASSES
            expected = fill_isolated_land_cells(remove_isolated_sea_cells(expected))
        end

        @test Array(interior(bottom_height, :, :, 1)) == expected

        # The passage stage is wired into smooth_bathymetry_gaps! and its effect survives the two
        # stages after it. The same canal shows why despiking cannot substitute: with the stage off,
        # `spike_ratio` still lifts the 2 m canal to its neighbours' depth — fixing the *depth* and
        # leaving the spurious connection wide open, which is the failure this stage exists for.
        canal_grid = LatitudeLongitudeGrid(
            CPU();
            size = (7, 7, 1),
            halo = (1, 1, 1),
            longitude = (10.0, 11.0),
            latitude = (59.0, 60.0),
            z = [-10.0, 0.0],
        )
        canal_field = Field{Center, Center, Nothing}(canal_grid)

        set!(canal_field, Float32.(loop))
        FjordSim.Bathymetry.smooth_bathymetry_gaps!(
            canal_field;
            close_narrow_passages = true,
            spike_ratio = 0.5,
            max_slope_factor = 0.5,
        )
        closed_canal = Array(interior(canal_field, :, :, 1))
        @test closed_canal[4, 4] >= 0                              # canal closed
        @test remove_narrow_passages(closed_canal) == closed_canal  # and nothing reopened it
        @test closed_canal[4, 3] == -9.0f0                          # neighbours untouched

        set!(canal_field, Float32.(loop))
        FjordSim.Bathymetry.smooth_bathymetry_gaps!(
            canal_field;
            close_narrow_passages = false,
            spike_ratio = 0.5,
            max_slope_factor = 0.5,
        )
        open_canal = Array(interior(canal_field, :, :, 1))
        @test open_canal[4, 4] == -9.0f0  # despiked to its neighbours' depth, still wet

        # The island stage is wired in and is opt-in the same way. Left at its default the whole
        # pipeline leaves the island standing, which is what keeps every source that configures
        # nothing unchanged: the topological cleanup cannot take it, because
        # `fill_isolated_land_cells` only ever fires on a single cell.
        island_field = Field{Center, Center, Nothing}(canal_grid)

        set!(island_field, Float32.(island))
        FjordSim.Bathymetry.smooth_bathymetry_gaps!(island_field)
        @test Array(interior(island_field, :, :, 1)) == Float32.(island)

        set!(island_field, Float32.(island))
        FjordSim.Bathymetry.smooth_bathymetry_gaps!(island_field; max_island_cells = 3)
        @test all(<(0), Array(interior(island_field, :, :, 1)))
    end

    @testset "point sampling" begin
        # Synthetic dybdepunkt/dybdekurve-like layers in EPSG:25833, built on a Memory
        # dataset so this test needs no real Geonorge FileGDB.
        source_points = [(10000.0, 6_600_000.0), (10500.0, 6_600_500.0), (11000.0, 6_601_000.0)]
        source_depths = [5.0, -12.5, 30.0]

        contour_vertices =
            [(20000.0, 6_610_000.0), (20100.0, 6_610_100.0), (20200.0, 6_610_200.0), (20300.0, 6_610_300.0)]
        contour_depth = 42.0
        contour_stride = 2

        # `collect_depth_layer_coordinates!` reads the spatial filter and the contour stride off
        # the config, so set the derived `filter_bounds` directly instead of via `native_region!`.
        bathymetry_config = test_bathymetry_config(
            contour_stride = contour_stride,
            filter_bounds = (0.0, 6_500_000.0, 100_000.0, 6_700_000.0),
        )
        samples = FjordSim.Bathymetry.DepthSamples()
        xs, ys, bottom_heights = samples.xs, samples.ys, samples.bottom_heights

        ArchGDAL.importEPSG(25833; order = :trad) do source_srs
            ArchGDAL.create(ArchGDAL.getdriver("Memory")) do dataset
                ArchGDAL.createlayer(
                    name = "dybdepunkt",
                    dataset = dataset,
                    geom = ArchGDAL.wkbPoint,
                    spatialref = source_srs,
                ) do layer
                    ArchGDAL.addfielddefn!(layer, "dybde", ArchGDAL.OFTReal)
                    depth_index = ArchGDAL.findfieldindex(layer, "dybde", false)
                    for ((x, y), depth) in zip(source_points, source_depths)
                        ArchGDAL.createfeature(layer) do feature
                            ArchGDAL.setfield!(feature, depth_index, depth)
                            ArchGDAL.setgeom!(feature, 0, ArchGDAL.createpoint(x, y))
                            return nothing
                        end
                    end
                end

                ArchGDAL.createlayer(
                    name = "dybdekurve",
                    dataset = dataset,
                    geom = ArchGDAL.wkbMultiLineString,
                    spatialref = source_srs,
                ) do layer
                    ArchGDAL.addfielddefn!(layer, "dybde", ArchGDAL.OFTReal)
                    depth_index = ArchGDAL.findfieldindex(layer, "dybde", false)
                    wkt_points = join(("$x $y" for (x, y) in contour_vertices), ", ")
                    ArchGDAL.createfeature(layer) do feature
                        ArchGDAL.setfield!(feature, depth_index, contour_depth)
                        ArchGDAL.setgeom!(feature, 0, ArchGDAL.fromWKT("MULTILINESTRING (($wkt_points))"))
                        return nothing
                    end
                end

                FjordSim.Bathymetry.collect_depth_layer_coordinates!(
                    samples,
                    dataset,
                    "dybdepunkt",
                    bathymetry_config;
                    geometry = :point,
                )
                FjordSim.Bathymetry.collect_depth_layer_coordinates!(
                    samples,
                    dataset,
                    "dybdekurve",
                    bathymetry_config;
                    geometry = :line,
                )
            end
        end

        expected_contour_indices = FjordSim.Bathymetry.contour_point_indices(length(contour_vertices), contour_stride)
        n_points = length(source_points)

        @test length(xs) == n_points + length(expected_contour_indices)
        @test xs[1:n_points] == first.(source_points)
        @test ys[1:n_points] == last.(source_points)
        @test bottom_heights[1:n_points] == -abs.(source_depths)
        @test xs[n_points+1:end] == first.(contour_vertices[expected_contour_indices.+1])
        @test ys[n_points+1:end] == last.(contour_vertices[expected_contour_indices.+1])
        @test all(==(-abs(contour_depth)), bottom_heights[n_points+1:end])

        # The single bulk transform now used by `sample_bathymetry_points!` must match
        # transforming each point individually — the per-point behavior it replaces.
        n = length(xs)
        bulk_xs, bulk_ys = copy(xs), copy(ys)

        ArchGDAL.importEPSG(25833; order = :trad) do source_srs
            ArchGDAL.importEPSG(4326; order = :trad) do target_srs
                ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                    ArchGDAL.transform!(bulk_xs, bulk_ys, zeros(Float64, n), transform)

                    for i in 1:n
                        point = ArchGDAL.createpoint(xs[i], ys[i])
                        ArchGDAL.transform!(point, transform)
                        px, py, _ = ArchGDAL.getpoint(point, 0)
                        @test bulk_xs[i] ≈ px
                        @test bulk_ys[i] ≈ py
                    end
                end
            end
        end
    end
end

@testset "Forcing" begin
    @testset "preparation helpers" begin
        forcing_dimension_names = FjordSim.Forcing.forcing_dimension_names
        source_fill = FjordSim.Forcing.source_fill
        fill_source! = FjordSim.Forcing.fill_source!
        solve_vertical_faces = FjordSim.Forcing.solve_vertical_faces
        nearest_valid_map = FjordSim.Forcing.nearest_valid_map
        daily_time_steps = FjordSim.Forcing.daily_time_steps
        SourceRecord = FjordSim.Forcing.SourceRecord

        # Staggered variables are written on face dimensions, matching the layout
        # `forcing_from_file` builds its FieldTimeSeries with.
        @test forcing_dimension_names("T") == ("Nx", "Ny", "Nz", "time")
        @test forcing_dimension_names("S") == ("Nx", "Ny", "Nz", "time")
        @test forcing_dimension_names("u") == ("Nx_faces", "Ny", "Nz", "time")
        @test forcing_dimension_names("v") == ("Nx", "Ny_faces", "Nz", "time")

        # An Oceananigans grid is specified by faces with centres at face midpoints, so putting
        # centres on NorKyst's non-uniform depth levels needs the faces solved for. The feasible
        # interval for the deepest face is narrow, so this is checked exactly.
        norkyst_depths = [0.0, 3, 10, 15, 25, 50, 75, 100, 150, 200, 250, 300, 500, 1000, 2000, 3000]
        faces = solve_vertical_faces(norkyst_depths)
        @test length(faces) == length(norkyst_depths) + 1
        @test all(>(0), diff(faces))  # monotonic, so Oceananigans accepts it
        @test (faces[1:end-1] .+ faces[2:end]) ./ 2 ≈ -reverse(norkyst_depths)  # centres land on the levels
        # A depth list whose implied faces oscillate for every seed must fail loudly rather than
        # silently building a grid with misplaced levels. Tightly spaced deep levels under a wide
        # shallow gap do it: the fourth level needs `2(c3 - c2) < c4 - c1`.
        @test_throws ErrorException solve_vertical_faces([0.0, 1.0, 9.0, 10.0])

        # The source mask is filled before interpolation, because `interpolate` would otherwise
        # propagate NaN out of every stencil touching land. A masked cell takes its nearest valid
        # neighbour; a level with no data at all inherits the deepest level that has some, which is
        # what the old vertical clamping achieved. Nothing may stay NaN — `interpolate` has no
        # NaN-awareness, and NumericalEarth's inpaint_mask! would have turned such a level into zeros.
        slab = Float32[
            1.0 2.0
            NaN 4.0
        ]
        slab = reshape(vcat(vec(slab), fill(NaN32, 4)), 2, 2, 2)   # level 2 entirely masked
        mask_fill = source_fill(isfinite.(slab))
        @test mask_fill.level_valid == [true, false]
        filled = fill_source!(copy(slab), mask_fill)
        @test count(isnan, filled) == 0
        @test filled[2, 1, 1] == 1.0f0             # nearest valid neighbour of the masked cell
        @test filled[:, :, 2] == filled[:, :, 1]   # dry level inherits the deepest level with data
        @test filled[1, 1, 1] == 1.0f0 && filled[2, 2, 1] == 4.0f0  # valid cells untouched
        @test_throws ErrorException source_fill(falses(2, 2, 2))  # no data anywhere
        @test_throws DimensionMismatch fill_source!(zeros(Float32, 3, 3, 3), mask_fill)

        # `prepare_forcing` writes zero lambdas everywhere: the interior relaxation band it used to
        # carry along the open edge is gone, because the boundary itself now nudges towards hourly
        # exterior data and a band relaxing the same variables a few cells in would fight it.
        # `add_rivers` writes the only nonzero lambdas the file ever carries.

        # Every cell resolves to the nearest valid source cell, so a fjord arm NorKyst treats as
        # land still receives data. Valid cells resolve to themselves.
        valid = falses(3, 3, 1)
        valid[2, 2, 1] = true
        @test all(==(Int32(5)), nearest_valid_map(valid, 1))  # 2 + (2 - 1) * 3
        valid[1, 1, 1] = true
        nearest = nearest_valid_map(valid, 1)
        @test nearest[1, 1] == 1
        @test nearest[2, 2] == 5
        @test nearest[3, 3] == 5  # closer to (2, 2) than to (1, 1)

        # A short month leaves a hole in the cyclical forcing period, so a missing day is the linear
        # blend of the records bracketing it rather than a copy of the previous one.
        records = [
            SourceRecord(DateTime(2020, 1, 1, 12), "a.nc", 1),
            SourceRecord(DateTime(2020, 1, 2, 12), "a.nc", 2),
            SourceRecord(DateTime(2020, 1, 5, 12), "b.nc", 7),
        ]
        steps = @test_logs (:warn,) daily_time_steps(records)
        @test [step.date for step in steps] == [DateTime(2020, 1, d, 12) for d in 1:5]
        # Present days read a single record.
        @test all(iszero(steps[n].weight) for n in (1, 2, 5))
        @test steps[2].lower === steps[2].upper === records[2]
        # The two missing days blend 2020-01-02 with 2020-01-05 at 1/3 and 2/3.
        @test steps[3].lower === records[2] && steps[3].upper === records[3]
        @test steps[3].weight ≈ 1 / 3
        @test steps[4].weight ≈ 2 / 3
        # A gap-free run needs no blending and emits no warning.
        contiguous = @test_logs daily_time_steps(records[1:2])
        @test [step.date for step in contiguous] == [record.date for record in records[1:2]]
        @test all(iszero(step.weight) for step in contiguous)

        # A gap longer than max_gap is still filled, but warned about separately: MET's archive is
        # missing whole weeks in 2017 and 2018, and papering over those silently would mislead.
        long_gap = [records[1], SourceRecord(DateTime(2020, 1, 20, 12), "c.nc", 3)]
        @test_logs (:warn,) (:warn,) daily_time_steps(long_gap; max_gap = 6)
        @test_logs (:warn,) daily_time_steps(long_gap; max_gap = 30)  # only the summary warning

        # Sub-daily cadences are left alone rather than being collapsed onto a daily axis.
        hourly = [
            SourceRecord(DateTime(2020, 1, 1, 0), "a.nc", 1),
            SourceRecord(DateTime(2020, 1, 1, 1), "a.nc", 2),
        ]
        hourly_steps = @test_logs (:warn,) daily_time_steps(hourly)
        @test [step.date for step in hourly_steps] == [record.date for record in hourly]
        @test all(iszero(step.weight) for step in hourly_steps)
    end

    @testset "forcing terms" begin
        forcing_term_x_flux = FjordSim.Forcing.forcing_term_x_flux
        forcing_term_y_flux = FjordSim.Forcing.forcing_term_y_flux
        forcing_term_relax = FjordSim.Forcing.forcing_term_relax
        Ax = Oceananigans.Operators.Ax
        Ay = Oceananigans.Operators.Ay
        volume = Oceananigans.Operators.volume

        mktempdir() do tmp
            (; grid) = immersed_test_grid(joinpath(tmp, "bathymetry.nc"); size = (2, 3, 2))
            i, j, k = 1, 1, 1
            Ax_over_V = Ax(i, j, k, grid, Center(), Center(), Center()) /
                        volume(i, j, k, grid, Center(), Center(), Center())
            Ay_over_V = Ay(i, j, k, grid, Center(), Center(), Center()) /
                        volume(i, j, k, grid, Center(), Center(), Center())

            # Each formula is pinned against Oceananigans' own area/volume operators directly, rather
            # than re-deriving the expected numbers by hand.
            λ, flux = 3.7, 2.5
            @test forcing_term_x_flux(λ, flux, i, j, k, grid) == flux * Ax_over_V
            @test forcing_term_y_flux(λ, flux, i, j, k, grid) == flux * Ay_over_V

            field, value = 12.0, 5.0
            @test forcing_term_relax(0.3, value, i, j, k, grid, field) == -0.3 * (field - value)

            # The `ForcingFromFile` kernel routes to the matching term by λ's sign/magnitude. No
            # producer in this repo ever writes `|λ| > 1` (`prepare_forcing` writes zeros and the
            # river forcing stays inside `(-1, 1)`), so this is the only coverage of that routing.
            forcing_config = test_forcing_config(data_root = tmp)

            function forced_value(lambda)
                write_prepared_forcing(
                    forcing_path(forcing_config);
                    size = (2, 3, 2),
                    names = ("T",),
                    value = (name, index) -> 4.0f0,
                    lambda,
                )
                forcing = forcing_from_file(forcing_config; grid, tracers = (:T,))
                clock = (; time = 0.0)
                fields = (; T = zeros(Float64, 2, 3, 2))
                return forcing.T(i, j, k, grid, clock, fields)
            end

            @test forced_value(2.0f0) ≈ 4.0 * Ax_over_V     # λ > 1 -> x-flux
            @test forced_value(-2.0f0) ≈ 4.0 * Ay_over_V    # λ < -1 -> y-flux
            # λ read back from the (Float32-written) file is promoted to Float64 from the *rounded*
            # Float32 value, so the expected value is derived the same way rather than from a Float64
            # literal that would differ in the last bit.
            @test forced_value(0.3f0) ≈ -Float64(0.3f0) * (0.0 - 4.0)  # |λ| < 1 -> relaxation, field = 0
        end
    end

    @testset "water mask" begin
        water_mask = FjordSim.Forcing.water_mask

        mktempdir() do tmp
            # A 4x2 channel deepening to the east. Column 1 is land, column 2's bottom sits inside
            # the deep cell (a partial cell), columns 3-4 are fully wet. Both rows are identical, so
            # the staggered assertions below can compare against the tracer row directly.
            (; grid) = immersed_test_grid(
                joinpath(tmp, "bathymetry.nc");
                size = (4, 2, 2),
                longitude = (10.0, 10.4),
                latitude = (59.0, 59.2),
                bottom_height = [0.0 0.0; -13.0 -13.0; -20.0 -20.0; -20.0 -20.0],
            )

            # The mask must agree with what the model treats as wet, which for PartialCellBottom is
            # decided by minimum_fractional_cell_height (0.2), not by a cell-centre test. Column 2's
            # deep cell keeps 7 m of its 10 m, so it is wet even though its centre (-15) is below the
            # bottom height (-13) — the old hand-rolled mask marked it land.
            tracer = water_mask(grid, Center, Center, :south)
            @test size(tracer) == (4, 2, 2)
            @test !tracer[1, 1, 1]      # land column, both levels
            @test !tracer[1, 1, 2]
            @test tracer[2, 1, 1]       # partial cell, still wet
            @test tracer[3, 1, 1]       # fully wet
            @test tracer[4, 1, 1]

            # A velocity face is land when either tracer cell it separates is land, matching
            # Oceananigans: a face against land is a wall. The old mask called it wet.
            u = water_mask(grid, Face, Center, :south)
            @test size(u) == (5, 2, 2)
            @test !u[2, 1, 1]     # between land column 1 and wet column 2 -> wall
            @test u[3, 1, 1]      # between two wet columns
            @test !u[1, 1, 1]     # closed west wall
            @test !u[5, 1, 1]     # closed east wall

            # ...but an open boundary the config names must survive, because that is where the
            # setup puts its radiating NormalFlowBoundaryCondition and where `prepare_boundaries`
            # samples the exterior velocity. peripheral_node alone marks it land, since the tracer
            # cell outside the domain is an inactive halo cell.
            v_south = water_mask(grid, Center, Face, :south)
            @test size(v_south) == (4, 3, 2)
            @test v_south[:, 1, :] == tracer[:, 1, :]      # southern row restored from the tracer row
            @test !any(v_south[:, 3, :])                   # northern row still a closed wall
            v_north = water_mask(grid, Center, Face, :north)
            @test !any(v_north[:, 1, :])                   # southern row now closed
            @test v_north[:, 3, :] == tracer[:, 2, :]      # northern row restored
            # A south edge must not touch the u mask, and a west edge must open its western column.
            @test water_mask(grid, Face, Center, :south)[1, 1, :] == [false, false]
            @test water_mask(grid, Face, Center, :west)[1, 1, :] == tracer[1, 1, :]

            # A setup naming no open-boundary data has no open edge: every lateral boundary stays the
            # closed wall `peripheral_node` calls it, which is the right reading with no exterior
            # state to read there.
            @test !any(water_mask(grid, Center, Face, nothing)[:, 1, :])
            @test !any(water_mask(grid, Center, Face, nothing)[:, 3, :])
            # ...and the tracer mask is edge-independent, which is why `add_rivers` needs no edge:
            # every `open_boundary_water!` method returns early unless the location is staggered
            # across its own edge.
            @test water_mask(grid, Center, Center, nothing) == water_mask(grid, Center, Center, :south)
            @test water_mask(grid, Center, Center, nothing) == water_mask(grid, Center, Center, :west)
        end
    end

    @testset "file round-trip" begin
        mktempdir() do tmp
            (; grid) = immersed_test_grid(joinpath(tmp, "bathymetry.nc"); size = (2, 3, 2))

            forcing_config = test_forcing_config(data_root = tmp)

            # A file in exactly the layout prepare_forcing writes: staggered variables on
            # face dimensions, land as the NaN fill value, times decodable to DateTime.
            write_prepared_forcing(forcing_path(forcing_config); size = (2, 3, 2))

            # The config method resolves the path itself; the filepath method still works.
            forcing = forcing_from_file(forcing_config; grid, tracers = (:T, :S))
            @test forcing isa NamedTuple
            @test Set(keys(forcing)) == Set((:T, :S, :u, :v))
            @test forcing_from_file(; grid, filepath = forcing_path(forcing_config), tracers = (:T, :S)) isa NamedTuple

            # Only requested tracers are picked up, alongside the velocities.
            @test Set(keys(forcing_from_file(forcing_config; grid, tracers = (:T,)))) == Set((:T, :u, :v))

            # `simulation_forcing` is the hook `build_simulation` calls instead of a hardcoded
            # `forcing_from_file` — dispatched on the forcing config, taking the resolved filepath
            # rather than resolving it, so it also serves the rivers-augmented copy. Every
            # built-in source shares the default, which just forwards to `forcing_from_file`.
            @test Set(keys(simulation_forcing(
                forcing_config, grid, forcing_path(forcing_config), (:T, :S), nothing,
            ))) == Set((:T, :S, :u, :v))

            # A grid the file was not written for must be rejected rather than silently misread.
            other = immersed_test_grid(
                joinpath(tmp, "bathymetry_other.nc");
                size = (4, 3, 2),
                longitude = (10.0, 12.0),
                latitude = (59.0, 62.0),
            ).grid
            @test_throws DimensionMismatch forcing_from_file(
                forcing_config;
                grid = other,
                tracers = (:T, :S),
            )
        end
    end

    @testset "OF800 rivers config" begin
        rivers = OF800RiversConfig(data_root = "/data/oslofjord")

        @test rivers isa AbstractRiverConfig
        @test all(isconcretetype, fieldtypes(typeof(rivers)))
        @test river_forcing_path(rivers) == "/data/oslofjord/forcing_rivers.nc"
        @test FjordSim.Forcing.river_locations_path(rivers) == "/data/oslofjord/OF800_rivers.csv"
        @test FjordSim.Forcing.river_series_path(rivers) == "/data/oslofjord/of800_rivers_v9_1990_2022_RA1.nc"
        @test river_search_radius(rivers) == 10
        @test rivers.relaxation_timescale == 3600.0

        # Both source files download from per-file links by default; a folder link cannot work.
        for url in (rivers.locations_url, rivers.series_url)
            @test occursin("/scl/fi/", url)   # a per-file link, not a /scl/fo/ folder archive
            @test occursin("dl=1", url)       # serves the file rather than the web app
        end

        # An expired Dropbox `st` token answers with a login page and HTTP 200, so the download
        # reports success. That page must be rejected and removed, not left to masquerade as data.
        mktempdir() do tmp
            page = joinpath(tmp, "OF800_rivers.csv")
            write(page, "<!DOCTYPE html>\n<html><title>Log in to Dropbox</title></html>")
            @test_throws ErrorException FjordSim.Forcing.validate_river_download(page, "https://example.invalid")
            @test !isfile(page)

            empty_file = joinpath(tmp, "empty.nc")
            write(empty_file, "")
            @test_throws ErrorException FjordSim.Forcing.validate_river_download(empty_file, "https://example.invalid")

            # Real data is passed through untouched: the CSV header, and the HDF5 magic the
            # NetCDF-4 series file starts with.
            good_csv = joinpath(tmp, "good.csv")
            write(good_csv, "River number,Name,LatOutlet,LonOutlet,Zero\n1,A,59.0,10.5,0\n")
            @test FjordSim.Forcing.validate_river_download(good_csv, "https://example.invalid") == good_csv
            @test isfile(good_csv)

            good_nc = joinpath(tmp, "good.nc")
            write(good_nc, UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a])
            @test FjordSim.Forcing.validate_river_download(good_nc, "https://example.invalid") == good_nc
        end
    end

    @testset "boundary preparation helpers" begin
        boundary_dimension_names = FjordSim.Forcing.boundary_dimension_names
        boundary_location = FjordSim.Forcing.boundary_location
        boundary_variable_name = FjordSim.Forcing.boundary_variable_name
        boundary_mask_slice = FjordSim.Forcing.boundary_mask_slice
        boundary_domain = FjordSim.Forcing.boundary_domain
        fill_boundary_gaps! = FjordSim.Forcing.fill_boundary_gaps!
        hourly_time_steps = FjordSim.Forcing.hourly_time_steps
        SourceRecord = FjordSim.Forcing.SourceRecord

        boundary_source_slab = FjordSim.Forcing.boundary_source_slab
        blended_slab = FjordSim.Forcing.blended_slab

        # `ubar`/`vbar` are the one pair NorKyst publishes along its own curvilinear axes rather than
        # derotated, so they are rotated to eastward/northward on the source grid. Every other
        # variable, and every other source, keeps the plain `blended_slab` read.
        mktempdir() do tmp
            dates = [DateTime(2020, 1, 1), DateTime(2020, 1, 1, 1)]
            config = test_boundaries_config(data_root = tmp, architecture = :cpu)
            mkpath(boundary_data_directory(config))
            stub = write_hourly_source_stub(
                joinpath(
                    boundary_data_directory(config),
                    FjordSim.Forcing.boundary_monthly_filename(config, 2020, 1),
                );
                dates,
                longitude = (10.0, 11.0),
                latitude = (59.0, 60.0),
                # A quarter turn makes the rotation unmistakable: east becomes north exactly.
                grid_angle = π / 2,
                value = (name, level, index) -> name == "ubar" ? 3.0f0 :
                                                name == "vbar" ? 4.0f0 : Float32(level + index),
            )

            steps = hourly_time_steps([SourceRecord(date, stub, index)
                                       for (index, date) in enumerate(dates)])
            reader = FjordSim.Forcing.SourceReader(stub)
            try
                step = first(steps)

                # angle = pi/2: (3, 4) -> (-4, 3).
                @test all(≈(-4.0f0), boundary_source_slab(config, reader, step, "ubar"))
                @test all(≈(3.0f0), boundary_source_slab(config, reader, step, "vbar"))

                # A rotation preserves speed, which the raw pair does not have in common with a
                # component-wise rescale — this is what catches a sign or transpose slip.
                rotated = hypot.(
                    boundary_source_slab(config, reader, step, "ubar"),
                    boundary_source_slab(config, reader, step, "vbar"),
                )
                @test all(≈(5.0f0), rotated)

                # Everything else is the untouched default, byte for byte.
                for name in ("temperature", "salinity", "u_eastward", "v_northward", "zeta")
                    @test boundary_source_slab(config, reader, step, name) ==
                          blended_slab(reader, step, name)
                end

                # And the supertype default is what a source that overrides nothing gets, including
                # for the two velocity names.
                other = test_boundaries_config(data_root = tmp)
                invoke_default = FjordSim.Configs.AbstractBoundaryDataConfig
                for name in ("ubar", "vbar", "zeta")
                    @test invoke(
                        boundary_source_slab,
                        Tuple{invoke_default,Any,Any,Any},
                        other, reader, step, name,
                    ) == blended_slab(reader, step, name)
                end
            finally
                close(reader)
            end
        end

        # The side is in the variable name, not the filename, so one file holds several boundaries.
        @test boundary_variable_name(:south, "T") == "south_T"
        @test boundary_variable_name(:west, "eta") == "west_eta"
        @test_throws ArgumentError boundary_variable_name(:middle, "T")

        # A boundary variable drops the dimension across the edge, and a surface variable drops the
        # vertical one too. The staggering follows the variable, as in the forcing file.
        @test boundary_dimension_names(Val(:south), "T") == ("Nx", "Nz", "time")
        @test boundary_dimension_names(Val(:south), "u") == ("Nx_faces", "Nz", "time")
        @test boundary_dimension_names(Val(:south), "v") == ("Nx", "Nz", "time")
        @test boundary_dimension_names(Val(:south), "eta") == ("Nx", "time")
        @test boundary_dimension_names(Val(:south), "vbar") == ("Nx", "time")
        @test boundary_dimension_names(Val(:west), "T") == ("Ny", "Nz", "time")
        @test boundary_dimension_names(Val(:west), "v") == ("Ny_faces", "Nz", "time")
        @test boundary_dimension_names(Val(:east), "ubar") == ("Ny", "time")
        @test_throws ArgumentError boundary_dimension_names(Val(:middle), "T")

        # The reduced locations are exactly the flavours `getbc` accepts as a boundary condition —
        # `XZFTS` on a south or north edge, `YZFTS` on a west or east one — which is what lets a
        # series be passed straight through with no discrete-form wrapper.
        @test boundary_location(Val(:south), "T") == (Center, Nothing, Center)
        @test boundary_location(Val(:south), "u") == (Face, Nothing, Center)
        @test boundary_location(Val(:south), "eta") == (Center, Nothing, Nothing)
        @test boundary_location(Val(:west), "v") == (Nothing, Face, Center)
        @test boundary_location(Val(:west), "ubar") == (Nothing, Center, Nothing)

        # The boundary row is one cell thick across the edge, so every three-dimensional structure
        # in the forcing core applies to it unchanged. A surface variable also takes the top level,
        # which is the surface because `znodes` increases upwards.
        mask = reshape(1:24, 4, 3, 2) .> 0
        @test size(boundary_mask_slice(Val(:south), mask, false)) == (4, 1, 2)
        @test size(boundary_mask_slice(Val(:south), mask, true)) == (4, 1, 1)
        @test size(boundary_mask_slice(Val(:west), mask, false)) == (1, 3, 2)
        marked = trues(4, 3, 2)
        marked[:, 3, :] .= false
        @test all(boundary_mask_slice(Val(:south), marked, false))     # the low end of the axis
        @test !any(boundary_mask_slice(Val(:north), marked, false))    # the high end

        # The download covers the full extent along the edge and only a margin across it — a thin
        # band, because an hourly whole-domain download is more than twenty times the interior one.
        grid = LatitudeLongitudeGrid(
            CPU(); size = (2, 2, 1), halo = (1, 1, 1),
            longitude = (10.0, 11.0), latitude = (59.0, 60.0), z = (-10.0, 0.0),
        )
        longitude, latitude = boundary_domain(Val(:south), grid, 0.1)
        @test all(isapprox.(longitude, (9.9, 11.1)))
        @test all(isapprox.(latitude, (58.9, 59.1)))
        longitude, latitude = boundary_domain(Val(:east), grid, 0.1)
        @test all(isapprox.(longitude, (10.9, 11.1)))
        @test all(isapprox.(latitude, (58.9, 60.1)))
        @test_throws ArgumentError boundary_domain(Val(:middle), grid, 0.1)

        # A dry boundary cell is `NaN` in the file and is filled from its nearest wet neighbour on
        # read. Not cosmetic: the tangential velocity's faces at both ends of a south edge are
        # peripheral, so the file always has holes there, and `NormalRadiation` does not zero them —
        # `immersed_peripheral_node` is false at a wet domain-edge face.
        # The middle gap takes whichever side is nearer, not whichever comes first.
        line = reshape(Float64[NaN, 2.0, NaN, NaN, 5.0], 5, 1, 1, 1)
        @test vec(fill_boundary_gaps!(copy(line))) == [2.0, 2.0, 2.0, 5.0, 5.0]
        @test all(iszero, fill_boundary_gaps!(fill(NaN, 3, 1, 1, 1)))  # a wholly dry line
        # A west edge is reduced in x instead, so the line lies along the second axis.
        west_line = reshape(Float64[NaN, 4.0, NaN], 1, 3, 1, 1)
        @test vec(fill_boundary_gaps!(copy(west_line))) == [4.0, 4.0, 4.0]

        # The boundary axis is hourly rather than daily, which is the whole reason it is a separate
        # file. Gap filling is the same argument as for the interior forcing, one unit down.
        hours = [
            SourceRecord(DateTime(2020, 1, 1, hour), "one.nc", hour + 1) for hour in (0, 1, 4)
        ]
        steps = @test_logs (:warn,) hourly_time_steps(hours)
        @test [step.date for step in steps] ==
              [DateTime(2020, 1, 1, hour) for hour in 0:4]
        @test steps[3].weight ≈ 1 / 3          # 02:00 is one third of the way across the gap
        @test iszero(steps[2].weight)          # a downloaded hour reads a single record
    end

    @testset "boundary round-trip" begin
        # The only test in this suite that drives a real regrid-and-write pipeline end to end:
        # `write_prepared_forcing` hand-rolls the *output* layout, so a rank or staggering bug in
        # `write_boundaries_file` would otherwise be invisible. It also pins the two things the
        # boundary file does that the forcing file does not — a one-cell-thick target slab, and
        # surface variables interpolated on a single-level source grid.
        mktempdir() do tmp
            longitude, latitude = (10.0, 11.0), (59.0, 60.0)
            dates = [DateTime(2020, 1, 1) + Hour(hour) for hour = 0:2]

            # A land column on the southern boundary row, so a dry boundary cell is exercised.
            bottom_height = fill(-20.0, 4, 3)
            bottom_height[2, 1] = 0.0
            (; grid) = immersed_test_grid(
                joinpath(tmp, "bathymetry.nc");
                size = (4, 3, 2),
                longitude,
                latitude,
                z_faces = [-20.0, -10.0, 0.0],
                bottom_height,
            )

            boundaries_config = test_boundaries_config(data_root = tmp, architecture = :cpu)
            mkpath(boundary_data_directory(boundaries_config))
            write_hourly_source_stub(
                joinpath(
                    boundary_data_directory(boundaries_config),
                    FjordSim.Forcing.boundary_monthly_filename(boundaries_config, 2020, 1),
                );
                dates,
                # Wider than the target grid, so every target node has source cells around it.
                longitude = (longitude[1] - 0.2, longitude[2] + 0.2),
                latitude = (latitude[1] - 0.2, latitude[2] + 0.2),
            )

            result = prepare_boundaries(grid, boundaries_config)
            @test result.output_file == boundary_data_path(boundaries_config)
            @test result.times == dates
            @test Set(result.variables) == Set(
                "south_" * name for name in ("T", "S", "u", "v", "eta", "ubar", "vbar")
            )

            NCDataset(result.output_file) do ds
                # The file states its own side as well as its variable names do.
                @test ds.attrib["open_edge"] == "south"

                # Shapes and staggering, straight from `boundary_dimension_names`.
                @test size(ds["south_T"]) == (4, 2, length(dates))
                @test size(ds["south_u"]) == (5, 2, length(dates))
                @test size(ds["south_eta"]) == (4, length(dates))
                @test size(ds["south_ubar"]) == (5, length(dates))
                @test NCDatasets.dimnames(ds["south_v"]) == ("Nx", "Nz", "time")
                @test NCDatasets.dimnames(ds["south_vbar"]) == ("Nx", "time")

                # No `_lambda` twins: this file is boundary data, not a relaxation forcing. The
                # timescales are the boundary *condition*'s, not the dataset's.
                @test !any(occursin("_lambda", name) for name in keys(ds))

                # The land column is `NaN` in the file, exactly as `water_mask` marks it.
                @test all(ismissing, ds["south_T"][2, :, 1])
                @test ismissing(ds["south_eta"][2, 1])
                @test !any(ismissing, ds["south_T"][1, :, 1])

                # A surface variable is interpolated at the single-level source grid's own cell
                # centre, so it comes back exactly rather than interpolated in the vertical.
                @test ds["south_eta"][1, 2] ≈ 3.0f0      # value(name, level = 1, index = 2)

                # A full-depth variable lands on the target cell centres: -5 m is the source's own
                # second level, -15 m is two thirds of the way to its third.
                @test ds["south_T"][1, 2, 1] ≈ 3.0f0
                @test ds["south_T"][1, 1, 1] ≈ 3.0f0 + 2 / 3
            end

            # Reading it back gives reduced series on the model grid, keyed by the bare names, with
            # the dry column filled and the time axis zeroed at the instant the run starts.
            all_series = FjordSim.boundary_series(boundaries_config, grid, DateTime(2020, 1, 1))
            @test keys(all_series) == (:south,)
            series = all_series.south
            @test Set(keys(series)) == Set((:T, :S, :u, :v, :eta, :ubar, :vbar))
            @test size(series.T) == (4, 1, 2, length(dates))
            @test size(series.u) == (5, 1, 2, length(dates))
            @test size(series.eta) == (4, 1, 1, length(dates))
            @test series.T.times == [0, 3600, 7200]
            @test all(isfinite, interior(series.T, :, :, :, :))
            @test interior(series.eta, 2, 1, 1, 1) ≈ interior(series.eta, 1, 1, 1, 1)  # filled

            # An `Hour(1)` axis zeroed at a later instant is negative before it, which is what keeps
            # the boundary in phase with the forcing and the atmosphere rather than with its own
            # first record.
            shifted = FjordSim.boundary_series(boundaries_config, grid, DateTime(2020, 1, 1, 2))
            @test shifted.south.T.times == [-7200, -3600, 0]

            # A file written for one edge cannot be read as a config that opens another: the
            # variables would be missing, and the stated edge says so first.
            @test_throws ErrorException FjordSim.boundary_series(
                test_boundaries_config(data_root = tmp, open_edges = :north), grid, DateTime(2020, 1, 1),
            )

            # The plot reads both panel shapes out of the same file, taking the side from the config.
            @test isfile(plot_boundaries(boundaries_config))
        end
    end

    @testset "open ocean round-trip" begin
        # A domain open on all four sides, through the real regrid-and-write path: one file with four
        # sides in it, which is what the layout was built for — variables named per side, all six
        # spatial dimensions defined, one time axis. The single-edge round-trip above cannot catch a
        # writer or reader that assumed one side.
        mktempdir() do tmp
            longitude, latitude = (10.0, 11.0), (59.0, 60.0)
            dates = [DateTime(2020, 1, 1) + Hour(hour) for hour = 0:2]
            edges = (:south, :north, :west, :east)

            (; grid) = immersed_test_grid(
                joinpath(tmp, "bathymetry.nc");
                size = (4, 3, 2),
                halo = (1, 1, 1),
                longitude,
                latitude,
                z_faces = [-20.0, -10.0, 0.0],
                bottom_height = fill(-20.0, 4, 3),
            )

            config = test_boundaries_config(data_root = tmp, architecture = :cpu, open_edges = edges)
            mkpath(boundary_data_directory(config))
            write_hourly_source_stub(
                joinpath(
                    boundary_data_directory(config),
                    FjordSim.Forcing.boundary_monthly_filename(config, 2020, 1),
                );
                dates,
                longitude = (longitude[1] - 0.2, longitude[2] + 0.2),
                latitude = (latitude[1] - 0.2, latitude[2] + 0.2),
            )

            # The download box is the bounding box of the four bands, which for opposite edges is the
            # whole domain — the thin-band saving is a single-edge one, and this is where that shows.
            box_longitude, box_latitude = FjordSim.Forcing.boundary_domain(edges, grid, 0.05)
            @test box_longitude == (longitude[1] - 0.05, longitude[2] + 0.05)
            @test box_latitude == (latitude[1] - 0.05, latitude[2] + 0.05)

            result = prepare_boundaries(grid, config)
            @test length(result.variables) == 4 * 7
            @test Set(result.variables) ==
                  Set("$(edge)_$name" for edge in edges
                      for name in ("T", "S", "u", "v", "eta", "ubar", "vbar"))

            NCDataset(result.output_file) do ds
                # One file, four sides, one time axis, and the attribute lists every side.
                @test Set(Symbol.(split(ds.attrib["open_edge"], ","))) == Set(edges)
                @test ds.dim["time"] == length(dates)
                # The along-boundary axis and the staggering follow the edge: a south section runs
                # along x, a west one along y, and each normal velocity is on its own face axis.
                @test size(ds["south_T"]) == (4, 2, length(dates))
                @test size(ds["west_T"]) == (3, 2, length(dates))
                @test size(ds["south_v"]) == (4, 2, length(dates))
                @test size(ds["west_u"]) == (3, 2, length(dates))
                @test size(ds["north_eta"]) == (4, length(dates))
                # ...and a *tangential* component is staggered along the edge it runs down, so `vbar`
                # on an east edge sits on y-faces: one node more than the tracer row beside it.
                @test size(ds["east_vbar"]) == (4, length(dates))
                @test size(ds["south_ubar"]) == (5, length(dates))
            end

            series = FjordSim.boundary_series(config, grid, DateTime(2020, 1, 1))
            @test Set(keys(series)) == Set(edges)
            for edge in edges
                @test Set(keys(getproperty(series, edge))) ==
                      Set((:T, :S, :u, :v, :eta, :ubar, :vbar))
                @test all(isfinite, interior(getproperty(series, edge).T, :, :, :, :))
            end
            # Each group is reduced across its own edge, so a south series is one cell thick in y and
            # a west series one cell thick in x.
            @test size(series.south.T) == (4, 1, 2, length(dates))
            @test size(series.west.T) == (1, 3, 2, length(dates))

            # The forcing land mask opens all four velocity face rows at once.
            u_mask = FjordSim.Forcing.water_mask(grid, Face, Center, edges)
            v_mask = FjordSim.Forcing.water_mask(grid, Center, Face, edges)
            @test any(u_mask[1, :, :]) && any(u_mask[5, :, :])
            @test any(v_mask[:, 1, :]) && any(v_mask[:, 4, :])

            # One plot, one block of rows per edge.
            @test isfile(plot_boundaries(config))
        end
    end

    @testset "river cell snapping" begin
        is_coastal_cell = FjordSim.Forcing.is_coastal_cell
        nearest_coastal_cell = FjordSim.Forcing.nearest_coastal_cell
        coastal_water_mask = FjordSim.Forcing.coastal_water_mask
        river_cells = FjordSim.Forcing.river_cells

        # A hand-built mask exercises the two rules on their own: a cell is coastal when it is water
        # and touches land, so open water in the middle is not a valid river mouth.
        mask = trues(5, 5)
        mask[1, :] .= false          # a land column along the western edge
        @test !is_coastal_cell(mask, 1, 3)   # land itself
        @test is_coastal_cell(mask, 2, 3)    # water touching the land column
        @test !is_coastal_cell(mask, 4, 3)   # open water, no land neighbour
        @test !is_coastal_cell(mask, 0, 3)   # outside the grid

        # An open-water cell is pulled to the coast — column 2 is the only coastal column here, so
        # the search has to walk out to it rather than stopping at the first water cell it meets.
        @test nearest_coastal_cell(mask, 2, 3, 10) == (2, 3, 0.0)
        @test nearest_coastal_cell(mask, 4, 3, 10) == (2, 3, 2.0)
        @test isnothing(nearest_coastal_cell(mask, 5, 3, 1))  # coast is 3 cells away, radius is 1

        # Ties are broken by iteration order — latitude offset outermost, ascending — so of the two
        # cells one step away from (3, 3) the one with the lower j wins.
        tie_mask = trues(5, 5)
        tie_mask[3, 1] = false
        tie_mask[3, 5] = false
        @test nearest_coastal_cell(tie_mask, 3, 3, 10) == (3, 2, 1.0)

        mktempdir() do tmp
            # A 5x4 basin with a land column at i = 2, so columns 1 and 3 are coastal while columns
            # 4 and 5 are open water. The land is off the domain edge because an outlet has to sit
            # strictly inside the grid to be accepted at all.
            (; grid) = land_column_test_grid(joinpath(tmp, "bathymetry.nc"))

            # The mask comes from the same water_mask prepare_forcing uses, taken at the surface.
            mask = coastal_water_mask(grid)
            @test size(mask) == (5, 4)
            @test !mask[2, 1]
            @test mask[1, 1]
            @test mask[3, 1]

            longitudes = Array(Oceananigans.Grids.λnodes(grid, Center()))
            latitudes = Array(Oceananigans.Grids.φnodes(grid, Center()))

            # An outlet on the land column relocates to the coast; one in open water does too; one
            # outside the domain is dropped rather than clamped to the nearest edge cell.
            locations = [
                FjordSim.Forcing.RiverLocation(1, "on land", longitudes[2], latitudes[2]),
                FjordSim.Forcing.RiverLocation(2, "open water", longitudes[4], latitudes[2]),
                FjordSim.Forcing.RiverLocation(3, "outside", 20.0, latitudes[2]),
            ]
            cells = river_cells(grid, locations, 10)
            @test length(cells) == 2
            @test [cell.location.id for cell in cells] == [1, 2]
            @test (cells[1].i, cells[1].j, cells[1].distance) == (1, 2, 1.0)
            @test (cells[2].i, cells[2].j, cells[2].distance) == (3, 2, 1.0)
            @test all(mask[cell.i, cell.j] for cell in cells)

            # An outlet sitting exactly on the outermost node counts as outside, matching the
            # reference's strict bounds test.
            edge = FjordSim.Forcing.RiverLocation(4, "on the edge", longitudes[1], latitudes[2])
            @test isempty(river_cells(grid, [edge], 10))

            # With no reach, the on-land outlet is dropped too rather than written into land.
            @test isempty(river_cells(grid, [locations[1]], 0))
        end
    end

    @testset "add rivers round-trip" begin
        mktempdir() do tmp
            (; grid) = land_column_test_grid(joinpath(tmp, "bathymetry.nc"))

            longitudes = Array(Oceananigans.Grids.λnodes(grid, Center()))
            latitudes = Array(Oceananigans.Grids.φnodes(grid, Center()))

            # One outlet on the land column, which relocates to the coastal cell (1, 2). "N" is a
            # variable this forcing file does not carry, so it must be skipped, not fail.
            rivers = StubRivers(
                tmp,
                "forcing_rivers.nc",
                3600.0,
                10,
                [FjordSim.Forcing.RiverLocation(7, "test river", longitudes[2], latitudes[2])],
                Dict(
                    "T" => Float32[3.0 4.0],
                    "S" => Float32[0.0 0.0],
                    "N" => Float32[9.0 9.0],
                ),
            )
            forcing_config = test_forcing_config(data_root = tmp, rivers = rivers)
            @test forcing_config.rivers isa AbstractRiverConfig
            @test all(isconcretetype, fieldtypes(typeof(forcing_config)))

            # The boundary lambda is nonzero, so the assertions below can tell a river cell's own
            # lambda from the band's.
            write_prepared_forcing(forcing_path(forcing_config); size = (5, 4, 2), lambda = 2.0f-5)

            result = add_rivers(grid, forcing_config)
            @test result.output_file == joinpath(tmp, "forcing_rivers.nc")
            @test result.variables == ["S", "T"]  # "N" is not in the file, so it is skipped
            @test length(result.cells) == 1
            @test (result.cells[1].i, result.cells[1].j) == (1, 2)

            NCDataset(result.output_file) do written
                # The river cell carries the river values and the river lambda at the surface, for
                # every time step.
                @test written["T"][1, 2, 2, :] == Float32[3.0, 4.0]
                @test written["S"][1, 2, 2, :] == Float32[0.0, 0.0]
                @test all(written["T_lambda"][1, 2, 2, :] .≈ Float32(1 / 3600))
                @test all(written["S_lambda"][1, 2, 2, :] .≈ Float32(1 / 3600))

                # ...and nothing else moves: not the level below, not a neighbouring cell, not the
                # variables the river dataset did not name.
                @test written["T"][1, 2, 1, :] == Float32[1.0, 1.0]
                @test all(written["T_lambda"][1, 2, 1, :] .== 2.0f-5)
                @test written["T"][1, 3, 2, :] == Float32[1.0, 1.0]
                @test all(written["T_lambda"][1, 3, 2, :] .== 2.0f-5)
                @test all(written["u"][:, :, :, :] .== 1.0f0)
                @test !haskey(written, "N")
            end

            # The prepared forcing itself is left alone, so the step can be re-run.
            NCDataset(forcing_path(forcing_config)) do original
                @test all(original["T"][:, :, :, :] .== 1.0f0)
                @test all(original["T_lambda"][:, :, :, :] .== 2.0f-5)
            end
        end
    end

    @testset "standalone rivers" begin
        mktempdir() do tmp
            (; grid) = land_column_test_grid(joinpath(tmp, "bathymetry.nc"))

            longitudes = Array(Oceananigans.Grids.λnodes(grid, Center()))
            latitudes = Array(Oceananigans.Grids.φnodes(grid, Center()))
            location = FjordSim.Forcing.RiverLocation(7, "test river", longitudes[2], latitudes[2])

            # Two days plus the step past the window's end: a daily axis anchored at `start_date`
            # covering `start_date + 1day`, which is what `river_forcing_times` builds.
            coverage = (DateTime(2020, 1, 1), DateTime(2020, 1, 2))
            dates = [DateTime(2020, 1, 1), DateTime(2020, 1, 2)]
            rivers(; standalone) = StubRivers(
                tmp,
                "forcing_rivers.nc",
                3600.0,
                10,
                [location],
                Dict("T" => Float32[3.0 4.0], "S" => Float32[0.0 0.0]),
                standalone,
            )
            forcing_config = test_forcing_config(data_root = tmp, rivers = rivers(standalone = true))
            @test all(isconcretetype, fieldtypes(typeof(forcing_config)))

            # No prepared forcing exists, and none is needed: the whole point of `standalone`.
            @test !isfile(forcing_path(forcing_config))
            result = add_rivers(grid, forcing_config; coverage)
            @test result.output_file == joinpath(tmp, "forcing_rivers.nc")
            @test result.times == dates
            @test result.variables == ["S", "T"]
            @test !isfile(forcing_path(forcing_config))

            NCDataset(result.output_file) do written
                # The layout is the prepared-forcing contract, so `forcing_from_file`'s dimension
                # check against the grid holds — that is what makes the file readable at all.
                for (name, expected) in
                    ("Nx" => 5, "Ny" => 4, "Nz" => 2, "Nx_faces" => 6, "Ny_faces" => 5, "time" => 2)
                    @test written.dim[name] == expected
                end
                @test DateTime.(written["time"][:]) == dates
                @test written.attrib["rivers_only"] == "true"

                # The river cell carries its values and lambda at the surface, for every step...
                @test written["T"][1, 2, 2, :] == Float32[3.0, 4.0]
                @test all(written["T_lambda"][1, 2, 2, :] .≈ Float32(1 / 3600))

                # ...and every other cell is empty rather than zero-valued, which the reader turns
                # into the same land sentinel a dry cell gets: no forcing anywhere but the river.
                @test all(ismissing, written["T"][1, 3, 2, :])
                @test all(iszero, written["T_lambda"][1, 3, 2, :])

                # Nothing but the surface level is written at all, and only the variables the river
                # dataset named exist.
                @test all(ismissing, written["T"][1, 2, 1, :])
                @test !haskey(written, "u")
            end

            # The axis reaches past the window, so `validate_time_coverage` accepts a run that ends
            # inside the last day rather than exactly on a record.
            @test first(result.times) <= first(coverage)
            @test last(result.times) >= last(coverage)

            # A standalone file has no ocean state in it, so starting a run from one is refused
            # rather than silently beginning at T = S = 0 everywhere.
            @test_throws ErrorException FjordSim.Simulations.forcing_state(
                result.output_file, grid, first(dates), (:T, :S),
            )

            # Without a run window there is no axis to build, and no prepared file to fall back on.
            @test_throws ErrorException add_rivers(grid, forcing_config)

            # ...while a non-standalone config still demands the prepared forcing, naming the step
            # that writes it — a forgotten `prepare_forcing` must not silently become a river-only run.
            plain = test_forcing_config(data_root = tmp, rivers = rivers(standalone = false))
            @test_throws ErrorException add_rivers(grid, plain; coverage)
        end
    end
end

@testset "Atmosphere" begin
    atmospheres = FjordSim.Atmospheres

    @testset "preparation helpers" begin
        uniform_centers = atmospheres.uniform_centers
        nora3_runs = atmospheres.nora3_runs
        nora3_month_dates = atmospheres.nora3_month_dates
        grid_rotation_angle = atmospheres.grid_rotation_angle
        rotate_to_east_north = atmospheres.rotate_to_east_north
        interpolate_to_target! = atmospheres.interpolate_to_target!
        projected_atmosphere_nodes = atmospheres.projected_atmosphere_nodes
        validate_target_coverage = atmospheres.validate_target_coverage
        validate_atmosphere_records = atmospheres.validate_atmosphere_records

        # The prepared axis starts on a whole multiple of the resolution, is uniform, and covers the
        # domain plus the padding on both sides.
        centers = uniform_centers((10.2, 11.02), 0.1, 0.02)
        # Both ends are computed products, so they are compared with a tolerance rather than exactly —
        # `floor(10.1 / 0.02) * 0.02` lands on 10.1 only by luck.
        @test first(centers) ≈ 10.1
        @test all(isapprox(0.02), diff(centers))
        @test last(centers) > 11.12 - 1e-9
        @test length(centers) == 52

        config = test_atmosphere_config()
        @test_throws ArgumentError atmosphere_target_axes(
            LatitudeLongitudeGrid(CPU(); size = (1, 1, 1), longitude = (10.0, 11.0), latitude = (59.0, 60.0), z = (-1.0, 0.0)),
            test_atmosphere_config(resolution = 0.0),
        )

        # One run covers six hours, so a day needs its four runs plus the previous day's 18Z run for
        # hours 00:00 to 03:00.
        runs = nora3_runs(2020, 1)
        @test first(runs) == DateTime(2019, 12, 31, 18)
        @test last(runs) == DateTime(2020, 1, 31, 18)
        @test length(runs) == 1 + 31 * 4

        # Every hour of the month is supplied exactly once, with no gap and no duplicate — the
        # property that lets all eight variables share one time axis.
        january = nora3_month_dates(2020, 1)
        @test length(january) == 31 * 24
        @test first(january) == DateTime(2020, 1, 1, 0)
        @test last(january) == DateTime(2020, 1, 31, 23)
        @test allunique(january)
        @test all(==(Hour(1)), diff(january))

        # ...including across month and year boundaries, and in a leap February.
        @test length(nora3_month_dates(2020, 2)) == 29 * 24
        @test length(nora3_month_dates(2021, 2)) == 28 * 24
        year = vcat((nora3_month_dates(2020, month) for month in 1:12)...)
        @test length(year) == 366 * 24
        @test allunique(year)
        @test all(==(Hour(1)), diff(year))
        @test first(year) == DateTime(2020, 1, 1, 0)
        @test last(year) == DateTime(2020, 12, 31, 23)

        # The URL layout is deterministic, so no catalog listing is needed.
        @test atmospheres.nora3_url(config, DateTime(2020, 1, 1, 0), 4) ==
              "https://thredds.met.no/thredds/dodsC/nora3/2020/01/01/00/fc2020010100_004_fp.nc"
        @test atmospheres.nora3_monthly_filename(config, 2020, 3) == "NORA3_202003.nc"

        # The rotation angle is measured from east to the source grid's local x axis.
        east_aligned_longitude = [10.0 10.0; 10.1 10.1; 10.2 10.2]
        east_aligned_latitude = [59.0 59.5; 59.0 59.5; 59.0 59.5]
        @test all(isapprox(0.0), grid_rotation_angle(east_aligned_longitude, east_aligned_latitude))

        north_aligned_longitude = [10.0 10.5; 10.0 10.5; 10.0 10.5]
        north_aligned_latitude = [59.0 59.0; 59.1 59.1; 59.2 59.2]
        @test all(isapprox(π / 2), grid_rotation_angle(north_aligned_longitude, north_aligned_latitude))

        # The last column repeats the previous difference rather than being left undefined.
        angle = grid_rotation_angle(east_aligned_longitude, east_aligned_latitude)
        @test size(angle) == size(east_aligned_longitude)
        @test angle[end, 1] == angle[end-1, 1]

        # Rotating by zero is the identity; by π/2 sends (u, v) to (-v, u).
        eastward, northward = rotate_to_east_north(zeros(2, 2), fill(3.0, 2, 2), fill(4.0, 2, 2))
        @test all(eastward .== 3.0) && all(northward .== 4.0)
        eastward, northward = rotate_to_east_north(fill(π / 2, 2, 2), fill(3.0, 2, 2), fill(4.0, 2, 2))
        @test all(isapprox(-4.0), eastward) && all(isapprox(3.0), northward)

        # Rotating and unrotating returns the original components.
        angle = [0.3 -1.2; 2.0 0.7]
        u = [1.0 -2.0; 3.0 0.5]
        v = [-1.5 2.5; 0.0 4.0]
        eastward, northward = rotate_to_east_north(angle, u, v)
        back_u, back_v = rotate_to_east_north(-angle, eastward, northward)
        @test all(isapprox.(back_u, u; atol = 1e-12))
        @test all(isapprox.(back_v, v; atol = 1e-12))

        # Bilinear interpolation is exact for a field linear in the projected coordinates.
        source = ProjectedAtmosphereGrid(collect(0.0:100.0:400.0), collect(0.0:100.0:300.0), "+proj=longlat +datum=WGS84")
        slab = Float32[1 + 2 * x + 3 * y for x in source.x, y in source.y]
        x = [50.0 150.0; 250.0 375.0]
        y = [25.0 175.0; 100.0 290.0]
        output = Matrix{Float32}(undef, 2, 2)
        interpolate_to_target!(output, slab, x, y, source)
        expected = Float32[1 + 2 * x[index] + 3 * y[index] for index in eachindex(x)]
        @test all(isapprox.(vec(output), expected; rtol = 1e-5))

        # A target node outside the source window is refused rather than extrapolated, because the
        # usual cause is a download that predates a change to `resolution` or `padding`.
        @test isnothing(validate_target_coverage(x, y, source))
        @test_throws ErrorException validate_target_coverage([500.0;;], [0.0;;], source)

        # An empty or unsorted record axis is an error; a hole in it is a warning.
        @test_throws ErrorException validate_atmosphere_records(AtmosphereRecord[])
        contiguous = [AtmosphereRecord(DateTime(2020, 1, 1, hour), "a.nc", hour + 1) for hour in 0:3]
        @test isnothing(validate_atmosphere_records(contiguous))
        gapped = [
            AtmosphereRecord(DateTime(2020, 1, 1, 0), "a.nc", 1),
            AtmosphereRecord(DateTime(2020, 1, 1, 5), "a.nc", 2),
        ]
        @test_logs (:warn,) validate_atmosphere_records(gapped)
    end

    @testset "file round-trip" begin
        # The real NORA3 projection, so the coordinate transform is exercised for real.
        proj4 = "+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +no_defs +R=6.371e+06"

        mktempdir() do tmp
            config = test_atmosphere_config(data_root = tmp)
            grid = LatitudeLongitudeGrid(
                CPU();
                size = (8, 10, 2),
                longitude = (10.2, 11.02),
                latitude = (59.0, 59.93),
                z = (-10.0, 0.0),
            )
            longitude, latitude = atmosphere_target_axes(grid, config)

            # A source window in projected meters that comfortably covers the prepared grid, found by
            # projecting the prepared nodes and padding the result by two source cells.
            probe = ProjectedAtmosphereGrid([0.0, 1.0], [0.0, 1.0], proj4)
            probe_x, probe_y = atmospheres.projected_atmosphere_nodes(longitude, latitude, probe)
            spacing = 3000.0
            source_x = collect(
                range(
                    (floor(minimum(probe_x) / spacing) - 2) * spacing,
                    step = spacing,
                    length = ceil(Int, (maximum(probe_x) - minimum(probe_x)) / spacing) + 5,
                ),
            )
            source_y = collect(
                range(
                    (floor(minimum(probe_y) / spacing) - 2) * spacing,
                    step = spacing,
                    length = ceil(Int, (maximum(probe_y) - minimum(probe_y)) / spacing) + 5,
                ),
            )
            Nx, Ny = length(source_x), length(source_y)
            x_matrix = repeat(source_x, 1, Ny)
            y_matrix = repeat(reshape(source_y, 1, Ny), Nx, 1)

            # Linear in the projected coordinates, so bilinear interpolation must reproduce it.
            analytic(x, y) = 1.0 + 2.0e-5 * x + 3.0e-5 * y

            # `MultiYearNORA3`'s default backend keeps ten time indices in memory, so the fixture
            # needs at least that many steps.
            dates = atmospheres.nora3_month_dates(2020, 1)[1:12]
            mkpath(atmosphere_directory(config))
            source_file = joinpath(atmosphere_directory(config), "NORA3_202001.nc")

            NCDataset(source_file, "c") do ds
                defDim(ds, "x", Nx)
                defDim(ds, "y", Ny)
                defDim(ds, "time", length(dates))
                defVar(ds, "x", Float64, ("x",))[:] = source_x
                defVar(ds, "y", Float64, ("y",))[:] = source_y
                defVar(ds, "time", dates, ("time",))
                defVar(ds, "projection_lambert", Int32, (); attrib = ["proj4" => proj4])
                for variable in atmospheres.ATMOSPHERE_VARIABLES
                    written = defVar(ds, variable.name, Float32, ("x", "y", "time"))
                    for step in 1:length(dates)
                        written[:, :, step] = Float32.(analytic.(x_matrix, y_matrix) .+ step)
                    end
                end
            end

            result = prepare_atmosphere(grid, config)
            @test result.output_file == atmosphere_path(config)
            @test length(result.times) == length(dates)
            @test result.times == dates
            @test sort(result.variables) == sort([v.name for v in atmospheres.ATMOSPHERE_VARIABLES])

            NCDataset(result.output_file) do ds
                # The layout the simulation-time reader requires: (longitude, latitude, time) in
                # Julia order, uniform 1D centers, CF time.
                @test NCDatasets.dimnames(ds["air_temperature_2m"]) == ("lon", "lat", "time")
                @test size(ds["air_temperature_2m"]) == (length(longitude), length(latitude), length(dates))
                @test all(isapprox(ds["lon"][2] - ds["lon"][1]), diff(ds["lon"][:]))
                @test all(isapprox(ds["lat"][2] - ds["lat"][1]), diff(ds["lat"][:]))
                @test ds["time"][:] == dates
                @test ds["air_temperature_2m"].attrib["units"] == "K"
                @test ds["swrad"].attrib["units"] == "W m-2"

                # Every variable is present and finite — an atmospheric field has no land mask to
                # excuse a NaN.
                for variable in atmospheres.ATMOSPHERE_VARIABLES
                    @test haskey(ds, variable.name)
                    @test all(isfinite, ds[variable.name][:, :, :])
                end

                # The regridded values match the analytic field at the prepared nodes.
                source = ProjectedAtmosphereGrid(source_x, source_y, proj4)
                node_x, node_y = atmospheres.projected_atmosphere_nodes(ds["lon"][:], ds["lat"][:], source)
                for step in (1, 7, 12)
                    expected = Float32.(analytic.(node_x, node_y) .+ step)
                    @test all(isapprox.(ds["swrad"][:, :, step], expected; rtol = 1e-5))
                end
            end

            # A second month cut to a different source window is refused rather than read with the
            # first month's coordinates — the hazard when `padding` changes and the download skips
            # months already present.
            february = atmospheres.nora3_month_dates(2020, 2)[1:3]
            NCDataset(joinpath(atmosphere_directory(config), "NORA3_202002.nc"), "c") do ds
                defDim(ds, "x", Nx - 1)
                defDim(ds, "y", Ny)
                defDim(ds, "time", length(february))
                defVar(ds, "x", Float64, ("x",))[:] = source_x[1:end-1]
                defVar(ds, "y", Float64, ("y",))[:] = source_y
                defVar(ds, "time", february, ("time",))
                defVar(ds, "projection_lambert", Int32, (); attrib = ["proj4" => proj4])
                for variable in atmospheres.ATMOSPHERE_VARIABLES
                    defVar(ds, variable.name, Float32, ("x", "y", "time"))[:, :, :] .= 1.0f0
                end
            end
            @test_throws ErrorException prepare_atmosphere(grid, config)
            rm(joinpath(atmosphere_directory(config), "NORA3_202002.nc"))

            # ...and the produced file satisfies the reader contract end to end.
            dataset = MultiYearNORA3(config)
            @test dataset.size == (length(longitude), length(latitude))
            @test dataset.all_dates == dates
            @test_nowarn NORA3PrescribedAtmosphere(CPU(), Float32; dataset)
            @test_nowarn NORA3PrescribedRadiation(CPU(), Float32; dataset)

            nora3 = FjordSim.Atmospheres.NORA3
            series(; kw...) = nora3.NORA3FieldTimeSeries(
                :temperature,
                CPU(),
                Float32;
                dataset,
                start_date = first(dates),
                end_date = last(dates),
                kw...,
            )
            record(step) = NCDataset(result.output_file) do ds
                Float32.(coalesce.(ds["air_temperature_2m"][:, :, step], NaN32))
            end

            full = series()
            # The time axis must carry the grid's float type. `native_times` gives `Float64` seconds,
            # and against a `Float32` grid that makes the interpolation weight — and so the
            # interpolated value — `Float64`, which boxes inside GPU kernels.
            @test eltype(full.times) === eltype(full.grid)
            @test full.times[1] == 0
            @test full.times[2] - full.times[1] == 3600
            @test interior(full, :, :, 1, 1) ≈ record(1)

            # A sub-range must read the records it is timestamped for. `time_indices` counts slots of
            # the series, which coincide with records of the file only when the series spans it, so
            # reading the file by slot shifted the whole atmosphere by the selection's offset.
            offset = 4
            shifted = nora3.NORA3FieldTimeSeries(
                :temperature,
                CPU(),
                Float32;
                dataset,
                start_date = dates[offset],
                end_date = last(dates),
            )
            @test interior(shifted, :, :, 1, 1) ≈ record(offset)
            @test interior(shifted, :, :, 1, 1) ≉ record(1)
            @test length(shifted.times) == length(dates) - offset + 1
            @test shifted.times[1] == 0

            # `reference_date` moves where t = 0 sits without touching which records are read, which
            # is what lets the atmosphere and the forcing share an origin.
            anchored = series(reference_date = first(dates) - Hour(3))
            @test anchored.times[1] == 3 * 3600
            @test interior(anchored, :, :, 1, 1) ≈ record(1)

            # A date the file does not carry is refused. Dropping it would leave the series one slot
            # short of the axis it was built for, so every later slot would read the wrong record.
            @test_throws ErrorException nora3.nora3_time_indices(
                dataset,
                [first(dates), first(dates) - Hour(1)],
                :temperature,
            )
        end
    end

    @testset "reader window" begin
        # A prepared file with fewer steps than the backend's default window. Written by hand rather
        # than through `prepare_atmosphere`, because the point is the reader, not the regrid.
        nora3 = atmospheres.NORA3

        mktempdir() do tmp
            config = test_atmosphere_config(data_root = tmp)
            dates = [DateTime(2020, 1, 1) + Hour(step - 1) for step in 1:3]
            longitude = collect(10.0:0.02:10.1)
            latitude = collect(59.0:0.02:59.1)

            NCDataset(atmosphere_path(config), "c") do ds
                defDim(ds, "lon", length(longitude))
                defDim(ds, "lat", length(latitude))
                defDim(ds, "time", length(dates))
                defVar(ds, "lon", Float64, ("lon",))[:] = longitude
                defVar(ds, "lat", Float64, ("lat",))[:] = latitude
                defVar(ds, "time", dates, ("time",))
                for variable in atmospheres.ATMOSPHERE_VARIABLES
                    written = defVar(ds, variable.name, Float32, ("lon", "lat", "time"))
                    for step in eachindex(dates)
                        written[:, :, step] = fill(Float32(step), length(longitude), length(latitude))
                    end
                end
            end

            dataset = MultiYearNORA3(config)
            @test dataset.size == (length(longitude), length(latitude))
            @test atmosphere_date_range(config) == (first(dates), last(dates))

            # The default backend asks for ten slots. Left unclamped that wraps `time_indices` past 1
            # more than once, which the wrap split cannot express: the first chunk comes out empty and
            # `copyto!` accepts the short read, leaving the tail of the window stale.
            fts = nora3.NORA3FieldTimeSeries(
                :temperature,
                CPU(),
                Float32;
                dataset,
                start_date = first(dates),
                end_date = last(dates),
            )
            @test length(fts.times) == length(dates)
            for step in eachindex(dates)
                @test all(isapprox(Float32(step)), interior(fts, :, :, 1, step))
            end
        end
    end
end

@testset "Simulations" begin
    @testset "config" begin
        simulation_forcing_path = FjordSim.Simulations.simulation_forcing_path
        checkpoints = FjordSim.Simulations.checkpoints
        loop_output_paths = FjordSim.Simulations.loop_output_paths
        validate_writers = FjordSim.Simulations.validate_writers

        # No field has a default, so omitting any is an UndefKeywordError rather than a
        # plausible-looking run. Every field is swept, not a hand-picked subset, so adding one to
        # `SimulationConfig` without a default is covered the moment `test_simulation_fields` names it.
        fields = test_simulation_fields()
        @test Set(keys(fields)) == Set(fieldnames(SimulationConfig))
        for name in fieldnames(SimulationConfig)
            @test_throws UndefKeywordError SimulationConfig(; delete!(copy(fields), name)...)
        end
        @test_throws UndefKeywordError SimulationConfig()

        # The same rule holds one level down, and that is where it now matters most: the model config
        # is where the scientific choices live, so a default there would be one fjord's physics
        # silently applied to another's.
        model_fields = test_model_fields()
        @test Set(keys(model_fields)) == Set(fieldnames(CoupledHydrostaticSimulation))
        for name in fieldnames(CoupledHydrostaticSimulation)
            @test_throws UndefKeywordError CoupledHydrostaticSimulation(;
                delete!(copy(model_fields), name)...
            )
        end

        # `Base.@kwdef` on a *parametric* struct dispatches on the declared field types instead of
        # converting. So a duration written as an integer is a MethodError, not a silent promotion —
        # which is why every setup writes `SimulationConfig`'s durations with `Oceananigans.Units`
        # constants, all of which are already Float64.
        @test_throws MethodError test_simulation_config(stop_time = 3600)
        @test test_simulation_config(stop_time = 1hour).stop_time === 3600.0
        # `loops` is the one field this does not apply to, being an Int already.
        @test test_simulation_config(loops = 2).loops === 2
        # The nested configs are the other exception, and deliberately so: each has a hand-written
        # keyword constructor that converts, so the rule a reader learns from `SimulationConfig` does
        # not silently mislead them about the four types they will write far more often.
        @test test_checkpoint_writer(interval = 3600).interval === 3600.0
        @test test_snapshot_writer(interval = 3600).interval === 3600.0
        @test test_time_stepping(cfl = 1//10).cfl === 0.1

        # A writer with no schedule, or with nothing to write, is rejected at construction — that is
        # a struct invariant, not something `build_simulation` should be discovering.
        @test_throws ArgumentError test_snapshot_writer(interval = 0.0)
        @test_throws ArgumentError test_snapshot_writer(variables = ())
        @test_throws ArgumentError test_checkpoint_writer(interval = 0.0)
        @test_throws ArgumentError test_progress_callback(interval = 0.0)

        # `extra_kwargs` is the escape hatch, and it can only *add*. A slot that reaches no
        # constructor, and a key that would win the splat over a keyword `coupled_simulation`
        # already passes, are both rejected at construction — where the setup file that made the
        # mistake is still what the error is about. Without the second check, a stray `tracers`
        # would silently disagree with `model_tracers`, which `build_simulation` already used to
        # pick the forcing terms and the open tracer boundaries before the model existed.
        @test_throws ArgumentError test_model_config(extra_kwargs = (ocean_mdoel = (;),))
        @test_throws ArgumentError test_model_config(
            extra_kwargs = (ocean_model = (tracers = (:T,),),),
        )
        @test_throws ArgumentError test_model_config(
            extra_kwargs = (ocean_model = (free_surface = nothing,),),
        )
        @test_throws ArgumentError test_model_config(extra_kwargs = (coupled_model = (radiation = nothing,),))
        @test_throws ArgumentError test_model_config(extra_kwargs = (ocean_simulation = (Δt = 1.0,),))
        # A slot holding something that is not a keyword set says so, rather than failing at the
        # splat inside `coupled_simulation` after the grid has been built.
        @test_throws ArgumentError test_model_config(extra_kwargs = (ocean_model = 1,))
        # An omitted slot is empty rather than missing, so `coupled_simulation` splats all four
        # unconditionally.
        empty_slots = test_model_config()
        for slot in keys(FjordSim.Simulations.EXTRA_KWARG_SLOTS)
            @test FjordSim.Simulations.extra_kwargs(empty_slots, slot) === (;)
        end
        # A keyword the built-in method does not already pass goes through untouched.
        timestepped = test_model_config(extra_kwargs = (ocean_model = (timestepper = :SplitRungeKutta3,),))
        @test FjordSim.Simulations.extra_kwargs(timestepped, :ocean_model) ==
              (timestepper = :SplitRungeKutta3,)

        # Two callbacks under one name replace each other in `simulation.callbacks`, so only the
        # last would ever fire — the same failure `validate_writers` rejects for writers.
        @test_throws ArgumentError FjordSim.Simulations.validate_callbacks(
            test_simulation_config(
                callbacks = (test_progress_callback(), test_progress_callback()),
            ),
        )
        @test isnothing(
            FjordSim.Simulations.validate_callbacks(test_simulation_config(callbacks = ())),
        )

        config = test_simulation_config()
        @test config isa AbstractSimulationConfig
        # Every field is concrete or parameterized — `isconcretetype(typeof(config))` would be vacuous,
        # since `typeof` never returns an abstract type.
        @test all(isconcretetype, fieldtypes(typeof(config)))
        @test all(isconcretetype, fieldtypes(typeof(test_snapshot_writer())))
        @test all(isconcretetype, fieldtypes(CheckpointWriter))
        @test all(isconcretetype, fieldtypes(AdaptiveTimeStep))
        # Six structural parameters, one per nested config. They grow only when a genuinely new
        # *kind* of nested config appears — `callbacks` was the sixth, replacing a bare
        # `progress_interval::Float64` that could only change how often the one hardcoded report
        # fired. A new subtype of an existing supertype is never a new field.
        @test length(typeof(config).parameters) == 6
        # Nine components plus `extra_kwargs`, which *is* the collapse the ninth was headed for:
        # rather than a field per remaining constructor keyword, one `NamedTuple` carries all of
        # them. Ten is the ceiling — anything further goes inside `extra_kwargs`.
        # `free_surface` counts here too: it is itself a config, not a bare `Float64`, dispatched
        # by its own `free_surface(config, grid)` hook.
        @test length(typeof(test_model_config()).parameters) == 10

        # The tracer list reaches `build_simulation` through a hook rather than a field read, so a
        # model config that names its tracers differently still composes.
        @test model_tracers(test_model_config()) == (:T, :S)

        # `free_surface` is dispatched, not a plain scalar: a fixture grid is enough to build the
        # object `coupled_simulation` would hand `HydrostaticFreeSurfaceModel`.
        free_surface_grid = LatitudeLongitudeGrid(CPU(), test_grid_config())
        @test free_surface(SplitExplicitFreeSurfaceConfig(cfl = 0.7), free_surface_grid) isa
              SplitExplicitFreeSurface

        # The device is a Symbol, not a live CPU()/GPU(), so a setup file loads without a GPU.
        @test fieldtype(typeof(config), :architecture) === Symbol

        # `run_tag` is the wall-clock launch instant, so every invocation names its own files, and moving
        # `start_date` — which used to be the tag — leaves it alone. It is constant for the life of the
        # process, since `results_path`, `log_path` and the checkpointer all have to agree.
        @test occursin(r"^\d{8}T\d{6}$", run_tag(config))
        @test run_tag(config) == FjordSim.Configs.LAUNCH_TAG[]
        @test run_tag(test_simulation_config(start_date = DateTime(1999, 12, 31))) ==
              FjordSim.Configs.LAUNCH_TAG[]

        # The tag is inserted before the extension either way, and the loop index goes on top of it so a
        # looped run gives each repetition its own file. The filename is the writer's; only the root
        # is the simulation config's.
        writer = test_snapshot_writer()
        with_run_tag() do
            @test results_path(writer, config) ==
                  joinpath(TEST_RESULTS_ROOT, "snapshots_test_$PINNED_RUN_TAG.nc")
            @test results_path(writer, config, 7) ==
                  joinpath(TEST_RESULTS_ROOT, "snapshots_test_$(PINNED_RUN_TAG)_loop07.nc")

            # An absolute path relocates just that file.
            relocated = test_snapshot_writer(output_file = "/elsewhere/run.nc")
            @test results_path(relocated, config) == "/elsewhere/run_$PINNED_RUN_TAG.nc"
        end

        # A single run keeps the shorter name; only a looped one gains the index.
        @test FjordSim.Simulations.loop_output_path(writer, config, 1) ==
              results_path(writer, config)
        looped = test_simulation_config(loops = 3)
        @test FjordSim.Simulations.loop_output_path(writer, looped, 1) ==
              results_path(writer, looped, 1)

        # Which writers checkpoint, and which name a file the run should report, are traits rather
        # than `isa` tests — so a setup naming no `CheckpointWriter` writes none, and the reported
        # list is the snapshots alone even though the checkpointer also writes files.
        @test !checkpoints(test_snapshot_writer())
        @test checkpoints(test_checkpoint_writer())
        @test checkpoints(config)
        snapshots_only = test_simulation_config(writers = (test_snapshot_writer(),))
        @test !checkpoints(snapshots_only)
        @test !checkpoints(test_simulation_config(writers = ()))
        with_run_tag() do
            @test loop_output_paths(config, 1) == [results_path(writer, config)]
            @test loop_output_paths(test_simulation_config(writers = ()), 1) == String[]
        end

        # Two checkpointing writers, or two writers under one name, are rejected before anything is
        # read: `run!` cannot tell which checkpoint to resume from, and the second writer under a
        # name silently replaces the first so one of the requested files is never written.
        @test isnothing(validate_writers(config))
        @test_throws ArgumentError validate_writers(test_simulation_config(
            writers = (test_checkpoint_writer(), test_checkpoint_writer()),
        ))
        @test_throws ArgumentError validate_writers(test_simulation_config(
            writers = (test_snapshot_writer(), test_snapshot_writer(output_file = "other.nc")),
        ))

        # The window every loop replays, which is what the prepare steps pad their axes to. It does not
        # grow with `loops`, since each repetition replays the same interval.
        window = test_simulation_config(start_date = DateTime(2020, 1, 1), stop_time = 365days)
        @test coverage_window(window) == (DateTime(2020, 1, 1), DateTime(2020, 12, 31))
        @test coverage_window(test_simulation_config(
            start_date = DateTime(2020, 1, 1),
            stop_time = 365days,
            loops = 5,
        )) == coverage_window(window)
        @test isnothing(coverage_window(nothing))     # a setup naming no simulation config

        # The checkpoint prefix carries the loop index because the checkpoint itself does not: the state
        # records the clock but not which repetition produced it. It carries no run tag, so a later launch
        # can name — and so resume — the checkpoints of the one before it.
        @test FjordSim.Simulations.checkpoint_prefix(3) == "checkpoint_loop03"

        # Architecture resolution shares `interpolation_architecture`'s Val methods, so :cpu pins the
        # CPU and an unknown selector is rejected rather than silently defaulted.
        @test simulation_architecture(test_simulation_config(architecture = :cpu)) isa CPU
        @test_throws ArgumentError simulation_architecture(
            test_simulation_config(architecture = :tpu),
        )

        # Which edge the domain is open on is stated once, on the open-boundary data config, and read
        # from there by everything that acts on it — see the "open lateral boundary" testset for the
        # per-edge dispatch. It is a `Configs` accessor rather than a pair of one-liners reaching
        # through the forcing config, which is what a `boundary_config` independent of forcing needs.
        @test open_edges(test_boundaries_config(open_edges = :west)) == [:west]
        @test isempty(open_edges(nothing))

        # A setup that configures rivers simulates the rivers-augmented copy, never the pre-rivers
        # file. Both registered setups do; the no-rivers branch is covered by `no_rivers` below.
        oslo = oslofjorden()
        @test simulation_forcing_path(oslo) == river_forcing_path(oslo.forcing_config.rivers)
        drammen = drammensfjorden()
        @test simulation_forcing_path(drammen) == river_forcing_path(drammen.forcing_config.rivers)
        @test simulation_forcing_path(drammen) != forcing_path(drammen.forcing_config)

        # Every prerequisite is reported as the command that produces it, and every check runs before
        # anything is read or allocated, so the temporary root stays empty apart from the stub.
        mktempdir() do tmp
            bathymetry_config = test_bathymetry_config(data_root = tmp)
            fjord(forcing_config) = FjordConfig(;
                grid_config = test_grid_config(),
                bathymetry_config,
                forcing_config,
                simulation_config = test_simulation_config(results_root = tmp),
            )

            # `rivers` is a type parameter of NorKystConfig, so whether a config has any is a
            # construction-time choice — hence two configs rather than one mutated.
            with_rivers = fjord(test_forcing_config(
                data_root = tmp,
                rivers = OF800RiversConfig(data_root = tmp),
            ))
            without_rivers = fjord(test_forcing_config(data_root = tmp))

            @test_throws "prepare_bathymetry" build_simulation(with_rivers)

            touch(bathymetry_path(bathymetry_config))
            # A config naming rivers is told to run add_rivers, since that is what writes the file
            # `simulation_forcing_path` picked; one naming none is told to run prepare_forcing.
            @test_throws "add_rivers" build_simulation(with_rivers)
            @test isnothing(without_rivers.forcing_config.rivers)
            @test simulation_forcing_path(without_rivers) ==
                  forcing_path(without_rivers.forcing_config)
            @test_throws "prepare_forcing" build_simulation(without_rivers)
        end
    end

    @testset "time coverage" begin
        validate_time_coverage = FjordSim.Simulations.validate_time_coverage
        covered = (DateTime(2020, 1, 1, 12), DateTime(2020, 12, 31, 12))

        # Both readers use `Cyclical()` time indexing, which wraps rather than failing outside its
        # data, so a run that outlasts its forcing would quietly replay the beginning.
        @test isnothing(validate_time_coverage("forcing", covered, DateTime(2020, 1, 1, 12), 86400.0, "x"))
        @test isnothing(validate_time_coverage("forcing", covered, DateTime(2020, 6, 1), 86400.0, "x"))
        # A source that cannot report its dates is skipped rather than blocking the run.
        @test isnothing(validate_time_coverage("atmosphere", nothing, DateTime(2019, 1, 1), 1e9, "x"))

        @test_throws ErrorException validate_time_coverage(
            "forcing", covered, DateTime(2020, 1, 1), 86400.0, "x",   # starts before the data
        )
        @test_throws ErrorException validate_time_coverage(
            "forcing", covered, DateTime(2020, 12, 31), 86400.0, "x", # runs past the data
        )

        mktempdir() do tmp
            filepath = joinpath(tmp, "forcing.nc")
            dates = [DateTime(2020, 1, 1, 12) + Hour(24 * (step - 1)) for step in 1:5]
            NCDataset(filepath, "c") do ds
                defDim(ds, "time", length(dates))
                defVar(ds, "time", dates, ("time",))
            end
            # A hook on the forcing config now, the mirror of `atmosphere_date_range`, rather
            # than a bare read of a path dispatched on nothing — so a source whose prepared files
            # are not this NetCDF layout is not validated by a reader it does not use.
            forcing_config = test_forcing_config(data_root = tmp)
            @test forcing_date_range(forcing_config, filepath) == (first(dates), last(dates))
            @test isnothing(forcing_date_range(nothing, filepath))
        end
    end

    @testset "time axis padding" begin
        pad_time_steps = FjordSim.Forcing.pad_time_steps
        pad_atmosphere_records = FjordSim.Atmospheres.pad_atmosphere_records
        ForcingTimeStep = FjordSim.Forcing.ForcingTimeStep
        SourceRecord = FjordSim.Forcing.SourceRecord

        step(date, index) =
            ForcingTimeStep(date, SourceRecord(date, "a.nc", index), SourceRecord(date, "a.nc", index), 0.0f0)
        steps = [step(DateTime(2020, 1, 1, 12), 1), step(DateTime(2020, 1, 2, 12), 2)]

        # A setup naming no simulation config prepares exactly the downloaded range.
        @test pad_time_steps(steps, nothing) === steps
        # So does a window the axis already spans.
        @test length(pad_time_steps(steps, (DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)))) == 2

        # A window reaching past either end gets one replicated record there, carrying the same source
        # records and blend weight, so the written field is identical to its neighbour.
        padded = pad_time_steps(steps, (DateTime(2020, 1, 1), DateTime(2020, 1, 3)))
        @test [entry.date for entry in padded] == [
            DateTime(2020, 1, 1),
            DateTime(2020, 1, 1, 12),
            DateTime(2020, 1, 2, 12),
            DateTime(2020, 1, 3),
        ]
        @test first(padded).lower.index == 1 && first(padded).upper.index == 1
        @test last(padded).lower.index == 2 && last(padded).weight == 0.0f0

        # A pad may reach at most one record spacing. Unbounded it would manufacture the very coverage
        # `validate_time_coverage` exists to verify: one replicated December would "cover" a window a
        # year past the data, and the run would interpolate twelve months between identical records.
        @test_throws ErrorException pad_time_steps(
            steps, (DateTime(2019, 12, 31, 11), DateTime(2020, 1, 2, 12)),
        )
        @test_throws ErrorException pad_time_steps(
            steps, (DateTime(2020, 1, 1, 12), DateTime(2020, 1, 3, 13)),
        )
        # Exactly one spacing is allowed at each end.
        @test length(pad_time_steps(steps, (DateTime(2019, 12, 31, 12), DateTime(2020, 1, 3, 12)))) == 4
        # A one-record axis has no spacing to bound a pad by, and says so instead of guessing.
        @test_throws ErrorException pad_time_steps(
            [first(steps)], (DateTime(2020, 1, 1), DateTime(2020, 1, 1, 12)),
        )

        records = [
            AtmosphereRecord(DateTime(2020, 1, 1, 0), "n.nc", 1),
            AtmosphereRecord(DateTime(2020, 1, 1, 1), "n.nc", 2),
        ]
        @test pad_atmosphere_records(records, nothing) === records

        # A pad points at the same file and index, so `validate_source_consistency` sees no new file.
        padded = pad_atmosphere_records(records, (DateTime(2019, 12, 31, 23), DateTime(2020, 1, 1, 2)))
        @test [entry.date for entry in padded] == [
            DateTime(2019, 12, 31, 23),
            DateTime(2020, 1, 1, 0),
            DateTime(2020, 1, 1, 1),
            DateTime(2020, 1, 1, 2),
        ]
        @test [entry.index for entry in padded] == [1, 1, 2, 2]
        @test all(entry -> entry.filepath == "n.nc", padded)

        @test_throws ErrorException pad_atmosphere_records(
            records, (DateTime(2019, 12, 31, 22), DateTime(2020, 1, 1, 1)),
        )
        @test_throws ErrorException pad_atmosphere_records(
            records, (DateTime(2020, 1, 1, 0), DateTime(2020, 1, 1, 3)),
        )
    end

    @testset "loop restart clocks" begin
        rewind_clock! = FjordSim.Simulations.rewind_clock!

        # NumericalEarth's own `reset_clock!(::EarthSystemModel)` cannot be used: its per-component
        # fallback is `reset!(getproperty(component, :clock))` and `components` includes `sea_ice`, so a
        # `FreezingLimitedOceanTemperature` — a liquidus and nothing else — makes it throw. This is the
        # regression test for that, and the reason `rewind_clock!` dispatches on having a clock at all.
        @test_throws Exception Oceananigans.Simulations.reset_clock!(FreezingLimitedOceanTemperature())
        @test isnothing(rewind_clock!(FreezingLimitedOceanTemperature()))
        @test isnothing(rewind_clock!(nothing))

        # A component that does have a clock is sent back to zero, state untouched.
        clock = Oceananigans.TimeSteppers.Clock{Float64}(time = 1234.5, iteration = 42)
        rewind_clock!(ClockHolder(clock))
        @test clock.time == 0
        @test clock.iteration == 0
    end

    @testset "initial conditions" begin
        resolve_initial_conditions = FjordSim.Simulations.resolve_initial_conditions
        initial_conditions_date = FjordSim.Simulations.initial_conditions_date
        forcing_state = FjordSim.Simulations.forcing_state
        results_state = FjordSim.Simulations.results_state

        # An unnamed date reads the run's own start_date; a named one is taken as given.
        @test initial_conditions_date(FromForcing(), DateTime(2020, 1, 1, 12)) == DateTime(2020, 1, 1, 12)
        @test initial_conditions_date(FromForcing(DateTime(2020, 6, 1)), DateTime(2020, 1, 1, 12)) ==
              DateTime(2020, 6, 1)

        # A literal NamedTuple is already what `set!` wants, so it passes through by identity — the
        # behaviour every setup had before the other two shapes existed.
        constants = (T = 5.0, S = 33.0)
        @test resolve_initial_conditions(constants, nothing, "unused.nc", nothing) === constants

        Nx, Ny, Nz = 6, 5, 4

        mktempdir() do tmp
            bathymetry_config = test_bathymetry_config(data_root = tmp)
            (; grid_config, grid) = immersed_test_grid(
                bathymetry_path(bathymetry_config);
                size = (Nx, Ny, Nz),
                halo = (3, 3, 3),
            )

            # A prepared forcing file: land is NaN, and every variable has a `_lambda` twin the state
            # reader must ignore. Each record is filled with its own index, so a state read from the
            # wrong record is visible.
            dates = [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)]
            forcing_file = write_prepared_forcing(
                joinpath(tmp, "forcing.nc");
                size = (Nx, Ny, Nz),
                dates,
                value = (name, index) -> index,
                land = (1, 1, 1),
            )

            state = forcing_state(forcing_file, grid, DateTime(2020, 1, 2, 12), (:T, :S))
            @test keys(state) == (:T, :S, :u, :v)
            @test size(state.T) == (Nx, Ny, Nz)
            @test size(state.u) == (Nx + 1, Ny, Nz)     # u is on x faces
            @test size(state.v) == (Nx, Ny + 1, Nz)     # v is on y faces
            # Converted to the grid's element type, since a Union{Missing,Float32} array cannot be
            # moved to a GPU and a mismatched element type need not convert on the device.
            @test all(array -> eltype(array) === eltype(grid), values(state))
            @test all(array -> all(isfinite, array), values(state))
            @test state.T[1, 1, 1] == 0                 # land zeroed
            @test state.T[2, 2, 2] == 2                 # the second record, not the first

            # The date must be on the axis: a nearest-record fallback would silently start the run
            # somewhere else.
            @test_throws ErrorException forcing_state(forcing_file, grid, DateTime(2020, 1, 3, 12), (:T, :S))

            # A results file: time is seconds from the writing run's model zero, so a calendar date
            # needs the `start_date` attribute `build_simulation` records.
            function write_results(filepath, attributes)
                NCDataset(filepath, "c"; attrib = attributes) do ds
                    defDim(ds, "λ_caa", Nx)
                    defDim(ds, "φ_aca", Ny)
                    defDim(ds, "z_aac", Nz)
                    defDim(ds, "λ_faa", Nx + 1)
                    defDim(ds, "φ_afa", Ny + 1)
                    defDim(ds, "time", 3)
                    defVar(ds, "time", [0.0, 3600.0, 7200.0], ("time",))
                    for (name, xdim, ydim, nx, ny) in (
                        ("T", "λ_caa", "φ_aca", Nx, Ny),
                        ("S", "λ_caa", "φ_aca", Nx, Ny),
                        ("u", "λ_faa", "φ_aca", Nx + 1, Ny),
                        ("v", "λ_caa", "φ_afa", Nx, Ny + 1),
                    )
                        defVar(ds, name, Float32, (xdim, ydim, "z_aac", "time"))
                        for index in 1:3
                            ds[name][:, :, :, index] = fill(Float32(10 * index), (nx, ny, Nz))
                        end
                    end
                end

                return filepath
            end

            tagged = joinpath(tmp, "snapshots_tagged.nc")
            write_results(tagged, ["start_date" => string(DateTime(2020, 1, 1, 12))])

            # No date takes the last record, which is the only thing an untagged file can offer.
            untagged = joinpath(tmp, "snapshots_untagged.nc")
            write_results(untagged, Pair{String,String}[])
            @test results_state(untagged, grid, nothing, (:T, :S)).T[1, 1, 1] == 30
            @test_throws ErrorException results_state(untagged, grid, DateTime(2020, 1, 1, 13), (:T, :S))

            # With the attribute, a date resolves to its record.
            @test results_state(tagged, grid, DateTime(2020, 1, 1, 13), (:T, :S)).T[1, 1, 1] == 20
            @test results_state(tagged, grid, DateTime(2020, 1, 1, 12), (:T, :S)).T[1, 1, 1] == 10
            @test_throws ErrorException results_state(tagged, grid, DateTime(2020, 6, 1), (:T, :S))

            # A state file from another grid is a clear error rather than a silent misread.
            wrong = joinpath(tmp, "wrong_grid.nc")
            NCDataset(wrong, "c") do ds
                defDim(ds, "λ_caa", Nx + 1)
                defDim(ds, "φ_aca", Ny)
                defDim(ds, "z_aac", Nz)
                defDim(ds, "time", 1)
                defVar(ds, "time", [0.0], ("time",))
            end
            @test_throws DimensionMismatch results_state(wrong, grid, nothing, (:T, :S))

            # `FromResults` resolves a relative path against `results_root`, like `output_file` does.
            simulation_config = test_simulation_config(results_root = tmp)
            @test resolve_initial_conditions(
                FromResults("snapshots_tagged.nc"),
                grid,
                forcing_file,
                simulation_config,
            ).T[1, 1, 1] == 30
            @test_throws ErrorException resolve_initial_conditions(
                FromResults("absent.nc"),
                grid,
                forcing_file,
                simulation_config,
            )

            # And `FromForcing` with no date reads the run's own start_date, which the fixture puts at
            # noon to match the axis `write_prepared_forcing` writes.
            @test resolve_initial_conditions(
                FromForcing(),
                grid,
                forcing_file,
                simulation_config,
            ).T[2, 2, 2] == 1
        end
    end

    @testset "checkpoint resume" begin
        resume_loop = FjordSim.Simulations.resume_loop
        checkpoint_prefix = FjordSim.Simulations.checkpoint_prefix

        mktempdir() do tmp
            resuming(; kwargs...) =
                test_simulation_config(; results_root = tmp, loops = 4, pickup = true, kwargs...)

            # Without `pickup` the run always starts at the first repetition.
            @test resume_loop(resuming(pickup = false)) == 1

            # With `pickup` but no checkpoints, it warns and starts from the beginning rather than
            # failing — `run!` itself tolerates a missing checkpoint the same way.
            @test (@test_logs (:warn,) resume_loop(resuming())) == 1

            # The loop index is read back out of the checkpoint filename, because the checkpoint state
            # records the clock but not which repetition produced it.
            touch(joinpath(tmp, checkpoint_prefix(1) * "_iteration50.jld2"))
            touch(joinpath(tmp, checkpoint_prefix(3) * "_iteration7.jld2"))
            touch(joinpath(tmp, checkpoint_prefix(3) * "_iteration90.jld2"))
            @test resume_loop(resuming()) == 3

            # Checkpoints carry no run tag, so they are resumed whichever launch wrote them and whatever
            # `start_date` it ran — there is one resumable run per results directory. That is what keeps
            # `pickup` working now that the snapshots are named per launch.
            with_run_tag() do
                @test resume_loop(resuming(start_date = DateTime(2021, 5, 6, 7))) == 3
            end

            # `pickup` with no checkpointing writer is a config contradicting itself — nothing this
            # run writes could ever be resumed from — so it is rejected rather than left to fail
            # inside `run!`. Only expressible now that checkpointing is a writer, not a threshold.
            @test_throws ArgumentError resume_loop(
                resuming(writers = (test_snapshot_writer(),)),
            )
        end
    end

    @testset "build simulation" begin
        # The "Simulation config" testset only reaches build_simulation's two prerequisite errors, so
        # everything past them — the HydrostaticFreeSurfaceModel, the coupled OceanSeaIceModel and the
        # tracer top boundary conditions NumericalEarth reads its freshwater exchange out of — is
        # otherwise never executed. Synthetic inputs keep it on the CPU with no ~/FjordSim_data.
        mktempdir() do tmp
            # Small, but the halo is a real setup's: WENO plus a biharmonic closure sets the floor, and
            # an ImmersedBoundaryGrid needs one point more than that in every direction.
            Nx, Ny, Nz = 16, 16, 8
            bathymetry_config = test_bathymetry_config(data_root = tmp)
            (; grid_config) = immersed_test_grid(
                bathymetry_path(bathymetry_config);
                size = (Nx, Ny, Nz),
                halo = (7, 7, 7),
            )

            # The open lateral boundary reads its exterior state from the prepared boundary file, so
            # this build needs one — the file is what makes the boundary open rather than a wall.
            boundaries_config = test_boundaries_config(data_root = tmp)
            forcing_config = test_forcing_config(data_root = tmp)
            write_prepared_forcing(forcing_path(forcing_config); size = (Nx, Ny, Nz))
            write_prepared_boundaries(
                boundary_data_path(boundaries_config);
                size = (Nx, Ny, Nz),
                dates = [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)],
            )

            # `atmosphere_config` stays unnamed, so both prescribed_atmosphere hooks yield nothing and the
            # coupled model runs ocean-only — which is what keeps this off the network and off a GPU.
            #
            # Every build below varies only String/Float64/Int/Bool fields, so all of them share one
            # `SimulationConfig` instantiation and one `build_simulation` specialization. Overriding
            # any of the five nested-config fields here — `writers` above all — would change the
            # config's type and recompile the whole coupled model.
            fjord(; kwargs...) = FjordConfig(;
                grid_config,
                bathymetry_config,
                forcing_config,
                boundary_config = boundaries_config,
                simulation_config = test_simulation_config(; results_root = tmp, kwargs...),
            )

            config = fjord()
            simulation_config = config.simulation_config
            snapshot, checkpointer = simulation_config.writers
            @test isnothing(config.atmosphere_config)
            @test isnothing(config.forcing_config.rivers)

            simulation = build_simulation(config)
            @test simulation isa Simulation
            @test simulation.stop_time == simulation_config.stop_time

            ocean_model = simulation.model.ocean.model
            @test Set(keys(ocean_model.tracers)) ⊇ Set((:T, :S))

            # The coupled model must have accepted (ocean, sea_ice) in that order — the flipped call
            # throws rather than mis-assembling, so reaching here at all pins the argument order.
            @test simulation.model.ocean isa Simulation
            @test !isnothing(simulation.model.sea_ice)

            # And the model's own tracer top conditions still expose the freshwater exchange after
            # passing through HydrostaticFreeSurfaceModel's boundary-condition materialization.
            S_top = ocean_model.tracers.S.boundary_conditions.top
            @test NumericalEarth.Oceans.extract_freshwater_flux(S_top.condition) isa AbstractField

            # The progress callback is attached under the name its config names, rather than the
            # fixed `:progress` key `build_simulation` used to hardcode along with the function and
            # the schedule type.
            @test haskey(simulation.callbacks, :progress)

            # `extra_kwargs` reaches the constructors it names. `verbose` is a stored `Simulation`
            # field and compiles no kernels, so this proves the splat arrives without paying for a
            # second time-stepper specialization. Only one such case: overriding `model` changes
            # `SimulationConfig`'s type and forces a fresh `build_simulation` specialization.
            passthrough = build_simulation(
                fjord(
                    model = test_model_config(
                        extra_kwargs = (
                            ocean_simulation = (verbose = false,),
                            coupled_simulation = (verbose = false,),
                        ),
                    ),
                ),
            )
            @test passthrough.verbose == false
            @test passthrough.model.ocean.verbose == false

            # The open boundary landed on the velocity normal to the boundary config's open edge, and
            # it radiates rather than being a wall — the whole point of the change, and something
            # only a real model build can show.
            edge = only(open_edges(config.boundary_config))
            normal_velocity = edge in (:south, :north) ? :v : :u
            edge_bc = getproperty(
                getproperty(ocean_model.velocities, normal_velocity).boundary_conditions, edge,
            )
            @test edge_bc.classification isa
                  Oceananigans.BoundaryConditions.NormalFlow{<:Oceananigans.BoundaryConditions.NormalRadiation}

            # The barotropic transport carries the Flather condition, and Oceananigans paired the
            # Chapman condition onto η by itself — which is why nothing in FjordSim states it.
            barotropic = edge in (:south, :north) ? :V : :U
            free_surface = ocean_model.free_surface
            barotropic_bc = getproperty(
                getproperty(free_surface.barotropic_velocities, barotropic).boundary_conditions, edge,
            )
            @test barotropic_bc.classification isa
                  Oceananigans.BoundaryConditions.NormalFlow{<:Oceananigans.BoundaryConditions.GravityWaveRadiation}
            @test getproperty(free_surface.displacement.boundary_conditions, edge).classification isa
                  Oceananigans.BoundaryConditions.Value{<:Oceananigans.BoundaryConditions.SurfaceWaveRadiation}

            # The tracer open boundaries landed on the same edge, on T and S, with a real
            # NormalRadiation scheme — proving the wiring reaches a real model, not just the
            # isolated `boundary_condition_sides` function.
            for name in (:T, :S)
                tracer_edge_bc =
                    getproperty(getproperty(ocean_model.tracers, name).boundary_conditions, edge)
                @test tracer_edge_bc.classification.scheme isa
                      Oceananigans.BoundaryConditions.NormalRadiation
            end

            # The snapshot writer lands under its own name, writes exactly the variables the config
            # asked for — the check that closes the gap between "the writer is attached" and
            # "`snapshot_outputs` resolved those Symbols to the right model fields" — and records
            # `start_date`, which is the only thing that lets a later `FromResults(path, date)` turn
            # a date into a record, a snapshot's own time axis being seconds from model zero.
            #
            # As a `Set` of `Symbol`, because `NetCDFWriter` re-keys the `NamedTuple` it is given
            # into a `Dict{String}`: neither the order nor the key type survives.
            ocean_writer = simulation.model.ocean.output_writers[snapshot.name]
            @test ocean_writer.filepath == results_path(snapshot, simulation_config)
            @test Set(Symbol.(keys(ocean_writer.outputs))) == Set(snapshot.variables)
            @test ocean_writer.global_attributes["start_date"] == string(simulation_config.start_date)

            # A variable the model does not have is an error naming it, not a silently dropped
            # column discovered a run later.
            @test_throws ArgumentError FjordSim.Simulations.snapshot_outputs(
                test_snapshot_writer(variables = (:T, :not_a_field)),
                ocean_model,
            )

            # The checkpointer goes on the *coupled* simulation: `prognostic_state` of the coupled model
            # is what a resumable state is, and `run!(…; pickup)` looks for it there. Exactly one, which
            # is what `checkpoint_path` requires.
            checkpointers =
                filter(writer -> writer isa Checkpointer, collect(values(simulation.output_writers)))
            @test length(checkpointers) == 1
            @test !any(
                writer -> writer isa Checkpointer,
                values(simulation.model.ocean.output_writers),
            )

            # Naming no `CheckpointWriter` attaches none at all. Asserted on the config rather than
            # by building a second simulation: dropping the writer changes `SimulationConfig`'s type
            # and would recompile the whole coupled model for one boolean.
            @test !FjordSim.Simulations.checkpoints(
                test_simulation_config(results_root = tmp, writers = (snapshot,)),
            )
            # ...and `attach_writers!` over that tuple leaves the coupled simulation's writers alone,
            # which is the behavioral half of the same claim, on the simulation already built.
            empty!(simulation.output_writers)
            FjordSim.Simulations.attach_writers!(
                simulation,
                test_simulation_config(results_root = tmp, writers = (snapshot,)),
                1,
            )
            @test !any(writer -> writer isa Checkpointer, values(simulation.output_writers))
            # Put the checkpointer back, since later assertions in this testset expect it.
            attach_writer!(simulation, checkpointer, simulation_config, 1)

            # Looping does not widen what the prepared files must span, so a two-loop run of the same
            # window still builds; the loop index reaches the output filename.
            looped_config = fjord(loops = 2)
            looped = build_simulation(looped_config)
            @test looped.stop_time == looped_config.simulation_config.stop_time
            @test looped.model.ocean.output_writers[snapshot.name].filepath ==
                  results_path(snapshot, looped_config.simulation_config, 1)

            # Restarting a loop over the real object graph: every clock goes back to zero, the ocean
            # state does not, the writer is replaced with the next loop's file, and the ocean
            # sub-simulation is marked uninitialized so its fresh writer gets a schedule and a t = 0
            # record. Exercised without a run!, since the loop's correctness lives here rather than in
            # the time stepping.
            looped.model.clock.time = 42.0
            looped.model.clock.iteration = 7
            looped.model.ocean.model.clock.time = 42.0
            looped.model.ocean.model.clock.iteration = 7
            set!(looped.model.ocean.model, T = 9.0)
            FjordSim.Simulations.restart_loop!(looped, looped_config.simulation_config, 2)

            @test looped.model.clock.time == 0
            @test looped.model.clock.iteration == 0
            @test looped.model.ocean.model.clock.time == 0
            @test looped.model.ocean.model.clock.iteration == 0
            @test !looped.model.ocean.initialized
            @test looped.model.ocean.output_writers[snapshot.name].filepath ==
                  results_path(snapshot, looped_config.simulation_config, 2)
            # The state carried over — that is the whole point of resetting only the clocks.
            @test maximum(interior(looped.model.ocean.model.tracers.T)) ≈ 9.0
            # stop_time is untouched, so the next repetition runs the same window.
            @test looped.stop_time == looped_config.simulation_config.stop_time

            # A loop count below one is rejected before anything is read or allocated.
            @test_throws ArgumentError build_simulation(fjord(loops = 0))

            # Building into a `results_root` that does not exist yet creates it. Nothing else would:
            # `NetCDFWriter` creates only the `dir` keyword it is not given, and a setup naming no
            # `CheckpointWriter` has nothing else that mkpaths. The old model-assembly function
            # did it; `coupled_simulation` is about assembling a model, so `build_simulation` does.
            fresh_root = joinpath(tmp, "not_created_yet", "nested")
            @test !isdir(fresh_root)
            @test build_simulation(fjord(results_root = fresh_root)) isa Simulation
            @test isdir(fresh_root)
        end
    end

    @testset "build simulation open all round" begin
        # A coupled model on a domain open on all four sides, which is what an open-ocean region needs
        # and what a single `open_edge` Symbol made impossible. Only a real build shows it: the four
        # edges' conditions have to merge without colliding, `HydrostaticFreeSurfaceModel` has to
        # regularize a `u` carrying a normal condition on two sides and a tangential one on the other
        # two, and Oceananigans has to pair its Chapman condition onto all four.
        mktempdir() do tmp
            Nx, Ny, Nz = 16, 16, 8
            edges = (:south, :north, :west, :east)
            bathymetry_config = test_bathymetry_config(data_root = tmp)
            (; grid_config) = immersed_test_grid(
                bathymetry_path(bathymetry_config);
                size = (Nx, Ny, Nz),
                halo = (7, 7, 7),
            )

            boundaries_config = test_boundaries_config(data_root = tmp, open_edges = edges)
            write_prepared_boundaries(
                boundary_data_path(boundaries_config);
                size = (Nx, Ny, Nz),
                edges,
                dates = [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)],
            )

            config = FjordConfig(;
                grid_config,
                bathymetry_config,
                boundary_config = boundaries_config,
                simulation_config = test_simulation_config(results_root = tmp),
            )

            simulation = build_simulation(config)
            @test simulation isa Simulation

            ocean_model = simulation.model.ocean.model
            free_surface = ocean_model.free_surface

            # Every side of both velocities is open, and each is the right kind: normal where the
            # component crosses the edge, tangential where it runs along it.
            for (component, normal_sides, tangential_sides) in (
                (:u, (:west, :east), (:south, :north)),
                (:v, (:south, :north), (:west, :east)),
            )
                conditions = getproperty(ocean_model.velocities, component).boundary_conditions
                for side in normal_sides
                    @test getproperty(conditions, side).classification isa
                          Oceananigans.BoundaryConditions.NormalFlow{
                        <:Oceananigans.BoundaryConditions.NormalRadiation,
                    }
                end
                for side in tangential_sides
                    @test getproperty(conditions, side).classification.scheme isa
                          Oceananigans.BoundaryConditions.NormalRadiation
                end
            end

            # Both barotropic transports carry the Flather condition on the two sides they are normal
            # to, and Oceananigans paired the Chapman condition onto all four sides of η by itself.
            for (barotropic, sides) in ((:U, (:west, :east)), (:V, (:south, :north)))
                conditions =
                    getproperty(free_surface.barotropic_velocities, barotropic).boundary_conditions
                for side in sides
                    @test getproperty(conditions, side).classification isa
                          Oceananigans.BoundaryConditions.NormalFlow{
                        <:Oceananigans.BoundaryConditions.GravityWaveRadiation,
                    }
                end
            end
            for side in edges
                @test getproperty(free_surface.displacement.boundary_conditions, side).classification.scheme isa
                      Oceananigans.BoundaryConditions.SurfaceWaveRadiation
            end

            # ...and every tracer is open on all four sides.
            for tracer in (:T, :S)
                conditions = getproperty(ocean_model.tracers, tracer).boundary_conditions
                for side in edges
                    @test getproperty(conditions, side).classification.scheme isa
                          Oceananigans.BoundaryConditions.NormalRadiation
                end
            end
        end
    end

    @testset "build simulation without forcing" begin
        # A run with no interior forcing at all, but a data-driven open boundary — the combination
        # that was impossible while the boundary dataset and the open edge both hung off the forcing
        # config, since a `nothing` forcing config took them down with it. `examples/oslofjorden.jl`
        # is the real-world shape of this.
        mktempdir() do tmp
            Nx, Ny, Nz = 16, 16, 8
            bathymetry_config = test_bathymetry_config(data_root = tmp)
            (; grid_config) = immersed_test_grid(
                bathymetry_path(bathymetry_config);
                size = (Nx, Ny, Nz),
                halo = (7, 7, 7),
            )

            boundaries_config = test_boundaries_config(data_root = tmp)
            write_prepared_boundaries(
                boundary_data_path(boundaries_config);
                size = (Nx, Ny, Nz),
                dates = [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)],
            )

            config = FjordConfig(;
                grid_config,
                bathymetry_config,
                boundary_config = boundaries_config,
                simulation_config = test_simulation_config(results_root = tmp),
            )
            @test isnothing(config.forcing_config)

            # No prepared forcing is looked for, so no prerequisite error fires and no forcing terms
            # are built — the model gets Oceananigans' own empty forcing.
            @test isnothing(FjordSim.Simulations.resolve_forcing_file(config))
            @test simulation_forcing(config.forcing_config, nothing, nothing, (:T, :S), nothing) ==
                  NamedTuple()

            simulation = build_simulation(config)
            @test simulation isa Simulation

            # ...and the open boundary is still there, on the edge the boundary config names.
            ocean_model = simulation.model.ocean.model
            @test ocean_model.velocities.v.boundary_conditions.south.classification isa
                  Oceananigans.BoundaryConditions.NormalFlow{
                <:Oceananigans.BoundaryConditions.NormalRadiation,
            }

            # Reading a state from a forcing file that does not exist is a stated error rather than a
            # failure inside `NCDataset(nothing)` after the grid has been built.
            @test_throws ErrorException FjordSim.Simulations.resolve_initial_conditions(
                FromForcing(), nothing, nothing, config.simulation_config,
            )
        end
    end

    @testset "run! smoke test" begin
        # No other test in this suite calls `run!`/`time_step!` for real — every other
        # `build_simulation` call either stops at a prerequisite check or inspects the assembled
        # model without stepping it. This is the check that would actually catch a `NaN` entering
        # through the new tracer open boundary condition (or any other boundary-condition wiring),
        # rather than just asserting object types.
        mktempdir() do tmp
            Nx, Ny, Nz = 16, 16, 8
            bathymetry_config = test_bathymetry_config(data_root = tmp)
            (; grid_config) = immersed_test_grid(
                bathymetry_path(bathymetry_config);
                size = (Nx, Ny, Nz),
                halo = (7, 7, 7),
            )

            boundaries_config = test_boundaries_config(data_root = tmp)
            forcing_config = test_forcing_config(data_root = tmp)
            write_prepared_forcing(forcing_path(forcing_config); size = (Nx, Ny, Nz))
            write_prepared_boundaries(
                boundary_data_path(boundaries_config);
                size = (Nx, Ny, Nz),
                dates = [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)],
            )

            config = FjordConfig(;
                grid_config,
                bathymetry_config,
                forcing_config,
                boundary_config = boundaries_config,
                simulation_config = test_simulation_config(results_root = tmp),
            )

            simulation = build_simulation(config)
            run!(simulation)

            ocean_model = simulation.model.ocean.model
            @test all(isfinite, interior(ocean_model.velocities.u))
            @test all(isfinite, interior(ocean_model.velocities.v))
            @test all(isfinite, interior(ocean_model.tracers.T))
            @test all(isfinite, interior(ocean_model.tracers.S))

            # The free surface and the barotropic transport are the two fields only the open
            # boundary drives, and the two the split-explicit substepping fills every substep. A
            # boundary that walled itself off, or a `NaN` from the Flather condition's own
            # arithmetic, would show up here and nowhere else in the suite.
            free_surface = ocean_model.free_surface
            @test all(isfinite, interior(free_surface.displacement))
            @test all(isfinite, interior(free_surface.barotropic_velocities.V))
        end
    end
end
