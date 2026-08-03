using FjordSim
using FjordSim.Bathymetry: write_bathymetry_file
using Dates: DateTime, Hour
using Test
using ArchGDAL
using NCDatasets
using Oceananigans
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: AbstractField
using NumericalEarth

# Alternative config types used by the "Config extensibility" testset. They exist only to
# check that `FjordConfig` accepts any subtype of the abstract config supertypes and that new
# behavior is added by overloading, without editing FjordSim. Defined at top level because a
# @testset body is a function scope, where `struct` is not allowed.
struct SingleColumnGrid <: AbstractGridConfig
    depth::Float64
end

struct MinimalBathymetry <: AbstractBathymetryConfig
    data_root::String
    output_file::String
    plot_file::String
end

struct ConstantForcing <: AbstractForcingConfig
    data_root::String
    output_file::String
    plot_file::String
    temperature::Float64
    rivers::Nothing
end

struct MinimalRivers <: AbstractRiverConfig
    data_root::String
    output_file::String
    relaxation_timescale::Float64
    search_radius::Int
end

struct MinimalAtmosphere <: AbstractAtmosphereConfig
    data_root::String
    output_file::String
    plot_file::String
    output_directory::String
    resolution::Float64
    padding::Float64
end

struct MinimalSimulation <: AbstractSimulationConfig
    results_root::String
    output_file::String
end

# A river config the "Add rivers round-trip" testset drives `add_rivers` with, standing in for a
# real dataset: fixed outlets and a value that counts up in time, so a misplaced write is
# visible.
struct StubRivers <: AbstractRiverConfig
    data_root::String
    output_file::String
    relaxation_timescale::Float64
    search_radius::Int
    locations::Vector{FjordSim.Forcing.RiverLocation}
    series::Dict{String,Matrix{Float32}}
end

FjordSim.Forcing.river_locations(config::StubRivers) = config.locations
FjordSim.Forcing.river_series(config::StubRivers, times) = config.series

# New behavior for a new grid config, added without touching Grids.jl.
Oceananigans.LatitudeLongitudeGrid(architecture, config::SingleColumnGrid) = LatitudeLongitudeGrid(
    architecture;
    size = (1, 1, 2),
    halo = (1, 1, 1),
    longitude = (10.0, 11.0),
    latitude = (59.0, 60.0),
    z = [-config.depth, -config.depth / 2, 0.0],
)

# One forcing hook overloaded on the new forcing config, without touching Forcing.jl.
FjordSim.Forcing.forcing_variable_names(config::ConstantForcing) = Dict("temperature" => "T")

@testset "Backward Compatibility — API Exports" begin
    # Verify all exported symbols are present in the public interface
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
        :ProjectedSourceGrid,
        :RiverLocation,
        :geodatabase_path,
        :top_bottom_boundary_conditions,
        :coupled_hydrostatic_simulation,
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
        :build_simulation,
        :run_simulation,
        :simulation_architecture,
        :results_path,
        :prescribed_atmosphere,
        :prescribed_radiation,
    ]

    for sym in exported_symbols
        @test isdefined(FjordSim, sym)  # All core exports must be defined
    end

    # Verify submodule exports
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
        (:Simulations, :coupled_hydrostatic_simulation),
        (:CLI, :parse_arguments),
    ]

    for (module_name, sym) in submodule_exports
        mod = getfield(FjordSim, module_name)
        @test isdefined(mod, sym)  # All submodule exports must be defined
    end

    @test isdefined(FjordSim, :Bathymetry)
    @test isdefined(FjordSim, :Setups)
    @test isdefined(FjordSim, :Simulations)
    @test isdefined(FjordSim, :CLI)
    @test parentmodule(coupled_hydrostatic_simulation) === FjordSim.Simulations
end

@testset "Backward Compatibility — Function Signatures" begin
    mktempdir() do tmp
        # Create minimal test files
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

        architecture = CPU()
        grid = ImmersedBoundaryGrid(bathymetry_file, architecture, (1, 1, 1))

        # Test top_bottom_boundary_conditions signature
        @test_nowarn top_bottom_boundary_conditions(; grid, bottom_drag_coefficient = 0.003)
        bcs = top_bottom_boundary_conditions(; grid, bottom_drag_coefficient = 0.003)
        
        # Verify return type: should be a NamedTuple with u, v, T, S
        @test haskey(bcs, :u)  # Return must include u
        @test haskey(bcs, :v)  # Return must include v
        @test haskey(bcs, :T)  # Return must include T
        @test haskey(bcs, :S)  # Return must include S
        
        # Verify nested structure for u and v: (top, bottom)
        @test haskey(bcs.u, :top)  # u must have top
        @test haskey(bcs.u, :bottom)  # u must have bottom
        @test haskey(bcs.v, :top)  # v must have top
        @test haskey(bcs.v, :bottom)  # v must have bottom
        
        # Verify return types are callable boundary conditions (not checking exact type, just that they exist and are not nothing)
        @test !isnothing(bcs.u.top)  # u.top must exist
        @test !isnothing(bcs.u.bottom)  # u.bottom must exist
        @test !isnothing(bcs.v.top)  # v.top must exist
        @test !isnothing(bcs.v.bottom)  # v.bottom must exist
        @test !isnothing(bcs.T.top)  # T.top must exist
        @test !isnothing(bcs.S.top)  # S.top must exist

        # The tracer top conditions must carry a freshwater exchange NumericalEarth can read back
        # out of them: `net_fluxes` pulls the volume flux and its heat content straight from the
        # boundary condition, so a bare FluxBoundaryCondition raises a MethodError as soon as the
        # coupled model assembles fluxes — which nothing else here would catch.
        @test NumericalEarth.Oceans.extract_freshwater_flux(bcs.S.top.condition) isa AbstractField
        @test NumericalEarth.Oceans.extract_freshwater_flux(bcs.T.top.condition) isa AbstractField
        @test NumericalEarth.Oceans.freshwater_exchange(bcs.T.top.condition).content_flux isa AbstractField
    end
end

@testset "Backward Compatibility — Struct Fields" begin
    mktempdir() do tmp
        # Create minimal NORA3 file. The dimensions are `lon`/`lat` because that is the prepared
        # file's contract — `NORA3FieldTimeSeries` reads the axes by those names, so a file
        # dimensioned otherwise could never be read back.
        nora3_filename = "NORA3_test.nc"
        nora3_path = joinpath(tmp, nora3_filename)
        ds = NCDataset(nora3_path, "c")
        defDim(ds, "lon", 2)
        defDim(ds, "lat", 2)
        defDim(ds, "time", 2)
        air_temperature_2m = defVar(ds, "air_temperature_2m", Float64, ("lon", "lat", "time"))
        time_var = defVar(ds, "time", Float64, ("time",))
        air_temperature_2m[:, :, :] .= 273.15  # Use broadcasting assignment
        time_var[:] = [0.0, 3600.0]
        close(ds)

        # Test MultiYearNORA3 struct fields
        nora3 = MultiYearNORA3(nora3_filename, tmp)
        @test hasfield(typeof(nora3), :metadata_filename)  # Field required
        @test hasfield(typeof(nora3), :default_download_directory)  # Field required
        @test hasfield(typeof(nora3), :size)  # Field required
        @test hasfield(typeof(nora3), :all_dates)  # Field required
        
        # Verify field types
        @test isa(nora3.metadata_filename, String)  # metadata_filename is String
        @test isa(nora3.default_download_directory, String)  # default_download_directory is String
        @test isa(nora3.size, Tuple)  # size is Tuple
        @test isa(nora3.all_dates, Vector)  # all_dates is Vector
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

        @test DW.is_three_dimensional(vector_metadata) == false     # upstream default is `true`
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

@testset "Setup configs" begin
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
    shared = DybdedataConfig(
        data_root = data_root,
        output_file = "bathymetry.nc",
        plot_file = "bathymetry.png",
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

    # The prepared forcing file, its plot and the relaxation zone are defaulted, since there is
    # one prepared forcing file per setup and the notebook-derived relaxation matches both
    # current setups.
    @test forcing_path(forcing_config) == joinpath(data_root, "forcing.nc")
    @test plot_path(forcing_config) == joinpath(data_root, "forcing.png")
    @test forcing_config.relaxation_edge === :south
    @test forcing_config.relaxation_cells == 10
    @test forcing_config.relaxation_timescale == 86400.0

    # Where the interpolation runs is a config field, not a command-line flag. `:auto` is the
    # default so one setup runs on a GPU machine and a laptop alike; `:cpu` is honoured
    # regardless of the hardware present, which is what makes this assertion machine-independent.
    @test forcing_config.architecture === :auto
    @test interpolation_architecture(
        NorKystConfig(
            data_root = data_root,
            output_directory = "norkyst",
            architecture = :cpu,
            parameters = ["temperature"],
            years = [2020],
        ),
    ) == CPU()
    # An unknown selector is rejected up front rather than falling back to some default.
    @test_throws ArgumentError interpolation_architecture(
        NorKystConfig(
            data_root = data_root,
            output_directory = "norkyst",
            architecture = :tpu,
            parameters = ["temperature"],
            years = [2020],
        ),
    )

    # ...and an absolute path relocates just that file, as for the bathymetry config.
    relocated = NorKystConfig(
        data_root = data_root,
        output_directory = "norkyst",
        output_file = "/shared/forcing.nc",
        parameters = ["temperature"],
        years = [2020],
    )
    @test forcing_path(relocated) == "/shared/forcing.nc"
    @test plot_path(relocated) == joinpath(data_root, "forcing.png")
    @test forcing_directory(
        NorKystConfig(
            data_root = data_root,
            output_directory = "/nk",
            parameters = ["temperature"],
            years = [2020],
        ),
    ) == "/nk"
    @test forcing_monthly_filename(forcing_config, 2020, 3) == "NorKyst-800m_ZDEPTHS_avg_202003.nc"
    @test occursin("thredds.met.no", forcing_config.catalog_url)

    # A setup states its own atmosphere directory and years; only the NORA3 endpoint, the
    # prepared file names and the target resolution are defaulted.
    @test_throws UndefKeywordError NORA3Config(data_root = data_root, years = [2020])
    @test_throws UndefKeywordError NORA3Config(data_root = data_root, output_directory = "nora3")
    @test_throws UndefKeywordError NORA3Config(output_directory = "nora3", years = [2020])

    atmosphere_config = NORA3Config(data_root = data_root, output_directory = "nora3", years = [2020])
    @test atmosphere_config isa AbstractAtmosphereConfig
    @test isconcretetype(typeof(atmosphere_config))
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
        NORA3Config(data_root = data_root, output_directory = "/nora3", years = [2020]),
    ) == "/nora3"

    grid_config = EvenGrid(
        size = (2, 3, 2),
        halo = (1, 1, 1),
        longitude = (10.0, 12.0),
        latitude = (59.0, 62.0),
        z_faces = [-20.0, -10.0, 0.0],
    )

    config = FjordConfig(; grid_config, bathymetry_config, forcing_config, atmosphere_config)
    grid = LatitudeLongitudeGrid(CPU(), config.grid_config)

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

@testset "Setup registry" begin
    @test setup_names() == ["drammensfjorden", "oslofjorden"]  # sorted: Dict order is unspecified
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
        @test isconcretetype(typeof(config))  # every instantiation stays concretely typed
        @test config.bathymetry_config.data_root == data_root
        @test config.forcing_config.data_root == data_root
        @test bathymetry_path(config.bathymetry_config) == joinpath(data_root, "bathymetry.nc")
        @test forcing_path(config.forcing_config) == joinpath(data_root, "forcing.nc")
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

    oslo = oslofjorden()
    @test oslo.forcing_config.rivers isa OF800RiversConfig
    @test oslo.atmosphere_config isa NORA3Config
    # The river data downloads under the setup's own root rather than being shared from elsewhere.
    @test oslo.forcing_config.rivers.data_root == oslo.forcing_config.data_root
    @test atmosphere_path(oslo.atmosphere_config) ==
          joinpath(oslo.atmosphere_config.data_root, "atmosphere.nc")
    # Results are rooted separately from the input data, so `results_root` is its own path.
    @test oslo.simulation_config isa SimulationConfig
    @test results_path(oslo.simulation_config) ==
          joinpath(homedir(), "FjordSim_results", "oslofjorden", "snapshots_ocean.nc")

    drammen = drammensfjorden()
    @test isnothing(drammen.forcing_config.rivers)     # so add_rivers is a no-op for it
    @test isnothing(drammen.atmosphere_config)         # ...and so are both atmosphere steps
    @test isnothing(drammen.simulation_config)         # ...and so is run_simulation
end

@testset "CLI" begin
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

    subcommands = FjordSim.CLI.subcommand_names()
    @test subcommands == [
        "prepare_bathymetry",
        "download_forcing",
        "prepare_forcing",
        "add_rivers",
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

    @test isdefined(FjordSim, :main)
    @test hasmethod(FjordSim.main, Tuple{Vector{String}})
    # Exporting `main` would make Julia run the CLI after the body of any script that says
    # `using FjordSim` — including this test file.
    @test !(:main in names(FjordSim))
    @test FjordSim.main(["--help"]) === 0            # must be an exit code, not a driver's result
    @test FjordSim.main(["frobnicate", "--config", "oslofjorden"]) == 2
    @test FjordSim.main(["prepare_forcing", "--config", "nordfjorden"]) == 2
    @test FjordSim.main(String[]) == 2

    # A step the setup opts out of is a no-op through the FjordConfig arity too, matching
    # `add_rivers(grid, config, ::Nothing)` and `prepare_atmosphere(grid, ::Nothing)`. Rooted in a
    # temporary directory so nothing can touch the real ~/FjordSim_data.
    mktempdir() do tmp
        config = drammensfjorden()
        config.bathymetry_config.data_root = tmp
        config.forcing_config.data_root = tmp

        @test isnothing(add_rivers(config))
        @test isnothing(prepare_atmosphere(config))
        @test isnothing(download_atmosphere(config))
        # Both return on the `nothing` simulation config, before any prerequisite is checked, so
        # neither touches the empty temporary root.
        @test isnothing(build_simulation(config))
        @test isnothing(run_simulation(config))
        @test FjordSim.main(["add_rivers", "--config", joinpath(tmp, "unused.jl")]) == 2
    end

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

    # A step that fails is exit code 1, and its error and backtrace are in the log rather than only
    # in the terminal's scrollback — which is the whole point of the tee, since a stacktrace through
    # `SimulationConfig` is long enough to scroll the error message itself away. An out-of-tree
    # config rooted in `tmp` rather than the registered name, so the assertion does not depend on
    # whether ~/FjordSim_data/drammensfjorden happens to exist.
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

            logged = read(FjordSim.CLI.LOG_FILE, String)
            @test occursin("Processed bathymetry", logged)
            @test occursin("Stacktrace", logged)
        end
    end

    # The log goes beside the output it describes. `results_root` is a simulation-config field, so a
    # setup naming no simulation config has no results directory and falls back to the cwd — the
    # `drammensfjorden` case the test above relies on.
    oslo = oslofjorden()
    @test FjordSim.CLI.log_path(oslo) ==
          joinpath(oslo.simulation_config.results_root, FjordSim.CLI.LOG_FILE)
    @test FjordSim.CLI.log_path(drammensfjorden()) == FjordSim.CLI.LOG_FILE
    @test FjordSim.CLI.log_path(nothing) == FjordSim.CLI.LOG_FILE

    # End to end: a failing step on a setup that does name a results root writes its log there and
    # not into the working directory, creating the directory if it does not exist yet.
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
            @test FjordSim.main(["prepare_forcing", "--config", config_file]) == 1
            @test !isfile(FjordSim.CLI.LOG_FILE)

            logged = read(joinpath(results, FjordSim.CLI.LOG_FILE), String)
            @test occursin("Processed bathymetry", logged)
        end
    end
end

@testset "OF800 rivers config" begin
    rivers = OF800RiversConfig(data_root = "/data/oslofjord")

    @test rivers isa AbstractRiverConfig
    @test isconcretetype(typeof(rivers))
    @test river_forcing_path(rivers) == "/data/oslofjord/forcing_rivers.nc"
    @test FjordSim.Forcing.river_locations_path(rivers) == "/data/oslofjord/OF800_rivers.csv"
    @test FjordSim.Forcing.river_series_path(rivers) == "/data/oslofjord/of800_rivers_v9_1990_2022_RA1.nc"
    @test river_search_radius(rivers) == 10
    @test rivers.relaxation_timescale == 3600.0

    # Both source files download from per-file links by default; a folder link cannot work.
    @test occursin("/scl/fi/", rivers.locations_url) && occursin("dl=1", rivers.locations_url)
    @test occursin("/scl/fi/", rivers.series_url) && occursin("dl=1", rivers.series_url)

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

@testset "Simulation config" begin
    open_boundary_conditions = FjordSim.Simulations.open_boundary_conditions
    simulation_forcing_path = FjordSim.Simulations.simulation_forcing_path

    # No field has a default, so every one of them is named by the setup and omitting any is an
    # UndefKeywordError rather than a plausible-looking run. `results_root` is the last one left,
    # so dropping it from an otherwise complete argument set still fails.
    fields = Dict(name => getfield(oslofjorden().simulation_config, name)
                  for name in fieldnames(SimulationConfig))
    @test_throws UndefKeywordError SimulationConfig(; delete!(copy(fields), :results_root)...)
    @test_throws UndefKeywordError SimulationConfig(; delete!(copy(fields), :coriolis)...)
    @test_throws UndefKeywordError SimulationConfig(; delete!(copy(fields), :tracer_advection)...)
    @test_throws UndefKeywordError SimulationConfig(; delete!(copy(fields), :stop_time)...)
    @test_throws UndefKeywordError SimulationConfig()

    config = oslofjorden().simulation_config
    @test config isa AbstractSimulationConfig
    @test isconcretetype(typeof(config))
    @test all(isconcretetype, fieldtypes(typeof(config)))
    # One parameter per prebuilt Oceananigans component. A tenth is the signal to collapse them
    # into a single `NamedTuple` field instead of adding another.
    @test length(typeof(config).parameters) == 9

    # The device is a Symbol, not a live CPU()/GPU(), so a setup file loads without a GPU.
    @test fieldtype(typeof(config), :architecture) === Symbol
    @test results_path(config) ==
          joinpath(homedir(), "FjordSim_results", "oslofjorden", "snapshots_ocean.nc")
    config.output_file = "/elsewhere/run.nc"      # an absolute path relocates just that file
    @test results_path(config) == "/elsewhere/run.nc"

    # The values Oslofjord runs with, stated in the setup rather than inherited. Durations are
    # seconds, and must be Float64 — @kwdef on a parametric struct does not convert.
    @test config.tracers == (:T, :S)
    @test config.initial_conditions == (T = 5.0, S = 33.0)
    @test isnothing(config.biogeochemistry)
    @test config.architecture == :auto
    @test config.stop_time == 365 * 24 * 3600.0
    @test config.output_interval == 3600.0
    @test config.max_time_step == 180.0

    # Architecture resolution shares `interpolation_architecture`'s Val methods, so :cpu pins the
    # CPU and an unknown selector is rejected rather than silently defaulted.
    config.architecture = :cpu
    @test simulation_architecture(config) isa CPU
    config.architecture = :tpu
    @test_throws ArgumentError simulation_architecture(config)

    # The open boundary is derived from the forcing config's `relaxation_edge`, on the velocity
    # normal to that edge, so the two cannot drift apart.
    @test keys(open_boundary_conditions(Val(:south))) == (:v,)
    @test keys(open_boundary_conditions(Val(:south)).v) == (:south,)
    @test keys(open_boundary_conditions(Val(:north))) == (:v,)
    @test keys(open_boundary_conditions(Val(:west))) == (:u,)
    @test keys(open_boundary_conditions(Val(:east)).u) == (:east,)
    @test_throws ArgumentError open_boundary_conditions(Val(:middle))

    # A setup that configures rivers simulates the rivers-augmented copy, never the pre-rivers
    # file; one with no rivers reads what prepare_forcing wrote.
    oslo = oslofjorden()
    @test simulation_forcing_path(oslo) == river_forcing_path(oslo.forcing_config.rivers)
    drammen = drammensfjorden()
    @test simulation_forcing_path(drammen) == forcing_path(drammen.forcing_config)

    # Every prerequisite is reported as the command that produces it, and both checks run before
    # anything is read or allocated, so the temporary root stays empty apart from the stub.
    mktempdir() do tmp
        config = oslofjorden()
        config.bathymetry_config.data_root = tmp
        config.forcing_config.data_root = tmp
        config.forcing_config.rivers.data_root = tmp

        @test_throws "prepare_bathymetry" build_simulation(config)

        touch(bathymetry_path(config.bathymetry_config))
        @test_throws "add_rivers" build_simulation(config)

        # ...and a setup with no rivers names prepare_forcing instead. Assembled rather than
        # mutated, because `rivers` is a type parameter of NorKystConfig and cannot be unset.
        without_rivers = drammensfjorden()
        no_rivers = FjordConfig(
            grid_config = without_rivers.grid_config,
            bathymetry_config = config.bathymetry_config,
            forcing_config = without_rivers.forcing_config,
            simulation_config = config.simulation_config,
        )
        no_rivers.forcing_config.data_root = tmp
        @test isnothing(no_rivers.forcing_config.rivers)
        @test_throws "prepare_forcing" build_simulation(no_rivers)
    end
end

@testset "Config extensibility" begin
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
        simulation_config = MinimalSimulation(data_root, "column_snapshots.nc"),
    )
    rivers = MinimalRivers(data_root, "column_rivers.nc", 3600.0, 10)

    @test config.grid_config isa AbstractGridConfig
    @test config.bathymetry_config isa AbstractBathymetryConfig
    @test config.forcing_config isa AbstractForcingConfig
    @test config.atmosphere_config isa AbstractAtmosphereConfig
    @test config.simulation_config isa AbstractSimulationConfig
    @test rivers isa AbstractRiverConfig

    # Field types are still concrete, so the struct stays type-stable per instantiation.
    @test isconcretetype(typeof(config))
    @test fieldtype(typeof(config), :grid_config) === SingleColumnGrid
    @test fieldtype(typeof(config), :atmosphere_config) === MinimalAtmosphere
    @test fieldtype(typeof(config), :simulation_config) === MinimalSimulation

    # ...and the built-in types remain valid, i.e. FjordConfig is genuinely generic.
    @test FjordConfig(
        grid_config = EvenGrid(
            size = (2, 3, 2),
            halo = (1, 1, 1),
            longitude = (10.0, 12.0),
            latitude = (59.0, 62.0),
            z_faces = [-20.0, -10.0, 0.0],
        ),
        bathymetry_config = DybdedataConfig(
            data_root = data_root,
            output_file = "bathymetry.nc",
            plot_file = "bathymetry.png",
        ),
        forcing_config = NorKystConfig(
            data_root = data_root,
            output_directory = "norkyst",
            parameters = ["temperature"],
            years = [2020],
        ),
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
    # The simulation config resolves against `results_root`, being the one config that writes.
    @test results_path(config.simulation_config) == joinpath(data_root, "column_snapshots.nc")

    # `river_search_radius` is the river pipeline's one optional hook.
    @test river_search_radius(rivers) == 10

    # A method overloaded on the new grid config is picked up by existing call sites.
    grid = LatitudeLongitudeGrid(CPU(), config.grid_config)
    @test size(grid) == (1, 1, 2)
    @test collect(Oceananigans.Grids.znodes(grid, Face())) == [-120.0, -60.0, 0.0]

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

@testset "FjordSim.jl" begin
    mktempdir() do tmp
        bathymetry_file = joinpath(tmp, "bathymetry.nc")
        nora3_filename = "NORA3_test.nc"
        nora3_path = joinpath(tmp, nora3_filename)

        # Minimal bathymetry file for Grids.ImmersedBoundaryGrid.
        ds = NCDataset(bathymetry_file, "c")
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "zf", 3)
        z_faces = defVar(ds, "z_faces", Float64, ("zf",))
        h = defVar(ds, "h", Float64, ("x", "y"))
        lat = defVar(ds, "lat", Float64, ("x",))
        lon = defVar(ds, "lon", Float64, ("y",))
        z_faces[:] = [-20.0, -10.0, 0.0]
        h[:, :] = fill(10.0, (2, 2))  # Create and assign a 2x2 matrix of 10.0
        lat[:] = [59.0, 60.0]
        lon[:] = [10.0, 11.0]
        close(ds)

        # Minimal NORA3 file for MultiYearNORA3 dataset construction, on the prepared file's own
        # `lon`/`lat`/`time` dimensions.
        ds = NCDataset(nora3_path, "c")
        defDim(ds, "lon", 2)
        defDim(ds, "lat", 2)
        defDim(ds, "time", 2)
        air_temperature_2m = defVar(ds, "air_temperature_2m", Float64, ("lon", "lat", "time"))
        time = defVar(ds, "time", Float64, ("time",))
        air_temperature_2m[:, :, :] .= 273.15
        time[:] = [0.0, 3600.0]
        close(ds)

        architecture = CPU()
        grid = @test_nowarn ImmersedBoundaryGrid(bathymetry_file, architecture, (1, 1, 1))
        @test_nowarn top_bottom_boundary_conditions(; grid, bottom_drag_coefficient = 0.003)
        @test_nowarn MultiYearNORA3(nora3_filename, tmp)
    end
end

@testset "Bathymetry writer" begin
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

@testset "Bathymetry gap filling" begin
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
    for i = 1:3, j = 1:3
        smoothed[i, j] < 0 || continue
        for (di, dj) in ((1, 0), (0, 1))
            ii, jj = i + di, j + dj
            (ii <= 3 && jj <= 3) || continue
            smoothed[ii, jj] < 0 || continue
            here, there = -smoothed[i, j], -smoothed[ii, jj]
            @test abs(here - there) / (here + there) <= 0.4 + 1e-9
        end
    end

    # Both stages are opt-in: the default keeps the topological cleanup and nothing else, which is
    # what leaves every other source's bathymetry unchanged.
    @test FjordSim.Bathymetry.smoothing_options(drammensfjorden().bathymetry_config) ==
          (spike_ratio = 0.0, max_slope_factor = 0.0, minimum_depth = 0.0)
    oslo_smoothing = FjordSim.Bathymetry.smoothing_options(oslofjorden().bathymetry_config)
    @test oslo_smoothing.spike_ratio == 0.5
    @test oslo_smoothing.max_slope_factor == 0.5
    @test oslo_smoothing.minimum_depth == 2.0

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
    for _ = 1:FjordSim.Bathymetry.BATHYMETRY_GAP_FILL_PASSES
        expected = fill_isolated_land_cells(remove_isolated_sea_cells(expected))
    end

    @test Array(interior(bottom_height, :, :, 1)) == expected
end

@testset "Bathymetry point sampling (bulk coordinate transform)" begin
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
    bathymetry_config = DybdedataConfig(
        data_root = tempdir(),
        output_file = "bathymetry.nc",
        plot_file = "bathymetry.png",
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

                for i = 1:n
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

@testset "Forcing preparation helpers" begin
    forcing_dimension_names = FjordSim.Forcing.forcing_dimension_names
    source_fill = FjordSim.Forcing.source_fill
    fill_source! = FjordSim.Forcing.fill_source!
    solve_vertical_faces = FjordSim.Forcing.solve_vertical_faces
    relaxation_lambda = FjordSim.Forcing.relaxation_lambda
    nearest_valid_map = FjordSim.Forcing.nearest_valid_map
    daily_time_steps = FjordSim.Forcing.daily_time_steps
    ForcingTimeStep = FjordSim.Forcing.ForcingTimeStep
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

    # The relaxation band hugs the configured edge and only covers water cells.
    wet = trues(3, 4, 1)
    band(edge) = begin
        config = NorKystConfig(
            data_root = tempdir(),
            output_directory = "norkyst",
            relaxation_edge = edge,
            relaxation_cells = 2,
            parameters = ["temperature"],
            years = [2020],
        )
        findall(!iszero, relaxation_lambda(wet, config)[:, :, 1])
    end
    @test Set(index[2] for index in band(:south)) == Set([1, 2])
    @test Set(index[2] for index in band(:north)) == Set([3, 4])
    @test Set(index[1] for index in band(:west)) == Set([1, 2])
    @test Set(index[1] for index in band(:east)) == Set([2, 3])

    dry = falses(3, 4, 1)
    dry_config = NorKystConfig(
        data_root = tempdir(),
        output_directory = "norkyst",
        relaxation_cells = 2,
        parameters = ["temperature"],
        years = [2020],
    )
    @test all(iszero, relaxation_lambda(dry, dry_config))  # land never relaxes
    @test maximum(relaxation_lambda(wet, dry_config)) == Float32(1 / 86400)
    @test_throws ArgumentError relaxation_lambda(
        wet,
        NorKystConfig(
            data_root = tempdir(),
            output_directory = "norkyst",
            relaxation_edge = :up,
            parameters = ["temperature"],
            years = [2020],
        ),
    )

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
    @test [step.date for step in steps] == [DateTime(2020, 1, d, 12) for d = 1:5]
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

@testset "Forcing water mask" begin
    water_mask = FjordSim.Forcing.water_mask

    mktempdir() do tmp
        architecture = CPU()
        # A 4x2 channel deepening to the east. Column 1 is land, column 2's bottom sits inside
        # the deep cell (a partial cell), columns 3-4 are fully wet. Both rows are identical, so
        # the staggered assertions below can compare against the tracer row directly.
        grid_config = EvenGrid(
            size = (4, 2, 2),
            halo = (1, 1, 1),
            longitude = (10.0, 10.4),
            latitude = (59.0, 59.2),
            z_faces = [-20.0, -10.0, 0.0],
        )
        underlying_grid = LatitudeLongitudeGrid(architecture, grid_config)
        bottom_height = Field{Center, Center, Nothing}(underlying_grid)
        set!(bottom_height, [0.0 0.0; -13.0 -13.0; -20.0 -20.0; -20.0 -20.0])

        bathymetry_file = joinpath(tmp, "bathymetry.nc")
        write_bathymetry_file(bathymetry_file, underlying_grid, bottom_height)
        grid = ImmersedBoundaryGrid(bathymetry_file, architecture, grid_config.halo)

        # The mask must agree with what the model treats as wet, which for PartialCellBottom is
        # decided by minimum_fractional_cell_height (0.2), not by a cell-centre test. Column 2's
        # deep cell keeps 7 m of its 10 m, so it is wet even though its centre (-15) is below the
        # bottom height (-13) — the old hand-rolled mask marked it land.
        tracer = water_mask(grid, Center, Center, :south)
        @test size(tracer) == (4, 2, 2)
        @test tracer[1, 1, 1] == false && tracer[1, 1, 2] == false  # land column
        @test tracer[2, 1, 1] == true                               # partial cell, still wet
        @test tracer[3, 1, 1] == true && tracer[4, 1, 1] == true

        # A velocity face is land when either tracer cell it separates is land, matching
        # Oceananigans: a face against land is a wall. The old mask called it wet.
        u = water_mask(grid, Face, Center, :south)
        @test size(u) == (5, 2, 2)
        @test u[2, 1, 1] == false  # between land column 1 and wet column 2 -> wall
        @test u[3, 1, 1] == true   # between two wet columns
        @test u[1, 1, 1] == false && u[5, 1, 1] == false  # closed east/west walls

        # ...but the open boundary named by relaxation_edge must survive, because that is where
        # the forcing relaxes and the setup puts a NormalFlowBoundaryCondition. peripheral_node alone
        # marks it land, since the tracer cell outside the domain is an inactive halo cell.
        v_south = water_mask(grid, Center, Face, :south)
        @test size(v_south) == (4, 3, 2)
        @test v_south[:, 1, :] == tracer[:, 1, :]      # southern row restored from the tracer row
        @test all(v_south[:, 3, :] .== false)          # northern row still a closed wall
        v_north = water_mask(grid, Center, Face, :north)
        @test all(v_north[:, 1, :] .== false)          # southern row now closed
        @test v_north[:, 3, :] == tracer[:, 2, :]      # northern row restored
        # A south edge must not touch the u mask, and a west edge must open its western column.
        @test water_mask(grid, Face, Center, :south)[1, 1, :] == [false, false]
        @test water_mask(grid, Face, Center, :west)[1, 1, :] == tracer[1, 1, :]
    end
end

@testset "Forcing file round-trip" begin
    mktempdir() do tmp
        architecture = CPU()
        grid_config = EvenGrid(
            size = (2, 3, 2),
            halo = (1, 1, 1),
            longitude = (10.0, 12.0),
            latitude = (59.0, 62.0),
            z_faces = [-20.0, -10.0, 0.0],
        )
        underlying_grid = LatitudeLongitudeGrid(architecture, grid_config)
        bottom_height = Field{Center, Center, Nothing}(underlying_grid)
        set!(bottom_height, fill(-20.0, (2, 3)))

        bathymetry_file = joinpath(tmp, "bathymetry.nc")
        write_bathymetry_file(bathymetry_file, underlying_grid, bottom_height)
        grid = ImmersedBoundaryGrid(bathymetry_file, architecture, grid_config.halo)

        forcing_config = NorKystConfig(
            data_root = tmp,
            output_directory = "norkyst",
            parameters = ["temperature", "salinity", "u_eastward", "v_northward"],
            years = [2020],
        )

        # A file in exactly the layout prepare_forcing writes: staggered variables on
        # face dimensions, land as the NaN fill value, times decodable to DateTime.
        ds = NCDataset(forcing_path(forcing_config), "c")
        defDim(ds, "Nx", 2)
        defDim(ds, "Ny", 3)
        defDim(ds, "Nz", 2)
        defDim(ds, "Nx_faces", 3)
        defDim(ds, "Ny_faces", 4)
        defDim(ds, "time", 2)
        defVar(ds, "time", [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)], ("time",))
        for name in ("T", "S", "u", "v")
            dimensions = FjordSim.Forcing.forcing_dimension_names(name)
            shape = (ds.dim[dimensions[1]], ds.dim[dimensions[2]], ds.dim[dimensions[3]], ds.dim[dimensions[4]])
            for variable_name in (name, name * "_lambda")
                variable = defVar(ds, variable_name, Float32, dimensions; attrib = ["_FillValue" => NaN32])
                variable[:, :, :, :] = fill(1.0f0, shape)
            end
        end
        close(ds)

        # The config method resolves the path itself; the filepath method still works.
        forcing = forcing_from_file(forcing_config; grid, tracers = (:T, :S))
        @test forcing isa NamedTuple
        @test Set(keys(forcing)) == Set((:T, :S, :u, :v))
        @test forcing_from_file(; grid, filepath = forcing_path(forcing_config), tracers = (:T, :S)) isa NamedTuple

        # Only requested tracers are picked up, alongside the velocities.
        @test Set(keys(forcing_from_file(forcing_config; grid, tracers = (:T,)))) == Set((:T, :u, :v))

        # A grid the file was not written for must be rejected rather than silently misread.
        other_grid_config = EvenGrid(
            size = (4, 3, 2),
            halo = (1, 1, 1),
            longitude = (10.0, 12.0),
            latitude = (59.0, 62.0),
            z_faces = [-20.0, -10.0, 0.0],
        )
        other_underlying = LatitudeLongitudeGrid(architecture, other_grid_config)
        other_bottom = Field{Center, Center, Nothing}(other_underlying)
        set!(other_bottom, fill(-20.0, (4, 3)))
        other_file = joinpath(tmp, "bathymetry_other.nc")
        write_bathymetry_file(other_file, other_underlying, other_bottom)
        other = ImmersedBoundaryGrid(other_file, architecture, other_grid_config.halo)
        @test_throws DimensionMismatch forcing_from_file(forcing_config; grid = other, tracers = (:T, :S))
    end
end

@testset "River cell snapping" begin
    is_coastal_cell = FjordSim.Forcing.is_coastal_cell
    nearest_coastal_cell = FjordSim.Forcing.nearest_coastal_cell
    coastal_water_mask = FjordSim.Forcing.coastal_water_mask
    river_cells = FjordSim.Forcing.river_cells

    # A hand-built mask exercises the two rules on their own: a cell is coastal when it is water
    # and touches land, so open water in the middle is not a valid river mouth.
    mask = trues(5, 5)
    mask[1, :] .= false          # a land column along the western edge
    @test is_coastal_cell(mask, 1, 3) == false   # land itself
    @test is_coastal_cell(mask, 2, 3) == true    # water touching the land column
    @test is_coastal_cell(mask, 4, 3) == false   # open water, no land neighbour
    @test is_coastal_cell(mask, 0, 3) == false   # outside the grid

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
        architecture = CPU()
        # A 5x4 basin with a land column at i = 2, so columns 1 and 3 are coastal while columns
        # 4 and 5 are open water. The land is off the domain edge because an outlet has to sit
        # strictly inside the grid to be accepted at all.
        grid_config = EvenGrid(
            size = (5, 4, 2),
            halo = (1, 1, 1),
            longitude = (10.0, 10.5),
            latitude = (59.0, 59.4),
            z_faces = [-20.0, -10.0, 0.0],
        )
        underlying_grid = LatitudeLongitudeGrid(architecture, grid_config)
        bottom_height = Field{Center, Center, Nothing}(underlying_grid)
        depths = fill(-20.0, (5, 4))
        depths[2, :] .= 0.0
        set!(bottom_height, depths)

        bathymetry_file = joinpath(tmp, "bathymetry.nc")
        write_bathymetry_file(bathymetry_file, underlying_grid, bottom_height)
        grid = ImmersedBoundaryGrid(bathymetry_file, architecture, grid_config.halo)

        # The mask comes from the same water_mask prepare_forcing uses, taken at the surface.
        mask = coastal_water_mask(grid, :south)
        @test size(mask) == (5, 4)
        @test mask[2, 1] == false
        @test mask[1, 1] == true && mask[3, 1] == true

        longitudes = Array(Oceananigans.Grids.λnodes(grid, Center()))
        latitudes = Array(Oceananigans.Grids.φnodes(grid, Center()))

        # An outlet on the land column relocates to the coast; one in open water does too; one
        # outside the domain is dropped rather than clamped to the nearest edge cell.
        locations = [
            FjordSim.Forcing.RiverLocation(1, "on land", longitudes[2], latitudes[2]),
            FjordSim.Forcing.RiverLocation(2, "open water", longitudes[4], latitudes[2]),
            FjordSim.Forcing.RiverLocation(3, "outside", 20.0, latitudes[2]),
        ]
        cells = river_cells(grid, locations, :south, 10)
        @test length(cells) == 2
        @test [cell.location.id for cell in cells] == [1, 2]
        @test (cells[1].i, cells[1].j, cells[1].distance) == (1, 2, 1.0)
        @test (cells[2].i, cells[2].j, cells[2].distance) == (3, 2, 1.0)
        @test all(mask[cell.i, cell.j] for cell in cells)

        # An outlet sitting exactly on the outermost node counts as outside, matching the
        # reference's strict bounds test.
        edge = FjordSim.Forcing.RiverLocation(4, "on the edge", longitudes[1], latitudes[2])
        @test isempty(river_cells(grid, [edge], :south, 10))

        # With no reach, the on-land outlet is dropped too rather than written into land.
        @test isempty(river_cells(grid, [locations[1]], :south, 0))
    end
end

@testset "Add rivers round-trip" begin
    mktempdir() do tmp
        architecture = CPU()
        grid_config = EvenGrid(
            size = (5, 4, 2),
            halo = (1, 1, 1),
            longitude = (10.0, 10.5),
            latitude = (59.0, 59.4),
            z_faces = [-20.0, -10.0, 0.0],
        )
        underlying_grid = LatitudeLongitudeGrid(architecture, grid_config)
        bottom_height = Field{Center, Center, Nothing}(underlying_grid)
        depths = fill(-20.0, (5, 4))
        depths[2, :] .= 0.0
        set!(bottom_height, depths)

        bathymetry_file = joinpath(tmp, "bathymetry.nc")
        write_bathymetry_file(bathymetry_file, underlying_grid, bottom_height)
        grid = ImmersedBoundaryGrid(bathymetry_file, architecture, grid_config.halo)

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
        forcing_config = NorKystConfig(
            data_root = tmp,
            output_directory = "norkyst",
            parameters = ["temperature", "salinity"],
            years = [2020],
            rivers = rivers,
        )
        @test forcing_config.rivers isa AbstractRiverConfig
        @test isconcretetype(typeof(forcing_config))

        ds = NCDataset(forcing_path(forcing_config), "c")
        defDim(ds, "Nx", 5)
        defDim(ds, "Ny", 4)
        defDim(ds, "Nz", 2)
        defDim(ds, "Nx_faces", 6)
        defDim(ds, "Ny_faces", 5)
        defDim(ds, "time", 2)
        defVar(ds, "time", [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)], ("time",))
        for name in ("T", "S", "u", "v")
            dimensions = FjordSim.Forcing.forcing_dimension_names(name)
            shape = (ds.dim[dimensions[1]], ds.dim[dimensions[2]], ds.dim[dimensions[3]], ds.dim[dimensions[4]])
            defVar(ds, name, Float32, dimensions; attrib = ["_FillValue" => NaN32])[:, :, :, :] =
                fill(1.0f0, shape)
            defVar(ds, name * "_lambda", Float32, dimensions; attrib = ["_FillValue" => NaN32])[:, :, :, :] =
                fill(2.0f-5, shape)
        end
        close(ds)

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

@testset "Atmosphere preparation helpers" begin
    atmospheres = FjordSim.Atmospheres
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
    @test first(centers) == 10.1
    @test all(isapprox(0.02), diff(centers))
    @test last(centers) >= 11.12
    @test length(centers) == 52

    config = NORA3Config(data_root = tempdir(), output_directory = "nora3", years = [2020])
    @test_throws ArgumentError atmosphere_target_axes(
        LatitudeLongitudeGrid(CPU(); size = (1, 1, 1), longitude = (10.0, 11.0), latitude = (59.0, 60.0), z = (-1.0, 0.0)),
        NORA3Config(data_root = tempdir(), output_directory = "n", resolution = 0.0, years = [2020]),
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
    year = vcat((nora3_month_dates(2020, month) for month = 1:12)...)
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
    contiguous = [AtmosphereRecord(DateTime(2020, 1, 1, hour), "a.nc", hour + 1) for hour = 0:3]
    @test isnothing(validate_atmosphere_records(contiguous))
    gapped = [
        AtmosphereRecord(DateTime(2020, 1, 1, 0), "a.nc", 1),
        AtmosphereRecord(DateTime(2020, 1, 1, 5), "a.nc", 2),
    ]
    @test_logs (:warn,) validate_atmosphere_records(gapped)
end

@testset "Atmosphere file round-trip" begin
    atmospheres = FjordSim.Atmospheres
    # The real NORA3 projection, so the coordinate transform is exercised for real.
    proj4 = "+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +no_defs +R=6.371e+06"

    mktempdir() do tmp
        config = NORA3Config(data_root = tmp, output_directory = "nora3", years = [2020])
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
                for step = 1:length(dates)
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

@testset "Atmosphere reader window" begin
    # A prepared file with fewer steps than the backend's default window. Written by hand rather
    # than through `prepare_atmosphere`, because the point is the reader, not the regrid.
    atmospheres = FjordSim.Atmospheres
    nora3 = atmospheres.NORA3

    mktempdir() do tmp
        config = NORA3Config(data_root = tmp, output_directory = "nora3", years = [2020])
        dates = [DateTime(2020, 1, 1) + Hour(step - 1) for step = 1:3]
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

@testset "Simulation time coverage" begin
    validate_time_coverage = FjordSim.Simulations.validate_time_coverage
    forcing_date_range = FjordSim.Simulations.forcing_date_range
    range = (DateTime(2020, 1, 1, 12), DateTime(2020, 12, 31, 12))

    # Both readers use `Cyclical()` time indexing, which wraps rather than failing outside its
    # data, so a run that outlasts its forcing would quietly replay the beginning.
    @test isnothing(validate_time_coverage("forcing", range, DateTime(2020, 1, 1, 12), 86400.0, "x"))
    @test isnothing(validate_time_coverage("forcing", range, DateTime(2020, 6, 1), 86400.0, "x"))
    # A source that cannot report its dates is skipped rather than blocking the run.
    @test isnothing(validate_time_coverage("atmosphere", nothing, DateTime(2019, 1, 1), 1e9, "x"))

    @test_throws ErrorException validate_time_coverage(
        "forcing", range, DateTime(2020, 1, 1), 86400.0, "x",   # starts before the data
    )
    @test_throws ErrorException validate_time_coverage(
        "forcing", range, DateTime(2020, 12, 31), 86400.0, "x", # runs past the data
    )

    mktempdir() do tmp
        filepath = joinpath(tmp, "forcing.nc")
        dates = [DateTime(2020, 1, 1, 12) + Hour(24 * (step - 1)) for step = 1:5]
        NCDataset(filepath, "c") do ds
            defDim(ds, "time", length(dates))
            defVar(ds, "time", dates, ("time",))
        end
        @test forcing_date_range(filepath) == (first(dates), last(dates))
    end
end

@testset "Build simulation" begin
    # The "Simulation config" testset only reaches build_simulation's two prerequisite errors, so
    # everything past them — the HydrostaticFreeSurfaceModel, the coupled OceanSeaIceModel and the
    # tracer top boundary conditions NumericalEarth reads its freshwater exchange out of — is
    # otherwise never executed. Synthetic inputs keep it on the CPU with no ~/FjordSim_data.
    mktempdir() do tmp
        # Small, but the halo is oslofjorden's: the setup's WENO and biharmonic closure set the
        # floor, and an ImmersedBoundaryGrid needs one point more than that in every direction.
        architecture = CPU()
        Nx, Ny, Nz = 16, 16, 8
        grid_config = EvenGrid(
            size = (Nx, Ny, Nz),
            halo = (7, 7, 7),
            longitude = (10.0, 11.0),
            latitude = (59.0, 60.0),
            z_faces = collect(range(-80.0, 0.0, length = Nz + 1)),
        )

        bathymetry_config = DybdedataConfig(
            data_root = tmp,
            output_file = "bathymetry.nc",
            plot_file = "bathymetry.png",
        )
        underlying_grid = LatitudeLongitudeGrid(architecture, grid_config)
        bottom_height = Field{Center, Center, Nothing}(underlying_grid)
        set!(bottom_height, fill(-80.0, (Nx, Ny)))
        write_bathymetry_file(bathymetry_path(bathymetry_config), underlying_grid, bottom_height)

        forcing_config = NorKystConfig(
            data_root = tmp,
            output_directory = "norkyst",
            parameters = ["temperature", "salinity", "u_eastward", "v_northward"],
            years = [2020],
        )
        ds = NCDataset(forcing_path(forcing_config), "c")
        defDim(ds, "Nx", Nx)
        defDim(ds, "Ny", Ny)
        defDim(ds, "Nz", Nz)
        defDim(ds, "Nx_faces", Nx + 1)
        defDim(ds, "Ny_faces", Ny + 1)
        defDim(ds, "time", 2)
        defVar(ds, "time", [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)], ("time",))
        for name in ("T", "S", "u", "v")
            dimensions = FjordSim.Forcing.forcing_dimension_names(name)
            shape = ntuple(i -> ds.dim[dimensions[i]], 4)
            for variable_name in (name, name * "_lambda")
                variable = defVar(ds, variable_name, Float32, dimensions; attrib = ["_FillValue" => NaN32])
                variable[:, :, :, :] = fill(1.0f0, shape)
            end
        end
        close(ds)

        # Oslofjord's physics, on the CPU and stopped almost immediately. `atmosphere_config` stays
        # unnamed, so both prescribed_atmosphere hooks yield nothing and the coupled model runs
        # ocean-only — which is what keeps this off the network and off a GPU.
        simulation_config = oslofjorden().simulation_config
        simulation_config.architecture = :cpu
        simulation_config.results_root = tmp
        simulation_config.stop_time = 60.0

        config = FjordConfig(
            grid_config = grid_config,
            bathymetry_config = bathymetry_config,
            forcing_config = forcing_config,
            simulation_config = simulation_config,
        )
        @test isnothing(config.atmosphere_config)
        @test isnothing(config.forcing_config.rivers)

        simulation = build_simulation(config)
        @test simulation isa Simulation
        @test simulation.stop_time == 60.0

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

        # The open boundary landed on the velocity normal to the forcing config's relaxation edge.
        normal_velocity = config.forcing_config.relaxation_edge in (:south, :north) ? :v : :u
        edge_bc = getproperty(
            getproperty(ocean_model.velocities, normal_velocity).boundary_conditions,
            config.forcing_config.relaxation_edge,
        )
        @test edge_bc isa BoundaryCondition{<:Oceananigans.BoundaryConditions.NormalFlow}
    end
end
