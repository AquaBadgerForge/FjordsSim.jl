using FjordSim
using FjordSim.Bathymetry: write_bathymetry_file
using Dates: DateTime
using Test
using ArchGDAL
using NCDatasets
using Oceananigans
using Oceananigans.BoundaryConditions: FluxBoundaryCondition

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
    ]

    for (module_name, sym) in submodule_exports
        mod = getfield(FjordSim, module_name)
        @test isdefined(mod, sym)  # All submodule exports must be defined
    end

    @test isdefined(FjordSim, :Bathymetry)
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
    end
end

@testset "Backward Compatibility — Struct Fields" begin
    mktempdir() do tmp
        # Create minimal NORA3 file
        nora3_filename = "NORA3_test.nc"
        nora3_path = joinpath(tmp, nora3_filename)
        ds = NCDataset(nora3_path, "c")
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "time", 2)
        air_temperature_2m = defVar(ds, "air_temperature_2m", Float64, ("x", "y", "time"))
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

    grid_config = EvenGrid(
        size = (2, 3, 2),
        halo = (1, 1, 1),
        longitude = (10.0, 12.0),
        latitude = (59.0, 62.0),
        z_faces = [-20.0, -10.0, 0.0],
    )

    config = FjordConfig(; grid_config, bathymetry_config, forcing_config)
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
    @test config.grid_config.size == (4, 6, 2)
    @test config.bathymetry_config.padding_cells == 0
    @test config.forcing_config.years == [2020, 2021]
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

@testset "Config extensibility" begin
    data_root = joinpath(tempdir(), "fjordsim_extensibility_test")

    # FjordConfig accepts the alternative subtypes without any change to its definition.
    config = FjordConfig(
        grid_config = SingleColumnGrid(120.0),
        bathymetry_config = MinimalBathymetry(data_root, "column.nc", "column.png"),
        forcing_config = ConstantForcing(data_root, "column_forcing.nc", "column_forcing.png", 8.0, nothing),
    )
    rivers = MinimalRivers(data_root, "column_rivers.nc", 3600.0, 10)

    @test config.grid_config isa AbstractGridConfig
    @test config.bathymetry_config isa AbstractBathymetryConfig
    @test config.forcing_config isa AbstractForcingConfig
    @test rivers isa AbstractRiverConfig

    # Field types are still concrete, so the struct stays type-stable per instantiation.
    @test isconcretetype(typeof(config))
    @test fieldtype(typeof(config), :grid_config) === SingleColumnGrid

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

    # A forcing config carrying no rivers skips the step rather than needing a river dataset.
    @test isnothing(add_rivers(grid, config.forcing_config))
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

        # Minimal NORA3 file for MultiYearNORA3 dataset construction.
        ds = NCDataset(nora3_path, "c")
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "time", 2)
        air_temperature_2m = defVar(ds, "air_temperature_2m", Float64, ("x", "y", "time"))
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
        # the forcing relaxes and the setup puts an OpenBoundaryCondition. peripheral_node alone
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
