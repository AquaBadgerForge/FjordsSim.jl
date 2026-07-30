using FjordSim
using FjordSim.Bathymetry: write_bathymetry_file
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
    temperature::Float64
end

# New behavior for a new grid config, added without touching Grids.jl.
Oceananigans.LatitudeLongitudeGrid(arch, config::SingleColumnGrid) = LatitudeLongitudeGrid(
    arch;
    size = (1, 1, 2),
    halo = (1, 1, 1),
    longitude = (10.0, 11.0),
    latitude = (59.0, 60.0),
    z = [-config.depth, -config.depth / 2, 0.0],
)

@testset "Backward Compatibility — API Exports" begin
    # Verify all exported symbols are present in the public interface
    exported_symbols = [
        :ImmersedBoundaryGrid,
        :FjordConfig,
        :AbstractGridConfig,
        :AbstractBathymetryConfig,
        :AbstractForcingConfig,
        :EvenGrid,
        :DybdedataConfig,
        :NorKystConfig,
        :forcing_from_file,
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
        (:FDatasets, :DSForcing),
        (:FDatasets, :DSResults),
        (:FDatasets, :last_date),
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

        arch = CPU()
        grid = ImmersedBoundaryGrid(bathymetry_file, arch, (1, 1, 1))

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

    forcing_config = NorKystConfig(data_root = data_root, years = [2020])
    @test norkyst_directory(forcing_config) == joinpath(data_root, "norkyst")
    @test norkyst_directory(NorKystConfig(data_root = data_root, output_directory = "/nk", years = [2020])) == "/nk"
    @test norkyst_monthly_filename(2020, 3) == "NorKyst-800m_ZDEPTHS_avg_202003.nc"
    @test occursin("thredds.met.no", forcing_config.catalog_url)
    # Each config gets its own parameter vector, not the shared module-level default.
    @test forcing_config.parameters !== FjordSim.Forcing.NORKYST_PARAMETERS

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

@testset "Config extensibility" begin
    data_root = joinpath(tempdir(), "fjordsim_extensibility_test")

    # FjordConfig accepts the alternative subtypes without any change to its definition.
    config = FjordConfig(
        grid_config = SingleColumnGrid(120.0),
        bathymetry_config = MinimalBathymetry(data_root, "column.nc", "column.png"),
        forcing_config = ConstantForcing(8.0),
    )

    @test config.grid_config isa AbstractGridConfig
    @test config.bathymetry_config isa AbstractBathymetryConfig
    @test config.forcing_config isa AbstractForcingConfig

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
        forcing_config = NorKystConfig(data_root = data_root, years = [2020]),
    ) isa FjordConfig

    # Path resolution is inherited from AbstractBathymetryConfig — no new methods needed.
    @test bathymetry_path(config.bathymetry_config) == joinpath(data_root, "column.nc")
    @test plot_path(config.bathymetry_config) == joinpath(data_root, "column.png")

    # A method overloaded on the new grid config is picked up by existing call sites.
    grid = LatitudeLongitudeGrid(CPU(), config.grid_config)
    @test size(grid) == (1, 1, 2)
    @test collect(Oceananigans.Grids.znodes(grid, Face())) == [-120.0, -60.0, 0.0]
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

        arch = CPU()
        grid = @test_nowarn ImmersedBoundaryGrid(bathymetry_file, arch, (1, 1, 1))
        @test_nowarn top_bottom_boundary_conditions(; grid, bottom_drag_coefficient = 0.003)
        @test_nowarn MultiYearNORA3(nora3_filename, tmp)
    end
end

@testset "Bathymetry writer" begin
    mktempdir() do tmp
        arch = CPU()
        z_faces = [-20.0, -10.0, 0.0]
        grid = LatitudeLongitudeGrid(
            arch;
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
        reloaded = ImmersedBoundaryGrid(bathymetry_file, arch, (1, 1, 1))
        @test collect(Oceananigans.Grids.znodes(reloaded.underlying_grid, Face())) == z_faces

        @test FjordSim.Bathymetry.vertical_faces(grid) == z_faces
        # Independent of the halo size, which is what the old indexing got wrong.
        for halo_size in (1, 2, 3)
            deep_faces = [-450.0, -200.0, -50.0, 0.0]
            deep_grid = LatitudeLongitudeGrid(
                arch;
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
    arch = CPU()
    grid = LatitudeLongitudeGrid(
        arch;
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
    npoints = length(source_points)

    @test length(xs) == npoints + length(expected_contour_indices)
    @test xs[1:npoints] == first.(source_points)
    @test ys[1:npoints] == last.(source_points)
    @test bottom_heights[1:npoints] == -abs.(source_depths)
    @test xs[npoints+1:end] == first.(contour_vertices[expected_contour_indices.+1])
    @test ys[npoints+1:end] == last.(contour_vertices[expected_contour_indices.+1])
    @test all(==(-abs(contour_depth)), bottom_heights[npoints+1:end])

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
