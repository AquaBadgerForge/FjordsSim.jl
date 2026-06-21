using FjordSim
using FjordSim.Bathymetry: write_bathymetry_file
using Test
using NCDatasets
using Oceananigans
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
import Oceananigans: iteration
import Oceananigans.Fields: interior
import Oceananigans.Utils: prettytime

struct MockProgressField
    values
end

interior(field::MockProgressField) = field.values

struct MockProgressSimulation
    model
    Δt
end

iteration(::MockProgressSimulation) = 7
prettytime(::MockProgressSimulation) = "12 hours"

function mock_progress_simulation()
    u = MockProgressField([1.0, -2.0])
    v = MockProgressField([0.5, -0.25])
    w = MockProgressField([0.01, -0.02])
    T = MockProgressField([3.0, 5.0])
    NUT = MockProgressField([1.0, 2.5])
    HET = MockProgressField([0.1, 0.4])

    ocean_model = (
        velocities = (u, v, w),
        tracers = (T = T, NUT = NUT, HET = HET),
    )

    return MockProgressSimulation((ocean = (model = ocean_model,),), 1.0)
end

@testset "Backward Compatibility — API Exports" begin
    # Verify all exported symbols are present in the public interface
    exported_symbols = [
        :ImmersedBoundaryGrid,
        :forcing_from_file,
        :top_bottom_boundary_conditions,
        :coupled_hydrostatic_simulation,
        :recursive_merge,
        :progress,
        :progress_callback,
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
        (:Utils, :progress_callback),
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

@testset "Configurable progress logging" begin
    sim = mock_progress_simulation()

    @test_logs (:info, r"max\|u\|: .*extrema\(T\): \(5\.00, 3\.00\) ᵒC, extrema\(NUT\): \(2\.50, 1\.00\), wall time:") progress(
        sim;
        tracers = (:T, :NUT),
        wall_time_reference = Ref(time_ns()),
    )

    callback = progress_callback(; tracers = (:HET,))
    @test callback isa Function
    @test_logs (:info, r"extrema\(HET\): \(0\.40, 0\.10\), wall time:") callback(sim)

    @test_throws ArgumentError progress(sim; tracers = (:POM,), wall_time_reference = Ref(time_ns()))
end

@testset "Backward Compatibility — Function Signatures" begin
    mktempdir() do tmp
        # Create minimal test files
        bathymetry_path = joinpath(tmp, "bathymetry.nc")
        ds = NCDataset(bathymetry_path, "c")
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
        grid = ImmersedBoundaryGrid(bathymetry_path, arch, (1, 1, 1))

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

@testset "FjordSim.jl" begin
    mktempdir() do tmp
        bathymetry_path = joinpath(tmp, "bathymetry.nc")
        nora3_filename = "NORA3_test.nc"
        nora3_path = joinpath(tmp, nora3_filename)

        # Minimal bathymetry file for Grids.ImmersedBoundaryGrid.
        ds = NCDataset(bathymetry_path, "c")
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
        grid = @test_nowarn ImmersedBoundaryGrid(bathymetry_path, arch, (1, 1, 1))
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

        bathymetry_path = joinpath(tmp, "bathymetry_written.nc")
        @test_nowarn write_bathymetry_file(bathymetry_path, grid, bottom_height)

        ds = NCDataset(bathymetry_path)
        @test ds["lon"][:] == [10.5, 11.5]
        @test ds["lat"][:] == [59.5, 60.5, 61.5]
        @test ds["h"][:, :] == Float32[-15.0 -16.0 -17.0; -18.0 -19.0 -20.0]
        close(ds)

        @test_nowarn ImmersedBoundaryGrid(bathymetry_path, arch, (1, 1, 1))
        @test FjordSim.Bathymetry.contour_point_indices(10, 3) == [0, 3, 6, 9]
        @test FjordSim.Bathymetry.contour_point_indices(5, 3) == [0, 3, 4]
    end
end
