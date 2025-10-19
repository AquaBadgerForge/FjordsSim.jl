using Oceananigans.Grids: LatitudeLongitudeGrid, ImmersedBoundaryGrid, Center
using Oceananigans.Fields: Field, set!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: GridFittedBottom
using JLD2: @load
using NCDatasets: Dataset

function compute_faces(centers)
    spacing = diff(centers)[1]  # Assuming uniform spacing
    faces = vcat([centers[1] - spacing / 2], (centers[1:end-1] .+ centers[2:end]) / 2, [centers[end] + spacing / 2])
    return faces
end

"""
This grid maker uses a netcdf file to create a grid.
A netcdf file should have 4 variables:
- z_faces - 1d
- h - 2d
- lat and lon - 1d
"""
function grid_from_nc(arch, halo, filepath)
    ds = Dataset(filepath)
    z_faces = ds["z_faces"][:]
    # depths from an nc file are for grid centers
    depth = ds["h"][:, :]
    latitude = compute_faces(ds["lat"][:])
    longitude = compute_faces(ds["lon"][:])

    Nx, Ny = size(depth)
    Nz = length(z_faces)
    # Size should be for grid centers,
    # but z, latitude and langitude should be for faces
    underlying_grid =
        LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz - 1), halo = halo, z = z_faces, latitude, longitude)
    bathymetry = Field{Center, Center, Nothing}(underlying_grid)
    set!(bathymetry, coalesce.(depth, 0.0))
    fill_halo_regions!(bathymetry)
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry); active_cells_map = true)
    return grid
end

"""Takes depths from a jld2 file and preset lats and lon, return a grid."""
function grid_from_bathymetry_file(arch, halo, filepath, latitude, longitude)
    @load filepath depth z_faces
    Nx, Ny = size(depth)
    # Sanitize z_faces: finite, sorted, strictly increasing
    zf = collect(Float64.(z_faces))
    zf = filter(isfinite, zf)
    sort!(zf)
    # Enforce strict monotonicity by nudging duplicates by a tiny epsilon
    eps = 1e-6
    for n in 2:length(zf)
        if zf[n] <= zf[n-1]
            zf[n] = zf[n-1] + eps
        end
    end
    Nz = length(zf) - 1
    underlying_grid = LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz), halo = halo, z = zf, latitude, longitude)
    bathymetry = Field{Center, Center, Nothing}(underlying_grid)
    # Ensure bathymetry has no NaNs/Infs and is non-negative
    _depth = map(x -> (isfinite(x) && x > 0) ? x : 0.0, depth)
    set!(bathymetry, _depth)
    fill_halo_regions!(bathymetry)
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry); active_cells_map = true)
    return grid
end
