module Bathymetry

export DybdedataConfig,
    prepare_bathymetry,
    bathymetry_dataset,
    regrid_options,
    geodatabase_path

using NCDatasets
using NumericalEarth
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain, znodes

using NumericalEarth.DataWrangling: Metadatum, metadata_path

using ..Configs: AbstractBathymetryConfig, bathymetry_path

# Matches the fixed loop count and neighbor threshold used in the Oslofjord notebook's
# post-regrid gap-filling pass.
const BATHYMETRY_GAP_FILL_PASSES = 10
const ISOLATED_SEA_CELL_LAND_SIDES = 3
# NumericalEarth bathymetry regridding constructs a native grid with halo = (10, 10, 1).
# Keep a generated raw dataset comfortably larger than that minimum.
const MIN_NATIVE_BATHYMETRY_SIZE = 24

# --- Extension hooks ---

"""
    bathymetry_dataset(target_grid, config)

The NumericalEarth dataset `prepare_bathymetry` regrids onto `target_grid`, materializing it
locally if needed. This is the one source-specific step of the pipeline: everything after it
— regridding, smoothing and writing — is the same for every source.

A new bathymetry source implements this method on its `AbstractBathymetryConfig` subtype; see
`bathymetry_dataset(target_grid, config::DybdedataConfig)` in `src/Bathymetry/geonorge.jl`
for the built-in one.
"""
function bathymetry_dataset end

"""
    regrid_options(config)

Named tuple of keyword arguments forwarded to `NumericalEarth.regrid_bathymetry` by
`prepare_bathymetry`, letting a source expose its own regridding knobs as config fields.
Defaults to no options, so implementing it is optional.
"""
regrid_options(config::AbstractBathymetryConfig) = (;)

# --- Preparation pipeline ---

"""
    prepare_bathymetry(target_grid, config::AbstractBathymetryConfig; regrid_kw...)

Regrid the bathymetry source described by `config` onto `target_grid` with
`NumericalEarth.regrid_bathymetry`, and write a processed NetCDF file at
`bathymetry_path(config)` compatible with `FjordSim.Grids.ImmersedBoundaryGrid`.

The source enters only through `bathymetry_dataset(target_grid, config)`; regridding options
come from `regrid_options(config)`, with `regrid_kw...` overriding them per call.

# Returns
A named tuple with `dataset`, `raw_file`, `output_file`, and `bottom_height`.
"""
function prepare_bathymetry(target_grid, config::AbstractBathymetryConfig; regrid_kw...)
    dataset = bathymetry_dataset(target_grid, config)
    metadata = Metadatum(:bottom_height; dataset)

    @info "Regridding bathymetry onto target grid"
    options = (; regrid_options(config)..., regrid_kw...)
    bottom_height = NumericalEarth.regrid_bathymetry(target_grid, metadata; options...)

    output_file = bathymetry_path(config)
    @info "Smoothing small-scale bathymetry gaps"
    smooth_bathymetry_gaps!(bottom_height)
    @info "Writing processed bathymetry file to $output_file"
    write_bathymetry_file(output_file, target_grid, bottom_height)
    @info "Finished preparing bathymetry"

    return (; dataset, raw_file = metadata_path(metadata), output_file, bottom_height)
end

"""
    write_bathymetry_file(filepath, target_grid, bottom_height)

Write a processed NetCDF bathymetry file compatible with
`FjordSim.Grids.ImmersedBoundaryGrid`.
"""
function write_bathymetry_file(filepath::String, target_grid, bottom_height)
    isdir(dirname(filepath)) || mkpath(dirname(filepath))

    Nx, Ny, _ = size(target_grid)
    longitude = center_coordinates(x_domain(target_grid), Nx)
    latitude = center_coordinates(y_domain(target_grid), Ny)
    z_faces = vertical_faces(target_grid)

    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    h = Array(interior(cpu_bottom_height, :, :, 1))

    isfile(filepath) && rm(filepath; force = true)

    ds = NCDataset(filepath, "c")
    try
        defDim(ds, "lon", Nx)
        defDim(ds, "lat", Ny)
        defDim(ds, "zf", length(z_faces))

        longitude_variable = defVar(ds, "lon", Float64, ("lon",))
        latitude_variable = defVar(ds, "lat", Float64, ("lat",))
        z_faces_variable = defVar(ds, "z_faces", Float64, ("zf",))
        bottom_height_variable = defVar(ds, "h", Float32, ("lon", "lat"))

        longitude_variable[:] = longitude
        latitude_variable[:] = latitude
        z_faces_variable[:] = z_faces
        bottom_height_variable[:, :] = h
    finally
        close(ds)
    end

    return filepath
end

"""
    smooth_bathymetry_gaps!(bottom_height)

Remove small-scale diagonal checkerboard artifacts and single-cell sea/land noise left
by regridding, following the fixed cleanup pass used for the Oslofjord ROMS-based
bathymetry: one diagonal-pair fill, then `BATHYMETRY_GAP_FILL_PASSES` rounds of isolated
sea/land cell cleanup. Mutates `bottom_height` in place.
"""
function smooth_bathymetry_gaps!(bottom_height)
    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    h = Array(interior(cpu_bottom_height, :, :, 1))

    h = fill_secondary_diagonal_pairs(fill_diagonal_pairs(h))
    for _ = 1:BATHYMETRY_GAP_FILL_PASSES
        h = fill_isolated_land_cells(remove_isolated_sea_cells(h))
    end

    set!(bottom_height, h)
    return bottom_height
end

"""
    fill_diagonal_pairs(h)

For every 2x2 block where the top-left/bottom-right cells are sea (`h < 0`) and the
top-right/bottom-left cells are land (`h >= 0`), fill the top-right cell with the mean of
the two sea corners.
"""
function fill_diagonal_pairs(h)
    filled = copy(h)
    Nx, Ny = size(h)

    for i = 1:Nx-1, j = 1:Ny-1
        top_left = h[i, j]
        top_right = h[i, j+1]
        bottom_left = h[i+1, j]
        bottom_right = h[i+1, j+1]

        if top_left < 0 && bottom_right < 0 && top_right >= 0 && bottom_left >= 0
            filled[i, j+1] = (top_left + bottom_right) / 2
        end
    end

    return filled
end

"""
    fill_secondary_diagonal_pairs(h)

For every 2x2 block where the top-right/bottom-left cells are sea (`h < 0`) and the
top-left/bottom-right cells are land (`h >= 0`), fill the top-left cell with the mean of
the two sea corners.
"""
function fill_secondary_diagonal_pairs(h)
    filled = copy(h)
    Nx, Ny = size(h)

    for i = 1:Nx-1, j = 1:Ny-1
        top_left = h[i, j]
        top_right = h[i, j+1]
        bottom_left = h[i+1, j]
        bottom_right = h[i+1, j+1]

        if top_right < 0 && bottom_left < 0 && top_left >= 0 && bottom_right >= 0
            filled[i, j] = (top_right + bottom_left) / 2
        end
    end

    return filled
end

"""
    remove_isolated_sea_cells(h; sides = ISOLATED_SEA_CELL_LAND_SIDES)

Turn interior sea cells (`h < 0`) with at least `sides` land neighbors (out of
north/south/east/west) into land, by setting them to zero.
"""
function remove_isolated_sea_cells(h; sides = ISOLATED_SEA_CELL_LAND_SIDES)
    replaced = copy(h)
    Nx, Ny = size(h)

    for i = 2:Nx-1, j = 2:Ny-1
        if h[i, j] < 0
            land_neighbors = count((h[i-1, j] >= 0, h[i+1, j] >= 0, h[i, j-1] >= 0, h[i, j+1] >= 0))
            land_neighbors >= sides && (replaced[i, j] = zero(eltype(h)))
        end
    end

    return replaced
end

"""
    fill_isolated_land_cells(h)

Fill interior land cells (`h >= 0`) whose north/south/east/west neighbors are all sea
(`h < 0`) with the mean of those 4 neighbor depths.
"""
function fill_isolated_land_cells(h)
    filled = copy(h)
    Nx, Ny = size(h)

    for i = 2:Nx-1, j = 2:Ny-1
        if h[i, j] >= 0
            west = h[i, j-1]
            north = h[i-1, j]
            east = h[i, j+1]
            south = h[i+1, j]

            if west < 0 && north < 0 && east < 0 && south < 0
                filled[i, j] = (west + north + east + south) / 4
            end
        end
    end

    return filled
end

# --- Grid domain helpers ---

function center_coordinates(domain, N)
    delta = domain_step(domain, N)
    return collect(range(domain[1] + delta / 2, step = delta, length = N))
end

domain_step(domain, N) = (domain[2] - domain[1]) / N

"""
    vertical_faces(grid)

The `Nz + 1` interior vertical face coordinates of `grid`, bottom to top.

Uses `znodes` rather than slicing `grid.z.cᵃᵃᶠ` by hand: that array is an `OffsetArray` whose
indices already account for the halo, so offsetting the start by `grid.Hz` skipped the
deepest face and appended one above the surface.
"""
vertical_faces(grid) = collect(znodes(grid, Face()))

function expand_domain(domain, N, padding_cells)
    delta = domain_step(domain, N)
    lower = domain[1] - padding_cells * delta
    upper = domain[2] + padding_cells * delta
    return (lower, upper)
end

include("geonorge.jl")

end  # module Bathymetry
