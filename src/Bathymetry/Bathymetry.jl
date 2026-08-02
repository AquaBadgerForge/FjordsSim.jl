module Bathymetry

export DybdedataConfig,
    prepare_bathymetry,
    bathymetry_dataset,
    regrid_options,
    smoothing_options,
    geodatabase_path

using NCDatasets
using NumericalEarth
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain, znodes
using Printf: @printf
using Statistics: mean, median

using NumericalEarth.DataWrangling: Metadatum, metadata_path

using ..Configs: AbstractBathymetryConfig, FjordConfig, bathymetry_path
using ..Plotting: plot_bathymetry

# Matches the fixed loop count and neighbor threshold used in the Oslofjord notebook's
# post-regrid gap-filling pass.
const BATHYMETRY_GAP_FILL_PASSES = 10
const ISOLATED_SEA_CELL_LAND_SIDES = 3
# A spike needs neighbors to be compared against; two is the minimum that gives a meaningful
# median and keeps single-cell inlets, which have one wet neighbor, out of it.
const SPIKE_MIN_WET_NEIGHBOURS = 2
# `limit_bottom_slope` converges in tens of passes at the slope factors a fjord needs; the cap
# only exists so a pathological input cannot spin forever.
const SLOPE_LIMIT_MAX_PASSES = 400
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

"""
    smoothing_options(config)

Named tuple of keyword arguments forwarded to `smooth_bathymetry_gaps!` by `prepare_bathymetry`,
letting a source expose the post-regrid smoothing knobs as config fields. Defaults to no options,
which leaves only the topological cleanup every source gets — so implementing it is optional and
an existing source keeps its behavior.
"""
smoothing_options(config::AbstractBathymetryConfig) = (;)

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
    smooth_bathymetry_gaps!(bottom_height; smoothing_options(config)...)
    @info "Writing processed bathymetry file to $output_file"
    write_bathymetry_file(output_file, target_grid, bottom_height)
    @info "Finished preparing bathymetry"

    return (; dataset, raw_file = metadata_path(metadata), output_file, bottom_height)
end

"""
    prepare_bathymetry(config::FjordConfig)

Prepare the bathymetry a whole setup names: build its grid, regrid the source onto it, write the
processed NetCDF and the diagnostic plot.

This is the setup-level driver, the same shape as `download_forcing(config::FjordConfig)`. The
grid is built on the CPU because that is what the rest of the preparation pipeline reads, and
`bathymetry_path`'s directory is created here since this is the first step of a setup and nothing
has written into `data_root` yet.

# Returns
The `prepare_bathymetry(target_grid, config)` named tuple with `plot_file` added.
"""
function prepare_bathymetry(config::FjordConfig)
    grid = LatitudeLongitudeGrid(CPU(), config.grid_config)
    mkpath(dirname(bathymetry_path(config.bathymetry_config)))
    print_grid_extents(grid)

    result = prepare_bathymetry(grid, config.bathymetry_config)
    plot_file = plot_bathymetry(grid, result.bottom_height, config.bathymetry_config)

    @info "Raw bathymetry saved to $(result.raw_file)"
    @info "Processed FjordSim bathymetry saved to $(result.output_file)"
    @info "Bathymetry plot saved to $plot_file"

    return (; result..., plot_file)
end

"""
    print_grid_extents(grid)

Print the cell side lengths of `grid` in meters. The grid config states its extent in degrees, so
this is the only place the resulting resolution becomes visible before a long regrid starts.
"""
function print_grid_extents(grid)
    dx = xspacings(grid)
    dy = yspacings(grid)
    dz = zspacings(grid)

    println("Grid cell side extents (m):")
    @printf("  N = (%d, %d, %d)\n", grid.Nx, grid.Ny, grid.Nz)
    @printf("  Δx min/max/mean = %.3f / %.3f / %.3f\n", minimum(dx), maximum(dx), mean(dx))
    @printf("  Δy min/max/mean = %.3f / %.3f / %.3f\n", minimum(dy), maximum(dy), mean(dy))
    @printf("  Δz min/max/mean = %.3f / %.3f / %.3f\n", minimum(dz), maximum(dz), mean(dz))

    return nothing
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
    smooth_bathymetry_gaps!(bottom_height; spike_ratio = 0, max_slope_factor = 0, minimum_depth = 0)

Clean up the regridded bathymetry in place, in three stages.

First the topological pass every source gets: one diagonal-pair fill, then
`BATHYMETRY_GAP_FILL_PASSES` rounds of isolated sea/land cell cleanup, following the fixed cleanup
used for the Oslofjord ROMS-based bathymetry. This removes the checkerboard noise that would
otherwise skew the neighbor medians the next stage takes.

Then two optional stages, each skipped when its parameter is zero so a source that configures
neither keeps exactly the behavior above: `fill_shallow_spikes` with `spike_ratio`, and
`limit_bottom_slope` with `max_slope_factor` and `minimum_depth`. Despiking runs first, because a
spike is precisely the kind of one-cell feature that slope limiting would otherwise smear into its
neighbors instead of removing.
"""
function smooth_bathymetry_gaps!(
    bottom_height;
    spike_ratio = 0,
    max_slope_factor = 0,
    minimum_depth = 0,
)
    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    h = Array(interior(cpu_bottom_height, :, :, 1))

    h = fill_secondary_diagonal_pairs(fill_diagonal_pairs(h))
    for _ = 1:BATHYMETRY_GAP_FILL_PASSES
        h = fill_isolated_land_cells(remove_isolated_sea_cells(h))
    end

    spike_ratio > 0 && (h = fill_shallow_spikes(h; ratio = spike_ratio))
    max_slope_factor > 0 && (h = limit_bottom_slope(h; max_slope_factor, minimum_depth))

    set!(bottom_height, h)
    return bottom_height
end

"""
    fill_shallow_spikes(h; ratio, min_neighbours = SPIKE_MIN_WET_NEIGHBOURS)

Replace every sea cell shallower than `ratio` times the median depth of its sea neighbors with
that median, leaving land (`h >= 0`) untouched.

This targets the isolated shallow spikes regridding leaves along a coastline — and the ones a
`minimum_depth` floor creates, by lifting a sub-metre sliver to a constant depth that is still far
shallower than the water around it. Such a cell carries the same transport as its neighbors in a
fraction of the water column, so the velocity there grows until the run goes unstable. A genuinely
shallow *region* is not affected, because its neighbors are shallow too and the median moves with
it: this only fires where a cell disagrees with its surroundings.

A cell with fewer than `min_neighbours` sea neighbors is skipped, since a median over one value
would turn every single-cell inlet into its neighbor.
"""
function fill_shallow_spikes(h; ratio, min_neighbours = SPIKE_MIN_WET_NEIGHBOURS)
    filled = copy(h)
    Nx, Ny = size(h)

    for i = 1:Nx, j = 1:Ny
        h[i, j] < 0 || continue

        depths = eltype(h)[]
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = i + di, j + dj
            (1 <= ii <= Nx && 1 <= jj <= Ny) || continue
            h[ii, jj] < 0 && push!(depths, -h[ii, jj])
        end
        length(depths) >= min_neighbours || continue

        neighbour_depth = median(depths)
        -h[i, j] < ratio * neighbour_depth && (filled[i, j] = -neighbour_depth)
    end

    return filled
end

"""
    limit_bottom_slope(h; max_slope_factor, minimum_depth = 0, max_passes = SLOPE_LIMIT_MAX_PASSES)

Limit the bathymetry's steepness so that no two adjacent sea cells differ by more than
`max_slope_factor` in the Beckmann–Haidvogel slope parameter

```
r = |d₁ - d₂| / (d₁ + d₂)
```

where `d` is positive depth. Land is never touched, and no cell is made shallower than
`minimum_depth`.

`PartialCellBottom` bounds how *thin* a cell may be but says nothing about how much the depth may
change between adjacent columns, which is the quantity that destabilizes a regional run: a shallow
cell beside a deep one carries the same transport in a fraction of the water column. Bounding `r`
is the standard remedy.

Each offending pair is moved symmetrically — the deeper cell up and the shallower one down by the
same `Δ = (|d₁ - d₂| - max_slope_factor (d₁ + d₂)) / 2` — so the pair's depth sum is unchanged and
the domain's water volume is conserved rather than being quietly shifted. Passes repeat until no
pair exceeds the limit, since fixing one pair can push a neighbouring pair over it.
"""
function limit_bottom_slope(
    h;
    max_slope_factor,
    minimum_depth = 0,
    max_passes = SLOPE_LIMIT_MAX_PASSES,
)
    depth = map(value -> value < 0 ? -value : convert(eltype(h), NaN), h)
    Nx, Ny = size(h)

    for _ = 1:max_passes
        steepest = zero(eltype(h))

        for i = 1:Nx, j = 1:Ny
            isnan(depth[i, j]) && continue

            for (di, dj) in ((1, 0), (0, 1))
                ii, jj = i + di, j + dj
                (ii <= Nx && jj <= Ny) || continue
                isnan(depth[ii, jj]) && continue

                here, there = depth[i, j], depth[ii, jj]
                slope = abs(here - there) / (here + there)
                steepest = max(steepest, slope)
                slope > max_slope_factor || continue

                shift = (abs(here - there) - max_slope_factor * (here + there)) / 2
                deeper, shallower = here > there ? ((i, j), (ii, jj)) : ((ii, jj), (i, j))
                depth[deeper...] = max(depth[deeper...] - shift, minimum_depth)
                depth[shallower...] = max(depth[shallower...] + shift, minimum_depth)
            end
        end

        steepest <= max_slope_factor && break
    end

    return map((original, d) -> isnan(d) ? original : -d, h, depth)
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
