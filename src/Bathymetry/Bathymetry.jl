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

using ..Configs:
    AbstractBathymetryConfig,
    FjordConfig,
    bathymetry_path,
    domain_grid,
    open_edges,
    lateral_edges,
    LATERAL_EDGES
using ..Plotting: plot_bathymetry

# Matches the fixed loop count and neighbor threshold used in the Oslofjord notebook's
# post-regrid gap-filling pass.
const BATHYMETRY_GAP_FILL_PASSES = 10
const ISOLATED_SEA_CELL_LAND_SIDES = 3
# A spike needs neighbors to be compared against; two is the minimum that gives a meaningful
# median and keeps single-cell inlets, which have one wet neighbor, out of it.
const SPIKE_MIN_WET_NEIGHBOURS = 2
# `limit_bottom_slope` reaches its target in tens of passes at the slope factors a fjord needs, but
# its exact convergence test can keep sweeping after that (see the function), so the cap is what
# stops it in practice and is set high enough that the extra sweeps are affordable.
const SLOPE_LIMIT_MAX_PASSES = 1000
# How far above `max_slope_factor` the achieved slope may sit before it is a real failure rather than
# `Float32` rounding on a pair the limiter has driven exactly onto the limit.
const SLOPE_LIMIT_TOLERANCE = 1e-3
# The element type `write_bathymetry_file` stores bottom heights as, and therefore the precision the
# simulation actually runs on: `Grids.ImmersedBoundaryGrid` reads the file back and immerses what it
# finds there. Stated once because `snap_partial_bottom_cells` has to reason in it — see
# `snap_to_face` — and a snap that only holds in the pipeline's own precision does not hold.
const BATHYMETRY_ELTYPE = Float32
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

`edges` names the domain's open lateral boundaries, and reaches `smooth_bathymetry_gaps!` for the
one stage that acts on them. It is a keyword rather than a field of the bathymetry config for the
same reason `prepare_forcing` takes one: which edges are open is a property of the domain and its
exterior data, not of a bathymetry source. `nothing` means none, and that stage is then a no-op.

# Returns
A named tuple with `dataset`, `raw_file`, `output_file`, and `bottom_height`.
"""
function prepare_bathymetry(target_grid, config::AbstractBathymetryConfig; edges = nothing, regrid_kw...)
    dataset = bathymetry_dataset(target_grid, config)
    metadata = Metadatum(:bottom_height; dataset)

    @info "Regridding bathymetry onto target grid"
    options = (; regrid_options(config)..., regrid_kw...)
    bottom_height = NumericalEarth.regrid_bathymetry(target_grid, metadata; options...)

    output_file = bathymetry_path(config)
    @info "Smoothing small-scale bathymetry gaps"
    smooth_bathymetry_gaps!(bottom_height; smoothing_options(config)..., edges)
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

The open edges come from `open_edges(config.boundary_config)` — the one place a setup states them —
so a setup naming no boundary data config prepares its bathymetry with every lateral boundary
treated as a wall, which is the right reading when there is no exterior state to admit.

# Returns
The `prepare_bathymetry(target_grid, config)` named tuple with `plot_file` added.
"""
function prepare_bathymetry(config::FjordConfig)
    grid = domain_grid(config.grid_config, CPU())
    mkpath(dirname(bathymetry_path(config.bathymetry_config)))
    print_grid_extents(grid)

    result = prepare_bathymetry(grid, config.bathymetry_config; edges = open_edges(config.boundary_config))
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
        bottom_height_variable = defVar(ds, "h", BATHYMETRY_ELTYPE, ("lon", "lat"))

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
    smooth_bathymetry_gaps!(bottom_height; open_boundary_land_cells = 0, max_island_cells = 0,
                            close_narrow_passages = false, spike_ratio = 0,
                            minimum_cell_fraction = 0, max_slope_factor = 0, minimum_depth = 0,
                            edges = nothing)

Clean up the regridded bathymetry in place, in seven stages.

First the topological pass every source gets: one diagonal-pair fill, then
`BATHYMETRY_GAP_FILL_PASSES` rounds of isolated sea/land cell cleanup, following the fixed cleanup
used for the Oslofjord ROMS-based bathymetry. This removes the checkerboard noise that would
otherwise skew the neighbor medians a later stage takes.

Then six optional stages — seven runs, since `snap_partial_bottom_cells` runs once each side of
slope limiting — each skipped when its parameter is `false` or zero so a source that configures none
keeps exactly the behavior above: `clear_open_boundary_land` with
`open_boundary_land_cells` and `edges`, `fill_small_islands` with `max_island_cells`,
`remove_narrow_passages` with `close_narrow_passages`, `fill_shallow_spikes` with `spike_ratio`,
`snap_partial_bottom_cells` with `minimum_cell_fraction`, and `limit_bottom_slope` with
`max_slope_factor` and `minimum_depth`.

`clear_open_boundary_land` is first of the six, and belongs there for the same reason
`fill_small_islands` precedes `remove_narrow_passages`: it only ever turns land into sea, so it can
neither manufacture a narrow passage nor strand a cell, while running it first lets the island and
passage stages judge the band as it will actually be. A land component straddling the band's inner
edge loses its in-band cells here and is then correctly measured — as whatever is left of it — by
`fill_small_islands`. It is also the only stage that reads `edges`, and it is a no-op when the setup
names none.

The order of the remaining five is forced. The two stages that change the *land mask* run before the two that
change depths: run the other way round, the cells they are about to turn into land or water would
already have contributed to a neighbor median and to a slope pair. Of the two land-mask stages,
`fill_small_islands` runs first, because flooding an island only ever widens water and so cannot
create a passage, while closing a passage adds land exactly where a spurious loop closes — which is
where these islands live — and can bridge one to the mainland, hiding it from the size test for
good. Despiking then precedes slope limiting, because a spike is precisely the kind of one-cell
feature that slope limiting would otherwise smear into its neighbors instead of removing.

`snap_partial_bottom_cells` runs **twice**, once each side of slope limiting, and that is not
redundancy. It is the only stage that knows where the vertical faces are, and slope limiting moves
depths without regard for them, so the limiter is guaranteed to put back slivers the snap removed —
the tighter the limit, the more of them. Measured on `oslofjorden` at `max_slope_factor = 0.25`, one
snap before the limiter leaves 2628 slivers (6.0% of bottom cells) in the finished field; snapping
again afterwards leaves **none**. The price is that the second snap breaks the limit it was just
handed: the steepest slope goes from 0.250 to 0.307, on 16 pairs out of 84 000. That is the whole
trade — a slope parameter 23% over its target on a handful of pairs, against six percent of the
domain's bottom cells being slivers.

The first snap is not made redundant by the second. It hands the limiter a field already aligned to
the faces, so the limiter has less to move and the second snap has less to correct; and where the
limiter does not run at all (`max_slope_factor = 0`) it is the only one there is.

Closing a passage leaves a dead-end stub behind — the cell on the far side of it now has three land
neighbors — so the isolated-cell cleanup loop runs again afterwards. That loop cannot re-open a
closed passage: a passage cell has two land neighbors by definition and they stay land, so
`fill_isolated_land_cells`, which needs all four neighbors wet, can never fire on it.

Flooding an island needs no such cleanup, which is why it has none: every neighbour of a flooded
component is sea by maximality, so no flooded cell borders land and neither isolated-cell stage has
anything new to fire on.
"""
function smooth_bathymetry_gaps!(
    bottom_height;
    open_boundary_land_cells = 0,
    max_island_cells = 0,
    close_narrow_passages = false,
    spike_ratio = 0,
    minimum_cell_fraction = 0,
    max_slope_factor = 0,
    minimum_depth = 0,
    edges = nothing,
)
    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    h = Array(interior(cpu_bottom_height, :, :, 1))

    h = fill_secondary_diagonal_pairs(fill_diagonal_pairs(h))
    for _ = 1:BATHYMETRY_GAP_FILL_PASSES
        h = fill_isolated_land_cells(remove_isolated_sea_cells(h))
    end

    if open_boundary_land_cells > 0
        wet_before = count(<(0), h)
        for edge in lateral_edges(edges)
            h = clear_open_boundary_land(h, Val(edge); cells = open_boundary_land_cells)
        end
        @info "Flooded $(count(<(0), h) - wet_before) land cells within $open_boundary_land_cells cells of the open boundary"
    end

    if max_island_cells > 0
        wet_before = count(<(0), h)
        h = fill_small_islands(h; max_cells = max_island_cells)
        @info "Flooded $(count(<(0), h) - wet_before) cells of interior land islands of $max_island_cells cells or fewer"
    end

    if close_narrow_passages
        wet_before = count(<(0), h)
        h = remove_narrow_passages(h)
        for _ = 1:BATHYMETRY_GAP_FILL_PASSES
            h = fill_isolated_land_cells(remove_isolated_sea_cells(h))
        end
        @info "Closed $(wet_before - count(<(0), h)) cells of one-cell-wide passages and their stubs"
    end

    spike_ratio > 0 && (h = fill_shallow_spikes(h; ratio = spike_ratio))

    z_faces = minimum_cell_fraction > 0 ? vertical_faces(bottom_height.grid) : nothing

    minimum_cell_fraction > 0 &&
        (h = snap_and_report(h, z_faces, minimum_cell_fraction, minimum_depth, "before slope limiting"))

    if max_slope_factor > 0
        h = limit_bottom_slope(h; max_slope_factor, minimum_depth)
        # Slope limiting moves depths without regard for where the vertical faces are, so it puts
        # back the very slivers the snap above removed — and the tighter the limit, the more of
        # them. Snapping again is what makes the two stages compose instead of one undoing the
        # other; see the docstring for what it costs.
        minimum_cell_fraction > 0 &&
            (h = snap_and_report(h, z_faces, minimum_cell_fraction, minimum_depth, "after slope limiting"))
    end

    set!(bottom_height, h)
    return bottom_height
end

"""
    snap_and_report(h, z_faces, minimum_cell_fraction, minimum_depth, when)

`snap_partial_bottom_cells`, plus the log line saying how many columns it moved and `when` it ran.

A named helper because the stage runs twice — once before slope limiting and once after — and the
two counts together are what tells a reader how much of the first pass the limiter undid.
"""
function snap_and_report(h, z_faces, minimum_cell_fraction, minimum_depth, when)
    snapped = snap_partial_bottom_cells(h, z_faces; minimum_cell_fraction, minimum_depth)
    @info "Raised $(count(!=(0), snapped - h)) columns off a bottom cell thinner than $minimum_cell_fraction of its layer ($when)"
    return snapped
end

"""
    clear_open_boundary_land(h, edge::Val; cells)

Flood every land cell within `cells` rows of one open lateral boundary, giving each the mean depth
of its wet neighbours.

A coastline that runs into an open boundary is the worst place in the domain for one to be. The
boundary condition there has to reconcile a radiation scheme, a prescribed exterior state and a wall
within a cell or two of each other, and a headland poking through the boundary row splits the
prescribed inflow around an obstacle the exterior dataset never saw. Clearing a band of it leaves the
open edge a clean, uninterrupted channel, which is the geometry every open-boundary scheme is derived
for.

The cost is that the newly wet cells have no exterior profile of their own — the prepared boundary
file marks them dry, and `FjordSim.Forcing.fill_boundary_gaps!` fills them from the nearest wet cell
along the boundary. That is the same treatment a boundary column already gets wherever the model and
the source disagree about the coastline, so the trade is a few laterally interpolated columns against
a headland in the boundary row. Keep `cells` small enough that it stays a few.

Depths are assigned by repeated relaxation rather than in one pass: each round fills the cells that
have at least one wet neighbour *now*, so a cell in the middle of a cleared headland inherits from
the coast through the cells between it, and the band ramps rather than stepping. A cell no round can
reach is left as land, which can only happen if the band holds a land component touching nothing wet.

# Arguments
- `h`: bottom height, `h < 0` sea.
- `edge`: `Val(:south)`, `Val(:north)`, `Val(:west)` or `Val(:east)`.
- `cells`: how many rows in from that edge to clear.
"""
function clear_open_boundary_land(h, edge::Val; cells)
    Nx, Ny = size(h)
    cleared = copy(h)
    band = open_boundary_band(edge, Nx, Ny, cells)

    pending = [(i, j) for j in band[2], i in band[1] if h[i, j] >= 0]
    isempty(pending) && return cleared

    # One round per row of the band is enough for the innermost cell to reach the coast, and the
    # loop breaks as soon as a round fills nothing, so the bound only caps a band that cannot fill.
    for _ = 1:(cells+1)
        isempty(pending) && break

        filled = Tuple{Int,Int}[]
        depths = eltype(h)[]
        for (i, j) in pending
            neighbours = eltype(h)[]
            for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                ii, jj = i + di, j + dj
                (1 <= ii <= Nx && 1 <= jj <= Ny) || continue
                cleared[ii, jj] < 0 && push!(neighbours, cleared[ii, jj])
            end
            isempty(neighbours) && continue
            push!(filled, (i, j))
            push!(depths, mean(neighbours))
        end
        isempty(filled) && break

        # Written after the whole round is computed, so every cell in one round sees the same
        # field and the fill does not depend on the order the band is walked.
        for (n, cell) in enumerate(filled)
            cleared[cell...] = depths[n]
        end
        pending = filter(cell -> cleared[cell...] >= 0, pending)
    end

    return cleared
end

"""
    open_boundary_band(edge::Val, Nx, Ny, cells)

The `(i_range, j_range)` of the `cells` rows nearest one lateral edge, clamped to the domain.

One method per edge, so the four are independent statements rather than branches, and an edge that
is not one of `LATERAL_EDGES` raises rather than silently selecting nothing.
"""
open_boundary_band(::Val{:south}, Nx, Ny, cells) = (1:Nx, 1:min(cells, Ny))
open_boundary_band(::Val{:north}, Nx, Ny, cells) = (1:Nx, max(1, Ny - cells + 1):Ny)
open_boundary_band(::Val{:west}, Nx, Ny, cells) = (1:min(cells, Nx), 1:Ny)
open_boundary_band(::Val{:east}, Nx, Ny, cells) = (max(1, Nx - cells + 1):Nx, 1:Ny)

open_boundary_band(::Val{edge}, Nx, Ny, cells) where {edge} =
    throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))

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
    snap_to_face(FT, face)

`face` as an `FT`, rounded *away* from the water: the smallest `FT` value at or above it.

`FT` here is `BATHYMETRY_ELTYPE`, the precision the *file* stores, not the precision the pipeline
computes in — and that distinction is the whole point. `smooth_bathymetry_gaps!` works in the grid's
float type, usually `Float64`, where `face_above` lands on the face exactly; the narrowing to
`Float32` happens later, in `write_bathymetry_file`, on a value the snap has already stopped looking
at. So a snap that is exact when it is made can be undone by the write.

It is undone whenever a face is not representable in `Float32`. Four of `oslofjorden`'s are not —
−10.8, −7.9, −3.7 and −2.2 all narrow to a hair *below* the face they name. The stored column then
ends an infinitesimal distance under the face, so the layer beneath is its bottom cell at about
5 × 10⁻⁸ of full thickness: the pipeline writes out, in the worst possible form, exactly the sliver
it removed. Measured on `oslofjorden`, 1265 columns, 2.9 % of the domain, every one of them at those
four faces.

Rounding up ends the column a nanometre *above* the face instead, which survives the narrowing and
reads back as a full bottom cell. The error introduced is at the seventh significant figure of a
depth in metres.

The old 18-level `z_faces` were every one exactly representable in `Float32`, which is why this never
bit before, and why it is stated here rather than left as a silent guard: the next setup to write a
face like −10.8 would hit it again.
"""
@inline function snap_to_face(FT, face)
    rounded = convert(FT, face)
    return rounded < face ? nextfloat(rounded) : rounded
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

    passes = 0
    for pass = 1:max_passes
        passes = pass
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

    # The loop's own test is measured before each pass's fixes and compares exactly, so a field whose
    # steepest pair has been driven *onto* the limit can still round marginally above it in `Float32`
    # and never register a clean pass. Report what the field actually achieved rather than what the
    # loop concluded, and say so when the cap was what stopped it.
    achieved = steepest_slope(depth)
    if passes == max_passes
        @info "limit_bottom_slope stopped at the $max_passes-pass cap with steepest slope $achieved (limit $max_slope_factor)"
        achieved > max_slope_factor * (1 + SLOPE_LIMIT_TOLERANCE) && @warn(
            "Bathymetry slope limiting did not converge: steepest slope is $achieved against a " *
            "limit of $max_slope_factor. Raise `max_passes` or loosen `max_slope_factor`."
        )
    else
        @info "limit_bottom_slope converged in $passes passes with steepest slope $achieved (limit $max_slope_factor)"
    end

    return map((original, d) -> isnan(d) ? original : -d, h, depth)
end

"""
    steepest_slope(depth)

The largest `r = |d1 - d2| / (d1 + d2)` over adjacent wet pairs of a `depth` field whose land is
`NaN`. What `limit_bottom_slope` actually achieved, as opposed to what it measured on its way there.
"""
function steepest_slope(depth)
    Nx, Ny = size(depth)
    steepest = zero(eltype(depth))

    for i = 1:Nx, j = 1:Ny
        isnan(depth[i, j]) && continue

        for (di, dj) in ((1, 0), (0, 1))
            ii, jj = i + di, j + dj
            (ii <= Nx && jj <= Ny) || continue
            isnan(depth[ii, jj]) && continue

            here, there = depth[i, j], depth[ii, jj]
            steepest = max(steepest, abs(here - there) / (here + there))
        end
    end

    return steepest
end

"""
    snap_partial_bottom_cells(h, z_faces; minimum_cell_fraction, minimum_depth = 0)

Raise every wet column whose bottom partial cell would be thinner than `minimum_cell_fraction` of
its layer to the vertical face above, so the sliver is never created.

`PartialCellBottom` refuses to make a bottom cell thinner than `minimum_fractional_cell_height`
times the layer thickness: its kernel takes `min(z⁺ - ϵ Δz, zb)`, pushing the bottom *down* until
the cell is exactly that thick. So a sounding lying just below a face does not give a thin cell —
it gives a cell of exactly `ϵ Δz` whose floor is somewhere the sounding never was. This stage moves
the bottom the other way instead, up to the face, which ends the column one layer higher and leaves
every remaining bottom cell a full one. `minimum_cell_fraction` must therefore be the *same* number
the grid gives `PartialCellBottom` — 0.2, Oceananigans' default, for the grids `Grids.jl` builds.

Such a cell is dangerous out of proportion to its size, and for tracers rather than for momentum.
It holds a fraction of the water its neighbours do, so any flux into it moves its concentration by a
large amount, and because it reaches into a layer its neighbours may not reach at all it can be
nearly cut off horizontally as well. On `oslofjorden` the two worst were at the *open southern
boundary*: soundings of 51.3 m and 54.5 m against a layer spanning 50 to 75 m, floored to fractions
of 0.052 and 0.179, with one and two lateral neighbours at their own level because the columns
beside them (44.1 m, 48.5 m) stop a whole layer higher. Salinity there climbed from 33 to 65 psu and
temperature from 5 to 28 °C over four and a half days, while the prepared boundary file was asking
for 33.2 psu and 8.1 °C and the cell was inflowing 70% of the time — the nudging was active
throughout and simply could not keep up. Across the domain the tracer overshoots were six times
over-represented in bottom cells with at most one lateral neighbour.

Refining `z_faces` is the obvious alternative and it does not substitute: a finer grid gives the
seabed more faces to cross, so the count of laterally isolated bottom cells *rises* — measured on
`oslofjorden`, 3.9% of columns on the 18-level grid against 5.0% on the 24-level one that replaced
it. It is the sliver that has to go, not the layer that has to shrink. The two are complementary,
and the finer grid makes this stage matter slightly more rather than less.

Called **twice** by `smooth_bathymetry_gaps!`, once each side of slope limiting; see there for why.

The value written is `snap_to_face(BATHYMETRY_ELTYPE, face_above)`, not `face_above` itself, because
the snap has to survive being narrowed to the file's precision — see `snap_to_face`.

Two guards. A column is left alone when the face above is the surface or when snapping to it would
breach `minimum_depth`, so the stage never dries a cell out or undercuts the depth floor; such a
column keeps its sliver and `PartialCellBottom` floors it as before. On `oslofjorden` neither guard
fires, and with both passes in place `PartialCellBottom`'s clamp does not fire anywhere in the
finished field.
"""
function snap_partial_bottom_cells(h, z_faces; minimum_cell_fraction, minimum_depth = 0)
    snapped = copy(h)
    Nz = length(z_faces) - 1
    Nx, Ny = size(h)

    for i = 1:Nx, j = 1:Ny
        h[i, j] < 0 || continue

        layer = findlast(face -> face <= h[i, j], z_faces)
        (layer === nothing || layer > Nz) && continue

        face_above = z_faces[layer+1]
        fraction = (face_above - h[i, j]) / (face_above - z_faces[layer])
        fraction < minimum_cell_fraction || continue
        (face_above < 0 && -face_above >= minimum_depth) || continue

        snapped[i, j] = oftype(h[i, j], snap_to_face(BATHYMETRY_ELTYPE, face_above))
    end

    return snapped
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
    wet_component(h, seed, blocked)

The set of sea cells (`h < 0`) reachable from `seed` by north/south/east/west steps, treating
`blocked` as land, as a `BitMatrix` the size of `h`.

A hand-rolled flood fill rather than `ImageMorphology.label_components`, which
`NumericalEarth.remove_minor_basins!` uses: that package reaches FjordSim only transitively, and
one breadth-first walk is all `remove_narrow_passages` needs.
"""
function wet_component(h, seed, blocked)
    Nx, Ny = size(h)
    seen = falses(Nx, Ny)
    seen[seed...] = true
    stack = [seed]

    while !isempty(stack)
        i, j = pop!(stack)
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = i + di, j + dj
            (1 <= ii <= Nx && 1 <= jj <= Ny) || continue
            (seen[ii, jj] || (ii, jj) == blocked) && continue
            h[ii, jj] < 0 || continue
            seen[ii, jj] = true
            push!(stack, (ii, jj))
        end
    end

    return seen
end

"""
    remove_narrow_passages(h)

Turn every one-cell-wide sea passage whose removal leaves both of its sides connected into land.

A one-cell-wide passage is a sea cell that is sea on both sides along one axis and land on both
sides along the other — regridding leaves these where a channel too narrow to resolve cuts through
a peninsula. They are what destabilizes a regional run once the boundary admits a tide: the two
basins such a cell joins are usually already connected elsewhere, so it closes a loop, and the
barotropic head difference around that loop is pushed through a cross-section of one cell width and
a few metres depth. Velocity there grows until the time-step wizard collapses. Neither
`fill_shallow_spikes` nor `limit_bottom_slope` can help, because both bound depth *contrast* and
the passage agrees with its neighbors — it is its *width* that is wrong.

A passage that is the sole link to a basin is kept: closing it would delete that water from the
domain. Candidates are therefore tested one at a time against the partially closed field rather than
all at once against the input, which is what makes the stage unable to sever a basin — two passages
that are each redundant while the other is open would together disconnect one if closed as a batch.

One pass, deliberately: closing a passage can leave a neighbouring cell one-cell-wide in turn, and
iterating to convergence would erode a genuine narrow arm cell by cell.
"""
function remove_narrow_passages(h)
    closed = copy(h)
    Nx, Ny = size(h)

    for i = 2:Nx-1, j = 2:Ny-1
        closed[i, j] < 0 || continue

        west, east = closed[i-1, j] < 0, closed[i+1, j] < 0
        south, north = closed[i, j-1] < 0, closed[i, j+1] < 0

        sides = if west && east && !south && !north
            ((i - 1, j), (i + 1, j))
        elseif south && north && !west && !east
            ((i, j - 1), (i, j + 1))
        else
            continue
        end

        # Redundant only if the far side is still reachable with this cell blocked.
        wet_component(closed, sides[1], (i, j))[sides[2]...] || continue
        closed[i, j] = zero(eltype(h))
    end

    return closed
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

"""
    land_component!(visited, h, seed)

The land cells (`h >= 0`) reachable from `seed` by north/south/east/west steps, as a vector of
`(i, j)`, marking each of them in `visited`.

The companion of `wet_component` for the other phase, kept separate rather than folded into it:
that one answers a connectivity question about one candidate cell, takes a `blocked` cell for it
and returns a full-domain mask, while this one is swept over the whole field. So it needs a
`visited` buffer shared across components — the thing that makes the sweep walk every land cell
exactly once, whatever the component count — and the membership list itself, which is what the
caller fills and measures.
"""
function land_component!(visited, h, seed)
    Nx, Ny = size(h)
    visited[seed...] = true
    cells = [seed]
    stack = [seed]

    while !isempty(stack)
        i, j = pop!(stack)
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = i + di, j + dj
            (1 <= ii <= Nx && 1 <= jj <= Ny) || continue
            (visited[ii, jj] || h[ii, jj] < 0) && continue
            visited[ii, jj] = true
            push!(cells, (ii, jj))
            push!(stack, (ii, jj))
        end
    end

    return cells
end

"""
    fill_small_islands(h; max_cells)

Flood every 4-connected patch of land (`h >= 0`) of at most `max_cells` cells that does not touch
the domain edge, filling all of its cells with the mean depth of the sea around it.

This is the dual of `remove_narrow_passages` — an unresolved *land* island rather than an
unresolved water channel — and it destabilizes a run the same way. Flow splits around a rock a few
cells across, so the head difference between its two ends is worked out around a closed loop whose
arms are each a handful of cells wide: the same loop geometry, drawn by the land instead of by the
water. Nothing about the flow such an island obstructs is resolved by one cell of it, so keeping it
buys no fidelity.

On `oslofjorden` the domain's velocity maximum — 2.67 m s⁻¹, near uniform over the whole water
column — sat on the circulation around a three-cell, one-cell-wide island at i = 83, j = 104-106.
Of that bathymetry's 85 land components, 46 are interior clusters of 2 to 6 cells, 153 cells in
all. Neither `fill_shallow_spikes` nor `limit_bottom_slope` can see them, for the same reason
neither sees a narrow passage: both bound depth *contrast*, and the defect is geometry.

`fill_isolated_land_cells` is the `max_cells = 1` case of this and runs on every source; this
generalizes it to the clusters that survive it. A three-cell ridge survives because every one of
its cells has a land neighbour, so none of them has the four wet neighbours that stage needs.

A component touching the domain edge is kept whatever its size, because land continuing outside the
domain may have only a few cells inside it, and flooding those would open the domain into water
that is not there. This is the component-wise form of the `2:Nx-1, 2:Ny-1` bound the single-cell
stages use.

The patch gets one depth, the mean over the sea cells orthogonally adjacent to it: a cell in the
middle of a 2x3 patch has no sea neighbour of its own, and a flat floor over a patch this small is
what slope limiting would produce anyway. A sea cell touching the patch on two sides counts twice,
weighting the mean by contact length. For a single cell this is exactly the value
`fill_isolated_land_cells` computes.

No cleanup pass is needed afterwards. Every neighbour of a flooded component is sea — a land
neighbour would by maximality be part of the component — so no flooded cell borders land, and
neither `remove_isolated_sea_cells` nor `fill_isolated_land_cells` has anything new to fire on.
Flooding cannot create a one-cell-wide passage either: it only ever turns land into sea, so the
count of land neighbours of any sea cell can only fall, and the passage predicate needs it to rise.
"""
function fill_small_islands(h; max_cells)
    filled = copy(h)
    Nx, Ny = size(h)
    visited = falses(Nx, Ny)

    for i = 1:Nx, j = 1:Ny
        (h[i, j] >= 0 && !visited[i, j]) || continue

        # Walked in full even when it is far over the limit: stopping early would leave the rest of
        # the mainland unvisited, and the remainder would then be picked up as small components of
        # its own and flooded.
        cells = land_component!(visited, h, (i, j))
        length(cells) <= max_cells || continue
        any(cell -> cell[1] in (1, Nx) || cell[2] in (1, Ny), cells) && continue

        depths = eltype(h)[]
        for (ci, cj) in cells, (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = ci + di, cj + dj
            (1 <= ii <= Nx && 1 <= jj <= Ny) || continue
            h[ii, jj] < 0 && push!(depths, h[ii, jj])
        end
        # An interior component always has one, since every neighbour of it is in-domain sea; the
        # guard keeps the helper total for any other input.
        isempty(depths) && continue

        island_depth = mean(depths)
        for cell in cells
            filled[cell...] = island_depth
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
