# Rivers: the generic half of the river pipeline. `add_rivers` copies a prepared forcing file
# and writes river relaxation into the surface level of the copy, one grid cell per river. Only
# the three hooks `river_locations`, `river_series` and `download_rivers` are dataset-specific —
# `src/Forcing/of800_rivers.jl` is the template to copy for a new river dataset.

const RIVER_NEIGHBOUR_OFFSETS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))

"""
    RiverLocation(id, name, longitude, latitude)

One river outlet as the source dataset describes it, before it is snapped to a grid cell.
"""
struct RiverLocation
    id::Int
    name::String
    longitude::Float64
    latitude::Float64
end

"""
    RiverCell(location, i, j, distance)

A `RiverLocation` snapped to the grid cell `(i, j)` that will receive its relaxation.
`distance` is how many cells the outlet had to move to reach a coastal water cell.
"""
struct RiverCell
    location::RiverLocation
    i::Int
    j::Int
    distance::Float64
end

"""
    river_locations(config)

The river outlets of `config`, as a `Vector{RiverLocation}`. Required hook for every
`AbstractRiverConfig`.
"""
function river_locations end

"""
    river_series(config, times)

The river values to write, as a `Dict` mapping a FjordSim forcing variable name (`"T"`, `"S"`,
...) to a `(river, time)` matrix with one row per `river_locations` entry and one column per
entry of `times`. Required hook for every `AbstractRiverConfig`.
"""
function river_series end

"""
    download_rivers(config)

Fetch the source data of a river dataset. Only needed by a dataset that downloads.
"""
function download_rivers end

"""
    river_search_radius(config)

How far, in grid cells, `add_rivers` looks for a coastal water cell when an outlet does not
land on one. Defaults to `config.search_radius`.
"""
river_search_radius(config::AbstractRiverConfig) = config.search_radius

"""
    coastal_water_mask(target_grid, edge)

Surface-level water mask of `target_grid`, as a `Ny`-by-`Nx`-shaped `(i, j)` matrix. Built from
the same `water_mask` that `prepare_forcing` uses, so "water" means exactly what it means when
the forcing file is written, including `PartialCellBottom` cells.
"""
function coastal_water_mask(target_grid, edge)
    mask = water_mask(target_grid, Center, Center, edge)
    return mask[:, :, size(mask, 3)]
end

"""
    is_coastal_cell(mask, i, j)

Whether `(i, j)` is a water cell with at least one land cell among its eight neighbours. A
river has to enter the domain at the coastline: an outlet snapped to open water would inject
freshwater in the middle of the fjord.
"""
function is_coastal_cell(mask, i, j)
    checkbounds(Bool, mask, i, j) || return false
    mask[i, j] || return false

    for (di, dj) in RIVER_NEIGHBOUR_OFFSETS
        neighbour_i = i + di
        neighbour_j = j + dj
        checkbounds(Bool, mask, neighbour_i, neighbour_j) || continue
        mask[neighbour_i, neighbour_j] || return true
    end

    return false
end

"""
    nearest_coastal_cell(mask, i, j, radius)

The coastal water cell closest to `(i, j)`, as `(i, j, distance)`, searching rings of growing
radius out to `radius` cells. Returns `nothing` when no ring holds one.

Ties are broken by the iteration order — latitude offset ascending outermost, longitude offset
ascending inside it — which is the order the reference implementation uses, so both pick the
same cell.
"""
function nearest_coastal_cell(mask, i, j, radius)
    is_coastal_cell(mask, i, j) && return (i, j, 0.0)

    for ring = 1:radius
        best = nothing

        for dj = -ring:ring, di = -ring:ring
            distance = sqrt(di^2 + dj^2)
            ring - 1//2 <= distance <= ring + 1//2 || continue
            is_coastal_cell(mask, i + di, j + dj) || continue
            if isnothing(best) || distance < best[3]
                best = (i + di, j + dj, distance)
            end
        end

        isnothing(best) || return best
    end

    return nothing
end

"""
    river_cells(target_grid, locations, edge, radius)

Snap each `RiverLocation` to the coastal water cell that will carry it. An outlet is located by
independent nearest-node lookups in longitude and latitude, then moved to the nearest coastal
water cell within `radius`.

Outlets outside the grid, and outlets with no coastal cell in reach, are dropped with a
warning — writing a river into a land cell would be silently lost when the model reads the
forcing back.
"""
function river_cells(target_grid, locations, edge, radius)
    mask = coastal_water_mask(target_grid, edge)
    longitudes = Array(λnodes(target_grid, Center()))
    latitudes = Array(φnodes(target_grid, Center()))

    cells = RiverCell[]
    for location in locations
        label = "river $(location.id) ($(location.name))"

        if !(first(longitudes) < location.longitude < last(longitudes)) ||
           !(first(latitudes) < location.latitude < last(latitudes))
            @info "Skipping $label: outlet is outside the grid"
            continue
        end

        i = argmin(abs.(longitudes .- location.longitude))
        j = argmin(abs.(latitudes .- location.latitude))

        nearest = nearest_coastal_cell(mask, i, j, radius)
        if isnothing(nearest)
            @warn "Skipping $label: no coastal water cell within $radius cells of ($i, $j)"
            continue
        end

        push!(cells, RiverCell(location, nearest[1], nearest[2], nearest[3]))
    end

    return cells
end

"""
    add_rivers(config::FjordConfig)
    add_rivers(target_grid, config::AbstractForcingConfig)

Write river relaxation on top of the forcing file prepared by `prepare_forcing`, into the copy
at `river_forcing_path(config.rivers)`. Returns `nothing` when the setup has no rivers, or no
forcing at all.

Rivers enter as relaxation, not as a mass flux: each river cell gets its value and a lambda of
`1 / relaxation_timescale` at the surface level, for every time step. That lambda is well
inside the `|λ| < 1` regime of `ForcingFromFile`, so the river cells relax toward the river
values while the rest of the domain keeps whatever `prepare_forcing` wrote — including the
boundary relaxation band, which a river cell landing inside it would override.

The `FjordConfig` method is the setup-level driver: it reads the grid from the processed
bathymetry, so the river cells are snapped against the same land mask the forcing was written
with, and downloads the river data first. The original forcing file is never modified, so the
step can be re-run without redoing `prepare_forcing`.
"""
function add_rivers(config::FjordConfig)
    isnothing(config.forcing_config) && return nothing
    rivers = config.forcing_config.rivers
    isnothing(rivers) && return nothing

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    download_rivers(rivers)

    # The grid stays on the CPU: building the land masks walks `peripheral_node` cell by cell.
    grid = ImmersedBoundaryGrid(bathymetry_file, CPU(), config.grid_config.halo)
    result = add_rivers(grid, config.forcing_config)

    @info "Placed $(length(result.cells)) of $(length(river_locations(rivers))) rivers"
    @info "Patched variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing with rivers saved to $(result.output_file)"

    return result
end

add_rivers(target_grid, config::AbstractForcingConfig) = add_rivers(target_grid, config, config.rivers)

add_rivers(target_grid, config::AbstractForcingConfig, ::Nothing) = nothing

function add_rivers(target_grid, config::AbstractForcingConfig, rivers::AbstractRiverConfig)
    forcing_file = forcing_path(config)
    isfile(forcing_file) || error(
        "Forcing file $forcing_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_forcing` for this setup first.",
    )

    times = NCDataset(forcing_file) do ds
        DateTime.(ds["time"][:])
    end

    locations = river_locations(rivers)
    cells = river_cells(target_grid, locations, config.open_edge, river_search_radius(rivers))
    isempty(cells) && error("None of the $(length(locations)) river outlets landed on the grid.")

    report_river_cells(cells)

    series = river_series(rivers, times)
    output_file = river_forcing_path(rivers)
    abspath(output_file) == abspath(forcing_file) && error(
        "The river output would overwrite the prepared forcing at $forcing_file. " *
        "Give the river config an `output_file` that differs from the forcing config's.",
    )

    @info "Copying $forcing_file to $output_file"
    isdir(dirname(output_file)) || mkpath(dirname(output_file))
    cp(forcing_file, output_file; force = true)

    lambda = Float32(1 / rivers.relaxation_timescale)
    written = write_rivers(output_file, cells, locations, series, lambda, length(times))

    @info "Finished adding rivers"

    return (; output_file, cells, variables = written, times)
end

"""
    report_river_cells(cells)

Log where each river landed, and warn about rivers sharing a cell — the last one written wins.
"""
function report_river_cells(cells)
    for cell in cells
        moved = iszero(cell.distance) ? "" : " (moved $(round(cell.distance, digits = 2)) cells)"
        @info "River $(cell.location.id) ($(cell.location.name)) at cell ($(cell.i), $(cell.j))$moved"
    end

    seen = Dict{Tuple{Int,Int},Int}()
    for cell in cells
        previous = get(seen, (cell.i, cell.j), nothing)
        isnothing(previous) ||
            @warn "River $(cell.location.id) shares cell ($(cell.i), $(cell.j)) with river $previous; " *
                  "the later river wins"
        seen[(cell.i, cell.j)] = cell.location.id
    end

    return cells
end

"""
    write_rivers(output_file, cells, locations, series, lambda, n_times)

Patch the surface level of each river cell in `output_file`, in place. Variables the forcing
file does not carry are skipped with a warning, so a river dataset may offer more than a given
setup prepares.

The file is chunked one horizontal slab per `(level, time)`, so this reads, patches and writes
back whole surface slabs rather than individual cells — a point write would rewrite the entire
slab anyway.
"""
function write_rivers(output_file, cells, locations, series, lambda, n_times)
    rows = Dict(location.id => row for (row, location) in enumerate(locations))
    written = String[]

    NCDataset(output_file, "a") do ds
        surface = ds.dim["Nz"]

        for (name, values) in series
            if !haskey(ds, name) || !haskey(ds, name * "_lambda")
                @warn "Skipping river variable $name: not in $(basename(output_file))"
                continue
            end

            size(values, 2) == n_times || error(
                "River series for $name has $(size(values, 2)) time steps, " *
                "but the forcing file has $n_times.",
            )

            value_variable = ds[name]
            lambda_variable = ds[name*"_lambda"]

            for index = 1:n_times
                value_slab = value_variable[:, :, surface, index]
                lambda_slab = lambda_variable[:, :, surface, index]

                for cell in cells
                    value_slab[cell.i, cell.j] = values[rows[cell.location.id], index]
                    lambda_slab[cell.i, cell.j] = lambda
                end

                value_variable[:, :, surface, index] = value_slab
                lambda_variable[:, :, surface, index] = lambda_slab
            end

            push!(written, name)
        end
    end

    isempty(written) && error("None of the river variables $(sort(collect(keys(series)))) are in $output_file.")

    @info "Wrote rivers into $(join(sort(written), ", ")) at $(length(cells)) cell(s)"

    return sort(written)
end
