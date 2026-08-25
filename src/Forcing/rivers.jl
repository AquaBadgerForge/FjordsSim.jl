# Rivers: the generic half of the river pipeline. `add_rivers` writes river relaxation into the
# surface level of a forcing file, one grid cell per river — either a copy of the one
# `prepare_forcing` wrote, or a river-only file of its own. Only the three hooks `river_locations`,
# `river_series` and `download_rivers` are dataset-specific — `src/Forcing/of800_rivers.jl` is the
# template to copy for a new river dataset.

const RIVER_NEIGHBOUR_OFFSETS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))

# Spacing of a standalone river file's time axis. River records are daily and `river_series` matches
# them by calendar date, so a finer axis would only repeat records.
const RIVER_TIME_SPACING = Day(1)

"""
Global attribute marking a forcing file that carries only rivers, so a reader wanting ocean state
rather than forcing terms can refuse it rather than reading `NaN` everywhere. Set by
`create_river_forcing_file`, read by `FjordSim.Simulations.forcing_state`.
"""
const RIVERS_ONLY_ATTRIBUTE = "rivers_only"

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
    river_minimum_levels(config)

How many wet levels a column must have before `add_rivers` will put an outlet in it. Defaults to
`0`, which accepts any water cell and is what the placement did before this hook existed.

Unlike `river_search_radius` just above, the fallback does **not** read a field of the same name.
That one is essential to placement and every river config has always had to state it; this one is an
opt-in refinement, so a config written before it existed — or one that simply does not care — must
keep working rather than meeting a `FieldError` from inside `add_rivers`. A config that wants the
rule overloads this hook, as `OF800RiversConfig` does.

A river is written into the *surface* level alone, so the column beneath it is what has to carry
the exchange the freshening drives. Give it one cell and it cannot: the fresh surface cell sets up
an estuarine circulation, and the salty inflow at depth concentrates in the single cell below
instead of spreading through a column. On `oslofjorden` the four outlets that landed on the
`minimum_depth` floor — two 1 m cells — held 32 to 64 psu below a surface cell at 0, while all
fifteen outlets with four levels or more stayed between 29 and 35. Column salt was conserved
throughout, so nothing was created; the column simply could not resolve the redistribution.

This is the depth counterpart of the rule `is_coastal_cell` already applies: a river may not enter
open water, and it may not enter a column too shallow to carry it either.
"""
river_minimum_levels(::AbstractRiverConfig) = 0

"""
    coastal_water_mask(target_grid, minimum_levels)

Surface-level water mask of `target_grid`, as an `(i, j)` matrix, with every column of fewer than
`minimum_levels` wet cells masked out. Built from the same `water_mask` that `prepare_forcing`
uses, so "water" means exactly what it means when the forcing file is written, including
`PartialCellBottom` cells.

Masking the shallow columns out of the mask itself, rather than testing depth separately, is what
keeps `is_coastal_cell` and `nearest_coastal_cell` unchanged: a column too shallow to carry a river
simply counts as shore for this purpose, so the nearest acceptable cell is by construction both
coastal and deep enough. `minimum_levels` of `0` or `1` masks nothing.

The open edge does not enter: this asks for the tracer location, and every `open_boundary_water!`
method returns the mask untouched unless the location is staggered across its own edge. A river
enters at a tracer cell, so the edge could only ever have been a no-op here.
"""
function coastal_water_mask(target_grid, minimum_levels)
    mask = water_mask(target_grid, Center, Center, nothing)
    surface = mask[:, :, size(mask, 3)]
    levels = dropdims(sum(mask; dims = 3); dims = 3)
    return surface .& (levels .>= minimum_levels)
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
    river_cells(target_grid, locations, radius, minimum_levels)

Snap each `RiverLocation` to the coastal water cell that will carry it. An outlet is located by
independent nearest-node lookups in longitude and latitude, then moved to the nearest coastal
water cell of at least `minimum_levels` wet levels within `radius`.

Outlets outside the grid, and outlets with no such cell in reach, are dropped with a
warning — writing a river into a land cell would be silently lost when the model reads the
forcing back.
"""
function river_cells(target_grid, locations, radius, minimum_levels)
    mask = coastal_water_mask(target_grid, minimum_levels)
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
            @warn "Skipping $label: no coastal water cell of at least $minimum_levels levels " *
                  "within $radius cells of ($i, $j)"
            continue
        end

        push!(cells, RiverCell(location, nearest[1], nearest[2], nearest[3]))
    end

    return cells
end

"""
    add_rivers(config::FjordConfig)
    add_rivers(target_grid, config::AbstractForcingConfig; coverage = nothing)

Write river relaxation into the forcing file at `river_forcing_path(config.rivers)`. Returns
`nothing` when the setup has no rivers, or no forcing config at all.

Rivers enter as relaxation, not as a mass flux: each river cell gets its value and a lambda of
`1 / relaxation_timescale` at the surface level, for every time step. That lambda is well
inside the `|λ| < 1` regime of `ForcingFromFile`, so the river cells relax toward the river
values while the rest of the domain keeps whatever the rest of the file says.

Where that file comes from is the river config's `standalone`:

- `false` copies the file `prepare_forcing` wrote and patches the copy, taking its time axis from
  it. The original is never modified, so the step re-runs without redoing `prepare_forcing`, and a
  missing prepared forcing is an error naming the step that produces it.
- `true` writes the file from scratch, carrying *only* rivers — no interior forcing to download or
  regrid. Its time axis comes from the run window instead of from a prepared file, so it needs
  `coverage`. See `create_river_forcing_file`.

The `FjordConfig` method is the setup-level driver: it reads the grid from the processed
bathymetry, so the river cells are snapped against the same land mask the forcing would be written
with, downloads the river data first, and takes `coverage` from the simulation config.
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
    grid = simulation_grid(config.grid_config, bathymetry_file, CPU())
    result = add_rivers(grid, config.forcing_config; coverage = coverage_window(config.simulation_config))

    @info "Placed $(length(result.cells)) of $(length(river_locations(rivers))) rivers"
    @info "Patched variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing with rivers saved to $(result.output_file)"

    return result
end

add_rivers(target_grid, config::AbstractForcingConfig; coverage = nothing) =
    add_rivers(target_grid, config, config.rivers; coverage)

add_rivers(target_grid, config::AbstractForcingConfig, ::Nothing; coverage = nothing) = nothing

function add_rivers(
    target_grid,
    config::AbstractForcingConfig,
    rivers::AbstractRiverConfig;
    coverage = nothing,
)
    forcing_file = forcing_path(config)
    times = river_forcing_times(rivers, forcing_file, coverage)

    locations = river_locations(rivers)
    cells = river_cells(
        target_grid, locations, river_search_radius(rivers), river_minimum_levels(rivers),
    )
    isempty(cells) && error("None of the $(length(locations)) river outlets landed on the grid.")

    report_river_cells(cells)

    series = river_series(rivers, times)
    output_file = river_forcing_path(rivers)
    abspath(output_file) == abspath(forcing_file) && error(
        "The river output would overwrite the prepared forcing at $forcing_file. " *
        "Give the river config an `output_file` that differs from the forcing config's.",
    )

    isdir(dirname(output_file)) || mkpath(dirname(output_file))
    if rivers.standalone
        create_river_forcing_file(output_file, target_grid, sort(collect(keys(series))), times)
    else
        @info "Copying $forcing_file to $output_file"
        cp(forcing_file, output_file; force = true)
    end

    lambda = Float32(1 / rivers.relaxation_timescale)
    written = write_rivers(output_file, cells, locations, series, lambda, length(times))

    @info "Finished adding rivers"

    return (; output_file, cells, variables = written, times)
end

"""
    river_forcing_times(rivers, forcing_file, coverage)

The time axis the river forcing file is written on: the prepared forcing's own axis when patching a
copy of it, and a daily axis spanning the run window when `standalone`.

The standalone axis is anchored at the window's first instant and extended until it passes the last,
so `forcing_date_range` reports a span `validate_time_coverage` accepts — every reader indexes time
`Cyclical()`, which wraps rather than failing outside its data, and the coverage check is what keeps
that wrap unreachable.

Daily because a river dataset's records are daily and `river_series` matches them by *calendar
date*, so a finer cadence would ask for the same record several times over. A sub-daily dataset would
want a hook here.
"""
function river_forcing_times(rivers::AbstractRiverConfig, forcing_file, coverage)
    rivers.standalone || return prepared_forcing_times(forcing_file)

    isnothing(coverage) && error(
        "A standalone river forcing file takes its time axis from the run window, but this setup " *
        "names no `simulation_config` to read `start_date` and `stop_time` from. Name one, or set " *
        "the river config's `standalone` to false and prepare interior forcing first.",
    )

    first_date, last_date = coverage
    dates = collect(first_date:RIVER_TIME_SPACING:last_date)
    # Extend until the axis reaches past the window's end, and to at least two records: the reader
    # keeps two time indices in memory, so a single-record axis reads past the end of the file.
    while last(dates) < last_date || length(dates) < 2
        push!(dates, last(dates) + RIVER_TIME_SPACING)
    end

    return dates
end

function prepared_forcing_times(forcing_file)
    isfile(forcing_file) || error(
        "Forcing file $forcing_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_forcing` for this setup first, or set the river " *
        "config's `standalone` to true for a forcing file carrying only rivers.",
    )

    return NCDataset(forcing_file) do ds
        DateTime.(ds["time"][:])
    end
end

"""
    create_river_forcing_file(output_file, target_grid, names, dates)

Write an otherwise-empty forcing file for `write_rivers` to patch: the layout `forcing_from_file`
reads, on `dates`, with a `value`/`_lambda` pair per name in `names` and no data in it.

The dimensions and coordinates are `define_forcing_dimensions!`'s, so this is the same contract
`prepare_forcing` writes and the reader's dimension check against the grid holds unchanged.

Only the **surface level** is filled in, with `NaN32` values and zero lambdas — that is the level
`write_rivers` patches, and the levels below it carry no river forcing at all, so they are left
unwritten and read back as the file's `_FillValue`. Both are inert twice over: a non-finite value
becomes the `-999.0` sentinel every branch of `ForcingFromFile` gates on with `value > -990`, and a
non-finite lambda fails each of `λ > 1`, `λ < -1` and `-1 < λ < 1`. It is also what dry cells already
read as in a file `prepare_forcing` wrote.

The `rivers_only` global attribute marks what the file is, so a reader that needs actual ocean state
— `FromForcing` initial conditions, which would otherwise start the run from zeros everywhere — can
refuse it instead of silently succeeding.
"""
function create_river_forcing_file(output_file, target_grid, names, dates)
    @info "Creating river-only forcing file $output_file with $(length(dates)) time steps: " *
          "$(first(dates)) to $(last(dates))"
    isfile(output_file) && rm(output_file; force = true)

    NCDataset(output_file, "c") do ds
        define_forcing_dimensions!(ds, target_grid, dates)
        ds.attrib[RIVERS_ONLY_ATTRIBUTE] = "true"

        for name in names
            dimensions = forcing_dimension_names(name)
            nx, ny = ds.dim[dimensions[1]], ds.dim[dimensions[2]]
            surface = ds.dim["Nz"]

            for (variable_name, fill_value) in ((name, NaN32), (name * "_lambda", 0.0f0))
                variable = defVar(ds, variable_name, Float32, dimensions;
                    chunksizes = [nx, ny, 1, 1], deflatelevel = FORCING_DEFLATE_LEVEL,
                    attrib = ["_FillValue" => NaN32])
                slab = fill(fill_value, nx, ny)
                for index in eachindex(dates)
                    variable[:, :, surface, index] = slab
                end
            end
        end
    end

    return output_file
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
