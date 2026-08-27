# OF800 rivers adapter: the source-specific half of the river pipeline. Copy this file's shape
# to add a new river dataset — it owns a config subtype and the hook methods `river_locations`
# and `river_series` (plus `download_rivers` if it fetches data). Everything else comes from the
# generic core in rivers.jl.

using Downloads

# Per-file Dropbox links for the two source files, each shared as "anyone with the link".
# Two things they must keep: `dl=1`, so the link serves the file rather than the web app, and
# the per-*file* `/scl/fi/` form — a `/scl/fo/` folder link cannot be used, because Dropbox
# renders folder listings client-side and serves the whole folder as one archive. The `st`
# parameter a browser adds when copying a link is not needed and is left off, so these do not
# go stale. A link that loses public access answers HTTP 200 with a login page instead of
# failing, which is what `validate_river_download` catches.
const OF800_LOCATIONS_URL = "https://www.dropbox.com/scl/fi/bzjwsnm97laa5w1xwa83u/OF800_rivers.csv?rlkey=0z3l30vlbvjippfy0lrymybqg&dl=1"
const OF800_SERIES_URL = "https://www.dropbox.com/scl/fi/3bfcu0wjje2uhdzbcso91/of800_rivers_v9_1990_2022_RA1.nc?rlkey=dxqgwj325ogx5odir5wnurquo&dl=1"

# The `river` coordinate of the source file runs 5…29 while the outlet CSV numbers the same
# rivers 1…25.
const OF800_RIVER_NUMBER_OFFSET = 4

"""
    OF800RiversConfig(; data_root, kwargs...)

Rivers of the OF800 Oslofjord setup: outlet coordinates from a CSV and daily time series from a
ROMS river forcing NetCDF.

# Fields
- `data_root`: directory holding the two source files and receiving the output.
- `output_file`: the rivers-augmented copy of the forcing file, relative to `data_root`.
- `locations_file`, `series_file`: source file names, relative to `data_root`.
- `locations_url`, `series_url`: per-file download URLs, defaulting to the dataset's own Dropbox
  links. See the note at the top of this file for the form they have to take.
- `variables`: source variable name => FjordSim forcing variable name.
- `constants`: FjordSim forcing variable name => a value held constant in time, for a variable
  the source does not carry. River salinity is 0 by definition of fresh water.
- `relaxation_timescale`: seconds; the river cells relax this fast toward the river values.
- `search_radius`: how far to look for a coastal water cell, in grid cells.
- `minimum_levels`: how many wet levels a column must have before an outlet may be placed in it.
  `0` accepts any water cell, which is what the placement did before.

  A river is relaxed into the *surface* level alone, so the column beneath it is what carries the
  exchange the freshening drives, and one cell cannot. On `oslofjorden` the four outlets that landed
  on the `minimum_depth` floor — a column of two 1 m cells — ran to 32 and 64 psu in the cell below a
  surface cell held at 0, while all fifteen outlets with four levels or more stayed between 29 and
  35 psu. Column salt stayed conserved, so this is a redistribution the column cannot resolve rather
  than anything being created. See `river_minimum_levels`.
- `standalone`: whether `add_rivers` writes a forcing file carrying only rivers rather than patching
  a copy of the prepared interior forcing. `false` by default, so the step keeps needing
  `prepare_forcing` unless a setup says otherwise.
"""
Base.@kwdef mutable struct OF800RiversConfig <: AbstractRiverConfig
    data_root::String
    output_file::String = "forcing_rivers.nc"
    plot_file::String = "forcing_rivers.png"
    locations_file::String = "OF800_rivers.csv"
    series_file::String = "of800_rivers_v9_1990_2022_RA1.nc"
    locations_url::String = OF800_LOCATIONS_URL
    series_url::String = OF800_SERIES_URL
    variables::Dict{String,String} = Dict("river_temp" => "T")
    constants::Dict{String,Float64} = Dict("S" => 0.0)
    relaxation_timescale::Float64 = 3600.0
    search_radius::Int = 10
    minimum_levels::Int = 0
    standalone::Bool = false
end

"""
    river_minimum_levels(config::OF800RiversConfig)

The wet-level floor for outlet placement, from `config.minimum_levels`. See the generic
`river_minimum_levels` for why the supertype fallback is a plain `0` rather than this field read.
"""
river_minimum_levels(config::OF800RiversConfig) = config.minimum_levels

"""
    river_locations_path(config)
    river_series_path(config)

Resolve the two source files against `config.data_root`. An absolute name is returned
unchanged, so a single copy of the river data can be shared across setups.
"""
river_locations_path(config::OF800RiversConfig) = joinpath(config.data_root, config.locations_file)
river_series_path(config::OF800RiversConfig) = joinpath(config.data_root, config.series_file)

"""
    download_rivers(target_grid, config::OF800RiversConfig)

Fetch the two source files if they are not present. A file that is already there is left alone,
so the step is cheap to re-run — the series file is ~176 MB.

`target_grid` goes unused: this dataset's outlets are a column of the locations CSV, so it is told
where its rivers are rather than having to find them in the domain.
"""
function download_rivers(target_grid, config::OF800RiversConfig)
    ensure_river_file(river_locations_path(config), config.locations_url)
    ensure_river_file(river_series_path(config), config.series_url)
    return config.data_root
end

function ensure_river_file(filepath, url)
    isfile(filepath) && return filepath

    mkpath(dirname(filepath))
    @info "Downloading $(basename(filepath)) from $url"
    Downloads.download(url, filepath)

    return validate_river_download(filepath, url)
end

"""
    validate_river_download(filepath, url)

Reject a download that is a web page rather than the file, and delete it. A Dropbox link that is
not publicly shared answers with a login page and HTTP 200, so the download reports success and
writes that page over the data file — which would otherwise surface much later as a confusing
parse error.
"""
function validate_river_download(filepath, url)
    first_byte = open(io -> read(io, 1), filepath)

    if isempty(first_byte) || first_byte[1] == UInt8('<')
        rm(filepath; force = true)
        error(
            "Downloading $(basename(filepath)) returned a web page instead of the file. " *
            "Check that $url is still shared with anyone who has the link, and that it is a " *
            "per-file link ending in `dl=1`.",
        )
    end

    return filepath
end

"""
    river_locations(config::OF800RiversConfig)

Read the outlet CSV: `River number,Name,LatOutlet,LonOutlet,Zero`. Parsed by hand because the
file is a couple of dozen rows and the project carries no CSV reader.
"""
function river_locations(config::OF800RiversConfig)
    filepath = river_locations_path(config)
    isfile(filepath) || error("River locations file $filepath does not exist. Run download_rivers first.")

    locations = RiverLocation[]
    for (number, line) in enumerate(eachline(filepath))
        number == 1 && continue  # header
        isempty(strip(line)) && continue

        fields = split(line, ',')
        length(fields) >= 4 || error("Malformed river locations line $number in $filepath: $line")

        push!(
            locations,
            RiverLocation(
                parse(Int, fields[1]),
                String(strip(fields[2])),
                parse(Float64, fields[4]),
                parse(Float64, fields[3]),
            ),
        )
    end

    isempty(locations) && error("No rivers found in $filepath.")

    return locations
end

"""
    river_series(config::OF800RiversConfig, times)

Read the daily river series for `times`. The source is a ROMS river file whose values are
identical across its `s_rho` levels, so only the topmost level is read, and whose records are
stamped at midday, so a forcing date matches the record whose date it shares.
"""
function river_series(config::OF800RiversConfig, times)
    filepath = river_series_path(config)
    isfile(filepath) || error("River series file $filepath does not exist. Run download_rivers first.")

    locations = river_locations(config)
    series = Dict{String,Matrix{Float32}}()

    NCDataset(filepath) do ds
        indices = river_time_indices(ds, times, filepath)
        columns = river_columns(ds, locations, filepath)
        levels = ds.dim["s_rho"]

        for (source_name, name) in config.variables
            haskey(ds, source_name) || error("River variable $source_name is not in $filepath.")
            values = Float32.(ds[source_name][:, levels, :])  # (river, s_rho, river_time)
            series[name] = values[columns, indices]
        end
    end

    for (name, value) in config.constants
        series[name] = fill(Float32(value), length(locations), length(times))
    end

    return series
end

"""
    river_time_indices(ds, times, filepath)

Position of each forcing date in the river file's time axis, matched by calendar date so the
river file's midday stamps line up with whatever time of day the forcing carries.
"""
function river_time_indices(ds, times, filepath)
    river_times = DateTime.(ds["river_time"][:])
    positions = Dict(Dates.Date(date) => index for (index, date) in enumerate(river_times))

    indices = Vector{Int}(undef, length(times))
    for (index, date) in enumerate(times)
        position = get(positions, Dates.Date(date), nothing)
        isnothing(position) && error(
            "Forcing date $date has no river record in $filepath, " *
            "which covers $(first(river_times)) to $(last(river_times)).",
        )
        indices[index] = position
    end

    return indices
end

"""
    river_columns(ds, locations, filepath)

Position of each outlet in the river file's `river` dimension. The file numbers its rivers from
`OF800_RIVER_NUMBER_OFFSET + 1` while the outlet CSV numbers them from 1.
"""
function river_columns(ds, locations, filepath)
    numbers = round.(Int, ds["river"][:]) .- OF800_RIVER_NUMBER_OFFSET
    positions = Dict(number => index for (index, number) in enumerate(numbers))

    columns = Vector{Int}(undef, length(locations))
    for (index, location) in enumerate(locations)
        position = get(positions, location.id, nothing)
        isnothing(position) && error("River $(location.id) ($(location.name)) is not in $filepath.")
        columns[index] = position
    end

    return columns
end
