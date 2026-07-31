# NORA3 adapter: the source-specific half of the atmosphere pipeline. Copy this file's shape to
# add a new atmosphere dataset — it owns a config subtype, its own download, and the three hook
# methods `atmosphere_time_steps`, `atmosphere_source_grid` and `atmosphere_variable_names`
# (plus `download_atmosphere` if it fetches data). Everything else comes from the generic core in
# Atmospheres.jl.

using Dates: Date, Day, lastdayofmonth
import Dates

const NORA3_OPENDAP_URL = "https://thredds.met.no/thredds/dodsC/nora3"

# NORA3 publishes one file per forecast lead hour of a 6-hourly run, each holding a single time
# step. Leads 4..9 of one run supply six consecutive hours, so the four daily runs tile a day
# exactly with no overlap and no dedup pass. Lead 3 is read only to de-accumulate lead 4.
const NORA3_RUN_HOURS = (0, 6, 12, 18)
const NORA3_FIRST_LEAD = 3
const NORA3_LAST_LEAD = 9
const NORA3_SECONDS_PER_LEAD = 3600

# Extra source cells kept around the prepared grid's bounding box, so bilinear interpolation at
# the outermost prepared cell still has a full stencil.
const NORA3_SUBSET_MARGIN_CELLS = 2

const NORA3_OPEN_RETRIES = 4
const NORA3_RETRY_BACKOFF_SECONDS = 5

# Prepared name => source variable, for the fields taken straight from the file.
#
# NORA3 publishes air temperature in Kelvin and the prepared file wants Kelvin, so nothing is
# converted. The `/ NORA3_SECONDS_PER_LEAD` in de-accumulation is the pipeline's only unit change.
const NORA3_INSTANTANEOUS_VARIABLES = Dict(
    "air_temperature_2m" => "air_temperature_2m",
    "specific_humidity_2m" => "specific_humidity_2m",
    "air_pressure_at_sea_level" => "air_pressure_at_sea_level",
)

# Prepared name => source variable, for the fluxes accumulated since the start of each forecast
# run.
#
# Shortwave is the *downwelling* flux, where the reference Python pipeline this ports uses
# `integral_of_surface_net_downward_shortwave_flux_wrt_time`: `NORA3PrescribedRadiation` treats
# `swrad` as downwelling and applies its own surface albedo, so a net flux counts the albedo
# twice. Longwave was already the downwelling flux there.
const NORA3_ACCUMULATED_VARIABLES = Dict(
    "swrad" => "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time",
    "lwrad" => "integral_of_surface_downwelling_longwave_flux_in_air_wrt_time",
    "precipitation" => "precipitation_amount_acc",
)

# Wind components relative to the Lambert grid axes, rotated to eastward/northward on the source
# grid by `rotate_to_east_north` before anything is interpolated.
const NORA3_EASTWARD_SOURCE = "x_wind_10m"
const NORA3_NORTHWARD_SOURCE = "y_wind_10m"
const NORA3_EASTWARD_NAME = "u_wind_10m"
const NORA3_NORTHWARD_NAME = "v_wind_10m"

"""
    NORA3Config(; data_root, output_directory, years, kw...)

Atmosphere configuration for the MET Norway NORA3 reanalysis, served per forecast lead hour over
OPeNDAP from `thredds.met.no`.

Only `opendap_url` is defaulted among the source's own settings, being NORA3's public endpoint.
`data_root`, `output_directory` and `years` have no defaults, so a setup that forgets one gets an
`UndefKeywordError` rather than a silent fallback.

The download normalizes NORA3's forecast structure away: it selects the forecast lead supplying
each hour, de-accumulates the accumulated fluxes *within* a run, rotates the winds to
eastward/northward, and writes one gap-free hourly file per month already carrying the eight
prepared variable names. `prepare_atmosphere` is then a pure regrid.

Two deliberate differences from the reference Python pipeline this ports, and therefore from a
file it produced: `swrad` here is the downwelling shortwave flux rather than the net one, and the
hourly axis is built from the forecast leads themselves rather than padded at the start by
repeating the first step.

# Fields
- `data_root`: Directory the other paths resolve against.
- `output_directory`: Where the downloaded monthly files go, via `atmosphere_directory`.
- `output_file`: The prepared file, via `atmosphere_path`.
- `plot_file`: The diagnostic plot, via `plot_path`.
- `opendap_url`: Base OPeNDAP URL of the NORA3 file archive.
- `resolution`: Prepared grid spacing in degrees.
- `padding`: Degrees of margin added around the ocean domain.
- `years`: Calendar years to prepare. Hours 00:00 to 03:00 of 1 January come from the 18Z run of
  31 December of the preceding year, which is fetched automatically.

A relative `output_file`, `plot_file` or `output_directory` resolves against `data_root`; an
absolute one relocates just that file.
"""
Base.@kwdef mutable struct NORA3Config <: AbstractAtmosphereConfig
    data_root::String
    output_directory::String
    output_file::String = "atmosphere.nc"
    plot_file::String = "atmosphere.png"
    opendap_url::String = NORA3_OPENDAP_URL
    resolution::Float64 = 0.02
    padding::Float64 = 0.1
    years::Vector{Int}
end

"""
    nora3_monthly_filename(config, year, month)

Name of the downloaded file holding one month of hourly NORA3 data.
"""
nora3_monthly_filename(config::NORA3Config, year, month) = "NORA3_$(year)$(lpad(month, 2, '0')).nc"

nora3_monthly_path(config::NORA3Config, year, month) =
    joinpath(atmosphere_directory(config), nora3_monthly_filename(config, year, month))

"""
    nora3_url(config, run_time, lead)

OPeNDAP URL of one NORA3 forecast lead. The archive's layout is fully determined by the run time
and the lead, so no catalog listing is needed.
"""
function nora3_url(config::NORA3Config, run_time, lead)
    directory = Dates.format(run_time, "yyyy/mm/dd/HH")
    stamp = Dates.format(run_time, "yyyymmddHH")

    return string(config.opendap_url, "/", directory, "/fc", stamp, "_", lpad(lead, 3, '0'), "_fp.nc")
end

"""
    atmosphere_variable_names(config::NORA3Config)

The prepared variables, as downloaded name => prepared name.

This is the identity for NORA3 because the download already writes the prepared names: five of the
eight variables are *derived* rather than copied (the two rotated winds and the three
de-accumulated fluxes), so there is no source variable left to rename at this stage. A dataset
whose download carries the source's own names through would return a real mapping here.
"""
atmosphere_variable_names(config::NORA3Config) =
    Dict(variable.name => variable.name for variable in ATMOSPHERE_VARIABLES)

"""
    atmosphere_time_steps(config::NORA3Config)

The hourly records present in the downloaded monthly files, sorted by date.
"""
function atmosphere_time_steps(config::NORA3Config)
    output_directory = atmosphere_directory(config)
    isdir(output_directory) || error(
        "NORA3 download directory $output_directory does not exist. " *
        "Run scripts/atmosphere_download.jl for this config first.",
    )

    records = AtmosphereRecord[]
    for year in config.years, month = 1:12
        filepath = nora3_monthly_path(config, year, month)
        isfile(filepath) || continue

        NCDataset(filepath) do ds
            for (index, date) in enumerate(ds["time"][:])
                push!(records, AtmosphereRecord(DateTime(date), filepath, index))
            end
        end
    end

    isempty(records) && error(
        "No downloaded NORA3 files found in $output_directory for years " *
        "$(join(config.years, ", ")). Run scripts/atmosphere_download.jl for this config first.",
    )
    sort!(records, by = record -> record.date)

    return unique(record -> record.date, records)
end

"""
    atmosphere_source_grid(config::NORA3Config, filepath)

Geometry of a downloaded monthly file: its Lambert conformal conic subset, as a
`ProjectedAtmosphereGrid`.
"""
function atmosphere_source_grid(config::NORA3Config, filepath)
    return NCDataset(filepath) do ds
        x = Float64.(ds["x"][:])
        y = Float64.(ds["y"][:])
        proj4 = variable(ds, "projection_lambert").attrib["proj4"]

        validate_regular_axis(x, "x")
        validate_regular_axis(y, "y")

        return ProjectedAtmosphereGrid(x, y, proj4)
    end
end

"""
    validate_regular_axis(values, name)

Check a projected axis is long enough and regularly spaced, which `interpolate_to_target!`
assumes when it turns a coordinate into a fractional index.
"""
function validate_regular_axis(values, name)
    length(values) >= 2 || error("The NORA3 $name axis needs at least two points, got $(length(values))")
    spacing = values[2] - values[1]
    all(isapprox(spacing), diff(values)) ||
        error("The NORA3 $name axis is not regularly spaced")

    return nothing
end

# --- Download ---

"""
    download_atmosphere(target_grid, config::NORA3Config)

Download NORA3 for `config.years` into `atmosphere_directory(config)`, one gap-free hourly file
per month.

A month whose file already exists is skipped, so an interrupted download resumes. Each month is
written to a `.tmp` path and renamed on success — a month is close to 900 OPeNDAP opens, so a
half-written file left behind by an interrupted run must never be mistaken for a complete one.

Returns the output directory.
"""
function download_atmosphere(target_grid, config::NORA3Config)
    output_directory = atmosphere_directory(config)
    mkpath(output_directory)

    longitude, latitude = atmosphere_target_axes(target_grid, config)
    pending = [
        (year, month) for year in config.years, month = 1:12 if
        !isfile(nora3_monthly_path(config, year, month))
    ]

    if isempty(pending)
        @info "All requested NORA3 months are already downloaded in $output_directory"
        return output_directory
    end

    @info "Downloading NORA3 years $(join(config.years, ", ")) to $output_directory"
    @info "Subsetting to longitude $(first(longitude)) to $(last(longitude)), " *
          "latitude $(first(latitude)) to $(last(latitude))"

    reference_year, reference_month = first(pending)
    reference_url = nora3_url(config, DateTime(reference_year, reference_month, 1), NORA3_LAST_LEAD)
    subset = nora3_subset(longitude, latitude, reference_url)
    @info "NORA3 subset: $(length(subset.x)) x $(length(subset.y)) source cells"

    for (year, month) in pending
        process_month(year, month, subset, config)
    end

    @info "Finished downloading NORA3 atmosphere"

    return output_directory
end

"""
    NORA3Subset

The NORA3 index window covering the prepared grid, plus everything about the source grid that is
constant across time steps.

`angle` is precomputed once because the rotation depends only on the grid geometry. It is derived
from the subset rather than the full grid, so its last column repeats the second-to-last
difference; that column sits in the `NORA3_SUBSET_MARGIN_CELLS` margin, outside the prepared
grid.

# Fields
- `x_range`, `y_range`: Index ranges into the native grid.
- `x`, `y`: Projected coordinate centers of the subset, in meters.
- `longitude`, `latitude`: Geographic coordinates of the subset.
- `angle`: `grid_rotation_angle` of the subset.
- `proj4`: PROJ.4 definition of the NORA3 projection.
"""
struct NORA3Subset
    x_range::UnitRange{Int}
    y_range::UnitRange{Int}
    x::Vector{Float64}
    y::Vector{Float64}
    longitude::Matrix{Float64}
    latitude::Matrix{Float64}
    angle::Matrix{Float64}
    proj4::String
end

"""
    nora3_subset(longitude, latitude, reference_url)

The NORA3 index window covering the prepared grid's longitude/latitude box, grown by
`NORA3_SUBSET_MARGIN_CELLS`.

The box is the *prepared* grid's, not the ocean grid's: the prepared grid is deliberately padded
beyond the ocean domain, and a subset cut to the ocean domain alone would leave its outer cells
with nothing to interpolate from.
"""
function nora3_subset(longitude, latitude, reference_url)
    return retry_open(reference_url) do ds
        source_longitude = Float64.(Array(ds["longitude"][:, :]))
        source_latitude = Float64.(Array(ds["latitude"][:, :]))

        inside =
            (source_longitude .>= first(longitude)) .&
            (source_longitude .<= last(longitude)) .&
            (source_latitude .>= first(latitude)) .&
            (source_latitude .<= last(latitude))
        any(inside) || error(
            "The prepared grid's box (longitude $(first(longitude)) to $(last(longitude)), " *
            "latitude $(first(latitude)) to $(last(latitude))) does not overlap the NORA3 domain",
        )

        x_range = grown_range(bounding_range(inside, 1), size(inside, 1))
        y_range = grown_range(bounding_range(inside, 2), size(inside, 2))

        x = Float64.(ds["x"][x_range])
        y = Float64.(ds["y"][y_range])
        validate_regular_axis(x, "x")
        validate_regular_axis(y, "y")

        subset_longitude = source_longitude[x_range, y_range]
        subset_latitude = source_latitude[x_range, y_range]
        proj4 = variable(ds, "projection_lambert").attrib["proj4"]

        return NORA3Subset(
            x_range,
            y_range,
            x,
            y,
            subset_longitude,
            subset_latitude,
            grid_rotation_angle(subset_longitude, subset_latitude),
            proj4,
        )
    end
end

"""
    bounding_range(mask, dimension)

The index range along `dimension` holding every true entry of a 2D `mask`.
"""
function bounding_range(mask, dimension)
    reduced = vec(any(mask, dims = dimension == 1 ? 2 : 1))
    indices = findall(reduced)

    return first(indices):last(indices)
end

"""
    grown_range(range, limit)

`range` widened by `NORA3_SUBSET_MARGIN_CELLS` on both sides, clamped to `1:limit`.
"""
grown_range(range, limit) =
    max(1, first(range) - NORA3_SUBSET_MARGIN_CELLS):min(limit, last(range) + NORA3_SUBSET_MARGIN_CELLS)

"""
    retry_open(body, url)

Open an OPeNDAP dataset, run `body` on it and close it, retrying `NORA3_OPEN_RETRIES` times with
exponential backoff. The THREDDS server intermittently refuses connections under load, and a
year's download is thousands of opens, so a single transient failure must not lose a month.
"""
function retry_open(body, url)
    for attempt = 1:NORA3_OPEN_RETRIES
        try
            ds = NCDataset(url)
            try
                return body(ds)
            finally
                close(ds)
            end
        catch exception
            attempt == NORA3_OPEN_RETRIES && rethrow()
            wait_seconds = NORA3_RETRY_BACKOFF_SECONDS * 2^(attempt - 1)
            @warn "Reading $url failed, retrying in $(wait_seconds)s" exception
            sleep(wait_seconds)
        end
    end

    error("Unreachable: retry_open exhausted its attempts without returning or rethrowing")
end

"""
    nora3_runs(year, month)

The forecast runs whose leads 4..9 cover every hour of `year`-`month`.

Each run covers six hours, so a day needs its four runs plus, for hours 00:00 to 03:00, the 18Z
run of the previous day. The 18Z run of the last day of the month is still read: it supplies
22:00 and 23:00.
"""
function nora3_runs(year, month)
    first_day = Date(year, month, 1)
    runs = [DateTime(first_day - Day(1)) + Hour(18)]

    for day in first_day:Day(1):lastdayofmonth(first_day), hour in NORA3_RUN_HOURS
        push!(runs, DateTime(day) + Hour(hour))
    end

    return runs
end

"""
    nora3_month_dates(year, month)

The hourly time axis of one month, gap-free and sorted: every valid time that `nora3_runs`
supplies and that falls inside the month.
"""
function nora3_month_dates(year, month)
    dates = DateTime[]

    for run_time in nora3_runs(year, month), lead = (NORA3_FIRST_LEAD+1):NORA3_LAST_LEAD
        date = run_time + Hour(lead)
        Dates.year(date) == year && Dates.month(date) == month && push!(dates, date)
    end

    return sort!(dates)
end

"""
    process_month(year, month, subset::NORA3Subset, config::NORA3Config)

Download one month into `nora3_monthly_path(config, year, month)`.
"""
function process_month(year, month, subset::NORA3Subset, config::NORA3Config)
    output_path = nora3_monthly_path(config, year, month)
    temporary_path = output_path * ".tmp"
    dates = nora3_month_dates(year, month)
    time_index = Dict(date => index for (index, date) in enumerate(dates))

    @info "Downloading NORA3 $(year)-$(lpad(month, 2, '0')) ($(length(dates)) hourly steps)"

    isfile(temporary_path) && rm(temporary_path)
    output = NCDataset(temporary_path, "c")
    try
        define_monthly_file!(output, subset, dates)

        for run_time in nora3_runs(year, month)
            process_run(output, run_time, subset, config, time_index)
        end
    catch
        close(output)
        isfile(temporary_path) && rm(temporary_path)
        rethrow()
    end
    close(output)

    mv(temporary_path, output_path; force = true)
    @info "Saved $output_path"

    return output_path
end

"""
    define_monthly_file!(ds, subset::NORA3Subset, dates)

Define the downloaded file: the subset's projected and geographic coordinates, its projection,
the hourly time axis, and the eight prepared variables.
"""
function define_monthly_file!(ds, subset::NORA3Subset, dates)
    Nx = length(subset.x)
    Ny = length(subset.y)

    defDim(ds, "x", Nx)
    defDim(ds, "y", Ny)
    defDim(ds, "time", length(dates))

    defVar(ds, "x", Float64, ("x",); attrib = ["units" => "m"])[:] = subset.x
    defVar(ds, "y", Float64, ("y",); attrib = ["units" => "m"])[:] = subset.y
    defVar(ds, "longitude", Float64, ("x", "y"); attrib = ["units" => "degrees_east"])[:, :] =
        subset.longitude
    defVar(ds, "latitude", Float64, ("x", "y"); attrib = ["units" => "degrees_north"])[:, :] =
        subset.latitude
    defVar(ds, "time", dates, ("time",))
    defVar(
        ds,
        "projection_lambert",
        Int32,
        ();
        attrib = ["grid_mapping_name" => "lambert_conformal_conic", "proj4" => subset.proj4],
    )

    for variable in ATMOSPHERE_VARIABLES
        defVar(
            ds,
            variable.name,
            Float32,
            ("x", "y", "time");
            chunksizes = [Nx, Ny, 1],
            deflatelevel = ATMOSPHERE_DEFLATE_LEVEL,
            attrib = ["_FillValue" => NaN32; atmosphere_variable_attributes(variable.name)],
        )
    end

    return ds
end

"""
    process_run(output, run_time, subset::NORA3Subset, config::NORA3Config, time_index)

Read one forecast run's leads and write the hours of it that belong to the month being built.

Leads 3..9 are walked in order, carrying the previous lead's accumulated fluxes forward, so every
de-accumulation is a difference *within* one run. That matters because the accumulators restart at
every run: downwelling longwave for valid time 09:00 reads 11_181_925 J/m² as lead 9 of the 00Z run
but 3_734_820 J/m² as lead 3 of the 06Z run. Differencing a naive concatenation of the runs would
therefore report about -2000 W/m² at each 6-hour boundary.

Each increment is labelled at the *end* of its interval, which is what puts the fluxes on the same
hourly axis as the instantaneous fields. The reference Python pipeline labels them at the interval
midpoint instead and so needs a second `time_acc` axis, which the simulation-time reader cannot
consume — it reads one `time`.
"""
function process_run(output, run_time, subset::NORA3Subset, config::NORA3Config, time_index)
    previous_accumulated = read_lead(nora3_url(config, run_time, NORA3_FIRST_LEAD), subset; instantaneous = false)

    for lead = (NORA3_FIRST_LEAD+1):NORA3_LAST_LEAD
        date = run_time + Hour(lead)
        index = get(time_index, date, 0)
        slabs = read_lead(nora3_url(config, run_time, lead), subset; instantaneous = index != 0)

        if index != 0
            write_step!(output, index, slabs, previous_accumulated, subset)
        end

        previous_accumulated = slabs
    end

    return nothing
end

"""
    read_lead(url, subset::NORA3Subset; instantaneous = true)

One forecast lead's slabs over the subset, as source variable name => `Float32` matrix.

`instantaneous = false` reads only the accumulated fluxes, which is all that is needed of a lead
whose own hour lies outside the month being written.
"""
function read_lead(url, subset::NORA3Subset; instantaneous = true)
    return retry_open(url) do ds
        slabs = Dict{String,Matrix{Float32}}()

        for source_name in values(NORA3_ACCUMULATED_VARIABLES)
            slabs[source_name] = read_slab(ds, source_name, subset)
        end

        if instantaneous
            for source_name in values(NORA3_INSTANTANEOUS_VARIABLES)
                slabs[source_name] = read_slab(ds, source_name, subset)
            end
            slabs[NORA3_EASTWARD_SOURCE] = read_slab(ds, NORA3_EASTWARD_SOURCE, subset)
            slabs[NORA3_NORTHWARD_SOURCE] = read_slab(ds, NORA3_NORTHWARD_SOURCE, subset)
        end

        return slabs
    end
end

"""
    read_slab(ds, source_name, subset::NORA3Subset)

One horizontal slab of `source_name`, squeezing out the singleton height and time dimensions.
Every NORA3 surface variable is dimensioned `(x, y, height, time)` in Julia order, whatever its
height coordinate is called.
"""
function read_slab(ds, source_name, subset::NORA3Subset)
    data = ds[source_name][subset.x_range, subset.y_range, 1, 1]

    return Float32.(coalesce.(data, NaN32))
end

"""
    write_step!(output, index, slabs, previous_accumulated, subset::NORA3Subset)

Write one hourly step: the instantaneous fields as they come, the accumulated fluxes as the mean
rate over the hour ending at this step, and the winds rotated to eastward/northward.
"""
function write_step!(output, index, slabs, previous_accumulated, subset::NORA3Subset)
    for (name, source_name) in NORA3_INSTANTANEOUS_VARIABLES
        output[name][:, :, index] = slabs[source_name]
    end

    for (name, source_name) in NORA3_ACCUMULATED_VARIABLES
        output[name][:, :, index] =
            (slabs[source_name] .- previous_accumulated[source_name]) ./ NORA3_SECONDS_PER_LEAD
    end

    eastward, northward = rotate_to_east_north(
        subset.angle,
        slabs[NORA3_EASTWARD_SOURCE],
        slabs[NORA3_NORTHWARD_SOURCE],
    )
    output[NORA3_EASTWARD_NAME][:, :, index] = eastward
    output[NORA3_NORTHWARD_NAME][:, :, index] = northward

    return nothing
end
