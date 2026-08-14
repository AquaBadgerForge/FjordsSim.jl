module Atmospheres

export prepare_atmosphere,
    download_atmosphere,
    atmosphere_time_steps,
    atmosphere_source_grid,
    atmosphere_variable_names,
    atmosphere_target_axes,
    prescribed_atmosphere,
    prescribed_radiation,
    atmosphere_date_range,
    rotate_to_east_north,
    AtmosphereRecord,
    ProjectedAtmosphereGrid,
    NORA3Config,
    NORA3PrescribedAtmosphere,
    NORA3PrescribedRadiation,
    MultiYearNORA3

using Oceananigans
using Oceananigans.Grids: x_domain, y_domain
using ArchGDAL
using Dates: DateTime, Hour
import Dates
using NCDatasets
using NumericalEarth.DataWrangling: first_date, last_date

using ..Configs:
    AbstractAtmosphereConfig,
    FjordConfig,
    atmosphere_path,
    atmosphere_directory,
    coverage_window,
    domain_grid
using ..Plotting: plot_atmosphere

include("NORA3.jl")

using .NORA3: NORA3PrescribedAtmosphere, NORA3PrescribedRadiation, MultiYearNORA3

# Compression level for the prepared file, matching the forcing writer.
const ATMOSPHERE_DEFLATE_LEVEL = 5

# The prepared file's variables, with the units every source must deliver them in.
#
# This is a contract, not a style choice: the simulation-time reader resolves variables through
# `NORA3.NORA3_dataset_variable_names` and NumericalEarth consumes the values directly, so both the
# names and the units are fixed for every atmosphere source. In particular air temperature is
# Kelvin — NumericalEarth's `default_atmosphere_tracers` initializes it to `273.15 + 20` and the
# value feeds `Thermodynamics.air_density` and `latent_heat_vapor`, which are Kelvin-only — and
# both radiative fluxes are *downwelling*, because `NORA3PrescribedRadiation` applies its own
# surface albedo and would count it twice given a net flux.
const ATMOSPHERE_VARIABLES = (
    (name = "u_wind_10m", units = "m s-1", long_name = "Eastward component of 10 metre wind"),
    (name = "v_wind_10m", units = "m s-1", long_name = "Northward component of 10 metre wind"),
    (name = "air_temperature_2m", units = "K", long_name = "Screen level air temperature"),
    (name = "specific_humidity_2m", units = "kg kg-1", long_name = "Screen level specific humidity"),
    (name = "air_pressure_at_sea_level", units = "Pa", long_name = "Mean sea level pressure"),
    (name = "precipitation", units = "kg m-2 s-1", long_name = "Precipitation flux"),
    (name = "swrad", units = "W m-2", long_name = "Surface downwelling shortwave radiation"),
    (name = "lwrad", units = "W m-2", long_name = "Surface downwelling longwave radiation"),
)

"""
    atmosphere_variable_attributes(name)

The `units` and `long_name` the prepared file records for `name`, or no attributes for a variable
outside `ATMOSPHERE_VARIABLES`.
"""
function atmosphere_variable_attributes(name)
    index = findfirst(variable -> variable.name == name, ATMOSPHERE_VARIABLES)
    isnothing(index) && return Pair{String,String}[]
    variable = ATMOSPHERE_VARIABLES[index]

    return ["units" => variable.units, "long_name" => variable.long_name]
end

"""
    AtmosphereRecord

One time record present in a downloaded atmosphere file.

The download normalizes the dataset's forecast structure away, so by the time
`prepare_atmosphere` sees the data every record is one hour of one file and nothing about the
source's run/lead layout survives.

# Fields
- `date`: Time of the record.
- `filepath`: Downloaded file holding it.
- `index`: Time index within that file.
"""
struct AtmosphereRecord
    date::DateTime
    filepath::String
    index::Int
end

"""
    ProjectedAtmosphereGrid

A downloaded atmosphere subset that lives on a regular grid in a projected coordinate system,
plus the projection it is defined in.

Nothing here is dataset-specific, so any source on a regular projected grid can return one of
these from `atmosphere_source_grid` and inherit `projected_atmosphere_nodes` and
`interpolate_to_target!` unchanged. There is no vertical axis: the prepared variables are all
surface or screen-level fields.

# Fields
- `x`, `y`: Projected coordinate centers in meters, regularly spaced.
- `proj4`: PROJ.4 definition of the projection.
"""
struct ProjectedAtmosphereGrid
    x::Vector{Float64}
    y::Vector{Float64}
    proj4::String
end

# --- Extension hooks ---

"""
    atmosphere_time_steps(config)

The time records present in the downloaded files, as a sorted `Vector{AtmosphereRecord}` with
one entry per prepared time step.

Required for every `AbstractAtmosphereConfig`. See the `NORA3Config` method in
`src/Atmospheres/nora3_source.jl`.
"""
function atmosphere_time_steps end

"""
    atmosphere_source_grid(config, filepath)

Geometry of the downloaded data in `filepath`, e.g. a `ProjectedAtmosphereGrid`.

Required for every `AbstractAtmosphereConfig`. See the `NORA3Config` method in
`src/Atmospheres/nora3_source.jl`.
"""
function atmosphere_source_grid end

"""
    atmosphere_variable_names(config)

`Dict` mapping downloaded variable name to prepared variable name. The prepared names are the
ones `FjordSim.Atmospheres.NORA3.MultiYearNORA3` reads, so they are a fixed contract rather than
a free choice.

Required for every `AbstractAtmosphereConfig`. See the `NORA3Config` method in
`src/Atmospheres/nora3_source.jl`.
"""
function atmosphere_variable_names end

"""
    download_atmosphere(config::FjordConfig)
    download_atmosphere(target_grid, config)

Fetch and subset the atmosphere dataset a setup names, into `atmosphere_directory(config)`.

The `FjordConfig` method builds the setup's grid on the CPU and dispatches on the atmosphere
config, so a dataset only implements `download_atmosphere(target_grid, config)`. The grid rather
than the grid config is passed on, because a dataset needs the domain bounds and `x_domain` and
`y_domain` give those for any `AbstractGridConfig`.

The bounds it needs are the *padded* ones from `atmosphere_target_axes`, not the ocean domain:
the prepared grid deliberately overhangs the ocean grid, so a download subset to the ocean
domain alone would leave the outer prepared cells with nothing to interpolate from.

A setup naming no atmosphere config is a no-op.
"""
download_atmosphere(config::FjordConfig) =
    download_atmosphere(domain_grid(config.grid_config, CPU()), config.atmosphere_config)

download_atmosphere(target_grid, ::Nothing) = nothing

"""
    prescribed_atmosphere(config, architecture; reference_date = nothing)
    prescribed_radiation(config, architecture; reference_date = nothing)

Read the prepared atmosphere file back at simulation time, as the NumericalEarth
`PrescribedAtmosphere` and `PrescribedRadiation` that `coupled_simulation` consumes.

These are the only atmosphere hooks the simulation side needs, and they are what keeps
`FjordSim.Simulations` from naming a dataset: which reader, which file and which backend are all
the source's business. A setup naming no atmosphere config gets `nothing` for both, which is what
`OceanSeaIceModel` takes to mean an uncoupled run.

`reference_date` is the instant the returned time axes are zeroed at, `nothing` meaning the
prepared file's own first record. `build_simulation` passes the simulation config's `start_date`,
so the atmosphere and the ocean forcing agree on what model time zero stands for rather than each
zeroing at its own first record.

Deliberately no float-type argument: both NumericalEarth constructors default to `Float32`, and
passing `Oceananigans.defaults.FloatType` would silently promote the atmosphere to `Float64`.

Required for every `AbstractAtmosphereConfig` belonging to a setup that is simulated. See the
`NORA3Config` methods in `src/Atmospheres/nora3_source.jl`.
"""
prescribed_atmosphere(::Nothing, architecture; reference_date = nothing) = nothing
prescribed_radiation(::Nothing, architecture; reference_date = nothing) = nothing

"""
    atmosphere_date_range(config)

First and last date of the prepared atmosphere file, as a tuple, for `build_simulation`'s check
that the run does not outlast its atmosphere.

Optional: the default is `nothing`, meaning "this source cannot report its dates", which skips
the check rather than blocking a simulation. It exists as a hook rather than a direct read in
`FjordSim.Simulations` for the same reason `prescribed_atmosphere` does — resolving the prepared
file is the source's business, so the simulation module names no dataset.
"""
atmosphere_date_range(::Nothing) = nothing
atmosphere_date_range(::AbstractAtmosphereConfig) = nothing

# --- Target grid ---

"""
    atmosphere_target_axes(target_grid, config)

Longitude and latitude cell centers of the prepared atmosphere grid: the longitude/latitude
domain of `target_grid` grown by `config.padding` degrees and sampled at `config.resolution`,
snapped to whole multiples of the resolution.

The bounds come from `x_domain`/`y_domain` rather than a grid config's own fields, so any
`AbstractGridConfig` works. The atmosphere grid is deliberately coarser and larger than the
ocean grid: NumericalEarth interpolates from it, so it only has to cover the ocean domain, and
the padding is what keeps every ocean cell strictly inside it.

Both axes are uniformly spaced by construction, which `Utils.compute_faces` relies on when the
prepared file is read back.
"""
function atmosphere_target_axes(target_grid, config::AbstractAtmosphereConfig)
    config.resolution > 0 || throw(ArgumentError("resolution must be positive, got $(config.resolution)"))
    config.padding >= 0 || throw(ArgumentError("padding must not be negative, got $(config.padding)"))

    longitude = uniform_centers(x_domain(target_grid), config.padding, config.resolution)
    latitude = uniform_centers(y_domain(target_grid), config.padding, config.resolution)

    return longitude, latitude
end

"""
    uniform_centers(domain, padding, resolution)

Uniformly spaced centers at `resolution` covering `domain` grown by `padding` on both sides,
starting at a whole multiple of `resolution`.
"""
function uniform_centers(domain, padding, resolution)
    lower = domain[1] - padding
    upper = domain[2] + padding
    first_center = floor(lower / resolution) * resolution
    count = ceil(Int, (upper - first_center) / resolution) + 1

    return collect(range(first_center, step = resolution, length = count))
end

"""
    projected_atmosphere_nodes(longitude, latitude, source::ProjectedAtmosphereGrid)

The prepared grid's nodes projected into the source projection, as
`(length(longitude), length(latitude))`-shaped matrices. The whole point cloud is transformed in
a single GDAL call, as in `FjordSim.Forcing.projected_target_nodes`.
"""
function projected_atmosphere_nodes(longitude, latitude, source::ProjectedAtmosphereGrid)
    point_count = length(longitude) * length(latitude)
    x = Vector{Float64}(undef, point_count)
    y = Vector{Float64}(undef, point_count)

    point = 1
    for φ in latitude, λ in longitude
        x[point] = λ
        y[point] = φ
        point += 1
    end

    ArchGDAL.importEPSG(4326; order = :trad) do source_srs
        ArchGDAL.importPROJ4(source.proj4) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                ArchGDAL.transform!(x, y, zeros(Float64, point_count), transform)
            end
        end
    end

    shape = (length(longitude), length(latitude))
    return reshape(x, shape), reshape(y, shape)
end

"""
    validate_target_coverage(x, y, source::ProjectedAtmosphereGrid)

Check that every projected target node lies inside the downloaded subset, so the interpolation
never has to extrapolate. Fails loudly rather than silently clamping, because the usual cause is
a download subset that predates a change to `resolution` or `padding`.
"""
function validate_target_coverage(x, y, source::ProjectedAtmosphereGrid)
    x_range = extrema(source.x)
    y_range = extrema(source.y)
    inside = all(x_range[1] .<= x .<= x_range[2]) && all(y_range[1] .<= y .<= y_range[2])

    inside || error(
        "The prepared atmosphere grid reaches outside the downloaded subset " *
        "(target x $(extrema(x)), y $(extrema(y)); source x $x_range, y $y_range). " *
        "Re-run the atmosphere download for this config.",
    )

    return nothing
end

# --- Interpolation ---

"""
    interpolate_to_target!(output, slab, x, y, source::ProjectedAtmosphereGrid)

Bilinearly interpolate one source slab onto the projected target nodes `x`, `y`.

The source is regular in projected meters, so a target node's fractional source index is
arithmetic and no search is needed. This is deliberately not
`FjordSim.Forcing.interpolate_to_target!`: that one is trilinear and mask-driven, built for 3D
ocean fields on a `RectilinearGrid` with land to exclude, whereas atmospheric fields are 2D and
defined everywhere, including over land. The target grid here is also small enough (order
50 x 60) that there is nothing worth moving to a GPU.
"""
function interpolate_to_target!(output, slab, x, y, source::ProjectedAtmosphereGrid)
    Δx = source.x[2] - source.x[1]
    Δy = source.y[2] - source.y[1]
    Nx = length(source.x)
    Ny = length(source.y)
    x_origin = source.x[1]
    y_origin = source.y[1]

    for j in axes(output, 2), i in axes(output, 1)
        fractional_i = (x[i, j] - x_origin) / Δx + 1
        fractional_j = (y[i, j] - y_origin) / Δy + 1
        lower_i = clamp(floor(Int, fractional_i), 1, Nx - 1)
        lower_j = clamp(floor(Int, fractional_j), 1, Ny - 1)
        weight_i = fractional_i - lower_i
        weight_j = fractional_j - lower_j

        value =
            (1 - weight_i) * (1 - weight_j) * slab[lower_i, lower_j] +
            weight_i * (1 - weight_j) * slab[lower_i+1, lower_j] +
            (1 - weight_i) * weight_j * slab[lower_i, lower_j+1] +
            weight_i * weight_j * slab[lower_i+1, lower_j+1]
        output[i, j] = convert(Float32, value)
    end

    return output
end

# --- Wind rotation ---

"""
    grid_rotation_angle(longitude, latitude)

Angle from east to the source grid's local x axis, one per source cell, from finite differences
of the 2D `longitude` and `latitude` fields.

Source datasets on a projected grid give wind components relative to their own axes, which have
to be rotated to eastward/northward before the fields mean anything geographically. Deriving the
angle from the coordinate fields keeps this source-agnostic: it works for any curvilinear grid
publishing 2D longitude and latitude, without knowing the projection. The analytic alternative
for a Lambert conformal conic, `(longitude - longitude_0) * sin(latitude_1)`, needs the
projection parameters and would not carry over to another dataset.

The difference is taken along the first (x) axis, and the last column repeats the previous
difference.
"""
function grid_rotation_angle(longitude, latitude)
    Nx, Ny = size(longitude)
    Nx >= 2 || throw(ArgumentError("need at least two source columns to derive a rotation angle"))
    angle = Matrix{Float64}(undef, Nx, Ny)

    for j = 1:Ny, i = 1:Nx
        upper = min(i + 1, Nx)
        lower = upper - 1
        Δlongitude = (longitude[upper, j] - longitude[lower, j]) * cosd(latitude[i, j])
        Δlatitude = latitude[upper, j] - latitude[lower, j]
        angle[i, j] = atan(Δlatitude, Δlongitude)
    end

    return angle
end

"""
    rotate_to_east_north(angle, u, v)

Rotate grid-relative velocity components `u`, `v` into eastward and northward components, given
the `grid_rotation_angle` of each cell.
"""
function rotate_to_east_north(angle, u, v)
    eastward = similar(u)
    northward = similar(v)

    for index in eachindex(u, v, angle)
        cos_angle = cos(angle[index])
        sin_angle = sin(angle[index])
        eastward[index] = u[index] * cos_angle - v[index] * sin_angle
        northward[index] = v[index] * cos_angle + u[index] * sin_angle
    end

    return eastward, northward
end

# --- Preparation ---

"""
    prepare_atmosphere(target_grid, config::AbstractAtmosphereConfig)

Regrid the downloaded atmosphere files onto a regular longitude/latitude grid covering
`target_grid`, and write the single NetCDF file the simulation reads.

The pipeline: `atmosphere_variable_names(config)` picks the variables →
`atmosphere_time_steps(config)` gives the hourly records → `atmosphere_source_grid(config, …)`
gives the downloaded geometry → `atmosphere_target_axes` derives the prepared grid →
`projected_atmosphere_nodes` projects its nodes into the source projection → one bilinear
`interpolate_to_target!` per variable per step → streaming NetCDF write. Only the three hooks
are dataset-specific.

Everything depending on the source's forecast structure — which forecast lead supplies which
hour, de-accumulating fluxes within a run, unit conversion, wind rotation — happens in the
download, so this step is a pure regrid and re-running it at a different `resolution` costs no
re-download.

# The prepared file

Its layout is what `FjordSim.Atmospheres.NORA3.MultiYearNORA3` reads and is therefore fixed:
`Float32` variables of shape `(longitude, latitude, time)`, uniformly spaced `Float64` `lon` and
`lat` cell centers, and a CF-encoded `time`.

# Coverage

`coverage` is the `(first, last)` calendar interval the run needs, from `coverage_window`. The time
axis is padded to reach both ends by replicating the nearest record; `nothing` prepares exactly the
downloaded range. See `pad_atmosphere_records`.

# Returns

A named tuple `(; output_file, times, variables)`.
"""
function prepare_atmosphere(target_grid, config::AbstractAtmosphereConfig; coverage = nothing)
    variable_names = atmosphere_variable_names(config)
    isempty(variable_names) && error("$(nameof(typeof(config))) names no atmosphere variables to prepare")
    # Sorted so the prepared file's variable order does not depend on Dict iteration order.
    variables = sort!([source_name => name for (source_name, name) in variable_names], by = pair -> pair.second)

    records = atmosphere_time_steps(config)
    validate_atmosphere_records(records)
    records = pad_atmosphere_records(records, coverage)
    source = atmosphere_source_grid(config, first(records).filepath)
    validate_source_consistency(config, records, source)

    longitude, latitude = atmosphere_target_axes(target_grid, config)
    x, y = projected_atmosphere_nodes(longitude, latitude, source)
    validate_target_coverage(x, y, source)

    @info "Regridding $(length(records)) hourly steps from $(first(records).date) to $(last(records).date)"
    @info "Prepared grid: $(length(longitude)) x $(length(latitude)) cells at $(config.resolution) degrees"

    output_file = atmosphere_path(config)
    write_atmosphere_file(output_file, longitude, latitude, variables, records, source, x, y)

    return (;
        output_file,
        times = [record.date for record in records],
        variables = [pair.second for pair in variables],
    )
end

prepare_atmosphere(target_grid, ::Nothing; coverage = nothing) = nothing

"""
    prepare_atmosphere(config::FjordConfig)

Regrid the atmosphere a whole setup names onto a regular longitude/latitude grid, and write the
diagnostic plot. Returns `nothing` when the setup names no atmosphere.

This is the setup-level driver, the same shape as `download_atmosphere(config::FjordConfig)`. It
needs no bathymetry: the prepared grid is derived from the setup's longitude/latitude domain
rather than from the ocean grid's land mask, so unlike `prepare_forcing` it does not depend on
`prepare_bathymetry`. It does read the files `download_atmosphere` wrote.

The coverage window comes from the setup's simulation config, so the prepared file spans the run the
setup describes; a setup naming no simulation config gets `nothing` and the downloaded range.

# Returns
The `prepare_atmosphere(target_grid, config)` named tuple with `plot_file` added.
"""
function prepare_atmosphere(config::FjordConfig)
    atmosphere_config = config.atmosphere_config
    isnothing(atmosphere_config) && return nothing

    source_directory = atmosphere_directory(atmosphere_config)
    isdir(source_directory) || error(
        "Atmosphere source directory $source_directory does not exist. " *
        "Run `julia --project -m FjordSim download_atmosphere` for this setup first.",
    )

    grid = domain_grid(config.grid_config, CPU())
    result = prepare_atmosphere(
        grid,
        atmosphere_config;
        coverage = coverage_window(config.simulation_config),
    )
    plot_file = plot_atmosphere(atmosphere_config)

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Atmosphere file saved to $(result.output_file)"
    @info "Atmosphere plot saved to $plot_file"

    return (; result..., plot_file)
end

"""
    validate_atmosphere_records(records)

Check the record axis is non-empty, sorted and gap-free hourly, warning about any gap. A gap is
not fatal — `MultiYearNORA3` reads whatever axis it finds — but it means a forecast run failed to
download and the simulation will interpolate across the hole.
"""
function validate_atmosphere_records(records)
    isempty(records) &&
        error("No atmosphere time records found. Run the atmosphere download for this config first.")
    issorted(records, by = record -> record.date) || error("Atmosphere records must be sorted by date")

    for index = 2:length(records)
        span = records[index].date - records[index-1].date
        span == Hour(1) && continue
        @warn "Gap in the atmosphere time axis" from = records[index-1].date to = records[index].date
    end

    return nothing
end

"""
    pad_atmosphere_records(records, coverage)

Extend `records` to span `coverage`, replicating the nearest record at each end that falls short. A
`nothing` coverage leaves `records` alone.

Applied after `validate_atmosphere_records`, whose hourly-gap warning would otherwise fire on a pad
that is deliberately off the hourly cadence — the run's window need not land on an hour of the
downloaded axis.

A pad points at the same file and index as the record it copies, so the written field is identical
to its neighbour and `validate_source_consistency`, which walks the records' distinct files, sees
nothing new. See `FjordSim.Forcing.pad_time_steps` for why both prepared files are padded rather
than the readers being made tolerant: both read with `Cyclical()` time indexing, which fills a
shortfall by wrapping to the far end of the year instead of failing.

A pad may reach at most one record spacing past the downloaded axis, for the reason spelled out
there: unbounded, it would manufacture the very coverage `Simulations.validate_time_coverage` exists
to verify.
"""
pad_atmosphere_records(records, ::Nothing) = records

function pad_atmosphere_records(records, coverage)
    first_needed, last_needed = coverage
    padded = collect(records)

    if first(padded).date > first_needed
        head = first(padded)
        spacing = atmosphere_record_spacing(padded)
        head.date - first_needed <= spacing || error(
            "The run starts at $first_needed but the downloaded atmosphere only begins at " *
            "$(head.date), more than one record spacing " *
            "($(Dates.canonicalize(spacing))) later. Download the preceding period, or move the " *
            "simulation config's `start_date` to " *
            "$(head.date - spacing) or later.",
        )
        pushfirst!(padded, AtmosphereRecord(first_needed, head.filepath, head.index))
        @info "Padded atmosphere start: replicated $(head.date) at $first_needed"
    end

    if last(padded).date < last_needed
        tail = last(padded)
        spacing = atmosphere_record_spacing(padded)
        last_needed - tail.date <= spacing || error(
            "The run ends at $last_needed but the downloaded atmosphere already stops at " *
            "$(tail.date), more than one record spacing " *
            "($(Dates.canonicalize(spacing))) earlier. Download the following period, or shorten " *
            "`stop_time` so the run ends by " *
            "$(tail.date + spacing).",
        )
        push!(padded, AtmosphereRecord(last_needed, tail.filepath, tail.index))
        @info "Padded atmosphere end: replicated $(tail.date) at $last_needed"
    end

    return padded
end

"""
    atmosphere_record_spacing(records)

The axis's record spacing, taken from its first pair — `validate_atmosphere_records` has already
reported any departure from hourly, so one pair speaks for the whole axis. A single-record axis has
no spacing to bound a pad by, and is an error rather than a guess.
"""
function atmosphere_record_spacing(records)
    length(records) >= 2 || error(
        "Cannot pad an atmosphere axis of $(length(records)) record(s): its record spacing, " *
        "which bounds how far a pad may reach, is unknown.",
    )
    return records[2].date - records[1].date
end

"""
    validate_source_consistency(config, records, source)

Check every downloaded file shares the geometry of the first one.

The interpolation is set up once from the first file, so a file cut to a different window would be
read with the wrong coordinates and produce a plausible-looking wrong answer. That is a live hazard
because the download skips months already present: changing `resolution` or `padding` and re-running
it leaves the earlier months on the old window.
"""
function validate_source_consistency(config, records, source)
    for filepath in unique(record.filepath for record in records)
        other = atmosphere_source_grid(config, filepath)
        (other.x == source.x && other.y == source.y && other.proj4 == source.proj4) || error(
            "$filepath was downloaded on a different source window than " *
            "$(first(records).filepath) ($(length(other.x))x$(length(other.y)) against " *
            "$(length(source.x))x$(length(source.y))). Delete $(atmosphere_directory(config)) " *
            "and re-run the atmosphere download.",
        )
    end

    return nothing
end

"""
    AtmosphereReader(filepath)

Holds the one downloaded file currently open, reopening on demand. Records are written in date
order, so a month's worth of steps costs one open rather than one per step. Opened eagerly so the
`dataset` field stays concretely typed.
"""
mutable struct AtmosphereReader{D}
    filepath::String
    dataset::D
end

AtmosphereReader(filepath::String) = AtmosphereReader(filepath, NCDataset(filepath))

Base.close(reader::AtmosphereReader) = close(reader.dataset)

"""
    source_slab(reader::AtmosphereReader, record::AtmosphereRecord, source_name)

One horizontal slab of `source_name` at `record`, as `Float32` with missing values turned into
`NaN`.
"""
function source_slab(reader::AtmosphereReader, record::AtmosphereRecord, source_name)
    if reader.filepath != record.filepath
        close(reader.dataset)
        reader.dataset = NCDataset(record.filepath)
        reader.filepath = record.filepath
    end

    return Float32.(coalesce.(reader.dataset[source_name][:, :, record.index], NaN32))
end

"""
    write_atmosphere_file(filepath, longitude, latitude, variables, records, source, x, y)

Stream the prepared file: one interpolated slab per variable per record.

Chunked one horizontal slab per time step, as in `FjordSim.Forcing.write_forcing_file` — the
default chunking spreads a single time step across the whole file and made that writer an order
of magnitude slower.
"""
function write_atmosphere_file(filepath, longitude, latitude, variables, records, source, x, y)
    isdir(dirname(filepath)) || mkpath(dirname(filepath))
    isfile(filepath) && rm(filepath)

    Nx = length(longitude)
    Ny = length(latitude)
    output = Matrix{Float32}(undef, Nx, Ny)
    reader = AtmosphereReader(first(records).filepath)

    ds = NCDataset(filepath, "c")
    try
        define_atmosphere_dimensions!(ds, longitude, latitude, records)

        for (_, name) in variables
            defVar(
                ds,
                name,
                Float32,
                ("lon", "lat", "time");
                chunksizes = [Nx, Ny, 1],
                deflatelevel = ATMOSPHERE_DEFLATE_LEVEL,
                attrib = ["_FillValue" => NaN32; atmosphere_variable_attributes(name)],
            )
        end

        announced = ""
        for (index, record) in enumerate(records)
            if record.filepath != announced
                @info "Regridding $(basename(record.filepath))"
                announced = record.filepath
            end

            for (source_name, name) in variables
                slab = source_slab(reader, record, source_name)
                interpolate_to_target!(output, slab, x, y, source)
                # An atmospheric field is defined everywhere, so unlike the forcing writer there
                # is no land mask here to excuse a NaN. One reaching `PrescribedAtmosphere` would
                # poison the flux computation silently, so refuse to write it.
                all(isfinite, output) || error(
                    "Non-finite values in $name at $(record.date) " *
                    "($(count(!isfinite, output)) of $(length(output)) cells). " *
                    "The downloaded file $(record.filepath) is incomplete or corrupt.",
                )
                ds[name][:, :, index] = output
            end
        end
    finally
        close(ds)
        close(reader)
    end

    return filepath
end

"""
    define_atmosphere_dimensions!(ds, longitude, latitude, records)

Define the dimensions and coordinate variables the prepared file's reader expects: 1D `lon` and
`lat` cell centers, and `time` written as `DateTime`s so NCDatasets encodes CF units.
"""
function define_atmosphere_dimensions!(ds, longitude, latitude, records)
    defDim(ds, "lon", length(longitude))
    defDim(ds, "lat", length(latitude))
    defDim(ds, "time", length(records))

    defVar(ds, "lon", Float64, ("lon",))[:] = longitude
    defVar(ds, "lat", Float64, ("lat",))[:] = latitude
    defVar(ds, "time", [record.date for record in records], ("time",))

    return ds
end

include("nora3_source.jl")

end  # module Atmospheres
