# Hourly NorKyst-800m adapter: the source-specific half of the open-boundary pipeline. Copy this
# file's shape to add a new boundary dataset — it owns a config subtype, its own download, and the
# three hook methods `boundary_time_steps`, `boundary_source_grid` and `boundary_variable_names`.
# Everything else comes from the generic core in `boundaries.jl`.
#
# A *different* NorKyst collection from `norkyst.jl`'s. That one reads `fou-hi/norkyst800m`,
# whose files are `..._avg.an.*` — daily averages, one record each. This one reads
# `fou-hi/norkyst800m-1h`, whose files are `..._his.an.*` — 24 hourly instants each. The interior
# forcing wants the daily means; an open boundary wants the hourly instants, because a daily mean has
# the tide averaged out of it and elevation with no tide is not worth prescribing.
#
# NorKyst's own lateral boundary conditions, from its global attributes, are the same scheme set
# FjordSim now applies at its own edge: `zeta: Che` (Chapman), `ubar/vbar: Shc` (Flather-family),
# `u/v/temp/salt: RadNud` (Orlanski radiation plus nudging).

const NORKYST_HOURLY_CATALOG_URL = "https://thredds.met.no/thredds/catalog/fou-hi/norkyst800m-1h/catalog.xml"
const NORKYST_HOURLY_OPENDAP_URL = "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/"

# NorKyst variable names and the FjordSim boundary names they become. `zeta`, `ubar` and `vbar` are
# the three surface variables; the other four are full-depth.
const NORKYST_BOUNDARY_VARIABLE_NAMES = Dict(
    "temperature" => "T",
    "salinity" => "S",
    "u_eastward" => "u",
    "v_northward" => "v",
    "zeta" => "eta",
    "ubar" => "ubar",
    "vbar" => "vbar",
)

"""
    NorKystBoundariesConfig

Configuration for downloading and subsetting hourly NorKyst-800m data along one open lateral
boundary.

`output_directory`, `output_file` and `plot_file` are names relative to `data_root`, resolved by
`boundary_data_directory`, `boundary_data_path` and `plot_path`. Setting one to an absolute path
overrides `data_root` for that entry only.

# Fields
- `data_root`: Directory holding this setup's boundary files. Required.
- `output_directory`: Name of the directory the monthly downloads are written to. Required.
- `output_file`: Name of the prepared boundary NetCDF written by `prepare_boundaries`.
- `plot_file`: Name of the diagnostic boundary plot.
- `margin`: Degrees of latitude/longitude the download reaches past the boundary row, on both
  sides. Must leave room for the source cells the bilinear interpolation reads around the row —
  NorKyst's spacing is 800 m, so a few hundredths of a degree. Kept small because the download is
  hourly: the whole domain box at 24 records a day is more than twenty times the interior forcing,
  while a thin band is a tenth of the box.
- `architecture`: Where `prepare_boundaries` interpolates — `:auto`, `:cpu` or `:gpu`, resolved by
  `interpolation_architecture`.
- `catalog_url`, `opendap_url`: The hourly NorKyst collection's own public endpoints; the only
  defaulted fields.
- `parameters`: Source variable names to extract. Required.
- `years`: Calendar years to download. Required.
"""
Base.@kwdef mutable struct NorKystBoundariesConfig <: AbstractBoundaryDataConfig
    data_root::String
    output_directory::String
    output_file::String = "boundaries.nc"
    plot_file::String = "boundaries.png"
    margin::Float64 = 0.05
    architecture::Symbol = :auto
    catalog_url::String = NORKYST_HOURLY_CATALOG_URL
    opendap_url::String = NORKYST_HOURLY_OPENDAP_URL
    parameters::Vector{String}
    years::Vector{Int}
end

"""
    boundary_monthly_filename(config::NorKystBoundariesConfig, year, month)

Name of the combined monthly hourly-NorKyst NetCDF written for `year` and `month`.
"""
boundary_monthly_filename(config::NorKystBoundariesConfig, year, month) =
    "NorKyst-800m_ZDEPTHS_his_$(year)$(lpad(month, 2, '0')).nc"

"""
    boundary_variable_names(config::NorKystBoundariesConfig)

The NorKyst variables this dataset can supply along a boundary and the FjordSim names they become.
"""
boundary_variable_names(config::NorKystBoundariesConfig) = NORKYST_BOUNDARY_VARIABLE_NAMES

"""
    boundary_time_steps(config::NorKystBoundariesConfig)

Every time record of every downloaded monthly file for `config.years`, sorted by date with
duplicates dropped. Errors if the directory or the files are missing.
"""
function boundary_time_steps(config::NorKystBoundariesConfig)
    directory = boundary_data_directory(config)
    isdir(directory) || error(
        "Hourly NorKyst directory $directory does not exist. " *
        "Run `julia --project -m FjordSim download_boundaries` for this setup first.",
    )

    records = SourceRecord[]
    for year in config.years, month = 1:12
        filepath = joinpath(directory, boundary_monthly_filename(config, year, month))
        isfile(filepath) || continue
        NCDataset(filepath) do ds
            for (index, date) in enumerate(ds["time"][:])
                push!(records, SourceRecord(DateTime(date), filepath, index))
            end
        end
    end

    isempty(records) &&
        error("No hourly NorKyst monthly files for years $(config.years) found in $directory.")
    sort!(records; by = record -> record.date)

    return unique(record -> record.date, records)
end

"""
    boundary_source_grid(config::NorKystBoundariesConfig, filepath)

Read the projected coordinates, depth levels and projection of a downloaded hourly NorKyst subset.
Identical in shape to `forcing_source_grid(config::NorKystConfig, filepath)` — the two collections
share their grid, their projection variable and their depth axis.
"""
function boundary_source_grid(config::NorKystBoundariesConfig, filepath)
    return NCDataset(filepath) do ds
        x = Array{Float64}(ds["X"][:])
        y = Array{Float64}(ds["Y"][:])
        depths = Array{Float64}(ds["depth"][:])
        proj4 = NCDatasets.variable(ds, NORKYST_PROJECTION_VARIABLE).attrib["proj4"]

        length(x) >= 2 && length(y) >= 2 || error(
            "Hourly NorKyst subset in $filepath is too small to interpolate from: " *
            "$(length(x))x$(length(y)). Raise the boundary config's `margin`.",
        )
        all(difference -> isapprox(difference, x[2] - x[1]), diff(x)) &&
            all(difference -> isapprox(difference, y[2] - y[1]), diff(y)) ||
            error("Hourly NorKyst projected coordinates in $filepath are not regularly spaced.")

        return ProjectedSourceGrid(x, y, depths, proj4)
    end
end

# --- Download ---

"""
    download_boundaries(target_grid, edge, config::NorKystBoundariesConfig)

Download the hourly NorKyst-800m months covering `config.years`, each combined into one NetCDF file
in `boundary_data_directory(config)` subset to a `config.margin`-wide lon/lat band along `edge` of
`target_grid`. A month whose file already exists is skipped, so an interrupted download resumes.

The band, rather than the whole domain box, is what keeps an hourly download affordable — see
`boundary_domain`.
"""
function download_boundaries(target_grid, edge, config::NorKystBoundariesConfig)
    output_directory = boundary_data_directory(config)
    mkpath(output_directory)

    longitude, latitude = boundary_domain(Val(validate_open_edge(edge)), target_grid, config.margin)
    @info "Downloading hourly NorKyst-800m years $(join(config.years, ", ")) along the :$edge " *
          "boundary to $output_directory"
    @info "  Band: longitude $longitude, latitude $latitude"

    files = list_opendap_files(catalog_url = config.catalog_url)
    for year in config.years, month = 1:12
        process_boundary_month(year, month, longitude, latitude, config; files)
    end

    @info "Finished downloading hourly NorKyst-800m boundary data"
    return output_directory
end

"""
    boundary_subset_ranges(ds, longitude, latitude, config::NorKystBoundariesConfig)

The NorKyst index window covering the `longitude`/`latitude` band, as the same `NorKystSubset` the
interior download uses — so `define_output_file`, `write_parameter_chunk!` and
`write_time_dependent_coordinates!` are reused verbatim.

The band is thin in geographic space but NorKyst's grid is rotated about 59 degrees from east here,
so its bounding index box is not thin in either X or Y. It is still an order of magnitude smaller
than the whole-domain box, which is the point.
"""
function boundary_subset_ranges(ds, longitude, latitude, config::NorKystBoundariesConfig)
    latitude_variable = variable(ds, "lat")
    longitude_variable = variable(ds, "lon")
    source_latitude = Array(latitude_variable[ntuple(_ -> :, ndims(latitude_variable))...])
    source_longitude = Array(longitude_variable[ntuple(_ -> :, ndims(longitude_variable))...])

    mask = (source_latitude .>= latitude[1]) .&
           (source_latitude .<= latitude[2]) .&
           (source_longitude .>= longitude[1]) .&
           (source_longitude .<= longitude[2])

    ranges = Dict{String,UnitRange{Int}}()
    for (index, dimension) in enumerate(dimnames(latitude_variable))
        ranges[dimension] = bounding_range(mask, index)
    end

    spatial_dimensions = dimnames(latitude_variable)
    subset_mask = mask[(ranges[dimension] for dimension in spatial_dimensions)...]
    return NorKystSubset(ranges, subset_mask, spatial_dimensions, config.parameters)
end

"""
    process_boundary_month(year, month, longitude, latitude, config; files)

Download one month of hourly NorKyst into a single NetCDF, subset to the boundary band.

The interior download's `process_month` in shape, and it shares every helper below the driver, but a
separate function: this collection has its own catalog, its own filename pattern and 24 records per
file rather than one, and folding both into one driver would mean a config-shaped branch inside it.
"""
function process_boundary_month(year, month, longitude, latitude, config::NorKystBoundariesConfig; files = nothing)
    files = isnothing(files) ? list_opendap_files(catalog_url = config.catalog_url) : files
    month_string = ".$(year)$(lpad(month, 2, '0'))"
    year_month = month_string[2:end]
    output_path = joinpath(boundary_data_directory(config), boundary_monthly_filename(config, year, month))

    if isfile(output_path)
        @info "Skipping $year_month (already exists)"
        return output_path
    end

    @info "Processing $year_month..."
    monthly_files = sort([file for file in files if occursin(month_string, file)])

    if isempty(monthly_files)
        @info "  No files found for $year_month, skipping."
        return nothing
    end

    urls = [joinpath(config.opendap_url, file) for file in monthly_files]
    @info "  Opening $(length(urls)) datasets..."

    datasets = NCDataset[]
    try
        for url in urls
            push!(datasets, NCDataset(url))
        end

        subset = boundary_subset_ranges(first(datasets), longitude, latitude, config)
        total_time = sum(time_length, datasets)
        @info "  Subset $(join(("$dimension=$(length(range))" for (dimension, range) in sort(collect(subset.ranges))), ", ")), $total_time records"

        @info "  Writing output to: $output_path"
        output = define_output_file(output_path, first(datasets), subset, total_time)
        try
            time_start = 1
            for ds in datasets
                write_time_dependent_coordinates!(output, ds, subset, time_start)
                for name in subset.parameters
                    write_parameter_chunk!(output, ds, name, subset, time_start)
                end
                time_start += time_length(ds)
            end
        finally
            close(output)
        end
    finally
        foreach(close, datasets)
    end

    @info "Finished $year_month"
    return output_path
end
