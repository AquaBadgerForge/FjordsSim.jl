# NorKyst-800m adapter: the source-specific half of the forcing pipeline. Copy this file's
# shape to add a new forcing dataset — it owns a config subtype, its own download, and the
# three hook methods `forcing_time_steps`, `forcing_source_grid` and `forcing_variable_names`
# (plus `download_forcing` if it fetches data). Everything else comes from the generic core
# in Forcing.jl.

using Downloads

const NORKYST_CATALOG_URL = "https://thredds.met.no/thredds/catalog/fou-hi/norkyst800m/catalog.xml"
const NORKYST_OPENDAP_URL = "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m/"

# NorKyst variable names and the FjordSim forcing names they become. Only the intersection of
# this mapping with `config.parameters` is prepared.
const NORKYST_VARIABLE_NAMES = Dict(
    "temperature" => "T",
    "salinity" => "S",
    "u_eastward" => "u",
    "v_northward" => "v",
)
# Scalar NorKyst variable whose `proj4` attribute defines the projected X/Y coordinates.
const NORKYST_PROJECTION_VARIABLE = "projection_stere"

"""
    NorKystConfig

Configuration for downloading and subsetting NorKyst-800m reanalysis data.

A setup states where its forcing goes, which variables to extract and which years to cover.
Only `catalog_url` and `opendap_url` are defaulted, being the NorKyst dataset's own public
endpoints.

`output_directory`, `output_file` and `plot_file` are names relative to `data_root`, resolved
by `forcing_directory`, `forcing_path` and `plot_path`. Setting one to an absolute path
overrides `data_root` for that entry only.

# Fields
- `data_root`: Directory holding this setup's forcing files. Required.
- `output_directory`: Name of the directory where monthly NetCDF files are written. Required.
- `output_file`: Name of the prepared forcing NetCDF written by `prepare_forcing`.
- `plot_file`: Name of the diagnostic forcing plot.
- `relaxation_edge`: Lateral boundary the forcing relaxes towards NorKyst on, one of
  `:south`, `:north`, `:west` or `:east`.
- `relaxation_cells`: Width of the relaxation band in grid cells.
- `relaxation_timescale`: Relaxation timescale in seconds; the written lambda is its inverse.
- `architecture`: Where `prepare_forcing` interpolates — `:auto`, `:cpu` or `:gpu`, resolved by
  `interpolation_architecture`. `:auto` keeps one config usable on a GPU machine and a laptop.
- `catalog_url`: THREDDS catalog URL listing available files.
- `opendap_url`: OPeNDAP base URL for streaming data.
- `parameters`: Variable names to extract (e.g. `["temperature", "salinity"]`). Required.
- `years`: Calendar years to download. Required.
"""
Base.@kwdef mutable struct NorKystConfig <: AbstractForcingConfig
    data_root::String
    output_directory::String
    output_file::String = "forcing.nc"
    plot_file::String = "forcing.png"
    relaxation_edge::Symbol = :south
    relaxation_cells::Int = 10
    relaxation_timescale::Float64 = 86400.0
    architecture::Symbol = :auto
    catalog_url::String = NORKYST_CATALOG_URL
    opendap_url::String = NORKYST_OPENDAP_URL
    parameters::Vector{String}
    years::Vector{Int}
end

"""
    forcing_monthly_filename(config::NorKystConfig, year, month)

Name of the combined monthly NorKyst-800m NetCDF file written for `year` and `month`.
"""
forcing_monthly_filename(config::NorKystConfig, year, month) =
    "NorKyst-800m_ZDEPTHS_avg_$(year)$(lpad(month, 2, '0')).nc"

"""
    forcing_variable_names(config::NorKystConfig)

The NorKyst variables this dataset can supply and the FjordSim forcing names they become.
"""
forcing_variable_names(config::NorKystConfig) = NORKYST_VARIABLE_NAMES

"""
    forcing_time_steps(config::NorKystConfig)

Every time record of every downloaded monthly file for `config.years`, sorted by date with
duplicates dropped. Errors if the directory or the files are missing.
"""
function forcing_time_steps(config::NorKystConfig)
    directory = forcing_directory(config)
    isdir(directory) || error(
        "NorKyst directory $directory does not exist. " *
        "Run scripts/forcing_download.jl for this config first.",
    )

    records = SourceRecord[]
    for year in config.years, month = 1:12
        filepath = joinpath(directory, forcing_monthly_filename(config, year, month))
        isfile(filepath) || continue
        NCDataset(filepath) do ds
            for (index, date) in enumerate(ds["time"][:])
                push!(records, SourceRecord(DateTime(date), filepath, index))
            end
        end
    end

    isempty(records) && error("No NorKyst monthly files for years $(config.years) found in $directory.")
    sort!(records; by = record -> record.date)

    return unique(record -> record.date, records)
end

"""
    forcing_source_grid(config::NorKystConfig, filepath)

Read the projected coordinates, depth levels and projection of a downloaded NorKyst subset.
Errors unless the projected coordinates are regularly spaced, which `source_field_grid` needs in
order to express them as a `RectilinearGrid`.
"""
function forcing_source_grid(config::NorKystConfig, filepath)
    return NCDataset(filepath) do ds
        x = Array{Float64}(ds["X"][:])
        y = Array{Float64}(ds["Y"][:])
        depths = Array{Float64}(ds["depth"][:])
        proj4 = NCDatasets.variable(ds, NORKYST_PROJECTION_VARIABLE).attrib["proj4"]

        length(x) >= 2 && length(y) >= 2 ||
            error("NorKyst subset in $filepath is too small to interpolate from: $(length(x))x$(length(y)).")
        all(difference -> isapprox(difference, x[2] - x[1]), diff(x)) &&
            all(difference -> isapprox(difference, y[2] - y[1]), diff(y)) ||
            error("NorKyst projected coordinates in $filepath are not regularly spaced.")

        return ProjectedSourceGrid(x, y, depths, proj4)
    end
end

# --- Download ---

"""
    download_forcing(target_grid, config::NorKystConfig)

Download the NorKyst-800m months covering `config.years`, each combined into one NetCDF file
in `forcing_directory(config)` subset to the lon/lat box of `target_grid`. A month whose file
already exists is skipped, so an interrupted download resumes.
"""
function download_forcing(target_grid, config::NorKystConfig)
    output_directory = forcing_directory(config)
    mkpath(output_directory)

    @info "Downloading NorKyst-800m years $(join(config.years, ", ")) to $output_directory"

    files = list_opendap_files(catalog_url = config.catalog_url)
    for year in config.years, month = 1:12
        process_month(year, month, target_grid, config; files)
    end

    @info "Finished downloading NorKyst-800m forcing"
    return output_directory
end

function list_opendap_files(; catalog_url)
    catalog_path = Downloads.download(catalog_url)
    catalog = read(catalog_path, String)

    files = String[]
    for match in eachmatch(r"<dataset\b[^>]*\bname=\"([^\"]+\.nc)\""i, catalog)
        push!(files, match.captures[1])
    end

    return files
end

function bounding_range(mask, dimension)
    reduced_dims = Tuple(i for i in 1:ndims(mask) if i != dimension)
    axis_mask = vec(any(mask; dims = reduced_dims))
    indices = findall(axis_mask)
    isempty(indices) && error("No NorKyst points found inside the requested lon/lat range.")
    return first(indices):last(indices)
end

"""
    NorKystSubset

The spatial window of a NorKyst dataset selected by a setup, plus the variables to extract.
Derived from a target grid and a `NorKystConfig` by `subset_ranges` and passed to every
function that writes the subset out.

# Fields
- `ranges`: Index range per dataset dimension covering the requested lon/lat box.
- `mask`: In-box flag per point of the subset, used to blank out points outside the box.
- `spatial_dimensions`: Names of the dimensions `mask` is indexed by, in order.
- `parameters`: Variable names to extract.
"""
struct NorKystSubset{M<:AbstractArray,S<:Tuple}
    ranges::Dict{String,UnitRange{Int}}
    mask::M
    spatial_dimensions::S
    parameters::Vector{String}
end

"""
    subset_ranges(ds, target_grid, config::NorKystConfig)

The NorKyst index window covering the lon/lat domain of `target_grid`.

The bounds come from `x_domain`/`y_domain` rather than a grid config's own fields, so any
`AbstractGridConfig` works; they are the same face bounds an `EvenGrid` is built from.
"""
function subset_ranges(ds, target_grid, config::NorKystConfig)
    longitude_range = x_domain(target_grid)
    latitude_range = y_domain(target_grid)

    latitude_variable = variable(ds, "lat")
    longitude_variable = variable(ds, "lon")
    latitude = Array(latitude_variable[ntuple(_ -> :, ndims(latitude_variable))...])
    longitude = Array(longitude_variable[ntuple(_ -> :, ndims(longitude_variable))...])
    mask = (latitude .>= latitude_range[1]) .&
           (latitude .<= latitude_range[2]) .&
           (longitude .>= longitude_range[1]) .&
           (longitude .<= longitude_range[2])

    ranges = Dict{String,UnitRange{Int}}()
    for (i, dimension) in enumerate(dimnames(latitude_variable))
        ranges[dimension] = bounding_range(mask, i)
    end

    spatial_dimensions = dimnames(latitude_variable)
    subset_mask = mask[(ranges[dimension] for dimension in spatial_dimensions)...]
    return NorKystSubset(ranges, subset_mask, spatial_dimensions, config.parameters)
end

function variable_indices(variable, ranges)
    return ntuple(ndims(variable)) do i
        dimension = dimnames(variable)[i]
        get(ranges, dimension, :)
    end
end

function time_dimension(ds)
    time_dimensions = [name for name in dimnames(ds) if occursin("time", lowercase(name))]
    !isempty(time_dimensions) && return first(time_dimensions)
    haskey(ds, "time") && return first(dimnames(ds["time"]))
    error("Could not find a time dimension in NorKyst dataset.")
end

function time_length(ds)
    return NCDatasets.dim(ds, time_dimension(ds))
end

function copy_attributes!(dest, source)
    for (key, value) in source.attrib
        key == "_FillValue" && continue
        dest.attrib[key] = value
    end
end

filtered_attributes(source) = [(key, value) for (key, value) in source.attrib if key != "_FillValue"]
decoded_attributes(source) = [
    (key, value) for (key, value) in source.attrib
    if key ∉ ("_FillValue", "missing_value", "scale_factor", "add_offset")
]

function concrete_float_data(data)
    T = nonmissingtype(eltype(data))
    T <: AbstractFloat || error("Expected floating-point NorKyst data, got $(eltype(data)).")
    output = Array{T}(undef, size(data))

    for index in eachindex(data)
        value = data[index]
        output[index] = ismissing(value) ? convert(T, NaN) : value
    end

    return output
end

function define_subset_variable(output, source, name, subset::NorKystSubset; deflatelevel = 5)
    variable = NCDatasets.variable(source, name)
    dimensions = dimnames(variable)
    indices = variable_indices(variable, subset.ranges)
    data = variable[indices...]
    variable_type = data isa AbstractArray ? eltype(data) : typeof(data)

    output_variable = defVar(
        output,
        name,
        variable_type,
        dimensions;
        deflatelevel,
        attrib = filtered_attributes(variable),
    )

    if time_dimension(source) ∉ dimensions
        if ndims(output_variable) == 0
            output_variable[] = data
        else
            output_variable[ntuple(_ -> :, ndims(output_variable))...] = data
        end
    end

    return output_variable
end

function copy_auxiliary_variable(name, variable, time_dim, subset::NorKystSubset)
    name in subset.parameters && return false
    dimensions = dimnames(variable)
    return name == time_dim || time_dim ∉ dimensions
end

function define_output_file(output_path, template, subset::NorKystSubset, total_time)
    time_dim = time_dimension(template)
    isfile(output_path) && rm(output_path; force = true)

    output = NCDataset(output_path, "c")
    try
        for dimension in dimnames(template)
            dimension_length = if dimension == time_dim
                total_time
            elseif haskey(subset.ranges, dimension)
                length(subset.ranges[dimension])
            else
                NCDatasets.dim(template, dimension)
            end
            defDim(output, dimension, dimension_length)
        end

        for name in keys(template)
            variable = NCDatasets.variable(template, name)
            copy_auxiliary_variable(name, variable, time_dim, subset) || continue
            all(dimension -> dimension in dimnames(output), dimnames(variable)) || continue
            define_subset_variable(output, template, name, subset; deflatelevel = 0)
        end

        for name in subset.parameters
            variable = template[name]
            variable_type = nonmissingtype(eltype(variable))
            defVar(
                output,
                name,
                variable_type,
                dimnames(variable);
                deflatelevel = 5,
                attrib = decoded_attributes(variable),
            )
        end

        copy_attributes!(output, template)
    catch
        close(output)
        rethrow()
    end

    return output
end

function masked_fill_value(data)
    return if eltype(data) <: AbstractFloat
        convert(eltype(data), NaN)
    else
        error("Cannot mask data with element type $(eltype(data)).")
    end
end

function apply_spatial_mask(data, variable, subset::NorKystSubset)
    spatial_indices = [findfirst(==(dimension), dimnames(variable)) for dimension in subset.spatial_dimensions]
    any(isnothing, spatial_indices) && return data

    masked_data = copy(data)
    fill_value = masked_fill_value(masked_data)

    for mask_index in CartesianIndices(subset.mask)
        subset.mask[mask_index] && continue
        data_indices = ntuple(ndims(masked_data)) do i
            mask_dimension = findfirst(==(i), spatial_indices)
            mask_dimension === nothing ? (:) : mask_index[mask_dimension]
        end
        masked_data[data_indices...] .= fill_value
    end

    return masked_data
end

function write_parameter_chunk!(output, source, name, subset::NorKystSubset, time_start)
    variable = source[name]
    input_indices = variable_indices(variable, subset.ranges)
    data = concrete_float_data(variable[input_indices...])
    data = apply_spatial_mask(data, variable, subset)

    output_variable = output[name]
    output_indices = ntuple(ndims(output_variable)) do i
        dimension = dimnames(output_variable)[i]
        if dimension == time_dimension(source)
            time_start:(time_start + size(data, i) - 1)
        else
            (:)
        end
    end

    output_variable[output_indices...] = data
    return size(data, findfirst(==(time_dimension(source)), dimnames(variable)))
end

function write_time_dependent_coordinates!(output, source, subset::NorKystSubset, time_start)
    time_dim = time_dimension(source)

    for name in keys(source)
        name in subset.parameters && continue
        variable = NCDatasets.variable(source, name)
        time_index = findfirst(==(time_dim), dimnames(variable))
        time_index === nothing && continue
        haskey(output, name) || continue

        input_indices = variable_indices(variable, subset.ranges)
        data = variable[input_indices...]
        output_variable = NCDatasets.variable(output, name)
        output_indices = ntuple(ndims(output_variable)) do i
            if i == time_index
                time_start:(time_start + size(data, i) - 1)
            else
                (:)
            end
        end
        output_variable[output_indices...] = data
    end
end

function process_month(year, month, target_grid, config::NorKystConfig; files = nothing)
    files = isnothing(files) ? list_opendap_files(catalog_url = config.catalog_url) : files
    month_string = ".$(year)$(lpad(month, 2, '0'))"
    year_month = month_string[2:end]
    output_path = joinpath(forcing_directory(config), forcing_monthly_filename(config, year, month))

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
            @info "    Opened: $url"
        end

        subset = subset_ranges(first(datasets), target_grid, config)
        total_time = sum(time_length, datasets)

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
