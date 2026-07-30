#!/usr/bin/env julia

using Downloads
using FjordSim
using NCDatasets

const DEFAULT_CONFIG_PATH = joinpath(@__DIR__, "..", "configs", "drammensfjorden.jl")

function parse_args(args = ARGS)
    config_path = DEFAULT_CONFIG_PATH

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i == length(args) && error("--config requires a value")
            config_path = args[i + 1]
            i += 2
        elseif startswith(arg, "--config=")
            config_path = split(arg, "=", limit = 2)[2]
            i += 1
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    config = include(abspath(config_path))
    return (
        config = config,
        config_path = abspath(config_path),
    )
end

function print_usage()
    println("""
    Download and combine NorKyst-800m monthly data for configured years.

    Usage:
      julia --project scripts/forcing_download_norkyst.jl [--config PATH]

    Options:
      --config PATH      Setup config (Julia file). Default: configs/drammensfjorden.jl
    """)
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

function subset_ranges(ds; latitude_range, longitude_range)
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
    return ranges, subset_mask, spatial_dimensions
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

function define_subset_variable(output, source, name, ranges; deflatelevel = 5)
    variable = NCDatasets.variable(source, name)
    dimensions = dimnames(variable)
    indices = variable_indices(variable, ranges)
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

function copy_auxiliary_variable(name, variable, time_dim, parameters)
    name in parameters && return false
    dimensions = dimnames(variable)
    return name == time_dim || time_dim ∉ dimensions
end

function define_output_file(output_path, template, ranges, total_time; parameters)
    time_dim = time_dimension(template)
    isfile(output_path) && rm(output_path; force = true)

    output = NCDataset(output_path, "c")
    try
        for dimension in dimnames(template)
            dimension_length = if dimension == time_dim
                total_time
            elseif haskey(ranges, dimension)
                length(ranges[dimension])
            else
                NCDatasets.dim(template, dimension)
            end
            defDim(output, dimension, dimension_length)
        end

        for name in keys(template)
            variable = NCDatasets.variable(template, name)
            copy_auxiliary_variable(name, variable, time_dim, parameters) || continue
            all(dimension -> dimension in dimnames(output), dimnames(variable)) || continue
            define_subset_variable(output, template, name, ranges; deflatelevel = 0)
        end

        for name in parameters
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

function apply_spatial_mask(data, variable, mask, spatial_dimensions)
    spatial_indices = [findfirst(==(dimension), dimnames(variable)) for dimension in spatial_dimensions]
    any(isnothing, spatial_indices) && return data

    masked_data = copy(data)
    fill_value = masked_fill_value(masked_data)

    for mask_index in CartesianIndices(mask)
        mask[mask_index] && continue
        data_indices = ntuple(ndims(masked_data)) do i
            mask_dimension = findfirst(==(i), spatial_indices)
            mask_dimension === nothing ? (:) : mask_index[mask_dimension]
        end
        masked_data[data_indices...] .= fill_value
    end

    return masked_data
end

function write_parameter_chunk!(output, source, name, ranges, mask, spatial_dimensions, time_start)
    variable = source[name]
    input_indices = variable_indices(variable, ranges)
    data = concrete_float_data(variable[input_indices...])
    data = apply_spatial_mask(data, variable, mask, spatial_dimensions)

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

function write_time_dependent_coordinates!(output, source, ranges, time_start; parameters)
    time_dim = time_dimension(source)

    for name in keys(source)
        name in parameters && continue
        variable = NCDatasets.variable(source, name)
        time_index = findfirst(==(time_dim), dimnames(variable))
        time_index === nothing && continue
        haskey(output, name) || continue

        input_indices = variable_indices(variable, ranges)
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

function process_month(year, month, config::FjordConfig; files = nothing)
    forcing_config = config.forcing_config
    parameters = forcing_config.parameters
    latitude_range = config.grid_config.latitude
    longitude_range = config.grid_config.longitude

    files = isnothing(files) ? list_opendap_files(catalog_url = forcing_config.catalog_url) : files
    month_string = ".$(year)$(lpad(month, 2, '0'))"
    year_month = month_string[2:end]
    output_path = joinpath(norkyst_directory(forcing_config), norkyst_monthly_filename(year, month))

    if isfile(output_path)
        println("Skipping $year_month (already exists)")
        return output_path
    end

    println("Processing $year_month...")
    monthly_files = sort([file for file in files if occursin(month_string, file)])

    if isempty(monthly_files)
        println("  No files found for $year_month, skipping.")
        return nothing
    end

    urls = [joinpath(forcing_config.opendap_url, file) for file in monthly_files]
    println("  Opening $(length(urls)) datasets...")

    datasets = NCDataset[]
    try
        for url in urls
            push!(datasets, NCDataset(url))
            println("    Opened: $url")
        end

        ranges, mask, spatial_dimensions = subset_ranges(first(datasets); latitude_range, longitude_range)
        total_time = sum(time_length, datasets)

        println("  Writing output to: $output_path")
        output = define_output_file(output_path, first(datasets), ranges, total_time; parameters)
        try
            time_start = 1
            for ds in datasets
                write_time_dependent_coordinates!(output, ds, ranges, time_start; parameters)
                for name in parameters
                    write_parameter_chunk!(output, ds, name, ranges, mask, spatial_dimensions, time_start)
                end
                time_start += time_length(ds)
            end
        finally
            close(output)
        end
    finally
        foreach(close, datasets)
    end

    println("Finished $year_month\n")
    return output_path
end

function main()
    args = parse_args()
    config = args.config
    forcing_config = config.forcing_config
    output_directory = norkyst_directory(forcing_config)
    mkpath(output_directory)

    println("Processing years: $(join(forcing_config.years, ", "))")
    println("Config: $(args.config_path)")
    println("Output directory: $output_directory\n")

    files = list_opendap_files(catalog_url = forcing_config.catalog_url)
    for year in forcing_config.years
        for month in 1:12
            process_month(year, month, config; files)
        end
    end

    println("All done.")
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
