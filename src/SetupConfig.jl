module SetupConfig

using TOML

const DEFAULT_SETUP_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "configs", "drammensfjorden.toml"))

expand_user(path::AbstractString) = path == "~" ? homedir() : startswith(path, "~/") ? joinpath(homedir(), path[3:end]) : String(path)

function require_key(table, key, context)
    haskey(table, key) || error("Missing required config key `$key` in $context.")
    return table[key]
end

function tuple_value(table, key, context, ::Val{N}, T) where {N}
    values = require_key(table, key, context)
    values isa AbstractVector || error("Config key `$key` in $context must be an array.")
    length(values) == N || error("Config key `$key` in $context must contain $N values.")
    return Tuple(convert(T, value) for value in values)
end

function vector_value(table, key, context, T)
    values = require_key(table, key, context)
    values isa AbstractVector || error("Config key `$key` in $context must be an array.")
    return [convert(T, value) for value in values]
end

function path_from_config(path::AbstractString, data_root::AbstractString)
    expanded = expand_user(path)
    return isabspath(expanded) ? expanded : joinpath(data_root, expanded)
end

function load_setup_config(path::AbstractString = DEFAULT_SETUP_CONFIG_PATH)
    config_path = normpath(expand_user(path))
    raw = TOML.parsefile(config_path)

    name = String(require_key(raw, "name", "top-level config"))
    data_root = expand_user(String(require_key(raw, "data_root", "top-level config")))

    grid = require_key(raw, "grid", "top-level config")
    bathymetry = require_key(raw, "bathymetry", "top-level config")
    forcing = require_key(raw, "forcing", "top-level config")
    norkyst = require_key(forcing, "norkyst", "forcing config")

    grid_config = (
        size = tuple_value(grid, "size", "grid config", Val(3), Int),
        halo = tuple_value(grid, "halo", "grid config", Val(3), Int),
        longitude = tuple_value(grid, "longitude", "grid config", Val(2), Float64),
        latitude = tuple_value(grid, "latitude", "grid config", Val(2), Float64),
        z_faces = vector_value(grid, "z_faces", "grid config", Float64),
    )

    setup_data_dir = joinpath(data_root, name)
    bathymetry_name = "bathymetry_$name"
    bathymetry_config = (
        name = bathymetry_name,
        output_dir = setup_data_dir,
        output_path = joinpath(setup_data_dir, "$bathymetry_name.nc"),
        plot_path = joinpath(setup_data_dir, "$bathymetry_name.png"),
        geodatabase_path = path_from_config(String(require_key(bathymetry, "geodatabase_path", "bathymetry config")), data_root),
        raw_resolution_factor = convert(Int, require_key(bathymetry, "raw_resolution_factor", "bathymetry config")),
        padding_cells = convert(Int, require_key(bathymetry, "padding_cells", "bathymetry config")),
        include_contours = convert(Bool, require_key(bathymetry, "include_contours", "bathymetry config")),
        interpolation_passes = convert(Int, require_key(bathymetry, "interpolation_passes", "bathymetry config")),
        major_basins = convert(Int, require_key(bathymetry, "major_basins", "bathymetry config")),
        cache = convert(Bool, require_key(bathymetry, "cache", "bathymetry config")),
    )

    norkyst_config = (
        name = "forcing_$name",
        catalog_url = String(require_key(norkyst, "catalog_url", "NorKyst config")),
        opendap_url = String(require_key(norkyst, "opendap_url", "NorKyst config")),
        parameters = Tuple(vector_value(norkyst, "parameters", "NorKyst config", String)),
        years = Tuple(vector_value(norkyst, "years", "NorKyst config", Int)),
        output_dir = setup_data_dir,
    )

    return (
        path = config_path,
        name = name,
        data_root = data_root,
        grid = grid_config,
        bathymetry = bathymetry_config,
        norkyst = norkyst_config,
    )
end

end # module
