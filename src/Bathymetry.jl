module Bathymetry

export DybdedataConfig, prepare_geonorge_bathymetry, bathymetry_path, plot_path, geodatabase_path

using Scratch
using ArchGDAL
using Downloads
using NCDatasets
using p7zip_jll
using NumericalEarth
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain, znodes

using NumericalEarth.DataWrangling: AbstractStaticBathymetry, Metadatum, metadata_path

using ..Configs: AbstractBathymetryConfig

import NumericalEarth.DataWrangling:
    dataset_variable_name,
    default_download_directory,
    download_dataset,
    latitude_interfaces,
    longitude_interfaces,
    metadata_filename,
    reversed_vertical_axis

const GEONORGE_LAND_LAYERS = ("landareal", "skjer")
# Geonorge serves complete-dataset FileGDB archives as static zip files. The archive
# basename matches the `.gdb` directory basename with a `.zip` extension, so the download
# URL is derivable from the configured `geodatabase_file`.
const GEONORGE_FGDB_BASE_URL = "https://nedlasting.geonorge.no/geonorge/Basisdata/Dybdedata/FGDB"
const GEONORGE_DYBDEDATA_GDB = "Basisdata_0000_Norge_25833_Dybdedata_FGDB.gdb"
# NumericalEarth bathymetry regridding constructs a native grid with halo = (10, 10, 1).
# Keep the generated raw dataset comfortably larger than that minimum.
const MIN_NATIVE_BATHYMETRY_SIZE = 24
# Matches the fixed loop count and neighbor threshold used in the Oslofjord notebook's
# post-regrid gap-filling pass.
const BATHYMETRY_GAP_FILL_PASSES = 10
const ISOLATED_SEA_CELL_LAND_SIDES = 3

download_bathymetry_cache::String = ""

function __init__()
    global download_bathymetry_cache = @get_scratch!("Bathymetry")
end

# --- Bathymetry config ---

"""
    DybdedataConfig

Configuration for `prepare_geonorge_bathymetry`, driven by the Geonorge Sjøkart Dybdedata
FileGDB dataset.

`data_root`, `output_file` and `plot_file` are required: a setup names its own output files.
The FileGDB name is a property of the Geonorge dataset, so it defaults to
`GEONORGE_DYBDEDATA_GDB` and is downloaded into `data_root` on first use.

`output_file`, `plot_file` and `geodatabase_file` are names relative to `data_root`, resolved
by `bathymetry_path`, `plot_path` and `geodatabase_path`. Setting one to an absolute path
overrides `data_root` for that file, which is how a single shared FileGDB copy is reused
across fjords.

# User-facing fields
- `data_root`: Directory holding this setup's bathymetry files. Required.
- `output_file`: Name of the processed FjordSim bathymetry NetCDF. Required.
- `plot_file`: Name of the diagnostic bathymetry plot. Required.
- `geodatabase_file`: Name of the Geonorge Sjøkart FileGDB database.
- `raw_directory`: Scratch directory holding the intermediate regional raw dataset.
- `raw_resolution_factor`: Native raw-grid refinement relative to the target grid.
- `padding_cells`: Number of target-grid cell widths added around the region.
- `include_contours`: Sample `dybdekurve` contour lines in addition to depth points.
- `contour_stride`: Sample every n-th contour vertex.
- `interpolation_passes`: Passed to `NumericalEarth.regrid_bathymetry`.
- `major_basins`: Passed to `NumericalEarth.regrid_bathymetry`.
- `geonorge_cache`: Reuse the cached regional raw NetCDF when available.
- `regrid_cache`: Use NumericalEarth's on-disk bathymetry cache.

# Fields derived from the target grid by `native_region!`
- `raw_file`: Path of the intermediate regional raw NetCDF.
- `longitude`: Padded native longitude bounds in degrees.
- `latitude`: Padded native latitude bounds in degrees.
- `size`: Native `(Nx, Ny, Nz)` size consumed by `NumericalEarth.regrid_bathymetry`.
- `filter_bounds`: `(xmin, ymin, xmax, ymax)` spatial filter in EPSG:25833 meters.
"""
Base.@kwdef mutable struct DybdedataConfig <: AbstractBathymetryConfig
    data_root::String
    output_file::String
    plot_file::String
    geodatabase_file::String = GEONORGE_DYBDEDATA_GDB
    raw_directory::String = download_bathymetry_cache
    raw_resolution_factor::Int = 4
    padding_cells::Int = 0
    include_contours::Bool = true
    contour_stride::Int = 1
    interpolation_passes::Int = 1
    major_basins::Int = 1
    geonorge_cache::Bool = true
    regrid_cache::Bool = true
    raw_file::String = ""
    longitude::NTuple{2,Float64} = (0.0, 0.0)
    latitude::NTuple{2,Float64} = (0.0, 0.0)
    size::NTuple{3,Int} = (0, 0, 1)
    filter_bounds::NTuple{4,Float64} = (0.0, 0.0, 0.0, 0.0)
end

"""
    bathymetry_path(config)
    plot_path(config)

Resolve `config.output_file` and `config.plot_file` against `config.data_root`. Defined for
every `AbstractBathymetryConfig`, so a new bathymetry source inherits path resolution. A
field holding an absolute path is returned unchanged, relocating that file outside
`data_root`.
"""
bathymetry_path(config::AbstractBathymetryConfig) = joinpath(config.data_root, config.output_file)
plot_path(config::AbstractBathymetryConfig) = joinpath(config.data_root, config.plot_file)

"""
    geodatabase_path(config::DybdedataConfig)

Resolve `config.geodatabase_file` against `config.data_root`. Point it at an absolute path to
share one FileGDB copy across setups.
"""
geodatabase_path(config::DybdedataConfig) = joinpath(config.data_root, config.geodatabase_file)

"""
    DepthSamples

Accumulator for a sampled bathymetry point cloud: `xs` and `ys` in the coordinate reference
system currently being read, and one bottom height per coordinate pair.
"""
struct DepthSamples
    xs::Vector{Float64}
    ys::Vector{Float64}
    bottom_heights::Vector{Float64}
end

DepthSamples() = DepthSamples(Float64[], Float64[], Float64[])

Base.length(samples::DepthSamples) = length(samples.xs)

"""
    GeonorgeBathymetry

Static NumericalEarth-compatible metadata wrapper for a regional raw bathymetry
NetCDF derived from Geonorge Sjøkart bathymetry data.

# Fields
- `metadata_filename`: Filename of the generated raw NetCDF.
- `default_download_directory`: Directory containing the raw NetCDF.
- `longitude_interfaces`: Native longitude bounds of the raw dataset.
- `latitude_interfaces`: Native latitude bounds of the raw dataset.
- `size`: Native `(Nx, Ny, Nz)` size consumed by `NumericalEarth.regrid_bathymetry`.
"""
struct GeonorgeBathymetry <: AbstractStaticBathymetry
    metadata_filename::String
    default_download_directory::String
    longitude_interfaces::NTuple{2,Float64}
    latitude_interfaces::NTuple{2,Float64}
    size::NTuple{3,Int}
end

GeonorgeBathymetry(config::DybdedataConfig) = GeonorgeBathymetry(
    basename(config.raw_file),
    dirname(config.raw_file),
    config.longitude,
    config.latitude,
    config.size,
)

const GeonorgeBathymetryMetadatum = Metadatum{<:GeonorgeBathymetry}

default_download_directory(dataset::GeonorgeBathymetry) = dataset.default_download_directory
metadata_filename(dataset::GeonorgeBathymetry, args...) = dataset.metadata_filename
longitude_interfaces(dataset::GeonorgeBathymetry) = dataset.longitude_interfaces
latitude_interfaces(dataset::GeonorgeBathymetry) = dataset.latitude_interfaces
reversed_vertical_axis(::GeonorgeBathymetry) = false
Base.size(dataset::GeonorgeBathymetry) = dataset.size

dataset_variable_name(::GeonorgeBathymetryMetadatum) = "z"

"""
    download_dataset(metadata::GeonorgeBathymetryMetadatum)

Return the generated regional raw bathymetry file.
This dataset is materialized locally by `prepare_geonorge_bathymetry`.
"""
function download_dataset(metadata::GeonorgeBathymetryMetadatum)
    filepath = metadata_path(metadata)
    isfile(filepath) || error("Raw bathymetry file $filepath does not exist. Run prepare_geonorge_bathymetry first.")
    return filepath
end

"""
    geonorge_geodatabase_url(geodatabase_file)

Derive the Geonorge complete-dataset download URL for a local FileGDB path. The remote
zip archive shares the `.gdb` directory basename with a `.zip` extension.
"""
function geonorge_geodatabase_url(geodatabase_file)
    stem = replace(basename(geodatabase_file), r"\.gdb$" => "")
    return "$(GEONORGE_FGDB_BASE_URL)/$(stem).zip"
end

"""
    ensure_geodatabase(config::DybdedataConfig)

Ensure a Geonorge FileGDB directory exists at `geodatabase_path(config)`, downloading and
extracting it from Geonorge's public file server when absent. Returns the FileGDB path.
"""
function ensure_geodatabase(config::DybdedataConfig)
    geodatabase_file = geodatabase_path(config)
    isdir(geodatabase_file) && return geodatabase_file

    url = geonorge_geodatabase_url(geodatabase_file)
    destination_directory = dirname(geodatabase_file)
    mkpath(destination_directory)

    stem = replace(basename(geodatabase_file), r"\.gdb$" => "")
    zip_path = joinpath(destination_directory, "$(stem).zip")
    staging_directory = joinpath(destination_directory, "$(stem)_staging")

    @info "Local Geonorge geodatabase not found at $geodatabase_file"
    @info "Downloading Geonorge geodatabase from $url (this is a large file, ~2.3 GB)"

    try
        Downloads.download(url, zip_path)

        isdir(staging_directory) && rm(staging_directory; recursive = true, force = true)
        mkpath(staging_directory)

        @info "Extracting Geonorge geodatabase to $geodatabase_file"
        run(`$(p7zip()) x $zip_path -o$staging_directory -y`)

        extracted = find_geodatabase_directory(staging_directory)
        isnothing(extracted) && error("No .gdb directory found in Geonorge archive $url")

        isdir(geodatabase_file) && rm(geodatabase_file; recursive = true, force = true)
        mv(extracted, geodatabase_file)
    finally
        rm(zip_path; force = true)
        rm(staging_directory; recursive = true, force = true)
    end

    @info "Finished preparing Geonorge geodatabase at $geodatabase_file"
    return geodatabase_file
end

function find_geodatabase_directory(root)
    for (directory, subdirectories, _) in walkdir(root)
        for subdirectory in subdirectories
            endswith(subdirectory, ".gdb") && return joinpath(directory, subdirectory)
        end
    end
    return nothing
end

"""
    native_region!(config::DybdedataConfig, target_grid)

Resolve the native raw-dataset region implied by `config` and `target_grid`, storing the
padded longitude/latitude bounds, the native size, the EPSG:25833 spatial filter bounds and
the raw NetCDF path on `config`. Returns `config`.
"""
function native_region!(config::DybdedataConfig, target_grid)
    config.raw_resolution_factor >= 1 || throw(ArgumentError("raw_resolution_factor must be >= 1"))
    config.padding_cells >= 0 || throw(ArgumentError("padding_cells must be >= 0"))
    config.contour_stride >= 1 || throw(ArgumentError("contour_stride must be >= 1"))

    Nx, Ny, _ = size(target_grid)
    config.longitude = expand_domain(x_domain(target_grid), Nx, config.padding_cells)
    config.latitude = expand_domain(y_domain(target_grid), Ny, config.padding_cells)
    config.size = (
        max(config.raw_resolution_factor * (Nx + 2 * config.padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        max(config.raw_resolution_factor * (Ny + 2 * config.padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        1,
    )
    config.filter_bounds = transformed_filter_bounds(config.longitude, config.latitude)
    config.raw_file = joinpath(config.raw_directory, geonorge_raw_filename(config))

    return config
end

"""
    prepare_geonorge_bathymetry(target_grid, config::DybdedataConfig; regrid_kw...)

Read the Geonorge Sjøkart FileGDB bathymetry dataset described by `config`, build a regional
NumericalEarth-style raw bathymetry dataset in scratch storage, regrid it onto `target_grid`
with `NumericalEarth.regrid_bathymetry`, and write a processed NetCDF file compatible with
`FjordSim.Grids.ImmersedBoundaryGrid`.

`config` is mutated by `native_region!` to record the native region derived from
`target_grid`. `regrid_kw...` is forwarded to `NumericalEarth.regrid_bathymetry`.

# Returns
A named tuple with `dataset`, `raw_file`, `output_file`, and `bottom_height`.
"""
function prepare_geonorge_bathymetry(target_grid, config::DybdedataConfig; regrid_kw...)
    native_region!(config, target_grid)
    ensure_geodatabase(config)

    @info "Preparing Geonorge bathymetry"
    dataset = geonorge_dataset(config)

    metadata = Metadatum(:bottom_height; dataset)
    @info "Regridding Geonorge bathymetry onto target grid"
    bottom_height = NumericalEarth.regrid_bathymetry(
        target_grid,
        metadata;
        cache = config.regrid_cache,
        interpolation_passes = config.interpolation_passes,
        major_basins = config.major_basins,
        regrid_kw...,
    )
    output_file = bathymetry_path(config)
    @info "Smoothing small-scale bathymetry gaps"
    smooth_bathymetry_gaps!(bottom_height)
    @info "Writing processed bathymetry file to $output_file"
    write_bathymetry_file(output_file, target_grid, bottom_height)
    @info "Finished preparing Geonorge bathymetry"

    return (; dataset, raw_file = metadata_path(metadata), output_file, bottom_height)
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

        lon = defVar(ds, "lon", Float64, ("lon",))
        lat = defVar(ds, "lat", Float64, ("lat",))
        zf = defVar(ds, "z_faces", Float64, ("zf",))
        hvar = defVar(ds, "h", Float32, ("lon", "lat"))

        lon[:] = longitude
        lat[:] = latitude
        zf[:] = z_faces
        hvar[:, :] = h
    finally
        close(ds)
    end

    return filepath
end

"""
    smooth_bathymetry_gaps!(bottom_height)

Remove small-scale diagonal checkerboard artifacts and single-cell sea/land noise left
by regridding, following the fixed cleanup pass used for the Oslofjord ROMS-based
bathymetry: one diagonal-pair fill, then `BATHYMETRY_GAP_FILL_PASSES` rounds of isolated
sea/land cell cleanup. Mutates `bottom_height` in place.
"""
function smooth_bathymetry_gaps!(bottom_height)
    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    h = Array(interior(cpu_bottom_height, :, :, 1))

    h = fill_secondary_diagonal_pairs(fill_diagonal_pairs(h))
    for _ = 1:BATHYMETRY_GAP_FILL_PASSES
        h = fill_isolated_land_cells(remove_isolated_sea_cells(h))
    end

    set!(bottom_height, h)
    return bottom_height
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

function geonorge_dataset(config::DybdedataConfig)
    isdir(config.raw_directory) || mkpath(config.raw_directory)

    if !config.geonorge_cache || !isfile(config.raw_file)
        @info "Building raw Geonorge bathymetry at $(config.raw_file)"
        write_native_bathymetry(config)
    else
        @info "Using cached raw Geonorge bathymetry at $(config.raw_file)"
    end

    return GeonorgeBathymetry(config)
end

function write_native_bathymetry(config::DybdedataConfig)
    Nx, Ny, _ = config.size
    @info "Writing native bathymetry grid with size ($Nx, $Ny)"
    longitude_centers = center_coordinates(config.longitude, Nx)
    latitude_centers = center_coordinates(config.latitude, Ny)
    z_data = build_native_bathymetry_data(config)

    isfile(config.raw_file) && rm(config.raw_file; force = true)

    ds = NCDataset(config.raw_file, "c")
    try
        defDim(ds, "lon", Nx)
        defDim(ds, "lat", Ny)

        lon = defVar(ds, "lon", Float64, ("lon",))
        lat = defVar(ds, "lat", Float64, ("lat",))
        z = defVar(ds, "z", Float32, ("lon", "lat"))

        lon[:] = longitude_centers
        lat[:] = latitude_centers
        z[:, :] = z_data
    finally
        close(ds)
    end

    @info "Finished writing native bathymetry file to $(config.raw_file)"
    return config.raw_file
end

function build_native_bathymetry_data(config::DybdedataConfig)
    return ArchGDAL.importEPSG(25833; order = :trad) do source_srs
        ArchGDAL.importEPSG(4326; order = :trad) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                bathymetry = create_point_dataset(config, transform, target_srs) do point_dataset
                    @info "Gridding sampled bathymetry depths"
                    grid_point_dataset(point_dataset, config)
                end

                land_mask = create_land_dataset(config, transform, target_srs) do land_dataset
                    @info "Rasterizing land features"
                    rasterize_land_dataset(land_dataset, config)
                end

                bathymetry[land_mask] .= 0.0f0
                bathymetry
            end
        end
    end
end

function create_point_dataset(f::Function, config::DybdedataConfig, transform, target_srs)
    ArchGDAL.create(ArchGDAL.getdriver("Memory")) do point_dataset
        ArchGDAL.createlayer(
            name = "bathymetry_points",
            dataset = point_dataset,
            geom = ArchGDAL.wkbPoint,
            spatialref = target_srs,
        ) do point_layer
            ArchGDAL.addfielddefn!(point_layer, "z", ArchGDAL.OFTReal)
            point_count = sample_bathymetry_points!(point_layer, config, transform)
            point_count > 0 || error("No Geonorge bathymetry features intersect the requested region.")
            return f(point_dataset)
        end
    end
end

function create_land_dataset(f::Function, config::DybdedataConfig, transform, target_srs)
    ArchGDAL.create(ArchGDAL.getdriver("Memory")) do land_dataset
        ArchGDAL.createlayer(
            name = "land_features",
            dataset = land_dataset,
            geom = ArchGDAL.wkbUnknown,
            spatialref = target_srs,
        ) do land_layer
            sample_land_features!(land_layer, config, transform)
            return f(land_dataset)
        end
    end
end

function grid_point_dataset(point_dataset, config::DybdedataConfig)
    Nx, Ny, _ = config.size
    options = [
        "-of",
        "MEM",
        "-a",
        "invdistnn:power=2:smoothing=0.2:max_points=16:min_points=1:nodata=0",
        "-zfield",
        "z",
        "-txe",
        string(config.longitude[1]),
        string(config.longitude[2]),
        "-tye",
        string(config.latitude[1]),
        string(config.latitude[2]),
        "-outsize",
        string(Nx),
        string(Ny),
        "-a_srs",
        "EPSG:4326",
        "-l",
        "bathymetry_points",
    ]

    ArchGDAL.gdalgrid(point_dataset, options) do raster_dataset
        band = ArchGDAL.getband(raster_dataset, 1)
        raster_to_bottom_height(band, Nx, Ny)
    end
end

function rasterize_land_dataset(land_dataset, config::DybdedataConfig)
    Nx, Ny, _ = config.size
    options = [
        "-of",
        "MEM",
        "-burn",
        "1",
        "-ot",
        "Byte",
        "-a_srs",
        "EPSG:4326",
        "-te",
        string(config.longitude[1]),
        string(config.latitude[1]),
        string(config.longitude[2]),
        string(config.latitude[2]),
        "-outsize",
        string(Nx),
        string(Ny),
        "-l",
        "land_features",
    ]

    ArchGDAL.gdalrasterize(land_dataset, options) do raster_dataset
        band = ArchGDAL.getband(raster_dataset, 1)
        raster_to_mask(band, Nx, Ny)
    end
end

function sample_bathymetry_points!(point_layer, config::DybdedataConfig, transform)
    samples = DepthSamples()

    ArchGDAL.read(geodatabase_path(config)) do dataset
        collect_depth_layer_coordinates!(samples, dataset, "dybdepunkt", config; geometry = :point)
        config.include_contours &&
            collect_depth_layer_coordinates!(samples, dataset, "dybdekurve", config; geometry = :line)
    end

    point_count = length(samples)
    point_count > 0 || return 0

    # Transform the whole point cloud in a single GDAL call instead of one call per
    # point/vertex, which otherwise dominates wall-clock time for dense sounding data.
    ArchGDAL.transform!(samples.xs, samples.ys, zeros(Float64, point_count), transform)

    z_field_index = ArchGDAL.findfieldindex(point_layer, "z", false)
    for i = 1:point_count
        ArchGDAL.createfeature(point_layer) do feature
            ArchGDAL.setfield!(feature, z_field_index, samples.bottom_heights[i])
            ArchGDAL.setgeom!(feature, 0, ArchGDAL.createpoint(samples.xs[i], samples.ys[i]))
            return nothing
        end
    end

    return point_count
end

function sample_land_features!(land_layer, config::DybdedataConfig, transform)
    xmin, ymin, xmax, ymax = config.filter_bounds
    @info "Sampling Geonorge land features"

    ArchGDAL.read(geodatabase_path(config)) do dataset
        for layer_name in GEONORGE_LAND_LAYERS
            layer = find_layer(dataset, layer_name)
            isnothing(layer) && continue

            ArchGDAL.setspatialfilter!(layer, xmin, ymin, xmax, ymax)
            feature_count = ArchGDAL.nfeature(layer, true)
            @info "Sampling $layer_name: $feature_count land features"

            for source_feature in layer
                geometry = ArchGDAL.clone(ArchGDAL.getgeom(source_feature))
                ArchGDAL.transform!(geometry, transform)

                ArchGDAL.createfeature(land_layer) do target_feature
                    ArchGDAL.setgeom!(target_feature, 0, geometry)
                    return nothing
                end
            end
        end
    end

    @info "Finished sampling Geonorge land features"
    return nothing
end

function collect_depth_layer_coordinates!(
    samples::DepthSamples,
    dataset,
    layer_name,
    config::DybdedataConfig;
    geometry,
)
    layer = find_layer(dataset, layer_name)
    isnothing(layer) && return nothing

    xmin, ymin, xmax, ymax = config.filter_bounds
    ArchGDAL.setspatialfilter!(layer, xmin, ymin, xmax, ymax)
    depth_index = ArchGDAL.findfieldindex(layer, "dybde", false)

    for feature in layer
        depth = ArchGDAL.getfield(feature, depth_index)
        ismissing(depth) && continue

        bottom_height = -abs(Float64(depth))
        geometry == :point && collect_point_coordinates!(samples, ArchGDAL.getgeom(feature), bottom_height)
        geometry == :line && collect_linestring_coordinates!(samples, ArchGDAL.getgeom(feature), bottom_height, config)
    end

    return nothing
end

function collect_point_coordinates!(samples::DepthSamples, point_geometry, bottom_height)
    x, y, _ = ArchGDAL.getpoint(point_geometry, 0)
    push!(samples.xs, x)
    push!(samples.ys, y)
    push!(samples.bottom_heights, bottom_height)
    return nothing
end

function collect_linestring_coordinates!(samples::DepthSamples, line, bottom_height, config::DybdedataConfig)
    if ArchGDAL.geomname(line) == "MULTILINESTRING"
        for geometry_index = 0:ArchGDAL.ngeom(line)-1
            collect_linestring_coordinates!(samples, ArchGDAL.getgeom(line, geometry_index), bottom_height, config)
        end
        return nothing
    end

    npoints = ArchGDAL.ngeom(line)

    for point_index in contour_point_indices(npoints, config.contour_stride)
        x, y, _ = ArchGDAL.getpoint(line, point_index)
        push!(samples.xs, x)
        push!(samples.ys, y)
        push!(samples.bottom_heights, bottom_height)
    end

    return nothing
end

function contour_point_indices(npoints, stride)
    last_index = npoints - 1
    indices = collect(0:stride:last_index)
    isempty(indices) || last(indices) == last_index || push!(indices, last_index)
    return indices
end

function raster_to_bottom_height(band, Nx, Ny)
    data = Array(ArchGDAL.read(band))
    if size(data) == (Ny, Nx)
        data = permutedims(data, (2, 1))
    elseif size(data) != (Nx, Ny)
        error("Unexpected raster shape $(size(data)); expected ($Nx, $Ny) or ($Ny, $Nx).")
    end
    return reverse(Float32.(data), dims = 2)
end

function raster_to_mask(band, Nx, Ny)
    data = Array(ArchGDAL.read(band))
    if size(data) == (Ny, Nx)
        data = permutedims(data, (2, 1))
    elseif size(data) != (Nx, Ny)
        error("Unexpected raster mask shape $(size(data)); expected ($Nx, $Ny) or ($Ny, $Nx).")
    end
    return reverse(data .> 0, dims = 2)
end

function transformed_filter_bounds(longitude, latitude)
    corners = (
        (longitude[1], latitude[1]),
        (longitude[1], latitude[2]),
        (longitude[2], latitude[1]),
        (longitude[2], latitude[2]),
    )

    xs = Float64[]
    ys = Float64[]

    ArchGDAL.importEPSG(4326; order = :trad) do source_srs
        ArchGDAL.importEPSG(25833; order = :trad) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                for (lon, lat) in corners
                    point = ArchGDAL.createpoint(lon, lat)
                    ArchGDAL.transform!(point, transform)
                    x, y, _ = ArchGDAL.getpoint(point, 0)
                    push!(xs, x)
                    push!(ys, y)
                end
            end
        end
    end

    return minimum(xs), minimum(ys), maximum(xs), maximum(ys)
end

function find_layer(dataset, layer_name)
    for index = 0:ArchGDAL.nlayer(dataset)-1
        layer = ArchGDAL.getlayer(dataset, index)
        ArchGDAL.getname(layer) == layer_name && return layer
    end
    return nothing
end

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

function geonorge_raw_filename(config::DybdedataConfig)
    Nx, Ny, _ = config.size
    lon1 = replace(string(round(config.longitude[1], digits = 3)), '.' => 'p')
    lon2 = replace(string(round(config.longitude[2], digits = 3)), '.' => 'p')
    lat1 = replace(string(round(config.latitude[1], digits = 3)), '.' => 'p')
    lat2 = replace(string(round(config.latitude[2], digits = 3)), '.' => 'p')
    contour_suffix = config.include_contours ? "_contour_stride_$(config.contour_stride)" : "_points_only"
    return "geonorge_sjokart_bathymetry_$(lon1)_$(lon2)_$(lat1)_$(lat2)_$(Nx)x$(Ny)$(contour_suffix).nc"
end

end  # module Bathymetry
