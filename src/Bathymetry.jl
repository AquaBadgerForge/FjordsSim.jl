module Bathymetry

export BathymetryConfig, prepare_geonorge_bathymetry

using Scratch
using ArchGDAL
using Downloads
using NCDatasets
using p7zip_jll
using NumericalEarth
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain

using NumericalEarth.DataWrangling: AbstractStaticBathymetry, Metadatum, metadata_path

import NumericalEarth.DataWrangling:
    dataset_variable_name,
    default_download_directory,
    download_dataset,
    latitude_interfaces,
    longitude_interfaces,
    metadata_filename,
    reversed_vertical_axis

# --- Bathymetry config ---

"""
    BathymetryConfig

Configuration for `prepare_geonorge_bathymetry`.

# Fields
- `output_path`: Destination for the processed FjordSim bathymetry NetCDF.
- `plot_path`: Destination for the diagnostic bathymetry plot.
- `geodatabase_path`: Path to the local Geonorge Sjøkart FileGDB database.
- `raw_resolution_factor`: Native raw-grid refinement relative to the target grid.
- `padding_cells`: Number of target-grid cell widths added around the region.
- `include_contours`: Sample `dybdekurve` contour lines in addition to depth points.
- `contour_stride`: Sample every n-th contour vertex.
- `interpolation_passes`: Passed to `NumericalEarth.regrid_bathymetry`.
- `major_basins`: Passed to `NumericalEarth.regrid_bathymetry`.
- `geonorge_cache`: Reuse the cached regional raw NetCDF when available.
- `regrid_cache`: Use NumericalEarth's on-disk bathymetry cache.
"""
Base.@kwdef struct BathymetryConfig
    output_path::String
    plot_path::String
    geodatabase_path::String
    raw_resolution_factor::Int
    padding_cells::Int
    include_contours::Bool
    contour_stride::Int
    interpolation_passes::Int
    major_basins::Int
    geonorge_cache::Bool
    regrid_cache::Bool
end

const GEONORGE_LAND_LAYERS = ("landareal", "skjer")
# Geonorge serves complete-dataset FileGDB archives as static zip files. The archive
# basename matches the `.gdb` directory basename with a `.zip` extension, so the download
# URL is derivable from the configured `geodatabase_path`.
const GEONORGE_FGDB_BASE_URL = "https://nedlasting.geonorge.no/geonorge/Basisdata/Dybdedata/FGDB"
# NumericalEarth bathymetry regridding constructs a native grid with halo = (10, 10, 1).
# Keep the generated raw dataset comfortably larger than that minimum.
const MIN_NATIVE_BATHYMETRY_SIZE = 24

download_bathymetry_cache::String = ""

function __init__()
    global download_bathymetry_cache = @get_scratch!("Bathymetry")
end

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
    geonorge_geodatabase_url(geodatabase_path)

Derive the Geonorge complete-dataset download URL for a local FileGDB path. The remote
zip archive shares the `.gdb` directory basename with a `.zip` extension.
"""
function geonorge_geodatabase_url(geodatabase_path)
    stem = replace(basename(geodatabase_path), r"\.gdb$" => "")
    return "$(GEONORGE_FGDB_BASE_URL)/$(stem).zip"
end

"""
    ensure_geodatabase(geodatabase_path)

Ensure a Geonorge FileGDB directory exists at `geodatabase_path`, downloading and
extracting it from Geonorge's public file server when absent. Returns `geodatabase_path`.
"""
function ensure_geodatabase(geodatabase_path)
    isdir(geodatabase_path) && return geodatabase_path

    url = geonorge_geodatabase_url(geodatabase_path)
    destination_directory = dirname(geodatabase_path)
    mkpath(destination_directory)

    stem = replace(basename(geodatabase_path), r"\.gdb$" => "")
    zip_path = joinpath(destination_directory, "$(stem).zip")
    staging_directory = joinpath(destination_directory, "$(stem)_staging")

    @info "Local Geonorge geodatabase not found at $geodatabase_path"
    @info "Downloading Geonorge geodatabase from $url (this is a large file, ~2.3 GB)"

    try
        Downloads.download(url, zip_path)

        isdir(staging_directory) && rm(staging_directory; recursive = true, force = true)
        mkpath(staging_directory)

        @info "Extracting Geonorge geodatabase to $geodatabase_path"
        run(`$(p7zip()) x $zip_path -o$staging_directory -y`)

        extracted = find_geodatabase_directory(staging_directory)
        isnothing(extracted) && error("No .gdb directory found in Geonorge archive $url")

        isdir(geodatabase_path) && rm(geodatabase_path; recursive = true, force = true)
        mv(extracted, geodatabase_path)
    finally
        rm(zip_path; force = true)
        rm(staging_directory; recursive = true, force = true)
    end

    @info "Finished preparing Geonorge geodatabase at $geodatabase_path"
    return geodatabase_path
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
    prepare_geonorge_bathymetry(target_grid; output_path, geodatabase_path, raw_dir=download_bathymetry_cache,
                                raw_resolution_factor=4, padding_cells=0,
                                include_contours=true, contour_stride=1,
                                geonorge_cache=true, regrid_cache=true, reporter, regrid_kw...)

Read the local Geonorge Sjøkart FileGDB bathymetry dataset, build a regional
NumericalEarth-style raw bathymetry dataset in scratch storage, regrid it onto
`target_grid` with `NumericalEarth.regrid_bathymetry`, and write a processed
NetCDF file compatible with `FjordSim.Grids.ImmersedBoundaryGrid`.

# Keyword arguments
- `output_path`: Destination for the processed FjordSim bathymetry NetCDF.
- `geodatabase_path`: Path to the local Geonorge Sjøkart FileGDB database.
- `raw_dir`: Scratch directory used for the intermediate regional raw dataset.
- `raw_resolution_factor`: Native raw-grid refinement relative to `target_grid`.
- `padding_cells`: Number of target-grid cell widths added around the requested region.
- `include_contours`: If `true` (default), sample both `dybdepunkt` and `dybdekurve`.
- `contour_stride`: Sample every `contour_stride`-th contour vertex when `include_contours=true`.
- `geonorge_cache`: If `true` (default), reuse the generated regional raw NetCDF.
- `regrid_cache`: If `true` (default), use NumericalEarth's on-disk bathymetry cache.
- `reporter`: Progress reporter; defaults to `LoggingProgressReporter()`.
- `regrid_kw...`: Forwarded directly to `NumericalEarth.regrid_bathymetry`.

# Returns
A named tuple with `dataset`, `raw_path`, `output_path`, and `bottom_height`.
"""
function prepare_geonorge_bathymetry(
    target_grid;
    output_path::String,
    geodatabase_path::String,
    raw_dir::String = download_bathymetry_cache,
    raw_resolution_factor::Int = 4,
    padding_cells::Int = 0,
    include_contours::Bool = true,
    contour_stride::Int = 1,
    geonorge_cache::Bool = true,
    regrid_cache::Bool = true,
    regrid_kw...,
)
    raw_resolution_factor >= 1 || throw(ArgumentError("raw_resolution_factor must be >= 1"))
    padding_cells >= 0 || throw(ArgumentError("padding_cells must be >= 0"))
    contour_stride >= 1 || throw(ArgumentError("contour_stride must be >= 1"))
    ensure_geodatabase(geodatabase_path)

    @info "Preparing Geonorge bathymetry"
    dataset = geonorge_dataset(
        target_grid;
        raw_dir,
        raw_resolution_factor,
        padding_cells,
        include_contours,
        contour_stride,
        geonorge_cache,
        geodatabase_path,
    )

    metadata = Metadatum(:bottom_height; dataset)
    @info "Regridding Geonorge bathymetry onto target grid"
    bottom_height = NumericalEarth.regrid_bathymetry(target_grid, metadata; cache = regrid_cache, regrid_kw...)
    @info "Writing processed bathymetry file to $output_path"
    write_bathymetry_file(output_path, target_grid, bottom_height)
    @info "Finished preparing Geonorge bathymetry"

    return (; dataset, raw_path = metadata_path(metadata), output_path, bottom_height)
end

function prepare_geonorge_bathymetry(target_grid, config::BathymetryConfig; kw...)
    return prepare_geonorge_bathymetry(
        target_grid;
        output_path = config.output_path,
        geodatabase_path = config.geodatabase_path,
        raw_resolution_factor = config.raw_resolution_factor,
        padding_cells = config.padding_cells,
        include_contours = config.include_contours,
        contour_stride = config.contour_stride,
        geonorge_cache = config.geonorge_cache,
        regrid_cache = config.regrid_cache,
        interpolation_passes = config.interpolation_passes,
        major_basins = config.major_basins,
        kw...,
    )
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

function geonorge_dataset(
    target_grid;
    raw_dir,
    raw_resolution_factor,
    padding_cells,
    include_contours,
    contour_stride,
    geonorge_cache,
    geodatabase_path,
)
    isdir(raw_dir) || mkpath(raw_dir)

    Nx, Ny, _ = size(target_grid)
    longitude = expand_domain(x_domain(target_grid), Nx, padding_cells)
    latitude = expand_domain(y_domain(target_grid), Ny, padding_cells)
    raw_size = (
        max(raw_resolution_factor * (Nx + 2 * padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        max(raw_resolution_factor * (Ny + 2 * padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        1,
    )

    raw_filename = geonorge_raw_filename(longitude, latitude, raw_size; include_contours, contour_stride)
    raw_path = joinpath(raw_dir, raw_filename)

    if !geonorge_cache || !isfile(raw_path)
        @info "Building raw Geonorge bathymetry at $raw_path"
        write_native_bathymetry(
            raw_path,
            geodatabase_path;
            longitude,
            latitude,
            size = raw_size,
            include_contours,
            contour_stride,
        )
    else
        @info "Using cached raw Geonorge bathymetry at $raw_path"
    end

    return GeonorgeBathymetry(raw_filename, raw_dir, longitude, latitude, raw_size)
end

function write_native_bathymetry(
    filepath,
    geodatabase_path;
    longitude,
    latitude,
    size,
    include_contours::Bool = true,
    contour_stride::Int = 1,
)
    Nx, Ny, _ = size
    @info "Writing native bathymetry grid with size ($Nx, $Ny)"
    longitude_centers = center_coordinates(longitude, Nx)
    latitude_centers = center_coordinates(latitude, Ny)
    z_data =
        build_native_bathymetry_data(geodatabase_path; longitude, latitude, Nx, Ny, include_contours, contour_stride)

    isfile(filepath) && rm(filepath; force = true)

    ds = NCDataset(filepath, "c")
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

    @info "Finished writing native bathymetry file to $filepath"
    return filepath
end

function build_native_bathymetry_data(
    geodatabase_path;
    longitude,
    latitude,
    Nx,
    Ny,
    include_contours::Bool = true,
    contour_stride::Int = 1,
)
    filter_bounds = transformed_filter_bounds(longitude, latitude)

    return ArchGDAL.importEPSG(25833; order = :trad) do source_srs
        ArchGDAL.importEPSG(4326; order = :trad) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                bathymetry = create_point_dataset(
                    geodatabase_path,
                    transform,
                    filter_bounds,
                    target_srs;
                    include_contours,
                    contour_stride,
                ) do point_dataset
                    @info "Gridding sampled bathymetry depths"
                    grid_point_dataset(point_dataset; longitude, latitude, Nx, Ny)
                end

                land_mask = create_land_dataset(geodatabase_path, transform, filter_bounds, target_srs) do land_dataset
                    @info "Rasterizing land features"
                    rasterize_land_dataset(land_dataset; longitude, latitude, Nx, Ny)
                end

                bathymetry[land_mask] .= 0.0f0
                bathymetry
            end
        end
    end
end

function create_point_dataset(
    f::Function,
    geodatabase_path,
    transform,
    filter_bounds,
    target_srs;
    include_contours::Bool = true,
    contour_stride::Int = 1,
)
    ArchGDAL.create(ArchGDAL.getdriver("Memory")) do point_dataset
        ArchGDAL.createlayer(
            name = "bathymetry_points",
            dataset = point_dataset,
            geom = ArchGDAL.wkbPoint,
            spatialref = target_srs,
        ) do point_layer
            ArchGDAL.addfielddefn!(point_layer, "z", ArchGDAL.OFTReal)
            point_count = sample_bathymetry_points!(
                point_layer,
                geodatabase_path,
                transform,
                filter_bounds;
                include_contours,
                contour_stride,
            )
            point_count > 0 || error("No Geonorge bathymetry features intersect the requested region.")
            return f(point_dataset)
        end
    end
end

function create_land_dataset(f::Function, geodatabase_path, transform, filter_bounds, target_srs)
    ArchGDAL.create(ArchGDAL.getdriver("Memory")) do land_dataset
        ArchGDAL.createlayer(
            name = "land_features",
            dataset = land_dataset,
            geom = ArchGDAL.wkbUnknown,
            spatialref = target_srs,
        ) do land_layer
            sample_land_features!(land_layer, geodatabase_path, transform, filter_bounds)
            return f(land_dataset)
        end
    end
end

function grid_point_dataset(point_dataset; longitude, latitude, Nx, Ny)
    options = [
        "-of",
        "MEM",
        "-a",
        "invdistnn:power=2:smoothing=0.2:max_points=16:min_points=1:nodata=0",
        "-zfield",
        "z",
        "-txe",
        string(longitude[1]),
        string(longitude[2]),
        "-tye",
        string(latitude[1]),
        string(latitude[2]),
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

function rasterize_land_dataset(land_dataset; longitude, latitude, Nx, Ny)
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
        string(longitude[1]),
        string(latitude[1]),
        string(longitude[2]),
        string(latitude[2]),
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

function sample_bathymetry_points!(
    point_layer,
    geodatabase_path,
    transform,
    filter_bounds;
    include_contours::Bool = true,
    contour_stride::Int = 1,
)
    xmin, ymin, xmax, ymax = filter_bounds
    xs = Float64[]
    ys = Float64[]
    bottom_heights = Float64[]

    ArchGDAL.read(geodatabase_path) do dataset
        collect_depth_layer_coordinates!(
            xs,
            ys,
            bottom_heights,
            dataset,
            "dybdepunkt",
            xmin,
            ymin,
            xmax,
            ymax;
            geometry = :point,
        )
        include_contours && collect_depth_layer_coordinates!(
            xs,
            ys,
            bottom_heights,
            dataset,
            "dybdekurve",
            xmin,
            ymin,
            xmax,
            ymax;
            geometry = :line,
            contour_stride,
        )
    end

    point_count = length(xs)
    point_count > 0 || return 0

    # Transform the whole point cloud in a single GDAL call instead of one call per
    # point/vertex, which otherwise dominates wall-clock time for dense sounding data.
    ArchGDAL.transform!(xs, ys, zeros(Float64, point_count), transform)

    z_field_index = ArchGDAL.findfieldindex(point_layer, "z", false)
    for i = 1:point_count
        ArchGDAL.createfeature(point_layer) do feature
            ArchGDAL.setfield!(feature, z_field_index, bottom_heights[i])
            ArchGDAL.setgeom!(feature, 0, ArchGDAL.createpoint(xs[i], ys[i]))
            return nothing
        end
    end

    return point_count
end

function sample_land_features!(land_layer, geodatabase_path, transform, filter_bounds)
    xmin, ymin, xmax, ymax = filter_bounds
    @info "Sampling Geonorge land features"

    ArchGDAL.read(geodatabase_path) do dataset
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
    xs,
    ys,
    bottom_heights,
    dataset,
    layer_name,
    xmin,
    ymin,
    xmax,
    ymax;
    geometry,
    contour_stride::Int = 1,
)
    layer = find_layer(dataset, layer_name)
    isnothing(layer) && return nothing

    ArchGDAL.setspatialfilter!(layer, xmin, ymin, xmax, ymax)
    depth_index = ArchGDAL.findfieldindex(layer, "dybde", false)

    for feature in layer
        depth = ArchGDAL.getfield(feature, depth_index)
        ismissing(depth) && continue

        bottom_height = -abs(Float64(depth))
        geometry == :point &&
            collect_point_coordinates!(xs, ys, bottom_heights, ArchGDAL.getgeom(feature), bottom_height)
        geometry == :line && collect_linestring_coordinates!(
            xs,
            ys,
            bottom_heights,
            ArchGDAL.getgeom(feature),
            bottom_height;
            stride = contour_stride,
        )
    end

    return nothing
end

function collect_point_coordinates!(xs, ys, bottom_heights, point_geometry, bottom_height)
    x, y, _ = ArchGDAL.getpoint(point_geometry, 0)
    push!(xs, x)
    push!(ys, y)
    push!(bottom_heights, bottom_height)
    return nothing
end

function collect_linestring_coordinates!(xs, ys, bottom_heights, line, bottom_height; stride::Int = 1)
    if ArchGDAL.geomname(line) == "MULTILINESTRING"
        for geometry_index = 0:ArchGDAL.ngeom(line)-1
            collect_linestring_coordinates!(
                xs,
                ys,
                bottom_heights,
                ArchGDAL.getgeom(line, geometry_index),
                bottom_height;
                stride,
            )
        end
        return nothing
    end

    npoints = ArchGDAL.ngeom(line)

    for point_index in contour_point_indices(npoints, stride)
        x, y, _ = ArchGDAL.getpoint(line, point_index)
        push!(xs, x)
        push!(ys, y)
        push!(bottom_heights, bottom_height)
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

function vertical_faces(grid)
    z_faces_with_halo = getproperty(grid.z, Symbol("cᵃᵃᶠ"))
    first_index = grid.Hz + 1
    last_index = first_index + grid.Nz
    return collect(z_faces_with_halo[first_index:last_index])
end

function expand_domain(domain, N, padding_cells)
    delta = domain_step(domain, N)
    lower = domain[1] - padding_cells * delta
    upper = domain[2] + padding_cells * delta
    return (lower, upper)
end

function geonorge_raw_filename(longitude, latitude, size; include_contours::Bool = true, contour_stride::Int = 1)
    Nx, Ny, _ = size
    lon1 = replace(string(round(longitude[1], digits = 3)), '.' => 'p')
    lon2 = replace(string(round(longitude[2], digits = 3)), '.' => 'p')
    lat1 = replace(string(round(latitude[1], digits = 3)), '.' => 'p')
    lat2 = replace(string(round(latitude[2], digits = 3)), '.' => 'p')
    contour_suffix = include_contours ? "_contour_stride_$(contour_stride)" : "_points_only"
    return "geonorge_sjokart_bathymetry_$(lon1)_$(lon2)_$(lat1)_$(lat2)_$(Nx)x$(Ny)$(contour_suffix).nc"
end

end  # module Bathymetry
