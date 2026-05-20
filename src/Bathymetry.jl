module Bathymetry

export GeonorgeBathymetry, prepare_geonorge_bathymetry, write_bathymetry_file

using Scratch
using ArchGDAL
using NCDatasets
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

const GEONORGE_LAND_LAYERS = ("landareal", "skjer")
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
    prepare_geonorge_bathymetry(target_grid; output_path, geodatabase_path, raw_dir=download_bathymetry_cache,
                                raw_resolution_factor=4, padding_cells=2,
                                include_contours=true, cache=true, regrid_kw...)

Read the local Geonorge Sjøkart FileGDB bathymetry dataset, build a regional
NumericalEarth-style raw bathymetry dataset in scratch storage, regrid it onto
`target_grid` with `NumericalEarth.regrid_bathymetry`, and write a processed
NetCDF file compatible with `FjordSim.Grids.ImmersedBoundaryGrid`.

# Keyword arguments
- `output_path`: Destination for the processed FjordSim bathymetry NetCDF.
- `geodatabase_path`: Path to the local Geonorge Sjøkart FileGDB database.
- `raw_dir`: Scratch directory used for the intermediate regional raw dataset.
- `raw_resolution_factor`: Native raw-grid refinement relative to `target_grid`.
    This controls the resolution of the intermediate regional raw bathymetry NetCDF
    built from the local Geonorge FileGDB before NumericalEarth regrids it onto
    `target_grid`. Larger values produce a finer temporary source grid but do not
    change the size of the final FjordSim grid.
- `padding_cells`: Number of target-grid cell widths added around the requested region.
- `include_contours`: If `true` (default), sample both `dybdepunkt` and
    `dybdekurve`. Set to `false` to grid only depth points, which is usually much
    faster for dense local datasets.
- `cache`: If `true` (default), reuse both the generated regional raw NetCDF and
    NumericalEarth's on-disk bathymetry cache. Set to `false` to rebuild the raw
    NetCDF from the local FileGDB and force regridding.
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
    padding_cells::Int = 2,
    include_contours::Bool = true,
    cache::Bool = true,
    regrid_kw...,
)
    raw_resolution_factor >= 1 || throw(ArgumentError("raw_resolution_factor must be >= 1"))
    padding_cells >= 0 || throw(ArgumentError("padding_cells must be >= 0"))
    :cache in keys(regrid_kw) &&
        throw(ArgumentError("Pass `cache` directly to prepare_geonorge_bathymetry, not via `regrid_kw...`."))
    isdir(geodatabase_path) ||
        error("Local Geonorge bathymetry geodatabase not found at $geodatabase_path.")

    @info "Preparing Geonorge bathymetry"
    dataset = geonorge_dataset(target_grid; raw_dir, raw_resolution_factor, padding_cells, include_contours, cache, geodatabase_path)

    metadata = Metadatum(:bottom_height; dataset)
    @info "Regridding Geonorge bathymetry onto target grid"
    bottom_height = NumericalEarth.regrid_bathymetry(target_grid, metadata; cache, regrid_kw...)
    try
        validate_land_representation(
            Array(interior(on_architecture(CPU(), bottom_height), :, :, 1));
            context = "regridded",
        )
        @info "Validated regridded Geonorge bathymetry"
    catch err
        if cache
            @info "Cached Geonorge bathymetry lacks land cells; recomputing without NumericalEarth cache."
            bottom_height = NumericalEarth.regrid_bathymetry(target_grid, metadata; cache = false, regrid_kw...)
            validate_land_representation(
                Array(interior(on_architecture(CPU(), bottom_height), :, :, 1));
                context = "regridded",
            )
        else
            rethrow(err)
        end
    end
    @info "Writing processed bathymetry file to $output_path"
    write_bathymetry_file(output_path, target_grid, bottom_height)
    @info "Finished preparing Geonorge bathymetry"

    return (; dataset, raw_path = metadata_path(metadata), output_path, bottom_height)
end

"""
    write_bathymetry_file(filepath, target_grid, bottom_height)

Write a processed NetCDF bathymetry file compatible with
`FjordSim.Grids.ImmersedBoundaryGrid`.

# Arguments
- `filepath`: Output NetCDF path.
- `target_grid`: Oceananigans target grid.
- `bottom_height`: Bottom-height field already regridded onto `target_grid`.
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
    geonorge_dataset(target_grid; raw_dir, raw_resolution_factor, padding_cells, include_contours, cache, geodatabase_path)

Construct the regional raw bathymetry dataset wrapper used as input to
`NumericalEarth.regrid_bathymetry`.
"""
function geonorge_dataset(target_grid; raw_dir, raw_resolution_factor, padding_cells, include_contours, cache, geodatabase_path)
    isdir(raw_dir) || mkpath(raw_dir)

    Nx, Ny, _ = size(target_grid)
    longitude = expand_domain(x_domain(target_grid), Nx, padding_cells)
    latitude = expand_domain(y_domain(target_grid), Ny, padding_cells)
    raw_size = (
        max(raw_resolution_factor * (Nx + 2 * padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        max(raw_resolution_factor * (Ny + 2 * padding_cells), MIN_NATIVE_BATHYMETRY_SIZE),
        1,
    )

    raw_filename = geonorge_raw_filename(longitude, latitude, raw_size; include_contours)
    raw_path = joinpath(raw_dir, raw_filename)

    if !cache || !isfile(raw_path)
        @info "Building raw Geonorge bathymetry at $raw_path"
        write_native_bathymetry(raw_path, geodatabase_path; longitude, latitude, size = raw_size, include_contours)
    else
        @info "Using cached raw Geonorge bathymetry at $raw_path"
    end

    return GeonorgeBathymetry(raw_filename, raw_dir, longitude, latitude, raw_size)
end



"""
    write_native_bathymetry(filepath, geodatabase_path; longitude, latitude, size, include_contours=true)

Create the regional raw bathymetry NetCDF consumed by `GeonorgeBathymetry`.
The native raster is built from Sjøkart depth points and depth contours,
combined with land polygons and skerries so that land cells remain `h >= 0`.
"""
function write_native_bathymetry(filepath, geodatabase_path; longitude, latitude, size, include_contours::Bool = true)
    Nx, Ny, _ = size
    @info "Writing native bathymetry grid with size ($Nx, $Ny)"
    longitude_centers = center_coordinates(longitude, Nx)
    latitude_centers = center_coordinates(latitude, Ny)
    z_data = build_native_bathymetry_data(geodatabase_path; longitude, latitude, Nx, Ny, include_contours)
    validate_land_representation(z_data; context = "raw")
    @info "Validated raw Geonorge bathymetry"

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

"""
    build_native_bathymetry_data(geodatabase_path; longitude, latitude, Nx, Ny, include_contours=true)

Create the raw regional bathymetry array by sampling Sjøkart depth features,
gridding them onto a regular WGS84 longitude-latitude raster, and then burning
land features back to `0 m`.
"""
function build_native_bathymetry_data(geodatabase_path; longitude, latitude, Nx, Ny, include_contours::Bool = true)
    @info "Native bathymetry build: transforming bounds"
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
                ) do point_dataset
                    @info "Native bathymetry build: gridding sampled depths"
                    grid_point_dataset(point_dataset; longitude, latitude, Nx, Ny)
                end

                land_mask = create_land_dataset(geodatabase_path, transform, filter_bounds, target_srs) do land_dataset
                    @info "Native bathymetry build: rasterizing land"
                    rasterize_land_dataset(land_dataset; longitude, latitude, Nx, Ny)
                end

                bathymetry[land_mask] .= 0.0f0
                @info "Native bathymetry build complete"
                bathymetry
            end
        end
    end
end

"""
    validate_land_representation(z_data; context = "bathymetry")

Reject bathymetry rasters that contain only underwater values.

Oceananigans and NumericalEarth interpret `h >= 0` as land, so a raster whose
maximum value is still below zero cannot encode a usable land-sea boundary.
"""
function validate_land_representation(z_data; context = "bathymetry")
    maximum_height = maximum(z_data)
    maximum_height >= 0 && return nothing

    error(
        "The source Geonorge bathymetry produced no land cells for this region " *
        "(maximum $(context) value = $(maximum_height) m). Oceananigans requires land cells with h >= 0 " *
        "to distinguish land from sea.",
    )
end

"""
    create_point_dataset(geodatabase_path, transform, filter_bounds, target_srs; include_contours=true, f)

Create an in-memory point dataset in WGS84 containing sampled Sjøkart depth
points and contour vertices,
then invoke `f(point_dataset)`.
"""
function create_point_dataset(
    f::Function,
    geodatabase_path,
    transform,
    filter_bounds,
    target_srs;
    include_contours::Bool = true,
)
    @info "Creating native bathymetry point dataset"
    ArchGDAL.create(ArchGDAL.getdriver("Memory")) do point_dataset
        ArchGDAL.createlayer(
            name = "bathymetry_points",
            dataset = point_dataset,
            geom = ArchGDAL.wkbPoint,
            spatialref = target_srs,
        ) do point_layer
            ArchGDAL.addfielddefn!(point_layer, "z", ArchGDAL.OFTReal)
            point_count =
                sample_bathymetry_points!(point_layer, geodatabase_path, transform, filter_bounds; include_contours)
            point_count > 0 || error("No Geonorge bathymetry features intersect the requested region.")
            @info "Created native bathymetry point dataset with $point_count sampled points"
            return f(point_dataset)
        end
    end
end

"""
    create_land_dataset(geodatabase_path, transform, filter_bounds, target_srs, f)

Create an in-memory WGS84 vector dataset containing transformed land features,
then invoke `f(land_dataset)`.
"""
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

"""
    grid_point_dataset(point_dataset; longitude, latitude, Nx, Ny)

Grid the sampled point dataset onto a regular native bathymetry raster.
"""
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

"""
    rasterize_land_dataset(land_dataset; longitude, latitude, Nx, Ny)

Rasterize transformed land features onto the same native WGS84 grid used for
the raw bathymetry raster.
"""
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

"""
    sample_bathymetry_points!(point_layer, geodatabase_path, transform, filter_bounds; include_contours=true)

Sample depth points and contour vertices from the clipped Sjøkart geodatabase,
transform them to WGS84, and append them as point features with a `z` value.
Returns the number of appended points.
"""
function sample_bathymetry_points!(
    point_layer,
    geodatabase_path,
    transform,
    filter_bounds;
    include_contours::Bool = true,
)
    xmin, ymin, xmax, ymax = filter_bounds
    point_count = 0
    @info "Sampling Geonorge bathymetry points"

    ArchGDAL.read(geodatabase_path) do dataset
        point_count += sample_depth_layer!(
            point_layer,
            dataset,
            "dybdepunkt",
            transform,
            xmin,
            ymin,
            xmax,
            ymax;
            geometry = :point,
        )
        include_contours && (
            point_count += sample_depth_layer!(
                point_layer,
                dataset,
                "dybdekurve",
                transform,
                xmin,
                ymin,
                xmax,
                ymax;
                geometry = :line,
            )
        )
    end

    @info "Finished sampling Geonorge bathymetry points: $point_count points"
    return point_count
end

"""
    sample_land_features!(land_layer, geodatabase_path, transform, filter_bounds)

Append transformed land-related features to `land_layer`.
"""
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

"""
    sample_depth_layer!(point_layer, dataset, layer_name, transform, xmin, ymin, xmax, ymax; geometry)

Sample one depth-bearing layer from the Sjøkart geodatabase.
"""
function sample_depth_layer!(point_layer, dataset, layer_name, transform, xmin, ymin, xmax, ymax; geometry)
    layer = find_layer(dataset, layer_name)
    isnothing(layer) && return 0

    ArchGDAL.setspatialfilter!(layer, xmin, ymin, xmax, ymax)
    depth_index = ArchGDAL.findfieldindex(layer, "dybde", false)
    total_features = ArchGDAL.nfeature(layer, true)
    processed_features = 0
    progress_interval = max(cld(total_features, 20), 1)
    point_count = 0
    @info "Sampling $layer_name: 0 / $total_features features processed"

    for feature in layer
        processed_features += 1
        depth = ArchGDAL.getfield(feature, depth_index)
        if ismissing(depth)
            processed_features % progress_interval == 0 &&
                @info "Sampling $layer_name: $processed_features / $total_features features processed, $(total_features - processed_features) remaining"
            continue
        end

        bottom_height = -abs(Float64(depth))
        geometry == :point &&
            (point_count += add_point_geometry!(point_layer, ArchGDAL.getgeom(feature), transform, bottom_height))
        geometry == :line &&
            (point_count += add_linestring_points!(point_layer, ArchGDAL.getgeom(feature), transform, bottom_height))

        processed_features % progress_interval == 0 &&
            @info "Sampling $layer_name: $processed_features / $total_features features processed, $(total_features - processed_features) remaining"
    end

    @info "Finished sampling $layer_name: $processed_features / $total_features features processed, $point_count points added"
    return point_count
end

"""
    add_point_geometry!(point_layer, point_geometry, transform, bottom_height)

Append one transformed depth point to `point_layer`.
"""
function add_point_geometry!(point_layer, point_geometry, transform, bottom_height)
    point = ArchGDAL.clone(point_geometry)
    ArchGDAL.transform!(point, transform)
    longitude, latitude, _ = ArchGDAL.getpoint(point, 0)

    ArchGDAL.createfeature(point_layer) do feature
        ArchGDAL.setfield!(feature, ArchGDAL.findfieldindex(feature, "z"), bottom_height)
        ArchGDAL.setgeom!(feature, 0, ArchGDAL.createpoint(longitude, latitude))
        return nothing
    end

    return 1
end

"""
    add_linestring_points!(point_layer, line, transform, bottom_height)

Append each vertex of a contour line to `point_layer` with a constant
`bottom_height` attribute.
"""
function add_linestring_points!(point_layer, line, transform, bottom_height)
    if ArchGDAL.geomname(line) == "MULTILINESTRING"
        added = 0
        for geometry_index = 0:ArchGDAL.ngeom(line)-1
            added +=
                add_linestring_points!(point_layer, ArchGDAL.getgeom(line, geometry_index), transform, bottom_height)
        end
        return added
    end

    npoints = ArchGDAL.ngeom(line)
    added = 0

    for point_index = 0:npoints-1
        x, y, _ = ArchGDAL.getpoint(line, point_index)
        point = ArchGDAL.createpoint(x, y)
        ArchGDAL.transform!(point, transform)
        longitude, latitude, _ = ArchGDAL.getpoint(point, 0)

        ArchGDAL.createfeature(point_layer) do feature
            ArchGDAL.setfield!(feature, ArchGDAL.findfieldindex(feature, "z"), bottom_height)
            ArchGDAL.setgeom!(feature, 0, ArchGDAL.createpoint(longitude, latitude))
            return nothing
        end

        added += 1
    end

    return added
end

"""
    raster_to_bottom_height(band, Nx, Ny)

Read a GDAL raster band into the `(Nx, Ny)` orientation used by FjordSim's
bathymetry NetCDF files.
"""
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

"""
    transformed_filter_bounds(longitude, latitude)

Transform a WGS84 longitude/latitude bounding box into EPSG:25833 bounds for
spatial filtering of the raw Geonorge Sjøkart geodatabase.
"""
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

"""
    center_coordinates(domain, N)

Return `N` uniformly spaced cell-center coordinates for the interval
`domain = (lower, upper)`.
"""
function center_coordinates(domain, N)
    delta = domain_step(domain, N)
    return collect(range(domain[1] + delta / 2, step = delta, length = N))
end

domain_step(domain, N) = (domain[2] - domain[1]) / N

"""
    vertical_faces(grid)

Extract the physical vertical face coordinates from a live Oceananigans grid,
excluding halo entries.
"""
function vertical_faces(grid)
    z_faces_with_halo = getproperty(grid.z, Symbol("cᵃᵃᶠ"))
    first_index = grid.Hz + 1
    last_index = first_index + grid.Nz
    return collect(z_faces_with_halo[first_index:last_index])
end

"""
    expand_domain(domain, N, padding_cells)

Expand a uniformly spaced domain by `padding_cells` cell widths on both sides.
"""
function expand_domain(domain, N, padding_cells)
    delta = domain_step(domain, N)
    lower = domain[1] - padding_cells * delta
    upper = domain[2] + padding_cells * delta
    return (lower, upper)
end

"""
    geonorge_raw_filename(longitude, latitude, size)

Build a deterministic filename for a regional raw bathymetry cache file.
"""
function geonorge_raw_filename(longitude, latitude, size; include_contours::Bool = true)
    Nx, Ny, _ = size
    lon1 = replace(string(round(longitude[1], digits = 3)), '.' => 'p')
    lon2 = replace(string(round(longitude[2], digits = 3)), '.' => 'p')
    lat1 = replace(string(round(latitude[1], digits = 3)), '.' => 'p')
    lat2 = replace(string(round(latitude[2], digits = 3)), '.' => 'p')
    contour_suffix = include_contours ? "" : "_points_only"

    return "geonorge_sjokart_bathymetry_$(lon1)_$(lon2)_$(lat1)_$(lat2)_$(Nx)x$(Ny)$(contour_suffix).nc"
end

end  # module Bathymetry
