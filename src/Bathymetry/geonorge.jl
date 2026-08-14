# Geonorge Sjøkart Dybdedata adapter: the source-specific half of the bathymetry pipeline.
# Copy this file's shape to add a new bathymetry source — it owns a config subtype, its own
# download and sampling, and the two hook methods `bathymetry_dataset` and `regrid_options`.
# Everything else comes from the generic core in Bathymetry.jl.

using Scratch
using ArchGDAL
using Downloads
using p7zip_jll

using NumericalEarth.DataWrangling: AbstractStaticBathymetry

import NumericalEarth.DataWrangling:
    dataset_variable_name,
    default_download_directory,
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

download_bathymetry_cache::String = ""

function __init__()
    global download_bathymetry_cache = @get_scratch!("Bathymetry")
end

# --- Bathymetry config ---

"""
    DybdedataConfig

Configuration for `prepare_bathymetry` driven by the Geonorge Sjøkart Dybdedata FileGDB
dataset, which it reaches through `bathymetry_dataset` and `regrid_options`.

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
- `minimum_depth`: Passed to `NumericalEarth.regrid_bathymetry`, as a positive depth in metres.
  Every wet cell shallower than this is *deepened* to exactly it (`min(h, -minimum_depth)`), so
  no water area is lost and the coastline does not move — note that `regrid_bathymetry`'s own
  docstring claims such cells become land, which its kernel does not do.

  Defaults to `0.0`, which keeps every cell the source resolves — including the sub-millimetre
  slivers interpolation leaves along a coastline, which `PartialCellBottom` then floors at
  `minimum_fractional_cell_height * Δz`. Such a sliver destabilizes its *neighbours*, not itself:
  the horizontal spacings are `Δrᶠᶜᶜ = min(Δrᶜᶜᶜ(i-1,j,k), Δrᶜᶜᶜ(i,j,k))`, so the face between a
  floored 0.2 m cell and the ordinary cell beside it is 0.2 m thick while the volume it feeds is
  metres thick, and `cell_advection_timescale` builds its vertical term from `Δzᶜᶜᶠ` and never
  sees the thin cell at all. The time-step wizard therefore under-constrains `Δt` there, velocity
  grows monotonically in the neighbouring cell over ~1500 iterations, and the run dies with
  `InexactError: Int64(NaN)` from `calculate_substeps` once the timescale reduction goes NaN — the
  last symptom rather than the cause. On `oslofjorden` a 9.4 cm sliver did this to the 4.9 m cell
  beside it. Set it to a depth that leaves no sliver.

  A floor alone is not enough, and on its own makes a second problem: lifting a sliver to a
  constant depth in much deeper water leaves a shallow *spike*, which destabilizes its neighbours
  the same way. Pair it with `spike_ratio`.
- `close_narrow_passages`: Run `remove_narrow_passages`, which turns every one-cell-wide sea passage
  whose removal leaves both of its sides connected into land. `false` disables the stage.

  This is a *width* problem, and the only one of the four stages that moves the coastline. Regridding
  leaves a one-cell channel wherever a strait too narrow to resolve cuts through a peninsula; the two
  basins it joins are usually already connected elsewhere, so it closes a loop, and the barotropic
  head difference around that loop is forced through a cross-section one cell wide and a few metres
  deep. On `oslofjorden` there were 66 such passages, 23 of them at exactly the `minimum_depth` floor,
  and one of them — a 2.35 m canal through a peninsula at 10.43°E, 59.09°N — carried a coherent
  47 m s⁻¹ barotropic jet that collapsed the time step from three minutes to 0.3 s. Neither
  `spike_ratio` nor `max_slope_factor` touches them, because both bound depth *contrast* and such a
  cell agrees with its neighbours.

  A passage that is the sole link to a basin is kept, so no water is deleted from the domain.
- `spike_ratio`: `fill_shallow_spikes` threshold — a sea cell shallower than this fraction of its
  neighbours' median depth is replaced by that median. `0.0` disables despiking. This is what
  removes the interpolation spikes, and the ones `minimum_depth` itself creates.
- `max_slope_factor`: `limit_bottom_slope` limit on the Beckmann–Haidvogel slope parameter
  `r = |d₁ - d₂| / (d₁ + d₂)` between adjacent sea cells. `0.0` disables slope limiting. Despiking
  handles the one-cell artifacts; this bounds what is left, including real topography. Lower is
  more stable and less faithful — 0.2 is the textbook limit but flattens genuine fjord sills, so
  prefer the largest value that keeps the run stable.
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
    minimum_depth::Float64 = 0.0
    close_narrow_passages::Bool = false
    spike_ratio::Float64 = 0.0
    max_slope_factor::Float64 = 0.0
    geonorge_cache::Bool = true
    regrid_cache::Bool = true
    raw_file::String = ""
    longitude::NTuple{2,Float64} = (0.0, 0.0)
    latitude::NTuple{2,Float64} = (0.0, 0.0)
    size::NTuple{3,Int} = (0, 0, 1)
    filter_bounds::NTuple{4,Float64} = (0.0, 0.0, 0.0, 0.0)
end

"""
    bathymetry_dataset(target_grid, config::DybdedataConfig)

Materialize the regional raw Geonorge dataset for `target_grid`: derive the native region,
download and extract the FileGDB if absent, and build (or reuse) the cached raw NetCDF.

`config` is mutated by `native_region!` to record the native region derived from `target_grid`.
"""
function bathymetry_dataset(target_grid, config::DybdedataConfig)
    native_region!(config, target_grid)
    ensure_geodatabase(config)

    @info "Preparing Geonorge bathymetry"
    return geonorge_dataset(config)
end

"""
    regrid_options(config::DybdedataConfig)

The `NumericalEarth.regrid_bathymetry` options this setup configures.
"""
regrid_options(config::DybdedataConfig) = (;
    cache = config.regrid_cache,
    interpolation_passes = config.interpolation_passes,
    major_basins = config.major_basins,
    minimum_depth = config.minimum_depth,
)

"""
    smoothing_options(config::DybdedataConfig)

The `smooth_bathymetry_gaps!` options this setup configures.

`minimum_depth` is passed on as well as being a `regrid_bathymetry` option, because
`limit_bottom_slope` must not undo the floor while flattening a slope.
"""
smoothing_options(config::DybdedataConfig) = (;
    close_narrow_passages = config.close_narrow_passages,
    spike_ratio = config.spike_ratio,
    max_slope_factor = config.max_slope_factor,
    minimum_depth = config.minimum_depth,
)

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
    Downloads.download(metadata::GeonorgeBathymetryMetadatum)

Return the generated regional raw bathymetry file.
This dataset is materialized locally by `bathymetry_dataset`.
"""
function Downloads.download(metadata::GeonorgeBathymetryMetadatum)
    filepath = metadata_path(metadata)
    isfile(filepath) || error("Raw bathymetry file $filepath does not exist. Run prepare_bathymetry first.")
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

    n_points = ArchGDAL.ngeom(line)

    for point_index in contour_point_indices(n_points, config.contour_stride)
        x, y, _ = ArchGDAL.getpoint(line, point_index)
        push!(samples.xs, x)
        push!(samples.ys, y)
        push!(samples.bottom_heights, bottom_height)
    end

    return nothing
end

function contour_point_indices(n_points, stride)
    last_index = n_points - 1
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

function geonorge_raw_filename(config::DybdedataConfig)
    Nx, Ny, _ = config.size
    lon1 = replace(string(round(config.longitude[1], digits = 3)), '.' => 'p')
    lon2 = replace(string(round(config.longitude[2], digits = 3)), '.' => 'p')
    lat1 = replace(string(round(config.latitude[1], digits = 3)), '.' => 'p')
    lat2 = replace(string(round(config.latitude[2], digits = 3)), '.' => 'p')
    contour_suffix = config.include_contours ? "_contour_stride_$(config.contour_stride)" : "_points_only"
    return "geonorge_sjokart_bathymetry_$(lon1)_$(lon2)_$(lat1)_$(lat2)_$(Nx)x$(Ny)$(contour_suffix).nc"
end
