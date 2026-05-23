module FDatasets

export DSForcing, DSResults, last_date

using Dates

using NCDatasets
using Oceananigans
using NumericalEarth.DataWrangling: Metadatum, metadata_path

import Oceananigans: location
import NumericalEarth.DataWrangling:
    metadata_filename,
    default_download_directory,
    all_dates,
    first_date,
    last_date,
    dataset_variable_name,
    download_dataset,
    longitude_interfaces,
    latitude_interfaces,
    z_interfaces,
    reversed_vertical_axis,
    inpainted_metadata_path,
    retrieve_data,
    default_inpainting

struct DSForcing{D, Z}
    metadata_filename::String
    default_download_directory::String
    reversed_vertical_axis::Bool
    longitude_interfaces::Tuple{Float64, Float64}
    latitude_interfaces::Tuple{Float64, Float64}
    size::NTuple{3, Int}
    all_dates::Vector{D}
    first_date::D
    last_date::D
    z_interfaces::Vector{Z}
end

function DSForcing(metadata_filename, default_download_directory)
    filepath = joinpath(default_download_directory, metadata_filename)
    ds = NCDataset(filepath)
    longitude_interfaces = (Float64(ds["Nx"][1]), Float64(ds["Nx"][end]))
    latitude_interfaces = (Float64(ds["Ny"][1]), Float64(ds["Ny"][end]))
    sz = NTuple{3, Int}((ds.dim["Nx"], ds.dim["Ny"], ds.dim["Nz"]))
    all_dates = ds["time"][:]
    first_date = ds["time"][1]
    last_date = ds["time"][end]
    z_interfaces = ds["Nz_faces"][:]
    close(ds)

    return DSForcing(
        metadata_filename,
        default_download_directory,
        false,
        longitude_interfaces,
        latitude_interfaces,
        sz,
        all_dates,
        first_date,
        last_date,
        z_interfaces,
    )
end  # function

struct DSResults{D, Z}
    metadata_filename::String
    default_download_directory::String
    reversed_vertical_axis::Bool
    longitude_interfaces::Tuple{Float64, Float64}
    latitude_interfaces::Tuple{Float64, Float64}
    size::NTuple{3, Int}
    all_dates::Vector{D}
    first_date::D
    last_date::D
    z_interfaces::Vector{Z}
end

"""
    DSResults(metadata_filename, default_download_directory; start_date_time)

Load and save info from a NetCDF file.

# Arguments
- `metadata_filename::String`
    Name of the NetCDF file.

- `default_download_directory::String`
    Directory where the file is stored.

# Keyword Arguments
- `start_date_time::DateTime`
    DateTime of the first record in the file at
    `joinpath(default_download_directory, metadata_filename)`.

"""
function DSResults(
    metadata_filename::String,
    default_download_directory::String;
    start_date_time::DateTime
)
    filepath = joinpath(default_download_directory, metadata_filename)
    ds = NCDataset(filepath)
    longitude_interfaces = (Float64(ds["λ_faa"][1]), Float64(ds["λ_faa"][end]))
    latitude_interfaces = (Float64(ds["φ_afa"][1]), Float64(ds["φ_afa"][end]))
    sz = NTuple{3, Int}(size(ds["T"])[1:3])  # that is for a grid, so size should be for tracers
    all_dates = start_date_time .+ Millisecond.(round.(Int, ds["time"][:] .* 1000))
    first_date = all_dates[1]
    last_date = all_dates[end]
    z_interfaces = ds["z_aaf"][:]
    close(ds)

    return DSResults(
        metadata_filename,
        default_download_directory,
        false,
        longitude_interfaces,
        latitude_interfaces,
        sz,
        all_dates,
        first_date,
        last_date,
        z_interfaces,
    )
end  # function

Variable_names = Dict(:temperature => "T", :salinity => "S", :u_velocity => "u", :v_velocity => "v")
Variable_location = Dict(
    :temperature => (Center, Center, Center),
    :salinity => (Center, Center, Center),
    :free_surface => (Center, Center, Nothing),
    :sea_ice_thickness => (Center, Center, Nothing),
    :sea_ice_concentration => (Center, Center, Nothing),
    :net_heat_flux => (Center, Center, Nothing),
    :u_velocity => (Face, Center, Center),
    :v_velocity => (Center, Face, Center),
    :sensible_heat_flux => (Center, Center, Nothing),
    :latent_heat_flux => (Center, Center, Nothing),
    :net_longwave => (Center, Center, Nothing),
    :downwelling_shortwave => (Center, Center, Nothing),
    :downwelling_longwave => (Center, Center, Nothing),
)

const MetadatumForcing = Metadatum{<:DSForcing,<:Any,<:Any}
const MetadatumResults = Metadatum{<:DSResults,<:Any,<:Any}

# Metadatum constructors query metadata_filename(dataset, ...) with extra context arguments.
metadata_filename(ds::Union{DSResults, DSForcing}, args...) = ds.metadata_filename
metadata_filename(metadata::Union{MetadatumForcing, MetadatumResults}) = metadata.dataset.metadata_filename
default_download_directory(ds::Union{DSResults, DSForcing}) = ds.default_download_directory
reversed_vertical_axis(ds::Union{DSResults, DSForcing}) = ds.reversed_vertical_axis
longitude_interfaces(ds::Union{DSResults, DSForcing}) = ds.longitude_interfaces
latitude_interfaces(ds::Union{DSResults, DSForcing}) = ds.latitude_interfaces

Base.size(ds::Union{DSResults, DSForcing}) = ds.size
Base.size(ds::Union{DSResults, DSForcing}, variable) = size(ds)

all_dates(ds::Union{DSResults, DSForcing}, args...) = ds.all_dates
first_date(ds::Union{DSResults, DSForcing}, args...) = ds.first_date
last_date(ds::Union{DSResults, DSForcing}, args...) = ds.last_date

z_interfaces(metadata::Union{MetadatumResults, MetadatumForcing}) = metadata.dataset.z_interfaces

function dataset_variable_name(metadata::Union{MetadatumResults, MetadatumForcing})
    if haskey(Variable_names, metadata.name)
        return Variable_names[metadata.name]
    end
    return metadata.name
end

function location(metadata::Union{MetadatumResults, MetadatumForcing})
    if haskey(Variable_location, metadata.name)
        return Variable_location[metadata.name]
    end
    # assume a tracer
    @info string("Assuming ", metadata.name, " is at (Center, Center, Center).")
    return (Center, Center, Center)
end

function download_dataset(metadata::Union{MetadatumResults, MetadatumForcing})
    filepath = metadata_path(metadata)
    return filepath
end

function inpainted_metadata_filename(metadata::Union{MetadatumResults, MetadatumForcing})
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::Union{MetadatumResults, MetadatumForcing}) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

"""
    retrieve_data(metadata)

Retrieve data from netcdf file according to `metadata`.
"""
function retrieve_data(metadata::Union{Metadatum{V} where V<:DSResults, Metadatum{V} where V<:DSForcing})
    path = metadata_path(metadata)
    name = dataset_variable_name(metadata)
    # NetCDF shenanigans
    ds = Dataset(path)
    # find a proper datetime index - that is a difference
    i = findfirst(==(metadata.dates), metadata.dataset.all_dates)
    data = ds[name][:, :, :, i]
    close(ds)
    return data
end

function default_inpainting(metadata::Union{MetadatumResults})
    # no inpainting
    return nothing
end

end  # module
