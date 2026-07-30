module NORA3

export NORA3PrescribedAtmosphere, NORA3PrescribedRadiation, MultiYearNORA3

using ...Utils: compute_faces

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices, FieldTimeSeries
using NumericalEarth
using NumericalEarth: PrescribedAtmosphere
using NumericalEarth.Radiations: PrescribedRadiation, SurfaceRadiationProperties, default_stefan_boltzmann_constant
using NumericalEarth.DataWrangling: compute_native_date_range, Metadata, metadata_path, native_times
using Adapt
using NCDatasets
using JLD2
using Dates

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!
import NumericalEarth: all_dates
import NumericalEarth.DataWrangling: metadata_filename

const NORA3_FILE = "NORA3.nc"

"""
    default_nora3_dataset()

`MultiYearNORA3` for the default NORA3 reanalysis file at
`~/FjordSim_data/NORA3/NORA3.nc`.
"""
default_nora3_dataset() = MultiYearNORA3(NORA3_FILE, joinpath(homedir(), "FjordSim_data", "NORA3"))

NORA3_variable_names = (
    :freshwater_flux,
    :specific_humidity,
    :sea_level_pressure,
    :downwelling_longwave_radiation,
    :downwelling_shortwave_radiation,
    :temperature,
    :eastward_velocity,
    :northward_velocity,
)

struct MultiYearNORA3{D}
    metadata_filename::String
    default_download_directory::String
    size::NTuple{2, Int}
    all_dates::Vector{D}
end

function MultiYearNORA3(metadata_filename::String, default_download_directory::String)
    filepath = joinpath(default_download_directory, metadata_filename)
    ds = NCDataset(filepath)
    sz = NTuple{2, Int}(size(ds["air_temperature_2m"])[1:2])
    all_dates = ds["time"][:]
    close(ds)

    return MultiYearNORA3(metadata_filename, default_download_directory, sz, all_dates)
end  # function

available_variables(::MultiYearNORA3) = NORA3_variable_names

const NORA3Metadata{D} = Metadata{<:MultiYearNORA3,D}
const NORA3Metadatum = Metadatum{<:MultiYearNORA3}
Base.size(metadata::NORA3Metadata) = (metadata.dataset.size..., length(metadata.dates))
Base.size(::NORA3Metadatum) = (metadata.dataset.size..., 1)

is_three_dimensional(data::NORA3Metadata) = false
location(::NORA3Metadata) = (Center, Center, Center)
dataset_variable_name(data::NORA3Metadata) = NORA3_dataset_variable_names[data.name]

NORA3_dataset_variable_names = Dict(
    :freshwater_flux => "precipitation",   # Freshwater fluxes
    :specific_humidity => "specific_humidity_2m",     # Surface specific humidity
    :sea_level_pressure => "air_pressure_at_sea_level",      # Sea level pressure
    :downwelling_longwave_radiation => "lwrad",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "swrad",     # Downwelling shortwave radiation
    :temperature => "air_temperature_2m",      # Near-surface air temperature
    :eastward_velocity => "u_wind_10m",      # Eastward near-surface wind
    :northward_velocity => "v_wind_10m",      # Northward near-surface wind
)

all_dates(ds::MultiYearNORA3, name) = ds.all_dates
all_dates(ds::MultiYearNORA3) = ds.all_dates
first_date(ds::MultiYearNORA3) = first(all_dates(ds))
last_date(ds::MultiYearNORA3) = last(all_dates(ds))
metadata_filename(ds::MultiYearNORA3, args...) = ds.metadata_filename

# this is field data series stuff to get data during simulation
struct NORA3NetCDFBackend{M} <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
    metadata::M
end

Adapt.adapt_structure(to, b::NORA3NetCDFBackend) = NORA3NetCDFBackend(b.start, b.length, nothing)

NORA3NetCDFBackend(length, metadata::Metadata) = NORA3NetCDFBackend(1, length, metadata)
NORA3NetCDFBackend(start::Integer, length::Integer) = NORA3NetCDFBackend(start, length, nothing)
NORA3NetCDFBackend(length) = NORA3NetCDFBackend(1, length, nothing)

Base.length(backend::NORA3NetCDFBackend) = backend.length
Base.summary(backend::NORA3NetCDFBackend) = string("NORA3NetCDFBackend(", backend.start, ", ", backend.length, ")")

const NORA3NetCDFFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:NORA3NetCDFBackend}

new_backend(b::NORA3NetCDFBackend, start, length) = NORA3NetCDFBackend(start, length, b.metadata)

function NORA3_time_indices(ds::MultiYearNORA3, dates, name)
    nora3_all_dates = all_dates(ds, name)
    indices = sizehint!(Int[], length(dates))

    for date in dates
        index = findfirst(x -> x == date, nora3_all_dates)
        !isnothing(index) && push!(indices, index)
    end

    return indices
end

function NORA3FieldTimeSeries(variable_name::Symbol, architecture, FT; dataset, start_date, end_date, kw...)

    native_dates = all_dates(dataset, variable_name)
    dates = compute_native_date_range(native_dates, start_date, end_date)
    metadata = Metadata(variable_name; dataset, dates, dir = dataset.default_download_directory)

    return NORA3FieldTimeSeries(metadata, architecture, FT; kw...)
end

function set!(fts::NORA3NetCDFFTS, backend=fts.backend)

    metadata = backend.metadata
    filepath = joinpath(metadata.dataset.default_download_directory, metadata.dataset.metadata_filename)
    ds = Dataset(filepath)

    nn   = time_indices(fts)
    nn   = collect(nn)
    name = dataset_variable_name(fts.backend.metadata)

    if issorted(nn)
        data = ds[name][:, :, nn]
    else
        # The time indices may be cycling past 1; eg ti = [6, 7, 8, 1].
        # However, DiskArrays does not seem to support loading data with unsorted
        # indices. So to handle this, we load the data in chunks, where each chunk's
        # indices are sorted, and then glue the data together.
        m = findfirst(n -> n == 1, nn)
        n1 = nn[1:m-1]
        n2 = nn[m:end]

        data1 = ds[name][:, :, n1]
        data2 = ds[name][:, :, n2]
        data = cat(data1, data2, dims=3)
    end

    close(ds)

    copyto!(interior(fts, :, :, 1, :), data)
    fill_halo_regions!(fts)

    return nothing
end

function NORA3FieldTimeSeries(
    metadata::NORA3Metadata,
    architecture,
    FT;
    backend = NORA3NetCDFBackend(10),
    time_indexing = Cyclical(),
)

    dataset = metadata.dataset
    name = metadata.name

    # Change the metadata to reflect the actual time indices
    time_indices = NORA3_time_indices(dataset, metadata.dates, name)
    dates = all_dates(dataset, name)[time_indices]
    metadata = Metadata(metadata.name; dataset = metadata.dataset, dates, dir = metadata.dir)

    if backend.metadata isa Nothing
        backend = NORA3NetCDFBackend(backend.start, backend.length, metadata)
    end

    shortname = dataset_variable_name(metadata)
    filepath = joinpath(metadata.dataset.default_download_directory, metadata.dataset.metadata_filename)

    ds = Dataset(filepath)
    latitude = compute_faces(ds["lat"][:])
    longitude = compute_faces(ds["lon"][:])
    close(ds)

    Nrx = length(longitude) - 1
    Nry = length(latitude) - 1
    N = (Nrx, Nry)
    H = min.(N, (3, 3))

    grid = LatitudeLongitudeGrid(
        architecture,
        FT;
        halo = H,
        size = N,
        longitude = longitude,
        latitude = latitude,
        topology = (Bounded, Bounded, Flat),
    )
    boundary_conditions = FieldBoundaryConditions(grid, (Center(), Center(), nothing))
    start_time = first_date(metadata.dataset)
    times = native_times(metadata; start_time)
    fts = FieldTimeSeries{Center,Center,Nothing}(
        grid,
        times;
        backend,
        time_indexing,
        boundary_conditions,
        path = filepath,
        name = shortname,
    )
    set!(fts)

    return fts
end

function NORA3PrescribedAtmosphere(
    architecture = CPU(),
    FT = Float32;
    dataset = default_nora3_dataset(),
    start_date = first_date(dataset),
    end_date = last_date(dataset),
    backend = NORA3NetCDFBackend(10),
    time_indexing = Cyclical(),
    surface_layer_height = 10,  # meters
    other_kw...,
)

    kw = (; time_indexing, backend, start_date, end_date, dataset)
    kw = merge(kw, other_kw)

    ua = NORA3FieldTimeSeries(:eastward_velocity, architecture, FT; kw...)
    va = NORA3FieldTimeSeries(:northward_velocity, architecture, FT; kw...)
    Ta = NORA3FieldTimeSeries(:temperature, architecture, FT; kw...)
    qa = NORA3FieldTimeSeries(:specific_humidity, architecture, FT; kw...)
    pa = NORA3FieldTimeSeries(:sea_level_pressure, architecture, FT; kw...)
    Fra = NORA3FieldTimeSeries(:freshwater_flux, architecture, FT; kw...)

    times = ua.times
    grid = ua.grid

    freshwater_flux = (rain = Fra, snow = FieldTimeSeries{Center, Center, Nothing}(grid, times))

    velocities = (u = ua, v = va)

    tracers = (T = Ta, q = qa)

    pressure = pa

    FT = eltype(ua)
    surface_layer_height = convert(FT, surface_layer_height)

    atmosphere = PrescribedAtmosphere(
        grid,
        times;
        velocities,
        freshwater_flux,
        tracers,
        surface_layer_height,
        pressure,
    )

    return atmosphere
end # function

"""
    NORA3PrescribedRadiation([architecture = CPU(), FT = Float32];
                             dataset = default_nora3_dataset(),
                             start_date = first_date(dataset),
                             end_date = last_date(dataset),
                             backend = NORA3NetCDFBackend(10),
                             time_indexing = Cyclical(),
                             ocean_surface = SurfaceRadiationProperties(0.05, 0.97),
                             sea_ice_surface = SurfaceRadiationProperties(0.7, 1.0),
                             stefan_boltzmann_constant = default_stefan_boltzmann_constant,
                             other_kw...)

Return a `PrescribedRadiation` backed by NORA3 downwelling shortwave and
longwave `NORA3FieldTimeSeries`.
"""
function NORA3PrescribedRadiation(
    architecture = CPU(),
    FT = Float32;
    dataset = default_nora3_dataset(),
    start_date = first_date(dataset),
    end_date = last_date(dataset),
    backend = NORA3NetCDFBackend(10),
    time_indexing = Cyclical(),
    ocean_surface = SurfaceRadiationProperties(0.05, 0.97),
    sea_ice_surface = SurfaceRadiationProperties(0.7, 1.0),
    stefan_boltzmann_constant = default_stefan_boltzmann_constant,
    other_kw...,
)
    kw = (; time_indexing, backend, start_date, end_date, dataset)
    kw = merge(kw, other_kw)

    Qs = NORA3FieldTimeSeries(:downwelling_shortwave_radiation, architecture, FT; kw...)
    Ql = NORA3FieldTimeSeries(:downwelling_longwave_radiation,  architecture, FT; kw...)

    return PrescribedRadiation(Qs, Ql; ocean_surface, sea_ice_surface, stefan_boltzmann_constant)
end # function

end # module