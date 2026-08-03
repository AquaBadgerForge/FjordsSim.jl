module NORA3

export NORA3PrescribedAtmosphere, NORA3PrescribedRadiation, MultiYearNORA3

using ...Utils: compute_faces
using ...Configs: AbstractAtmosphereConfig, atmosphere_path

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.OutputReaders: Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices, FieldTimeSeries
using NumericalEarth.Atmospheres: PrescribedAtmosphere, PrescribedPrecipitationFlux
using NumericalEarth.Radiations: PrescribedRadiation, SurfaceRadiationProperties, default_stefan_boltzmann_constant
using NumericalEarth.DataWrangling: compute_native_date_range, Metadata, Metadatum, native_times
using Adapt
using NCDatasets

import Oceananigans: location
import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend
# The dataset interface below is *extended*, not redefined: these are NumericalEarth's own generic
# functions, so its machinery dispatches into this backend rather than silently falling back to the
# generic `Metadata` defaults. A bare `using NumericalEarth` would make each definition a new
# module-local binding that shadows the upstream one — see `Datasets.jl` for the same import block.
import NumericalEarth.DataWrangling:
    all_dates,
    available_variables,
    dataset_variable_name,
    first_date,
    is_three_dimensional,
    last_date,
    metadata_filename

const NORA3_FILE = "NORA3.nc"

"""
    default_nora3_dataset()

`MultiYearNORA3` for the default NORA3 reanalysis file at
`~/FjordSim_data/NORA3/NORA3.nc`.
"""
default_nora3_dataset() = MultiYearNORA3(NORA3_FILE, joinpath(homedir(), "FjordSim_data", "NORA3"))

const NORA3_variable_names = (
    :freshwater_flux,
    :specific_humidity,
    :sea_level_pressure,
    :downwelling_longwave_radiation,
    :downwelling_shortwave_radiation,
    :temperature,
    :eastward_velocity,
    :northward_velocity,
)

# NumericalEarth's variable symbol => the prepared file's variable name. The prepared names are
# fixed by `FjordSim.Atmospheres.ATMOSPHERE_VARIABLES`, which this module cannot reference: it is
# included *before* that definition, being the read side the contract is written against.
const NORA3_dataset_variable_names = Dict(
    :freshwater_flux => "precipitation",
    :specific_humidity => "specific_humidity_2m",
    :sea_level_pressure => "air_pressure_at_sea_level",
    :downwelling_longwave_radiation => "lwrad",
    :downwelling_shortwave_radiation => "swrad",
    :temperature => "air_temperature_2m",
    :eastward_velocity => "u_wind_10m",
    :northward_velocity => "v_wind_10m",
)

"""
    MultiYearNORA3

A prepared atmosphere NetCDF, as a NumericalEarth dataset.

Holds only what the `Metadata` interface is asked for without reopening the file: where the file
is, its horizontal size and its full date axis.

# Fields
- `metadata_filename`: File name, relative to `default_download_directory`.
- `default_download_directory`: Directory holding it.
- `size`: `(lon, lat)` cell counts.
- `all_dates`: Every date on the file's time axis, in file order.
"""
struct MultiYearNORA3{D}
    metadata_filename::String
    default_download_directory::String
    size::NTuple{2, Int}
    all_dates::Vector{D}
end

"""
    MultiYearNORA3(config::AbstractAtmosphereConfig)

The prepared atmosphere file of a setup, as a reader dataset.

`default_nora3_dataset()` points at the shared `~/FjordSim_data/NORA3/NORA3.nc`, but a file
written by `prepare_atmosphere` is specific to one setup — it is regridded onto that setup's own
box — so a setup resolves its own through `atmosphere_path`. Pass the result to
`NORA3PrescribedAtmosphere` or `NORA3PrescribedRadiation` as their `dataset` keyword.
"""
function MultiYearNORA3(config::AbstractAtmosphereConfig)
    filepath = atmosphere_path(config)
    isfile(filepath) || error(
        "Prepared atmosphere $filepath does not exist. " *
        "Run `julia --project -m FjordSim prepare_atmosphere` for this setup first.",
    )

    return MultiYearNORA3(basename(filepath), dirname(filepath))
end

function MultiYearNORA3(metadata_filename::String, default_download_directory::String)
    filepath = joinpath(default_download_directory, metadata_filename)

    # The size comes from the dimensions rather than from a variable: `lon`, `lat` and `time` are
    # the part of the prepared layout every variable shares, so nothing here has to name one.
    dataset_size, all_dates = NCDataset(filepath) do ds
        (NTuple{2, Int}((ds.dim["lon"], ds.dim["lat"])), ds["time"][:])
    end

    return MultiYearNORA3(metadata_filename, default_download_directory, dataset_size, all_dates)
end  # function

available_variables(::MultiYearNORA3) = NORA3_variable_names

const NORA3Metadata{D} = Metadata{<:MultiYearNORA3,D}
const NORA3Metadatum = Metadatum{<:MultiYearNORA3}
Base.size(metadata::NORA3Metadata) = (metadata.dataset.size..., length(metadata.dates))
Base.size(metadata::NORA3Metadatum) = (metadata.dataset.size..., 1)

is_three_dimensional(data::NORA3Metadata) = false
location(::NORA3Metadata) = (Center, Center, Center)
dataset_variable_name(data::NORA3Metadata) = NORA3_dataset_variable_names[data.name]

all_dates(ds::MultiYearNORA3, name) = ds.all_dates
all_dates(ds::MultiYearNORA3) = ds.all_dates
first_date(ds::MultiYearNORA3) = first(all_dates(ds))
last_date(ds::MultiYearNORA3) = last(all_dates(ds))
metadata_filename(ds::MultiYearNORA3, args...) = ds.metadata_filename

"""
    NORA3NetCDFBackend(length)
    NORA3NetCDFBackend(start, length)

Sliding-window in-memory backend for a `FieldTimeSeries` reading a prepared atmosphere file,
holding `length` time slots starting at slot `start`.

`metadata` and `file_indices` are filled in by `NORA3FieldTimeSeries` and are `nothing` until
then, which is also what `Adapt` strips them to: neither is usable inside a GPU kernel, matching
`NumericalEarth.DataWrangling.DatasetBackend`.

`file_indices` maps a slot of the series to its record in the file. It is what keeps `set!`
honest when `start_date`/`end_date` select a sub-range: the slot axis then starts at 1 while the
records it stands for do not.
"""
struct NORA3NetCDFBackend{M, I} <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
    metadata::M
    file_indices::I
end

Adapt.adapt_structure(to, b::NORA3NetCDFBackend) =
    NORA3NetCDFBackend(b.start, b.length, nothing, nothing)

NORA3NetCDFBackend(start::Integer, length::Integer) =
    NORA3NetCDFBackend(start, length, nothing, nothing)
NORA3NetCDFBackend(length) = NORA3NetCDFBackend(1, length)

Base.length(backend::NORA3NetCDFBackend) = backend.length
Base.summary(backend::NORA3NetCDFBackend) = string("NORA3NetCDFBackend(", backend.start, ", ", backend.length, ")")

const NORA3NetCDFFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:NORA3NetCDFBackend}

new_backend(b::NORA3NetCDFBackend, start, length) =
    NORA3NetCDFBackend(start, length, b.metadata, b.file_indices)

"""
    nora3_time_indices(dataset, dates, name)

The file record each of `dates` lives at, in the same order.

A date the file does not carry is an error rather than a silent omission: dropping it would leave
the series one slot short of the axis it was built for, so every later slot would read the wrong
record. `NumericalEarth`'s multi-year JRA55 reader refuses the same way.
"""
function nora3_time_indices(dataset::MultiYearNORA3, dates, name)
    file_dates = all_dates(dataset, name)
    file_index = Dict(date => index for (index, date) in enumerate(file_dates))
    indices = Vector{Int}(undef, length(dates))

    for (slot, date) in enumerate(dates)
        index = get(file_index, date, 0)
        index == 0 && error(
            "$date is absent from the time axis of $(metadata_filename(dataset)). The dates " *
            "requested for :$name must all be present in the prepared atmosphere file.",
        )
        indices[slot] = index
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
    name = dataset_variable_name(metadata)

    # `time_indices` are slots of the series, which coincide with records of the file only when
    # the series spans the whole file. `file_indices` is what turns one into the other.
    nn = collect(time_indices(fts))
    file_indices = backend.file_indices

    data = NCDataset(filepath) do ds
        if issorted(nn)
            ds[name][:, :, file_indices[nn]]
        else
            # The time indices may be cycling past 1; eg ti = [6, 7, 8, 1].
            # However, DiskArrays does not seem to support loading data with unsorted
            # indices. So to handle this, we load the data in chunks, where each chunk's
            # indices are sorted, and then glue the data together. The slot -> record map is
            # monotone, so splitting on the slots splits the records too.
            m = findfirst(n -> n == 1, nn)
            data1 = ds[name][:, :, file_indices[nn[1:m-1]]]
            data2 = ds[name][:, :, file_indices[nn[m:end]]]
            cat(data1, data2, dims=3)
        end
    end

    copyto!(interior(fts, :, :, 1, :), data)
    fill_halo_regions!(fts)

    return nothing
end

"""
    NORA3FieldTimeSeries(metadata, architecture, FT; backend, time_indexing, reference_date)

A `FieldTimeSeries` over one variable of a prepared atmosphere file.

`reference_date` is the instant the time axis is zeroed at, defaulting to the first date the
metadata selects. It is *not* `start_date`: that one picks which records to load, this one picks
where t = 0 sits, and a coupled run needs every component zeroed at the same instant.
"""
function NORA3FieldTimeSeries(
    metadata::NORA3Metadata,
    architecture,
    FT;
    backend = NORA3NetCDFBackend(10),
    time_indexing = Cyclical(),
    reference_date = nothing,
)

    dataset = metadata.dataset
    name = metadata.name

    # Change the metadata to reflect the actual time indices
    file_indices = nora3_time_indices(dataset, metadata.dates, name)
    dates = all_dates(dataset, name)[file_indices]
    metadata = Metadata(metadata.name; dataset = metadata.dataset, dates, dir = metadata.dir)

    # A window wider than the series would wrap `time_indices` past 1 more than once, which the
    # split in `set!` cannot express and `copyto!` would silently accept as a short read.
    window = min(backend.length, length(dates))
    backend = NORA3NetCDFBackend(backend.start, window, metadata, file_indices)

    shortname = dataset_variable_name(metadata)
    filepath = joinpath(metadata.dataset.default_download_directory, metadata.dataset.metadata_filename)

    longitude, latitude = NCDataset(filepath) do ds
        (compute_faces(ds["lon"][:]), compute_faces(ds["lat"][:]))
    end

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
    start_time = isnothing(reference_date) ? first(metadata.dates) : reference_date
    # Match the time axis to the grid's float type, as `NumericalEarth`'s own reader does:
    # `native_times` returns `Float64` seconds, and against a `Float32` grid that makes
    # `interpolate`'s time weight `Float64`, boxing the interpolated value inside GPU kernels.
    times = convert.(eltype(grid), native_times(metadata; start_time))
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

"""
    NORA3PrescribedAtmosphere([architecture = CPU(), FT = Float32];
                              dataset = default_nora3_dataset(),
                              start_date = first_date(dataset),
                              end_date = last_date(dataset),
                              backend = NORA3NetCDFBackend(10),
                              time_indexing = Cyclical(),
                              reference_date = nothing,
                              surface_layer_height = 10,
                              other_kw...)

Return a `PrescribedAtmosphere` backed by the six prepared NORA3 surface fields.

`start_date` and `end_date` select which records to load; `reference_date` sets the instant the
time axis is zeroed at, defaulting to the first selected record. `surface_layer_height` is the
elevation the prescribed variables are taken to sit at, which the similarity-theory solver reads
as the height of its lowest atmospheric level.
"""
function NORA3PrescribedAtmosphere(
    architecture = CPU(),
    FT = Float32;
    dataset = default_nora3_dataset(),
    start_date = first_date(dataset),
    end_date = last_date(dataset),
    backend = NORA3NetCDFBackend(10),
    time_indexing = Cyclical(),
    reference_date = nothing,
    surface_layer_height = 10,  # meters
    other_kw...,
)

    kw = (; time_indexing, backend, start_date, end_date, dataset, reference_date)
    kw = merge(kw, other_kw)

    ua = NORA3FieldTimeSeries(:eastward_velocity, architecture, FT; kw...)
    va = NORA3FieldTimeSeries(:northward_velocity, architecture, FT; kw...)
    Ta = NORA3FieldTimeSeries(:temperature, architecture, FT; kw...)
    qa = NORA3FieldTimeSeries(:specific_humidity, architecture, FT; kw...)
    pa = NORA3FieldTimeSeries(:sea_level_pressure, architecture, FT; kw...)
    Fra = NORA3FieldTimeSeries(:freshwater_flux, architecture, FT; kw...)

    times = ua.times
    grid = ua.grid

    velocities = (u = ua, v = va)

    # NORA3 carries rain only, and `PrescribedPrecipitationFlux` reads a `nothing` snow component as
    # "this dataset does not represent snowfall", so there is no placeholder field to build.
    precipitation_flux = PrescribedPrecipitationFlux(; rain = Fra)

    surface_layer_height = convert(FT, surface_layer_height)

    # Temperature and specific humidity are their own keywords: `tracers` is for gas species (CO₂
    # and the like). Passing them as tracers is accepted silently and leaves the atmosphere on its
    # constant placeholder temperature and humidity.
    atmosphere = PrescribedAtmosphere(
        grid,
        times;
        velocities,
        temperature = Ta,
        specific_humidity = qa,
        pressure = pa,
        precipitation_flux,
        surface_layer_height,
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
                             reference_date = nothing,
                             ocean_surface = SurfaceRadiationProperties(0.05, 0.97),
                             sea_ice_surface = SurfaceRadiationProperties(0.7, 1.0),
                             stefan_boltzmann_constant = default_stefan_boltzmann_constant,
                             other_kw...)

Return a `PrescribedRadiation` backed by NORA3 downwelling shortwave and
longwave `NORA3FieldTimeSeries`.

Both fluxes are downwelling, not net: `PrescribedRadiation` applies `ocean_surface`'s own albedo
and emissivity, so a net flux would count the albedo twice. `reference_date` zeroes the time axis
exactly as it does for `NORA3PrescribedAtmosphere`.
"""
function NORA3PrescribedRadiation(
    architecture = CPU(),
    FT = Float32;
    dataset = default_nora3_dataset(),
    start_date = first_date(dataset),
    end_date = last_date(dataset),
    backend = NORA3NetCDFBackend(10),
    time_indexing = Cyclical(),
    reference_date = nothing,
    ocean_surface = SurfaceRadiationProperties(0.05, 0.97),
    sea_ice_surface = SurfaceRadiationProperties(0.7, 1.0),
    stefan_boltzmann_constant = default_stefan_boltzmann_constant,
    other_kw...,
)
    kw = (; time_indexing, backend, start_date, end_date, dataset, reference_date)
    kw = merge(kw, other_kw)

    Qs = NORA3FieldTimeSeries(:downwelling_shortwave_radiation, architecture, FT; kw...)
    Ql = NORA3FieldTimeSeries(:downwelling_longwave_radiation,  architecture, FT; kw...)

    return PrescribedRadiation(Qs, Ql; ocean_surface, sea_ice_surface, stefan_boltzmann_constant)
end # function

end # module