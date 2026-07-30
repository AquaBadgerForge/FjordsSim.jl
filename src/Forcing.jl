module Forcing

using Oceananigans
using Oceananigans.Units: Time
using Oceananigans.OutputReaders: Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Operators: Ax, Ay, Az, volume
using Oceananigans: fill_halo_regions!, nodes, interior
using Dates: Second
using Adapt
using NCDatasets

import Oceananigans: on_architecture
import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend

using ..Configs: AbstractForcingConfig

export forcing_from_file, norkyst_directory, norkyst_monthly_filename, NorKystConfig

const NORKYST_CATALOG_URL = "https://thredds.met.no/thredds/catalog/fou-hi/norkyst800m/catalog.xml"
const NORKYST_OPENDAP_URL = "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m/"
const NORKYST_DIRECTORY = "norkyst"
const NORKYST_PARAMETERS = ["temperature", "salinity", "u_eastward", "v_northward"]

"""
    norkyst_monthly_filename(year, month)

Name of the combined monthly NorKyst-800m NetCDF file written for `year` and `month`.
"""
norkyst_monthly_filename(year, month) = "NorKyst-800m_ZDEPTHS_avg_$(year)$(lpad(month, 2, '0')).nc"

"""
    NorKystConfig

Configuration for downloading and subsetting NorKyst-800m reanalysis data.

Only `data_root` and `years` are required; the remaining fields default to the public
THREDDS endpoints and a `norkyst` subdirectory of `data_root`.

`output_directory` is a name relative to `data_root`, resolved by `norkyst_directory`.
Setting it to an absolute path overrides `data_root`.

# Fields
- `data_root`: Directory holding this setup's forcing files.
- `output_directory`: Name of the directory where monthly NetCDF files are written.
- `catalog_url`: THREDDS catalog URL listing available files.
- `opendap_url`: OPeNDAP base URL for streaming data.
- `parameters`: Variable names to extract (e.g. `["temperature", "salinity"]`).
- `years`: Calendar years to download.
"""
Base.@kwdef mutable struct NorKystConfig <: AbstractForcingConfig
    data_root::String
    output_directory::String = NORKYST_DIRECTORY
    catalog_url::String = NORKYST_CATALOG_URL
    opendap_url::String = NORKYST_OPENDAP_URL
    parameters::Vector{String} = copy(NORKYST_PARAMETERS)
    years::Vector{Int}
end

"""
    norkyst_directory(config::NorKystConfig)

Resolve `config.output_directory` against `config.data_root`. An absolute
`output_directory` is returned unchanged.
"""
norkyst_directory(config::NorKystConfig) = joinpath(config.data_root, config.output_directory)

""" Custom backend for FieldTimeSeries """
struct NetCDFBackend <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
end

NetCDFBackend(length) = NetCDFBackend(1, length)
new_backend(::NetCDFBackend, start, length) = NetCDFBackend(start, length)

Base.length(backend::NetCDFBackend) = backend.length
Base.summary(backend::NetCDFBackend) = string("NetCDFBackend(", backend.start, ", ", backend.length, ")")

const NetCDFFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:NetCDFBackend}

function get_data_location(var::Symbol)
    if var in (:u,)
        return (Face, Center, Center)
    elseif var in (:v,)
        return (Center, Face, Center)
    else
        return (Center, Center, Center)
    end
end

# Runtime generation of dictionaries based on forcing_variables_names
function get_oceananigans_fieldname(vars)
    Dict(Symbol(var) => Val{Symbol(var)}() for var in vars)
end

function get_data_location_dict(vars)
    Dict(Val{Symbol(var)}() => get_data_location(Symbol(var)) for var in vars)
end

# Generic Base.getindex for any variable using Val types
# fields[:VARIABLE] doesn't work in CUDA kernels
@inline Base.getindex(fields, i, j, k, ::Val{name}) where {name} = @inbounds getfield(fields, name)[i, j, k]

""" A custom forcing callable structure """
struct ForcingFromFile{FTS,V}
    fts_value::FTS
    fts_λ::FTS
    fieldname::V
end

Adapt.adapt_structure(to, p::ForcingFromFile) =
    ForcingFromFile(Adapt.adapt(to, p.fts_value), Adapt.adapt(to, p.fts_λ), Adapt.adapt(to, p.fieldname))

on_architecture(to, forcing::ForcingFromFile) = ForcingFromFile(
    on_architecture(to, forcing.fts_value),
    on_architecture(to, forcing.fts_λ),
    on_architecture(to, forcing.fieldname),
)

""" x direction """
function forcing_term_u(λ, flux, i, j, k, grid, args...)
    return flux * Ax(i, j, k, grid, Center(), Center(), Center()) / volume(i, j, k, grid, Center(), Center(), Center())
end

""" y direction """
function forcing_term_v(λ, flux, i, j, k, grid, args...)
    return flux * Ay(i, j, k, grid, Center(), Center(), Center()) / volume(i, j, k, grid, Center(), Center(), Center())
end

""" relaxation term """
function forcing_term_relax(λ, value, i, j, k, grid, field)
    return -λ * (field - value)
end

""" 
Return a result to be added to the tendency contributions (a kernel function).
"""
@inline function (p::ForcingFromFile{FTS,V})(i, j, k, grid, clock, fields) where {FTS,V}
    value = @inbounds p.fts_value[i, j, k, Time(clock.time)]
    λ = @inbounds p.fts_λ[i, j, k, Time(clock.time)]
    result = 0.0
    result += @inbounds ifelse(
        λ > 1 && value > -990,
        forcing_term_u(λ, value, i, j, k, grid, fields[i, j, k, p.fieldname]),
        0,
    )
    result += @inbounds ifelse(
        λ < -1 && value > -990,
        forcing_term_v(λ, value, i, j, k, grid, fields[i, j, k, p.fieldname]),
        0,
    )
    result += @inbounds ifelse(
        -1 < λ < 1 && value > -990,
        forcing_term_relax(λ, value, i, j, k, grid, fields[i, j, k, p.fieldname]),
        0,
    )
    return result
end

regularize_forcing(forcing::ForcingFromFile, field, field_name, model_field_names) = forcing

function native_times_to_seconds(native_times, start_time = native_times[1])
    return [Second(t - start_time).value for t in native_times]
end

function load_from_netcdf(; path::String, var_name::String, grid_size::Tuple, time_indices_in_memory::Tuple)
    ds = NCDataset(path)
    var = coalesce.(ds[var_name], -999.0)
    native_times = ds["time"]

    data = zeros(Float64, (grid_size[1:end]..., length(time_indices_in_memory)))
    j = 1
    for i in time_indices_in_memory
        @views data[:, :, :, j] .= var[:, :, :, i]
        j += 1
    end
    times = convert.(Int64, native_times_to_seconds(native_times))

    close(ds)
    return data, times
end

""" Update data in the FieldTimeSeries, e.g. in ForcingFromFile forcing structures. """
function set!(fts::NetCDFFTS, path::String = fts.path, name::String = fts.name)
    ti = time_indices(fts)
    data, _ = load_from_netcdf(; path, var_name = name, grid_size = size(fts)[1:end-1], time_indices_in_memory = ti)

    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    return nothing
end

"""
Return a custom ForcingFromFile forcing, which is a structure keeping 2 FieldTimeSeries
with forcing values and forcing 'lambdas'.
By default both FieldTimeSeries keep only 2 times indices in memory.
"""
function forcing_get_tuple(
    filepath,
    var_name,
    grid,
    time_indices_in_memory,
    backend,
    oceananigans_fieldname,
    data_location,
)
    field_name = oceananigans_fieldname[Symbol(var_name)]
    LX, LY, LZ = data_location[field_name]
    grid_size_tupled = size.(nodes(grid, (LX(), LY(), LZ())))
    grid_size = Tuple(x[1] for x in grid_size_tupled)

    data, times = load_from_netcdf(; path = filepath, var_name, grid_size, time_indices_in_memory)
    dataλ, timesλ =
        load_from_netcdf(; path = filepath, var_name = var_name * "_lambda", grid_size, time_indices_in_memory)

    fts = FieldTimeSeries{LX,LY,LZ}(grid, times; backend, time_indexing = Cyclical(), path = filepath, name = var_name)
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    ftsλ = FieldTimeSeries{LX,LY,LZ}(
        grid,
        timesλ;
        backend,
        time_indexing = Cyclical(),
        path = filepath,
        name = var_name * "_lambda",
    )
    copyto!(interior(ftsλ, :, :, :, :), dataλ)
    fill_halo_regions!(ftsλ)

    _forcing = ForcingFromFile(fts, ftsλ, field_name)
    result = NamedTuple{(Symbol(var_name),)}((_forcing,))
    return result
end

"""
Return a named tuple of forcing functions for all available variables in a NetCDF file.
"""
function forcing_from_file(; grid, filepath, tracers)
    ds = NCDataset(filepath)
    grid.underlying_grid.Nx == ds.dim["Nx"] &&
        grid.underlying_grid.Ny == ds.dim["Ny"] &&
        grid.underlying_grid.Nz == ds.dim["Nz"] ||
        throw(DimensionMismatch("forcing file dimensions not equal to grid dimensions"))
    forcing_variables_names = (map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)
    close(ds)

    oceananigans_fieldname = get_oceananigans_fieldname(forcing_variables_names)
    data_location = get_data_location_dict(forcing_variables_names)

    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    result = mapreduce(
        var_name -> forcing_get_tuple(
            filepath,
            var_name,
            grid,
            time_indices_in_memory,
            backend,
            oceananigans_fieldname,
            data_location,
        ),
        merge,
        forcing_variables_names,
    )

    return result
end

end # module