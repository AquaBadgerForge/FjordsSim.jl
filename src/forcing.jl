using Oceananigans.OutputReaders: FieldTimeSeries, Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: interior
using Oceananigans.Forcings: Forcing
using Oceananigans.Grids: Center, Face, nodes
using Oceananigans.Units: hours
using Oceananigans.Operators: Ax, Ay, Az, volume
using Dates: DateTime, Year, Second
using Statistics: mean
using NCDatasets: Dataset
using Adapt

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend

""" Custom backend for FieldTimeSeries """
struct NetCDFBackend <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
end

NetCDFBackend(length) = NetCDFBackend(1, length)
new_backend(::NetCDFBackend, start, length) = NetCDFBackend(start, length)

Base.length(backend::NetCDFBackend) = backend.length
Base.summary(backend::NetCDFBackend) = string("JLD2Backend(", backend.start, ", ", backend.length, ")")

const NetCDFFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:NetCDFBackend}

# Variable names for restoreable data (to be used in CUDA kernels)
struct Temperature end
struct Salinity end
struct UVelocity end
struct VVelocity end
struct Contaminant end
struct OXYDEP_NUT end
struct OXYDEP_PHY end
struct OXYDEP_HET end
struct OXYDEP_O₂ end
struct OXYDEP_DOM end

const oceananigans_fieldname = Dict(
    :T => Temperature(),
    :S => Salinity(),
    :u => UVelocity(),
    :v => VVelocity(),
    :C => Contaminant(),
    :NUT => OXYDEP_NUT(),
    :P => OXYDEP_PHY(),
    :HET => OXYDEP_HET(),
    :O₂ => OXYDEP_O₂(),
    :DOM => OXYDEP_DOM(),
)

const DATA_LOCATION = Dict(
    Temperature() => (Center, Center, Center),
    Salinity() => (Center, Center, Center),
    UVelocity() => (Face, Center, Center),
    VVelocity() => (Center, Face, Center),
    Contaminant() => (Center, Center, Center),
    OXYDEP_NUT() => (Center, Center, Center),
    OXYDEP_PHY() => (Center, Center, Center),
    OXYDEP_HET() => (Center, Center, Center),
    OXYDEP_O₂() => (Center, Center, Center),
    OXYDEP_DOM() => (Center, Center, Center),
)

# fields[:VARIABLE] doesn't work in CUDA kernels
@inline Base.getindex(fields, i, j, k, ::Temperature) = @inbounds fields.T[i, j, k]
@inline Base.getindex(fields, i, j, k, ::Salinity) = @inbounds fields.S[i, j, k]
@inline Base.getindex(fields, i, j, k, ::UVelocity) = @inbounds fields.u[i, j, k]
@inline Base.getindex(fields, i, j, k, ::VVelocity) = @inbounds fields.v[i, j, k]
@inline Base.getindex(fields, i, j, k, ::Contaminant) = @inbounds fields.C[i, j, k]
@inline Base.getindex(fields, i, j, k, ::OXYDEP_NUT) = @inbounds fields.NUT[i, j, k]
@inline Base.getindex(fields, i, j, k, ::OXYDEP_PHY) = @inbounds fields.P[i, j, k]
@inline Base.getindex(fields, i, j, k, ::OXYDEP_HET) = @inbounds fields.HET[i, j, k]
@inline Base.getindex(fields, i, j, k, ::OXYDEP_O₂) = @inbounds fields.O₂[i, j, k]
@inline Base.getindex(fields, i, j, k, ::OXYDEP_DOM) = @inbounds fields.DOM[i, j, k]

Base.summary(::Temperature) = "temperature"
Base.summary(::Salinity) = "salinity"
Base.summary(::Contaminant) = "contaminant"
Base.summary(::UVelocity) = "u_velocity"
Base.summary(::VVelocity) = "v_velocity"

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
Result will be added to the tendency contributions (a kernel function).
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
    times = []
    for native_time in native_times
        time = native_time - start_time
        time = Second(time).value
        push!(times, time)
    end
    return times
end

function load_from_netcdf(; path::String, var_name::String, grid_size::Tuple, time_indices_in_memory::Tuple)
    ds = Dataset(path)
    var = ds[var_name]
    native_times = ds["time"]

    data = zeros(Float64, (grid_size[1:end]..., length(time_indices_in_memory)))
    j = 1
    for i in time_indices_in_memory
        data[:, :, :, j] .= var[:, :, :, i]
        j += 1
    end
    times = convert.(Int64, native_times_to_seconds(native_times))

    close(ds)
    return data, times
end

""" Updates data in the FieldTimeSeries, e.g. in ForcingFromFile forcing structures. """
function set!(fts::NetCDFFTS, path::String = fts.path, name::String = fts.name)
    ti = time_indices(fts)
    data, _ = load_from_netcdf(; path, var_name = name, grid_size = size(fts)[1:end-1], time_indices_in_memory = ti)

    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    return nothing
end

"""
Constructs a custom ForcingFromFile forcing, which is a structure keeping 2 FieldTimeSeries
with forcing values and forcing 'lambdas'.
By default both FieldTimeSeries keep only 2 times indices in memory.
"""
function forcing_get_tuple(filepath, var_name, grid, time_indices_in_memory, backend)
    field_name = oceananigans_fieldname[Symbol(var_name)]
    LX, LY, LZ = DATA_LOCATION[field_name]
    grid_size_tupled = size.(nodes(grid, (LX(), LY(), LZ())))
    grid_size = Tuple(x[1] for x in grid_size_tupled)

    data, times = load_from_netcdf(; path = filepath, var_name, grid_size, time_indices_in_memory)
    dataλ, timesλ =
        load_from_netcdf(; path = filepath, var_name = var_name * "_lambda", grid_size, time_indices_in_memory)
    boundary_conditions = FieldBoundaryConditions(grid, (LX(), LY(), LZ()))

    fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
        path = filepath,
        name = var_name,
    )
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    ftsλ = FieldTimeSeries{LX,LY,LZ}(
        grid,
        timesλ;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
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
Returns a named tuple of forcing functions for all available variables in a NetCDF file.
"""
function forcing_from_file(grid_ref, filepath, tracers)
    grid = grid_ref[]
    ds = Dataset(filepath)
    grid.underlying_grid.Nx == ds.dim["Nx"] &&
        grid.underlying_grid.Ny == ds.dim["Ny"] &&
        grid.underlying_grid.Nz == ds.dim["Nz"] ||
        throw(DimensionMismatch("forcing file dimensions not equal to grid dimensions"))
    forcing_variables_names = (map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)
    close(ds)

    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    result = mapreduce(
        var_name -> forcing_get_tuple(filepath, var_name, grid, time_indices_in_memory, backend),
        merge,
        forcing_variables_names,
    )

    return result
end
