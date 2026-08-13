module Forcing

using Oceananigans
using Oceananigans.Units: Time
using Oceananigans.OutputReaders: Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.Operators: Ax, Ay, Az, volume
using Oceananigans: fill_halo_regions!, nodes, interior
using Oceananigans.Grids: active_cell, peripheral_node, x_domain, y_domain
using Oceananigans.Fields: interpolate
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index, @Const
using Dates: DateTime, Day, Hour, Millisecond, Second
import Dates
using Adapt
using ArchGDAL
using CUDA  # needed for GPU() to be usable, and for CUDA.functional()
using NCDatasets

import Oceananigans: on_architecture
import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend

using ..Configs:
    AbstractForcingConfig,
    AbstractRiverConfig,
    AbstractBoundaryDataConfig,
    FjordConfig,
    bathymetry_path,
    forcing_path,
    forcing_directory,
    river_forcing_path,
    boundary_data_path,
    boundary_data_directory,
    coverage_window,
    domain_grid,
    simulation_grid
using ..Plotting: plot_forcing, plot_boundaries

export forcing_from_file,
    simulation_forcing,
    forcing_date_range,
    prepare_forcing,
    download_forcing,
    interpolation_architecture,
    forcing_time_steps,
    forcing_source_grid,
    forcing_variable_names,
    forcing_monthly_filename,
    ProjectedSourceGrid,
    NorKystConfig,
    add_rivers,
    download_rivers,
    river_locations,
    river_series,
    river_search_radius,
    RiverLocation,
    OF800RiversConfig,
    LATERAL_EDGES,
    validate_open_edge,
    download_boundaries,
    prepare_boundaries,
    boundary_series,
    boundary_time_steps,
    boundary_source_grid,
    boundary_variable_names,
    boundary_date_range,
    boundary_variable_name,
    NorKystBoundariesConfig

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

function data_location(variable::Symbol)
    if variable in (:u,)
        return (Face, Center, Center)
    elseif variable in (:v,)
        return (Center, Face, Center)
    else
        return (Center, Center, Center)
    end
end

# Runtime generation of dictionaries based on forcing_variables_names
function oceananigans_field_names(variables)
    Dict(Symbol(variable) => Val{Symbol(variable)}() for variable in variables)
end

function data_locations(variables)
    Dict(Val{Symbol(variable)}() => data_location(Symbol(variable)) for variable in variables)
end

# Generic Base.getindex for any variable using Val types
# fields[:VARIABLE] doesn't work in CUDA kernels
@inline Base.getindex(fields, i, j, k, ::Val{name}) where {name} = @inbounds getfield(fields, name)[i, j, k]

""" A custom forcing callable structure """
struct ForcingFromFile{FTS,V}
    fts_value::FTS
    fts_λ::FTS
    field_name::V
end

Adapt.adapt_structure(to, p::ForcingFromFile) =
    ForcingFromFile(Adapt.adapt(to, p.fts_value), Adapt.adapt(to, p.fts_λ), Adapt.adapt(to, p.field_name))

on_architecture(to, forcing::ForcingFromFile) = ForcingFromFile(
    on_architecture(to, forcing.fts_value),
    on_architecture(to, forcing.fts_λ),
    on_architecture(to, forcing.field_name),
)

"""
    forcing_term_x_flux(λ, flux, i, j, k, grid, args...)

Convert a face-crossing flux density `flux` (units `[C]·length/time`, e.g. a tracer flux in
`units(C)·m/s`) entering cell `(i, j, k)` through its x-face into a volumetric tendency for
whatever field the caller applies it to, mirroring one face-term of Oceananigans' own divergence
operator `∇·f = (1/V)[δx(Ax·fx) + ...]` — but one-sided (a source, not a two-face divergence).

Only dimensionally valid for a tracer-like field: a velocity tendency needs acceleration units,
which nothing here supplies. `ForcingFromFile` is generic in `field_name`, so this is reachable for
any field whose `λ` happens to fall in `λ > 1` — it is the caller's responsibility to only route a
field here when `flux` really is a flux density for that field. There is currently no producer of
this `λ` regime at all — `prepared_variable` writes zeros and the river forcing writes
`|λ| = 1/relaxation_timescale < 1` — so this function is unreachable in practice; it is kept, named
honestly, for the day a producer of `λ > 1` exists.
"""
function forcing_term_x_flux(λ, flux, i, j, k, grid, args...)
    return flux * Ax(i, j, k, grid, Center(), Center(), Center()) / volume(i, j, k, grid, Center(), Center(), Center())
end

"""
    forcing_term_y_flux(λ, flux, i, j, k, grid, args...)

As `forcing_term_x_flux`, for a flux entering through the y-face (`Ay` in place of `Ax`). Reached
when `λ < -1`.
"""
function forcing_term_y_flux(λ, flux, i, j, k, grid, args...)
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
        forcing_term_x_flux(λ, value, i, j, k, grid, fields[i, j, k, p.field_name]),
        0,
    )
    result += @inbounds ifelse(
        λ < -1 && value > -990,
        forcing_term_y_flux(λ, value, i, j, k, grid, fields[i, j, k, p.field_name]),
        0,
    )
    result += @inbounds ifelse(
        -1 < λ < 1 && value > -990,
        forcing_term_relax(λ, value, i, j, k, grid, fields[i, j, k, p.field_name]),
        0,
    )
    return result
end

regularize_forcing(forcing::ForcingFromFile, field, field_name, model_field_names) = forcing

function native_times_to_seconds(native_times, start_time = native_times[1])
    return [Second(t - start_time).value for t in native_times]
end

function load_from_netcdf(;
    path::String,
    variable_name::String,
    grid_size::Tuple,
    time_indices_in_memory::Tuple,
    reference_date = nothing,
)
    ds = NCDataset(path)
    variable = ds[variable_name]
    native_times = ds["time"]

    data = zeros(Float64, (grid_size[1:end]..., length(time_indices_in_memory)))
    j = 1
    for i in time_indices_in_memory
        @views data[:, :, :, j] .= coalesce.(variable[:, :, :, i], -999.0)
        j += 1
    end
    start_time = isnothing(reference_date) ? native_times[1] : reference_date
    times = convert.(Int64, native_times_to_seconds(native_times, start_time))

    close(ds)
    return data, times
end

""" Update data in the FieldTimeSeries, e.g. in ForcingFromFile forcing structures. """
function set!(fts::NetCDFFTS, path::String = fts.path, name::String = fts.name)
    ti = time_indices(fts)
    data, _ =
        load_from_netcdf(; path, variable_name = name, grid_size = size(fts)[1:end-1], time_indices_in_memory = ti)

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
    variable_name,
    grid,
    time_indices_in_memory,
    backend,
    field_names,
    locations,
    reference_date,
)
    field_name = field_names[Symbol(variable_name)]
    LX, LY, LZ = locations[field_name]
    grid_size_tupled = size.(nodes(grid, (LX(), LY(), LZ())))
    grid_size = Tuple(x[1] for x in grid_size_tupled)

    data, times =
        load_from_netcdf(; path = filepath, variable_name, grid_size, time_indices_in_memory, reference_date)
    dataλ, timesλ = load_from_netcdf(;
        path = filepath,
        variable_name = variable_name * "_lambda",
        grid_size,
        time_indices_in_memory,
        reference_date,
    )

    fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        backend,
        time_indexing = Cyclical(),
        path = filepath,
        name = variable_name,
    )
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    ftsλ = FieldTimeSeries{LX,LY,LZ}(
        grid,
        timesλ;
        backend,
        time_indexing = Cyclical(),
        path = filepath,
        name = variable_name * "_lambda",
    )
    copyto!(interior(ftsλ, :, :, :, :), dataλ)
    fill_halo_regions!(ftsλ)

    _forcing = ForcingFromFile(fts, ftsλ, field_name)
    result = NamedTuple{(Symbol(variable_name),)}((_forcing,))
    return result
end

"""
Return a named tuple of forcing functions for all available variables in a NetCDF file.

`reference_date` is the instant the time axis is zeroed at, `nothing` meaning the file's own
first record. A coupled run passes its simulation config's `start_date`, so the forcing and the
atmosphere agree on what model time zero stands for.
"""
function forcing_from_file(; grid, filepath, tracers, reference_date = nothing)
    ds = NCDataset(filepath)
    grid.underlying_grid.Nx == ds.dim["Nx"] &&
        grid.underlying_grid.Ny == ds.dim["Ny"] &&
        grid.underlying_grid.Nz == ds.dim["Nz"] ||
        throw(DimensionMismatch("forcing file dimensions not equal to grid dimensions"))
    forcing_variables_names = (map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)
    close(ds)

    field_names = oceananigans_field_names(forcing_variables_names)
    locations = data_locations(forcing_variables_names)

    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    result = mapreduce(
        variable_name -> forcing_get_tuple(
            filepath,
            variable_name,
            grid,
            time_indices_in_memory,
            backend,
            field_names,
            locations,
            reference_date,
        ),
        merge,
        forcing_variables_names,
    )

    return result
end

"""
    forcing_from_file(config::AbstractForcingConfig; grid, tracers, reference_date = nothing)

Return the forcing named tuple for the prepared forcing file this setup names, resolved by
`forcing_path`. `config` is positional so the method dispatches on it.
"""
forcing_from_file(config::AbstractForcingConfig; grid, tracers, reference_date = nothing) =
    forcing_from_file(; grid, filepath = forcing_path(config), tracers, reference_date)

"""
    simulation_forcing(config::AbstractForcingConfig, grid, filepath, tracers, reference_date)

The forcing term object `coupled_simulation` consumes, read from `filepath` — the file
`build_simulation` resolved via `simulation_forcing_path`, which may be the rivers-augmented copy
rather than `forcing_path(config)` itself, so this hook takes the path instead of resolving it.

A dispatched hook rather than a hardcoded `forcing_from_file` call, so a forcing source whose
prepared files are not the FjordSim NetCDF layout `forcing_from_file` reads could supply a
completely different reader. Every built-in source shares this one default, since that layout is a
contract fixed by the read side, not a per-source detail — see `AbstractForcingConfig`.
"""
simulation_forcing(config::AbstractForcingConfig, grid, filepath, tracers, reference_date) =
    forcing_from_file(; grid, filepath, tracers, reference_date)

"""
    simulation_forcing(::Nothing, grid, filepath, tracers, reference_date)

No forcing term at all, for a setup whose `forcing_config` is `nothing`: `coupled_simulation`
still receives a splattable `forcing` keyword, just an empty one, which is exactly
`HydrostaticFreeSurfaceModel`'s own default.
"""
simulation_forcing(::Nothing, grid, filepath, tracers, reference_date) = NamedTuple()

"""
    forcing_date_range(config::AbstractForcingConfig, filepath)

First and last date the prepared forcing file at `filepath` covers, as a `(first, last)` tuple, or
`nothing` for a source that cannot report its dates.

What `FjordSim.Simulations.validate_time_coverage` checks the run's window against. Both readers use
`Cyclical()` time indexing, which wraps rather than failing outside its data, so this is what keeps
a run that outlasts its forcing from quietly replaying the beginning.

The supertype default reads `ds["time"]` from the prepared NetCDF, which is the layout every
built-in source writes. A source that overrode `simulation_forcing` to read something else overrides
this too — otherwise its coverage would be validated by a reader it does not use. The exact mirror of
`atmosphere_date_range`, which the forcing side previously had no counterpart to: this was a bare
`forcing_date_range(filepath)` in `Simulations`, dispatched on nothing at all.
"""
forcing_date_range(config::AbstractForcingConfig, filepath) = NCDataset(filepath) do ds
    dates = ds["time"][:]
    (first(dates), last(dates))
end

"""
    forcing_date_range(::Nothing, filepath)

No dates to check, for a setup whose `forcing_config` is `nothing`.
"""
forcing_date_range(::Nothing, filepath) = nothing

# --- Forcing preparation ---

const FORCING_DEFLATE_LEVEL = 5

"""
The four lateral boundaries a regional domain can be open on, in the order every error message
lists them.

Defined here rather than in `BoundaryConditions` because that module is included after this one and
both need it: `validate_open_edge` checks a forcing config's `open_edge` before any regridding, and
every `Val{edge}` dispatch in the boundary conditions falls back to an `ArgumentError` naming this
tuple. There used to be one copy per module, with identical contents, an identical membership test
and an identical error string.
"""
const LATERAL_EDGES = (:south, :north, :west, :east)
# Longest run of missing days interpolated without a separate warning, matching the default of
# `NumericalEarth.DataWrangling.fill_gaps!`.
const FORCING_MAX_GAP = 6

"""
    interpolation_architecture(config)

Where `prepare_forcing` runs its interpolation kernel, from `config.architecture`:

- `:auto`: the GPU when one is usable, else the CPU. The default, so one setup runs unchanged on
  a GPU machine and on a laptop.
- `:gpu`: the GPU, erroring when none is usable rather than silently running ~12x slower.
- `:cpu`: the CPU.

Only the kernel moves; `target_grid` and the masks stay on the CPU because building them walks
`peripheral_node` cell by cell.

Defined for a boundary-data config too, which reads its `architecture` field the same way and runs
the same kernel. `Simulations.simulation_architecture` also reuses the `Val` methods below, which is
why the resolution and the field read are separate methods rather than one function.
"""
interpolation_architecture(config::AbstractForcingConfig) =
    interpolation_architecture(Val(config.architecture))

interpolation_architecture(config::AbstractBoundaryDataConfig) =
    interpolation_architecture(Val(config.architecture))

interpolation_architecture(::Val{:cpu}) = CPU()
interpolation_architecture(::Val{:auto}) = CUDA.functional() ? GPU() : CPU()

function interpolation_architecture(::Val{:gpu})
    CUDA.functional() || error(
        "architecture = :gpu was requested but no usable GPU was found. " *
        "Use :auto to fall back to the CPU, or :cpu to ask for it explicitly.",
    )
    return GPU()
end

interpolation_architecture(::Val{selector}) where {selector} =
    throw(ArgumentError("architecture must be one of (:auto, :cpu, :gpu), got :$selector"))

# --- Extension hooks ---

"""
    forcing_time_steps(config)

Every source time record available to `config`, as `SourceRecord`s sorted by date with
duplicates dropped. `prepare_forcing` completes them to a gap-free daily axis with
`daily_time_steps`.

A new forcing dataset implements this on its `AbstractForcingConfig` subtype; see
`forcing_time_steps(config::NorKystConfig)` in `src/Forcing/norkyst.jl`.
"""
function forcing_time_steps end

"""
    forcing_source_grid(config, filepath)

Geometry of the source data in `filepath`: the coordinates it lives on and the projection
they are defined in. Return a `ProjectedSourceGrid` for any source on a regular grid in
projected meters — that reuses `source_field_grid` and `projected_target_nodes` unchanged.

A new forcing dataset implements this on its `AbstractForcingConfig` subtype.
"""
function forcing_source_grid end

"""
    forcing_variable_names(config)

Mapping from source variable name to the FjordSim forcing name it becomes, e.g.
`"temperature" => "T"`. Only the intersection with `config.parameters` is prepared, so this
also declares which variables the dataset can supply.

A new forcing dataset implements this on its `AbstractForcingConfig` subtype.
"""
function forcing_variable_names end

"""
    download_forcing(config::FjordConfig)
    download_forcing(target_grid, config::AbstractForcingConfig)

Fetch and subset the source data a later `prepare_forcing` call reads, into
`forcing_directory(config)`.

The `FjordConfig` form is the generic driver: it builds the setup's grid on the CPU and
dispatches on the forcing config, so a new dataset only implements the second method. The
grid is passed rather than the grid config because a dataset needs the domain bounds, which
`x_domain`/`y_domain` provide for any `AbstractGridConfig`.
"""
download_forcing(config::FjordConfig) =
    download_forcing(domain_grid(config.grid_config, CPU()), config.forcing_config)

download_forcing(target_grid, ::Nothing) = nothing

"""
    validate_open_edge(edge)

Return `edge` if it names a lateral boundary, else throw. Checked up front so a typo fails
before any regridding rather than deep inside the variable loop.

`LATERAL_EDGES` lives in `BoundaryConditions`, which is included after this module, so the tuple is
stated here and that module imports it — previously each held its own copy under a different name,
with an identical membership test and an identical error string.
"""
function validate_open_edge(edge)
    edge in LATERAL_EDGES ||
        throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))
    return edge
end

"""
    ProjectedSourceGrid

Geometry of a source subset living on a regular grid in projected meters: the coordinates,
its depth levels and the projection they are defined in.

Nothing here is dataset-specific, so any source on a regular projected grid can return one of
these from `forcing_source_grid` and inherit `source_field_grid` and
`projected_target_nodes` unchanged.

# Fields
- `x`, `y`: Projected coordinate centers in meters, regularly spaced.
- `depths`: Depth levels in meters, positive down.
- `proj4`: PROJ.4 definition of the projection.
"""
struct ProjectedSourceGrid
    x::Vector{Float64}
    y::Vector{Float64}
    depths::Vector{Float64}
    proj4::String
end

"""
    SourceFill

Where every masked source cell takes its value from, so that a filled source slab can be handed
to `Oceananigans.Fields.interpolate` without a `NaN` anywhere in it.

Built once per variable from the first time step, because the mask is a property of the dataset
rather than of the step. `nearest` is a linear index into the `(NX, NY)` plane of the same level;
levels with no valid cell at all (source levels deeper than the fjord) are marked in
`level_valid` and filled from the deepest level that does have data.

This is not `NumericalEarth`'s `inpaint_mask!`: that re-`NaN`s every masked cell before
propagating horizontally, discarding its own `continue_downwards!` pass, so a fully masked level
has nothing to spread from. With the default `NearestNeighborInpainting(Inf)` it then never
terminates, and with a finite cap `_fill_nans!` turns the level into zeros — a plausible-looking
fabricated temperature or salinity.

# Fields
- `valid`: Whether each source cell holds data.
- `nearest`: `(NX, NY, levels)` linear index of the nearest valid cell in the same level.
- `level_valid`: Whether a level has any valid cell inside the subset.
"""
struct SourceFill
    valid::BitArray{3}
    nearest::Array{Int32,3}
    level_valid::Vector{Bool}
end

"""
    SourceRecord

One time record present in a downloaded source file.

# Fields
- `date`: Time of the record.
- `filepath`: Source file holding it.
- `index`: Time index within that file.
"""
struct SourceRecord
    date::DateTime
    filepath::String
    index::Int
end

"""
    ForcingTimeStep

One output time step and the source records it is built from.

A step present in the download reads a single record (`lower === upper`, `weight == 0`). A day
missing from the download is the linear interpolation of the records bracketing it, following
`NumericalEarth.DataWrangling.fill_gaps!`.

# Fields
- `date`: Time of the output step.
- `lower`, `upper`: Bracketing source records; equal when the step is present in the download.
- `weight`: Weight of `upper`; zero when the step is present in the download.
"""
struct ForcingTimeStep
    date::DateTime
    lower::SourceRecord
    upper::SourceRecord
    weight::Float32
end

"""
    PreparedVariable

Everything needed to write one forcing variable: where it comes from, the land mask and target
nodes of its Oceananigans location, how to fill its source mask, and its relaxation lambdas.

`x` and `y` are the target nodes projected into the source coordinate system, held as 2D arrays
so the interpolation kernel only indexes them.

Parametric in the number of NetCDF dimension names, which is the one thing an interior forcing
variable and an open-boundary one disagree about: four for the forcing file's
`(along, across, Nz, time)`, three or two for the boundary file, whose across-edge dimension is one
cell and whose surface variables have no depth axis. Everything else — the three-dimensional mask
with its one-cell slab, the projected nodes, the source fill — is shared, which is why
`prepare_boundaries` reuses this type rather than defining its own.
"""
struct PreparedVariable{N}
    source_name::String
    name::String
    dimensions::NTuple{N,String}
    mask::Array{Bool,3}
    x::Matrix{Float64}
    y::Matrix{Float64}
    z::Vector{Float64}
    mask_fill::SourceFill
    lambda::Array{Float32,3}
end

# A hand-written constructor because `Base.@kwdef`-style automatic conversion does not happen for a
# *parametric* struct: the generated constructor keeps the declared field types in its signature, so
# the `BitArray` `water_mask` returns would not match `Array{Bool,3}`.
PreparedVariable(source_name, name, dimensions::NTuple{N,String}, mask, x, y, z, mask_fill, lambda) where {N} =
    PreparedVariable{N}(
        String(source_name),
        String(name),
        dimensions,
        Array{Bool,3}(mask),
        Matrix{Float64}(x),
        Matrix{Float64}(y),
        Vector{Float64}(z),
        mask_fill,
        Array{Float32,3}(lambda),
    )

"""
    prepare_forcing(target_grid, config::AbstractForcingConfig)

Regrid the source files already downloaded into `forcing_directory(config)` onto `target_grid`,
and write a forcing NetCDF at `forcing_path(config)` that `forcing_from_file` reads directly.

The dataset enters only through the three hooks `forcing_variable_names`,
`forcing_time_steps` and `forcing_source_grid`; everything after them is the same for every
source.

`target_grid` should be the `ImmersedBoundaryGrid` the simulation runs on, so the output
dimensions and land mask match it exactly — the mask comes from `peripheral_node`, the same
predicate the model uses. Interpolation is one trilinear `Oceananigans.Fields.interpolate` call
per target cell, taken in the source projection and clamped at the vertical ends. Cells the model
treats as land are written as `NaN`, which `forcing_from_file` reads back as the `value > -990`
sentinel the forcing kernel skips. Days missing from the download are interpolated between their
neighbours.

`config.architecture` selects where the interpolation kernel runs, via
`interpolation_architecture`, and is independent of `target_grid`, which must stay on the CPU:
building the masks walks `peripheral_node` cell by cell, which needs scalar access to the
bathymetry. Only the target nodes, mask and output buffer are moved to the interpolation device,
and only the finished slab comes back for each write, so a GPU run keeps the same streaming
memory profile.

Relaxation lambdas are written as **zero everywhere**. The interior sponge band this file used to
carry along `open_edge` is gone: the open lateral boundary now nudges towards hourly exterior data at
the boundary itself (`NormalRadiation`, `GravityWaveRadiation`), and a band relaxing the same
variables a few cells inside would fight it. `add_rivers` writes the only nonzero lambdas the file
ever carries, so the `*_lambda` variables exist for that step and for
`forcing_from_file`'s own contract, and the prepared values are otherwise read only as initial
conditions (`FromForcing`).

# Why the interpolation is not `NumericalEarth`'s

A source grid in projected meters — NorKyst's is rotated about 59 degrees from east in the
Oslofjord region — cannot be a `LatitudeLongitudeGrid`; treating its 2D `lon`/`lat` as separable
is wrong by roughly 100 grid cells. `NumericalEarth`'s dataset path (`native_grid`,
`set!(field, metadata)`, `DatasetRestoring`) always builds a `LatitudeLongitudeGrid`, and
`Oceananigans`' field-to-field `interpolate!` cannot bridge lon/lat and projected meters, so
neither applies here. The NORA3 atmosphere avoids this only because its file has already been
regridded to a regular lon/lat grid upstream. What does work is projecting the target nodes and
interpolating against a `RectilinearGrid` in projected meters; see `source_field_grid` and
`interpolate_to_target!`.

# Coverage
`coverage` is the `(first, last)` calendar interval the run needs, from `coverage_window`. The time
axis is padded to reach both ends by replicating the nearest step; `nothing` prepares exactly the
downloaded range. See `pad_time_steps`.

# Returns
A named tuple with `output_file`, `times` and `variables` (the written variable names).
"""
function prepare_forcing(target_grid, config::AbstractForcingConfig; coverage = nothing)
    architecture = interpolation_architecture(config)
    variable_names = forcing_variable_names(config)
    source_names = [name for name in config.parameters if haskey(variable_names, name)]
    isempty(source_names) && error(
        "None of the configured parameters $(config.parameters) map to a forcing variable. " *
        "Known $(nameof(typeof(config))) variables: $(sort(collect(keys(variable_names)))).",
    )

    validate_open_edge(config.open_edge)

    steps = pad_time_steps(daily_time_steps(forcing_time_steps(config)), coverage)
    reference_file = first(steps).lower.filepath
    source = forcing_source_grid(config, reference_file)
    @info "Preparing forcing from $(length(unique(step -> step.lower.filepath, steps))) file(s), " *
          "$(length(steps)) time steps: $(first(steps).date) to $(last(steps).date)"

    variables = [prepared_variable(name, target_grid, source, reference_file, config) for name in source_names]

    output_file = forcing_path(config)
    @info "Writing forcing file to $output_file, interpolating on $(summary(architecture))"
    write_forcing_file(output_file, target_grid, variables, steps, source, architecture)
    @info "Finished preparing forcing"

    return (; output_file, times = [step.date for step in steps], variables = [variable.name for variable in variables])
end

prepare_forcing(target_grid, ::Nothing; coverage = nothing) = nothing

"""
    prepare_forcing(config::FjordConfig)

Regrid the forcing a whole setup names onto its simulation grid, and write the diagnostic plot.
Returns `nothing` when the setup names no forcing.

This is the setup-level driver, the same shape as `download_forcing(config::FjordConfig)`. The
grid comes from the processed bathymetry rather than the grid config, so the output's land mask
matches the model exactly, which means `prepare_bathymetry` must have run first. It stays on the
CPU regardless of `config.forcing_config.architecture` — that field selects where the
interpolation kernel runs, and building the masks needs scalar access to the bathymetry.

The coverage window comes from the setup's simulation config, so the prepared file spans the run
the setup describes; a setup naming no simulation config gets `nothing` and the downloaded range.

# Returns
The `prepare_forcing(target_grid, config)` named tuple with `plot_file` added.
"""
function prepare_forcing(config::FjordConfig)
    isnothing(config.forcing_config) && return nothing

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    grid = simulation_grid(config.grid_config, bathymetry_file, CPU())
    result = prepare_forcing(
        grid,
        config.forcing_config;
        coverage = coverage_window(config.simulation_config),
    )
    plot_file = plot_forcing(grid, config.forcing_config)

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing file saved to $(result.output_file)"
    @info "Forcing plot saved to $plot_file"

    return (; result..., plot_file)
end

"""
    uniform_time_steps(records, spacing; max_gap, unit)

Complete `records` to a gap-free axis at `spacing`, linearly interpolating any missing step between
the records bracketing it. Downloaded months are occasionally short a record, which would otherwise
leave a hole in the cyclical forcing period.

A gap longer than `max_gap` steps is still interpolated, but warned about separately: MET's
NorKyst archive is missing whole weeks in 2017 and 2018, and papering over those silently would
misrepresent the forcing. Returned as single-record steps if the times are not on a whole-`spacing`
cadence.

`unit` only names the spacing in the warnings. Written generically because the interior forcing is
daily and the open-boundary data is hourly, and the gap-filling argument is identical for both.
"""
function uniform_time_steps(records, spacing::Dates.Period; max_gap, unit)
    single(record) = ForcingTimeStep(record.date, record, record, 0)

    differences = diff([record.date for record in records])
    if !all(difference -> difference % spacing == Millisecond(0), differences)
        @warn "Source times are not on a whole-$unit cadence; using them as downloaded."
        return map(single, records)
    end

    steps = ForcingTimeStep[]
    interpolated = DateTime[]
    step_milliseconds = Dates.value(Millisecond(spacing))

    for n = 1:length(records)-1
        lower = records[n]
        upper = records[n+1]
        push!(steps, single(lower))

        span = Dates.value(upper.date - lower.date) ÷ step_milliseconds
        span > max_gap + 1 && @warn "Source gap of $(span - 1) $unit(s) after $(lower.date) " *
                                    "exceeds max_gap = $max_gap; interpolating anyway"

        for missing_step = 1:span-1
            date = lower.date + missing_step * spacing
            push!(interpolated, date)
            push!(steps, ForcingTimeStep(date, lower, upper, missing_step / span))
        end
    end
    push!(steps, single(last(records)))

    isempty(interpolated) ||
        @warn "Interpolated $(length(interpolated)) missing $unit(s): $(join(interpolated, ", "))"

    return steps
end

"""
    daily_time_steps(records; max_gap = FORCING_MAX_GAP)

`uniform_time_steps` at a one-day spacing: the cadence of the NorKyst daily-average interior
forcing.
"""
daily_time_steps(records; max_gap = FORCING_MAX_GAP) =
    uniform_time_steps(records, Day(1); max_gap, unit = "day")

"""
    pad_time_steps(steps, coverage)

Extend `steps` to span `coverage`, replicating the nearest step at each end that falls short. A
`nothing` coverage leaves `steps` alone.

Written after `daily_time_steps` so its whole-day cadence check sees only the downloaded times: a
pad step is deliberately off that cadence whenever the run's window is.

The pad carries the *same* source records and blend weight as the step it copies, so the written
field is identical to its neighbour rather than a re-derived approximation of it. That is what makes
the padding a statement about the time axis alone.

The reason this exists is that a run's window rarely lines up with a dataset's records — NorKyst's
first daily record here is at 12:00, so a run starting at 00:00 on the same day has no forcing for
its first half-day. `forcing_from_file` reads with `Cyclical()` time indexing, which does not fail
outside its data but wraps to the far end, so the shortfall would be filled with the following
December rather than reported. Padding puts a real record at both ends instead, which is what lets
`Simulations.validate_time_coverage` stay a hard check.

That last point is also why a pad may reach **at most one record spacing** past the downloaded axis.
Unbounded, this would defeat the very check it exists to serve: a window running a year past the
data would be "covered" by one replicated December, `forcing_date_range` would report the invented
span, `validate_time_coverage` would pass, and the run would interpolate twelve months between two
identical records — quietly worse than the wrap. Overshooting by more than a spacing is a missing
download, and is reported as one.
"""
pad_time_steps(steps, ::Nothing) = steps

function pad_time_steps(steps, coverage)
    first_needed, last_needed = coverage
    padded = collect(steps)

    if first(padded).date > first_needed
        head = first(padded)
        spacing = forcing_record_spacing(padded)
        head.date - first_needed <= spacing || error(
            "The run starts at $first_needed but the downloaded forcing only begins at " *
            "$(head.date), more than one record spacing " *
            "($(Dates.canonicalize(spacing))) later. Download the preceding period, or move the " *
            "simulation config's `start_date` to " *
            "$(head.date - spacing) or later.",
        )
        pushfirst!(padded, ForcingTimeStep(first_needed, head.lower, head.upper, head.weight))
        @info "Padded forcing start: replicated $(head.date) at $first_needed"
    end

    if last(padded).date < last_needed
        tail = last(padded)
        spacing = forcing_record_spacing(padded)
        last_needed - tail.date <= spacing || error(
            "The run ends at $last_needed but the downloaded forcing already stops at " *
            "$(tail.date), more than one record spacing " *
            "($(Dates.canonicalize(spacing))) earlier. Download the following period, or shorten " *
            "`stop_time` so the run ends by " *
            "$(tail.date + spacing).",
        )
        push!(padded, ForcingTimeStep(last_needed, tail.lower, tail.upper, tail.weight))
        @info "Padded forcing end: replicated $(tail.date) at $last_needed"
    end

    return padded
end

"""
    forcing_record_spacing(steps)

The axis's record spacing, taken from its first pair — `daily_time_steps` has already made it
uniform, so one pair speaks for the whole axis. A single-record axis has no spacing to bound a pad
by, and is an error rather than a guess.
"""
function forcing_record_spacing(steps)
    length(steps) >= 2 || error(
        "Cannot pad a forcing axis of $(length(steps)) record(s): its record spacing, which " *
        "bounds how far a pad may reach, is unknown.",
    )
    return steps[2].date - steps[1].date
end

"""
    water_mask(target_grid, LX, LY, edge)

Water flag per cell of `target_grid` at the location `(LX, LY, Center)`, taken from
`Oceananigans.Grids.peripheral_node` so that it is by construction the same predicate the model
uses. That covers the staggering: a velocity face is peripheral when either tracer cell it
separates is inactive, and `PartialCellBottom` decides "inactive" with its own
`minimum_fractional_cell_height` rather than a hand-rolled cell-center test.

`edge` names the open lateral boundary, whose outermost face row `open_boundary_water!` restores.
"""
function water_mask(target_grid, LX, LY, edge)
    Nx, Ny, Nz = size(target_grid)
    nx = LX === Face ? Nx + 1 : Nx
    ny = LY === Face ? Ny + 1 : Ny

    mask = falses(nx, ny, Nz)
    for k = 1:Nz, j = 1:ny, i = 1:nx
        mask[i, j, k] = !peripheral_node(i, j, k, target_grid, LX(), LY(), Center())
    end

    return open_boundary_water!(mask, target_grid, LX, LY, Val(validate_open_edge(edge)))
end

"""
    open_boundary_water!(mask, target_grid, LX, LY, ::Val{edge})

Mark the velocity faces lying on the open lateral boundary `edge` as water wherever the tracer
cell just inside them is active.

`peripheral_node` treats those faces as peripheral because the tracer cell outside the domain is
a halo cell, and halo cells are inactive in `Bounded` directions. For a closed wall that is
right, but `edge` is where the domain is open and the normal velocity carries a
`NormalFlowBoundaryCondition`, so leaving that row masked would drop the data exactly where it is
needed — both for the prepared forcing's initial conditions and for the boundary file. Only the
component staggered across `edge` is affected, which is why each method tests one location.

One method per edge rather than a four-branch `if`: the edge is a `Val` everywhere it decides
anything in this package.
"""
open_boundary_water!(mask, target_grid, LX, LY, ::Val{edge}) where {edge} = mask

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:south})
    LY === Face || return mask
    for k = 1:size(target_grid, 3), i in axes(mask, 1)
        mask[i, 1, k] = active_cell(i, 1, k, target_grid)
    end
    return mask
end

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:north})
    LY === Face || return mask
    Ny = size(target_grid, 2)
    for k = 1:size(target_grid, 3), i in axes(mask, 1)
        mask[i, Ny+1, k] = active_cell(i, Ny, k, target_grid)
    end
    return mask
end

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:west})
    LX === Face || return mask
    for k = 1:size(target_grid, 3), j in axes(mask, 2)
        mask[1, j, k] = active_cell(1, j, k, target_grid)
    end
    return mask
end

function open_boundary_water!(mask, target_grid, LX, LY, ::Val{:east})
    LX === Face || return mask
    Nx = size(target_grid, 1)
    for k = 1:size(target_grid, 3), j in axes(mask, 2)
        mask[Nx+1, j, k] = active_cell(Nx, j, k, target_grid)
    end
    return mask
end

"""
    forcing_dimension_names(name)

NetCDF dimension names for a forcing variable, staggered according to its Oceananigans
location. Written in Julia (fastest-first) order, which yields the `(time, Nz, Ny, Nx)` layout
in the file.
"""
function forcing_dimension_names(name)
    LX, LY, _ = data_location(Symbol(name))
    return (LX === Face ? "Nx_faces" : "Nx", LY === Face ? "Ny_faces" : "Ny", "Nz", "time")
end

"""
    prepared_variable(source_name, target_grid, source, filepath, config)

Build the mask, projected target nodes, source fill and lambdas for one source variable at its
target location.
"""
function prepared_variable(source_name, target_grid, source, filepath, config::AbstractForcingConfig)
    name = forcing_variable_names(config)[source_name]
    LX, LY, LZ = data_location(Symbol(name))
    @info "Preparing target nodes and source fill for $source_name -> $name"

    mask = water_mask(target_grid, LX, LY, config.open_edge)
    longitude = Array(λnodes(target_grid, LX()))
    latitude = Array(φnodes(target_grid, LY()))
    size(mask)[1:2] == (length(longitude), length(latitude)) || error(
        "Target node count $(length(longitude))x$(length(latitude)) does not match the mask " *
        "$(size(mask)[1:2]) for $name; the grid topology must be bounded in both directions.",
    )

    x, y = projected_target_nodes(longitude, latitude, source)
    shape = (length(longitude), length(latitude))
    mask_fill = source_fill(source_validity(filepath, source_name))
    # Zero everywhere: the interior relaxation band is gone, and `add_rivers` writes the only
    # nonzero lambdas this file ever carries. See `prepare_forcing`.
    lambda = zeros(Float32, size(mask))

    return PreparedVariable(
        source_name,
        name,
        forcing_dimension_names(name),
        mask,
        reshape(x, shape),
        reshape(y, shape),
        Array(znodes(target_grid, LZ())),
        mask_fill,
        lambda,
    )
end

"""
    projected_target_nodes(longitude, latitude, source::ProjectedSourceGrid)

Project the target node grid into the source projection, flattened with longitude varying
fastest. The whole point cloud is transformed in a single GDAL call, as in
`FjordSim.Bathymetry.sample_bathymetry_points!`.
"""
function projected_target_nodes(longitude, latitude, source::ProjectedSourceGrid)
    point_count = length(longitude) * length(latitude)
    x = Vector{Float64}(undef, point_count)
    y = Vector{Float64}(undef, point_count)

    point = 1
    for φ in latitude, λ in longitude
        x[point] = λ
        y[point] = φ
        point += 1
    end

    ArchGDAL.importEPSG(4326; order = :trad) do source_srs
        ArchGDAL.importPROJ4(source.proj4) do target_srs
            ArchGDAL.createcoordtrans(source_srs, target_srs) do transform
                ArchGDAL.transform!(x, y, zeros(Float64, point_count), transform)
            end
        end
    end

    return x, y
end

"""
    source_validity(filepath, source_name)

Whether each source cell holds data, taken from the first time step. Source land cells and
the points the download masked out as being outside the requested box are `NaN`, and that
pattern is a property of the dataset rather than of the time step, so one read serves every
step.

Always `(NX, NY, levels)`: a surface variable with no depth axis (`zeta`, `ubar`, `vbar`, which the
boundary pipeline prepares) reads as a plane and is reshaped to one level, so `SourceFill`,
`fill_source!` and `set_source_field!` stay three-dimensional for every source variable.
"""
function source_validity(filepath, source_name)
    return NCDataset(filepath) do ds
        return isfinite.(as_source_slab(coalesce.(read_last_dimension(ds[source_name], 1), NaN32)))
    end
end

"""
    read_last_dimension(variable, index)

`variable` at `index` along its last (time) dimension, with `:` for every dimension before it —
`(X, Y, depth)` for a 3D source variable, `(X, Y)` for a surface one. The trailing index cannot be
written literally because the two ranks differ.
"""
read_last_dimension(variable, index) =
    variable[ntuple(_ -> Colon(), ndims(variable) - 1)..., index]

"""
    as_source_slab(data)

`data` as a `(NX, NY, levels)` array: a 3D slab unchanged, a 2D surface plane reshaped to a single
level. Everything downstream of the source read — `SourceFill`, `fill_source!`,
`set_source_field!`, the interpolation kernel — is written for three dimensions, and one reshape
here is cheaper than a rank parameter threaded through all of them.
"""
as_source_slab(data::AbstractArray{<:Any,3}) = data
as_source_slab(data::AbstractArray{<:Any,2}) = reshape(data, size(data, 1), size(data, 2), 1)

"""
    source_fill(valid)

Build the `SourceFill` for a source mask: the nearest valid cell per masked cell per level.
"""
function source_fill(valid)
    NX, NY, ND = size(valid)
    level_valid = [any(@view valid[:, :, level]) for level = 1:ND]
    any(level_valid) || error("Source data has no valid cell inside the requested region.")

    nearest = zeros(Int32, NX, NY, ND)
    for level = 1:ND
        level_valid[level] || continue
        nearest[:, :, level] = nearest_valid_map(valid, level)
    end

    return SourceFill(BitArray(valid), nearest, level_valid)
end

"""
    fill_source!(slab, mask_fill::SourceFill)

Replace every masked value in `slab` with the nearest valid value in the same level, and every
fully masked level with the deepest level that has data. Mutates and returns `slab`, which is then
free of `NaN` so `interpolate` can never propagate one.

Filling a level from the deepest level with data reproduces what the previous hand-written
vertical interpolation achieved by clamping: source levels below the fjord's own bathymetry carry
no information, so the deepest real value is the best available.
"""
function fill_source!(slab, mask_fill::SourceFill)
    NX, NY, ND = size(slab)
    size(mask_fill.valid) == (NX, NY, ND) || throw(DimensionMismatch("source fill does not match the slab"))
    flat = reshape(slab, NX * NY, ND)

    for level = 1:ND
        mask_fill.level_valid[level] || continue
        for j = 1:NY, i = 1:NX
            @inbounds mask_fill.valid[i, j, level] || (flat[i+(j-1)*NX, level] = flat[mask_fill.nearest[i, j, level], level])
        end
    end

    deepest = findlast(mask_fill.level_valid)
    for level = 1:ND
        mask_fill.level_valid[level] && continue
        @inbounds flat[:, level] .= @view flat[:, deepest]
    end

    return slab
end

"""
    nearest_valid_map(valid, level)

For every cell of source `level`, the linear index of the nearest valid cell in that level.

Computed by relaxing squared Euclidean distance out from all valid cells at once, propagating
the originating cell over 8 neighbors, which costs one near-linear pass per level. This plays the
role of the `scipy.ndimage.distance_transform_edt` nearest-fill used to prepare the Oslofjord
forcing, and of `NumericalEarth`'s `inpaint_mask!`, which cannot be used here — see `SourceFill`.
"""
function nearest_valid_map(valid, level)
    NX, NY, _ = size(valid)
    nearest = zeros(Int32, NX, NY)
    origin_i = zeros(Int32, NX, NY)
    origin_j = zeros(Int32, NX, NY)
    distance = fill(typemax(Int32), NX, NY)
    queue = Tuple{Int32,Int32}[]

    for j = 1:NY, i = 1:NX
        valid[i, j, level] || continue
        nearest[i, j] = i + (j - 1) * NX
        origin_i[i, j] = i
        origin_j[i, j] = j
        distance[i, j] = 0
        push!(queue, (i, j))
    end

    head = 1
    while head <= length(queue)
        i, j = queue[head]
        head += 1
        source_i = origin_i[i, j]
        source_j = origin_j[i, j]

        for dj = -1:1, di = -1:1
            di == 0 && dj == 0 && continue
            neighbor_i = i + di
            neighbor_j = j + dj
            1 <= neighbor_i <= NX && 1 <= neighbor_j <= NY || continue

            candidate = (neighbor_i - source_i)^2 + (neighbor_j - source_j)^2
            candidate < distance[neighbor_i, neighbor_j] || continue

            distance[neighbor_i, neighbor_j] = candidate
            origin_i[neighbor_i, neighbor_j] = source_i
            origin_j[neighbor_i, neighbor_j] = source_j
            nearest[neighbor_i, neighbor_j] = source_i + (source_j - 1) * NX
            push!(queue, (neighbor_i, neighbor_j))
        end
    end

    return nearest
end

"""
    solve_vertical_faces(depths)

Vertical faces of an Oceananigans grid whose cell centers land exactly on the source depth
levels, so that `interpolate` sees the levels the data is actually defined at.

Oceananigans grids are specified by faces with centers at face midpoints, while a source's
`depth` values *are* centers of a strongly non-uniform axis. Solving `c_k = (f_k + f_{k+1})/2`
leaves one free parameter `s = f_1`, with `f_k = (-1)^(k-1) s + g_k`, and a positive cell
thickness requires `f_k < c_k`. Each level therefore bounds `s` from one side, alternating with
parity. For NorKyst's 16 levels the resulting interval is about one metre wide in a 3 km domain —
it exists, but it is a property of MET's depth list rather than of this code, so an empty
interval is an error rather than a silently misplaced grid.
"""
function solve_vertical_faces(depths)
    centers = -reverse(depths)
    n = length(centers)

    g = zeros(n + 1)
    for k = 2:n+1
        g[k] = 2 * centers[k-1] - g[k-1]
    end

    lower = -Inf
    upper = Inf
    for k = 1:n
        if isodd(k)
            upper = min(upper, centers[k] - g[k])
        else
            lower = max(lower, g[k] - centers[k])
        end
    end

    upper > lower || error(
        "Source depth levels $(-centers) cannot be represented as an Oceananigans grid: the " *
        "faces placing cell centers there are non-monotonic for every choice of the deepest " *
        "face (feasible interval would be ($lower, $upper)).",
    )

    faces = zeros(n + 1)
    faces[1] = (lower + upper) / 2
    for k = 1:n
        faces[k+1] = 2 * centers[k] - faces[k]
    end

    return faces
end

"""
    source_field_grid(source::ProjectedSourceGrid, architecture = CPU())

The source subset as a `RectilinearGrid` in projected meters, on `architecture`.

The projected coordinates are regular, so this is a legal `RectilinearGrid` and
`Oceananigans.Fields.interpolate` applies. It is *not* expressible as a `LatitudeLongitudeGrid`:
NorKyst's grid is rotated about 59 degrees from east in the Oslofjord region, which also rules
out the `NumericalEarth` dataset path, whose `native_grid` is always a `LatitudeLongitudeGrid`.

A source on a differently-shaped grid overloads this, and `projected_target_nodes`, on its own
source-grid type.
"""
function source_field_grid(source::ProjectedSourceGrid, architecture = CPU())
    Δx = source.x[2] - source.x[1]
    Δy = source.y[2] - source.y[1]

    return RectilinearGrid(
        architecture;
        size = (length(source.x), length(source.y), length(source.depths)),
        x = collect(range(source.x[1] - Δx / 2, step = Δx, length = length(source.x) + 1)),
        y = collect(range(source.y[1] - Δy / 2, step = Δy, length = length(source.y) + 1)),
        z = solve_vertical_faces(source.depths),
        topology = (Bounded, Bounded, Bounded),
    )
end

"""
    set_source_field!(source_field, slab, mask_fill::SourceFill)

Load one source slab into `source_field`, filling its mask via `fill_source!` and flipping the
depth axis so it increases upwards like the grid. `slab` is mutated.
"""
function set_source_field!(source_field, slab, mask_fill::SourceFill)
    fill_source!(slab, mask_fill)
    set!(source_field, reverse(slab, dims = 3))
    fill_halo_regions!(source_field)

    return source_field
end

@kernel function _interpolate_to_target!(output, @Const(x), @Const(y), @Const(z), @Const(mask),
                                         source_field, source_grid, source_location)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        node = (x[i, j], y[i, j], z[k])
        output[i, j, k] = ifelse(
            mask[i, j, k],
            convert(Float32, interpolate(node, source_field, source_location, source_grid)),
            NaN32,
        )
    end
end

"""
    interpolate_to_target!(output, source_field, x, y, z, mask, architecture)

Interpolate `source_field` onto the target nodes `x`, `y`, `z`, writing `NaN` where `mask` says the
target cell is land. Every array must live on `architecture`.

The target nodes were projected into the source coordinate system by `projected_target_nodes`,
which is what lets `Oceananigans.Fields.interpolate` do the work: the source is a
`RectilinearGrid` in the same projected meters, so one trilinear call covers both the horizontal
and the vertical.

`launch!` takes its work size from `size(output)` and its device from `architecture`; the grid
argument is unused for a size-tuple work spec, so the source grid is passed for consistency.
"""
function interpolate_to_target!(output, source_field, x, y, z, mask, architecture)
    launch!(
        architecture,
        source_field.grid,
        size(output),
        _interpolate_to_target!,
        output,
        x,
        y,
        z,
        mask,
        source_field,
        source_field.grid,
        (Center(), Center(), Center()),
    )

    return output
end

"""
    SourceReader(filepath)

Holds the one source file currently open, reopening on demand. Output steps are written
in date order so consecutive reads almost always hit the same file; only an interpolated day
straddling a month boundary reopens, which happens at most once per gap.
"""
mutable struct SourceReader{D}
    filepath::String
    dataset::D
end

SourceReader(filepath::String) = SourceReader(filepath, NCDataset(filepath))

Base.close(reader::SourceReader) = close(reader.dataset)

"""
    source_slab(reader::SourceReader, record::SourceRecord, source_name)

One `(X, Y, depth)` slab of `source_name`, with masked points as `NaN`. A surface variable with no
depth axis comes back as a single level, per `as_source_slab`.
"""
function source_slab(reader::SourceReader, record::SourceRecord, source_name)
    if reader.filepath != record.filepath
        close(reader.dataset)
        reader.dataset = NCDataset(record.filepath)
        reader.filepath = record.filepath
    end

    plane = read_last_dimension(reader.dataset[source_name], record.index)
    return as_source_slab(Float32.(coalesce.(plane, NaN32)))
end

"""
    blended_slab(reader::SourceReader, step::ForcingTimeStep, source_name)

The source slab for `step`: a single record when the step is present in the download, otherwise
the linear blend of the records bracketing an interpolated day. Land stays `NaN` because it is
`NaN` in both records.
"""
function blended_slab(reader::SourceReader, step::ForcingTimeStep, source_name)
    lower = source_slab(reader, step.lower, source_name)
    iszero(step.weight) && return lower

    upper = source_slab(reader, step.upper, source_name)
    @. lower = (1 - step.weight) * lower + step.weight * upper

    return lower
end

"""
    write_forcing_file(filepath, target_grid, variables, steps, source, architecture)

Write the forcing NetCDF, streaming one time step at a time so peak memory stays at a single
slab rather than a whole variable. Interpolation runs on `architecture`.
"""
function write_forcing_file(filepath, target_grid, variables, steps, source, architecture)
    isdir(dirname(filepath)) || mkpath(dirname(filepath))
    isfile(filepath) && rm(filepath; force = true)

    ds = NCDataset(filepath, "c")
    try
        define_forcing_dimensions!(ds, target_grid, steps)
        for variable in variables
            # One horizontal level of one time step per chunk. Both this writer and
            # `load_from_netcdf` touch exactly one time index at a time, so a chunk spanning
            # several time steps would be decompressed and recompressed once per step; the
            # netCDF default chunking does span them and makes the write an order of magnitude
            # slower.
            chunksizes = [size(variable.mask, 1), size(variable.mask, 2), 1, 1]
            for name in (variable.name, variable.name * "_lambda")
                defVar(ds, name, Float32, variable.dimensions; chunksizes,
                    deflatelevel = FORCING_DEFLATE_LEVEL, attrib = ["_FillValue" => NaN32])
            end
        end

        # All variables of a source live on the same source grid, so one field is enough.
        source_field = Field{Center,Center,Center}(source_field_grid(source, architecture))

        # Target nodes, mask and output buffer move to the interpolation device once; only the
        # finished slab comes back per step, for the NetCDF write.
        devices = [(
            x = on_architecture(architecture, variable.x),
            y = on_architecture(architecture, variable.y),
            z = on_architecture(architecture, variable.z),
            mask = on_architecture(architecture, variable.mask),
            output = on_architecture(architecture, zeros(Float32, size(variable.mask))),
            host = Array{Float32}(undef, size(variable.mask)),
        ) for variable in variables]

        reader = SourceReader(first(steps).lower.filepath)
        announced = ""

        try
            for (index, step) in enumerate(steps)
                if step.lower.filepath != announced
                    announced = step.lower.filepath
                    @info "Regridding $(basename(announced))"
                end

                for (variable, device) in zip(variables, devices)
                    slab = blended_slab(reader, step, variable.source_name)
                    set_source_field!(source_field, slab, variable.mask_fill)
                    interpolate_to_target!(
                        device.output, source_field, device.x, device.y, device.z, device.mask, architecture,
                    )
                    copyto!(device.host, device.output)
                    ds[variable.name][:, :, :, index] = device.host
                    ds[variable.name*"_lambda"][:, :, :, index] = variable.lambda
                end
            end
        finally
            close(reader)
        end
    finally
        close(ds)
    end

    return filepath
end

"""
    define_forcing_dimensions!(ds, target_grid, steps)

Define the dimensions and coordinate variables `forcing_from_file` checks against the grid.
"""
function define_forcing_dimensions!(ds, target_grid, steps)
    coordinates = (
        "Nx" => Array(λnodes(target_grid, Center())),
        "Ny" => Array(φnodes(target_grid, Center())),
        "Nz" => Array(znodes(target_grid, Center())),
        "Nx_faces" => Array(λnodes(target_grid, Face())),
        "Ny_faces" => Array(φnodes(target_grid, Face())),
        "Nz_faces" => Array(znodes(target_grid, Face())),
    )

    for (name, values) in coordinates
        defDim(ds, name, length(values))
    end
    defDim(ds, "time", length(steps))

    for (name, values) in coordinates
        defVar(ds, name, Float64, (name,))[:] = values
    end
    defVar(ds, "time", [step.date for step in steps], ("time",))

    return ds
end

include("rivers.jl")
include("boundaries.jl")
include("norkyst.jl")
include("norkyst_boundaries.jl")
include("of800_rivers.jl")

end # module
