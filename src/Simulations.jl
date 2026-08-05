module Simulations

export SimulationConfig,
    CoupledHydrostaticSimulation,
    SplitExplicitFreeSurfaceConfig,
    free_surface,
    SnapshotWriter,
    CheckpointWriter,
    AdaptiveTimeStep,
    FromForcing,
    FromResults,
    simulation_architecture,
    model_tracers,
    coupled_simulation,
    attach_writer!,
    attach_time_stepping!,
    initial_time_step,
    build_simulation,
    run_simulation

using Oceananigans
using Oceananigans: fields
using Oceananigans.Utils: prettytime
using Oceananigans.TimeSteppers: reset!, update_state!
using NumericalEarth
using Dates: DateTime, Second
using NCDatasets

using ..Configs:
    AbstractSimulationConfig,
    AbstractCoupledSimulationConfig,
    AbstractFreeSurfaceConfig,
    AbstractWriterConfig,
    AbstractTimeSteppingConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    FjordConfig,
    bathymetry_path,
    forcing_path,
    river_forcing_path,
    results_path,
    run_tag
using ..Utils: progress, cell_advection_timescale_coupled_model
using ..Atmospheres: prescribed_atmosphere, prescribed_radiation, atmosphere_date_range
using ..Forcing: simulation_forcing, interpolation_architecture
using ..BoundaryConditions: field_boundary_conditions
using ..Grids: ImmersedBoundaryGrid

"""
    SplitExplicitFreeSurfaceConfig(; cfl)

The one built-in `AbstractFreeSurfaceConfig`: a barotropic split-explicit solver at a given CFL.

Holds only the knob, not the object — `SplitExplicitFreeSurface(grid, cfl = ...)` needs the grid,
which does not exist until `coupled_simulation` builds it, so `free_surface` is what does that
building once the grid is in hand.
"""
struct SplitExplicitFreeSurfaceConfig <: AbstractFreeSurfaceConfig
    cfl::Float64
end

SplitExplicitFreeSurfaceConfig(; cfl) = SplitExplicitFreeSurfaceConfig(Float64(cfl))

"""
    free_surface(config::SplitExplicitFreeSurfaceConfig, grid)

Build the `SplitExplicitFreeSurface` `coupled_simulation` passes to `HydrostaticFreeSurfaceModel`.
"""
free_surface(config::SplitExplicitFreeSurfaceConfig, grid) =
    SplitExplicitFreeSurface(grid, cfl = config.cfl)

"""
    CoupledHydrostaticSimulation(; buoyancy, closure, tracer_advection, momentum_advection,
                                 tracers, coriolis, sea_ice, biogeochemistry, free_surface)

What the model *is*: a `HydrostaticFreeSurfaceModel` inside a NumericalEarth `OceanSeaIceModel`.

Holds the components that depend on neither the grid nor the device, so they are stored as the
objects `coupled_simulation` consumes rather than as knobs it has to interpret. `free_surface` is
the one exception: it is itself an `AbstractFreeSurfaceConfig`, built by its own `free_surface(config,
grid)` hook rather than stored ready-made, for exactly the reason `coupled_simulation` — and not this
struct — is what assembles the model: `SplitExplicitFreeSurface` needs the grid, which only exists
once `coupled_simulation` is called.

**No field has a default**, deliberately. Every one is a scientific choice about a particular fjord,
and a default would let the next setup silently inherit it, so a setup that forgets one gets an
`UndefKeywordError` naming it rather than a plausible-looking run.

Nine type parameters, one per component. A tenth is the signal to collapse them into a single
`NamedTuple` field instead of adding another; the "config" testset guards the count.

Note that `Base.@kwdef` on a *parametric* struct generates a constructor that does not convert, so
`free_surface = SplitExplicitFreeSurfaceConfig(cfl = 1)` is a `MethodError` where `cfl = 1.0` is
fine — though `SplitExplicitFreeSurfaceConfig`'s own keyword constructor converts, the same
exception every nested config gets.
"""
Base.@kwdef struct CoupledHydrostaticSimulation{B,C,TA,MA,TR,CO,SI,BG,FS} <:
                   AbstractCoupledSimulationConfig
    buoyancy::B
    closure::C
    tracer_advection::TA
    momentum_advection::MA
    tracers::TR
    coriolis::CO
    sea_ice::SI
    biogeochemistry::BG
    free_surface::FS
end

"""
    model_tracers(model)

The tracer names a coupled-model config carries.

A hook rather than a field read, because `build_simulation` needs them long before the model exists:
`simulation_forcing` builds one term per tracer, `OpenLateralBoundary` opens one lateral condition per
tracer, and `resolve_initial_conditions` reads one state variable per tracer. All three ask the model
config what it simulates rather than being told a second time.
"""
model_tracers(model::CoupledHydrostaticSimulation) = model.tracers

"""
    SnapshotWriter(; name, output_file, variables, interval, overwrite_existing)

A NetCDF snapshot of named ocean fields, on a `TimeInterval` schedule.

# Fields
- `name`: the key it lives under in the ocean sub-simulation's `output_writers`, and what makes two
  snapshot writers distinguishable. Stated rather than derived, so a test or a post-processing step
  can name the writer it wants.
- `output_file`: resolved by `results_path`, which inserts the run tag and, for a looped run, the
  loop index.
- `variables`: the field names to write, as `Symbol`s, resolved against the model by
  `snapshot_outputs`. Anything `Oceananigans.fields` exposes — velocities, tracers, the free
  surface, auxiliaries — so a setup that adds a biogeochemical tracer writes it by naming it.
- `interval`, `overwrite_existing`: the schedule and the clobber policy.
"""
struct SnapshotWriter{V<:Tuple{Vararg{Symbol}}} <: AbstractWriterConfig
    name::Symbol
    output_file::String
    variables::V
    interval::Float64
    overwrite_existing::Bool
end

function SnapshotWriter(; name, output_file, variables, interval, overwrite_existing)
    interval > 0 || throw(
        ArgumentError("A snapshot writer's `interval` must be positive, got $interval."),
    )
    isempty(variables) &&
        throw(ArgumentError("Snapshot writer :$name names no `variables` to write."))

    return SnapshotWriter(
        Symbol(name),
        String(output_file),
        Tuple(Symbol(variable) for variable in variables),
        Float64(interval),
        Bool(overwrite_existing),
    )
end

"""
    CheckpointWriter(; interval, overwrite_existing, cleanup)

A JLD2 checkpoint of the coupled model's prognostic state, on a `TimeInterval` schedule.

Naming one is what makes a run resumable; a setup that names none writes no checkpoints at all,
which is how "checkpointing off" is now spelled. `pickup` without one is a configuration error
rather than a run that fails on its first `run!`.

Oceananigans 0.110 checkpoints *only* `prognostic_state(simulation)`, so the `Checkpointer`
docstring's warning that objects containing functions cannot be serialized does not apply here:
`ForcingFromFile`, its `FieldTimeSeries` backend and the `FreshwaterExchange` in the tracer top
boundary conditions are all outside that state. CATKE's diffusivities and its `e` tracer are inside
it, which is why a checkpoint is a few hundred MB on a real grid and why `cleanup` exists.
"""
struct CheckpointWriter <: AbstractWriterConfig
    interval::Float64
    overwrite_existing::Bool
    cleanup::Bool
end

function CheckpointWriter(; interval, overwrite_existing, cleanup)
    interval > 0 || throw(
        ArgumentError(
            "A checkpoint writer's `interval` must be positive, got $interval. A setup that wants " *
            "no checkpoints names no `CheckpointWriter` at all.",
        ),
    )

    return CheckpointWriter(Float64(interval), Bool(overwrite_existing), Bool(cleanup))
end

"""
    AdaptiveTimeStep(; initial_time_step, cfl, max_time_step, max_time_step_change)

An Oceananigans time-step wizard: the step is chosen each iteration from the advective CFL, starting
from `initial_time_step` and growing by at most `max_time_step_change` per step.

The wizard uses `cell_advection_timescale_coupled_model`, which reaches through the coupled model to
the ocean — the coupled model has no advective timescale of its own.
"""
struct AdaptiveTimeStep <: AbstractTimeSteppingConfig
    initial_time_step::Float64
    cfl::Float64
    max_time_step::Float64
    max_time_step_change::Float64
end

AdaptiveTimeStep(; initial_time_step, cfl, max_time_step, max_time_step_change) = AdaptiveTimeStep(
    Float64(initial_time_step),
    Float64(cfl),
    Float64(max_time_step),
    Float64(max_time_step_change),
)

"""
    initial_time_step(config)

The `Δt` a simulation starts at, in seconds.

Separate from `attach_time_stepping!` because the two happen at different moments: the step is needed
when the `Simulation` is constructed, the policy only once it exists.
"""
initial_time_step(config::AdaptiveTimeStep) = config.initial_time_step

"""
    SimulationConfig(; results_root, architecture, model, boundary_conditions, writers, ...)

Simulation configuration: how a setup whose data is already prepared is run.

**No field has a default**, here and in each of the four nested configs. Every knob is a scientific
choice about a particular fjord, and a default would let one setup silently inherit another's — so a
setup that forgets one gets an `UndefKeywordError` naming it rather than a plausible-looking run.
The setup file is the complete statement of what the simulation is; nothing about it is hidden here.

The fields split into four nested configs and the run control that ties them together, and the split
is by *what dispatches on what*. `model` is assembled by `coupled_simulation`, `boundary_conditions`
by `field_boundary_conditions`, each writer by `attach_writer!`, and `time_stepping` by
`attach_time_stepping!` — four generic functions, so a different model, a new boundary-condition
piece or another kind of output is a new subtype rather than an edit to `build_simulation`.

# Fields
- `results_root`: the directory every writer, and the run log, resolves against.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `simulation_architecture`. A `Symbol`
  rather than a live `CPU()`/`GPU()` so this config's field types stay concrete and a setup file
  loads on a machine with no GPU.
- `model`: an `AbstractCoupledSimulationConfig`.
- `boundary_conditions`: a tuple of `AbstractBoundaryConditionConfig`s. Tuple order is merge
  precedence; `()` means the model's defaults everywhere.
- `writers`: a tuple of `AbstractWriterConfig`s. `()` writes nothing.
- `time_stepping`: an `AbstractTimeSteppingConfig`.
- `initial_conditions`: where the ocean state starts from — a `NamedTuple` of constants, fields or
  functions, a `FromForcing`, or a `FromResults`. `build_simulation` puts it through
  `resolve_initial_conditions`, so what reaches `coupled_simulation` is always a `NamedTuple` `set!`
  accepts.
- `start_date`: the calendar instant model time zero stands for.
- `stop_time`: seconds of simulated time one pass through the run window lasts.
- `loops`: how many times to run that window, carrying the ocean state across each restart. `1`
  runs it once. See `run_simulation`.
- `progress_interval`: how often the `progress` callback reports.
- `pickup`: resume from the newest checkpoint under `results_root` instead of starting fresh.
  Requires a `CheckpointWriter`.

Nothing that could be derived from another config is a field at all: which forcing file, which open
boundary and which atmosphere reader all come from `forcing_config` and `atmosphere_config`.

`start_date` is a field rather than something derived from an input file because every prepared
file has its *own* first record, and each reader used to zero its own axis there: the Oslofjord
forcing starts at 12:00 and its atmosphere at 00:00, which silently ran the two twelve hours out
of phase. One stated instant is the only thing they can all agree on, and
`validate_time_coverage` checks each file actually spans `[start_date, start_date + stop_time]`
rather than letting `Cyclical()` wrap the shortfall.

Every duration is in seconds. `Base.@kwdef` on a *parametric* struct generates a constructor that
does not convert, so `stop_time = 3600` is a `MethodError` where `stop_time = 1hour` is fine — write
durations with `Oceananigans.Units`, whose constants are already `Float64`. `loops` is the exception,
being an `Int` already, and the nested configs are the other one: each has a hand-written keyword
constructor that converts, so an integer duration is fine in a writer or a time-stepping config.
"""
Base.@kwdef mutable struct SimulationConfig{M,BC,W,TS,I} <: AbstractSimulationConfig
    results_root::String
    architecture::Symbol
    model::M
    boundary_conditions::BC
    writers::W
    time_stepping::TS
    initial_conditions::I
    start_date::DateTime
    stop_time::Float64
    loops::Int
    progress_interval::Float64
    pickup::Bool
end

"""
    simulation_architecture(config)

Resolve `config.architecture` to a live `CPU()` or `GPU()`.

Shares the `Val`-dispatched methods of `interpolation_architecture`, so `:auto` picks the GPU when
`CUDA.functional()` and `:gpu` errors rather than silently falling back — the same behavior, and
the same message, as the forcing interpolation.
"""
simulation_architecture(config::AbstractSimulationConfig) =
    interpolation_architecture(Val(config.architecture))

"""
The global attribute a snapshot file records its `start_date` in.

The snapshot writer's own time axis is seconds from model zero, which says nothing about the
calendar instant zero stood for — so a later `FromResults(path, date)` would have no way to turn a
date into a record. Writing the instant into the file it describes keeps that knowledge with the
data instead of making it a second config field that could disagree with the first.
"""
const RESULTS_START_DATE_ATTRIBUTE = "start_date"

"""
    FromForcing(date = nothing)

Initial conditions read from the prepared forcing file the simulation is already reading, at
`date` — or at the simulation's `start_date` when `date` is `nothing`.

The file is on the model grid by construction (`prepare_forcing` regrids onto it), so this is a
plain read: no regridding, no inpainting, no `NumericalEarth` dataset wrapper.
"""
struct FromForcing{D}
    date::D
end

FromForcing() = FromForcing(nothing)

"""
    FromResults(path, date = nothing)

Initial conditions read from a previous run's snapshot file, at `date` — or from its last record
when `date` is `nothing`.

A relative `path` resolves against `results_root`, like `output_file` does. Naming a `date` requires
the file to carry the `$RESULTS_START_DATE_ATTRIBUTE` attribute `build_simulation` writes, since a
snapshot's time axis is seconds from its own model zero; files written before that attribute existed
can only be read by their last record.
"""
struct FromResults{D}
    path::String
    date::D
end

FromResults(path::String) = FromResults(path, nothing)

"""
    initial_conditions_date(initial_conditions, start_date)

Which instant to read, defaulting an unnamed date to the run's own `start_date`.
"""
initial_conditions_date(::FromForcing{Nothing}, start_date) = start_date
initial_conditions_date(initial_conditions::FromForcing, start_date) = initial_conditions.date

"""
    resolve_initial_conditions(initial_conditions, grid, forcing_file, config)

Turn a setup's `initial_conditions` into the `NamedTuple` `set!(model; ...)` consumes.

Dispatched on the kind of source, so `coupled_simulation` never learns that there is more than one:
it still receives something splattable and still applies it with a single `set!`.
A `NamedTuple` — constants, functions or fields — passes straight through, which is what every
setup did before the other two existed.

What gets set is every tracer `model_tracers(config.model)` names plus `u` and `v`, intersected with
what the source file actually carries — see `state_variables`. Nothing is enumerated here, so adding
a biogeochemical tracer to a setup is enough to have it read back.

The free surface `η` and any closure-owned tracer the config does not name (CATKE's `e`) keep their
defaults, so reading a state in is a *warm start*, not a restart: the barotropic mode and the
turbulence field re-adjust over the first hours. Use `pickup` for an exact continuation.
"""
resolve_initial_conditions(initial_conditions::NamedTuple, grid, forcing_file, config) =
    initial_conditions

function resolve_initial_conditions(initial_conditions::FromForcing, grid, forcing_file, config)
    date = initial_conditions_date(initial_conditions, config.start_date)
    @info "Initial conditions from forcing $forcing_file at $date"
    return forcing_state(forcing_file, grid, date, model_tracers(config.model))
end

function resolve_initial_conditions(initial_conditions::FromResults, grid, forcing_file, config)
    filepath = results_state_path(initial_conditions, config)
    isfile(filepath) || error("Results file $filepath does not exist.")
    @info "Initial conditions from results $filepath"
    return results_state(filepath, grid, initial_conditions.date, model_tracers(config.model))
end

"""
    results_state_path(initial_conditions, config)

Where `FromResults` reads from: its `path` as given when absolute, resolved against `results_root`
when relative — the same rule `output_file` follows, so a previous run's output can be named by its
filename alone.
"""
results_state_path(initial_conditions::FromResults, config) =
    isabspath(initial_conditions.path) ? initial_conditions.path :
    joinpath(config.results_root, initial_conditions.path)

"""
    forcing_state(filepath, grid, date, tracers)

The state variables a prepared forcing file holds at `date`.

The `_lambda` twins are relaxation rates rather than state, so only the bare names are read.
"""
function forcing_state(filepath, grid, date, tracers)
    return NCDataset(filepath) do ds
        validate_state_dimensions(filepath, ds, grid, ("Nx", "Ny", "Nz"))
        dates = ds["time"][:]
        index = findfirst(==(date), dates)
        isnothing(index) && error(
            "No forcing record at $date in $filepath, whose axis runs $(first(dates)) to " *
            "$(last(dates)). Name a date on that axis, or prepare forcing covering $date.",
        )
        return state_variables(ds, index, eltype(grid), tracers)
    end
end

"""
    results_state(filepath, grid, date, tracers)

The state variables a previous run's snapshot file holds at `date`, or at its last record when
`date` is `nothing`.
"""
function results_state(filepath, grid, date, tracers)
    return NCDataset(filepath) do ds
        validate_state_dimensions(filepath, ds, grid, ("λ_caa", "φ_aca", "z_aac"))
        return state_variables(ds, results_record_index(filepath, ds, date), eltype(grid), tracers)
    end
end

"""
    results_record_index(filepath, ds, date)

Which record of a snapshot file `date` names. Its time axis is seconds from the run's own model
zero, so the instant that zero stood for has to come from the file's own
`$RESULTS_START_DATE_ATTRIBUTE` attribute; `nothing` takes the last record and needs no attribute.
"""
results_record_index(filepath, ds, ::Nothing) = ds.dim["time"]

function results_record_index(filepath, ds, date::DateTime)
    haskey(ds.attrib, RESULTS_START_DATE_ATTRIBUTE) || error(
        "$filepath carries no `$RESULTS_START_DATE_ATTRIBUTE` attribute, so its time axis " *
        "(seconds from that run's model zero) cannot be turned into a date. Read its last record " *
        "with `FromResults(path)` instead, or re-run the simulation that wrote it.",
    )

    source_start = DateTime(ds.attrib[RESULTS_START_DATE_ATTRIBUTE])
    seconds = Second(date - source_start).value
    times = ds["time"][:]
    index = argmin(abs.(times .- seconds))

    abs(times[index] - seconds) <= 1 || error(
        "No record within 1 s of $date in $filepath, whose axis runs $source_start to " *
        "$(source_start + Second(round(Int, last(times)))).",
    )

    return index
end

"""
    state_variables(ds, index, FT, tracers)

The state variables `ds` carries at time `index`, as a `NamedTuple` of `FT` arrays.

Which variables those are comes from the simulation config's `tracers` plus the two horizontal
velocities — never a list written out here, so a setup that adds a biogeochemical tracer gets it read
back without this module being touched. The rule is
`(map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)`, the same one `forcing_from_file` uses to decide
which forcing terms to build, so the two cannot disagree about what the state is.

Intersecting with the file's own variables is what lets one reader serve a forcing file and a
snapshot file: they need not carry the same set, and a tracer the source lacks is simply left at its
default rather than being an error.
"""
function state_variables(ds, index, FT, tracers)
    names = (map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)
    return (; (Symbol(name) => finite_slab(ds[name][:, :, :, index], FT) for name in names)...)
end

"""
    finite_slab(data, FT)

`data` as `FT` with every missing or non-finite cell replaced by zero.

Both sources mark land that way, and in both it is exactly the cells this grid immerses — the
forcing file's mask comes from `peripheral_node` on this very grid, and a snapshot was written from
a field on it — so zeroing them fills only cells the model never reads.

`FT` is the grid's element type, not the file's: `set!` moves the array to the model's architecture
and `copyto!`s it into the field, and a `Union{Missing,Float32}` array cannot become a `CuArray` at
all while a mismatched element type need not convert on the device.
"""
finite_slab(data, FT) =
    map(value -> ismissing(value) || !isfinite(value) ? zero(FT) : convert(FT, value), data)

"""
    validate_state_dimensions(filepath, ds, grid, names)

Check the file's horizontal and vertical extents are the grid's, so a state file from another setup
is a clear error rather than a silent misread.
"""
function validate_state_dimensions(filepath, ds, grid, names)
    expected = size(grid)
    found = ntuple(index -> ds.dim[names[index]], 3)
    found == expected || throw(
        DimensionMismatch("$filepath is $found but the simulation grid is $expected"),
    )
    return nothing
end

"""
    simulation_forcing_path(config)

Which prepared forcing file the simulation reads: the rivers-augmented copy `add_rivers` wrote
when the setup names rivers, so they are never silently dropped, and the file `prepare_forcing`
wrote when it does not. Dispatched on `config.rivers` exactly like `add_rivers` — or `nothing`,
for a setup naming no forcing at all.
"""
simulation_forcing_path(config::FjordConfig) = simulation_forcing_path(config.forcing_config)

simulation_forcing_path(::Nothing) = nothing
simulation_forcing_path(config::AbstractForcingConfig) =
    simulation_forcing_path(config, config.rivers)

simulation_forcing_path(config::AbstractForcingConfig, ::Nothing) = forcing_path(config)
simulation_forcing_path(::AbstractForcingConfig, rivers::AbstractRiverConfig) =
    river_forcing_path(rivers)

"""
    forcing_prerequisite(rivers)

The subcommand that writes the file `simulation_forcing_path` picked, for its error message.
"""
forcing_prerequisite(::Nothing) = "prepare_forcing"
forcing_prerequisite(::AbstractRiverConfig) = "add_rivers"

"""
    resolve_forcing_file(config::FjordConfig)

The prepared forcing file `build_simulation` reads: `simulation_forcing_path(config)`, checked to
exist and reported by the step that writes it if it does not — or `nothing`, for a setup naming
no forcing at all, which is not a prerequisite to check.
"""
resolve_forcing_file(config::FjordConfig) = resolve_forcing_file(config.forcing_config)

resolve_forcing_file(::Nothing) = nothing

function resolve_forcing_file(forcing_config::AbstractForcingConfig)
    forcing_file = simulation_forcing_path(forcing_config)
    isfile(forcing_file) || error(
        "Prepared forcing $forcing_file does not exist. Run `julia --project -m FjordSim " *
        "$(forcing_prerequisite(forcing_config.rivers))` for this setup first.",
    )
    return forcing_file
end

"""
    validate_time_coverage(label, range, start_date, stop_time, remedy)

Check a prepared file's date `range` contains the run's whole interval.

Both readers use `Cyclical()` time indexing, which does not fail outside the data it was given —
it wraps, so a run that outlasts its forcing quietly replays the beginning and a run that starts
before it quietly reads the end. This is what makes that unreachable instead of merely unlikely.
A `nothing` range is a source that cannot report its dates, and is skipped.
"""
validate_time_coverage(label, ::Nothing, start_date, stop_time, remedy) = nothing

function validate_time_coverage(label, range, start_date, stop_time, remedy)
    first_available, last_available = range
    end_date = start_date + Second(round(Int, stop_time))

    start_date >= first_available || error(
        "The run starts at $start_date but the $label only begins at $first_available. Move the " *
        "simulation config's `start_date` to $first_available or later, or $remedy.",
    )
    end_date <= last_available || error(
        "The run ends at $end_date but the $label stops at $last_available. Shorten `stop_time` " *
        "to $(Second(last_available - start_date).value) seconds or less, or $remedy.",
    )

    return nothing
end

"""
    forcing_date_range(filepath)

First and last date on a prepared forcing file's time axis.
"""
forcing_date_range(filepath) = NCDataset(filepath) do ds
    dates = ds["time"][:]
    (first(dates), last(dates))
end

forcing_date_range(::Nothing) = nothing

"""
    loop_output_path(writer, config, loop)

One writer's file for one repetition: the plain run-tagged name for a setup that runs its window
once, and the loop-indexed one when it repeats. A single run therefore keeps the shorter name
rather than gaining a `_loop01` nobody asked for.
"""
loop_output_path(writer::AbstractWriterConfig, config::AbstractSimulationConfig, loop) =
    config.loops == 1 ? results_path(writer, config) : results_path(writer, config, loop)

"""
    CheckpointTrait

Whether a writer config contributes a `Checkpointer` — `Checkpointing()` or `NotCheckpointing()`,
defaulting to the latter.

A trait rather than an `isa` test at each of the three sites that need the answer, because the
answer has to be exact. `run!(…; checkpoint_at_end)` with no checkpointer to find does not fail: it
writes `checkpoint_iteration<N>.jld2` into the working directory behind a `@warn`, and
`run!(…; pickup)` with two of them cannot tell which to resume. It is also what lets
`build_simulation` reject `pickup` without a checkpointing writer as a configuration error, which
was not expressible while checkpointing was a threshold on a float.
"""
abstract type CheckpointTrait end
struct Checkpointing <: CheckpointTrait end
struct NotCheckpointing <: CheckpointTrait end

checkpoint_trait(::AbstractWriterConfig) = NotCheckpointing()
checkpoint_trait(::CheckpointWriter) = Checkpointing()

"""
    checkpoints(writer)
    checkpoints(config)

Whether a writer contributes a `Checkpointer`, or whether any of a simulation config's does.
"""
checkpoints(writer::AbstractWriterConfig) = checkpoints(checkpoint_trait(writer))
checkpoints(::Checkpointing) = true
checkpoints(::NotCheckpointing) = false
checkpoints(config::AbstractSimulationConfig) = any(checkpoints, config.writers)

"""
    OutputPathTrait

Whether a writer names a file the run should report — `NamesOutputFile()` or
`NamesNoOutputFile()`, defaulting to the latter.

Not simply the negation of `CheckpointTrait`: a checkpointer does write files, but they are
scratch rather than product — they carry no run tag, `cleanup` prunes them, and a `pickup` finds
them by scanning rather than by being told a path. A writer could reasonably be both or neither.
"""
abstract type OutputPathTrait end
struct NamesOutputFile <: OutputPathTrait end
struct NamesNoOutputFile <: OutputPathTrait end

output_path_trait(::AbstractWriterConfig) = NamesNoOutputFile()
output_path_trait(::SnapshotWriter) = NamesOutputFile()

"""
    reported_paths(writer, config, loop)

The files `writer` writes for repetition `loop` that `run_simulation` should report, as a tuple —
empty for a writer that names none, so flattening over the whole tuple is total.
"""
reported_paths(writer::AbstractWriterConfig, config::AbstractSimulationConfig, loop) =
    reported_paths(output_path_trait(writer), writer, config, loop)

reported_paths(::NamesOutputFile, writer, config, loop) =
    (loop_output_path(writer, config, loop),)

reported_paths(::NamesNoOutputFile, writer, config, loop) = ()

"""
    loop_output_paths(config, loop)

Every product file one repetition writes.
"""
loop_output_paths(config::AbstractSimulationConfig, loop) = String[
    path for writer in config.writers for path in reported_paths(writer, config, loop)
]

"""
    writer_keys(writer)

The `output_writers` keys a writer occupies, as a tuple. Used to reject a setup naming two writers
that would silently replace one another in the same dictionary.
"""
writer_keys(writer::AbstractWriterConfig) = writer_keys(output_path_trait(writer), writer)
writer_keys(::NamesOutputFile, writer) = (writer.name,)
writer_keys(::NamesNoOutputFile, ::AbstractWriterConfig) = ()

"""
    checkpoint_prefix(loop)

Prefix for one repetition's checkpoints, under which `Checkpointer` writes
`<prefix>_iteration<N>.jld2`.

The loop index is in the *name* because it is not in the checkpoint: the state records the clock but
not which repetition produced it, so without this `pickup` could not tell a loop-3 checkpoint from a
loop-1 one and would replay the whole spin-up. `resume_loop` reads the index back out.

The run tag is deliberately *not* in the name, even though `results_path` carries it: the tag is the
launch instant, so a later launch could not name — and so could not resume — the checkpoints of the
one before it. The cost is that checkpoints are shared per `results_root`, which is spelled out in
`resume_loop`.
"""
checkpoint_prefix(loop) = string("checkpoint_loop", lpad(loop, 2, '0'))

"""
    checkpointed_loops(config)

Every loop index that has a checkpoint under `results_root`.

The pattern is the exact inverse of `checkpoint_prefix`, and the coupling is silent: a
`CheckpointWriter` that gained a `prefix` field would leave this finding nothing, so `resume_loop`
would warn and replay the whole spin-up. Change both or neither.
"""
function checkpointed_loops(config::AbstractSimulationConfig)
    isdir(config.results_root) || return Int[]
    pattern = r"^checkpoint_loop(\d+)_iteration\d+\.jld2$"
    found = filter(!isnothing, match.(pattern, readdir(config.results_root)))
    return [parse(Int, captured[1]) for captured in found]
end

"""
    resume_loop(config)

Which repetition to start at: the highest one that has a checkpoint when `pickup` is set, and `1`
otherwise.

`build_simulation` attaches its writers for this loop rather than unconditionally for the first, so a
resumed run writes into the loop it left off in.

`pickup` with no checkpointing writer is rejected here rather than left to fail inside `run!`,
because it is a statement the config contradicts: nothing this run writes could ever be resumed
from. That check only became expressible once checkpointing was a writer rather than a threshold on
a float.

Since checkpoints carry no run tag, this is whatever state `results_root` holds, whichever launch
wrote it and whatever `start_date` it was written under — there is one resumable run per results
directory, and a second launch overwrites the first's checkpoints. The snapshots are per launch, so a
resumed run's file starts at the resume point and the records before it stay in the previous launch's
file.
"""
function resume_loop(config::AbstractSimulationConfig)
    config.pickup || return 1

    checkpoints(config) || throw(
        ArgumentError(
            "`pickup` is set but no writer checkpoints, so there is nothing to resume from. Add a " *
            "`CheckpointWriter` to `writers`, or set `pickup = false`.",
        ),
    )

    loops = checkpointed_loops(config)
    if isempty(loops)
        @warn "`pickup` is set but no checkpoint was found in $(config.results_root); " *
              "starting from the beginning"
        return 1
    end

    return maximum(loops)
end

"""
    attach_writers!(simulation, config, loop)

Point every writer the config names at `loop`'s own files.

Called once by `build_simulation` and again by `restart_loop!` for each subsequent repetition. What
each writer does with that is `attach_writer!`'s business — including which simulation it attaches
to, which is not the same for all of them.
"""
function attach_writers!(simulation, config::AbstractSimulationConfig, loop)
    for writer in config.writers
        attach_writer!(simulation, writer, config, loop)
    end

    return simulation
end

"""
    snapshot_outputs(writer, ocean_model)

The fields `writer.variables` names, as the `NamedTuple` `NetCDFWriter` consumes.

Resolved through `Oceananigans.fields`, the model's own flattened view of its velocities, free
surface, tracers and auxiliary fields — so a setup that adds a biogeochemical tracer, or names `w`,
`η` or CATKE's `e`, gets it written without this module learning what those are. Nothing is
enumerated here, the same rule `state_variables` follows on the read side.

A name the model does not have is an **error**, unlike `state_variables`, which intersects and moves
on. The asymmetry is deliberate and worth keeping: that function serves two kinds of file whose
variable sets legitimately differ and which no config named, while here the setup wrote the name
down. An over-eager error costs a typo fix; a silently dropped variable costs a whole run.
"""
function snapshot_outputs(writer::SnapshotWriter, ocean_model)
    available = fields(ocean_model)
    unknown = filter(name -> !haskey(available, name), writer.variables)

    isempty(unknown) || throw(
        ArgumentError(
            "Snapshot writer :$(writer.name) names $(join(unknown, ", ")), which the ocean model " *
            "does not have. Available: $(join(keys(available), ", ")).",
        ),
    )

    return NamedTuple(name => available[name] for name in writer.variables)
end

"""
    attach_writer!(simulation, writer::SnapshotWriter, config, loop)

Attach a NetCDF snapshot writer to the *ocean* sub-simulation, under `writer.name`.

A *fresh* writer each loop rather than a renamed one, because a `TimeInterval` accumulates its
actuation count: reused after a clock reset it would believe it was thousands of records ahead and
never fire again. The old one is closed first — it holds an open `NCDataset`, so replacing it
silently would leak a handle and an unflushed tail per loop.

The file records `start_date` as a global attribute so `FromResults(path, date)` can later turn a
date into a record — and, since the name carries the launch instant rather than the simulated one,
it is also the only place the window a file covers is written down. See
`RESULTS_START_DATE_ATTRIBUTE`.

The `mkpath` is load-bearing: `NetCDFWriter` creates only its `dir` keyword, which is left at the
default while the whole path goes in as `filename`, so nothing else would create the directory an
absolute `output_file` points into.
"""
function attach_writer!(simulation, writer::SnapshotWriter, config::AbstractSimulationConfig, loop)
    ocean_sim = simulation.model.ocean

    haskey(ocean_sim.output_writers, writer.name) &&
        close(pop!(ocean_sim.output_writers, writer.name))

    filepath = loop_output_path(writer, config, loop)
    mkpath(dirname(filepath))

    ocean_sim.output_writers[writer.name] = NetCDFWriter(
        ocean_sim.model,
        snapshot_outputs(writer, ocean_sim.model);
        filename = filepath,
        schedule = TimeInterval(writer.interval),
        overwrite_existing = writer.overwrite_existing,
        global_attributes = Dict(RESULTS_START_DATE_ATTRIBUTE => string(config.start_date)),
    )

    return simulation
end

"""
    attach_writer!(simulation, writer::CheckpointWriter, config, loop)

Attach a `Checkpointer` to the *coupled* simulation, under `loop`'s own prefix.

The coupled one, not the ocean one, for two reasons that both fail otherwise: `prognostic_state` of
the `OceanSeaIceModel` is what a resumable state actually is, and `run!(…; pickup)` looks for its
checkpointer in `simulation.output_writers`.

Replaced rather than closed and popped, unlike the snapshot writer: a `Checkpointer` holds no open
file handle, and its schedule is rebuilt with it.
"""
function attach_writer!(simulation, writer::CheckpointWriter, config::AbstractSimulationConfig, loop)
    simulation.output_writers[:checkpointer] = Checkpointer(
        simulation.model;
        schedule = TimeInterval(writer.interval),
        dir = config.results_root,
        prefix = checkpoint_prefix(loop),
        overwrite_existing = writer.overwrite_existing,
        cleanup = writer.cleanup,
    )

    return simulation
end

"""
    attach_time_stepping!(simulation, config::AdaptiveTimeStep)

Install an Oceananigans time-step wizard on `simulation`.

Deliberately does not touch `simulation.Δt`: the starting step is `initial_time_step`, applied once
when the `Simulation` is constructed, and `restart_loop!` leaves `Δt` alone so a repetition keeps the
step the previous one converged on rather than re-ramping from one second.
"""
function attach_time_stepping!(simulation, config::AdaptiveTimeStep)
    conjure_time_step_wizard!(
        simulation;
        cfl = config.cfl,
        max_Δt = config.max_time_step,
        max_change = config.max_time_step_change,
        cell_advection_timescale = cell_advection_timescale_coupled_model,
    )

    return simulation
end

"""
    rewind_clock!(component)

Send one coupled-model component's clock back to zero, or do nothing for a component that has no
clock.

`NumericalEarth`'s own `reset_clock!(::EarthSystemModel)` cannot be used for this. Its per-component
fallback is `reset!(getproperty(component, :clock))` and `components` includes `sea_ice`, so a
`FreezingLimitedOceanTemperature` — which is a liquidus and nothing else — makes it throw
`has no field clock`. Dispatching on whether the component *has* a clock keeps this module free of
any named component, so a setup that later names a real sea-ice or land model works unchanged, and
`Val` resolves the question at compile time.
"""
rewind_clock!(::Nothing) = nothing
rewind_clock!(component::Simulation) = rewind_clock!(component.model)
rewind_clock!(component) = rewind_clock!(component, Val(hasfield(typeof(component), :clock)))
rewind_clock!(component, ::Val{true}) = reset!(component.clock)
rewind_clock!(component, ::Val{false}) = nothing

"""
    restart_loop!(simulation, config, loop)

Send `simulation` back to the start of its window for repetition `loop`, keeping the ocean state.

Every clock in the coupled model goes back to zero — the coupled model's own, and each component's,
found by sweeping its properties rather than naming them. Nothing else is reset: the ocean's
velocities, tracers, CATKE diffusivities and free surface all carry over, which is the whole point.
`stop_time` is untouched, and so is `simulation.Δt`, so the time-step wizard keeps the step it
converged on instead of re-ramping from one second.

Two things then have to be prodded by hand, and both are silent if they are not.

The prescribed atmosphere and radiation hold their `FieldTimeSeries` window, and rewinding a clock
does not refill it. `time_step!(::EarthSystemModel)` assembles surface fluxes in
`maybe_prepare_first_time_step!` *before* it steps the atmosphere, so without this the first step of
each loop would be forced by last December — which is exactly why `NumericalEarth`'s own
`reset_clock!` ends with an `update_state!` on the atmosphere.

And the ocean is a `Simulation` of its own, which `run!` never touches: `run!` clears `initialized`
on the coupled simulation only, so `time_step!(ocean_sim)` would skip `initialize!` and the fresh
snapshot writer would never have its schedule initialized or its `t = 0` record written.

`clock.last_Δt` comes back as `Inf`, which makes the first step of each loop a forward Euler step.
That is both harmless and right: `G⁻` refers to a step taken under forcing a year away.
"""
function restart_loop!(simulation, config::AbstractSimulationConfig, loop)
    model = simulation.model

    rewind_clock!(model)
    for name in propertynames(model)
        rewind_clock!(getproperty(model, name))
    end

    update_state!(model.atmosphere)
    update_state!(model.radiation)

    model.ocean.initialized = false

    return attach_writers!(simulation, config, loop)
end

"""
    build_simulation(config::FjordConfig)

Assemble the coupled simulation a setup describes, with its output writer, progress callback and
time-step wizard attached, and return it without running it.

Split from `run_simulation` so a run can be inspected or stepped by hand:

```julia
simulation = build_simulation(oslofjorden())
run!(simulation)
```

Every input comes from the setup's own prepared files — `bathymetry_path`, the forcing file
`simulation_forcing_path` picks, and the atmosphere the `prescribed_atmosphere` and
`prescribed_radiation` hooks read. A missing prerequisite is reported naming the step that writes
it, rather than as a read error from deep inside NetCDF.

The initial conditions go through `resolve_initial_conditions`, so a `FromForcing` or `FromResults`
is read here — where the grid and the forcing path are known — and what `coupled_simulation`
receives is always a plain `NamedTuple`. Note that `pickup` supersedes them: the checkpoint restores
the state that `set!` had just written.

Returns `nothing` for a setup naming no simulation config.
"""
function build_simulation(config::FjordConfig)
    simulation_config = config.simulation_config
    isnothing(simulation_config) && return nothing

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    forcing_file = resolve_forcing_file(config)

    simulation_config.loops >= 1 ||
        throw(ArgumentError("loops must be at least 1, got $(simulation_config.loops)"))
    validate_writers(simulation_config)

    start_date = simulation_config.start_date
    stop_time = simulation_config.stop_time
    # One window, not `loops` of them: every repetition replays the same interval, so what the
    # prepared files have to span does not grow with the loop count.
    validate_time_coverage(
        "forcing $forcing_file",
        forcing_date_range(forcing_file),
        start_date,
        stop_time,
        "prepare more years of forcing",
    )
    validate_time_coverage(
        "atmosphere",
        atmosphere_date_range(config.atmosphere_config),
        start_date,
        stop_time,
        "prepare more years of atmosphere",
    )

    architecture = simulation_architecture(simulation_config)
    grid = ImmersedBoundaryGrid(bathymetry_file, architecture, config.grid_config.halo)
    tracers = model_tracers(simulation_config.model)

    # The one place a results directory is created. `coupled_simulation` is about assembling a
    # model, and `NetCDFWriter` only creates the `dir` keyword it is not given.
    mkpath(simulation_config.results_root)

    # Every time axis is zeroed at the same instant, so the components stay in phase: each
    # prepared file has its own first record, and left to itself each reader would zero there.
    # Dispatched on the forcing config rather than a hardcoded `forcing_from_file` call, so a
    # source whose prepared files are not that NetCDF contract could read a different way.
    forcing = simulation_forcing(config.forcing_config, grid, forcing_file, tracers, start_date)
    # Named `model_boundary_conditions` rather than `boundary_conditions`, which is the hook each
    # piece of it was built by: assigning to that name here would make it a local and shadow the
    # function for the rest of this body.
    model_boundary_conditions = field_boundary_conditions(
        simulation_config.boundary_conditions,
        grid,
        forcing,
        config.forcing_config,
        tracers,
    )

    initial_conditions = resolve_initial_conditions(
        simulation_config.initial_conditions,
        grid,
        forcing_file,
        simulation_config,
    )

    simulation = coupled_simulation(
        simulation_config.model,
        grid;
        forcing = forcing,
        boundary_conditions = model_boundary_conditions,
        initial_conditions = initial_conditions,
        atmosphere = prescribed_atmosphere(
            config.atmosphere_config,
            architecture;
            reference_date = start_date,
        ),
        radiation = prescribed_radiation(
            config.atmosphere_config,
            architecture;
            reference_date = start_date,
        ),
        stop_time = stop_time,
        initial_time_step = initial_time_step(simulation_config.time_stepping),
    )

    simulation.callbacks[:progress] =
        Callback(progress, TimeInterval(simulation_config.progress_interval))

    # For the loop `run_simulation` will start at, not unconditionally the first: a resumed run has
    # to write into the loop it left off in, and to checkpoint under that loop's prefix.
    start_loop = resume_loop(simulation_config)
    attach_writers!(simulation, simulation_config, start_loop)

    @info "Run $(run_tag(simulation_config)): output to " *
          "$(join(loop_output_paths(simulation_config, start_loop), ", "))"

    attach_time_stepping!(simulation, simulation_config.time_stepping)

    return simulation
end

"""
    validate_writers(config)

Reject a writers tuple the simulation cannot honour, before anything is read or allocated.

Two failure modes, both silent otherwise. More than one checkpointing writer: `run!(…; pickup)`
requires exactly one to resume from, and `checkpoint_at_end` with several writes its file behind a
`@warn`. And two writers under the same key: the second simply replaces the first in the
`output_writers` dictionary, so one of the files the setup asked for is never written.
"""
function validate_writers(config::AbstractSimulationConfig)
    checkpointing = count(checkpoints, config.writers)
    checkpointing <= 1 || throw(
        ArgumentError(
            "A simulation may name at most one checkpointing writer, got $checkpointing. " *
            "`run!` cannot tell which of several to resume from.",
        ),
    )

    # Not named `keys`: a local of that name shadows `Base.keys` for the rest of the body.
    names = [name for writer in config.writers for name in writer_keys(writer)]
    allunique(names) || throw(
        ArgumentError(
            "Writer names must be unique, got $(join(names, ", ")). Two writers under one name " *
            "replace each other and only the last one writes.",
        ),
    )

    return nothing
end

"""
    run_simulation(config::FjordConfig)

Build the simulation a setup describes and run its window `loops` times.

This is the setup-level driver, the same shape as `prepare_forcing(config::FjordConfig)`, and the
last step of the pipeline: it needs `prepare_bathymetry`, `prepare_forcing`, `add_rivers` if the
setup names rivers, and `prepare_atmosphere` if it names an atmosphere.

# Looping

Each repetition replays `[start_date, start_date + stop_time]` with the ocean state carried over, so
`loops > 1` is a spin-up: one forcing year run again and again until the deep basins stop drifting.
`restart_loop!` does the carrying over — the clock goes back to zero, the state does not — and each
repetition writes its own file, so the loops can be compared rather than overwriting one another.

The clock is reset rather than left to run on for `loops * stop_time` because that keeps the run
inside `[t¹, tᴺ]` of every prepared file. Both readers use `Cyclical()` time indexing, and a
monotonic clock would lean on its wrap-around, whose period is inferred from the file's own axis and
so stops matching the loop the moment a file is padded to a different window.

Returns `nothing` for a setup naming no simulation config, which is how `FjordSim.CLI.main`
reports a step a setup opts out of; otherwise `(; simulation, output_files)`.
"""
function run_simulation(config::FjordConfig)
    simulation = build_simulation(config)
    isnothing(simulation) && return nothing

    simulation_config = config.simulation_config
    # Asked again rather than threaded out of `build_simulation`, whose contract is the bare
    # `Simulation`: it is a directory scan, and both entry points have to agree on the answer
    # anyway. The one visible cost is that a `pickup` with no checkpoints warns twice.
    start_loop = resume_loop(simulation_config)
    output_files = String[]

    for loop = start_loop:simulation_config.loops
        loop > start_loop && restart_loop!(simulation, simulation_config, loop)

        loop_files = loop_output_paths(simulation_config, loop)
        append!(output_files, loop_files)
        @info "Loop $loop of $(simulation_config.loops): output to $(join(loop_files, ", "))"

        # `pickup` applies to the loop we resume into and to no other: the loops after it start from
        # a reset clock, and picking up there would send them back to the checkpoint every time.
        run!(
            simulation;
            pickup = simulation_config.pickup && loop == start_loop,
            checkpoint_at_end = checkpoints(simulation_config),
        )

        @info "Loop $loop reached $(prettytime(simulation.model.clock.time)) of model time"
    end

    # A setup may legitimately name no writer at all, in which case there is nothing to report.
    isempty(output_files) || @info "Output saved to $(join(output_files, ", "))"

    return (; simulation, output_files)
end

"""
    coupled_simulation(model, grid; forcing, boundary_conditions, initial_conditions,
                       atmosphere, radiation, stop_time, initial_time_step)

Assemble the coupled `Simulation` a model config describes, on `grid`, and return it without
running it.

The low-level entry point, dispatched on the model config and taking every grid-dependent component
already built — `build_simulation` is what takes a `FjordConfig` and derives them. Adding a
different kind of model is a new `AbstractCoupledSimulationConfig` subtype and a new method here;
neither `build_simulation` nor anything above it changes.

The method below builds a `HydrostaticFreeSurfaceModel` inside a NumericalEarth `OceanSeaIceModel`.
`model.free_surface`'s own `free_surface(config, grid)` hook builds the free-surface object here
rather than it arriving prebuilt, because `SplitExplicitFreeSurface` needs the grid, which is
exactly what this function is given and the config is not.
"""
function coupled_simulation(
    model::CoupledHydrostaticSimulation,
    grid;
    forcing,
    boundary_conditions,
    initial_conditions,
    atmosphere,
    radiation,
    stop_time,
    initial_time_step,
)
    @info "Compiling HydrostaticFreeSurfaceModel"
    ocean_model = HydrostaticFreeSurfaceModel(
        grid;
        buoyancy = model.buoyancy,
        closure = model.closure,
        tracer_advection = model.tracer_advection,
        momentum_advection = model.momentum_advection,
        tracers = model.tracers,
        free_surface = free_surface(model.free_surface, grid),
        coriolis = model.coriolis,
        forcing = forcing,
        boundary_conditions = boundary_conditions,
        biogeochemistry = model.biogeochemistry,
    )
    @info "Compiled HydrostaticFreeSurfaceModel"

    set!(ocean_model; initial_conditions...)

    Δt = initial_time_step
    ocean_sim = Simulation(ocean_model; Δt, stop_time)
    coupled_model = OceanSeaIceModel(ocean_sim, model.sea_ice; atmosphere, radiation)
    @info "Initialized coupled model"

    return Simulation(coupled_model; Δt, stop_time)
end

end  # module Simulations
