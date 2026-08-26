module Configs

export AbstractGridConfig,
    AbstractBathymetryConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    AbstractBoundaryDataConfig,
    AbstractAtmosphereConfig,
    AbstractSimulationConfig,
    AbstractCoupledSimulationConfig,
    AbstractFreeSurfaceConfig,
    AbstractClosureConfig,
    AbstractBoundaryConditionConfig,
    AbstractBoundaryConditionSetConfig,
    AbstractWriterConfig,
    AbstractCallbackConfig,
    AbstractTimeSteppingConfig,
    FjordConfig,
    domain_grid,
    simulation_grid,
    model_closure,
    bathymetry_path,
    forcing_path,
    forcing_directory,
    river_forcing_path,
    boundary_data_path,
    boundary_data_directory,
    atmosphere_path,
    atmosphere_directory,
    results_path,
    plot_path,
    open_edges,
    LATERAL_EDGES,
    validate_open_edge,
    lateral_edges,
    run_tag,
    coverage_window

using Dates: DateTime, Second, format, now

"""
The wall-clock instant this process started, as the fragment `run_tag` names a run's files with.

Filled in by `__init__` rather than by a `const` initialized with `now()`, which would freeze the
moment the package was precompiled into every later run.
"""
const LAUNCH_TAG = Ref("")

__init__() = LAUNCH_TAG[] = format(now(), "yyyymmddTHHMMSS")

"""
    AbstractGridConfig

Supertype for grid configurations. A concrete subtype describes how to build the simulation
grid.

# Methods a subtype provides
- `domain_grid(config, architecture)`: the bare grid, with no bathymetry. Required. Every
  `prepare_*` pipeline builds its target geometry with this.
- `simulation_grid(config, bathymetry_file, architecture)`: the grid the simulation runs on, with
  the processed bathymetry immersed into it. Required. `build_simulation` and
  `prepare_forcing(config::FjordConfig)` both read the model grid back through this.

Both are FjordSim generics rather than methods on `Oceananigans.LatitudeLongitudeGrid` and
`Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid`, which is what the grid hook used to be.
Those names state a *result* type, so a config describing a rectilinear, single-column or
curvilinear grid could not implement the hook under either name without lying about what it
returns. These two name the *role* instead, leaving the return type to the subtype.

`FjordSim.Grids.EvenGrid` is the built-in implementation, and returns a `LatitudeLongitudeGrid`
and an `ImmersedBoundaryGrid` wrapping one.
"""
abstract type AbstractGridConfig end

"""
    domain_grid(config::AbstractGridConfig, architecture)

The bare grid a config describes, with no bathymetry immersed into it.

The geometry every `prepare_*` pipeline regrids onto, and the underlying grid `simulation_grid`
wraps. Declared with no method, so a subtype that does not implement it fails as a `MethodError`
naming the hook rather than silently inheriting some other config's idea of a grid.

Named for the role rather than for a result type — it replaced a method on
`Oceananigans.LatitudeLongitudeGrid`, under which a config describing a rectilinear, single-column
or curvilinear grid could not have implemented the hook without lying about what it returns.

Declared here rather than in `Grids`, which is where its `EvenGrid` method lives: `Bathymetry`,
`Atmospheres` and `Forcing` all call it and are all included *before* `Grids`, so the generic has to
exist in a module they can already see.
"""
function domain_grid end

"""
    simulation_grid(config::AbstractGridConfig, bathymetry_file, architecture)

The grid the simulation runs on: `domain_grid`'s geometry with the processed bathymetry immersed
into it.

Separate from `domain_grid` because the two are needed at different moments and from different
inputs — the prepare steps have no bathymetry file yet, while `build_simulation` must read back the
very grid the bathymetry was written on. Declared here for the same include-order reason.
"""
function simulation_grid end

"""
    AbstractBathymetryConfig

Supertype for bathymetry configurations. A concrete subtype describes one bathymetry source
and the pipeline that turns it into a processed FjordSim bathymetry NetCDF.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `bathymetry_path` and `plot_path`.

# Methods a subtype provides
- `bathymetry_dataset(target_grid, config)`: the NumericalEarth dataset to regrid from.
  Required; `prepare_bathymetry` calls it and everything after it is source-agnostic.
- `regrid_options(config)`: named tuple forwarded to `NumericalEarth.regrid_bathymetry`.
  Optional, defaults to `(;)`.
- `smoothing_options(config)`: named tuple forwarded to `smooth_bathymetry_gaps!`, which is how a
  source opts into the five post-regrid stages — `max_island_cells`, `close_narrow_passages`,
  `spike_ratio`, `minimum_cell_fraction` and `max_slope_factor` (plus `minimum_depth`, which the
  slope limiter must not undo). Optional, defaults to `(;)`, leaving only the topological cleanup
  every source gets.

`FjordSim.Bathymetry.DybdedataConfig` is the built-in implementation, for Geonorge Sjøkart
Dybdedata; `src/Bathymetry/geonorge.jl` is the template to copy for a new source.
"""
abstract type AbstractBathymetryConfig end

"""
    AbstractForcingConfig

Supertype for forcing configurations. A concrete subtype describes one forcing dataset and
how to download and subset it.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `forcing_path` and `plot_path`.
- `output_directory`: resolved by `forcing_directory`; only needed by a dataset that downloads.
- `parameters`: source variable names to prepare.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `interpolation_architecture` to decide
  where `prepare_forcing` runs its interpolation kernel.
- `rivers`: an `AbstractRiverConfig` read by `add_rivers`, or `nothing` for no rivers. Only
  needed by a setup that adds rivers.

The open edge is *not* a field here. It belongs to the open-boundary data config
(`AbstractBoundaryDataConfig`), which hangs off `FjordConfig` rather than off this config, and
reaches `prepare_forcing` as its `edge` keyword — a forcing dataset says nothing about which edge of
the domain is open. It was a field here while the edge carried an interior relaxation band this
config's `relaxation_timescale` set the rate of; that band is gone.

# Methods a subtype provides
- `forcing_time_steps(config)`: the downloaded time records, as `SourceRecord`s. Required.
- `forcing_source_grid(config, filepath)`: geometry of the source data, e.g. a
  `ProjectedSourceGrid`. Required.
- `forcing_variable_names(config)`: source variable name => FjordSim forcing name. Required.
- `download_forcing(target_grid, config)`: fetch the source data. Only if it downloads.
- `simulation_forcing(config, grid, filepath, tracers, reference_date)`: the forcing term object
  `coupled_simulation` consumes, read from the resolved forcing file at `filepath` — which may be
  the rivers-augmented copy rather than `forcing_path(config)` itself, so this takes the path
  rather than resolving it. Optional; defaults to `forcing_from_file`, the FjordSim NetCDF layout
  every built-in source already writes. A source whose prepared files are not that contract
  overrides it.
- `forcing_date_range(config, filepath)`: first and last date the prepared file covers, as a
  `(first, last)` tuple, or `nothing` to skip the check. Optional; defaults to reading `ds["time"]`
  from the prepared NetCDF. The mirror of `atmosphere_date_range`, and a source that overrode
  `simulation_forcing` to read a different layout overrides this too — otherwise
  `validate_time_coverage` would check its coverage with a reader it does not use.

`FjordSim.Forcing.NorKystConfig` is the built-in implementation, for NorKyst-800m;
`src/Forcing/norkyst.jl` is the template to copy for a new dataset.
"""
abstract type AbstractForcingConfig end

"""
    AbstractRiverConfig

Supertype for river configurations. A concrete subtype describes one river dataset and how to
turn it into river relaxation written on top of a forcing file prepared by `prepare_forcing`.

Rivers are optional: a forcing config whose `rivers` field is `nothing` skips the step
entirely.

# Fields a subtype provides
- `data_root`, `output_file`: resolved by `river_forcing_path`, the rivers-augmented copy of
  the forcing file.
- `relaxation_timescale`: seconds; its reciprocal is the lambda written at each river cell.
- `search_radius`: read by the default `river_search_radius`.
- `standalone`: whether `add_rivers` writes a forcing file carrying *only* rivers, instead of
  patching a copy of the one `prepare_forcing` wrote. `true` skips the interior forcing entirely —
  no download, no regrid — and takes the time axis from the run window instead of from the prepared
  file. `false` is today's behaviour, including the error when no prepared forcing exists, so a
  forgotten `prepare_forcing` still fails loudly.

# Methods a subtype provides
- `river_locations(config)`: the river outlets, as `RiverLocation`s. Required.
- `river_series(config, times)`: FjordSim forcing variable name => a `(river, time)` matrix of
  values, one row per `river_locations` entry. Required.
- `download_rivers(config)`: fetch the source data. Only if it downloads.
- `river_search_radius(config)`: how far to look for a coastal cell, in cells. Optional,
  defaults to `config.search_radius`.

`FjordSim.Forcing.OF800RiversConfig` is the built-in implementation, for the OF800 Oslofjord
river dataset; `src/Forcing/of800_rivers.jl` is the template to copy for a new dataset.
"""
abstract type AbstractRiverConfig end

"""
    AbstractBoundaryDataConfig

Supertype for open-boundary *data* configurations. A concrete subtype describes one dataset of
exterior state along a lateral boundary, and how to download and regrid it into the boundary NetCDF
the open lateral boundary condition reads.

Not to be confused with `AbstractBoundaryConditionConfig`, one word away: that one describes the
boundary *condition* — which scheme acts where, and how fast it nudges — while this one describes
the exterior *values* the condition relaxes towards. One is physics, the other is data.

Open-boundary data is optional: a `FjordConfig` whose `boundary_config` is `nothing` skips both
boundary steps entirely, and its lateral boundaries are closed walls — `water_mask` unmasks no
velocity face row, and no boundary condition has an exterior state to nudge towards.

It hangs off `FjordConfig` rather than off the forcing config, and states the open edge itself. The
exterior values were always a separate file read by a separate pipeline; being configured *inside*
the forcing config said otherwise, and meant a setup running with no interior forcing could not have
a data-driven open boundary at all.

The boundary file is separate from the prepared forcing file rather than a part of it, because the
exterior state a Flather/Orlanski boundary needs is **hourly** — a daily mean has the tide averaged
out of it, and elevation with no tide is not worth prescribing. Only the boundary row needs that
cadence, so it is a few hundred MB rather than the hundreds of GB an hourly interior forcing would
be. Every variable in the file is named for its side (`south_T`, `south_eta`, …), so a second open
edge is a set of new variables in the same file rather than a new file or a format change.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `boundary_data_path` and `plot_path`.
- `output_directory`: resolved by `boundary_data_directory`; the directory the download writes its
  intermediate files into.
- `open_edges`: which lateral boundaries the domain is open on, each one of `:south`, `:north`,
  `:west` or `:east`, read by `open_edges` as a `Vector{Symbol}`. A domain in the open ocean names all
  four; a fjord names the one it connects to the sea through. Stated once, here:
  `download_boundaries` and `prepare_boundaries` read it off this config rather than taking it beside
  it, `water_mask` unmasks each of those edges' velocity face rows for both this pipeline and the
  interior forcing's, and the open lateral boundary condition puts its schemes on every one of them.
- `parameters`: source variable names to prepare.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `interpolation_architecture`.

# Methods a subtype provides
- `boundary_time_steps(config)`: the downloaded time records, as `SourceRecord`s. Required.
- `boundary_source_grid(config, filepath)`: geometry of the source data, e.g. a
  `ProjectedSourceGrid`. Required.
- `boundary_variable_names(config)`: source variable name => FjordSim boundary name. Required.
- `download_boundaries(target_grid, config)`: fetch the source data, reading the edges it subsets
  along with `open_edges(config)`. Only if it downloads.
- `boundary_date_range(config)`: first and last date the prepared file covers, as a
  `(first, last)` tuple, or `nothing` to skip the check. Optional; defaults to reading `ds["time"]`
  from the prepared NetCDF, so `validate_time_coverage` also guards the `Cyclical()` wrap here.
- `boundary_source_slab(config, reader, step, source_name)`: one source slab as the writer
  interpolates from. Optional; defaults to `blended_slab`, a plain read. Override it only for a
  variable that must be *derived* from more than one of the source's own — which a velocity component
  is, when the source states it along its own grid axes rather than eastward/northward.

`FjordSim.Forcing.NorKystBoundariesConfig` is the built-in implementation, for MET Norway's hourly
NorKyst-800m collection; `src/Forcing/norkyst_boundaries.jl` is the template to copy for a new
dataset.
"""
abstract type AbstractBoundaryDataConfig end

"""
    AbstractAtmosphereConfig

Supertype for atmosphere configurations. A concrete subtype describes one atmospheric
reanalysis dataset, how to download it and how to regrid it onto the regular longitude/latitude
grid the simulation-time reader consumes.

The prepared file is what `FjordSim.Atmospheres.NORA3.MultiYearNORA3` reads, so its layout is a
fixed contract rather than a dataset detail: variables of shape `(longitude, latitude, time)`,
uniformly spaced 1D `lon` and `lat` cell centers, and a CF-encoded `time`. `prepare_atmosphere`
writes it; a subtype only supplies the dataset-specific parts.

Unlike the forcing grid, the atmosphere grid is independent of the ocean grid — NumericalEarth
interpolates between them — so it only has to cover the ocean domain with a margin.

# Fields a subtype provides
- `data_root`, `output_file`, `plot_file`: resolved by `atmosphere_path` and `plot_path`.
- `output_directory`: resolved by `atmosphere_directory`; the directory the download writes its
  intermediate files into. Only needed by a dataset that downloads.
- `resolution`: target grid spacing in degrees.
- `padding`: degrees of margin added around the ocean domain, read by
  `atmosphere_target_axes`.

# Methods a subtype provides
- `atmosphere_time_steps(config)`: the downloaded time records, as `AtmosphereRecord`s.
  Required.
- `atmosphere_source_grid(config, filepath)`: geometry of the downloaded data, e.g. a
  `ProjectedAtmosphereGrid`. Required.
- `atmosphere_variable_names(config)`: source variable name => prepared variable name.
  Required.
- `download_atmosphere(target_grid, config)`: fetch the source data. Only if it downloads.
- `prescribed_atmosphere(config, architecture; reference_date)`,
  `prescribed_radiation(config, architecture; reference_date)`: read the prepared file back at
  simulation time, as the NumericalEarth objects `coupled_simulation` consumes.
  Only if the setup is simulated.
- `atmosphere_date_range(config)`: first and last date of the prepared file, for
  `build_simulation`'s coverage check. Optional; the default is `nothing`, which skips the check.

`FjordSim.Atmospheres.NORA3Config` is the built-in implementation, for the MET Norway NORA3
reanalysis; `src/Atmospheres/nora3_source.jl` is the template to copy for a new dataset.
"""
abstract type AbstractAtmosphereConfig end

"""
    AbstractSimulationConfig

Supertype for simulation configurations. A concrete subtype describes how a setup whose data is
already prepared is turned into a running simulation: the model components that do not depend on
the grid, and the knobs for the ones that do.

Unlike the other supertypes this one declares **fields only, no hooks**. `build_simulation` is
generic over the whole `FjordConfig` rather than over a simulation source, because everything
dataset-specific — which forcing file, which atmosphere reader — already comes from the other
configs. What *is* dispatched on lives one level down, in the five nested configs below, each of
which has its own supertype and its own hook.

# Fields a subtype provides
- `results_root`: the directory the run writes into, and what a writer's `output_file` and
  `FjordSim.CLI.log_path` resolve against. The one config rooted at a results directory rather
  than a `data_root`, being the only one that writes rather than reads.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `simulation_architecture`.
- `model`: an `AbstractCoupledSimulationConfig` — what the model *is*, assembled by
  `coupled_simulation`.
- `boundary_conditions`: an `AbstractBoundaryConditionSetConfig`, materialized by
  `FjordSim.BoundaryConditions.field_boundary_conditions`.
- `writers`: a tuple of `AbstractWriterConfig`s, attached by `attach_writer!`.
- `callbacks`: a tuple of `AbstractCallbackConfig`s, attached by `attach_callback!`.
- `time_stepping`: an `AbstractTimeSteppingConfig`, attached by `attach_time_stepping!`.
- `initial_conditions`: where the ocean state starts from, resolved by
  `FjordSim.Simulations.resolve_initial_conditions` before it reaches `coupled_simulation`.
- `start_date`: the calendar instant model time zero stands for. Read by `coverage_window` as well
  as by `build_simulation`.
- `stop_time`, `loops`, `pickup`: run control, read by `build_simulation` and `run_simulation`.

`results_root`, `start_date` and `stop_time` are the supertype-level minimum: the first is what
`FjordSim.CLI.log_path` and `results_path` root against, and the other two are what
`coverage_window` reports. A subtype omitting one inherits neither helper.

`FjordSim.Simulations.SimulationConfig` is the built-in implementation.
"""
abstract type AbstractSimulationConfig end

"""
    AbstractCoupledSimulationConfig

Supertype for coupled-model configurations. A concrete subtype names the model components that do
not depend on the grid, and its `coupled_simulation` method assembles them into a `Simulation`.

This is the seam that lets a setup run something other than a hydrostatic free-surface ocean inside
an `OceanSeaIceModel`: `build_simulation` derives the grid, the forcing, the boundary conditions,
the initial state and the atmosphere and then hands all of it to one generic call, so a different
model assembly is a new subtype rather than an edit to an existing method body.

# Methods a subtype provides
- `coupled_simulation(model, grid; forcing, boundary_conditions, initial_conditions, atmosphere,
  radiation, boundary_config, stop_time, initial_time_step)`: build and return the coupled
  `Simulation`. Required. `boundary_config` is the setup's `AbstractBoundaryDataConfig` (or
  `nothing`), passed for the same reason `boundary_condition_sides` receives it — so a nested config
  that acts on the open edges reads them from the one place they are stated.
- `model_tracers(model)`: the tracer names the model carries, as a tuple of `Symbol`s. Required —
  `build_simulation` needs them before the model exists, to pick which forcing terms to build,
  which lateral tracer boundaries to open and which state variables to read back.

`FjordSim.Simulations.CoupledHydrostaticSimulation` is the built-in implementation.
"""
abstract type AbstractCoupledSimulationConfig end

"""
    AbstractFreeSurfaceConfig

Supertype for free-surface configurations. A concrete subtype names a free-surface solver's knobs,
and its `free_surface(config, grid)` method builds the object `coupled_simulation` passes to
`HydrostaticFreeSurfaceModel`.

Nested inside a coupled-simulation config as its own field rather than a plain scalar, for the same
reason the model itself is a config rather than the built object: `SplitExplicitFreeSurface(grid,
cfl = ...)` needs the grid, which does not exist until `coupled_simulation` is called.

# Methods a subtype provides
- `free_surface(config, grid)`: build the free-surface object. Required.

`FjordSim.Simulations.SplitExplicitFreeSurfaceConfig` is the built-in implementation.
"""
abstract type AbstractFreeSurfaceConfig end

"""
    AbstractClosureConfig

Supertype for *closure* configurations. A concrete subtype names a turbulence closure's knobs, and
its `model_closure(config, grid, boundary_config)` method builds the closure — or closure tuple —
that `coupled_simulation` passes to `HydrostaticFreeSurfaceModel`.

The `closure` field of a coupled-simulation config holds either a pre-built Oceananigans closure
(the usual case, and what every setup wrote before this supertype existed) or one of these. Which
one it is is resolved by dispatch on `model_closure`, the way `initial_conditions` is resolved by
`FjordSim.Simulations.resolve_initial_conditions` — three methods rather than a branch, and an
identity fallback so a pre-built closure passes through untouched and no existing setup changes.

It exists because a closure can need what only exists at build time. `FjordSim.Simulations.BoundarySponge`
is the built-in implementation and needs two such things: the grid, for the domain extent and cell
size its ramp is expressed in, and the *open edges*, which are stated once on the
`AbstractBoundaryDataConfig` and read back with `open_edges` — so a setup that opens a second edge
cannot end up sponging only the first.

# Methods a subtype provides
- `model_closure(config, grid, boundary_config)`: build the closure. Required.
"""
abstract type AbstractClosureConfig end

"""
    model_closure(closure, grid, boundary_config)

The closure `coupled_simulation` gives `HydrostaticFreeSurfaceModel`.

The fallback is the identity, so anything a setup already writes — a single Oceananigans closure, or
a tuple of them — passes through unchanged. An `AbstractClosureConfig` overloads it to build one
from the grid and the setup's open-boundary data config.
"""
model_closure(closure, grid, boundary_config) = closure

"""
    AbstractBoundaryConditionConfig

Supertype for boundary-condition configurations. A concrete subtype describes one *contribution* to
the model's boundary conditions, and its `boundary_conditions` method builds that contribution.

A setup composes its boundaries out of several of these — air-sea fluxes, bottom drag and an open
lateral edge are three independent statements — and opts out of a piece by leaving it out. They are
collected by an `AbstractBoundaryConditionSetConfig`, whose built-in `MergedBoundaryConditions`
merges them with `FjordSim.Utils.recursive_merge` and materializes the result. So each method here
returns the *nested named tuple* shape (`(u = (top = …, bottom = …), T = (top = …,))`), not
`FieldBoundaryConditions`.

# Methods a subtype provides
- `boundary_condition_sides(config, grid, forcing, boundary_config, tracers[, boundaries])`: the
  nested named tuple this piece contributes, keyed by field and then by side. Required. Implement the
  six-argument form only if the piece reads the prepared open-boundary state (`boundaries`, from
  `FjordSim.Forcing.boundary_series`, or `nothing`); the five-argument form is reached through a
  forwarding fallback, which is what lets a piece that carries its own exterior state ignore it.

`FjordSim.BoundaryConditions.AirSeaFluxes`, `FjordSim.BoundaryConditions.QuadraticBottomDrag` and
`FjordSim.BoundaryConditions.OpenLateralBoundaryFromData` are the built-in implementations. Not to
be confused with `AbstractBoundaryDataConfig`, which is the exterior *data* rather than the scheme —
and which is what the fourth argument holds, so a piece can read `open_edges` from it.
"""
abstract type AbstractBoundaryConditionConfig end

"""
    AbstractBoundaryConditionSetConfig

Supertype for boundary-condition *set* configurations. A concrete subtype collects the individual
`AbstractBoundaryConditionConfig` contributions and materializes them into the object the model
consumes.

Two supertypes rather than one because the two answer different questions, and only the second was
previously dispatched on anything: a contribution says *what* one piece of the boundary is, while a
set says *how* the pieces combine and *what shape* the result takes. That second question used to be
answered by `field_boundary_conditions(configs::Tuple, …)` — dispatch on a bare `Tuple`, which no
subtype can specialize and which states nothing about what the tuple means. A model whose boundary
objects are not Oceananigans `FieldBoundaryConditions`, or that wants something other than
last-wins merging, is now a new subtype here rather than an edit to that function.

# Methods a subtype provides
- `field_boundary_conditions(config, grid, forcing, boundary_config, tracers, boundaries)`: the materialized
  boundary conditions, as the `NamedTuple` the model constructor consumes. Required.

`FjordSim.BoundaryConditions.MergedBoundaryConditions` is the built-in implementation.
"""
abstract type AbstractBoundaryConditionSetConfig end

"""
    AbstractWriterConfig

Supertype for output-writer configurations. A concrete subtype describes one thing a run writes, and
its `attach_writer!` method attaches it to the simulation.

A tuple of these replaces the hardcoded snapshot writer and the `checkpoint_interval` threshold that
used to decide whether a `Checkpointer` existed: a setup that wants no checkpoints names no
checkpointing writer, and a setup that wants two snapshot files at different intervals names two
writers. Which simulation a writer attaches to — the coupled one or the ocean sub-simulation — is
the method's business, not the caller's.

# Fields a subtype provides
- `output_file`: resolved by `results_path`, **only** for a writer whose `output_path_trait` is
  `FjordSim.Simulations.NamesOutputFile`. A writer that names no product file (a checkpointer,
  whose files carry no run tag and are pruned) needs no such field, and `results_path` reports that
  as an error rather than reaching for a field that is not there.

# Methods a subtype provides
- `attach_writer!(simulation, writer, config, loop)`: attach this writer for repetition `loop`.
  Required.
- `FjordSim.Simulations.checkpoint_trait(writer)`: whether this writer contributes a `Checkpointer`.
  Defaults to `NotCheckpointing()`.
- `FjordSim.Simulations.output_path_trait(writer)`: whether this writer names a file the run should
  report. Defaults to `NamesNoOutputFile()`.

`FjordSim.Simulations.SnapshotWriter` and `FjordSim.Simulations.CheckpointWriter` are the built-in
implementations.
"""
abstract type AbstractWriterConfig end

"""
    AbstractCallbackConfig

Supertype for callback configurations. A concrete subtype describes one thing a run does *while* it
runs — report progress, accumulate a diagnostic, watch for a stopping condition — and its
`attach_callback!` method installs it on the simulation.

The same shape as `AbstractWriterConfig`, and for the same reason: what a run reports used to be a
single hardcoded `Callback(progress, TimeInterval(config.progress_interval))` under the fixed key
`:progress`, so changing the report, the schedule type, or adding a second diagnostic all meant
rewriting `build_simulation`. A tuple of these makes each of those a config choice, and a setup that
wants a silent run names an empty tuple.

# Fields a subtype provides
- `name`: the key it occupies in `simulation.callbacks`. Two callbacks under one name replace each
  other, which `FjordSim.Simulations.validate_callbacks` rejects.

# Methods a subtype provides
- `attach_callback!(simulation, callback, config)`: attach this callback. Required. Which simulation
  it attaches to — the coupled one or the ocean sub-simulation — is the method's business, exactly
  as it is for a writer.

`FjordSim.Simulations.ProgressCallback` is the built-in implementation.
"""
abstract type AbstractCallbackConfig end

"""
    AbstractTimeSteppingConfig

Supertype for time-stepping configurations. A concrete subtype describes how the model's time step
is chosen, and its `attach_time_stepping!` method installs that choice on the simulation.

# Methods a subtype provides
- `attach_time_stepping!(simulation, config)`: install the time-step policy. Required. It must not
  set `simulation.Δt` — the step a loop starts from is `initial_time_step`, applied once when the
  simulation is built, and `FjordSim.Simulations.restart_loop!` deliberately leaves `Δt` alone so a
  repetition keeps the step the previous one converged on.
- `initial_time_step(config)`: the `Δt` the simulation starts at, in seconds. Required.

`FjordSim.Simulations.AdaptiveTimeStep` is the built-in implementation.
"""
abstract type AbstractTimeSteppingConfig end

"""
    bathymetry_path(config)
    forcing_path(config)
    atmosphere_path(config)
    plot_path(config)

Resolve `config.output_file` and `config.plot_file` against `config.data_root`. Defined for
every `AbstractBathymetryConfig`, `AbstractForcingConfig` and `AbstractAtmosphereConfig`, so a
new bathymetry source, forcing dataset or atmosphere dataset inherits path resolution. A field
holding an absolute path is returned unchanged, relocating that file outside `data_root`.
"""
bathymetry_path(config::AbstractBathymetryConfig) = joinpath(config.data_root, config.output_file)
forcing_path(config::AbstractForcingConfig) = joinpath(config.data_root, config.output_file)
atmosphere_path(config::AbstractAtmosphereConfig) = joinpath(config.data_root, config.output_file)
plot_path(config::AbstractBathymetryConfig) = joinpath(config.data_root, config.plot_file)
plot_path(config::AbstractForcingConfig) = joinpath(config.data_root, config.plot_file)
plot_path(config::AbstractBoundaryDataConfig) = joinpath(config.data_root, config.plot_file)
plot_path(config::AbstractAtmosphereConfig) = joinpath(config.data_root, config.plot_file)

"""
    forcing_directory(config)

Resolve `config.output_directory` against `config.data_root`: the directory a forcing
dataset downloads its source files into. An absolute `output_directory` is returned
unchanged.
"""
forcing_directory(config::AbstractForcingConfig) = joinpath(config.data_root, config.output_directory)

"""
    river_forcing_path(config)

Resolve `config.output_file` against `config.data_root`: the rivers-augmented copy of the
forcing file written by `add_rivers`. Defined for every `AbstractRiverConfig`, so a new river
dataset inherits path resolution. An absolute `output_file` is returned unchanged.
"""
river_forcing_path(config::AbstractRiverConfig) = joinpath(config.data_root, config.output_file)

"""
    boundary_data_path(config)
    boundary_data_directory(config)

Resolve `config.output_file` and `config.output_directory` against `config.data_root`: the prepared
open-boundary NetCDF written by `prepare_boundaries`, and the directory `download_boundaries` writes
its intermediate files into. Defined for every `AbstractBoundaryDataConfig`, so a new dataset
inherits path resolution. An absolute field is returned unchanged.

The edge is *not* in the filename: every variable in the file is prefixed with its side, so one file
holds however many boundaries a setup opens.
"""
boundary_data_path(config::AbstractBoundaryDataConfig) = joinpath(config.data_root, config.output_file)

"""
The four lateral boundaries a regional domain can be open on, in the order every error message
lists them.

Stated here, beside `open_edges`, because every module that reasons about an edge is included after
this one and none of them can see each other: `Bathymetry` clears land along the open ones, `Forcing`
restores their velocity face rows and regrids along them, and every `Val{edge}` dispatch in
`BoundaryConditions` falls back to an `ArgumentError` naming this tuple. It has already been a
per-module copy twice — with identical contents, an identical membership test and an identical error
string each time — which is what settled the question of where it belongs.
"""
const LATERAL_EDGES = (:south, :north, :west, :east)

"""
    validate_open_edge(edge)

Return `edge` if it names a lateral boundary, else throw. Checked up front so a typo fails
before any regridding rather than deep inside the variable loop.

Checked here rather than in each pipeline so one typo produces one error message wherever it is
written — in a boundary data config's `open_edges`, or in the `edges` a prepare step is handed.
"""
function validate_open_edge(edge)
    edge in LATERAL_EDGES ||
        throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))
    return edge
end

"""
    lateral_edges(edges)

Normalize whatever a caller names its open edges as into a validated `Vector{Symbol}`: one `Symbol`,
an iterable of them, or `nothing` for none.

Every consumer iterates, so a domain open on one edge, on all four, or on none is the same code path
rather than a scalar case and a plural one. Duplicates are rejected — a repeated edge would prepare
and write the same variables twice, and mean nothing the second time.
"""
lateral_edges(::Nothing) = Symbol[]
lateral_edges(edge::Symbol) = [validate_open_edge(edge)]

function lateral_edges(edges)
    normalized = Symbol[validate_open_edge(edge) for edge in edges]
    allunique(normalized) ||
        throw(ArgumentError("open_edges names an edge more than once: $normalized"))
    return normalized
end

"""
    open_edges(config)

Which lateral boundaries the domain is open on, as named by an `AbstractBoundaryDataConfig` — a
`Vector{Symbol}`, empty for a setup naming no open-boundary data, whose lateral boundaries are
therefore all closed walls.

Always a collection, never a bare `Symbol`, so a domain open on one edge, on all four, or on none is
one code path in every consumer rather than a scalar case and a plural one. A built-in config's
keyword constructor accepts either shape and normalizes with `lateral_edges`.

Declared here rather than beside the pipelines that read it because all of them do: `Forcing` needs
it for `water_mask` and is included first, both boundary steps regrid and subset along the edges, and
the open lateral boundary condition puts its schemes on each. Each of them reads it through this
accessor, so the edges are stated in exactly one place — the config that describes the boundary.
"""
open_edges(::Nothing) = Symbol[]
open_edges(config::AbstractBoundaryDataConfig) = config.open_edges
boundary_data_directory(config::AbstractBoundaryDataConfig) =
    joinpath(config.data_root, config.output_directory)

"""
    atmosphere_directory(config)

Resolve `config.output_directory` against `config.data_root`: the directory an atmosphere
dataset downloads its intermediate files into. An absolute `output_directory` is returned
unchanged.
"""
atmosphere_directory(config::AbstractAtmosphereConfig) = joinpath(config.data_root, config.output_directory)

"""
    run_tag(config)

The run's identity as a filename fragment: the wall-clock instant this process started, as
`yyyymmddTHHMMSS`.

The launch instant rather than the simulated `start_date`, so it distinguishes *invocations* and a
re-run cannot overwrite an earlier one's output. Which simulated window a file covers is recorded in
the snapshot's `start_date` global attribute rather than in its name. Seconds are in the format so a
crash and an immediate relaunch do not collide.

Constant for the life of the process, so `results_path`, `FjordSim.CLI.log_path` and the checkpointer
all agree however often they ask. `config` is taken only so a simulation-config subtype can name its
runs differently.
"""
run_tag(::AbstractSimulationConfig) = LAUNCH_TAG[]

"""
    results_path(writer, config)
    results_path(writer, config, loop)

Resolve `writer.output_file` against `config.results_root`: the file that writer writes to. Results
are rooted separately from the input data, which is why this reads `results_root` rather than
`data_root`. An absolute `output_file` keeps its own directory.

`run_tag(config)` is inserted before the extension so separate runs do not overwrite each other, and
the three-argument form appends the loop index on top of it, which is how a looped run gives each
repetition its own file.

The filename is the *writer's* rather than the simulation config's, so a setup can name several
outputs at different intervals. A writer that names no product file has no `output_file` and gets an
error saying so, rather than a `getproperty` failure — see `AbstractWriterConfig`.
"""
results_path(writer::AbstractWriterConfig, config::AbstractSimulationConfig) =
    tagged_path(writer, config, run_tag(config))

results_path(writer::AbstractWriterConfig, config::AbstractSimulationConfig, loop::Int) =
    tagged_path(writer, config, string(run_tag(config), "_loop", lpad(loop, 2, '0')))

"""
    tagged_path(writer, config, tag)

`writer.output_file` with `_tag` inserted before its extension, resolved against `results_root`
unless it is absolute — so an absolute `output_file` still relocates just that file, and the tag
lands on the basename either way.
"""
function tagged_path(writer::AbstractWriterConfig, config::AbstractSimulationConfig, tag)
    hasproperty(writer, :output_file) || throw(
        ArgumentError(
            "$(typeof(writer)) names no `output_file`, so it has no `results_path`. Only a writer " *
            "whose `output_path_trait` is `NamesOutputFile` resolves to one.",
        ),
    )

    directory, file = splitdir(writer.output_file)
    stem, extension = splitext(file)
    tagged = string(stem, "_", tag, extension)
    isabspath(writer.output_file) && return joinpath(directory, tagged)
    return joinpath(config.results_root, directory, tagged)
end

"""
    coverage_window(config)

The calendar interval a run needs its prepared inputs to span, as `(first, last)` — or `nothing`
for a setup naming no simulation config.

Passed to `prepare_forcing` and `prepare_atmosphere` as their `coverage`, so each pads its own time
axis to reach both ends. A plain tuple of dates rather than the config itself, so neither pipeline
ever sees a simulation config.

A looped run reuses one window rather than extending it, so this is `stop_time` from `start_date`
and not `loops * stop_time`.
"""
coverage_window(::Nothing) = nothing

coverage_window(config::AbstractSimulationConfig) =
    (config.start_date, config.start_date + Second(round(Int, config.stop_time)))

"""
    FjordConfig

Complete setup configuration for one fjord, as returned by the setup functions in
`FjordSim.Setups`.

The fields are parametric over the abstract config supertypes rather than fixed to the
built-in types, so a new grid, bathymetry source or forcing dataset only needs a struct
subtyping the matching supertype plus methods overloaded on it — `FjordConfig` and the
functions taking it are untouched. Each instantiation still has concrete field types.

# Fields
- `grid_config`: any `AbstractGridConfig`.
- `bathymetry_config`: any `AbstractBathymetryConfig`.
- `forcing_config`: any `AbstractForcingConfig`, or `nothing` for a setup that runs with no
  forcing at all. Defaults to `nothing`, so a setup opts in by naming one.
- `boundary_config`: any `AbstractBoundaryDataConfig`, or `nothing` for a setup whose lateral
  boundaries are closed. Defaults to `nothing`, so a setup opts in by naming one. Independent of
  `forcing_config`: the exterior state along the open edge is its own hourly file from its own
  pipeline, and either config can be named without the other.
- `atmosphere_config`: any `AbstractAtmosphereConfig`, or `nothing` for a setup that prepares no
  atmosphere. Defaults to `nothing`, so a setup opts in by naming one.
- `simulation_config`: any `AbstractSimulationConfig`, or `nothing` for a setup that is only
  prepared and not run. Defaults to `nothing`, so a setup opts in by naming one.

# Example

```julia
struct SingleColumn <: AbstractGridConfig
    depth::Float64
end

FjordSim.domain_grid(config::SingleColumn, architecture) = ...      # new behavior
FjordSim.simulation_grid(config::SingleColumn, file, architecture) = ...
```
"""
Base.@kwdef mutable struct FjordConfig{
    G<:AbstractGridConfig,
    B<:AbstractBathymetryConfig,
    F,
    O,
    A,
    S,
}
    grid_config::G
    bathymetry_config::B
    forcing_config::F = nothing
    boundary_config::O = nothing
    atmosphere_config::A = nothing
    simulation_config::S = nothing
end

end  # module Configs
