module Configs

export AbstractGridConfig,
    AbstractBathymetryConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    AbstractAtmosphereConfig,
    AbstractSimulationConfig,
    FjordConfig,
    bathymetry_path,
    forcing_path,
    forcing_directory,
    river_forcing_path,
    atmosphere_path,
    atmosphere_directory,
    results_path,
    plot_path,
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
- `LatitudeLongitudeGrid(architecture, config)`: build the grid. Required.

`FjordSim.Grids.EvenGrid` is the built-in implementation.
"""
abstract type AbstractGridConfig end

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
- `relaxation_edge`, `relaxation_cells`, `relaxation_timescale`: read by the generic
  `water_mask` and `relaxation_lambda`.
- `parameters`: source variable names to prepare.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `interpolation_architecture` to decide
  where `prepare_forcing` runs its interpolation kernel.
- `rivers`: an `AbstractRiverConfig` read by `add_rivers`, or `nothing` for no rivers. Only
  needed by a setup that adds rivers.

# Methods a subtype provides
- `forcing_time_steps(config)`: the downloaded time records, as `SourceRecord`s. Required.
- `forcing_source_grid(config, filepath)`: geometry of the source data, e.g. a
  `ProjectedSourceGrid`. Required.
- `forcing_variable_names(config)`: source variable name => FjordSim forcing name. Required.
- `download_forcing(target_grid, config)`: fetch the source data. Only if it downloads.

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
  simulation time, as the NumericalEarth objects `coupled_hydrostatic_simulation` consumes.
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
configs.

# Fields a subtype provides
- `results_root`, `output_file`: resolved by `results_path`, the simulation's output file.
- `architecture`: `:auto`, `:cpu` or `:gpu`, resolved by `simulation_architecture`.
- `buoyancy`, `closure`, `tracer_advection`, `momentum_advection`, `tracers`, `coriolis`,
  `sea_ice`, `biogeochemistry`: passed to `coupled_hydrostatic_simulation` as they are.
- `initial_conditions`: where the ocean state starts from, resolved by
  `FjordSim.Simulations.resolve_initial_conditions` before it reaches
  `coupled_hydrostatic_simulation`.
- `free_surface_cfl`, `bottom_drag_coefficient`: the grid-dependent components are built from
  these, since neither can exist before the grid does.
- `start_date`: the calendar instant model time zero stands for. Read by `coverage_window` as well
  as by `build_simulation`.
- `stop_time`, `loops`, `output_interval`, `progress_interval`, `overwrite_existing`,
  `checkpoint_interval`, `pickup`, `time_step_cfl`, `max_time_step`, `max_time_step_change`: run
  control, read by `build_simulation` and `run_simulation`.

`FjordSim.Simulations.SimulationConfig` is the built-in implementation.
"""
abstract type AbstractSimulationConfig end

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
    results_path(config)
    results_path(config, loop)

Resolve `config.output_file` against `config.results_root`: the file the simulation writes its
snapshots to. Results are rooted separately from the input data, which is why this reads
`results_root` rather than `data_root`. An absolute `output_file` keeps its own directory.

`run_tag(config)` is inserted before the extension so separate runs do not overwrite each other, and
the two-argument form appends the loop index on top of it, which is how a looped run gives each
repetition its own file.
"""
results_path(config::AbstractSimulationConfig) = tagged_path(config, run_tag(config))

results_path(config::AbstractSimulationConfig, loop::Int) =
    tagged_path(config, string(run_tag(config), "_loop", lpad(loop, 2, '0')))

"""
    tagged_path(config, tag)

`config.output_file` with `_tag` inserted before its extension, resolved against `results_root`
unless it is absolute — so an absolute `output_file` still relocates just that file, and the tag
lands on the basename either way.
"""
function tagged_path(config::AbstractSimulationConfig, tag)
    directory, file = splitdir(config.output_file)
    stem, extension = splitext(file)
    tagged = string(stem, "_", tag, extension)
    isabspath(config.output_file) && return joinpath(directory, tagged)
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
- `forcing_config`: any `AbstractForcingConfig`.
- `atmosphere_config`: any `AbstractAtmosphereConfig`, or `nothing` for a setup that prepares no
  atmosphere. Defaults to `nothing`, so a setup opts in by naming one.
- `simulation_config`: any `AbstractSimulationConfig`, or `nothing` for a setup that is only
  prepared and not run. Defaults to `nothing`, so a setup opts in by naming one.

# Example

```julia
struct SingleColumn <: AbstractGridConfig
    depth::Float64
end

Oceananigans.LatitudeLongitudeGrid(architecture, config::SingleColumn) = ...  # new behavior
```
"""
Base.@kwdef mutable struct FjordConfig{
    G<:AbstractGridConfig,
    B<:AbstractBathymetryConfig,
    F<:AbstractForcingConfig,
    A,
    S,
}
    grid_config::G
    bathymetry_config::B
    forcing_config::F
    atmosphere_config::A = nothing
    simulation_config::S = nothing
end

end  # module Configs
