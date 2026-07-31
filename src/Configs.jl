module Configs

export AbstractGridConfig,
    AbstractBathymetryConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    FjordConfig,
    bathymetry_path,
    forcing_path,
    forcing_directory,
    river_forcing_path,
    plot_path

"""
    AbstractGridConfig

Supertype for grid configurations. A concrete subtype describes how to build the simulation
grid.

# Methods a subtype provides
- `LatitudeLongitudeGrid(arch, config)`: build the grid. Required.

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
Dybdedata; `src/Bathymetry/Geonorge.jl` is the template to copy for a new source.
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
`src/Forcing/NorKyst.jl` is the template to copy for a new dataset.
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
river dataset; `src/Forcing/OF800Rivers.jl` is the template to copy for a new dataset.
"""
abstract type AbstractRiverConfig end

"""
    bathymetry_path(config)
    forcing_path(config)
    plot_path(config)

Resolve `config.output_file` and `config.plot_file` against `config.data_root`. Defined for
every `AbstractBathymetryConfig` and `AbstractForcingConfig`, so a new bathymetry source or
forcing dataset inherits path resolution. A field holding an absolute path is returned
unchanged, relocating that file outside `data_root`.
"""
bathymetry_path(config::AbstractBathymetryConfig) = joinpath(config.data_root, config.output_file)
forcing_path(config::AbstractForcingConfig) = joinpath(config.data_root, config.output_file)
plot_path(config::AbstractBathymetryConfig) = joinpath(config.data_root, config.plot_file)
plot_path(config::AbstractForcingConfig) = joinpath(config.data_root, config.plot_file)

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
    FjordConfig

Complete setup configuration for one fjord, as returned by the files in `configs/`.

The fields are parametric over the abstract config supertypes rather than fixed to the
built-in types, so a new grid, bathymetry source or forcing dataset only needs a struct
subtyping the matching supertype plus methods overloaded on it — `FjordConfig` and the
functions taking it are untouched. Each instantiation still has concrete field types.

# Fields
- `grid_config`: any `AbstractGridConfig`.
- `bathymetry_config`: any `AbstractBathymetryConfig`.
- `forcing_config`: any `AbstractForcingConfig`.

# Example

```julia
struct SingleColumn <: AbstractGridConfig
    depth::Float64
end

Oceananigans.LatitudeLongitudeGrid(arch, config::SingleColumn) = ...  # new behavior
```
"""
Base.@kwdef mutable struct FjordConfig{
    G<:AbstractGridConfig,
    B<:AbstractBathymetryConfig,
    F<:AbstractForcingConfig,
}
    grid_config::G
    bathymetry_config::B
    forcing_config::F
end

end  # module Configs
