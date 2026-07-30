module Configs

export AbstractGridConfig, AbstractBathymetryConfig, AbstractForcingConfig, FjordConfig

"""
    AbstractGridConfig

Supertype for grid configurations. A concrete subtype describes how to build the simulation
grid; behavior is added by overloading `LatitudeLongitudeGrid(arch, config)` on it.

`FjordSim.Grids.EvenGrid` is the built-in implementation.
"""
abstract type AbstractGridConfig end

"""
    AbstractBathymetryConfig

Supertype for bathymetry configurations. A concrete subtype describes one bathymetry source
and the pipeline that turns it into a processed FjordSim bathymetry NetCDF.

Subtypes are expected to provide `data_root`, `output_file` and `plot_file`, which
`bathymetry_path` and `plot_path` resolve. Source-specific behavior is added by overloading
on the concrete subtype.

`FjordSim.Bathymetry.DybdedataConfig` is the built-in implementation, for Geonorge Sjøkart
Dybdedata.
"""
abstract type AbstractBathymetryConfig end

"""
    AbstractForcingConfig

Supertype for forcing configurations. A concrete subtype describes one forcing dataset and
how to download and subset it.

`FjordSim.Forcing.NorKystConfig` is the built-in implementation, for NorKyst-800m.
"""
abstract type AbstractForcingConfig end

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
