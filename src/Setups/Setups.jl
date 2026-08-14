module Setups

export fjord_config, setup_names, oslofjorden, drammensfjorden

import Oceananigans        # for `Oceananigans.defaults.FloatType`
using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: day, days, hour, minutes, second
using Dates: DateTime
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature

using ..Configs: FjordConfig
using ..Utils: progress
using ..Grids: EvenGrid
using ..Bathymetry: DybdedataConfig
using ..Atmospheres: NORA3Config
using ..Forcing: NorKystConfig, OF800RiversConfig, NorKystBoundariesConfig
using ..BoundaryConditions:
    AirSeaFluxes, QuadraticBottomDrag, OpenLateralBoundaryFromData, MergedBoundaryConditions
using ..Simulations:
    SimulationConfig,
    CoupledHydrostaticSimulation,
    SplitExplicitFreeSurfaceConfig,
    SnapshotWriter,
    CheckpointWriter,
    ProgressCallback,
    AdaptiveTimeStep,
    FromForcing,
    FromResults

include("oslofjorden.jl")
include("drammensfjorden.jl")

# Each setup is a function rather than a `const FjordConfig`, and building one at load time would
# break in two ways that are silent rather than loud:
#
# - `DybdedataConfig`'s `raw_directory` defaults to the scratch path that `Bathymetry`'s `__init__`
#   fills in, which is still `""` while the package precompiles.
# - The config structs are mutable and `native_region!` mutates the bathymetry config, so a shared
#   instance would leak the first subcommand's state into the next.
#
# `Function` as the value type is deliberate: `--config` arrives as a runtime string, so a `Val`
# would only move the dynamic dispatch up a level. The `driver(config::FjordConfig)` methods are
# the function barrier.
const SETUPS = Dict{String,Function}(
    "oslofjorden" => oslofjorden,
    "drammensfjorden" => drammensfjorden,
)

"""
    setup_names()

The registered setup names, sorted. `SETUPS` is a `Dict`, whose iteration order is unspecified,
so both the help text and the tests read this instead of the keys directly.
"""
setup_names() = sort!(collect(keys(SETUPS)))

"""
    fjord_config(name_or_path)

The `FjordConfig` for a registered setup name, e.g. `fjord_config("oslofjorden")`.

An argument ending in `.jl` is instead loaded as an out-of-tree config file whose last expression
is a `FjordConfig`, so a fjord can be defined without adding it to the package. Anything else is
looked up in `SETUPS` and errors listing `setup_names()` — a misspelled name must not fall through
to the file branch and report a missing file.
"""
function fjord_config(name_or_path)
    endswith(name_or_path, ".jl") && return load_config_file(name_or_path)

    haskey(SETUPS, name_or_path) || throw(
        ArgumentError(
            "Unknown setup \"$name_or_path\". Available setups: $(join(setup_names(), ", ")). " *
            "Pass a path ending in .jl to load a config file instead.",
        ),
    )

    return SETUPS[name_or_path]()
end

"""
    load_config_file(path)

Load an out-of-tree config file and return the `FjordConfig` it evaluates to.

Evaluated in `Main` rather than in this module: the file's own `using FjordSim` and any helper
bindings belong in the caller's namespace, and a bare `include` here would resolve `path`
relative to `src/Setups/`.
"""
function load_config_file(path)
    filepath = abspath(expanduser(path))
    isfile(filepath) || throw(ArgumentError("Config file $filepath does not exist"))

    config = Base.include(Main, filepath)
    config isa FjordConfig || throw(
        ArgumentError(
            "$filepath evaluated to a $(typeof(config)), not a FjordConfig. " *
            "A config file's last expression must be the FjordConfig it describes.",
        ),
    )

    return config
end

end  # module Setups
