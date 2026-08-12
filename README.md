# FjordSim.jl

A framework for regional ocean simulations built on top of
[Oceananigans](https://github.com/CliMA/Oceananigans.jl) and
[NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl).

A simulation is completely defined in a config file, for example ['oslofjordenl.jl'](src/Setups/oslofjorden.jl).
It contains both the simulation parameters and all the logic of a simulation.
It is done using Julia main feature/paradigm - multiple dispatch.
Providing new structs and overloading existing/new methods on them, a user can change anything.
The same time, since the entire simulation is a config type file, all the parts are in one place, compact, and readable.
For example, [`examples/oslofjorden.jl`](examples/oslofjorden.jl) is a modified version of the original `src/Setups/oslofjorden.jl` but with implicit sea surface and open velocity boundary conditions (in the original setup has free explicit sea surface and velocity boundary conditions are closed waiting for open bc fro explicit free surface to appear in oceananigans).

## Installation

```bash
git clone https://github.com/NIVANorge/FjordSim.jl.git
cd FjordSim.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

`instantiate` resolves a `Manifest.toml` from `Project.toml`, then installs and precompiles the dependencies.
To track the repository from a project of your own rather than working inside it,
`add https://github.com/NIVANorge/FjordSim.jl.git`.

Julia **1.12 or newer** is required: the `-m FjordSim` command form used throughout this README is
Julia 1.12's package entry point.

A GPU is not strictly required — `architecture = :auto` falls back to the CPU — but running a
simulation without one is impractically slow. Preparing the input data is fine on a laptop.

## Quick start: Drammensfjorden

`drammensfjorden()` in `src/Setups/drammensfjorden.jl` is the smallest complete setup — the same
physics as the built-in Oslofjord one on a smaller domain.
Though data downloading still takes a while.

Preparing the data and running it is seven commands, in this order:

```bash
julia --project -m FjordSim prepare_bathymetry  --config drammensfjorden
julia --project -m FjordSim download_forcing    --config drammensfjorden
julia --project -m FjordSim prepare_forcing     --config drammensfjorden
julia --project -m FjordSim add_rivers          --config drammensfjorden
julia --project -m FjordSim download_atmosphere --config drammensfjorden
julia --project -m FjordSim prepare_atmosphere  --config drammensfjorden
julia --project -m FjordSim run_simulation      --config drammensfjorden
```

Each subcommand is named after the function it calls, so the same steps run from the REPL:

```julia
using FjordSim
config = drammensfjorden()
prepare_bathymetry(config)
download_forcing(config)
prepare_forcing(config)
add_rivers(config)
download_atmosphere(config)
prepare_atmosphere(config)
run_simulation(config)
```

Snapshots land in `~/FjordSim_results/drammensfjorden/`. 

## What a run is made of

A simulation is assembled from three data components, each regridded onto — or around — the same
simulation grid:

1. **Bathymetry** — the domain's depth and land/sea mask. `prepare_bathymetry` regrids a bathymetry
   source onto the simulation grid (built from the grid config's bounds and resolution) and writes the
   result as a NetCDF file. That file is what actually defines the model's `ImmersedBoundaryGrid`. The
   built-in source is the public [Geonorge](https://www.geonorge.no/) Sjøkart Dybdedata dataset.

2. **Forcing** — boundary and initial conditions from a larger-scale ocean model. `prepare_forcing`
   regrids a regional reanalysis (built in: NorKyst-800m) onto the simulation grid, producing
   relaxation values and lambdas along one open edge. `add_rivers` then writes river relaxation into a
   copy of that file. The same prepared file can seed a run's initial conditions, so the ocean starts
   from a realistic state rather than a uniform water column.

3. **Atmospheric data** — the surface forcing a run needs but a regional ocean model does not carry.
   `prepare_atmosphere` regrids a reanalysis product (built in: MET Norway's
   [NORA3](https://thredds.met.no/thredds/projects/nora3.html)) onto its own regular longitude/latitude
   grid, independent of the ocean grid since NumericalEarth interpolates between them, producing the
   eight variables a `PrescribedAtmosphere`/`PrescribedRadiation` pair needs.

Bathymetry and forcing go onto the *same* grid, which is why `prepare_forcing` needs the processed
bathymetry: that is where its land mask comes from. The atmosphere is independent of both.

A setup can configure any subset. Omitting the rivers, the atmosphere or the simulation config makes
the corresponding steps no-ops rather than errors. The grid, bathymetry and forcing are always
required.

### A fjord outside the package

A fjord does not have to live in FjordSim. Put the `FjordConfig` in a standalone `.jl` file whose last
expression is the config, and pass its path to `--config`:

```bash
julia --project -m FjordSim run_simulation --config examples/oslofjorden.jl
```
