# FjordSim.jl

A framework for regional ocean simulations built on top of
[Oceananigans](https://github.com/CliMA/Oceananigans.jl) and
[NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl).

A simulation is completely defined in a config file, for example
[`src/Setups/oslofjorden.jl`](src/Setups/oslofjorden.jl).
It contains the simulation parameters and the programming objects that define the logic of a simulation.
It is done using Julia's main feature/paradigm — multiple dispatch.
Providing new structs and overloading existing/new methods on them, a user can change anything.
At the same time, since the entire simulation is a config type file, all the parts are in one place,
compact, and readable.

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

## Quick start: Oslofjorden

`oslofjorden()` in `src/Setups/oslofjorden.jl` is the reference setup: Geonorge bathymetry,
NorKyst-800m forcing with NVE rivers, hourly NorKyst open-boundary data along the southern edge,
and NORA3 atmosphere. 
Please note: It may take several hours to download the data.

`add_rivers` also needs a free NVE HydAPI key — sign up at https://hydapi.nve.no/Users and set it
as `NVE_API_KEY` in your environment before running that step.

Preparing the data and running it is nine commands, in this order:

```bash
julia --project -m FjordSim prepare_bathymetry   --config oslofjorden
julia --project -m FjordSim download_forcing     --config oslofjorden
julia --project -m FjordSim prepare_forcing      --config oslofjorden
julia --project -m FjordSim add_rivers           --config oslofjorden
julia --project -m FjordSim download_boundaries  --config oslofjorden
julia --project -m FjordSim prepare_boundaries   --config oslofjorden
julia --project -m FjordSim download_atmosphere  --config oslofjorden
julia --project -m FjordSim prepare_atmosphere   --config oslofjorden
julia --project -m FjordSim run_simulation       --config oslofjorden
```

Each subcommand is named after the function it calls, so the same steps run from the REPL:

```julia
using FjordSim
config = oslofjorden()
prepare_bathymetry(config)
download_forcing(config)
prepare_forcing(config)
add_rivers(config)
download_boundaries(config)
prepare_boundaries(config)
download_atmosphere(config)
prepare_atmosphere(config)
run_simulation(config)
```

`build_simulation(config)` returns the assembled simulation without starting it, which is the entry
point for the REPL and the debugger. Snapshots and a per-launch log land in
`~/FjordSim_results/oslofjorden/`; `julia --project -m FjordSim --help` lists every subcommand and
setup.

## What a run is made of

A simulation is assembled from four data components, each regridded onto — or around — the same
simulation grid:

1. **Bathymetry** — the domain's depth and land/sea mask. `prepare_bathymetry` regrids a bathymetry
   source onto the simulation grid (built from the grid config's bounds and resolution) and writes the
   result as a NetCDF file. That file is what actually defines the model's `ImmersedBoundaryGrid`. The
   built-in source is the public [Geonorge](https://www.geonorge.no/) Sjøkart Dybdedata dataset.

2. **Forcing** — the interior state from a larger-scale ocean model. `prepare_forcing` regrids a
   regional reanalysis (built in: NorKyst-800m daily means) onto the simulation grid; `add_rivers`
   then writes river relaxation into a copy of that file. The prepared file is read for river forcing
   and, optionally, for a run's initial conditions, so the ocean starts from a realistic state rather
   than a uniform water column.

3. **Open-boundary data** — the exterior state along the open lateral edges, which is what makes them
   open. `prepare_boundaries` regrids the *hourly* NorKyst collection onto the boundary rows,
   including surface elevation and depth-averaged velocities, so the tide survives. It is a separate
   file from the forcing precisely because of that cadence: only the boundary rows need it. The
   boundary condition itself is Flather on the transport plus Orlanski radiation with inflow nudging
   on everything else.
   Any subset of the four edges can be open — a fjord names one, a region in the open ocean names all
   four, and all of them share one prepared file.

4. **Atmospheric data** — the surface forcing a run needs but a regional ocean model does not carry.
   `prepare_atmosphere` regrids a reanalysis product (built in: MET Norway's
   [NORA3](https://thredds.met.no/thredds/projects/nora3.html)) onto its own regular longitude/latitude
   grid, independent of the ocean grid since NumericalEarth interpolates between them, producing the
   eight variables a `PrescribedAtmosphere`/`PrescribedRadiation` pair needs.

Bathymetry, forcing and open-boundary data go onto the *same* grid, which is why `prepare_forcing`
and `prepare_boundaries` need the processed bathymetry: that is where their land mask comes from. The
atmosphere is independent.

### Skipping steps

Only the grid and the bathymetry are always required. 
Every other component is a config field defaulting to `nothing`, and a step whose config is absent is a no-op rather than an error — so a setup names what it wants and skips the rest.

With `download_forcing` and `prepare_forcing` skipped rivers can be kept without the NorKyst download and regrid by giving the river config `standalone = true`, which makes `add_rivers` write a river-only forcing file from scratch. 
What does need a prepared forcing file is `FromForcing` initial conditions;
use constants, functions, or `FromResults` from an earlier run instead.

### A fjord outside the package

A fjord does not have to live in FjordSim. Put the `FjordConfig` in a standalone `.jl` file whose last
expression is the config, and pass its path to `--config`:

```bash
julia --project -m FjordSim run_simulation --config examples/oslofjorden_npzd.jl
```

## Some `run_simulation` details

`run_simulation` has no command-line options: everything is stated in the setup's
`SimulationConfig`. These fields control the run itself.

| Field | Meaning |
|---|---|
| `start_date` | the calendar instant model time zero stands for |
| `stop_time` | how long one run window is |
| `loops` | how many times to repeat that window, carrying the ocean state over |
| `initial_conditions` | where the first state comes from |
| `pickup` | resume from the newest checkpoint instead of starting fresh |


### Checkpointing

```julia
CheckpointWriter(interval = 30days, cleanup = true)
```

It writes the coupled model's full prognostic state as JLD2 — a few hundred MB each on a real grid.
The filename carries the iteration, so every fire writes a new file and `cleanup = true` is what
keeps only the newest; `false` accumulates one per interval. Only one checkpointing writer is
allowed.

### Continuing a run

Set `pickup = true` and the newest checkpoint is restored — clock, time
step, tracers, velocities, free surface, CATKE's diffusivities and the Adams-Bashforth tendencies — so
this is an *exact* continuation, and it supersedes `initial_conditions` (they are set at build time
and then overwritten). `pickup` without a `CheckpointWriter` is a configuration error, not a run that
fails later; `pickup` with no checkpoint on disk warns and starts from the beginning.

The resumed run gets its own run tag, so its snapshots start at the resume point and the records
before it stay in the previous launch's file.

### Looping

`loops > 1` runs the same window repeatedly, saving the ocean state. 
Each iteration resets all clocks to zero, saves the time step reached by the master, and writes its own file (`snapshots_ocean_<tag>_loop02.nc`) so that loops can be compared without overwriting each other.

### Starting from forcing or from a previous run

`initial_conditions` takes one of three shapes:

```julia
initial_conditions = (T = 5.0, S = 33.0)                       # constants, functions or fields
initial_conditions = FromForcing()                             # prepared forcing at `start_date`
initial_conditions = FromForcing(DateTime(2020, 3, 1))         # ...or at a named date
initial_conditions = FromResults("snapshots_ocean_20260814T143012.nc")   # a previous run's last record
initial_conditions = FromResults("snapshots_ocean_20260814T143012.nc", DateTime(2020, 6, 1))
```

All are a **lossy warm start**: the tracers from the config plus `u` and `v` are set, while the free surface, the turbulence state and the timestepper tendencies keep their defaults, so the first hours are a barotropic adjustment. 
That is the difference from `pickup`, which is exact — the two are complementary, not alternatives.
