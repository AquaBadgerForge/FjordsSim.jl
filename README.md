# FjordSim.jl

A framework for ocean simulations built on top of [Oceananigans](https://github.com/CliMA/Oceananigans.jl) and [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl).

FjordSim's main contribution is a streamlined way to set up regional simulations.
A simulation is assembled from the following components:

1. **Grid**
   Domain bounds, horizontal size, and vertical faces.

2. **Bathymetry file**
   Contains the domain coordinates (longitude, latitude, depth) along with the bathymetric data.
   It can be generated from the public [Geonorge](https://www.geonorge.no/) Sjøkart Dybdedata dataset, see below.

3. **Forcing file**
   Includes information about sinks and sources (e.g., rivers), boundary conditions, and custom forcings.
   It can also be used to load initial conditions.

4. **Atmospheric data**
   Supports JRA55 from [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl) and [NORA3](https://thredds.met.no/thredds/projects/nora3.html).

## Installation

There are several options:

- Install from github:
1. Clone the git repository `git clone https://github.com/NIVANorge/FjordSim.jl.git`.
2. Move to the downloaded folder `cd FjordSim.jl`.
3. Run Julia REPL and activate the FjordSim environment `julia --project`.
4. Enter the Pkg REPL by pressing `]` from Julia REPL.
5. Type `instantiate` to 'resolve' a `Manifest.toml` from a `Project.toml` to install and precompile dependency packages.

- Add the latest FjordSim to your Julia project: `add https://github.com/NIVANorge/FjordSim.jl.git`.

- Add from the Julia registry: `add FjordSim`.

## Setups (work in progress)

Each fjord is a function in `src/Setups/` (`oslofjorden.jl`, `drammensfjorden.jl`) returning a
`FjordConfig`, which groups a grid, a bathymetry, a forcing and an atmosphere configuration:

```julia
function drammensfjorden()
    data_root = joinpath(homedir(), "FjordSim_data", "drammensfjorden")

    return FjordConfig(
        grid_config       = EvenGrid(size = (150, 200, 11), longitude = (10.20, 10.45), ...),
        bathymetry_config = DybdedataConfig(data_root = data_root, output_file = "bathymetry.nc", ...),
        forcing_config    = NorKystConfig(data_root = data_root, years = [2020]),
    )
end
```

`setup_names()` lists the registered setups and `fjord_config(name)` returns one by name; that is
what `--config` resolves.

By default, input data goes to `$HOME/FjordSim_data/<fjord>/` and results to
`$HOME/FjordSim_results/<fjord>/`. The config fields naming individual files are relative to the
setup's `data_root`; setting one to an absolute path relocates just that file, which is how a
single copy of the Geonorge database is shared between fjords.

To add a fjord to the package, copy an existing setup file into `src/Setups/`, adjust it, and add it
to the `SETUPS` registry in `src/Setups/Setups.jl`.

### A fjord outside the package

A fjord does not have to live in FjordSim. Put the `FjordConfig` in a standalone `.jl` file whose
last expression is the config, and pass its path to `--config`:

```julia
# ~/fjords/hardangerfjorden.jl
using FjordSim

data_root = joinpath(homedir(), "FjordSim_data", "hardangerfjorden")

FjordConfig(
    grid_config = EvenGrid(
        size      = (200, 300, 12),
        halo      = (7, 7, 7),
        longitude = (5.3, 6.5),
        latitude  = (59.9, 60.5),
        # Nz + 1 faces, bottom to top, ending at 0.0.
        z_faces   = [
            -800.0, -600.0, -400.0, -250.0, -150.0, -100.0,
            -50.0, -25.0, -15.0, -10.0, -5.0, -2.0, 0.0,
        ],
    ),
    bathymetry_config = DybdedataConfig(
        data_root   = data_root,
        output_file = "bathymetry.nc",
        plot_file   = "bathymetry.png",
        # Reuse the 2.3 GB Geonorge FileGDB already downloaded for another fjord instead of
        # fetching a second copy.
        geodatabase_file = joinpath(
            homedir(), "FjordSim_data", "oslofjorden",
            "Basisdata_0000_Norge_25833_Dybdedata_FGDB.gdb",
        ),
    ),
    forcing_config = NorKystConfig(
        data_root            = data_root,
        output_directory     = "norkyst",
        output_file          = "forcing.nc",
        plot_file            = "forcing.png",
        relaxation_edge      = :west,
        relaxation_cells     = 10,
        relaxation_timescale = 86400.0,
        architecture         = :auto,
        parameters           = ["temperature", "salinity", "u_eastward", "v_northward"],
        years                = [2020],
    ),
)
```

```bash
julia --project -m FjordSim prepare_bathymetry --config ~/fjords/hardangerfjorden.jl
```

or, equivalently, from the REPL:

```julia
using FjordSim
config = fjord_config(expanduser("~/fjords/hardangerfjorden.jl"))
prepare_bathymetry(config)
```

Three things to know:

- **The `.jl` suffix decides how `--config` reads its argument.** Anything ending in `.jl` is loaded
  as a file; anything else is looked up in the registry, so a misspelled setup name reports the
  available setups rather than a missing file.
- **The file's last expression must be the `FjordConfig`** — no `return`, no trailing assignment. A
  file evaluating to anything else is rejected up front instead of failing inside a pipeline later.
- **The file is evaluated in `Main`**, so it needs its own `using FjordSim`. `~` and relative paths
  are expanded.

Registering a fjord only buys the shorter `--config <name>`, its appearance in `--help`, and
coverage by the tests that loop over `setup_names()`. For a one-off fjord, the standalone file is
the better choice.

To add a new kind of grid, bathymetry source, or forcing dataset, define a struct subtyping
`AbstractGridConfig`, `AbstractBathymetryConfig`, or `AbstractForcingConfig` and overload that
supertype's hooks — `FjordConfig`, the generic pipelines and the command line stay unchanged. Each
pipeline is one generic function plus a handful of dispatch points:

| Pipeline | Generic entry point | Hooks a new source implements |
|---|---|---|
| Grid | — | `LatitudeLongitudeGrid(architecture, config)` |
| Bathymetry | `prepare_bathymetry(target_grid, config)` | `bathymetry_dataset(target_grid, config)`; optionally `regrid_options(config)` |
| Forcing | `prepare_forcing(target_grid, config)`, `download_forcing(config)` | `forcing_time_steps`, `forcing_source_grid`, `forcing_variable_names`; `download_forcing(target_grid, config)` if it downloads |
| Simulation | `build_simulation(config)`, `run_simulation(config)` | none — `AbstractSimulationConfig` is fields only |

`AbstractAtmosphereConfig` additionally has two read-side hooks, `prescribed_atmosphere(config,
architecture)` and `prescribed_radiation(config, architecture)`, which is how the simulation
reads a prepared atmosphere without naming a dataset.

Path resolution (`bathymetry_path`, `forcing_path`, `forcing_directory`, `results_path`,
`plot_path`) and the diagnostic plots (`plot_bathymetry`, `plot_forcing`) are defined on the
supertypes, so a new source inherits them. `src/Bathymetry/geonorge.jl` and `src/Forcing/norkyst.jl` are the built-in
adapters and the templates to copy; each supertype's docstring in `src/Configs.jl` lists the
fields and hooks it expects.

## Prepare input data

Every step takes the setup to prepare via `--config`, which is required and is the only option
— everything else is stated in the setup:

```bash
# Bathymetry: downloads the Geonorge Sjøkart FileGDB on first use (~2.3 GB), regrids it onto the
# configured grid, and writes the bathymetry NetCDF plus a diagnostic plot
julia --project -m FjordSim prepare_bathymetry --config oslofjorden

# Forcing: downloads and subsets the configured dataset (NorKyst-800m) over the setup's region
# and years
julia --project -m FjordSim download_forcing --config oslofjorden

# Forcing: regrids the download onto the simulation grid, writing the forcing NetCDF and a plot
julia --project -m FjordSim prepare_forcing --config oslofjorden

# Rivers: writes river relaxation into a copy of the prepared forcing
julia --project -m FjordSim add_rivers --config oslofjorden

# Atmosphere: downloads and subsets NORA3, then regrids it onto a regular lon/lat grid.
# The download is by far the slowest step — a year is close to 10000 OPeNDAP reads — and skips
# any month already on disk, so an interrupted run resumes.
julia --project -m FjordSim download_atmosphere --config oslofjorden
julia --project -m FjordSim prepare_atmosphere --config oslofjorden

# Simulation: builds and runs the coupled model, writing NetCDF snapshots to the results directory
julia --project -m FjordSim run_simulation --config oslofjorden

# All of the above, with the setups they accept
julia --project -m FjordSim --help
```

`-m` is Julia 1.12's package entry point. Each subcommand is named after the function it calls, so
the same steps run from the REPL:

```julia
using FjordSim
config = oslofjorden()
prepare_bathymetry(config)
download_forcing(config)
prepare_forcing(config)
```

A step the setup does not configure — `add_rivers` on a setup with no rivers, the atmosphere steps
on a setup with no atmosphere, `run_simulation` on a setup with no simulation config — does
nothing rather than failing.

Every step that runs writes a full transcript to `fjordsim.log` in the setup's results directory,
truncated each run, while still printing live to the terminal. If a step fails, that log is where to
read the error: a stacktrace through the simulation config spells out every type parameter, so the
error message itself scrolls out of the terminal long before the trace ends.

```bash
julia --project -m FjordSim run_simulation --config oslofjorden
head -60 ~/FjordSim_results/oslofjorden/fjordsim.log   # the error, before the stacktrace buries it
```

The log lands beside the output it describes. `results_root` is a field of the simulation config, so
a setup that names none — `drammensfjorden` — has no results directory and writes `fjordsim.log`
into the working directory instead. Either way, the path is printed as the run starts.

The forcing config's `architecture` field decides where the regridding interpolation runs:
`:auto` (the default) uses the GPU when one is usable and the CPU otherwise, so the same config
works on a compute node and a laptop; `:cpu` and `:gpu` pin it, and `:gpu` errors rather than
silently falling back to a ~12x slower CPU run.

## Run a simulation

Once every step the setup configures has run, `run_simulation` builds the coupled model from the
setup's own prepared files and runs it:

```bash
julia --project -m FjordSim run_simulation --config oslofjorden
```

or, equivalently, `julia --project examples/oslofjord.jl`. Snapshots go to the simulation config's
`results_root`, `$HOME/FjordSim_results/oslofjorden/` for the built-in Oslofjord setup — separate
from the input data root.

Everything the run needs comes from the setup: the grid from the processed bathymetry, the forcing
from the rivers-augmented file when the setup names rivers and the plain prepared file otherwise,
and the atmosphere from the file `prepare_atmosphere` wrote. A missing prerequisite is reported as
the command that produces it. The simulation config's `architecture` field selects the device, the
same `:auto`/`:cpu`/`:gpu` as the forcing config's.

To inspect or step through the assembly instead of running it, build the simulation without
starting it:

```julia
using FjordSim
simulation = build_simulation(oslofjorden())
run!(simulation)
```

What the simulation is made of — the buoyancy, the closure, the advection schemes, the tracers,
the initial conditions, the coriolis, the sea ice, the run length, the output interval and the
time-step wizard settings — is the `SimulationConfig` in `src/Setups/oslofjorden.jl`.
`SimulationConfig` has no defaults, so that block is the whole story: every knob is stated there,
and omitting one is an `UndefKeywordError` naming it rather than a silently inherited value.

The grid, forcing and atmospheric forcing for the Oslofjord are also available
[here](https://www.dropbox.com/scl/fo/gc3yc155b5eohi7998wgh/AGN2Yt3HyQ0LlZGImpcca6o?rlkey=x6okc3uxe2avud6sbxgd00l14&st=093llyqp&dl=0)
if you would rather not run the preparation steps yourself.

Further preparation scripts for the Oslofjord are available in the following repository:
[https://github.com/NIVANorge/oslofjord-sim](https://github.com/NIVANorge/oslofjord-sim)

![example_result](./artifacts/phytoplankton_multi.png)
