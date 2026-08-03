# FjordSim.jl

A framework for regional ocean simulations built on top of
[Oceananigans](https://github.com/CliMA/Oceananigans.jl) and
[NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl).

## The idea

Regional ocean models are usually set up one of two ways: through config files — text or YAML,
the way [ROMS](https://www.myroms.org/) does it — or through a programmable script, the way
Oceananigans and NumericalEarth themselves work. A config file lets you see every choice a run
makes in one place, but only the choices its authors anticipated: adding genuinely new behavior
means leaving the config format and writing code. A script can do anything, but the setup is now
spread across however much of the script builds it, and nothing stops two runs of "the same
model" from silently diverging.

FjordSim does both at once. A setup is a `FjordConfig` — a single Julia value a user can read
top to bottom — but its fields are not limited to numbers and strings. A field can hold a live,
programmable object: a turbulence closure, an advection scheme, a whole dataset adapter with its
own download and regrid logic. Changing the *logic* of a run — swapping in a different closure,
adding a new atmosphere source — is done by writing a new object and putting it in the config,
not by writing a parallel script; and because it is still a config, the complete list of what a
run does is still visible in one file rather than scattered across a codebase. This is a
direction more than a finished state: today the programmable fields are grid, bathymetry,
forcing, river, atmosphere and simulation configs, but more of the pipeline is meant to become
config-visible this way over time, and configs are expected to nest more deeply as they do (a
forcing config already holds a river config, for instance).

## Main data components

A simulation is assembled from three data components, each regridded onto — or around — the
same simulation grid:

1. **Bathymetry** — the domain's depth and land/sea mask. `prepare_bathymetry` regrids a
   bathymetry source onto the simulation grid (built from the grid config's bounds and
   resolution) and writes the result as a NetCDF file, which is what actually defines the model's
   `ImmersedBoundaryGrid`. The built-in source is the public
   [Geonorge](https://www.geonorge.no/) Sjøkart Dybdedata dataset.

2. **Forcing** — boundary and initial conditions from a larger-scale ocean model.
   `prepare_forcing` regrids a regional reanalysis (built in: NorKyst-800m) onto the simulation
   grid, producing relaxation values and lambdas along one open edge, plus optional river
   relaxation (`add_rivers`) written on top. The same prepared file can also seed a run's initial
   conditions, so the ocean starts from a realistic state instead of a uniform water column.

3. **Atmospheric data** — the surface forcing a run needs but the regional ocean model does not
   carry. `prepare_atmosphere` regrids a reanalysis product (built in: MET Norway's
   [NORA3](https://thredds.met.no/thredds/projects/nora3.html)) onto its own regular
   longitude/latitude grid — independent of the ocean grid, since NumericalEarth interpolates
   between them — producing the eight variables a `PrescribedAtmosphere`/`PrescribedRadiation`
   pair needs.

Bathymetry and forcing are downloaded and regridded onto the *same* grid, so `prepare_forcing`
needs the processed bathymetry to build the model's land mask. The atmosphere is independent of
both. A setup can configure any subset of these — omitting an atmosphere or a simulation config
just makes the corresponding steps no-ops — but the grid, bathymetry and forcing are always
required.

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

## Quick start: Drammensfjorden

`drammensfjorden()` in `src/Setups/drammensfjorden.jl` is the shortest complete example — the
same physics as the built-in Oslofjord setup, but no rivers and a 30-day run instead of a full
year, so it is a fast way to try every pipeline step once:

```julia
function drammensfjorden()
    data_root = joinpath(homedir(), "FjordSim_data", "drammensfjorden")
    FT = Oceananigans.defaults.FloatType

    return FjordConfig(
        grid_config       = EvenGrid(size = (150, 200, 11), longitude = (10.20, 10.45), ...),
        bathymetry_config = DybdedataConfig(data_root = data_root, output_file = "bathymetry.nc", ...),
        forcing_config    = NorKystConfig(data_root = data_root, years = [2020], ...),
        atmosphere_config = NORA3Config(data_root = data_root, years = [2020], ...),
        simulation_config = SimulationConfig(
            buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
            closure            = (CATKEVerticalDiffusivity(minimum_tke = 7e-6), ...),
            initial_conditions = FromForcing(),   # the NorKyst state at start_date
            start_date         = DateTime(2020, 1, 1),
            stop_time          = 30days,          # a short window, not a full year
            ...
        ),
    )
end
```

Preprocessing the data and running it is the same six commands as any setup with rivers left out
— this one names none:

```bash
julia --project -m FjordSim prepare_bathymetry --config drammensfjorden
julia --project -m FjordSim download_forcing --config drammensfjorden
julia --project -m FjordSim prepare_forcing --config drammensfjorden
julia --project -m FjordSim download_atmosphere --config drammensfjorden
julia --project -m FjordSim prepare_atmosphere --config drammensfjorden
julia --project -m FjordSim run_simulation --config drammensfjorden
```

or, equivalently, from the REPL:

```julia
using FjordSim
config = drammensfjorden()
prepare_bathymetry(config)
download_forcing(config)
prepare_forcing(config)
download_atmosphere(config)
prepare_atmosphere(config)
run_simulation(config)
```

Snapshots land in `~/FjordSim_results/drammensfjorden/`. To inspect or step through the assembled
model instead of running it outright, use `build_simulation` in place of `run_simulation`:

```julia
using FjordSim
simulation = build_simulation(drammensfjorden())
run!(simulation)
```

The full reference for every subcommand — including the features this quick setup does not use,
like rivers, looping and checkpointing — is under "Prepare input data" and "Run a simulation"
below, using the built-in Oslofjord setup.

## Setups

Each fjord is a zero-argument function in `src/Setups/` (`oslofjorden.jl`, `drammensfjorden.jl`)
returning a `FjordConfig`, which groups a grid, a bathymetry, a forcing, and optionally an
atmosphere and a simulation configuration. `setup_names()` lists the registered setups and
`fjord_config(name)` returns one by name; that is what `--config` resolves.

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

### Adding a new source

To add a new kind of grid, bathymetry source, forcing dataset, river dataset or atmosphere
dataset, define a struct subtyping the matching supertype and overload its hooks —
`FjordConfig`, the generic pipelines and the command line stay unchanged. Each pipeline is one
generic function plus a handful of dispatch points; each supertype's docstring in
`src/Configs.jl` lists the full field and hook contract, and the adapter files below are the
templates to copy.

| Config | Generic entry point | Hooks a new source implements | Template |
|---|---|---|---|
| Grid (`AbstractGridConfig`) | — | `LatitudeLongitudeGrid(architecture, config)` | `src/Grids.jl` |
| Bathymetry (`AbstractBathymetryConfig`) | `prepare_bathymetry(target_grid, config)` | `bathymetry_dataset(target_grid, config)`; optionally `regrid_options`, `smoothing_options` | `src/Bathymetry/geonorge.jl` |
| Forcing (`AbstractForcingConfig`) | `prepare_forcing(target_grid, config)`, `download_forcing(config)` | `forcing_time_steps`, `forcing_source_grid`, `forcing_variable_names`; `download_forcing(target_grid, config)` if it downloads | `src/Forcing/norkyst.jl` |
| Rivers (`AbstractRiverConfig`) | `add_rivers(target_grid, config)` | `river_locations`, `river_series`; `download_rivers` if it downloads | `src/Forcing/of800_rivers.jl` |
| Atmosphere (`AbstractAtmosphereConfig`) | `prepare_atmosphere(target_grid, config)`, `download_atmosphere(config)` | `atmosphere_time_steps`, `atmosphere_source_grid`, `atmosphere_variable_names`; `download_atmosphere(target_grid, config)` if it downloads; `prescribed_atmosphere`, `prescribed_radiation` if the setup is simulated | `src/Atmospheres/nora3_source.jl` |
| Simulation (`AbstractSimulationConfig`) | `build_simulation(config)`, `run_simulation(config)` | none — fields only | `src/Setups/oslofjorden.jl` |

A river config is not a `FjordConfig` field on its own — it goes in the forcing config's `rivers`
field, `nothing` for a setup with no rivers. Path resolution (`bathymetry_path`, `forcing_path`,
`forcing_directory`, `river_forcing_path`, `atmosphere_path`, `atmosphere_directory`,
`results_path`, `plot_path`) and the diagnostic plots (`plot_bathymetry`, `plot_forcing`,
`plot_atmosphere`) are defined on the supertypes, so a new source inherits them for free.

## Prepare input data

Every step takes the setup to prepare via `--config`, which is required and is the only option
— everything else is stated in the setup. The Oslofjord setup is used below since it is the one
that configures every step, including rivers:

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

`run_simulation` writes a full transcript to `fjordsim_<launch time>.log` in the setup's results
directory, while still printing live to the terminal; the other steps print to the terminal only. If
the run fails, that log is where to read the error: a stacktrace through the coupled model scrolls the
error message itself out of the terminal long before the trace ends. Type parameters in the trace are
abbreviated to `{…}`, as they are in the REPL, so a frame stays one readable line.

```bash
julia --project -m FjordSim run_simulation --config oslofjorden
# the path is printed as the run starts; the error is at the top of it
head -60 ~/FjordSim_results/oslofjorden/fjordsim_20260803T141530.log
```

The log lands beside the output it describes, and carries the same launch-time tag as the snapshots,
so runs do not overwrite each other's transcripts.

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
`initial_conditions` can be a literal `NamedTuple`, `FromForcing()` (the prepared forcing's own
state at `start_date`), or `FromResults(path)` (a previous run's snapshot); `drammensfjorden()`
above uses `FromForcing()`.

A run can repeat its window several times with the ocean state carried over (`loops`, a spin-up
for basins that do not equilibrate in one forcing year), and can write periodic checkpoints
(`checkpoint_interval`) that a later run resumes from with `pickup = true`. The quick-start setup
above sets `loops = 1` and `checkpoint_interval = 0.0` (no checkpointing); `oslofjorden()` shows
checkpointing turned on (`checkpoint_interval = 30days`), and both fields are documented on
`SimulationConfig` in `src/Simulations.jl`.

The grid, forcing and atmospheric forcing for the Oslofjord are also available
[here](https://www.dropbox.com/scl/fo/gc3yc155b5eohi7998wgh/AGN2Yt3HyQ0LlZGImpcca6o?rlkey=x6okc3uxe2avud6sbxgd00l14&st=093llyqp&dl=0)
if you would rather not run the preparation steps yourself.

Further preparation scripts for the Oslofjord are available in the following repository:
[https://github.com/NIVANorge/oslofjord-sim](https://github.com/NIVANorge/oslofjord-sim)

![example_result](./artifacts/phytoplankton_multi.png)
