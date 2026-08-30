# CLAUDE.md

> This file holds what's needed on every task: commands, hard rules, and gotchas. Design rationale,
> the full per-module architecture writeup, the "adding a new source" hook contracts, and per-fjord
> tuning history live in `docs/` and are meant to be read on demand — see the pointers below.

## Commands

```bash
# Run all tests
julia --project test/runtests.jl

# Run the out-of-tree Oslofjord example config (requires GPU + data files)
julia --project -m FjordSim run_simulation --config examples/oslofjorden.jl

# Prepare bathymetry for a configured fjord (downloads the Geonorge GDB on first use)
julia --project -m FjordSim prepare_bathymetry --config drammensfjorden

# Download and subset the configured forcing dataset for a fjord
julia --project -m FjordSim download_forcing --config oslofjorden

# Regrid the downloaded forcing onto the fjord's grid (needs the bathymetry and download first)
julia --project -m FjordSim prepare_forcing --config oslofjorden

# Write river relaxation on top of the prepared forcing (needs prepare_forcing first, unless the
# river config is `standalone`, which writes a forcing file carrying only rivers). `oslofjorden`
# discovers its river mouths from NVE's open map services and reads its gauges from NVE's HydAPI,
# so this step needs `NVE_API_KEY` — free at https://hydapi.nve.no/Users.
julia --project -m FjordSim add_rivers --config oslofjorden

# Download the hourly exterior state along the open lateral boundary (a thin band, not the box)
julia --project -m FjordSim download_boundaries --config oslofjorden

# Regrid it onto the open edge of the fjord's grid (needs the bathymetry and download first)
julia --project -m FjordSim prepare_boundaries --config oslofjorden

# Download and subset the configured atmosphere dataset (slow: ~10000 OPeNDAP reads per year)
julia --project -m FjordSim download_atmosphere --config oslofjorden

# Regrid the downloaded atmosphere onto a regular lon/lat grid (needs the download first)
julia --project -m FjordSim prepare_atmosphere --config oslofjorden

# Build and run the coupled simulation (needs every prepare step the setup configures, plus a GPU)
julia --project -m FjordSim run_simulation --config oslofjorden

# The subcommands and the setups they accept
julia --project -m FjordSim --help
```

`--config` is the only option, and it is required — there is no default setup. It takes a
registered setup name (`FjordSim.Setups.SETUPS`) or a path to an out-of-tree `.jl` config file
whose last expression is a `FjordConfig`. Every other knob, including which device the forcing
interpolation runs on, is a config field.

### Changing `start_date` or `stop_time`

Both prepare steps pad their time axes to span the run window the simulation config names
(`coverage_window`), so moving either end invalidates the prepared files. Re-run, in this order:

```bash
julia --project -m FjordSim prepare_forcing     --config oslofjorden   # rewrites forcing.nc
julia --project -m FjordSim add_rivers          --config oslofjorden   # re-copies forcing_rivers_nve.nc
julia --project -m FjordSim prepare_boundaries  --config oslofjorden   # rewrites boundaries.nc
julia --project -m FjordSim prepare_atmosphere  --config oslofjorden   # rewrites atmosphere.nc
```

`add_rivers` is not optional here: it `cp`s the forcing file and patches the copy, so the river file
— which is what `simulation_forcing_path` gives the simulation, `forcing_rivers_nve.nc` on
`oslofjorden` and `forcing_rivers.nc` on `drammensfjorden` — still carries
the *old* axis until it is re-run. A `standalone` river config takes the window even more directly,
building its whole time axis from it. No download step is affected; all three prepares are pure
regrids.
`loops` is not part of the window, so changing it re-runs nothing.

### Changing `z_faces` or the grid `size`

Also a data change, and a larger one. `z_faces` is written into `bathymetry.nc`, and
`simulation_grid` reads the *file* rather than the config — so a config edit alone changes nothing
about the run, and every 3D prepared file was regridded onto the old geometry. Re-run, in this order:

```bash
julia --project -m FjordSim prepare_bathymetry  --config oslofjorden   # rewrites bathymetry.nc
julia --project -m FjordSim prepare_forcing     --config oslofjorden   # rewrites forcing.nc
julia --project -m FjordSim add_rivers          --config oslofjorden   # re-copies forcing_rivers_nve.nc
julia --project -m FjordSim prepare_boundaries  --config oslofjorden   # rewrites boundaries.nc
```

No download step is affected — all four are pure regrids of data already on disk. `prepare_atmosphere`
is *not* in the list: the atmosphere is 2D on its own regular lon/lat grid and knows nothing about the
model's vertical coordinate.

The same list applies to a change in any `bathymetry_config` smoothing knob, since those change the
land mask that `prepare_forcing` and `prepare_boundaries` build their masks from. See "The vertical
grid" in `docs/setups.md` for what the current faces are and why.

`-m` is Julia 1.12's package entry point, which is why `Project.toml` has `[compat] julia = "1.12"`.
The equivalent without `-m` is
`julia --project -e 'using FjordSim; FjordSim.main(ARGS)' -- prepare_forcing --config oslofjorden`.

Each subcommand is a `FjordConfig` method on the generic function of the same name, so the same
steps run from the REPL, which is also how to debug one:
```julia
using Pkg; Pkg.activate(".")
using FjordSim
config = oslofjorden()
prepare_bathymetry(config)

using Debugger        # step through a step
@enter prepare_forcing(config)
```

`run_simulation` is the exception with a second entry point: `build_simulation(config)` returns
the assembled `Simulation` without starting it, which is what to reach for in the REPL or the
debugger.

## Architecture

FjordSim wraps [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and
[NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl) to set up regional ocean
simulations of Norwegian fjords. A simulation is assembled from a grid, a bathymetry file, forcing,
and atmospheric data.

**A new grid, bathymetry source, forcing dataset, river dataset, open-boundary dataset or
atmosphere dataset is added by subtyping the matching `Abstract*Config` supertype and overloading
methods on it — never by editing the generic pipeline.** Before doing that, read
**`docs/adding-a-source.md`** for the hook tables (what's required, what has a fallback, and why).

The modules, in `include` order from `src/FjordSim.jl` — **`docs/architecture.md` has the full
per-module writeup** (design rationale, ordering constraints, and the measured numbers behind
non-obvious choices like the bathymetry smoothing pipeline order or the open-boundary scheme set);
read it before changing any of these:

1. **Configs** (`src/Configs.jl`) — the abstract supertypes, `FjordConfig`, the edge vocabulary
   (`LATERAL_EDGES`, `open_edges`), and the path helpers every source inherits.
2. **Dataset adapters** (`src/Datasets.jl`) — dead code (unused and broken against NumericalEarth
   0.6); do not build on it.
3. **Utils** (`src/Utils.jl`) — `progress` callback, `recursive_merge`, time-step helpers.
4. **Plotting** (`src/Plotting.jl`) — one `plot_*` function per pipeline, dispatching on config
   supertypes.
5. **Bathymetry** (`src/Bathymetry/`) — `prepare_bathymetry` and the seven-stage
   `smooth_bathymetry_gaps!` pipeline (order is forced; see the doc before reordering).
6. **Atmosphere** (`src/Atmospheres/`) — `prepare_atmosphere`, NORA3 read/download.
7. **Forcing** (`src/Forcing/`) — interior forcing, rivers, and open-boundary data pipelines
   (`prepare_forcing`, `add_rivers`, `prepare_boundaries`).
8. **Boundary conditions** (`src/BoundaryConditions.jl`) — `AirSeaFluxes`, `QuadraticBottomDrag`,
   `OpenLateralBoundaryFromData`.
9. **Grids** (`src/Grids.jl`) — `EvenGrid`, the `domain_grid`/`simulation_grid` hooks.
10. **Simulations** (`src/Simulations.jl`) — `SimulationConfig`, `build_simulation`/`run_simulation`,
    writers, callbacks, looping, checkpointing, initial conditions.
11. **Setups** (`src/Setups/`) — the built-in fjords; see `docs/setups.md`.
12. **CLI** (`src/CLI.jl`) — subcommand dispatch, logging.
13. **Top-level** (`src/FjordSim.jl`) — re-exports and `main`.

## Setups

`src/Setups/` holds one lowercase file per fjord, each a zero-arg function returning a fresh
`FjordConfig`, registered in `SETUPS`. Adding a fjord is a two-place edit: the new file, and its
entry in `SETUPS`. Data paths default to `~/FjordSim_data/<fjord>/`; results to
`~/FjordSim_results/<fjord>/`. `SimulationConfig` and everything it nests have **no defaults** —
every field is a scientific choice a setup must state.

Read **`docs/setups.md`** before adding a fjord or tuning an existing one — it has the vertical-grid
rationale and the measured justification for `oslofjorden()`'s current parameter values (the "2020
diagnosis").

## GPU & Kernel Compatibility

Forcing callables (e.g. `ForcingFromFile`) and any function called inside an Oceananigans
kernel must follow these rules:

- Mark with `@inline`
- Use `ifelse` — never short-circuiting `&&`/`||` or ternary `?:` inside kernels
- Must be type-stable and allocation-free
- Never loop over grid points outside kernels — use `launch!` via KernelAbstractions
- Use literal zeros: `max(0, a)` not `max(zero(FT), a)`
- Never hardcode `Float64` literals (`0.0`, `1.0`) in kernels or constructors — use `zero(grid)`,
  `one(grid)`, `convert(FT, val)`, or rational literals like `1//2`
- Use `on_architecture` for CPU↔GPU data transfers — never call `Array()` or `CuArray()` directly
- Always index fields with three indices: `field[i, j, k]` — 2D indexing appears to work on some
  fields by coincidence but is unsupported and will break

`src/Forcing/Forcing.jl` already demonstrates the correct pattern; match it when adding new
forcings.

## Type Stability

- All structs must have concrete field types — no `Any`, no abstract field types
- User-facing constructors may accept flexible arguments, but the returned struct must be
  fully typed (follow the materialization pattern from NumericalEarth if needed)
- Type annotations are for dispatch only, not documentation
- Profile with `@code_warntype` when touching GPU-executed paths

## Import Conventions

- Extend functions via `function Mod.foo(...) ... end` — never `import Mod: foo` to extend
- Exports go at the top of each module file, before other code
- Keep imports explicit so the dependency surface stays auditable

## Code Style

- Avoid abbreviations: `latitude` not `lat`, `temperature` not `temp` (exception: physics
  symbols `T`, `S`, `u`, `v`, `λ`, `φ` are conventional)
- Keyword args: no-space inline `f(x=1)`, single-space multiline `f(\n    a = 1,\n    b = 2\n)`
- No trailing whitespace, no trailing blank lines; files end with exactly one newline
- Always use explicit `return` in functions longer than one expression
- Delete commented-out code — git is the history; no debugging artifacts or stale copy-paste
- Prefer dispatch over conditionals: multiple dispatch instead of `if`/`else` branching on types

## Common Pitfalls

- Never extend `getproperty` to fix undefined-property bugs — fix the caller instead
- A variable named the same as a function produces "type is not callable" — rename the variable
- Never add/remove/change `[deps]` in `Project.toml` on your own initiative; only touch
  `[compat]` when explicitly asked. But if a dependency would meaningfully simplify the
  implementation, *ask* — do not silently contort the design to avoid it. This applies to packages
  already present transitively (e.g. `KernelAbstractions` via Oceananigans), where a direct
  `[deps]` entry costs nothing but is still required to `using` them.
- `Pkg.resolve()` only *checks* the manifest against `[compat]`; it will not move a version to
  satisfy a widened bound, and reports `empty intersection between X@old and project compatibility`
  instead. Use `Pkg.update()` after raising a compat bound. (`Pkg.up` is not a name in current Pkg.)
- `Base.@kwdef` on a *parametric* struct does not convert its arguments, unlike on a
  non-parametric one: the generated constructor keeps the declared field types in its signature,
  so `SimulationConfig(stop_time = 3600, ...)` is a `MethodError` where `stop_time = 1hour` is
  fine. Every `Oceananigans.Units` constant is already a `Float64`, so writing durations with them —
  `365days`, `1hour`, `3minutes` — sidesteps this. `loops` is the one field this does *not* apply to,
  being an `Int` already. `SimulationConfig` is now the **only** config the rule reaches: every other
  one is either non-parametric, or — like `SnapshotWriter`, `ProgressCallback` and
  `CoupledHydrostaticSimulation`, all of which *are* parametric — has a hand-written keyword
  constructor that converts. `Forcing.PreparedVariable` is the one internal type the rule reaches
  too, for the same reason.

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land.
- Data files default to `~/FjordSim_data/<fjord>/` and results to `~/FjordSim_results/<fjord>/`,
  the latter from the simulation config's `results_root`.
- Open-boundary convention: the domain is open on any subset of its four lateral edges — none, one,
  or all four for a region in the open ocean — named once by the boundary data config's `open_edges`
  and read through the `open_edges` accessor as a `Vector{Symbol}`, empty for a setup naming no
  boundary config, whose lateral boundaries are then all closed walls. Everything that acts on an edge
  dispatches on `Val(edge)` and every consumer *iterates*, so the four edges are independent; the
  prepared boundary file names every variable for its side (`south_T`), which is how they all fit in
  one file on one time axis.
