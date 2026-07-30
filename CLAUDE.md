# CLAUDE.md

## Commands

```bash
# Run all tests
julia --project test/runtests.jl

# Run the Oslofjord example simulation (requires GPU + data files)
julia --project examples/oslofjord.jl

# Prepare bathymetry for a configured fjord (downloads the Geonorge GDB on first use)
julia --project scripts/bathymetry_prepare.jl --config configs/drammensfjorden.jl

# Download and subset NorKyst-800m forcing for a configured fjord
julia --project scripts/forcing_download_norkyst.jl --config configs/oslofjorden.jl
```

`--config` is required by both scripts — there is no default setup.

Activate the environment before running scripts interactively:
```julia
using Pkg; Pkg.activate(".")
```

## Architecture

FjordSim is a Julia package that wraps [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl) to set up regional ocean simulations of Norwegian fjords.

A simulation is assembled from a grid, a bathymetry file, forcing, and atmospheric data. The
modules, in `include` order from `src/FjordSim.jl`:

1. **Configs** (`src/Configs.jl`) — the three abstract supertypes `AbstractGridConfig`,
   `AbstractBathymetryConfig`, `AbstractForcingConfig`, plus `FjordConfig`, which holds one of
   each (parametrically, so every instantiation stays concretely typed). A new grid,
   bathymetry source or forcing dataset is added by subtyping the matching supertype and
   overloading methods on it — `FjordConfig` and its callers are untouched.

2. **Dataset adapters** (`src/FDatasets.jl`) — `DSForcing` and `DSResults` are NumericalEarth
   dataset wrappers for local FjordSim NetCDF files (forcing inputs and simulation outputs),
   used for initial conditions and restart.

3. **Utils** (`src/Utils.jl`) — `progress` callback, `recursive_merge` for nested boundary-condition
   named tuples, `cell_advection_timescale_coupled_model` for the time-step wizard, plus
   `compute_faces` and NetCDF/JLD2 helpers.

4. **Bathymetry** (`src/Bathymetry.jl`) — `DybdedataConfig <: AbstractBathymetryConfig` describes
   one setup's Geonorge Sjøkart Dybdedata bathymetry. `prepare_geonorge_bathymetry(target_grid, config)`
   runs the pipeline: derive the native region from the target grid (`native_region!`) → download
   and extract the FileGDB if absent (`ensure_geodatabase`, ~2.3 GB) → sample `dybdepunkt` points and
   optionally `dybdekurve` contours in EPSG:25833 → transform the point cloud to WGS84 in one bulk
   GDAL call → grid to a raster and burn in land features (`landareal`, `skjer`) →
   `NumericalEarth.regrid_bathymetry` → `smooth_bathymetry_gaps!` → processed NetCDF. The
   intermediate raw dataset is cached in `Scratch` storage and wrapped by
   `GeonorgeBathymetry <: AbstractStaticBathymetry`, which implements the NumericalEarth dataset
   interface. `bathymetry_path`/`plot_path` are defined on `AbstractBathymetryConfig`, so a new
   bathymetry source inherits path resolution.

5. **Atmosphere** (`src/Atmospheres/Atmospheres.jl`, `src/Atmospheres/NORA3.jl`) — `MultiYearNORA3` is
   a dataset wrapper for NORA3 reanalysis NetCDF. `NORA3PrescribedAtmosphere` and
   `NORA3PrescribedRadiation` construct NumericalEarth `PrescribedAtmosphere`/`PrescribedRadiation`
   objects backed by `NORA3FieldTimeSeries`. Default data path: `~/FjordSim_data/NORA3/NORA3.nc`.

6. **Forcing** (`src/Forcing.jl`) — loads river/relaxation forcing from NetCDF via
   `forcing_from_file`. The `ForcingFromFile` struct carries two `FieldTimeSeries`
   (values + lambdas) and dispatches to flux, advection, or relaxation terms based on the sign of
   lambda: `λ > 1` → x-flux, `λ < -1` → y-flux, `|λ| < 1` → relaxation. Uses a custom
   `NetCDFBackend` keeping 2 time indices in memory. Also holds
   `NorKystConfig <: AbstractForcingConfig` (THREDDS endpoints, variables, years) and
   `norkyst_directory`/`norkyst_monthly_filename`; the download itself lives in
   `scripts/forcing_download_norkyst.jl`.

7. **Boundary conditions** (`src/BoundaryConditions.jl`) — `top_bottom_boundary_conditions` creates
   wind/heat/salt flux fields at the top and quadratic bottom drag, returning a named tuple
   `(u, v, T, S)`.

8. **Grids** (`src/Grids.jl`) — `EvenGrid <: AbstractGridConfig` (size, halo, longitude, latitude,
   `z_faces`) with `LatitudeLongitudeGrid(arch, config::EvenGrid)`, and a constructor
   `ImmersedBoundaryGrid(filepath, arch, halo)` that reads the processed bathymetry NetCDF and
   returns an `ImmersedBoundaryGrid` wrapping a `LatitudeLongitudeGrid` with `PartialCellBottom`.
   The loader still accepts legacy files with positive depths or swapped `lon`/`lat` axes.

9. **Top-level** (`src/FjordSim.jl`) — re-exports the public API and defines
   `coupled_hydrostatic_simulation`, which assembles a `HydrostaticFreeSurfaceModel` inside an
   `OceanSeaIceModel` (NumericalEarth) and returns a `Simulation`. Also patches
   `compute_bounding_indices` from NumericalEarth to prevent off-by-one errors with custom
   longitude/latitude grids.

## Setups

`configs/` holds one Julia file per fjord (`oslofjorden.jl`, `drammensfjorden.jl`). Each file is a
script whose last expression is a `FjordConfig`, so scripts load it with
`config = include(abspath(path))`. Data paths are built from a per-setup `data_root` under
`~/FjordSim_data/<fjord>/`; the config fields naming files (`output_file`, `plot_file`,
`geodatabase_file`, `output_directory`) are names relative to `data_root`, and setting one to an
absolute path relocates just that file — which is how a single FileGDB copy is shared across fjords.

`examples/oslofjord.jl` is the end-to-end simulation script and does *not* go through
`FjordConfig`; it builds the grid straight from a bathymetry NetCDF and wires the components
together by hand.

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

`src/Forcing.jl` already demonstrates the correct pattern; match it when adding new forcings.

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
- Never add/remove/change `[deps]` in `Project.toml` unless the task requires it; only touch
  `[compat]` when explicitly asked

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land.
- Data files default to `~/FjordSim_data/` and results to `~/FjordSim_results/`.
