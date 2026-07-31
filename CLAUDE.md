# CLAUDE.md

## Commands

```bash
# Run all tests
julia --project test/runtests.jl

# Run the Oslofjord example simulation (requires GPU + data files)
julia --project examples/oslofjord.jl

# Prepare bathymetry for a configured fjord (downloads the Geonorge GDB on first use)
julia --project scripts/bathymetry_prepare.jl --config configs/drammensfjorden.jl

# Download and subset the configured forcing dataset for a fjord
julia --project scripts/forcing_download.jl --config configs/oslofjorden.jl

# Regrid the downloaded forcing onto the fjord's grid (needs the bathymetry and download first)
julia --project scripts/forcing_prepare.jl --config configs/oslofjorden.jl

# Write river relaxation on top of the prepared forcing (needs forcing_prepare first)
julia --project scripts/forcing_add_rivers.jl --config configs/oslofjorden.jl
```

`--config` is the only option any script takes, and it is required — there is no default setup.
Every other knob, including which device the forcing interpolation runs on, is a config field.

Activate the environment before running scripts interactively:
```julia
using Pkg; Pkg.activate(".")
```

## Architecture

FjordSim is a Julia package that wraps [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl) to set up regional ocean simulations of Norwegian fjords.

A simulation is assembled from a grid, a bathymetry file, forcing, and atmospheric data. The
modules, in `include` order from `src/FjordSim.jl`:

1. **Configs** (`src/Configs.jl`) — the abstract supertypes `AbstractGridConfig`,
   `AbstractBathymetryConfig`, `AbstractForcingConfig` and `AbstractRiverConfig`, plus
   `FjordConfig`, which holds one each of the first three (parametrically, so every
   instantiation stays concretely typed). A river config hangs off the forcing config's
   `rivers` field instead, where `nothing` means the setup has no rivers. A new grid,
   bathymetry source, forcing dataset or river dataset is added by subtyping the matching
   supertype and overloading methods on it — `FjordConfig` and its callers are untouched. The
   path helpers defined on the supertypes live here too, so a new source inherits them without
   loading the built-in source's module: `bathymetry_path`, `forcing_path`,
   `forcing_directory`, `river_forcing_path`, and `plot_path` (one method per supertype). Each
   supertype's docstring lists the fields and hook methods a subtype must provide — see
   "Adding a new source" below.

2. **Dataset adapters** (`src/Datasets.jl`) — `ForcingDataset` and `ResultsDataset` are NumericalEarth
   dataset wrappers for local FjordSim NetCDF files (forcing inputs and simulation outputs),
   used for initial conditions and restart.

3. **Utils** (`src/Utils.jl`) — `progress` callback, `recursive_merge` for nested boundary-condition
   named tuples, `cell_advection_timescale_coupled_model` for the time-step wizard, plus
   `compute_faces` and NetCDF/JLD2 helpers.

4. **Bathymetry** (`src/Bathymetry/Bathymetry.jl` generic core, `src/Bathymetry/geonorge.jl`
   source adapter, included into the same `Bathymetry` module).
   `prepare_bathymetry(target_grid, config::AbstractBathymetryConfig)` is the generic pipeline:
   `bathymetry_dataset(target_grid, config)` (the one source-specific step) →
   `NumericalEarth.regrid_bathymetry` with `regrid_options(config)` → `smooth_bathymetry_gaps!` →
   `write_bathymetry_file`. The core also owns the smoothing kernels and the
   `center_coordinates`/`expand_domain`/`vertical_faces` domain helpers.

   `geonorge.jl` holds `DybdedataConfig <: AbstractBathymetryConfig` and the Geonorge Sjøkart
   Dybdedata implementation of the two hooks: derive the native region from the target grid
   (`native_region!`) → download and extract the FileGDB if absent (`ensure_geodatabase`,
   ~2.3 GB) → sample `dybdepunkt` points and optionally `dybdekurve` contours in EPSG:25833 →
   transform the point cloud to WGS84 in one bulk GDAL call → grid to a raster and burn in land
   features (`landareal`, `skjer`). The intermediate raw dataset is cached in `Scratch` storage
   and wrapped by `GeonorgeBathymetry <: AbstractStaticBathymetry`, which implements the
   NumericalEarth dataset interface.

5. **Atmosphere** (`src/Atmospheres/Atmospheres.jl`, `src/Atmospheres/NORA3.jl`) — `MultiYearNORA3` is
   a dataset wrapper for NORA3 reanalysis NetCDF. `NORA3PrescribedAtmosphere` and
   `NORA3PrescribedRadiation` construct NumericalEarth `PrescribedAtmosphere`/`PrescribedRadiation`
   objects backed by `NORA3FieldTimeSeries`. Default data path: `~/FjordSim_data/NORA3/NORA3.nc`.

6. **Forcing** (`src/Forcing/Forcing.jl` generic core, `src/Forcing/norkyst.jl` dataset adapter,
   included into the same `Forcing` module).

   The read side loads river/relaxation forcing from NetCDF via `forcing_from_file`. The
   `ForcingFromFile` struct carries two `FieldTimeSeries` (values + lambdas) and dispatches to
   flux, advection, or relaxation terms based on the sign of lambda: `λ > 1` → x-flux,
   `λ < -1` → y-flux, `|λ| < 1` → relaxation. Uses a custom `NetCDFBackend` keeping 2 time
   indices in memory. `forcing_from_file` takes either a `filepath` keyword or an
   `AbstractForcingConfig` positionally, resolved by `forcing_path`.

   `prepare_forcing(target_grid, config::AbstractForcingConfig)` regrids the downloaded source
   files onto the simulation grid: `forcing_variable_names(config)` picks the variables →
   `forcing_time_steps(config)` completed by `daily_time_steps` to a gap-free daily axis →
   `forcing_source_grid(config, filepath)` → land mask from `peripheral_node` (so it matches the
   model, including `PartialCellBottom` and the velocity-face wall convention), with the open
   `relaxation_edge` row restored → source mask filled from the nearest valid cell (`SourceFill`)
   → one trilinear `Oceananigans.Fields.interpolate` per target cell in a `launch!` kernel,
   against the source subset expressed as a `RectilinearGrid` in projected meters
   (`ProjectedSourceGrid`, `source_field_grid`) → relaxation lambdas along `relaxation_edge` →
   streaming NetCDF write. Only the three hooks are dataset-specific; the rest is shared.

   Where the interpolation kernel runs is the config's `architecture` field (`:auto`, `:cpu`,
   `:gpu`), resolved by `interpolation_architecture` — `:auto` picks the GPU when
   `CUDA.functional()`, which is ~12x faster than a single-threaded `CPU()`, and `:gpu` errors
   rather than silently falling back. `target_grid` must stay on the CPU regardless, because
   building the masks walks `peripheral_node` cell by cell. A `Symbol` rather than a live
   `CPU()`/`GPU()` keeps config field types concrete and the config files loadable on a machine
   with no GPU.

   `download_forcing(config::FjordConfig)` is the generic download driver: it builds the setup's
   grid on the CPU and dispatches on the forcing config, so a dataset only implements
   `download_forcing(target_grid, config)`. `scripts/forcing_download.jl` is a thin CLI over it.

   `norkyst.jl` holds `NorKystConfig <: AbstractForcingConfig` (THREDDS endpoints, variables,
   years, output names, relaxation zone), `forcing_monthly_filename`, the three hook methods, and
   the NorKyst download: list the THREDDS catalog → open the month's OPeNDAP datasets → subset to
   the target grid's lon/lat box (`subset_ranges`, `NorKystSubset`) → write one combined monthly
   NetCDF.

   Two library functions deliberately *not* used, documented in `prepare_forcing` and
   `SourceFill`: NumericalEarth's dataset path (`native_grid`, `set!(field, metadata)`,
   `DatasetRestoring`) always builds a `LatitudeLongitudeGrid`, but the NorKyst grid is rotated
   ~59° from east here; and `inpaint_mask!` cannot fill a fully masked depth level, so on this
   regional subset it either never terminates or silently writes zeros.

   `rivers.jl` (generic core) and `of800_rivers.jl` (dataset adapter) are included into the same
   `Forcing` module and hold the rivers step, which runs *after* `prepare_forcing`.
   `add_rivers(target_grid, config::AbstractForcingConfig)` dispatches on `config.rivers` — a
   `nothing` river config is a no-op, so a setup opts in by naming one. The pipeline:
   `river_locations(rivers)` → snap each outlet to a grid cell (`river_cells`) → river values
   for the forcing file's own time axis (`river_series`) → copy the forcing file to
   `river_forcing_path(rivers)` → patch the copy's surface level (`write_rivers`). The original
   forcing file is never modified, so the step is re-runnable without redoing `forcing_prepare`.

   Rivers enter as **relaxation, not as a mass flux**: each river cell gets its value and
   `λ = 1 / relaxation_timescale` (1 hour by default) at the surface level for every time step,
   which lands in the existing `|λ| < 1` regime — no new forcing term or λ convention. A river
   cell inside the boundary relaxation band overrides that band's λ. Outlets are located by
   independent nearest-node lookups in longitude and latitude, then moved to the nearest
   *coastal* water cell — water with at least one land neighbour, so a river cannot be injected
   into open water. An outlet outside the grid, or with no coastal cell within `search_radius`,
   is dropped with a warning rather than written into land. The water mask comes from the same
   `water_mask` that `prepare_forcing` uses, so "water" means the same thing in both.
   `write_rivers` reads, patches and writes back whole surface slabs because the file is
   chunked one horizontal slab per `(level, time)`.

   `of800_rivers.jl` holds `OF800RiversConfig <: AbstractRiverConfig` and reads two files: an
   outlet CSV (hand-parsed; the project carries no CSV reader and the file is a couple of dozen
   rows) and a ROMS river NetCDF whose values repeat across `s_rho`, so only the top level is
   read and whose `river` coordinate is offset by `OF800_RIVER_NUMBER_OFFSET` from the CSV's
   numbering. `river_transport`, `river_Vshape` and `river_direction` are present in that file
   but deliberately unused — the reference pipeline this ports does not do a discharge
   conversion. `download_rivers` fetches both files from the per-file Dropbox links defaulted on
   the config, skipping any file already present. Those links must be per-*file* `/scl/fi/` URLs
   ending in `dl=1`: Dropbox renders folder listings client-side and serves a `/scl/fo/` folder
   link as one whole-folder archive, and without `dl=1` it serves the web app. A link that is not
   publicly shared answers HTTP 200 with a login page rather than failing, so
   `validate_river_download` rejects and deletes a downloaded file that starts with `<` instead
   of letting an HTML page masquerade as the data.

7. **Boundary conditions** (`src/BoundaryConditions.jl`) — `top_bottom_boundary_conditions` creates
   wind/heat/salt flux fields at the top and quadratic bottom drag, returning a named tuple
   `(u, v, T, S)`.

8. **Grids** (`src/Grids.jl`) — `EvenGrid <: AbstractGridConfig` (size, halo, longitude, latitude,
   `z_faces`) with `LatitudeLongitudeGrid(architecture, config::EvenGrid)`, and a constructor
   `ImmersedBoundaryGrid(filepath, architecture, halo)` that reads the processed bathymetry NetCDF and
   returns an `ImmersedBoundaryGrid` wrapping a `LatitudeLongitudeGrid` with `PartialCellBottom`.
   The loader still accepts legacy files with positive depths or swapped `lon`/`lat` axes.

9. **Plotting** (`src/Plotting.jl`) — `plot_bathymetry(grid, bottom_height, config)` and
   `plot_forcing(grid, config)`, both dispatching on the config *supertypes* so a new source
   inherits them, writing to `plot_path(config)`. Shares `default_figure_size` and `plot_axes`.

10. **Top-level** (`src/FjordSim.jl`) — re-exports the public API and defines
    `coupled_hydrostatic_simulation`, which assembles a `HydrostaticFreeSurfaceModel` inside an
    `OceanSeaIceModel` (NumericalEarth) and returns a `Simulation`. Also patches
    `compute_bounding_indices` from NumericalEarth to prevent off-by-one errors with custom
    longitude/latitude grids.

## Adding a new source

Every pipeline is a generic function on the config supertype plus a small set of hooks. Add a
source by subtyping and overloading the hooks — never by editing the generic function. The
adapter files (`src/Bathymetry/geonorge.jl`, `src/Forcing/norkyst.jl`, `src/Forcing/of800_rivers.jl`)
are the templates.

Grid — `AbstractGridConfig`:

| Hook | Required |
|---|---|
| `LatitudeLongitudeGrid(architecture, config)` | yes |

Bathymetry — `AbstractBathymetryConfig`, consumed by `prepare_bathymetry`:

| Hook | Required | Default |
|---|---|---|
| `bathymetry_dataset(target_grid, config)` → NumericalEarth dataset | yes | none |
| `regrid_options(config)` → NamedTuple for `regrid_bathymetry` | no | `(;)` |

Forcing — `AbstractForcingConfig`, consumed by `prepare_forcing`:

| Hook | Required |
|---|---|
| `forcing_time_steps(config)` → `Vector{SourceRecord}` | yes |
| `forcing_source_grid(config, filepath)` → source grid | yes |
| `forcing_variable_names(config)` → `Dict` source name => FjordSim name | yes |
| `download_forcing(target_grid, config)` | only if it downloads |
| `source_field_grid(source, architecture)`, `projected_target_nodes(longitude, latitude, source)` | only for a source grid that is not a regular projected grid; dispatch on the source-grid type, not the config |

Rivers — `AbstractRiverConfig`, consumed by `add_rivers`:

| Hook | Required | Default |
|---|---|---|
| `river_locations(config)` → `Vector{RiverLocation}` | yes | none |
| `river_series(config, times)` → `Dict` FjordSim name => `(river, time)` matrix | yes | none |
| `download_rivers(config)` | only if it downloads | none |
| `river_search_radius(config)` → cells to search for a coastal cell | no | `config.search_radius` |

A river config is not a `FjordConfig` field — it goes in the forcing config's `rivers` field,
`nothing` for a setup with no rivers. A variable `river_series` returns that the forcing file
does not carry is skipped with a warning, so one river dataset can serve setups that prepare
different variables.

Required fields are listed in each supertype's docstring in `src/Configs.jl`. Path resolution
(`bathymetry_path`, `forcing_path`, `forcing_directory`, `river_forcing_path`, `plot_path`) and
the plots come for free. A missing required hook surfaces as a `MethodError` naming it, which
the "Config extensibility" testset asserts.

## Setups

`configs/` holds one Julia file per fjord (`oslofjorden.jl`, `drammensfjorden.jl`). Each file is a
script whose last expression is a `FjordConfig`, so scripts load it with
`config = include(abspath(path))`. Data paths are built from a per-setup `data_root` under
`~/FjordSim_data/<fjord>/`; the config fields naming files (`output_file`, `plot_file`,
`geodatabase_file`, `output_directory`) are names relative to `data_root`, and setting one to an
absolute path relocates just that file — which is how a single FileGDB copy is shared across fjords.
A nested config carries its own `data_root` too, so it can be relocated independently, but
`oslofjorden.jl` gives its `OF800RiversConfig` the same `_data_root` as the rest of the setup —
the river data downloads there rather than being shared from elsewhere. `drammensfjorden.jl`
names no `rivers`, so it defaults to `nothing` and the rivers step is a no-op.

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
- `Pkg.resolve()` currently fails in this environment on an unrelated pinned `CUDACore` version.
  Adding a `[deps]` entry for a package already in `Manifest.toml` works without resolving; do not
  try to fix the resolve failure as a side effect.

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land.
- Data files default to `~/FjordSim_data/` and results to `~/FjordSim_results/`.
