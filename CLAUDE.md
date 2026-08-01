# CLAUDE.md

## Commands

```bash
# Run all tests
julia --project test/runtests.jl

# Run the Oslofjord example simulation (requires GPU + data files); same as run_simulation below
julia --project examples/oslofjord.jl

# Prepare bathymetry for a configured fjord (downloads the Geonorge GDB on first use)
julia --project -m FjordSim prepare_bathymetry --config drammensfjorden

# Download and subset the configured forcing dataset for a fjord
julia --project -m FjordSim download_forcing --config oslofjorden

# Regrid the downloaded forcing onto the fjord's grid (needs the bathymetry and download first)
julia --project -m FjordSim prepare_forcing --config oslofjorden

# Write river relaxation on top of the prepared forcing (needs prepare_forcing first)
julia --project -m FjordSim add_rivers --config oslofjorden

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

FjordSim is a Julia package that wraps [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl) to set up regional ocean simulations of Norwegian fjords.

A simulation is assembled from a grid, a bathymetry file, forcing, and atmospheric data. The
modules, in `include` order from `src/FjordSim.jl`:

1. **Configs** (`src/Configs.jl`) — the abstract supertypes `AbstractGridConfig`,
   `AbstractBathymetryConfig`, `AbstractForcingConfig`, `AbstractRiverConfig`,
   `AbstractAtmosphereConfig` and `AbstractSimulationConfig`, plus `FjordConfig`, which holds a
   grid, bathymetry, forcing, atmosphere and simulation config (parametrically, so every
   instantiation stays concretely typed). `atmosphere_config` and `simulation_config` default to
   `nothing`, so a setup opts in by naming one. A river config hangs off the forcing config's
   `rivers` field instead, where `nothing` means the setup has no rivers. A new grid, bathymetry
   source, forcing dataset, river dataset or atmosphere dataset is
   added by subtyping the matching supertype and overloading methods on it — `FjordConfig` and its
   callers are untouched. The path helpers defined on the supertypes live here too, so a new
   source inherits them without loading the built-in source's module: `bathymetry_path`,
   `forcing_path`, `forcing_directory`, `river_forcing_path`, `atmosphere_path`,
   `atmosphere_directory`, `results_path`, and `plot_path` (one method per supertype). Each
   supertype's docstring lists the fields and hook methods a subtype must provide — see "Adding a
   new source" below.

   `AbstractSimulationConfig` is the one supertype with fields but no hooks: `build_simulation` is
   generic over the whole `FjordConfig`, because everything dataset-specific already comes from
   the other configs. It is also the only one rooted at a `results_root` rather than a
   `data_root`, being the only config that writes rather than reads.

2. **Dataset adapters** (`src/Datasets.jl`) — `ForcingDataset` and `ResultsDataset` are NumericalEarth
   dataset wrappers for local FjordSim NetCDF files (forcing inputs and simulation outputs),
   used for initial conditions and restart.

3. **Utils** (`src/Utils.jl`) — `progress` callback, `recursive_merge` for nested boundary-condition
   named tuples, `cell_advection_timescale_coupled_model` for the time-step wizard, plus
   `compute_faces` and NetCDF/JLD2 helpers.

4. **Plotting** (`src/Plotting.jl`) — `plot_bathymetry(grid, bottom_height, config)`,
   `plot_forcing(grid, config)` and `plot_atmosphere(config)`, all dispatching on the config
   *supertypes* so a new source inherits them, writing to `plot_path(config)`. Shares
   `default_figure_size` and `plot_axes`. It is included *before* the pipelines rather than after
   them because each pipeline's setup-level driver plots as its last step, and `Plotting` itself
   only needs `Configs`, so there is no cycle.

5. **Bathymetry** (`src/Bathymetry/Bathymetry.jl` generic core, `src/Bathymetry/geonorge.jl`
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

6. **Atmosphere** (`src/Atmospheres/Atmospheres.jl` generic core, `src/Atmospheres/nora3_source.jl`
   dataset adapter included flat into the same module, `src/Atmospheres/NORA3.jl` a nested
   `module NORA3` holding the read side).

   The read side: `MultiYearNORA3` is a dataset wrapper for a prepared atmosphere NetCDF.
   `NORA3PrescribedAtmosphere` and `NORA3PrescribedRadiation` construct NumericalEarth
   `PrescribedAtmosphere`/`PrescribedRadiation` objects backed by `NORA3FieldTimeSeries`.
   `MultiYearNORA3(config)` resolves a setup's own prepared file through `atmosphere_path`;
   `default_nora3_dataset()` is the legacy shared `~/FjordSim_data/NORA3/NORA3.nc`.

   `prepare_atmosphere(target_grid, config::AbstractAtmosphereConfig)` regrids the downloaded files
   onto a regular lon/lat grid: `atmosphere_variable_names(config)` picks the variables →
   `atmosphere_time_steps(config)` gives the hourly `AtmosphereRecord`s →
   `atmosphere_source_grid(config, filepath)` gives the downloaded geometry →
   `atmosphere_target_axes` derives the prepared axes from `x_domain`/`y_domain` grown by
   `config.padding` and sampled at `config.resolution` → `projected_atmosphere_nodes` projects them
   into the source projection in one bulk GDAL call → one bilinear `interpolate_to_target!` per
   variable per step → streaming NetCDF write. `download_atmosphere(config::FjordConfig)` is the
   generic download driver, mirroring `download_forcing`.

   The prepared file's layout is a **contract**, fixed by the read side: `Float32` variables of
   shape `(lon, lat, time)`, uniformly spaced 1D `lon`/`lat` centers (`compute_faces` infers the
   spacing from the first difference), a CF-encoded `time`, and exactly the eight names and units
   in `ATMOSPHERE_VARIABLES`. Air temperature is **Kelvin** and both radiative fluxes are
   **downwelling**, because that is what NumericalEarth consumes — `PrescribedRadiation` applies
   its own surface albedo, so a net flux would count it twice.

   Interpolation is deliberately *not* reused from `Forcing`: `Atmospheres` is included before
   `Forcing`, and `Forcing`'s machinery is 3D, depth-oriented and mask-driven, whereas atmospheric
   fields are 2D and defined everywhere including over land. The cost is ~15 duplicated lines of the
   ArchGDAL transform in `projected_target_nodes`; the benefit is a self-contained module with no
   dummy depth level, no all-true mask, and no `architecture` field — the prepared grid is ~50x60,
   so there is nothing worth a GPU.

   `nora3_source.jl` holds `NORA3Config <: AbstractAtmosphereConfig` and the MET Norway NORA3
   download. NORA3 is served one file per forecast lead hour of a 6-hourly run, and the archive
   layout is deterministic, so no catalog listing is needed. The download owns everything that
   depends on the forecast structure and writes one gap-free hourly file per month already carrying
   the eight prepared names, which makes `prepare_atmosphere` a pure regrid — re-running it at a
   different `resolution` costs no re-download.

   Two non-obvious pieces. **Which lead supplies which hour**: leads 4..9 of one run cover six
   consecutive hours, so the four daily runs tile a day exactly with no overlap and no dedup pass —
   hours 04..09 from 00Z, 10..15 from 06Z, 16..21 from 12Z, 22..23 from 18Z, and 00..03 from the
   *previous* day's 18Z run. **De-accumulation**: the flux accumulators restart at every run
   (downwelling longwave for 09:00 reads 11_181_925 J/m² as lead 9 of the 00Z run but 3_734_820 J/m²
   as lead 3 of the 06Z run), so `process_run` walks leads 3..9 carrying the previous lead's values
   forward and every difference stays inside one run. Each increment is labelled at the *end* of its
   interval, which is what puts the fluxes on the same hourly axis as the instantaneous fields —
   the reference Python pipeline this ports labels them at the midpoint and so needs a second
   `time_acc` axis, which `MultiYearNORA3` cannot read.

   Winds arrive relative to the Lambert grid axes and are rotated on the native grid before
   anything is interpolated. `grid_rotation_angle` derives the angle by finite differences of the 2D
   longitude/latitude fields, so it carries over to any curvilinear source that publishes them; the
   analytic Lambert alternative would need the projection parameters. In the Oslofjord region it
   gives about -48°, matching the meridian convergence `(10.6 - (-42)) * sin(66.3°)`.

   met.no's aggregated `nora3_subset_atmos` collections are **not** usable here: they lack
   `specific_humidity_2m` and carry only *net* radiation, and net longwave cannot be inverted to
   downwelling.

7. **Forcing** (`src/Forcing/Forcing.jl` generic core, `src/Forcing/norkyst.jl` dataset adapter,
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
   `download_forcing(target_grid, config)`. `prepare_forcing(config::FjordConfig)` and
   `add_rivers(config::FjordConfig)` are the matching setup-level drivers for the other two steps;
   both take the grid from the processed bathymetry so the land mask matches the model, which is
   why they require `prepare_bathymetry` to have run.

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

8. **Boundary conditions** (`src/BoundaryConditions.jl`) — `top_bottom_boundary_conditions` creates
   wind/heat/salt flux fields at the top and quadratic bottom drag, returning a named tuple
   `(u, v, T, S)`.

9. **Grids** (`src/Grids.jl`) — `EvenGrid <: AbstractGridConfig` (size, halo, longitude, latitude,
   `z_faces`) with `LatitudeLongitudeGrid(architecture, config::EvenGrid)`, and a constructor
   `ImmersedBoundaryGrid(filepath, architecture, halo)` that reads the processed bathymetry NetCDF and
   returns an `ImmersedBoundaryGrid` wrapping a `LatitudeLongitudeGrid` with `PartialCellBottom`.
   The loader still accepts legacy files with positive depths or swapped `lon`/`lat` axes.

10. **Simulations** (`src/Simulations.jl`) — `SimulationConfig <: AbstractSimulationConfig`, the
    `build_simulation`/`run_simulation` drivers, and `coupled_hydrostatic_simulation`, which
    assembles a `HydrostaticFreeSurfaceModel` inside an `OceanSeaIceModel` (NumericalEarth) and
    returns a `Simulation`. Included after `Grids`, which it reads the grid back through, and
    before `Setups`, which constructs a `SimulationConfig`.

    `coupled_hydrostatic_simulation` lives here rather than in `src/FjordSim.jl` for that ordering
    reason: a function defined after the `include` block cannot be imported by a module included
    inside it. Its 15-positional signature is unchanged by the move, and `FjordSim` still exports
    it.

    `SimulationConfig`'s fields split three ways, and the split is forced rather than stylistic.
    `buoyancy`, `closure`, `tracer_advection`, `momentum_advection`, `tracers`,
    `initial_conditions`, `coriolis`, `sea_ice` and `biogeochemistry` depend on neither the grid
    nor the device, so they are stored as the objects `coupled_hydrostatic_simulation` consumes —
    one type parameter each, nine in total, which is the cost of the no-abstract-fields rule. A
    tenth is the signal to collapse them into one `NamedTuple` field instead of adding another,
    which the "Simulation config" testset guards with an exact parameter count.
    `free_surface_cfl` and `bottom_drag_coefficient` are scalars because
    `SplitExplicitFreeSurface` and `top_bottom_boundary_conditions` both need the grid.
    `architecture` is a `Symbol` for the same reason as the forcing config's, and
    `simulation_architecture` resolves it by reusing `interpolation_architecture`'s `Val` methods.

    No field has a default, deliberately, which is the one place `SimulationConfig` departs from
    the other configs: a defaulted closure or coriolis would be one fjord's physics quietly
    applied to another's. The setup file is the complete statement of the run.

    Three things are deliberately *not* fields, because they are already stated elsewhere and a
    second copy could only disagree: the forcing file (`simulation_forcing_path` takes the
    rivers-augmented copy when the forcing config names rivers, the plain prepared file
    otherwise), the open boundary (`open_boundary_conditions` puts it on the velocity normal to
    the forcing config's `relaxation_edge`), and the atmosphere (the `prescribed_atmosphere` and
    `prescribed_radiation` hooks on `config.atmosphere_config`).

    `build_simulation` returns the instrumented `Simulation` without running it — the REPL and
    debugger entry point — and `run_simulation` calls it and `run!`. Both return `nothing` for a
    setup naming no simulation config. Each prerequisite is checked before anything is read or
    allocated and reported as the command that produces it; the atmosphere is the exception,
    because `MultiYearNORA3(config)` already owns that error, which is what keeps the module free
    of any named dataset.

11. **Setups** (`src/Setups/Setups.jl` registry, one lowercase file per fjord beside it) — the
    built-in fjords, each a zero-arg function returning a fresh `FjordConfig`: `oslofjorden()`,
    `drammensfjorden()`. `SETUPS` maps a name to its function, `setup_names()` lists them sorted,
    and `fjord_config(name_or_path)` resolves either a registered name or a path to an out-of-tree
    `.jl` config file (evaluated in `Main`). Included after `Grids` because it constructs every
    config type, `EvenGrid` included.

    A setup is a *function*, not a `const`, for two reasons that both fail silently otherwise:
    `DybdedataConfig`'s `raw_directory` defaults to the scratch path `Bathymetry.__init__` fills
    in, which is still `""` during precompilation; and the config structs are mutable and
    `native_region!` edits the bathymetry config, so a shared instance would leak state between
    steps.

12. **CLI** (`src/CLI.jl`) — `SUBCOMMANDS` maps each subcommand to the driver it calls, one `USAGE`
    string, a pure `parse_arguments` returning `(; subcommand, config, help)` and throwing
    `ArgumentError`, and `main(args)` returning a process exit code. Included last, since it names
    every driver and every setup. `parse_arguments` deliberately does not `exit`: printing and exit
    codes live in `main`, which keeps the help path testable. Exit codes are 0 success, 1 the step
    failed, 2 bad arguments.

    `main` runs the driver inside `tee_output(f, LOG_FILE)`, which mirrors `stdout` and `stderr` to
    `fjordsim.log` in the working directory while still printing live to the terminal — a stacktrace
    through `SimulationConfig` spells out every type parameter and is long enough to push the error
    message itself out of the scrollback. `redirect_stdio` only accepts fd-backed streams, so the tee
    is a `Pipe` (needing an explicit `Base.link_pipe!`; an uninitialized `Pipe` throws from `eof`)
    plus a task copying each chunk to both destinations. Both streams share one pipe, so the log
    interleaves them in write order.

    Two things about it are load-bearing. The failure is **caught inside** the redirect and reported
    with `showerror`: an exception left to propagate would be printed by `Base._start` after
    `tee_output`'s `finally` had torn the redirect down, so the error would be the one thing missing
    from the log. And only the driver runs inside the tee — parsing, `--help` and config resolution
    stay outside it, so a usage error leaves no log file behind, which is also why the existing
    `main` tests (all of which fail before the driver is reached) write nothing.

    The cost is that `stdout` is a `Pipe` during a run, so `displaysize` reports the 24x80 default
    and wide `show` output (Oceananigans' grid and model summaries) may wrap at 80 columns. Colour
    survives, since `Base.have_color` is fixed at startup and the raw bytes reach the real terminal.
    On a single-threaded Julia the reader task only runs when the main task yields, but writing to
    the pipe is libuv I/O and does yield, and the `finally` waits for the reader to drain, so output
    can arrive in bursts but is never lost.

13. **Top-level** (`src/FjordSim.jl`) — re-exports the public API and defines `main` plus a bare
    `@main` for `julia -m FjordSim`. Also patches `compute_bounding_indices` from NumericalEarth to
    prevent off-by-one errors with custom longitude/latitude grids.

    Two non-obvious things about the entry point. `main` is **not exported**: Julia's startup runs
    `Main.main` after a script's body whenever that binding resolves to an entry point, so
    exporting it would make every `using FjordSim` in a script — `test/runtests.jl`,
    `examples/oslofjord.jl` — run the CLI on the way out. And the `@main` is bare, *after* the
    definition: `@main function main(args) ... end` expands to a **call**, which would run the CLI
    while the package precompiles.

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

Atmosphere — `AbstractAtmosphereConfig`, consumed by `prepare_atmosphere`:

| Hook | Required |
|---|---|
| `atmosphere_time_steps(config)` → `Vector{AtmosphereRecord}` | yes |
| `atmosphere_source_grid(config, filepath)` → source grid, e.g. `ProjectedAtmosphereGrid` | yes |
| `atmosphere_variable_names(config)` → `Dict` downloaded name => prepared name | yes |
| `download_atmosphere(target_grid, config)` | only if it downloads |
| `prescribed_atmosphere(config, architecture)` → `PrescribedAtmosphere` | only if the setup is simulated |
| `prescribed_radiation(config, architecture)` → `PrescribedRadiation` | only if the setup is simulated |

The last two are the read side, consumed by `build_simulation` rather than `prepare_atmosphere`,
and are what keep `Simulations` from naming a dataset. They take no float type on purpose: both
NumericalEarth constructors default to `Float32`, and passing `Oceananigans.defaults.FloatType`
would silently promote the atmosphere to `Float64`. A `nothing` atmosphere config yields `nothing`
for both.

The prepared variable names and units are *not* a hook — they are fixed by the read side in
`ATMOSPHERE_VARIABLES`. A source whose download already normalizes names (as NORA3's does, since
five of its eight variables are derived rather than copied) returns the identity mapping. The
core's `grid_rotation_angle` and `rotate_to_east_north` are available to any adapter whose source
gives wind relative to its own grid axes.

Simulation — `AbstractSimulationConfig`, consumed by `build_simulation`: **no hooks**. It is data
only, read field by field, so an alternative simulation config is a subtype supplying the field
set its docstring lists. It inherits `results_path`.

Required fields are listed in each supertype's docstring in `src/Configs.jl`. Path resolution
(`bathymetry_path`, `forcing_path`, `forcing_directory`, `river_forcing_path`, `atmosphere_path`,
`atmosphere_directory`, `results_path`, `plot_path`) and the plots come for free. A missing
required hook surfaces as a `MethodError` naming it, which the "Config extensibility" testset
asserts.

## Setups

`src/Setups/` holds one lowercase file per fjord (`oslofjorden.jl`, `drammensfjorden.jl`), each
defining a zero-arg function that returns a fresh `FjordConfig`, plus `Setups.jl` with the `SETUPS`
registry. Adding a fjord is a **two-place** edit: the new file, and its entry in `SETUPS` — the
registry is keyed by a runtime string from `--config`, so there is no dispatch alternative. A fjord
that should not live in the package can instead be a standalone `.jl` file whose last expression is
a `FjordConfig`, passed as `--config path/to/it.jl` and loaded by `fjord_config`.

Data paths are built from a per-setup `data_root` under `~/FjordSim_data/<fjord>/`, computed inside
the setup function with `homedir()`; the config fields naming files (`output_file`, `plot_file`,
`geodatabase_file`, `output_directory`) are names relative to `data_root`, and setting one to an
absolute path relocates just that file — which is how a single FileGDB copy is shared across fjords.
A nested config carries its own `data_root` too, so it can be relocated independently, but
`oslofjorden()` gives its `OF800RiversConfig` the same `data_root` as the rest of the setup —
the river data downloads there rather than being shared from elsewhere. `drammensfjorden()`
names no `rivers`, so it defaults to `nothing` and the rivers step is a no-op — and it names
neither an `atmosphere_config` nor a `simulation_config` either, which default to `nothing` and
make both atmosphere steps and `run_simulation` no-ops the same way.

The simulation config is rooted separately, at `~/FjordSim_results/<fjord>/` rather than under
`data_root`, since it writes rather than reads. Unlike every other config, `SimulationConfig` has
**no defaults at all**: `oslofjorden()` names all 21 fields, because each is a scientific choice
about that fjord and a default would let the next setup silently inherit it. Adding a field to
`SimulationConfig` therefore breaks every setup until each names it, which is the intent.

Each step is a `FjordConfig` method on the generic function of the same name —
`prepare_bathymetry`, `download_forcing`, `prepare_forcing`, `add_rivers`, `download_atmosphere`,
`prepare_atmosphere`, `run_simulation` — living beside the pipeline it drives. Each builds the
grid, checks the step before it has run, calls the generic pipeline, plots, and logs where the
output went. A step the setup opts out of returns `nothing` rather than raising, matching the
`::Nothing` methods of the lower arities; "you asked for a step this setup does not configure" is
reported by `CLI.main`, because that is user input rather than a pipeline condition.

`examples/oslofjord.jl` is the end-to-end simulation script and is now just
`run_simulation(oslofjorden())` — everything it used to wire by hand is the setup's
`SimulationConfig`.

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
- `Base.@kwdef` on a *parametric* struct does not convert its arguments, unlike on a
  non-parametric one: the generated constructor keeps the declared field types in its signature,
  so `SimulationConfig(stop_time = 3600, ...)` is a `MethodError` where `stop_time = 1hour` is
  fine. Every `Oceananigans.Units` constant is already a `Float64`, so writing durations with
  them — `365days`, `1hour`, `3minutes` — sidesteps this. Every other config is non-parametric
  and so has never hit it.

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land.
- Data files default to `~/FjordSim_data/<fjord>/` and results to `~/FjordSim_results/<fjord>/`,
  the latter from the simulation config's `results_root`.
