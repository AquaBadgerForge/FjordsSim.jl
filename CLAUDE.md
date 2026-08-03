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

### Changing `start_date` or `stop_time`

Both prepare steps pad their time axes to span the run window the simulation config names
(`coverage_window`), so moving either end invalidates the prepared files. Re-run, in this order:

```bash
julia --project -m FjordSim prepare_forcing     --config oslofjorden   # rewrites forcing.nc
julia --project -m FjordSim add_rivers          --config oslofjorden   # re-copies forcing_rivers.nc
julia --project -m FjordSim prepare_atmosphere  --config oslofjorden   # rewrites atmosphere.nc
```

`add_rivers` is not optional here: it `cp`s the forcing file and patches the copy, so
`forcing_rivers.nc` — which is what `simulation_forcing_path` gives the simulation — still carries
the *old* axis until it is re-run. Neither download step is affected; both prepares are pure regrids.
`loops` is not part of the window, so changing it re-runs nothing.

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

   Three helpers on `AbstractSimulationConfig` read `start_date`, which is therefore as much part of
   the supertype's field set as `results_root` is. `run_tag` is the run's identity as a filename
   fragment (`yyyymmddTHHMM`); `results_path(config)` inserts it before `output_file`'s extension and
   `results_path(config, loop)` appends a zero-padded loop index on top; `coverage_window` is the
   `(first, last)` calendar interval the run needs its prepared inputs to span, which
   `prepare_forcing` and `prepare_atmosphere` take as their `coverage`.

   The tag is derived from the *simulated* start instant, not the wall clock, so a config always
   names the same files: post-processing can predict them, and — load-bearing — a `pickup` appends to
   the file it was already writing rather than to a new one. It distinguishes *configurations*, not
   invocations; two runs of the same config still collide, governed by `overwrite_existing`.

   `AbstractSimulationConfig` is the one supertype with fields but no hooks: `build_simulation` is
   generic over the whole `FjordConfig`, because everything dataset-specific already comes from
   the other configs. It is also the only one rooted at a `results_root` rather than a
   `data_root`, being the only config that writes rather than reads.

2. **Dataset adapters** (`src/Datasets.jl`) — `ForcingDataset` and `ResultsDataset` are NumericalEarth
   dataset wrappers for local FjordSim NetCDF files (forcing inputs and simulation outputs).

   **Both are dead code, and broken.** Neither is constructed anywhere in `src/`, `test/` or
   `examples/`; the only call site ever written was already commented out when it was deleted. They
   also no longer work: NumericalEarth 0.6 added a `read_file_coords` step
   (`DataWrangling/metadata_field.jl`) that reads `ds[longitude_name(metadatum)]` and
   `ds[latitude_name(metadatum)]`, defaulting to `"longitude"`/`"latitude"`, whereas the forcing file
   names its axes `Nx`/`Ny` and a results file names them `λ_faa`/`φ_afa`.

   Initial conditions do **not** go through them — see the `Simulations` section. Both prepared files
   are already on the model grid, so reading a state is an `NCDatasets` read and a type conversion,
   and none of NumericalEarth's `Metadatum`/`native_grid`/regrid/inpaint path applies. This is the
   same argument the `Forcing` section makes for not reusing NumericalEarth's dataset path there.

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
   `NumericalEarth.regrid_bathymetry` with `regrid_options(config)` → `smooth_bathymetry_gaps!`
   with `smoothing_options(config)` → `write_bathymetry_file`. The core also owns the smoothing
   kernels and the `center_coordinates`/`expand_domain`/`vertical_faces` domain helpers.

   `smooth_bathymetry_gaps!` runs three stages. The topological cleanup every source gets
   (diagonal-pair fills, isolated sea/land cells), then two stages a source opts into through
   `smoothing_options`, each skipped when its parameter is zero: `fill_shallow_spikes` and
   `limit_bottom_slope`. The order is forced — the topological pass first, because checkerboard
   noise would skew the neighbour medians; despiking before slope limiting, because a spike is
   exactly the one-cell feature slope limiting would smear into its neighbours instead of removing.

   Both exist because `PartialCellBottom` bounds how *thin* a cell may be but says nothing about
   how much depth may change between adjacent *columns*, and that is what destabilizes a regional
   run: a shallow cell beside a deep one carries the same transport in a fraction of the water
   column, so velocity grows there until the time-step wizard's timescale goes NaN and
   `calculate_substeps` throws `InexactError: Int64(NaN)` — the last symptom, not the cause. A
   `minimum_depth` floor alone makes this *worse* in one respect: lifting a sliver to a constant
   depth in much deeper water leaves a shallow spike, which is why the two are configured together.
   `limit_bottom_slope` moves each offending pair symmetrically, so water volume is conserved
   exactly rather than quietly shifted.

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
   `atmosphere_time_steps(config)` gives the hourly `AtmosphereRecord`s → `pad_atmosphere_records`
   extends that axis to the run's `coverage` window →
   `atmosphere_source_grid(config, filepath)` gives the downloaded geometry →
   `atmosphere_target_axes` derives the prepared axes from `x_domain`/`y_domain` grown by
   `config.padding` and sampled at `config.resolution` → `projected_atmosphere_nodes` projects them
   into the source projection in one bulk GDAL call → one bilinear `interpolate_to_target!` per
   variable per step → streaming NetCDF write. `download_atmosphere(config::FjordConfig)` is the
   generic download driver, mirroring `download_forcing`.

   `pad_atmosphere_records` runs *after* `validate_atmosphere_records`, whose hourly-gap warning is
   about download integrity and would otherwise fire on a pad that is deliberately off the hourly
   cadence. See the `Forcing` section for why both prepared files are padded rather than the readers
   made tolerant, and for the one-record-spacing bound.

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
   `forcing_time_steps(config)` completed by `daily_time_steps` to a gap-free daily axis, then
   extended to the run's `coverage` window by `pad_time_steps` →
   `forcing_source_grid(config, filepath)` → land mask from `peripheral_node` (so it matches the
   model, including `PartialCellBottom` and the velocity-face wall convention), with the open
   `relaxation_edge` row restored → source mask filled from the nearest valid cell (`SourceFill`)
   → one trilinear `Oceananigans.Fields.interpolate` per target cell in a `launch!` kernel,
   against the source subset expressed as a `RectilinearGrid` in projected meters
   (`ProjectedSourceGrid`, `source_field_grid`) → relaxation lambdas along `relaxation_edge` →
   streaming NetCDF write. Only the three hooks are dataset-specific; the rest is shared.

   **Time-axis padding.** `pad_time_steps(steps, coverage)` extends the axis to reach both ends of
   the run window, replicating the nearest step — same source records, same blend weight, so the
   written field is identical to its neighbour. `coverage` comes from `coverage_window` at the one
   place that holds both configs, `prepare_forcing(config::FjordConfig)`; `nothing` (a setup naming
   no simulation config) prepares exactly the downloaded range.

   It exists because a run's window rarely lines up with a dataset's records — NorKyst's first daily
   record in Oslofjord is at 12:00, so a run starting at 00:00 has no forcing for its first half-day
   — and because both readers use `Cyclical()` time indexing, which does not fail outside its data
   but wraps to the far end of the year. Padding is what lets `validate_time_coverage` stay a hard
   check instead of the wrap silently filling the shortfall with the following December.

   Two things about it are load-bearing. It runs **after** `daily_time_steps`, because a pad at the
   `SourceRecord` level would put a 12-hour difference into that function's whole-day cadence test
   and collapse every step to a single record, silently disabling gap-filling for the whole year. And
   a pad may reach **at most one record spacing** past the downloaded axis
   (`forcing_record_spacing`): unbounded, it would manufacture the very coverage the validator
   exists to verify — one replicated December would "cover" a window a year past the data,
   `forcing_date_range` would report the invented span, and the run would interpolate twelve months
   between two identical records, which is strictly worse than the wrap it was written to prevent.
   Overshooting by more than a spacing is a missing download and is reported as one.

   Padding changes `Cyclical`'s *inferred* period, which is `(tᴺ - t¹) + Δt` and after a pad equals
   neither the loop period nor anything else meaningful. That is harmless only because the clock-reset
   loop never steps outside `[t¹, tᴺ]`, so the cycling branch is dead code — do not "fix" the period.
   It is also the argument that rules out the monotonic-clock alternative; see `Simulations`.

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

   The two tracer top conditions are not plain `FluxBoundaryCondition`s: they go through
   NumericalEarth's `build_tracer_top_bc`, which wraps the flux field together with a
   `FreshwaterExchange` carrying a freshwater volume flux and its heat content. That is not a
   physics choice here but an interface requirement — NumericalEarth's `net_fluxes` reads the
   exchange back *out of* the boundary condition (`Oceans/Oceans.jl`), so a bare flux condition is a
   `MethodError` the first time the coupled model assembles fluxes. The volume flux field stays
   zero: NumericalEarth turns freshwater into a volume change through a free-surface forcing it adds
   only on a mutable-z grid, and these setups' z faces are static, so the surface T and S fluxes are
   exactly the flux fields the interface writes.

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

    `SimulationConfig` has 25 fields. They split three ways, and the split is forced rather than
    stylistic.
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

    `start_date` is the exception to that list, and it is a field precisely *because* it cannot be
    derived. Every prepared file carries its own first record, and each reader used to zero its own
    time axis there — the Oslofjord forcing starts at 12:00 and its NORA3 atmosphere at 00:00, so
    the two ran twelve hours out of phase with nothing reporting it. One stated instant is the only
    thing they can all agree on, so `build_simulation` passes it to `forcing_from_file` and to both
    atmosphere read hooks as their `reference_date`. Picking a `start_date` earlier than a file's
    first record is therefore an error rather than a shift.

    `validate_time_coverage` enforces that: each prepared file must span
    `[start_date, start_date + stop_time]`. Both readers use `Cyclical()` time indexing, which does
    not fail outside the data it was given — it wraps — so a run that outlasts its forcing would
    otherwise quietly replay the beginning, and one starting early would read the end. The check is
    what makes the wrap unreachable rather than merely unlikely, which is why neither reader's
    `time_indexing` had to change. The atmosphere half goes through the optional
    `atmosphere_date_range` hook, so the module still names no dataset.

    The window checked is **one** `stop_time`, not `loops` of them: every repetition replays the same
    interval, so what the prepared files must span does not grow with the loop count. Multiplying by
    `loops` is the obvious wrong move.

    `build_simulation` returns the instrumented `Simulation` without running it — the REPL and
    debugger entry point — and `run_simulation` calls it and `run!`. Both return `nothing` for a
    setup naming no simulation config. Each prerequisite is checked before anything is read or
    allocated and reported as the command that produces it; the atmosphere is the exception,
    because `MultiYearNORA3(config)` already owns that error, which is what keeps the module free
    of any named dataset.

    ### Initial conditions

    `initial_conditions` holds one of three shapes, resolved by `resolve_initial_conditions` in
    `build_simulation` — where the grid and the forcing path are known — into the `NamedTuple` that
    `coupled_hydrostatic_simulation`'s single `set!(ocean_model; initial_conditions...)` consumes.
    That function is therefore unchanged and learns nothing about the new shapes; the dispatch is
    three methods, not a branch.

    - a `NamedTuple` of constants, functions or fields — already what `set!` wants, so it passes
      through by identity. This is what every setup did before the others existed.
    - `FromForcing(date = nothing)` — a record of the prepared forcing file, `nothing` meaning the
      run's own `start_date`.
    - `FromResults(path, date = nothing)` — a record of a previous run's snapshot file, `nothing`
      meaning its last. A relative `path` resolves against `results_root`, like `output_file`.

    Both files are already on the model grid — `prepare_forcing` regridded onto it, and a snapshot was
    written from fields on it — so this is a read and a type conversion, not an interpolation. Hence
    `src/Datasets.jl` stays unused: none of NumericalEarth's `Metadatum`/`native_grid`/regrid/inpaint
    path applies. Three details are load-bearing. Land is `NaN`/`missing` in the forcing file and is
    zeroed, which is exact rather than approximate because `prepare_forcing` writes `NaN` at precisely
    the cells `peripheral_node` calls dry on this same grid. Arrays are converted to `eltype(grid)`,
    because a `Union{Missing,Float32}` array cannot become a `CuArray` at all. And the date must be on
    the axis exactly — a nearest-record fallback would silently start the run somewhere else.

    A snapshot's `time` is seconds from *its* run's model zero and carries no calendar, so
    `FromResults(path, date)` needs the instant that zero stood for. The snapshot writer records it as
    the `start_date` global attribute (`RESULTS_START_DATE_ATTRIBUTE`), keeping that knowledge with the
    data rather than making it a second config field that could disagree. A file predating the
    attribute can only be read by its last record, and says so.

    What gets set is every tracer the simulation config's `tracers` names plus `u` and `v`,
    intersected with what the source file carries. That list is **never written out** in
    `Simulations`: `state_variables` derives it as
    `(map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)`, the same rule `forcing_from_file` uses to decide
    which forcing terms to build, so adding a biogeochemical tracer to a setup is enough to have it
    read back and the two cannot disagree about what the state is. A tracer the source lacks is left
    at its default rather than being an error, which is what lets one reader serve both file kinds.

    The free surface `η` and any closure-owned tracer the config does not name (CATKE's `e`) keep
    their defaults. So both `FromForcing` and `FromResults` are a **lossy warm start** — the
    turbulence state, the free surface and the Adams-Bashforth tendencies are absent, and the first
    hours are a barotropic adjustment. `pickup` is the exact restart; the two are complementary, not
    alternatives, and this is the most likely thing for a later reader to conflate.

    ### Looping

    `loops` runs the window repeatedly with the ocean state carried over — a spin-up, since one
    forcing year does not equilibrate the deep basins. Each repetition writes its own file
    (`loop_output_path`: the plain run-tagged name when `loops == 1`, `_loopNN` otherwise), so the
    loops can be compared rather than overwriting one another.

    `restart_loop!` sends every clock in the coupled model back to zero and nothing else, then
    re-attaches the writers. `run!` then re-initializes schedules, the timestepper and the initial
    output record, because `Oceananigans.initialize!` does all of that exactly when
    `clock.iteration == 0`. `stop_time` is untouched, and so is `simulation.Δt`, so the wizard keeps
    the step it converged on instead of re-ramping from one second.

    Four things about it are non-obvious and all fail silently otherwise.

    **`NumericalEarth`'s `reset_clock!(::EarthSystemModel)` cannot be used.** Its per-component
    fallback is `reset!(getproperty(component, :clock))` and `components` includes `sea_ice`, so
    `oslofjorden`'s `FreezingLimitedOceanTemperature` — a liquidus and nothing else — makes it throw
    `has no field clock`. FjordSim's own `rewind_clock!` dispatches on `Val(hasfield(…, :clock))`
    instead, which names no component and so keeps working if a setup later names a real sea-ice or
    land model. Do not "simplify" it back; the "Loop restart clocks" testset asserts the upstream
    method still throws.

    **The prescribed atmosphere and radiation need an explicit `update_state!`.** Rewinding a clock
    does not refill a `FieldTimeSeries` window, and `time_step!(::EarthSystemModel)` assembles surface
    fluxes in `maybe_prepare_first_time_step!` *before* it steps the atmosphere — so without this the
    first step of each loop would be forced by last December. This is why NumericalEarth's own
    `reset_clock!` ends the same way.

    **`ocean_sim.initialized` must be cleared by hand.** The ocean is a `Simulation` of its own and
    `run!` clears `initialized` only on the coupled one, so `time_step!(ocean_sim)` would skip
    `initialize!` and the fresh snapshot writer would never have its schedule initialized or its
    `t = 0` record written.

    **The writer is replaced, not renamed, and the old one is closed.** A `TimeInterval` accumulates
    its actuation count, so a writer reused after a clock reset believes it is thousands of records
    ahead and never fires again; and a `NetCDFWriter` holds an open `NCDataset`, so replacing it
    silently leaks a handle and an unflushed tail per loop.

    `clock.last_Δt` comes back as `Inf`, making each loop's first step a forward Euler step. That is
    right rather than merely harmless: `G⁻` refers to a step taken under forcing a year away.

    The alternative — one monotonic clock over `loops * stop_time`, leaning on `Cyclical`'s wrap — was
    rejected. It works numerically, but the wrap period is inferred from each file's own axis and stops
    matching the loop the moment a file is padded, so an explicit period would have to be threaded
    through `forcing_from_file` *and* through the `prescribed_atmosphere`/`prescribed_radiation` hook
    signatures, which are a documented extension contract. Clock-reset leaves both readers untouched.

    ### Checkpointing

    `checkpoint_interval` attaches a `Checkpointer`; `0.0` attaches none at all rather than one that
    never fires. It goes on the **coupled** simulation, not the ocean one, for two reasons that both
    fail otherwise: `prognostic_state` of the coupled model is what a resumable state is, and
    `run!(…; pickup)` looks for its checkpointer in `simulation.output_writers` and requires exactly
    one there.

    Oceananigans 0.110 checkpoints *only* `prognostic_state(simulation)`, so the `Checkpointer`
    docstring's warning that "objects containing functions cannot be serialized" does not apply here:
    `model.forcing` and `model.boundary_conditions` are not in `HydrostaticFreeSurfaceModel`'s
    prognostic state at all, `PrescribedAtmosphere` and `PrescribedRadiation` contribute only their
    clock, and `ComponentInterfaces` and `FreezingLimitedOceanTemperature` contribute `nothing`. So no
    `ForcingFromFile`, no `FieldTimeSeries` backend and no `FreshwaterExchange` is ever written — while
    CATKE's diffusivities and its `e` tracer are. Expect a few hundred MB per checkpoint on the
    Oslofjord grid, which is why `cleanup = true`.

    `pickup` resumes from the newest checkpoint. Two couplings make it work. The checkpoint prefix
    carries the loop index (`checkpoint_prefix`), because the state records the clock but not which
    repetition produced it — without that, `pickup` could not tell a loop-3 checkpoint from a loop-1
    one and would replay the whole spin-up, and loop 2 would overwrite loop 1's files at the same
    iteration number. `resume_loop` reads the index back out of the filename, and
    `build_simulation` attaches its writers for *that* loop rather than unconditionally for the first.
    And `picking_up` suppresses `overwrite_existing`, because the checkpoint restores the writer's
    actuation count: clobbering the file it is meant to continue would leave a schedule thousands of
    records ahead of an empty one. This is also why the run tag must be deterministic — a wall-clock
    tag would put the resumed run in a different file from the one it is appending to.

    `pickup` supersedes `initial_conditions`: `set!` still runs at build time and is then overwritten.

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

    `main` runs the driver inside `tee_output(f, log_path(config))`, which mirrors `stdout` and
    `stderr` to `fjordsim.log` while still printing live to the terminal — a stacktrace through
    `SimulationConfig` spells out every type parameter and is long enough to push the error message
    itself out of the scrollback. `redirect_stdio` only accepts fd-backed streams, so the tee is a
    `Pipe` (needing an explicit `Base.link_pipe!`; an uninitialized `Pipe` throws from `eof`) plus a
    task copying each chunk to both destinations. Both streams share one pipe, so the log interleaves
    them in write order.

    `log_path` puts the transcript under the simulation config's `results_root`, beside the output it
    describes, and creates that directory if absent. Its name carries the run tag —
    `fjordsim_<run_tag>.log` — so runs of different windows do not overwrite each other's transcript,
    matching how `results_path` tags the output. It dispatches on `config.simulation_config` just
    like `simulation_forcing_path`, because `results_root` and `start_date` are *simulation*-config
    fields and a setup need not name one: `drammensfjorden` has no results directory and no run to
    tag, so every step of it falls back to `LOG_FILE` (`fjordsim.log`) in the working directory. That
    fallback is why the name stays in `.gitignore`.

    Every subcommand for a setup that does name a simulation config writes to the *same* tagged log,
    so a later step overwrites an earlier one's transcript. That is pre-existing behaviour plus a tag,
    not a regression, but it means a `run_simulation` transcript does not survive a subsequent
    `add_rivers`.

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
    `@main` for `julia -m FjordSim`. It used to also patch NumericalEarth's
    `compute_bounding_indices` to stop a regional read from running past the file's extent; that is
    fixed upstream as of NumericalEarth 0.6 (`DataWrangling/set_region_data.jl` clamps both the
    bounding-box offset and every index), so the patch is gone and nothing here monkey-patches a
    dependency.

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
| `smoothing_options(config)` → NamedTuple for `smooth_bathymetry_gaps!` | no | `(;)` |

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
| `prescribed_atmosphere(config, architecture; reference_date)` → `PrescribedAtmosphere` | only if the setup is simulated |
| `prescribed_radiation(config, architecture; reference_date)` → `PrescribedRadiation` | only if the setup is simulated |
| `atmosphere_date_range(config)` → `(first, last)` `DateTime`s | no; defaults to `nothing`, skipping the coverage check |

The last three are the read side, consumed by `build_simulation` rather than `prepare_atmosphere`,
and are what keep `Simulations` from naming a dataset. The first two take no float type on
purpose: both NumericalEarth constructors default to `Float32`, and passing
`Oceananigans.defaults.FloatType` would silently promote the atmosphere to `Float64`. A `nothing`
atmosphere config yields `nothing` for all three.

`reference_date` is the instant the returned time axes are zeroed at, and is *not*
`NORA3PrescribedAtmosphere`'s `start_date`: that one selects which records to load, this one
selects where t = 0 sits. `build_simulation` passes the simulation config's `start_date` to both
this hook and `forcing_from_file`, which is the only thing keeping the two in phase — see the
`Simulations` section.

The prepared variable names and units are *not* a hook — they are fixed by the read side in
`ATMOSPHERE_VARIABLES`. A source whose download already normalizes names (as NORA3's does, since
five of its eight variables are derived rather than copied) returns the identity mapping. The
core's `grid_rotation_angle` and `rotate_to_east_north` are available to any adapter whose source
gives wind relative to its own grid axes.

Simulation — `AbstractSimulationConfig`, consumed by `build_simulation`: **no hooks**. It is data
only, read field by field, so an alternative simulation config is a subtype supplying the field
set its docstring lists. It inherits `results_path`, `run_tag` and `coverage_window` — the last two
read `start_date`, so a subtype omitting it inherits none of the three.

The `coverage` keyword `prepare_forcing` and `prepare_atmosphere` now take is a pipeline *argument*,
not a hook: the `FjordConfig` drivers derive it with `coverage_window` and every source inherits the
padding unchanged. No source gains a hook for it.

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
**no defaults at all**: `oslofjorden()` names all 25 fields, because each is a scientific choice
about that fjord and a default would let the next setup silently inherit it. Adding a field to
`SimulationConfig` therefore breaks every setup until each names it, which is the intent.

`start_date` and `stop_time` also decide what the prepare steps write, since both pad their time
axes to that window — so changing either is a data change, not just a run change. See "Changing
`start_date` or `stop_time`" under Commands.

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
- `Pkg.resolve()` only *checks* the manifest against `[compat]`; it will not move a version to
  satisfy a widened bound, and reports `empty intersection between X@old and project compatibility`
  instead. Use `Pkg.update()` after raising a compat bound. (`Pkg.up` is not a name in current Pkg.)
  The `CUDACore` resolve failure this note used to warn about is gone — CUDA now resolves to 6.2.1.
- `Base.@kwdef` on a *parametric* struct does not convert its arguments, unlike on a
  non-parametric one: the generated constructor keeps the declared field types in its signature,
  so `SimulationConfig(stop_time = 3600, ...)` is a `MethodError` where `stop_time = 1hour` is
  fine, and `checkpoint_interval = 0` is one where `0.0` is fine. Every `Oceananigans.Units`
  constant is already a `Float64`, so writing durations with them — `365days`, `1hour`, `3minutes`
  — sidesteps this. `loops` is the one new field this does *not* apply to, being an `Int` already.
  Every other config is non-parametric and so has never hit it.

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land.
- Data files default to `~/FjordSim_data/<fjord>/` and results to `~/FjordSim_results/<fjord>/`,
  the latter from the simulation config's `results_root`.
