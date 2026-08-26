# CLAUDE.md

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
# river config is `standalone`, which writes a forcing file carrying only rivers)
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
julia --project -m FjordSim add_rivers          --config oslofjorden   # re-copies forcing_rivers.nc
julia --project -m FjordSim prepare_boundaries  --config oslofjorden   # rewrites boundaries.nc
julia --project -m FjordSim prepare_atmosphere  --config oslofjorden   # rewrites atmosphere.nc
```

`add_rivers` is not optional here: it `cp`s the forcing file and patches the copy, so
`forcing_rivers.nc` — which is what `simulation_forcing_path` gives the simulation — still carries
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
julia --project -m FjordSim add_rivers          --config oslofjorden   # re-copies forcing_rivers.nc
julia --project -m FjordSim prepare_boundaries  --config oslofjorden   # rewrites boundaries.nc
```

No download step is affected — all four are pure regrids of data already on disk. `prepare_atmosphere`
is *not* in the list: the atmosphere is 2D on its own regular lon/lat grid and knows nothing about the
model's vertical coordinate.

The same list applies to a change in any `bathymetry_config` smoothing knob, since those change the
land mask that `prepare_forcing` and `prepare_boundaries` build their masks from. See "The vertical
grid" under Setups for what the current faces are and why.

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
   `AbstractBoundaryDataConfig`, `AbstractAtmosphereConfig` and `AbstractSimulationConfig`, plus the
   five the simulation config nests — `AbstractCoupledSimulationConfig`, `AbstractBoundaryConditionSetConfig`,
   `AbstractWriterConfig`, `AbstractCallbackConfig` and `AbstractTimeSteppingConfig` — plus
   `AbstractBoundaryConditionConfig`, `AbstractFreeSurfaceConfig` and `AbstractClosureConfig` — one
   level further down again, inside the boundary-condition set and inside the coupled-simulation
   config — and `FjordConfig`, which holds a
   grid, bathymetry, forcing, open-boundary, atmosphere and simulation config (parametrically, so
   every instantiation stays concretely typed). `forcing_config`, `boundary_config`,
   `atmosphere_config` and `simulation_config` all default to
   `nothing`, so a setup opts in by naming one. A river config is the exception: it hangs off the
   forcing config's `rivers` field, where `nothing` means the setup has no rivers, because rivers
   really are interior forcing — they are written into the forcing file itself. Open-boundary data is
   *not*: `boundary_config` is a `FjordConfig` field of its own, since the exterior state is a
   separate hourly file from a separate collection read by a separate pipeline, and either config is
   nameable without the other. It hung off `NorKystConfig.boundaries` until the day a setup needed an
   open boundary with no interior forcing and could not have one. The five
   simulation-level supertypes likewise hang off `SimulationConfig`'s `model`,
   `boundary_conditions`, `writers`, `callbacks` and `time_stepping`. A new grid, bathymetry
   source, forcing dataset, river dataset, open-boundary dataset or atmosphere dataset is
   added by subtyping the matching supertype and overloading methods on it — `FjordConfig` and its
   callers are untouched.

   `AbstractBoundaryDataConfig` sits one word away from `AbstractBoundaryConditionConfig` and
   describes the other half of the same boundary: this one the exterior *values*, that one the
   *scheme* acting on them. Data and physics, configured separately because a setup can change
   either without touching the other.

   The two grid hooks, `domain_grid` and `simulation_grid`, are also declared here rather than in
   `Grids`, because `Bathymetry`, `Atmospheres` and `Forcing` all call `domain_grid` and are all
   included before `Grids`. `model_closure` is declared here for the same reason and one more: its
   fallback is the identity, and an identity fallback is only useful if it is the *first* method
   anyone sees.

   The edge vocabulary lives here too — `LATERAL_EDGES`, `validate_open_edge` and `lateral_edges`,
   beside the `open_edges` accessor that reads a config's own. Every module that reasons about an
   edge is included after this one and none of them can see each other: `Bathymetry` clears land
   along the open ones, `Forcing` restores their velocity face rows and regrids along them, and
   every `Val{edge}` dispatch in `BoundaryConditions` raises an `ArgumentError` naming the tuple.

   The path helpers defined on the supertypes live here too, so a new
   source inherits them without loading the built-in source's module: `bathymetry_path`,
   `forcing_path`, `forcing_directory`, `river_forcing_path`, `boundary_data_path`,
   `boundary_data_directory`, `atmosphere_path`,
   `atmosphere_directory`, `results_path`, and `plot_path` (one method per supertype). Each
   supertype's docstring lists the fields and hook methods a subtype must provide — see "Adding a
   new source" below.

   Two helpers on `AbstractSimulationConfig` name its files and its window. `run_tag` is the run's
   identity as a filename fragment; `results_path(writer, config)` inserts it before the *writer's*
   `output_file`'s extension, resolved against the config's `results_root`, and
   `results_path(writer, config, loop)` appends a zero-padded loop index on top — the filename is the
   writer's so a setup can name several outputs, and a writer naming none (a checkpointer) gets a
   stated `ArgumentError` rather than a missing-field failure. `coverage_window` reads
   `start_date` — which is therefore as much part of the supertype's field set as `results_root` is —
   and is the `(first, last)` calendar interval the run needs its prepared inputs to span, which
   `prepare_forcing` and `prepare_atmosphere` take as their `coverage`.

   The tag is the **wall-clock instant the process started** (`yyyymmddTHHMMSS`), held in
   `Configs.LAUNCH_TAG` and filled by `Configs.__init__` — not by a `const` initialized with `now()`,
   which would freeze the precompilation instant into every later run. So it distinguishes
   *invocations*, not configurations: re-running a config cannot overwrite the earlier run's snapshots
   or log, and `overwrite_existing` only ever bites within one process. Seconds are in the format
   because a crash and an immediate relaunch otherwise collide.

   The consequences of that are all elsewhere and all deliberate. Which simulated window a file covers
   is no longer in its name — it is the snapshot's `start_date` global attribute (see
   `Simulations`), which is why a post-processing step now reads the file rather than predicting the
   path. `run_tag` takes a config it does not read, so a simulation-config subtype can still name its
   runs differently. And a `pickup` cannot name the previous launch's files, which is why checkpoints
   carry no tag — see "Checkpointing".

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

   `progress` takes every reduction on the `Field`, never on `interior(field)`, and that is the whole
   difference between a useful report and a misleading one. A reduction over a `Field` on an
   `ImmersedBoundaryGrid` excludes the immersed periphery for free
   (`Oceananigans.ImmersedBoundaries.NotImmersed`); `interior` hands back a bare array and throws that
   away. Oceananigans writes `zero(eltype)` into every immersed peripheral tracer cell at the top of
   each `update_state!`, and most of a fjord grid is land — so the array form reported the land mask
   as the minimum temperature for the whole of every run, and did it under a label reading
   `extrema(T)` while printing `(max, min)`. It also keeps the reduction on the GPU instead of
   indexing a `CuArray` cell by cell.

4. **Plotting** (`src/Plotting.jl`) — `plot_bathymetry(grid, bottom_height, config)`,
   `plot_forcing(grid, config)`, `plot_boundaries(config)` and `plot_atmosphere(config)`, all
   dispatching on the config *supertypes* so a new source inherits them, writing to
   `plot_path(config)`. `plot_boundaries` reads its edge from the config, because a boundary file's
   variables are named for their side and it has to know which of them to read; it draws two
   panel shapes, a section for a full-depth variable and a line for a surface one. Shares
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

   `smooth_bathymetry_gaps!` runs seven stages in eight passes — `snap_partial_bottom_cells` runs
   twice. The topological cleanup every source gets (diagonal-pair fills, isolated sea/land cells),
   then six stages a source opts into through `smoothing_options`, each skipped when its parameter is
   `false` or zero:
   `clear_open_boundary_land`, `fill_small_islands`, `remove_narrow_passages`,
   `fill_shallow_spikes`, `snap_partial_bottom_cells` and `limit_bottom_slope`.

   The order is forced. The topological pass first, because checkerboard noise would skew the
   neighbour medians. Then `clear_open_boundary_land`, which is the third stage that changes the
   *land mask* and goes before the other two for the same reason `fill_small_islands` goes before
   `remove_narrow_passages`: it only turns land into sea, so it can neither manufacture a narrow
   passage nor strand a cell, while running it first lets the island and passage stages judge the
   band as it will actually be. A land component straddling the band's inner edge loses its in-band
   cells there and is then correctly measured — as whatever remains of it — by `fill_small_islands`.
   It is also the only stage that reads `edges`, and is a no-op for a setup that names none.

   Then the two remaining land-mask stages, before the two that change
   depths: run the other way round, the cells they are about to turn into land or water would already
   have contributed to a neighbour median and to a slope pair. Of those two, `fill_small_islands`
   runs first — flooding an island only widens water and so can never create a passage, while closing
   a passage adds land exactly where a spurious loop closes, which is where these islands live, and
   can bridge one to the mainland where the size test will never see it again. Then despiking before
   slope limiting, because a spike is exactly the one-cell feature slope limiting would smear into
   its neighbours instead of removing.

   `snap_partial_bottom_cells` runs **twice**, once each side of slope limiting, and that is not
   redundancy. It is the only stage that knows where the vertical faces are, and slope limiting moves
   depths without regard for them, so the limiter is guaranteed to put back slivers the snap removed —
   the tighter the limit, the more of them. This used to be a single pass before the limiter, on the
   argument that slope limiting should keep the last word on smoothness, and it cost "a few" slivers:
   178 of 9907 at `max_slope_factor = 0.5`. At 0.25 on the 24-level grid it costs **2628, six percent
   of every bottom cell in the domain**, because a tighter limit means the limiter moves far more
   depth. Snapping again afterwards leaves **none**, and the price is that the second snap breaks the
   limit it was just handed: the steepest slope goes from 0.250 to 0.307, on 16 pairs out of 84 000.
   A slope parameter 23 % over target on a handful of pairs is the better trade, and it is still far
   inside the 0.500 the setup used to accept outright.

   The first snap is not made redundant by the second. It hands the limiter a field already aligned to
   the faces, so the limiter has less to move and the second snap less to correct; and where the
   limiter does not run at all (`max_slope_factor = 0`) it is the only one there is.

   Unlike a passage closure, a flood leaves no stub and so needs no cleanup rounds after it: every
   neighbour of a flooded component is sea by maximality, so no flooded cell borders land.

   ### `clear_open_boundary_land`

   Floods every land cell within `open_boundary_land_cells` rows of an **open** lateral boundary,
   giving each the mean depth of its wet neighbours.

   A coastline that runs into an open boundary is the worst place in the domain for one to be. The
   boundary condition there is already reconciling a radiation scheme against a prescribed exterior
   state within a cell or two, and a headland poking through the boundary row splits the prescribed
   inflow around an obstacle the exterior dataset never saw. Clearing a band leaves the open edge the
   clean, uninterrupted channel every open-boundary scheme is derived for.

   Three things about it are load-bearing.

   **The edges come from `open_edges`, not from a field of the bathymetry config.** They reach
   `smooth_bathymetry_gaps!` as its `edges` keyword, threaded by
   `prepare_bathymetry(target_grid, config; edges)` and supplied by
   `prepare_bathymetry(config::FjordConfig)` from `open_edges(config.boundary_config)` — the same
   shape, and the same argument, as `prepare_forcing`'s `edges` keyword. Which edges are open is a
   property of the domain and its exterior data, not of a bathymetry source, and a second statement
   of them could only disagree. A setup naming no boundary config prepares its bathymetry with every
   lateral boundary a wall, which is the right reading when there is no exterior state to admit.

   **Depths are assigned by repeated relaxation, not in one pass.** Each round fills the cells that
   have at least one wet neighbour *now*, so a cell in the middle of a cleared headland inherits from
   the coast through the cells between it and the band ramps rather than stepping. Each round is
   computed in full before it is written, so the result does not depend on the order the band is
   walked. A cell no round can reach stays land, which needs a land component in the band touching
   nothing wet at all.

   **The price is exterior data.** A newly wet cell has no profile of its own in the prepared
   boundary file, and `fill_boundary_gaps!` fills it from the nearest wet cell along the boundary —
   the same treatment a boundary column already gets wherever the model and the source disagree about
   the coastline. So the trade is a few laterally interpolated columns against a headland in the
   boundary row, and the width should stay a few cells rather than a few tens. On `oslofjorden` the
   southern band is nearly clear already: rows 1–5 hold seven land cells, row 1 is entirely water, and
   `open_boundary_land_cells = 5` moves those seven and changes nothing in the boundary row itself.

   `open_boundary_band` is four `Val{edge}` methods returning the band's index ranges, with a
   catch-all that raises — the same shape as `open_boundary_water!` in `Forcing`, and never a
   four-branch `if`.

   `fill_shallow_spikes` and `limit_bottom_slope` exist because `PartialCellBottom` bounds how *thin*
   a cell may be but says nothing about
   how much depth may change between adjacent *columns*, and that is what destabilizes a regional
   run: a shallow cell beside a deep one carries the same transport in a fraction of the water
   column, so velocity grows there until the time-step wizard's timescale goes NaN and
   `calculate_substeps` throws `InexactError: Int64(NaN)` — the last symptom, not the cause. A
   `minimum_depth` floor alone makes this *worse* in one respect: lifting a sliver to a constant
   depth in much deeper water leaves a shallow spike, which is why the two are configured together.
   `limit_bottom_slope` moves each offending pair symmetrically, so water volume is conserved
   exactly rather than quietly shifted.

   ### `remove_narrow_passages`

   `remove_narrow_passages` bounds something neither of those two can see: channel **width**. Both of
   them bound depth *contrast*, so a one-cell-wide sea passage whose neighbours are equally shallow
   passes every check they make — and it is the width that is wrong.

   A one-cell-wide passage is a sea cell that is sea on both sides along one axis and land on both
   sides along the other. Regridding leaves one wherever a strait too narrow to resolve cuts through a
   peninsula. The two basins such a cell joins are usually *already* connected elsewhere, so it closes
   a loop, and the barotropic head difference around that loop is forced through a cross-section one
   cell wide and a few metres deep. With only a quadratic bottom drag and a biharmonic viscosity to
   resist it, velocity there grows without bound.

   On `oslofjorden` there were 66 such passages, 23 of them at exactly the 2 m `minimum_depth` floor.
   One — a 2.35 m canal through a peninsula at 10.43°E, 59.09°N (i = 67–68, j = 49) — carried a
   coherent, depth-independent 47 m s⁻¹ jet that held the domain maximum in 36 of 39 snapshot records
   and collapsed the adaptive time step from its 3-minute cap to 0.3 s. The same defect blew up the
   pre-open-boundary runs too, after 3.2 days rather than 4 hours: a genuinely open boundary admitting
   a tide only excites it sooner.

   Three things about the stage are load-bearing.

   **A passage that is the sole link to a basin is kept**, since closing it would delete that water
   from the domain. On `oslofjorden` 12 of the 78 one-cell passages are of that kind.

   **Candidates are tested one at a time against the partially closed field**, not all at once against
   the input. Two passages that are each redundant *while the other is open* would together sever a
   basin if closed as a batch; sequential testing makes the stage unable to disconnect anything. It is
   therefore order-dependent, but deterministically so — row-major.

   **One pass, not iterated to convergence.** Closing a passage can leave a neighbouring cell
   one-cell-wide in turn, and iterating would erode a genuine narrow arm cell by cell. The
   isolated-cell cleanup loop does run again afterwards, to take the dead-end stubs a closure leaves
   behind; it cannot re-open a closed passage, because a passage cell has two land neighbours by
   definition and they stay land, so `fill_isolated_land_cells` — which needs all four wet — can never
   fire on it. That fill is also the likely *origin* of some of these canals: it floods a one-cell-thick
   isthmus, which is why the test fixtures make their peninsulas two cells thick.

   `wet_component` is the flood fill it tests with, hand-rolled rather than
   `ImageMorphology.label_components` (which `NumericalEarth.remove_minor_basins!` uses): that package
   reaches FjordSim only transitively, and one breadth-first walk is all this needs.

   ### `fill_small_islands`

   The dual of the stage above, and the one the *land* mask needs. `remove_narrow_passages` bounds
   the width of a water channel; this bounds the size of a land obstruction, which closes the same
   kind of loop with the phases swapped. Flow splits around a rock a few cells across, so the head
   difference between its two ends is worked out around a closed loop whose arms are each a handful
   of cells wide. Nothing about the flow such an island obstructs is resolved by one cell of it, so
   keeping it buys no fidelity.

   It floods every 4-connected patch of land of at most `max_island_cells` cells that does not touch
   the domain edge, setting the whole patch to the mean depth of the sea orthogonally adjacent to it.

   On `oslofjorden` this was the domain's velocity maximum: 2.67 m s⁻¹, near uniform over the water
   column (2.28 m s⁻¹ at k = 18 falling to 1.35 at k = 14), circulating around a **three-cell,
   one-cell-wide island** at i = 83, j = 104–106 — 10.485°E, 59.184°N, near Horten. A six-cell
   one-cell-wide ridge sits beside it at i = 84–86, j = 99–102 and forms the same pinch. Of that
   bathymetry's 85 land components, 46 are interior clusters of 2 to 6 cells, 153 cells in all, and
   60 are of 10 cells or fewer. `max_island_cells = 6` takes all 46.

   Four things about the stage are load-bearing.

   **`fill_isolated_land_cells` is its `max_cells = 1` case**, and cannot be stretched to cover the
   rest: it needs a land cell's four neighbours all wet, and every cell of a three-cell ridge has a
   land neighbour, so it fires on none of them. Neither `fill_shallow_spikes` nor `limit_bottom_slope`
   sees them either, for the same reason neither sees a narrow passage — both bound depth *contrast*,
   and an island agrees with its surroundings. It is its *presence* that is wrong.

   **A component touching the domain edge is kept**, whatever its size, because land continuing
   outside the domain may have only a few cells inside it and flooding those would open the domain
   into water that is not there. This is the component-wise form of the `2:Nx-1, 2:Ny-1` bound the
   single-cell stages use.

   **The patch gets one flat depth.** A cell in the middle of a 2x3 patch has no sea neighbour of its
   own, and a flat floor over a patch this small is what slope limiting would produce anyway. A sea
   cell touching the patch on two sides contributes twice, which weights the mean by contact length.
   For a single cell the value is exactly what `fill_isolated_land_cells` computes.

   **No cleanup pass follows, and none is needed.** Every neighbour of a flooded component is sea by
   maximality — a land neighbour would be part of the component — so no flooded cell borders land and
   neither isolated-cell stage has anything new to fire on. Flooding cannot manufacture a narrow
   passage either: it only turns land into sea, so the land-neighbour count of any sea cell can only
   fall, and the passage predicate needs it to rise. Measured on the `oslofjorden` field, flooding at
   `max_island_cells = 6` leaves the one-cell passage count unchanged at 18.

   The order matters concretely here: applied *after* the two depth stages the flood leaves
   `max r = 0.70` and two new spikes, where in its actual position `fill_shallow_spikes` and
   `limit_bottom_slope` clean up behind it.

   `land_component!` is the walk it uses, the land counterpart of `wet_component` and deliberately a
   second function rather than a predicate parameter on the first. That one answers a connectivity
   question about a single candidate, takes a `blocked` cell for it and returns a full-domain mask;
   this one is swept over the whole field, so it needs a `visited` buffer shared across components —
   what makes the sweep walk every land cell exactly once — and the membership list itself, which is
   what the caller measures against the threshold. Generalising `wet_component` to take an external
   buffer *and* a predicate *and* return a vector would leave nothing of it.

   ### `snap_partial_bottom_cells`

   The one stage that knows about the *vertical* grid, and the fix for the worst defect the open
   boundary had.

   `PartialCellBottom` will not make a bottom cell thinner than `minimum_fractional_cell_height`
   times its layer — its kernel takes `min(z⁺ - ϵ Δz, zb)`, pushing the bottom *down* until the cell
   is exactly that thick. So a sounding lying just below a face does not produce a thin cell; it
   produces a cell of exactly `ϵ Δz` whose floor is somewhere the sounding never was. This stage
   moves the bottom the other way, up to the face, ending the column one layer higher and leaving
   every surviving bottom cell a full one. `minimum_cell_fraction` must be the *same* number the grid
   hands `PartialCellBottom` — 0.2, Oceananigans' default, for the grids `Grids.jl` builds.

   Such a cell is dangerous out of proportion to its size, and for **tracers** rather than for
   momentum, which is what distinguishes it from the `minimum_depth` sliver problem. It holds a
   fraction of the water its neighbours do, so any flux into it moves its concentration a long way;
   and because it reaches into a layer its neighbours may not reach at all, it can be nearly cut off
   horizontally too.

   On `oslofjorden`, after the islands were flooded, this was the domain's remaining pathology and it
   sat on the **open southern boundary**. Soundings of 51.3 m and 54.5 m at i = 38–39, j = 1–2 against
   a layer spanning −75 to −50 m, floored to fractional heights of 0.052 and 0.179, with one and two
   lateral neighbours at their own level because the columns beside them (44.1 m, 48.5 m) stop a whole
   layer higher. Salinity there climbed 33 → 65 psu and temperature 5 → 28 °C over four and a half
   days, while the prepared boundary file was asking for 33.2 psu and 8.1 °C and the cell was
   inflowing 70% of the time — so the nudging was active the whole way and simply could not keep up.
   Domain-wide the tracer overshoots were six times over-represented in bottom cells with at most one
   lateral neighbour (25% of offenders against 4% of columns).

   Three things about it are load-bearing.

   **Refining `z_faces` is not the fix for *this*, though it is the obvious guess.** A finer vertical
   grid gives the seabed more faces to cross, so laterally isolated bottom cells become *more* common,
   not less — measured, 3.9% of columns on the 18-level grid against 5.5% on the 24-level one that
   replaced it. It is the sliver that has to go, not the layer that has to shrink, and that is what
   this stage does.

   The vertical grid was separately and genuinely wrong, and has since been fixed — see "The vertical
   grid" under Setups. The two are complementary: this stage removes a bottom cell that is too *thin*,
   the regrading removes one that is too *thick*. Neither substitutes for the other, and the finer
   grid makes this stage matter slightly more rather than less.

   **It snaps up, not down.** Down is what the model already does silently, by up to 3.7 m at the
   cells above, and it keeps the sliver. Up ends the column at a face, which also rejoins it to the
   level its neighbours are already on: the three runaway columns all land on exactly 50.0 m, where
   the 44.1 m and 48.5 m columns beside them already were. The price is water volume — 0.95% of the
   domain, median 1.85 m per affected column, 10 m at worst in the 50 m layers where 1% of columns
   live.

   **Two guards keep it total.** A column is left alone when the face above is the surface, or when
   snapping to it would breach `minimum_depth`, so the stage can never dry a cell out or undercut the
   floor; such a column keeps its sliver and is floored as before. On `oslofjorden` neither fires.

   **It snaps in the file's precision, not the pipeline's**, and that is not a detail. Smoothing runs
   in the grid's float type — `Float64` — where the snapped value lands on the face exactly; the
   narrowing to `Float32` happens later, in `write_bathymetry_file`, on a value the stage has stopped
   looking at. Wherever a face is not representable in `Float32`, that narrowing puts the column back
   a hair *below* it, and the layer beneath becomes its bottom cell at about 5 × 10⁻⁸ of full
   thickness — the pipeline writing out, in the worst possible form, the very sliver it removed. Four
   of `oslofjorden`'s 24 faces are inexact (−10.8, −7.9, −3.7, −2.2) and this hit 1265 columns, 2.9 %
   of the domain, every one of them at those four. `snap_to_face` rounds *away* from the water so the
   value survives the round trip, using `BATHYMETRY_ELTYPE` — the one statement of what the file
   stores, shared with `write_bathymetry_file`. The old 18-level faces were all exactly representable,
   which is why this never bit before.

   With both the second pass and this fix in place, `PartialCellBottom`'s clamp does not fire anywhere
   in the finished `oslofjorden` field: zero slivers, against 292 on the 18-level grid.

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

   The read side loads river relaxation from NetCDF via `forcing_from_file`. The
   `ForcingFromFile` struct carries two `FieldTimeSeries` (values + lambdas) and dispatches to
   flux, advection, or relaxation terms based on the sign of lambda: `λ > 1` → x-flux,
   `λ < -1` → y-flux, `|λ| < 1` → relaxation. Uses a custom `NetCDFBackend` keeping 2 time
   indices in memory. `forcing_from_file` takes either a `filepath` keyword or an
   `AbstractForcingConfig` positionally, resolved by `forcing_path`.

   `build_simulation` reaches this through one more layer of dispatch, `simulation_forcing(config,
   grid, filepath, tracers, reference_date)`, rather than calling `forcing_from_file` directly: the
   hook is what a forcing source overloads if its prepared files are not this NetCDF layout, and it
   takes `filepath` rather than resolving it because `build_simulation` may pass the
   rivers-augmented copy (`simulation_forcing_path`, in `Simulations`) instead of `forcing_path`
   itself. Every built-in source shares one default, since the layout is a contract fixed by the
   read side, not a per-source choice.

   `prepare_forcing(target_grid, config::AbstractForcingConfig)` regrids the downloaded source
   files onto the simulation grid: `forcing_variable_names(config)` picks the variables →
   `forcing_time_steps(config)` completed by `daily_time_steps` to a gap-free daily axis, then
   extended to the run's `coverage` window by `pad_time_steps` →
   `forcing_source_grid(config, filepath)` → land mask from `peripheral_node` (so it matches the
   model, including `PartialCellBottom` and the velocity-face wall convention), with the open edge's
   velocity face row restored → source mask filled from the nearest valid cell (`SourceFill`)
   → one trilinear `Oceananigans.Fields.interpolate` per target cell in a `launch!` kernel,
   against the source subset expressed as a `RectilinearGrid` in projected meters
   (`ProjectedSourceGrid`, `source_field_grid`) → **zero** lambdas → streaming NetCDF write. Only the
   three hooks are dataset-specific; the rest is shared.

   **There is no interior relaxation band.** `prepare_forcing` writes `λ = 0` everywhere, and
   `add_rivers` writes the only nonzero lambdas the file ever carries. The band used to taper
   `1/relaxation_timescale` inward from `relaxation_edge` across `relaxation_cells`; it is gone
   because the open lateral boundary now nudges towards hourly exterior data *at* the boundary
   (`NormalRadiation`, `GravityWaveRadiation` — see `BoundaryConditions`), and a band relaxing the
   same variables a few cells inside would fight it. `relaxation_lambda` and `edge_distance` were
   deleted with it, and `NorKystConfig` lost `relaxation_cells` and `relaxation_timescale`. What
   `forcing.nc`'s *values* are still read for is initial conditions (`FromForcing`) and rivers.

   This is a physics change, not a refactor: nothing damps the interior near the boundary any more.
   If a run turns out to need it, the fix is a viscosity sponge through the closure, not a
   reinstated tracer band fighting the boundary nudging.

   `relaxation_edge` was renamed **`open_edge`**, since plural **`open_edges`**, in the same move — a field named for a relaxation
   that no longer happens there was exactly the defect this change set out to fix — and
   `validate_relaxation_edge` became `validate_open_edge`. The field itself has since left this config
   entirely: which edges are open is a property of the domain and its exterior data, not of a forcing
   dataset, so they live on the `AbstractBoundaryDataConfig` and reach `prepare_forcing` as its
   `edges` keyword. A setup naming no boundary config passes `nothing`, `water_mask` restores no face
   row, and every lateral boundary is a closed wall — which is the right reading, since with no
   exterior state there is nothing to read at that row. A domain in the open ocean names all four and
   gets every face row back; `lateral_edges` normalizes one `Symbol`, a collection or `nothing` to one
   shape, so restoring none, one or four rows is the same loop. `LATERAL_EDGES` is now stated once, here,
   and imported by `BoundaryConditions`; the two modules previously held identical copies under
   different names with identical error strings. `open_boundary_water!` is now four `Val{edge}`
   methods rather than a four-branch `if`.

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
   years, output names), `forcing_monthly_filename`, the three hook methods, and
   the NorKyst download: list the THREDDS catalog → open the month's OPeNDAP datasets → subset to
   the target grid's lon/lat box (`subset_ranges`, `NorKystSubset`) → write one combined monthly
   NetCDF.

   Two library functions deliberately *not* used, documented in `prepare_forcing` and
   `SourceFill`: NumericalEarth's dataset path (`native_grid`, `set!(field, metadata)`,
   `DatasetRestoring`) always builds a `LatitudeLongitudeGrid`, but the NorKyst grid is rotated
   ~59° from east here; and `inpaint_mask!` cannot fill a fully masked depth level, so on this
   regional subset it either never terminates or silently writes zeros.

   `rivers.jl` (generic core) and `of800_rivers.jl` (dataset adapter) are included into the same
   `Forcing` module and hold the rivers step, which runs *after* `prepare_forcing` unless the river
   config is `standalone`.
   `add_rivers(target_grid, config::AbstractForcingConfig)` dispatches on `config.rivers` — a
   `nothing` river config is a no-op, so a setup opts in by naming one. The pipeline:
   `river_forcing_times` → `river_locations(rivers)` → snap each outlet to a grid cell
   (`river_cells`) → river values for that time axis (`river_series`) → get the base file → patch its
   surface level (`write_rivers`).

   ### `standalone`

   The base file is where the river config's `standalone` decides. `false` copies the prepared forcing
   to `river_forcing_path(rivers)` and patches the copy, taking the time axis from it: the original is
   never modified, so the step is re-runnable without redoing `prepare_forcing`, and a missing prepared
   forcing is an error naming the step that produces it. `true` writes the file from scratch
   (`create_river_forcing_file`), carrying **only** rivers — so a setup can have rivers with no
   interior forcing at all, no NorKyst download and no regrid.

   It is opt-in rather than inferred from whether `forcing.nc` happens to exist, because the two
   failure modes are not symmetric: a setup that meant to prepare forcing and forgot gets an error
   naming the step, where an inferred fallback would have silently run a whole year forcing-free. A
   `standalone` config never copies either, even when a prepared forcing is there — copying it would
   pull interior forcing into a run that asked for none.

   Three things about the standalone file are load-bearing. Its **time axis** comes from the run window
   (`coverage`, hence a `simulation_config`) rather than from a prepared file, daily, anchored at
   `start_date` and extended past `start_date + stop_time` so `validate_time_coverage` accepts it;
   daily because `river_series` matches records by *calendar date*, so a finer axis would only ask for
   the same record repeatedly. Its **variables** are whatever `river_series` returns, since there is
   nothing else in the file. And only its **surface level** is written — `NaN32` values and zero
   lambdas, the level `write_rivers` patches — with the levels below left unwritten and read back as
   the file's `_FillValue`, which is inert twice over: a non-finite value becomes the `-999.0` sentinel
   every branch of `ForcingFromFile` gates on with `value > -990`, and a non-finite lambda fails each of
   `λ > 1`, `λ < -1` and `-1 < λ < 1`. That is the same treatment dry cells already get in a prepared
   file, and it keeps the write to a fraction of a full 4D one.

   The file carries a `rivers_only` global attribute so it states what it is, and
   `Simulations.forcing_state` refuses it: `FromForcing` initial conditions would otherwise start a run
   from `T = S = 0` everywhere, silently, since `finite_slab` zeroes every non-finite cell.

   Rivers enter as **relaxation, not as a mass flux**: each river cell gets its value and
   `λ = 1 / relaxation_timescale` (1 hour by default) at the surface level for every time step,
   which lands in the existing `|λ| < 1` regime — no new forcing term or λ convention. Outlets are
   located by
   independent nearest-node lookups in longitude and latitude, then moved to the nearest
   *coastal* water cell — water with at least one land neighbour, so a river cannot be injected
   into open water — of at least `minimum_levels` wet levels, so it cannot be injected into a column
   too shallow to carry it either. An outlet outside the grid, or with no such cell within
   `search_radius`, is dropped with a warning rather than written into land.

   That depth rule exists because the river is written into the **surface level alone**, so the
   column *beneath* it is what carries the exchange the freshening drives — and one cell cannot. The
   fresh surface cell sets up an estuarine circulation, and the salty inflow at depth concentrates in
   the single cell below instead of spreading through a column. On `oslofjorden` four outlets had
   snapped onto the 2 m `minimum_depth` floor, a column of two 1 m cells, and held 32 to 64 psu under
   a surface cell relaxed to 0, while all fifteen outlets with four levels or more stayed between 29
   and 35. Column salt was conserved throughout — nothing was created, the column simply could not
   resolve the redistribution — and it oscillated rather than diverging, which is what distinguishes
   it from the open-boundary runaway `snap_partial_bottom_cells` fixes. All four had a four-to-six
   level coastal cell within three cells, so `minimum_levels = 4` moves them ~500 m and no further.

   `minimum_levels` is applied by masking the too-shallow columns out of `coastal_water_mask` itself
   rather than by a separate depth test, which is what keeps `is_coastal_cell` and
   `nearest_coastal_cell` unchanged: a column a river cannot enter simply counts as shore for this
   purpose, so the nearest acceptable cell is by construction both coastal and deep enough. The water mask comes from the same
   `water_mask` that `prepare_forcing` uses, so "water" means the same thing in both.
   `write_rivers` reads, patches and writes back whole surface slabs because the file is
   chunked one horizontal slab per `(level, time)`.

   `boundaries.jl` (generic core) and `norkyst_boundaries.jl` (dataset adapter) are included into
   the same `Forcing` module and hold the **open-boundary data** step, which is what makes the
   lateral boundary open at all. `prepare_boundaries(target_grid, config)` dispatches on
   `FjordConfig`'s own `boundary_config` — a `nothing` boundary config is a no-op, so a setup opts in
   by naming one, the way it opts into rivers, but on a `FjordConfig` field rather than inside the
   forcing config. It says nothing about forcing in either direction: a setup can have an open
   boundary with no interior forcing, or interior forcing with closed walls. The edges are
   `open_edges(config)`, read from the config that states them rather than passed beside it, where the
   two could disagree, and the whole per-variable build below runs once per edge into one file. The
   pipeline mirrors `prepare_forcing`:
   `boundary_variable_names` → `boundary_time_steps` completed by `hourly_time_steps` and padded by
   `pad_time_steps` → `boundary_source_grid` → the edge slice of `water_mask` →
   `SourceFill` → the same trilinear kernel → `write_boundaries_file`.

   It exists as a **separate file from `forcing.nc` because it is hourly**. A Flather boundary
   compares the model's own `η` against an exterior one, and NorKyst's daily-average collection has
   the tide averaged out of it — elevation with no tide is not worth prescribing. Only the boundary
   row needs that cadence, so `boundaries.nc` is a few hundred MB where an hourly interior forcing
   would be hundreds of GB. Putting an hourly variable *inside* `forcing.nc` was the alternative and
   was rejected: it would mean a second time axis inside a contract fixed by the read side, touching
   the streaming write, the padding, `forcing_date_range`, the rivers copy loop, plotting and the
   test fixture — the same second-axis problem the atmosphere pipeline was built to avoid.

   Every variable is named for its side — `south_T`, `west_T`, … — which is how several open edges live
   in one file: `prepare_boundaries` loops over `open_edges(config)` and writes every side's variables
   into the same file, on one time axis, with no format change. A domain in the open ocean naming all
   four therefore has one boundary file with four sides in it, not four files. All six spatial dimensions are defined in every boundary file, each variable
   using the two it needs. Layout, fixed by the read side: a full-depth variable is
   `(along, Nz, time)`, a surface one `(along, time)`; the along-boundary dimension is staggered with
   the variable (`south_u` and `south_ubar` on `Nx_faces`); `Float32`, `NaN` on land, CF-encoded
   `time`, and an `open_edge` global attribute — a comma-separated list — so the file states its own
   sides. **No `_lambda`
   twins** — this is boundary data, and the nudging timescales belong to the boundary *condition*.

   Three things about it are load-bearing.

   **It reuses `prepare_forcing`'s machinery on a one-cell-thick slab.** `PreparedVariable`'s mask is
   `(Nx, 1, Nz)` for a south edge, which every three-dimensional structure in the forcing core
   already accepts — `SourceFill`, the 3-index `_interpolate_to_target!` kernel, `launch!` over
   `size(output)`. `PreparedVariable` became parametric in its dimension-name count, which is the one
   thing the two files disagree about; nothing else needed widening.

   **The three surface variables need their own single-level source grid.** `solve_vertical_faces`
   cannot build one — for a single level it leaves the lower bound at `-Inf` and puts the first face
   there — so `surface_source_field_grid` centres its one cell on `SURFACE_TARGET_DEPTH`, which is
   also the vertical node the target slab gets, making the interpolation land exactly on the cell
   centre. `source_validity` and `source_slab` became `ndims`-aware and reshape a 2D source plane to
   one level (`as_source_slab`), which keeps everything downstream three-dimensional.

   **`fill_boundary_gaps!` fills dry boundary cells on read.** The file marks them `NaN` like every
   FjordSim prepared file, but a boundary series is read as an exterior *value* by schemes that do
   arithmetic on it, and only some of the nodes they write are ones the grid calls closed. The
   tangential velocity forces this: on a south edge it is `u`, whose faces at `i = 1` and
   `i = Nx + 1` are peripheral so the file has holes there, yet `immersed_peripheral_node` — what
   `NormalRadiation` zeroes on — is false at a wet domain-edge face. So the hole would be nudged into
   the halo at both bottom corners. The fill takes the nearest finite value along the boundary, which
   is what `SourceFill` and `daily_time_steps` already choose, and zero for a wholly dry line.

   The read side is `boundary_series(config, grid, reference_date)` → a `NamedTuple` keyed by **edge**,
   whose entries are `NamedTuple`s of **reduced `FieldTimeSeries`** keyed by the *bare* names
   (`(; south = (; T, eta, …), west = …)`). Keyed by edge because a domain may be open on several, and
   by bare name inside because the side is already the outer key — so no boundary condition spells a
   side into a variable name, and the group a scheme consumes is the same shape whichever edge it came
   from. `validate_boundary_file_edges` checks the file states every side the setup opens, so a setup
   that has since opened another edge is told to re-run the step rather than failing on a missing
   variable. The reduced locations
   are exactly the flavours `Oceananigans.BoundaryConditions.getbc` accepts as a boundary condition —
   `XZFTS` (`{LX, Nothing, LZ}`) for a south or north edge, `YZFTS` for west or east — which is why
   `BoundaryConditions` can pass a series straight through with no discrete-form wrapper. `Cyclical()`
   indexing as everywhere else, and `validate_time_coverage` keeps the wrap unreachable.
   `BoundaryBackend` is a second streaming backend rather than a widening of `NetCDFBackend`: the two
   read different ranks, and putting a rank branch in `load_from_netcdf` would put it in the interior
   forcing's hot read path for no benefit there.

   `norkyst_boundaries.jl` holds `NorKystBoundariesConfig <: AbstractBoundaryDataConfig` and reads a
   *different* NorKyst collection from `norkyst.jl`'s: `fou-hi/norkyst800m-1h`, whose
   `..._his.an.*` files carry 24 hourly instants each rather than one daily average, and which
   publish `zeta`, `ubar` and `vbar` alongside the four 3D variables. Two things differ from the
   interior download, and it shares every helper below the driver
   (`define_output_file`, `write_parameter_chunk!`, `write_time_dependent_coordinates!`): it subsets
   to a `margin`-wide lon/lat **band** along the edge rather than the whole domain box, without which
   an hourly download would be ~30 GB instead of ~3, and its month loop reads 24 records per file.
   NorKyst's grid is rotated ~59° from east here, so a thin geographic band is not thin in either `X`
   or `Y` — its bounding index box is still about a tenth of the whole-domain one, which is the point.

   ### `ubar`/`vbar` are grid-relative and must be rotated

   The two velocity pairs this collection publishes are **not in the same frame**, and that asymmetry
   is load-bearing. `u_eastward`/`v_northward` are diagnostics met.no has already derotated to
   geographic axes. `ubar`/`vbar` are what ROMS itself writes — depth-averaged velocity along the
   *model's own curvilinear x and y axes* — and the collection publishes no derotated twin for them.

   Since NorKyst's grid is rotated ~59° from east here, `vbar` in the Oslofjord region is very nearly
   the **eastward** barotropic velocity, not the northward one. Used unrotated it reaches
   `GravityWaveRadiation` as the exterior transport normal to a south edge, so the Flather condition
   nudges the barotropic mode towards a signal of roughly the right magnitude and essentially no
   correlation with the true normal flow. That is a persistent spurious transport across the open
   boundary, not a small error — measured on a year of prepared data, `cor(vbar, depth-mean v)` was
   **-0.034** while `cor(vbar, depth-mean u)` was **-0.916**, and the best-fit rotation carrying
   `(ubar, vbar)` onto the depth average of the geographic pair was **+59.05°** at gain 1.028, i.e. a
   pure rotation. In a run it showed up as a band of ~0.5 m/s mean flow 4–7 km inside the southern
   boundary, some 50x what the domain's tidal prism can account for.

   `boundary_source_slab(config::NorKystBoundariesConfig, …)` fixes this by rotating the pair with
   `rotate_to_east_north` and NorKyst's own `angle` variable — 58.86–59.92° over this subset, which
   `define_output_file` already copies into every downloaded monthly file, pre-subset to the band. So
   the fix costs no re-download.

   Three things about it are load-bearing. It runs on the **native source grid before interpolation**,
   the same order `nora3_source.jl` rotates NORA3's winds in, because afterwards the two components sit
   on different edge rows (`Nx_faces` versus `Nx`) under different masks and can no longer be combined.
   It needs a **hook** rather than a line in the pipeline because `prepare_boundaries` prepares one
   variable at a time and a rotation needs both — `boundary_source_slab` is the one step of the writer's
   per-step loop a source can override, and both slabs are reachable there from the same `reader` at the
   same `step`. And `angle` is read from whichever file the reader has open rather than cached, which is
   safe only because it is grid geometry and identical in every monthly file of one subset.

   The 3D pair needs none of this, and neither does `zeta`, being a scalar — so `norkyst.jl`, which
   reads only `u_eastward`/`v_northward`, has no equivalent problem.

   Worth knowing: NorKyst's own global attributes name its lateral boundary conditions
   `zeta: Che`, `ubar/vbar: Shc`, `u/v/temp/salt: RadNud` — Chapman, Flather-family and
   Orlanski-plus-nudging. That is exactly the scheme set FjordSim now applies at its own edge, which
   is not a coincidence so much as a reason to trust the choice.

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

8. **Boundary conditions** (`src/BoundaryConditions.jl`) — the `boundary_condition_sides` hook on
   `AbstractBoundaryConditionConfig`, the `field_boundary_conditions` hook on
   `AbstractBoundaryConditionSetConfig`, their built-in configs, and the functions they wrap.

   **Two levels, two hooks.** A piece says *what* one part of the boundary is; a set says *how* the
   pieces combine and *what shape* the result takes. `MergedBoundaryConditions(pieces...)` is the
   built-in set: it merges every piece's `boundary_condition_sides` with `recursive_merge` (last
   wins, so argument order is precedence; naming no pieces yields `(;)`) and materializes the result
   into `FieldBoundaryConditions`. That second level used to be
   `field_boundary_conditions(configs::Tuple, …)` — dispatch on a bare `Tuple`, which no subtype can
   specialize, so the merge strategy and the `FieldBoundaryConditions` result shape were fixed for
   every model FjordSim could assemble.

   Three built-in pieces, so a setup composes its boundaries out of parts it can drop or swap one at
   a time. `AirSeaFluxes()` creates the wind/heat/salt flux fields at the top;
   `QuadraticBottomDrag(coefficient)` creates the bottom drag. Those two were one config,
   `TopBottomFluxes`, fused behind `top_bottom_boundary_conditions` — but they share no computation,
   and swapping a drag law should not mean restating the surface. `top_bottom_boundary_conditions`
   remains as the one-line merge of both, which is the shape the tests assert against.

   ### `QuadraticBottomDrag` needs both halves

   The drag contributes **two** conditions per velocity component: a `bottom` one and an `immersed`
   one, an `ImmersedBoundaryCondition` wrapping `u_immersed_bottom_drag`. The second is the one that
   does anything on a fjord, and the config produced no drag at all until it was added.

   A `bottom` condition acts at the *underlying* grid's floor, and `u_quadratic_bottom_drag` hardcodes
   that level's index — it reads `Φ.u[i, j, 1]`. These setups give `z_faces` a deepest face well below
   the deepest sounding so no column is clipped, which means level `k = 1` is immersed across the whole
   domain and the flux is discarded: on `oslofjorden` the deepest of 43 398 wet columns is 395.1 m
   against a floor at −450 m, so not one cell of `k = 1` is wet. The seabed a fjord actually has is the
   *immersed* boundary, where Oceananigans' default is free slip.

   The symptom was a domain with **no momentum sink below the surface**. Barotropic circulation had
   nothing to spin it down: a closed circulation around a one-cell-wide land ridge near Horten
   (i = 83, j = 104:106) ratcheted from 0 to 2957 m² s⁻¹ over four simulated days without once
   changing sign, driving a depth-uniform 2.84 m s⁻¹ jet through the two-cell channel beside it and
   collapsing the adaptive time step from 15.5 s to 3.1 s. It kept accelerating after the storm that
   started it had passed — the wind was down to 2.4 m s⁻¹ and the pressure back to 1025 hPa when the
   jet peaked — which is what distinguished it from the wind-driven hotspots elsewhere in the domain,
   all of which decayed. With drag acting, the spin-down time there is H/(C_d·u) ≈ 22 minutes.

   NumericalEarth's own `ocean_simulation` builds both halves
   (`Oceans/ocean_simulation.jl`), but `coupled_simulation` assembles `HydrostaticFreeSurfaceModel`
   directly and so inherits neither — which is why this had to be stated here. Only `bottom` goes
   inside the `ImmersedBoundaryCondition`, matching NumericalEarth: immersed *side* walls stay
   free-slip, which is a separate choice and a separate change.

   ### `OpenLateralBoundaryFromData`

   The third piece, and the one this module exists for. It puts a **genuinely open** boundary on the
   edge the setup's `boundary_config` names, with every exterior value read from the prepared boundary
   file.

   The first version was called `OpenLateralBoundary` and was **not open**: it put
   `NormalFlowBoundaryCondition(nothing)` on the normal velocity, and `getbc` of `nothing` is
   `zero(grid)` — a closed wall. Only the tracers there were open. So the fjord's one connection to
   the Skagerrak passed heat and salt but no water, and the interior relaxation band was left to fake
   the exchange. The first rename, to `OpenLateralBoundaryFromForcing`, said where the values came
   from, which is the other half of what that name got wrong; the second, to `…FromData`, says it
   *correctly* — they are a separate file from a separate collection read by a separate pipeline, and
   they were never forcing.

   That was not an oversight while Oceananigans was at 0.110.13. A data-driven open velocity boundary
   was unusable with `SplitExplicitFreeSurface` there: the barotropic halo fills passed no `clock` and
   no `model_fields`, so any function-, `FieldTimeSeries`- or scheme-valued condition either
   `MethodError`ed inside the kernel or was silently skipped; the barotropic velocities inherited the
   *baroclinic* condition verbatim though `U = Σ Δz·u` is a transport (the source carried a
   `TODO: When open boundary conditions are online … integrate the BC upwards`); and the barotropic
   corrector ran *after* the boundary fill and over the boundary face, undoing it. **Oceananigans
   0.110.15** (PR #5351) closes all three and adds the scheme set below, which is why the
   `[compat]`-free Oceananigans dependency must not be resolved back below it.

   Four groups **per open edge**, each its own generic function with one method per `Val{edge}` and a
   catch-all that raises — never an `if edge === :south` chain. `open_edge_boundary_conditions` merges
   one edge's four, and `boundary_condition_sides` `mapreduce`s that over `open_edges(boundary_config)`
   with the same `recursive_merge`. For a south edge:

   | group | field | condition | scheme |
   |---|---|---|---|
   | `open_normal_velocity_boundary_conditions` | `v` | the `v` series | `NormalFlowBoundaryCondition(…; scheme = NormalRadiation)` |
   | `open_transport_boundary_conditions` | `V` | `south_boundary_transport`, discrete form | `GravityWaveRadiationBoundaryCondition` |
   | `open_tracer_boundary_conditions` | every model tracer | that tracer's series | `ValueBoundaryCondition(…; scheme = NormalRadiation)` |
   | `open_tangential_velocity_boundary_conditions` | `u` | the `u` series | `ValueBoundaryCondition(…; scheme = NormalRadiation)` |

   `GravityWaveRadiation` is Flather (1976) on the barotropic transport and is the part that actually
   opens the boundary; `NormalRadiation` is Orlanski (1976) with Marchesiello et al. (2001) nudging —
   radiate outgoing signals at the locally diagnosed phase speed, nudge towards the data on inflow.

   Five things about this are load-bearing.

   **The merge across edges composes rather than collides.** On a domain open to the south and the
   west, `u` picks up a *tangential* condition at `.south` and a *normal* one at `.west`, and each
   tracer picks up both sides — different leaves of the same nested tuple, which is exactly what
   `recursive_merge` is for. A region in the open ocean naming all four edges needs nothing extra
   here or in Oceananigans, which pairs its Chapman condition onto *every* side carrying a
   gravity-wave transport condition. The one place two edges genuinely meet is a **corner** cell,
   written by both edges' radiation conditions; which of them fills it is Oceananigans' halo-fill
   order, not anything stated here, and it is the part of a multi-edge run to watch rather than to
   assume.

   **`η` is not stated here at all.** Oceananigans pairs a `SurfaceWaveRadiation` (Chapman 1985)
   condition onto every side where `U` or `V` carries a gravity-wave condition
   (`default_free_surface_boundary_conditions` in
   `SplitExplicitFreeSurfaces/split_explicit_free_surface.jl`). A second statement here could only
   disagree with it.

   **Naming a `U`/`V` condition also changes the substepping.** `substep_halo_filling` returns
   `LocalHaloFilling()` when `U` or `V` has a normal-flow condition, so the barotropic halo is filled
   every substep. That is what makes the Flather condition act at all, and it needs no change to
   `SplitExplicitFreeSurfaceConfig`.

   **The velocity and tracer conditions pass the series straight through.** A reduced
   `FieldTimeSeries` *is* a boundary condition Oceananigans indexes by the two tangential indices, so
   the discrete-form wrappers the tracer boundary used to need (`x_boundary_tracer_value`,
   `y_boundary_tracer_value`, `tracer_open_boundary_condition`) are gone. Only the barotropic
   transport needs a function, because `GravityWaveRadiation` wants a 2-tuple.

   **The transport uses the model's own column depth.** `vbar · column_depthᶜᶠᵃ(i, 1, Nz+1, grid,
   model_fields.η)` — the same operator the scheme's own kernel uses — so the exterior transport is
   expressed in the model's bathymetry rather than NorKyst's. Four functions, one per edge, because
   the boundary node is a different index of a different operator on each side; each obeys the kernel
   rules in this file (`@inline`, no short-circuits, literal zeros). No land guard: the series are
   finite everywhere because `fill_boundary_gaps!` already filled them.

   Both timescales are **fields** on the config, `inflow_timescale` and `outflow_timescale`, stated
   per setup with no defaults. The predecessor derived one `relaxation_timescale` from the forcing
   config and used it for both, which was defensible only while the velocity wall made
   `inflow_timescale` inert; a genuinely open boundary treats the two differently, and the outflow
   rate is information no other config states.

   `boundary_condition_sides` gained a sixth positional argument, `boundaries` — the prepared series,
   or `nothing`. It is passed in rather than read here because `build_simulation` is the one place
   `start_date`, which every time axis is zeroed at, is known; a config that opened the file itself
   would need telling that instant a second time, and two statements of it could disagree. The
   **five-argument form is a forwarding fallback** that drops it, so `AirSeaFluxes`,
   `QuadraticBottomDrag` and any out-of-tree config that carries its own exterior state implement the
   shorter signature and never see it. That is what keeps `examples/oslofjorden.jl` — which
   implements the five-argument hook — working untouched.

   The **fourth** argument is the setup's `boundary_config`, and was the *forcing* config while the
   boundary dataset and the open edge both hung off that one. `OpenLateralBoundaryFromData` is the only
   built-in piece that reads it, for `open_edges`. The arity did not change, so a piece that ignores the
   argument is unaffected — both other built-ins, and `examples/oslofjorden.jl`'s
   `RadiatingLateralBoundary`, which carries its own edge. A piece that *reads* it now receives an
   `AbstractBoundaryDataConfig` or `nothing`, which is a change no `MethodError` would have caught, so
   it is stated here and in the hook tables.

   The two hooks have **different names on purpose**: one name returning a nested named tuple for a
   piece and a materialized one for a set would be two shapes behind one name. The per-piece hook is
   `boundary_condition_sides` — named for what it returns, the per-field *sides*
   `field_boundary_conditions` binds as `sides` in its `map`. It was called `boundary_conditions`,
   which collided with `Oceananigans.Fields.boundary_conditions`: that made it the one hook FjordSim
   could not re-export, so an out-of-tree config had to reach it as
   `FjordSim.BoundaryConditions.boundary_conditions`, and any caller binding a local of that name —
   which `build_simulation` used to do — shadowed the hook. Renamed, it is exported like every other
   hook, and the old name is gone rather than deprecated.

   `LATERAL_EDGES` now comes from `Configs`, along with `validate_open_edge` and `lateral_edges`.
   They lived in `Forcing` while that was the first module needing them, and moved when `Bathymetry` —
   which is included *before* `Forcing` — gained `clear_open_boundary_land`. `Configs` is the only
   module every one of them can see, and it already owns `open_edges`, so the edge vocabulary and the
   accessor that reads a config's edges now sit together. This is the tuple's third home; each of the
   first two was a per-module copy with identical contents, an identical membership test and an
   identical error string.

   Only `T` and `S` get a top condition, deliberately: NumericalEarth assembles air-sea fluxes for
   heat and salt specifically, so a third tracer has no exchange to write into and no
   `build_tracer_top_bc` method to build one with. A biogeochemical tracer's surface condition
   belongs in its own `AbstractBoundaryConditionConfig` — and so does its open-boundary series, which
   means a source variable for it in the boundary config's `parameters`.

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
   `z_faces`) and its two grid hooks, plus the constructors they wrap:
   `LatitudeLongitudeGrid(architecture, config::EvenGrid)`, and
   `ImmersedBoundaryGrid(filepath, architecture, halo)`, which reads the processed bathymetry NetCDF
   and returns an `ImmersedBoundaryGrid` wrapping a `LatitudeLongitudeGrid` with `PartialCellBottom`.
   The loader still accepts legacy files with positive depths or swapped `lon`/`lat` axes.

   The grid contract is **two hooks**, `domain_grid(config, architecture)` (the bare geometry every
   `prepare_*` pipeline regrids onto) and `simulation_grid(config, bathymetry_file, architecture)`
   (that geometry with the bathymetry immersed, which `build_simulation` and
   `prepare_forcing(config::FjordConfig)` read the model grid back through). Two rather than one
   because they are needed at different moments from different inputs: the prepare steps have no
   bathymetry file yet.

   Both are declared in **`Configs`**, not here, even though their `EvenGrid` methods live here:
   `Bathymetry`, `Atmospheres` and `Forcing` all call `domain_grid` and are all included *before*
   `Grids`, so the generic has to exist in a module they can already see.

   They are FjordSim generics rather than methods on `Oceananigans.LatitudeLongitudeGrid` and
   `ImmersedBoundaryGrid`, which is what the grid hook used to be. Those names state a *result*
   type, so a config describing a rectilinear, single-column or curvilinear grid could not implement
   the hook without lying about what it returns. The immersed half was worse: it was not a hook at
   all — `build_simulation` reached into `config.grid_config.halo` and called
   `ImmersedBoundaryGrid(file, arch, halo)` itself, so the grid config's only say in the grid the
   simulation ran on was one field. Note the hook is **not** named `target_grid`: that is a
   parameter name in ~112 places across the pipelines, and a function of the same name would shadow
   every one of them.

10. **Simulations** (`src/Simulations.jl`) — `SimulationConfig <: AbstractSimulationConfig`, the
    `build_simulation`/`run_simulation` drivers, and the four of the five nested config
    supertypes' built-in implementations that are not boundary conditions:
    `CoupledHydrostaticSimulation` with `coupled_simulation`, which assembles a
    `HydrostaticFreeSurfaceModel` inside an `OceanSeaIceModel` (NumericalEarth) and returns a
    `Simulation`; `SnapshotWriter`/`CheckpointWriter` with `attach_writer!`; `ProgressCallback` with
    `attach_callback!`; and `AdaptiveTimeStep` with `attach_time_stepping!`. Included after `Grids`,
    whose `simulation_grid` hook it reads the grid back through, and before `Setups`, which
    constructs all of them.

    `coupled_simulation` lives here rather than in `src/FjordSim.jl` for that ordering reason: a
    function defined after the `include` block cannot be imported by a module included inside it.

    `SimulationConfig` has 12 fields and 6 type parameters, one per nested config, and the split
    between it and them is by *what dispatches on what* rather than by convenience. Everything
    grid-independent about the model — `buoyancy`, `closure`, `tracer_advection`,
    `momentum_advection`, `tracers`, `coriolis`, `sea_ice`, `biogeochemistry`, `free_surface` — is
    stored on `CoupledHydrostaticSimulation` as the objects (or, for `free_surface`, the nested
    config) `coupled_simulation` consumes, one type parameter each, nine in total, which is the cost
    of the no-abstract-fields rule. The **tenth** is `extra_kwargs`, and it *is* the collapse the
    ninth was headed for: rather than a field per remaining constructor keyword — of which
    `HydrostaticFreeSurfaceModel` alone has eight — one `NamedTuple` carries all of them. Ten is the
    ceiling; anything further goes inside `extra_kwargs`. The "config" testset guards that count,
    and `SimulationConfig`'s six are the structural ones, growing only when a genuinely new *kind*
    of nested config appears — `callbacks` was the sixth. A new subtype of an existing supertype is
    never a new field.

    ### `extra_kwargs`

    Up to four slots, one per constructor `coupled_simulation` calls, each a `NamedTuple` splatted
    into that call: `ocean_model` → `HydrostaticFreeSurfaceModel` (`clock`, `timestepper`,
    `particles`, `velocities`, `pressure`, `closure_fields`, `auxiliary_fields`,
    `vertical_coordinate`), `ocean_simulation` and `coupled_simulation` → the two `Simulation`s
    (`verbose`, `stop_iteration`, `wall_time_limit`, `align_time_step`, `minimum_relative_step`),
    and `coupled_model` → `OceanSeaIceModel` (`land`, `interfaces`, the reference densities and heat
    capacities, and its `interface_kw…`). An omitted slot is empty, so all four splat
    unconditionally. Before this, needing any one of those meant reimplementing `coupled_simulation`
    wholesale for a model that differed from the built-in one by a single keyword.

    The escape hatch can only **add**. `validate_extra_kwargs` rejects a slot that reaches no
    constructor, a slot holding something that is not a `NamedTuple`, and — the load-bearing one — a
    key duplicating a keyword `coupled_simulation` already passes. `extra_kwargs` splats *after* the
    explicit keywords, so a duplicate wins silently; a stray `tracers` in `ocean_model` would then
    disagree with `model_tracers`, which `build_simulation` already used to pick the forcing terms
    and the open tracer boundaries before the model existed. The lists live in `EXTRA_KWARG_SLOTS`
    rather than inside `coupled_simulation`, so the check runs at construction, where the setup file
    that made the mistake is still what the error is about.

    `free_surface` is itself an `AbstractFreeSurfaceConfig` rather than a bare `Float64` CFL, with
    its own `free_surface(config, grid)` hook — the same reason the model as a whole is a config
    rather than a built object: `SplitExplicitFreeSurface` needs the grid, which does not exist
    until `coupled_simulation` calls the hook.

    ### `closure` and `BoundarySponge`

    `closure` holds *either* a pre-built Oceananigans closure — which is what every setup wrote before
    this existed, and still the usual case — *or* an `AbstractClosureConfig`. Which one it is is
    resolved by `model_closure(closure, grid, boundary_config)`, whose fallback is the identity, so
    the two shapes are three methods rather than a branch and nothing that already worked changed.
    Same pattern as `resolve_initial_conditions`, and for the same reason as `free_surface`: a closure
    can need what only exists at build time.

    `BoundarySponge` is the built-in one, and it needs *two* such things: the grid, for the domain
    extent and cell size its ramp is expressed in, and the open edges. It wraps a `base` closure and
    appends a `HorizontalScalarDiffusivity` whose `ν` and `κ` ramp to zero over `width_cells` inward
    from every open lateral edge, so it needs no field, no allocation and no architecture — the same
    coefficient function runs on CPU and GPU.

    It exists because an open lateral boundary radiates as well as admits, and nothing was absorbing
    what it radiated. On the 2020 `oslofjorden` run the boundary row carried velocities of 0.42 m/s
    standard deviation and 2.04 m/s peak against 0.13 and 0.43 in the interior, and grid-scale
    salinity roughness of 2.18 against 0.03 seventy rows in. That accumulated for fifty days in the
    poorly ventilated near-boundary bottom cells and then ran away: max salinity went 54 → 86 → 192
    psu over the last thirteen days while the *domain mean* moved by 1 psu, which is what identifies
    it as redistribution rather than a source. The prepared boundary file was asking for 24–34 psu
    throughout, so the exterior data was never the problem.

    Four things about it are load-bearing.

    **It is viscous, not a tracer relaxation band.** A band relaxing T and S a few cells inside the
    domain would fight the open boundary condition, which is nudging the same variables at the
    boundary towards the same data — which is exactly why the old interior relaxation band was deleted
    when the open boundary arrived. Viscosity and diffusivity name no target and so cannot disagree
    with one. This is the sponge that section said would be needed "if a run turns out to need it".

    **The edges are not a field.** They come from `open_edges(boundary_config)`, which is why
    `coupled_simulation` gained a `boundary_config` keyword — the same argument, for the same reason,
    that `boundary_condition_sides` already receives. A setup that later opens a second edge sponges
    both with no second edit, and one naming no boundary config gets `base` back unchanged rather than
    a sponge of zero strength.

    **The coefficient is bounded by stability, not by taste.** Explicit horizontal diffusion needs
    `Δt ≤ Δx²/4ν`, which on `oslofjorden`'s 193 m cell is 310 s at `ν = 30` and 155 s at `ν = 60`. A
    value that pushes that bound under the time-stepping config's `max_time_step` caps the run's time
    step instead of the CFL doing it. Raise one and check the other.

    **The ramp is `max(0, 1 - d/w)^2`, squared rather than linear**, so its derivative vanishes at the
    inner end too: a ramp reaching zero with a kink is itself a discontinuity in the momentum equation,
    which is the sort of thing an open boundary reflects off. `sponge_strength` takes the largest ramp
    over the four edges with a `1.0`/`0.0` multiplier per edge rather than a branch, so a corner
    belongs to whichever edge it is nearer and the kernel is one straight-line expression whatever
    subset the setup opened. A Gaussian `exp(-(d/w)^2)` is the other common choice
    and is what the reference implementation below uses; the squared ramp is preferred here only
    because it is compactly supported, so no viscosity leaks into the interior at all.

    **It is a `discrete_form = true` diffusivity, and that is not a style choice.** `ScalarDiffusivity`
    only threads `parameters` through in its discrete form. Asked for the continuous one it stores the
    bare function and calls it as `ν(λ, φ, z, t)` — four arguments, `parameters` silently dropped,
    despite the constructor having accepted them. A parameterized method then never matches, and
    because the call happens inside a GPU kernel the `MethodError` surfaces as
    `InvalidIRError: unsupported call to an unknown function (call to jl_f_throw_methoderror)`
    reported against an unrelated frame in the momentum tendency kernel — eight minutes into a build,
    naming nothing that would lead you to the closure. The "boundary sponge" testset calls the
    coefficient back through Oceananigans' own `νᶜᶜᶜ`/`κᶜᶜᶜ` dispatch precisely so that mismatch is a
    test failure rather than a compile failure an hour into a run.

    ### The restoring sponge, and why this is not one

    The other standard shape is a `Relaxation` forcing over a near-boundary mask, pulling each field
    towards a prescribed exterior state — NumericalEarth's `DatasetRestoring` is the ready-made
    version. A regional Arctic configuration built that way nudges `T` and `S` at `1/1day` and `u` and
    `v` at `1/20minutes` over a Gaussian mask four cells wide: velocities some seventy times harder
    than tracers, on the reasoning that it is the *momentum* field near the boundary that has to stay
    matched to what is prescribed.

    That asymmetry agrees with the Oslofjord diagnosis from the other end — the boundary row's velocity
    was 3–7x the interior's and the tracer extremes followed from it — and is why `viscosity` sits
    above `diffusivity` here rather than equal to it.

    FjordSim cannot simply adopt the restoring form. A relaxation forcing needs a *target throughout
    the band*, and the exterior state a setup prepares exists only **at** the boundary row:
    `boundary_series` returns reduced `FieldTimeSeries`, one cell thick. The only 3D target on the
    model grid is `forcing.nc`, whose lambdas `prepare_forcing` deliberately writes as zero.
    Reinstating a band there **for `u` and `v` alone** is the faithful translation, and it escapes the
    objection that killed the old tracer band — that band fought the open boundary's own nudging of
    the same tracers towards the same data, and a velocity band touches no tracer. That is the next
    thing to try if this sponge proves too blunt.

    `coefficient` stays a plain scalar on
    `QuadraticBottomDrag` rather than following the same pattern, because
    `boundary_condition_sides` is already the dispatched hook — there is no second function for it
    to delegate to the way `coupled_simulation` delegates to `free_surface`.
    `architecture` is a `Symbol` for the same reason as the forcing config's, and
    `simulation_architecture` resolves it by reusing `interpolation_architecture`'s `Val` methods.

    `model_tracers(model)` is a hook rather than a field read because `build_simulation` needs the
    tracer list before the model exists: `simulation_forcing`, `OpenLateralBoundaryFromData` and
    `resolve_initial_conditions` each build one thing per tracer.

    No field has a default anywhere in this stack, deliberately, which is the one place these configs
    depart from the rest of FjordSim: a defaulted closure or coriolis would be one fjord's physics
    quietly applied to another's. The setup file is the complete statement of the run.

    Three things are deliberately *not* fields, because they are already stated elsewhere and a
    second copy could only disagree: the forcing file (`simulation_forcing_path` takes the
    rivers-augmented copy when the forcing config names rivers, the plain prepared file
    otherwise), the open boundary (`OpenLateralBoundaryFromData` puts its schemes on the edge
    `config.boundary_config` names, and reads its exterior state from that same config's prepared
    file), and the atmosphere (the `prescribed_atmosphere` and
    `prescribed_radiation` hooks on `config.atmosphere_config`).

    `build_simulation` reads the boundary series once, with
    `boundary_series(config.boundary_config, grid, start_date)`, and hands it to
    `field_boundary_conditions` as its sixth argument — the same shape as `forcing`, and for the same
    reason: this is the one place `start_date` is known, so it is the one place a time axis can be
    zeroed at it. Both the config and the edge come straight off `FjordConfig` now; the
    `boundary_data_config` and `open_edge` one-liners that used to reach through the forcing config are
    gone, and `open_edges` is a `Configs` accessor on the boundary config instead, returning a
    `Vector{Symbol}` — empty for a setup that names none, so every consumer iterates.

    `start_date` is the exception to that list, and it is a field precisely *because* it cannot be
    derived. Every prepared file carries its own first record, and each reader used to zero its own
    time axis there — the Oslofjord forcing starts at 12:00 and its NORA3 atmosphere at 00:00, so
    the two ran twelve hours out of phase with nothing reporting it. One stated instant is the only
    thing they can all agree on, so `build_simulation` passes it to `simulation_forcing` and to both
    atmosphere read hooks as their `reference_date`. Picking a `start_date` earlier than a file's
    first record is therefore an error rather than a shift.

    `validate_time_coverage` enforces that: each prepared file must span
    `[start_date, start_date + stop_time]` — the forcing, the atmosphere, and the open-boundary data,
    through `boundary_date_range`. Every reader uses `Cyclical()` time indexing, which does
    not fail outside the data it was given — it wraps — so a run that outlasts its forcing would
    otherwise quietly replay the beginning, and one starting early would read the end. The check is
    what makes the wrap unreachable rather than merely unlikely, which is why neither reader's
    `time_indexing` had to change. **Both** halves now go through an optional hook on their own
    config — `atmosphere_date_range(config)` and `forcing_date_range(config, filepath)` — so the
    module names no dataset and no file layout. The forcing half used to be a bare
    `forcing_date_range(filepath)` living here, dispatched on nothing, hardcoding an `NCDataset` read
    of `ds["time"]`: a source that overrode `simulation_forcing` to read a different layout still had
    its coverage validated by a reader it did not use. The supertype default is that same NetCDF
    read, so no built-in source changed.

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
    `coupled_simulation`'s single `set!(ocean_model; initial_conditions...)` consumes. That function
    therefore learns nothing about the three shapes; the dispatch is three methods, not a branch.

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

    ### Callbacks

    `callbacks` is a tuple of `AbstractCallbackConfig`s, and `attach_callbacks!` loops over it
    calling `attach_callback!` — the same shape as `writers`, and `validate_callbacks` rejects two
    under one name for the same reason `validate_writers` does: the second replaces the first in
    `simulation.callbacks`, so only the last would ever fire.

    `ProgressCallback(; name, interval, report)` is the built-in one. All three fields matter: what a
    run reported used to be a hardcoded `Callback(progress, TimeInterval(config.progress_interval))`
    under the fixed key `:progress`, so a setup could change how *often* it reported and nothing
    else — not the report, not the schedule, and not by adding a second diagnostic.
    `progress_interval` is gone from `SimulationConfig` as a result.

    `report` is the field that matters most, because `Utils.progress` reaches
    `sim.model.ocean.model.tracers.T`. A model whose `tracers` omits `:T` — perfectly legal, since
    `model_tracers` exists precisely so tracers are a free choice — crashes at the first fire, after
    the whole coupled model has compiled and the run has started. Naming its own `report` is the way
    out. `Utils.progress` itself is unchanged, and so is `Utils.WALL_TIME`, a module-global `Ref` two
    simulations in one session still interleave.

    ### Writers

    `writers` is a tuple of `AbstractWriterConfig`s, and `attach_writers!` just loops over it calling
    `attach_writer!`. Which simulation a writer attaches to is the method's business, not the
    caller's: a `SnapshotWriter` and a `FieldSnapshotWriter` go on `simulation.model.ocean`, a
    `CheckpointWriter` on the coupled simulation, and a `ProgressCallback` on the coupled one too.

    A `SnapshotWriter`'s `variables` are `Symbol`s resolved by `snapshot_outputs` through
    `Oceananigans.fields(ocean_model)` — velocities, free surface, tracers, auxiliaries — so naming
    `:w` or a biogeochemical tracer is enough to have it written and nothing is enumerated here.

    ### `FieldSnapshotWriter`, and the one field NetCDF will not take

    `FieldSnapshotWriter` is the same idea written to JLD2, and it exists for exactly one reason:
    **Oceananigans' NetCDF writer cannot emit a `(Center, Center, Nothing)` user output at all.** The
    free surface `η` asks for a singleton `z_aaf = [0.0]` while the grid's own vertical coordinate
    already owns that name with the real faces, and `create_spatial_dimensions!` raises rather than
    reconciling them. Measured on `oslofjorden`: it fails beside the 3D fields, alone in a file of its
    own, and with `include_grid_metrics = false`. `bottom_height` is the same location and *is*
    written, because grid metrics take a different code path — which is what identifies this as an
    upstream defect in the user-output path rather than something a setup can configure around.

    So `η` goes in a second file, and the two writers split by *what a format can hold* rather than by
    subject: 3D fields in NetCDF, which is what every reader here expects, and z-reduced fields in
    JLD2, which serializes an array and has no dimension table to collide with. They share
    `snapshot_outputs` rather than duplicating it — which names a model has, and the error naming the
    ones it does not, are the same question whatever the format.

    The JLD2 layout is Oceananigans', not FjordSim's: `timeseries/<name>/<iteration>` holds one
    `Float32` array per record and `timeseries/t/<iteration>` the model time in seconds, with
    `with_halos = false` so an array is exactly the interior — `(Nx, Ny, 1)` for `η`. No grid is
    stored, so a reader needs the bathymetry file for geometry.
    A name the model lacks is an **error**, unlike `state_variables` on the read side, which
    intersects and moves on. The asymmetry is deliberate: that function serves two kinds of file whose
    variable sets legitimately differ and which no config named, while here the setup wrote the name
    down. An over-eager error costs a typo fix; a silently dropped variable costs a whole run.

    Two Holy traits replace what used to be threshold tests on floats. `checkpoint_trait` answers
    whether a writer contributes a `Checkpointer` — needed at three sites, and needed *exactly*, since
    `run!(…; checkpoint_at_end)` with none to find writes into `pwd()` behind a `@warn`. It is also
    what lets `resume_loop` reject `pickup` without a checkpointing writer, which was not expressible
    while checkpointing was a number. `output_path_trait` answers whether a writer names a file the
    run should report, and is **not** the negation of the first: a checkpointer writes files, but
    scratch rather than product ones.

    `validate_writers` rejects two checkpointing writers (`run!` cannot tell which to resume from) and
    two writers under one key (the second silently replaces the first in `output_writers`), before
    anything is read or allocated.

    `mkpath` is load-bearing in two places and in neither of the obvious ones. `build_simulation`
    creates `results_root`, and the snapshot `attach_writer!` creates its own file's directory —
    `NetCDFWriter` creates only the `dir` keyword it is given, and FjordSim passes the whole path as
    `filename` instead. `coupled_simulation` deliberately does not, being about assembling a model.

    ### Checkpointing

    Naming a `CheckpointWriter` among `writers` is what makes a run resumable; naming none writes no
    checkpoints, which is how checkpointing is turned off. It goes on the **coupled** simulation, not
    the ocean one, for two reasons that both fail otherwise: `prognostic_state` of the coupled model
    is what a resumable state is, and `run!(…; pickup)` looks for its checkpointer in
    `simulation.output_writers` and requires exactly one there.

    Oceananigans 0.110 checkpoints *only* `prognostic_state(simulation)`, so the `Checkpointer`
    docstring's warning that "objects containing functions cannot be serialized" does not apply here:
    `model.forcing` and `model.boundary_conditions` are not in `HydrostaticFreeSurfaceModel`'s
    prognostic state at all, `PrescribedAtmosphere` and `PrescribedRadiation` contribute only their
    clock, and `ComponentInterfaces` and `FreezingLimitedOceanTemperature` contribute `nothing`. So no
    `ForcingFromFile`, no `FieldTimeSeries` backend and no `FreshwaterExchange` is ever written — while
    CATKE's diffusivities and its `e` tracer are. Expect a few hundred MB per checkpoint on the
    Oslofjord grid, which is why `cleanup = true`.

    `pickup` resumes from the newest checkpoint. The checkpoint prefix carries the loop index
    (`checkpoint_prefix`), because the state records the clock but not which repetition produced it —
    without that, `pickup` could not tell a loop-3 checkpoint from a loop-1 one and would replay the
    whole spin-up, and loop 2 would overwrite loop 1's files at the same iteration number.
    `resume_loop` reads the index back out of the filename, and `build_simulation` attaches its
    writers for *that* loop rather than unconditionally for the first.

    It carries **no run tag**, unlike everything else a run writes, and that is what makes `pickup`
    work at all now that the tag is the launch instant: a later launch could not name — and so could
    not find — the checkpoints of the one before it. Three things follow, and each fails silently.
    Checkpoints are shared per `results_root`, so there is **one resumable run per results directory**:
    a second launch overwrites the first's (`Checkpointer.write_output!` opens `jldopen(path, "w")`
    whatever `overwrite_existing` says — which is why `CheckpointWriter` does not name that keyword at
    all, `Checkpointer` storing it and never reading it — and `cleanup = true` prunes the rest,
    a filename carrying the iteration meaning every fire writes a new file). A `pickup` therefore
    resumes whatever state is there, even one written under a different `start_date`. And the snapshots
    *are* per launch, so a resumed run's file starts at the resume point and the records before it stay
    in the previous launch's file — which is why the snapshot `attach_writer!` does not suppress
    `overwrite_existing` on a pickup: the file it writes is always new.

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

12. **CLI** (`src/CLI.jl`) — `SUBCOMMANDS` maps each subcommand to the driver it calls — nine of
    them now, `download_boundaries` and `prepare_boundaries` sitting between `add_rivers` and the
    atmosphere pair, in the order a setup is prepared — one `USAGE`
    string, a pure `parse_arguments` returning `(; subcommand, config, help)` and throwing
    `ArgumentError`, and `main(args)` returning a process exit code. Included last, since it names
    every driver and every setup. `parse_arguments` deliberately does not `exit`: printing and exit
    codes live in `main`, which keeps the help path testable. Exit codes are 0 success, 1 the step
    failed, 2 bad arguments.

    `run_step(driver, config, subcommand, setup)` runs one driver and turns its outcome into an exit
    code; `main` calls it inside `tee_output(f, log_path(simulation_config))` for **`run_simulation`
    only**, and directly for every other subcommand. `tee_output` mirrors `stdout` and `stderr` to the
    log while still printing live to the terminal. `redirect_stdio` only accepts fd-backed streams, so
    the tee is a `Pipe` (needing an explicit `Base.link_pipe!`; an uninitialized `Pipe` throws from
    `eof`) plus a task copying each chunk to both destinations. Both streams share one pipe, so the log
    interleaves them in write order.

    Only the simulation is logged because it is the one step whose failure is a stacktrace through the
    whole coupled model — long enough to push the error message itself out of the scrollback — the one
    whose output is not reproducible by re-running it, and the only one whose setup is guaranteed to
    name a `results_root` to put a transcript in. The prepare steps print a few lines and one error.
    That is also what retired `LOG_FILE` (`fjordsim.log` in the working directory) and `log_path`'s
    `FjordConfig`/`Nothing` methods: the fallback existed for a setup with no simulation config, which
    is exactly the case where `run_simulation` is a no-op and nothing is logged at all. `log_path` now
    takes the simulation config, and `.gitignore` no longer needs the name.

    `log_path` puts the transcript under `results_root`, beside the output it describes, and creates
    that directory if absent. Its name carries the run tag — `fjordsim_<run_tag>.log` — so, the tag
    being the launch instant, each launch keeps its own transcript instead of truncating the one
    before it.

    Errors go through `show_compact_error`, which sets `:limit`, `:displaysize` (`STACKTRACE_WIDTH`,
    120 columns) and `:stacktrace_types_limited` on the IO before `showerror`. That last one is the
    whole trick: Julia abbreviates the type parameters of a frame to `{…}` only when the IO carries
    it, which is what `REPL.repl_display_error` does and what a script's bare `showerror` does not — so
    a frame that used to spell out every parameter of a `FjordConfig` or a `HydrostaticFreeSurfaceModel`
    now fits on one line. `Base.type_depth_limit` floors the width at 120, so a smaller
    `STACKTRACE_WIDTH` would change nothing; stating `:displaysize` at all is what keeps the log
    independent of terminal width. It shortens *frames*, not messages — a `GPUCompiler`
    `InvalidIRError` body is as long as it ever was.

    Two things about the tee are load-bearing. The failure is **caught inside** the redirect, by
    `run_step`: an exception left to propagate would be printed by `Base._start` after `tee_output`'s
    `finally` had torn the redirect down, so the error would be the one thing missing from the log. And
    only the driver runs inside the tee — parsing, `--help` and config resolution stay outside it, so a
    usage error leaves no log file behind, which is also why the existing `main` tests (all of which
    fail before the driver is reached) write nothing.

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
    exporting it would make every `using FjordSim` in a script — `test/runtests.jl`, or a config
    file run directly rather than through `--config` — run the CLI on the way out. And the `@main` is
    bare, *after* the definition: `@main function main(args) ... end` expands to a **call**, which
    would run the CLI while the package precompiles.

## Adding a new source

Every pipeline is a generic function on the config supertype plus a small set of hooks. Add a
source by subtyping and overloading the hooks — never by editing the generic function. The
adapter files (`src/Bathymetry/geonorge.jl`, `src/Forcing/norkyst.jl`, `src/Forcing/of800_rivers.jl`)
are the templates.

Grid — `AbstractGridConfig`:

| Hook | Required |
|---|---|
| `domain_grid(config, architecture)` → the bare grid, no bathymetry | yes |
| `simulation_grid(config, bathymetry_file, architecture)` → the grid the simulation runs on | yes |

Both are declared in `Configs`, since `Bathymetry`, `Atmospheres` and `Forcing` all call
`domain_grid` and are included before `Grids`. Neither names a result type, so a grid that is not a
`LatitudeLongitudeGrid` can implement them.

Bathymetry — `AbstractBathymetryConfig`, consumed by `prepare_bathymetry`:

| Hook | Required | Default |
|---|---|---|
| `bathymetry_dataset(target_grid, config)` → NumericalEarth dataset | yes | none |
| `regrid_options(config)` → NamedTuple for `regrid_bathymetry` | no | `(;)` |
| `smoothing_options(config)` → NamedTuple for `smooth_bathymetry_gaps!` — `open_boundary_land_cells`, `max_island_cells`, `close_narrow_passages`, `spike_ratio`, `minimum_cell_fraction`, `max_slope_factor`, `minimum_depth` | no | `(;)`, leaving only the topological cleanup |

The `edges` keyword `prepare_bathymetry` takes is a pipeline *argument*, not a hook, exactly as
`prepare_forcing`'s is: the `FjordConfig` driver reads it with `open_edges(config.boundary_config)`
and every source inherits `clear_open_boundary_land` unchanged.

Forcing — `AbstractForcingConfig`, consumed by `prepare_forcing`:

| Hook | Required |
|---|---|
| `forcing_time_steps(config)` → `Vector{SourceRecord}` | yes |
| `forcing_source_grid(config, filepath)` → source grid | yes |
| `forcing_variable_names(config)` → `Dict` source name => FjordSim name | yes |
| `download_forcing(target_grid, config)` | only if it downloads |
| `simulation_forcing(config, grid, filepath, tracers, reference_date)` → the forcing term object `coupled_simulation` consumes | no; defaults to `forcing_from_file`, the FjordSim NetCDF contract |
| `forcing_date_range(config, filepath)` → `(first, last)` `DateTime`s | no; defaults to reading `ds["time"]` from the prepared NetCDF. A source that overrides `simulation_forcing` overrides this too |
| `source_field_grid(source, architecture)`, `projected_target_nodes(longitude, latitude, source)` | only for a source grid that is not a regular projected grid; dispatch on the source-grid type, not the config |

Rivers — `AbstractRiverConfig`, consumed by `add_rivers`:

| Hook | Required | Default |
|---|---|---|
| `river_locations(config)` → `Vector{RiverLocation}` | yes | none |
| `river_series(config, times)` → `Dict` FjordSim name => `(river, time)` matrix | yes | none |
| `download_rivers(config)` | only if it downloads | none |
| `river_search_radius(config)` → cells to search for a coastal cell | no | `config.search_radius` |
| `river_minimum_levels(config)` → wet levels a column must have to receive an outlet | no | `0`, accepting any water cell. Unlike `river_search_radius` the fallback reads no field, so a river config written before the hook existed keeps working; `OF800RiversConfig` overloads it |

Its `standalone` field decides whether `add_rivers` patches a copy of the prepared forcing or writes a
river-only file of its own; see the `Forcing` section. A `standalone` config needs the setup to name a
`simulation_config`, since the file's time axis comes from the run window.

A river config is not a `FjordConfig` field — it goes in the forcing config's `rivers` field,
`nothing` for a setup with no rivers, because rivers are written into the forcing file itself. A
variable `river_series` returns that the forcing file does not carry is skipped with a warning, so one
river dataset can serve setups that prepare different variables — vacuously so for a `standalone`
config, whose file carries exactly what `river_series` returned.

Open-boundary data — `AbstractBoundaryDataConfig`, consumed by `prepare_boundaries`:

| Hook | Required | Default |
|---|---|---|
| `boundary_time_steps(config)` → `Vector{SourceRecord}` | yes | none |
| `boundary_source_grid(config, filepath)` → source grid | yes | none |
| `boundary_variable_names(config)` → `Dict` source name => FjordSim name | yes | none |
| `download_boundaries(target_grid, config)` | only if it downloads; reads its edges with `open_edges(config)`, and subsets the box `boundary_domain` derives from them | none |
| `boundary_date_range(config)` → `(first, last)` `DateTime`s | no | reads `ds["time"]` from the prepared NetCDF |
| `boundary_source_slab(config, reader, step, source_name)` → one source slab | no | `blended_slab`, a plain read. Override only for a variable *derived* from more than one of the source's own — a velocity component stated along the source's own grid axes is the case that needs it |

Unlike a river config, this **is** a `FjordConfig` field, `boundary_config`, `nothing` for a setup
whose lateral boundary is not data-driven. It hung off the forcing config's `boundaries` field until a
setup needed an open boundary with no interior forcing; the exterior state was always its own hourly
file from its own collection read by its own pipeline, and either config is nameable without the other.

The edges **are** a field of it, `open_edges`, read through the `open_edges` accessor and stated nowhere
else. Write one `Symbol` or a collection of them — a domain in the open ocean names all four — and the
built-in config's keyword constructor normalizes either to a `Vector{Symbol}` with `lateral_edges`, so
every consumer iterates one shape rather than branching on a scalar. The three pipeline entry points
take the config and read the edges off it rather than taking both, which is what keeps them from
disagreeing; `prepare_forcing` takes them as an `edges` keyword, since a forcing dataset has no
business naming them, and `nothing` there means a closed domain. Not to be
confused with `AbstractBoundaryConditionConfig`, which is the scheme rather than the data.

`prepare_boundaries` reuses the forcing core on a one-cell-thick target slab, so a boundary source on
a regular projected grid inherits `source_field_grid`, `projected_target_nodes`, `SourceFill` and the
interpolation kernel unchanged, exactly as a forcing source does. It loops the whole per-variable build
over the config's edges and writes them all into one file, so neither hook above becomes per-edge — the
cost of that is the download: `boundary_domain` returns the bounding box of the edges' bands, which for
two opposite edges is the whole domain plus the margin, so a multi-edge download is the full box at
hourly cadence and its interior is fetched and never read. Per-edge downloads are the optimization, at
the price of `boundary_time_steps` and `boundary_source_grid` becoming per-edge. What the prepared file's variables
are *named* and *shaped* is not a hook — `boundary_variable_name`, `boundary_dimension_names` and
`boundary_location` fix it, because the read side does.

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
this hook and `simulation_forcing`, which is the only thing keeping the two in phase — see the
`Simulations` section.

The prepared variable names and units are *not* a hook — they are fixed by the read side in
`ATMOSPHERE_VARIABLES`. A source whose download already normalizes names (as NORA3's does, since
five of its eight variables are derived rather than copied) returns the identity mapping. The
core's `grid_rotation_angle` and `rotate_to_east_north` are available to any adapter whose source
gives wind relative to its own grid axes.

Simulation — `AbstractSimulationConfig`, consumed by `build_simulation`: **no hooks**. It is data
only, read field by field, so an alternative simulation config is a subtype supplying the field
set its docstring lists. It inherits `results_path`, `run_tag` and `coverage_window` — the last reads
`start_date`, so a subtype omitting that field inherits only the first two. Everything that *is*
dispatched on comes from the supertypes it nests, one hook table each. `AbstractCallbackConfig` is
the newest row: it replaced a bare `progress_interval::Float64`, which could only change how often
the one hardcoded `Callback(progress, TimeInterval(...))` fired.

| Supertype | Field | Required hooks |
|---|---|---|
| `AbstractCoupledSimulationConfig` | `model` | `coupled_simulation(model, grid; forcing, boundary_conditions, initial_conditions, atmosphere, radiation, boundary_config, stop_time, initial_time_step)`, `model_tracers(model)` |
| `AbstractBoundaryConditionSetConfig` | `boundary_conditions` | `field_boundary_conditions(config, grid, forcing, boundary_config, tracers, boundaries)` |
| `AbstractBoundaryConditionConfig` | the pieces inside the set | `boundary_condition_sides(config, grid, forcing, boundary_config, tracers[, boundaries])` — implement the six-argument form only if the piece reads the prepared open-boundary state; the five-argument one is reached through a forwarding fallback. The fourth argument was the *forcing* config before the open edge moved, so a piece that reads it needs checking |
| `AbstractCallbackConfig` | `callbacks` (a tuple) | `attach_callback!(simulation, callback, config)` |
| `AbstractWriterConfig` | `writers` (a tuple) | `attach_writer!(simulation, writer, config, loop)`; optionally `checkpoint_trait`, `output_path_trait`, both defaulting to the negative. `SnapshotWriter` (NetCDF), `FieldSnapshotWriter` (JLD2, for z-reduced fields NetCDF rejects) and `CheckpointWriter` are the built-ins |
| `AbstractTimeSteppingConfig` | `time_stepping` | `attach_time_stepping!(simulation, config)`, `initial_time_step(config)` |

A writer needs an `output_file` field only when its `output_path_trait` is `NamesOutputFile`; on any
other, `results_path` raises a stated `ArgumentError` rather than failing on a missing field.

`CoupledHydrostaticSimulation`, the built-in `AbstractCoupledSimulationConfig`, can nest two further
supertypes of its own. Neither is one of the rows above — those are what `SimulationConfig` itself
nests — but the same reasoning applies one level down: a different free surface or closure is a new
subtype and a new method, not an edit to `coupled_simulation`.

| Supertype | Field | Required hook | Built-in |
|---|---|---|---|
| `AbstractFreeSurfaceConfig` | `free_surface` | `free_surface(config, grid)` → the free-surface object `coupled_simulation` passes to `HydrostaticFreeSurfaceModel` | `SplitExplicitFreeSurfaceConfig` |
| `AbstractClosureConfig` | `closure` | `model_closure(config, grid, boundary_config)` → the closure, or tuple of closures, it passes | `BoundarySponge` |

The two differ in one way that matters. `free_surface` is *always* a config — there is nothing else
it could hold. `closure` holds a config **or** a pre-built Oceananigans closure, and usually the
latter; `model_closure`'s fallback is the identity, so the two shapes are resolved by dispatch and a
setup naming a plain closure needs no change and never sees the hook. Same pattern as
`initial_conditions` and `resolve_initial_conditions`.

`model_closure` takes `boundary_config` as well as the grid because a closure may act on the open
edges, which are stated only on the `AbstractBoundaryDataConfig`. That is why `coupled_simulation`
carries `boundary_config` too.

`CoupledHydrostaticSimulation` also carries `extra_kwargs`, which is not a hook but the reason a
setup rarely needs a *new* `AbstractCoupledSimulationConfig` subtype: it splats one `NamedTuple` per
slot into each of the four constructors `coupled_simulation` calls, so a model differing from the
built-in one by a keyword is a config value rather than a new type and a new method. See the
`Simulations` section for the slot names and for why the constructor rejects keys that would shadow
the fields.

The `coverage` keyword `prepare_forcing`, `prepare_boundaries` and `prepare_atmosphere` take is a
pipeline *argument*, not a hook: the `FjordConfig` drivers derive it with `coverage_window` and every
source inherits the padding unchanged. No source gains a hook for it.

Required fields are listed in each supertype's docstring in `src/Configs.jl`. Path resolution
(`bathymetry_path`, `forcing_path`, `forcing_directory`, `river_forcing_path`, `boundary_data_path`,
`boundary_data_directory`, `atmosphere_path`,
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
both `oslofjorden()` and `drammensfjorden()` give their `OF800RiversConfig` the same `data_root` as
the rest of the setup — the river data downloads there rather than being shared from elsewhere, so
each setup carries its own copy of the ~176 MB series file. The same goes for their
`NorKystBoundariesConfig`, and there it is not merely tidiness: the downloaded band is derived from
*that* setup's own open edge, and Drammensfjord's southern edge is 20 km north of Oslofjord's, so the
two bands are different data. A setup wanting no rivers or no boundary data leaves the
field unnamed so it defaults to `nothing`, making that step a no-op; because `rivers` is a type
parameter of `NorKystConfig` and `boundary_config` one of `FjordConfig`, that is a construction-time
choice in both cases and cannot be undone on an existing
instance. `forcing_config`, `atmosphere_config` and `simulation_config` default to `nothing` the same
way, making both
atmosphere steps and `run_simulation` no-ops — though both built-in setups name all three.

The simulation config is rooted separately, at `~/FjordSim_results/<fjord>/` rather than under
`data_root`, since it writes rather than reads. Unlike every other config, `SimulationConfig` has
**no defaults at all**, and neither do the configs it nests: `oslofjorden()` names every field
of all of them — `extra_kwargs = (;)` included — because each is a scientific choice about that
fjord and a default would let the next
setup silently inherit it. Adding a field to any of them therefore breaks every setup until each
names it, which is the intent.

### The vertical grid

`oslofjorden`'s `z_faces` are a geometric stretch — ratio 1.25, capped at 33.5 m — from a 1 m surface
cell to a deepest face at −400 m, 24 levels. They replaced a hand-written 18-level set running to
−450 m in 25 m and 50 m steps, and the numbers are worth keeping because the replacement was measured
against the actual bathymetry rather than chosen:

| | 18 levels, to −450 m | 24 levels, to −400 m |
|---|---|---|
| median thickness of the layer holding a column's floor | **25.0 m** | **9.0 m** |
| columns floored in a layer thicker than 20 m | **61.0 %** | **24.0 %** |
| largest ratio between adjacent layers | 2.50 | 1.33 |
| levels that are never wet | 1 | 0 |
| median wet levels per column | 9 | 11 |
| bottom cells thinner than 0.4 of their layer | 20.9 % | 20.3 % |
| bottom cells `PartialCellBottom` has to clamp (true slivers) | 292 | **0** |

The old grid's failure was the first row. A 25 m bottom cell is the entire near-bottom water column
of that site in one cell, with no vertical structure in it, and that is precisely where the 2020
run's tracer extremes sat: 96 % of the temperature offenders and 67 % of the salinity offenders were
in a bottom cell, and those were overwhelmingly at k = 9 and k = 10, the two 25 m layers. The first
layer, −450 to −400 m, held no water at all against a deepest sounding of 395.1 m, and `grid.Lz`
feeds `sqrt(g·Lz)` in `SplitExplicitFreeSurface`, so the unused depth was also buying a shorter
barotropic substep for nothing.

Three things to know before changing it.

**−400 m is deliberate headroom, not a coincidence.** The deepest sounding is 395.1 m, and
`limit_bottom_slope` can only make the deepest column *shallower* (it is the deeper half of every
pair it corrects), so the margin cannot close. A sounding below the deepest face is not an error —
`snap_partial_bottom_cells` skips it and `PartialCellBottom` clips it — so the basin would simply be
silently truncated.

**More levels is not free, and not monotonically better.** A finer grid gives the seabed more faces
to cross, so laterally isolated bottom cells rise from 3.9 % to 5.5 % of columns — see
`snap_partial_bottom_cells`, which is what handles them. A 28-level ratio-1.20 variant scores better
still (median floor layer 8.4 m, *no* column floored in a layer over 30 m) at 56 % more cells than
the original rather than 33 %; it is the option to reach for if the deep basins still look
under-resolved.

**Changing it is a data change.** `z_faces` is written into `bathymetry.nc`, and `simulation_grid`
reads the file rather than the config — so the run silently uses whatever the file holds. Every 3D
prepared file is regridded onto that grid, so re-run `prepare_bathymetry` → `prepare_forcing` →
`add_rivers` → `prepare_boundaries`. No download step is affected. The atmosphere is 2D and is not.

### The 2020 diagnosis

Several of `oslofjorden()`'s current values are there because of one run, and are commented in the
setup file with the number that justifies them. Collected here so the reasoning is not spread across
six comments:

- **`BoundarySponge`** on `closure` — the open boundary radiated noise nothing absorbed. See
  "`closure` and `BoundarySponge`" under Simulations.
- **`tracer_advection = WENO()`**, a scalar rather than `(T = WENO(), S = WENO())` — Oceananigans
  gives any tracer a `NamedTuple` omits the `Centered()` default, and CATKE contributes an `e` the
  setup never names, so `e` was being advected by an unbounded centered scheme. A scalar covers
  whatever the closure adds, now and later. This is the shape NumericalEarth's own `ocean_simulation`
  uses.
- **`inflow_timescale = 3hours`**, down from `1day` — at the ~10 s steps the run actually takes, a
  one-day timescale relaxes by 1.2 × 10⁻⁴ per step against an advective rate into the boundary cell
  of ~2.5 × 10⁻³ s⁻¹. Two hundred times weaker than what it was competing with. Both timescales are
  tuning knobs; too strong an inflow nudge over-constrains an open boundary and reflects.
- **`max_slope_factor = 0.25`**, down from `0.5` — measured on this bathymetry, the Drøbak sill goes
  from 65.2 m to 64.6 m, the deepest point is untouched and water volume changes by 0.000 %, while
  adjacent pairs steeper than r = 0.3 go from 4620 to none. The standing worry that a tight limit
  flattens a genuine sill does not survive contact with this domain.
- **`max_island_cells = 6`, unchanged and deliberately so** — every remaining interior land component
  is 7 cells or larger with a minimum width of 2 to 4 cells and bounding boxes like 4×8 and 5×4.
  Those are compact skerries, not the one-cell ridges the stage exists to remove, and raising the
  threshold would delete real topography.

`start_date` and `stop_time` also decide what the prepare steps write, since both pad their time
axes to that window — so changing either is a data change, not just a run change. See "Changing
`start_date` or `stop_time`" under Commands.

Each step is a `FjordConfig` method on the generic function of the same name —
`prepare_bathymetry`, `download_forcing`, `prepare_forcing`, `add_rivers`, `download_boundaries`,
`prepare_boundaries`, `download_atmosphere`,
`prepare_atmosphere`, `run_simulation` — living beside the pipeline it drives. Each builds the
grid, checks the step before it has run, calls the generic pipeline, plots, and logs where the
output went. A step the setup opts out of returns `nothing` rather than raising, matching the
`::Nothing` methods of the lower arities; "you asked for a step this setup does not configure" is
reported by `CLI.main`, because that is user input rather than a pipeline condition.

`examples/oslofjorden.jl` is the worked example of an out-of-tree config: a variant of
`src/Setups/oslofjorden.jl` running an implicit free surface and a radiating open lateral boundary,
which it gets by subtyping `AbstractFreeSurfaceConfig`, `AbstractBoundaryConditionConfig` and
`AbstractGridConfig` in the file itself. It is not a runner — it is passed as
`--config examples/oslofjorden.jl`, and it shares `oslofjorden()`'s `data_root` so the atmosphere
prepare steps do not have to run twice.

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
  fine. Every `Oceananigans.Units` constant is already a `Float64`, so writing durations with them —
  `365days`, `1hour`, `3minutes` — sidesteps this. `loops` is the one field this does *not* apply to,
  being an `Int` already. `SimulationConfig` is now the **only** config the rule reaches: every other
  one is either non-parametric, or — like `SnapshotWriter`, `ProgressCallback` and
  `CoupledHydrostaticSimulation`, all of which *are* parametric — has a hand-written keyword
  constructor that converts, precisely so a reader who learned the rule from `SimulationConfig` is
  not misled about the types they write far more often. `CoupledHydrostaticSimulation` gained one
  when `extra_kwargs` arrived, since that field needs validating at construction anyway.
  `Forcing.PreparedVariable` is the one internal type the rule reaches: it became parametric in its
  dimension-name count when the boundary pipeline started reusing it, so its hand-written constructor
  is what converts the `BitArray` `water_mask` returns into the declared `Array{Bool,3}` — the
  generated constructor would have demanded an exact match.

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
