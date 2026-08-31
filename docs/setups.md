# Setups

> Split out of `CLAUDE.md`. Read this when adding or tuning a fjord setup — it has the registry
> mechanics, the vertical-grid rationale, and the numbers behind `oslofjorden()`'s current
> parameter values.

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
`drammensfjorden()` gives its `OF800RiversConfig`, and `oslofjorden()` its `NVERiversConfig`, the
same `data_root` as the rest of the setup — the river data downloads there rather than being shared from elsewhere, so
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

## The vertical grid

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
`snap_partial_bottom_cells` in `docs/architecture.md` (Bathymetry section), which is what handles
them. A 28-level ratio-1.20 variant scores better still (median floor layer 8.4 m, *no* column
floored in a layer over 30 m) at 56 % more cells than the original rather than 33 %; it is the option
to reach for if the deep basins still look under-resolved.

**Changing it is a data change.** `z_faces` is written into `bathymetry.nc`, and `simulation_grid`
reads the file rather than the config — so the run silently uses whatever the file holds. Every 3D
prepared file is regridded onto that grid, so re-run `prepare_bathymetry` → `prepare_forcing` →
`add_rivers` → `prepare_boundaries`. No download step is affected. The atmosphere is 2D and is not.

## The 2020 diagnosis

Several of `oslofjorden()`'s current values are there because of one run, and are commented in the
setup file with the number that justifies them. Collected here so the reasoning is not spread across
six comments:

- **`BoundarySponge`** on `closure` — the open boundary radiated noise nothing absorbed. See
  "`closure` and `BoundarySponge`" in `docs/architecture.md` (Simulations section).
- **`HorizontalScalarBiharmonicDiffusivity(ν = 2e4, κ = 2e3)`**, down from `1e5`/`1e4` — a
  biharmonic coefficient only means something against `Δx⁴`, and the commit that raised ν from 15
  justified 1e5 as a "~1.7 hours" e-folding of the 2Δx mode when the discrete rate `ν₄·16/Δx⁴` makes
  it 14.5 min on this 193 m cell, 7x stronger than intended. 2e4 restores the intended 72 min, still
  ~1300x the value it replaced, and leaves an 8 km baroclinic eddy alone (56 h at 8Δx). Norkyst-800,
  the ROMS system feeding this run, applies **no** explicit interior viscosity and a 10 m² s⁻¹
  harmonic tracer diffusivity, against 0.2 m² s⁻¹ at 2Δx for `κ = 2e3` here. It also restores a
  stability margin nothing measures: explicit biharmonic diffusion needs `Δt ≤ Δx⁴/32ν₄`, 434 s at
  1e5 and 2170 s at 2e4, while `AdaptiveTimeStep` only ever sees the advective CFL. The same
  correction is larger on `drammensfjorden()`, whose 94 m cell had a 24 s diffusive limit at 1e5
  against ~10 s steps; its ν is 1e3, the Δx⁴-scaled equivalent, not a copy of Oslofjord's.
- **`minimum_tke = 7e-6` kept, not lowered** — CATKE takes `w★ = sqrt(max(minimum_tke, e))`, and
  measured on the 2020 run the prognostic `e` is below that floor in 85 % of wet cells (95-100 %
  below 9 m), so the floor and not the TKE equation sets the vertical diffusivity, at the documented
  background `κ = 0.098·e_min/N`, `ν = 0.242·e_min/N`. The domain-median `κᵤ/κᶜ` of 2.38 against
  `Cʰⁱᵤ/Cʰⁱᶜ = 2.47` is the confirmation. The value is ROMS' `GLS_KMIN` default (7.6e-6) and does not
  transfer cleanly — ROMS floors its length scale too (`GLS_PMIN`) where CATKE's is diagnostic — but
  it lands in the right place for this fjord: median `κᶜ` at 94-118 m is 8-9e-5 m² s⁻¹ against
  basin-mean density-budget values of 1.0-1.8e-4 (Bunnefjorden), 4.9-7.6e-4 (Vestfjorden) and 2.5e-3
  just inside the Drøbak sill (Gade 1970; Staalstrøm et al. 2012, *Ocean Sci.* 8, 525). The basin
  water is at the low end already, so lowering the floor would slow deep-water renewal, which is what
  this domain is judged on. `Cᵇ` stays at Oceananigans' 0.28 rather than NumericalEarth's regional
  0.01 for the same reason: it scales the near-bottom mixing length, and the high basin-mean
  diffusivity here is attributed to exactly that boundary mixing.
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
- **`NVERiversConfig` in place of `OF800RiversConfig`, with `minimum_levels = 0`** — the shallow-outlet
  runaway that killed the run at day 11.5 is now handled by filling the column rather than by
  relocating the outlet, so the two mechanisms are not stacked. `river_lambdas` scales λ by each
  river's own mean discharge, where OF800 gave all nineteen outlets the same λ, and Glomma and
  Drammenselva take `plume_depth = Inf` because their cells are river bed.
- **`minimum_discharge = 0.5`, and no stated coordinate anywhere** — the nine outlets this setup
  first carried were OF800's, and several were kilometres from their rivers; see "Outlet discovery"
  in `docs/architecture.md` (Forcing section) for the measurements. Discovery gives 21 mouths,
  including both of Glomma's, sized from their own REGINE catchments. Eight overrides add the gauges
  that survive inspection; the other thirteen mouths are freshwater-only, which is what
  `river_lambdas`' catchment-normal rule is for. Only four of the 21 carry water temperature —
  Akerselva's gauge reads 21.2 °C mean, Lysakerelva's 19.6 °C mean and 37.2 °C peak, and
  Gjersjøelva's 2020 series is empty.
- **Glomma is two rivers here**, `002.A21` (Østerelva) and `002.2A` (Vesterelva), both on
  Solbergfoss's series at `discharge_fraction` 2/3 and 1/3. Both beds are resolved on this grid —
  cells (222, 99) and (198, 101) — and the split is a **stated assumption**: NVE's own Nedre Glomma
  flood report describes the division around Kråkerøy and publishes no fraction, and none was found
  elsewhere. It is one number in the setup file to change when one turns up.

`start_date` and `stop_time` also decide what the prepare steps write, since both pad their time
axes to that window — so changing either is a data change, not just a run change. See "Changing
`start_date` or `stop_time`" in `CLAUDE.md`.

## Per-step drivers and the example config

Each step is a `FjordConfig` method on the generic function of the same name —
`prepare_bathymetry`, `download_forcing`, `prepare_forcing`, `add_rivers`, `download_boundaries`,
`prepare_boundaries`, `download_atmosphere`,
`prepare_atmosphere`, `run_simulation` — living beside the pipeline it drives. Each builds the
grid, checks the step before it has run, calls the generic pipeline, plots, and logs where the
output went. A step the setup opts out of returns `nothing` rather than raising, matching the
`::Nothing` methods of the lower arities; "you asked for a step this setup does not configure" is
reported by `CLI.main`, because that is user input rather than a pipeline condition.

`examples/oslofjorden_npzd.jl` is the worked example of an out-of-tree config: a variant of
`src/Setups/oslofjorden.jl` running an implicit free surface and a radiating open lateral boundary,
which it gets by subtyping `AbstractFreeSurfaceConfig`, `AbstractBoundaryConditionConfig` and
`AbstractGridConfig` in the file itself. It is not a runner — it is passed as
`--config examples/oslofjorden_npzd.jl`, and it shares `oslofjorden()`'s `data_root` so the atmosphere
prepare steps do not have to run twice.
