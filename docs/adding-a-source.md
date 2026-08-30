# Adding a new source

> Split out of `CLAUDE.md`. Read this before subtyping any `Abstract*Config` supertype — it lists
> every hook a new source must or may implement, with defaults and rationale.

Every pipeline is a generic function on the config supertype plus a small set of hooks. Add a
source by subtyping and overloading the hooks — never by editing the generic function. The
adapter files (`src/Bathymetry/geonorge.jl`, `src/Forcing/norkyst.jl`, `src/Forcing/of800_rivers.jl`,
`src/Forcing/nve_rivers.jl`) are the templates.

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
| `download_rivers(target_grid, config)` | only if it downloads | none |
| `river_search_radius(config)` → cells to search for a coastal cell | no | `config.search_radius` |
| `river_minimum_levels(config)` → wet levels a column must have to receive an outlet | no | `0`, accepting any water cell. Unlike `river_search_radius` the fallback reads no field, so a river config written before the hook existed keeps working; `OF800RiversConfig` overloads it |
| `river_plume_depth(config, location)` → metres one river's relaxation reaches below the surface | no | `0.0`, the surface level alone. `Inf` asks for the whole wet column. Reads no field, for the same reason as the row above; `NVERiversConfig` overloads it |
| `river_lambdas(config, cells, target_grid)` → `Vector{Float32}`, the coefficient at each cell | no | `1 / config.relaxation_timescale` everywhere. Takes the grid because a coefficient derived from discharge needs the plume volume; `NVERiversConfig` overloads it |

The last two are the extension points for a source whose rivers are *not* interchangeable — a plume
depth per river and a coefficient scaled by river size. Both fallbacks reproduce what the pipeline did
before they existed, so a source that ignores them is unaffected, and `of800_rivers.jl` does ignore
both. See "`river_plume_depth` and `river_lambdas`" in `docs/architecture.md` (Forcing section) for
why a depth rather than a level count, and why the λ cap is a stability requirement rather than a
preference.

Its `standalone` field decides whether `add_rivers` patches a copy of the prepared forcing or writes a
river-only file of its own; see `docs/architecture.md` (Forcing section). A `standalone` config needs
the setup to name a `simulation_config`, since the file's time axis comes from the run window.

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
this hook and `simulation_forcing`, which is the only thing keeping the two in phase — see
`docs/architecture.md` (Simulations section).

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
built-in one by a keyword is a config value rather than a new type and a new method. See
`docs/architecture.md` (Simulations section) for the slot names and for why the constructor rejects
keys that would shadow the fields.

The `coverage` keyword `prepare_forcing`, `prepare_boundaries` and `prepare_atmosphere` take is a
pipeline *argument*, not a hook: the `FjordConfig` drivers derive it with `coverage_window` and every
source inherits the padding unchanged. No source gains a hook for it.

Required fields are listed in each supertype's docstring in `src/Configs.jl`. Path resolution
(`bathymetry_path`, `forcing_path`, `forcing_directory`, `river_forcing_path`, `boundary_data_path`,
`boundary_data_directory`, `atmosphere_path`,
`atmosphere_directory`, `results_path`, `plot_path`) and the plots come for free. A missing
required hook surfaces as a `MethodError` naming it, which the "Config extensibility" testset
asserts.
