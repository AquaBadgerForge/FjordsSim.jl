# Config methods reference

Every FjordSim pipeline is one generic function plus a handful of *hooks* — dispatch points that take
a config as an argument. A new source is therefore purely additive: define a struct subtyping the
matching abstract config type, add methods for that struct to the hook functions, and the generic
pipelines, the plots, the path helpers and the command line all work on it without being edited.

This file is the complete inventory. Each supertype gets three tables:

- **Required hooks** — no method exists on the supertype. A missing one surfaces as a `MethodError`
  naming the function.
- **Optional hooks** — the supertype supplies a default, given in the table. Override only to change
  it.
- **Inherited generics** — defined on the supertype, so a subtype gets them for free. These are what
  you get in exchange for implementing the hooks, and they are the reason a new source needs so
  little code.

The `Built-in` column points at the method to copy. Line numbers are accurate as of this writing;
`grep` the name if one has drifted.

See the "Extending FjordSim" section of [`../README.md`](../README.md) for the short version, and each
supertype's docstring in [`../src/Configs.jl`](../src/Configs.jl) for the field contract that goes
alongside these methods.

---

## Contents

| Supertype | Defined at | Built-in subtype | Required hooks |
|---|---|---|---|
| [`AbstractGridConfig`](#abstractgridconfig) | `src/Configs.jl:44` | `EvenGrid` | 1 |
| [`AbstractBathymetryConfig`](#abstractbathymetryconfig) | `src/Configs.jl:64` | `DybdedataConfig` | 1 |
| [`AbstractForcingConfig`](#abstractforcingconfig) | `src/Configs.jl:93` | `NorKystConfig` | 3 (+1 if it downloads) |
| [`AbstractRiverConfig`](#abstractriverconfig) | `src/Configs.jl:121` | `OF800RiversConfig` | 2 (+1 if it downloads) |
| [`AbstractAtmosphereConfig`](#abstractatmosphereconfig) | `src/Configs.jl:164` | `NORA3Config` | 3 (+1 if it downloads, +2 if simulated) |
| [`AbstractSimulationConfig`](#abstractsimulationconfig) | `src/Configs.jl:196` | `SimulationConfig` | 0 — fields only |
| [`FjordConfig`](#fjordconfig) | `src/Configs.jl:340` | — (not a supertype) | — |

---

## `AbstractGridConfig`

Built-in subtype: `EvenGrid` (`src/Grids.jl:27`).

### Required hooks

| Hook | Built-in |
|---|---|
| `LatitudeLongitudeGrid(architecture, config)` | `src/Grids.jl:35` |

### Optional hooks

None.

### Inherited generics

**None.** This is the one supertype with no generic surface at all — no method anywhere dispatches on
`AbstractGridConfig` itself, and it declares no required fields. It is the cheapest supertype to
extend and, for the same reason, the least discoverable from `src/Configs.jl`.

Two things follow. The single hook extends a **foreign** generic — `LatitudeLongitudeGrid` is
imported from Oceananigans (`src/Grids.jl:10`), not declared here — so there is no local stub to
`grep` for and no FjordSim-owned function name to guide you. And the only thing anything reads off a
grid config directly is the field `config.grid_config.halo` (`src/Forcing/Forcing.jl:602`,
`src/Forcing/rivers.jl:201`); everything else goes through Oceananigans' `x_domain`/`y_domain` on the
*built grid*. So a grid config needs `halo` plus whatever its own constructor uses.

---

## `AbstractBathymetryConfig`

Built-in subtype: `DybdedataConfig` (`src/Bathymetry/geonorge.jl:101`), the Geonorge Sjøkart
Dybdedata source. Template file: `src/Bathymetry/geonorge.jl`.

### Required hooks

| Hook | Stub | Built-in |
|---|---|---|
| `bathymetry_dataset(target_grid, config)` → a NumericalEarth dataset | `src/Bathymetry/Bathymetry.jl:51` | `src/Bathymetry/geonorge.jl:133` |

### Optional hooks

| Hook | Default | Default at | Built-in |
|---|---|---|---|
| `regrid_options(config)` → `NamedTuple` for `regrid_bathymetry` | `(;)` | `src/Bathymetry/Bathymetry.jl:60` | `src/Bathymetry/geonorge.jl:146` |
| `smoothing_options(config)` → `NamedTuple` for `smooth_bathymetry_gaps!` | `(;)` | `src/Bathymetry/Bathymetry.jl:70` | `src/Bathymetry/geonorge.jl:161` |

### Inherited generics

| Function | Defined at |
|---|---|
| `prepare_bathymetry(target_grid, config; regrid_kw...)` | `src/Bathymetry/Bathymetry.jl:87` |
| `bathymetry_path(config)` | `src/Configs.jl:209` |
| `plot_path(config)` | `src/Configs.jl:212` |
| `plot_bathymetry(grid, bottom_height, config; title, figure_size)` | `src/Plotting.jl:97` |

---

## `AbstractForcingConfig`

Built-in subtype: `NorKystConfig{R}` (`src/Forcing/norkyst.jl:55`), NorKyst-800m over THREDDS.
Template file: `src/Forcing/norkyst.jl`.

`NorKystConfig` is parametric solely so its `rivers` field stays concretely typed whether it holds a
river config or `nothing`. A consequence worth knowing: **`rivers` is a construction-time choice**
and cannot be unset on an existing instance.

### Required hooks

| Hook | Stub | Built-in |
|---|---|---|
| `forcing_time_steps(config)` → `Vector{SourceRecord}` | `src/Forcing/Forcing.jl:343` | `src/Forcing/norkyst.jl:92` |
| `forcing_source_grid(config, filepath)` → source grid | `src/Forcing/Forcing.jl:354` | `src/Forcing/norkyst.jl:123` |
| `forcing_variable_names(config)` → `Dict` source name => FjordSim name | `src/Forcing/Forcing.jl:365` | `src/Forcing/norkyst.jl:84` |
| `download_forcing(target_grid, config)` — only if it downloads | **none** (see below) | `src/Forcing/norkyst.jl:149` |

`download_forcing(target_grid, config)` has **neither a supertype fallback nor a stub**, so a source
that omits it fails with a raw `MethodError` from the `FjordConfig` driver rather than a stated
error. A source that needs no download must define a no-op.

### Optional hooks

None on the config. A source grid that is *not* a regular projected grid additionally needs
`source_field_grid(source, architecture)` and
`projected_target_nodes(longitude, latitude, source)` — those dispatch on the **source-grid type**,
not on the config.

### Inherited generics

| Function | Defined at |
|---|---|
| `prepare_forcing(target_grid, config; coverage)` | `src/Forcing/Forcing.jl:551` |
| `prepared_variable(source_name, target_grid, source, filepath, config)` | `src/Forcing/Forcing.jl:822` |
| `relaxation_lambda(mask, config)` | `src/Forcing/Forcing.jl:1145` |
| `forcing_from_file(config; grid, tracers, reference_date)` | `src/Forcing/Forcing.jl:289` |
| `interpolation_architecture(config)` | `src/Forcing/Forcing.jl:314` |
| `add_rivers(target_grid, config)` → dispatches on `config.rivers` | `src/Forcing/rivers.jl:212` |
| `simulation_forcing_path(config, rivers)` | `src/Simulations.jl:376-377` |
| `forcing_path(config)` | `src/Configs.jl:210` |
| `forcing_directory(config)` | `src/Configs.jl:223` |
| `plot_path(config)` | `src/Configs.jl:213` |
| `plot_forcing(grid, config)` | `src/Plotting.jl:168` |

---

## `AbstractRiverConfig`

Built-in subtype: `OF800RiversConfig` (`src/Forcing/of800_rivers.jl:40`). Template file:
`src/Forcing/of800_rivers.jl`.

A river config is **not** a `FjordConfig` field. It goes in the forcing config's `rivers` field, where
`nothing` means the setup has no rivers and `add_rivers` is a no-op.

### Required hooks

| Hook | Stub | Built-in |
|---|---|---|
| `river_locations(config)` → `Vector{RiverLocation}` | `src/Forcing/rivers.jl:39` | `src/Forcing/of800_rivers.jl:114` |
| `river_series(config, times)` → `Dict` FjordSim name => `(river, time)` matrix | `src/Forcing/rivers.jl:48` | `src/Forcing/of800_rivers.jl:149` |
| `download_rivers(config)` | `src/Forcing/rivers.jl:55` | `src/Forcing/of800_rivers.jl:69` |

`download_rivers` has a bare stub but **is not optional**: `add_rivers(config::FjordConfig)` calls it
unconditionally (`src/Forcing/rivers.jl:198`), so a source that ships its data with the setup still
has to define a no-op.

### Optional hooks

| Hook | Default | Default at | Built-in |
|---|---|---|---|
| `river_search_radius(config)` → cells to search for a coastal cell | `config.search_radius` | `src/Forcing/rivers.jl:63` | not overridden — uses the default |

### Inherited generics

| Function | Defined at |
|---|---|
| `add_rivers(target_grid, forcing_config, rivers)` | `src/Forcing/rivers.jl:216` |
| `river_forcing_path(config)` | `src/Configs.jl:232` |
| `simulation_forcing_path(forcing_config, rivers)` | `src/Simulations.jl:377` |
| `forcing_prerequisite(rivers)` → `"add_rivers"` | `src/Simulations.jl:386` |

The last two are why naming rivers changes the *run*, not just the rivers step: the simulation reads
`river_forcing_path` instead of `forcing_path`, so `add_rivers` becomes a prerequisite of
`run_simulation`.

A variable `river_series` returns that the forcing file does not carry is skipped with a warning, so
one river dataset can serve setups that prepare different variables. An outlet outside the grid, or
with no coastal water cell within the search radius, is dropped with a message rather than written
into land; `add_rivers` fails only when *no* outlet lands.

---

## `AbstractAtmosphereConfig`

Built-in subtype: `NORA3Config` (`src/Atmospheres/nora3_source.jl:91`), MET Norway NORA3. Template
file: `src/Atmospheres/nora3_source.jl`.

### Required hooks

| Hook | Stub | Built-in |
|---|---|---|
| `atmosphere_time_steps(config)` → `Vector{AtmosphereRecord}` | `src/Atmospheres/Atmospheres.jl:128` | `src/Atmospheres/nora3_source.jl:143` |
| `atmosphere_source_grid(config, filepath)` → source grid | `src/Atmospheres/Atmospheres.jl:138` | `src/Atmospheres/nora3_source.jl:178` |
| `atmosphere_variable_names(config)` → `Dict` downloaded name => prepared name | `src/Atmospheres/Atmospheres.jl:150` | `src/Atmospheres/nora3_source.jl:135` |
| `download_atmosphere(target_grid, config)` — only if it downloads | `::Nothing` only, `src/Atmospheres/Atmospheres.jl:172` | `src/Atmospheres/nora3_source.jl:247` |
| `prescribed_atmosphere(config, architecture; reference_date)` — only if simulated | `::Nothing` only, `src/Atmospheres/Atmospheres.jl:197` | `src/Atmospheres/nora3_source.jl:216` |
| `prescribed_radiation(config, architecture; reference_date)` — only if simulated | `::Nothing` only, `src/Atmospheres/Atmospheres.jl:198` | `src/Atmospheres/nora3_source.jl:219` |

The last three have **only** a `::Nothing` method — no `::AbstractAtmosphereConfig` fallback and no
stub — so an unimplemented subtype fails with a raw `MethodError` rather than a stated error. Those
`::Nothing` methods exist to make a setup naming *no* atmosphere a no-op, not to serve as defaults.

`prescribed_atmosphere` and `prescribed_radiation` are the read side, consumed by `build_simulation`
rather than by `prepare_atmosphere`. They are what keep `src/Simulations.jl` from naming any dataset.
Neither takes a float type on purpose: both NumericalEarth constructors default to `Float32`, and
passing `Oceananigans.defaults.FloatType` would silently promote the atmosphere to `Float64`.

`reference_date` is the instant the returned time axes are zeroed at — *not* the same thing as
`NORA3PrescribedAtmosphere`'s `start_date`, which selects which records to load. `build_simulation`
passes the simulation config's `start_date` to both these hooks and to `forcing_from_file`, which is
the only thing keeping the atmosphere and the forcing in phase.

### Optional hooks

| Hook | Default | Default at | Built-in |
|---|---|---|---|
| `atmosphere_date_range(config)` → `(first, last)` `DateTime`s | `nothing`, skipping the coverage check | `src/Atmospheres/Atmospheres.jl:212` | `src/Atmospheres/nora3_source.jl:228` |

### Inherited generics

| Function | Defined at |
|---|---|
| `prepare_atmosphere(target_grid, config; coverage)` | `src/Atmospheres/Atmospheres.jl:439` |
| `atmosphere_target_axes(target_grid, config)` | `src/Atmospheres/Atmospheres.jl:231` |
| `atmosphere_path(config)` | `src/Configs.jl:211` |
| `atmosphere_directory(config)` | `src/Configs.jl:241` |
| `plot_path(config)` | `src/Configs.jl:214` |
| `plot_atmosphere(config)` | `src/Plotting.jl:226` |
| `MultiYearNORA3(config)` — the reader for a prepared file | `src/Atmospheres/NORA3.jl:99` |

The prepared variable names and units are **not** a hook. They are fixed by the read side in
`ATMOSPHERE_VARIABLES`: eight `Float32` variables of shape `(lon, lat, time)`, air temperature in
Kelvin, both radiative fluxes **downwelling**. A source whose download already normalizes names
returns the identity mapping from `atmosphere_variable_names`.

---

## `AbstractSimulationConfig`

Built-in subtype: `SimulationConfig` (`src/Simulations.jl:98`). Template: the `SimulationConfig` block
in `src/Setups/oslofjorden.jl`.

### Required hooks

**None.** This supertype is data only, read field by field, because everything dataset-specific
already comes from the other configs. An alternative simulation config is a subtype supplying the
field set its docstring lists — `results_root`, `output_file`, `start_date` and `stop_time` are the
supertype-level minimum. No method anywhere dispatches on the concrete `SimulationConfig`.

### Optional hooks

| Hook | Default | Default at |
|---|---|---|
| `run_tag(config)` → the run's identity as a filename fragment | `LAUNCH_TAG[]`, the wall-clock instant the process started | `src/Configs.jl:258` |

`run_tag` is the one exception to "no hooks". It takes a config it does not read, precisely so a
subtype can name its runs differently.

### Inherited generics

| Function | Defined at |
|---|---|
| `results_path(config)` / `results_path(config, loop)` | `src/Configs.jl:272`, `:274` |
| `coverage_window(config)` → the interval the prepared inputs must span | `src/Configs.jl:307` |
| `simulation_architecture(config)` | `src/Simulations.jl:135` |
| `loop_output_path(config, loop)` | `src/Simulations.jl:433` |
| `checkpointed_loops(config)` | `src/Simulations.jl:458` |
| `resume_loop(config)` | `src/Simulations.jl:480` |
| `attach_writers!(simulation, config, loop)` | `src/Simulations.jl:509` |
| `attach_checkpointer!(simulation, config, loop)` | `src/Simulations.jl:547` |
| `restart_loop!(simulation, config, loop)` | `src/Simulations.jl:607` |
| `log_path(config)` | `src/CLI.jl:102` |

`coverage_window` reads `start_date`, so a subtype omitting that field inherits only the other two
path helpers.

---

## `FjordConfig`

`FjordConfig` (`src/Configs.jl:340`) is not a supertype — it is the concrete container holding one of
each config, parametrically so every instantiation stays concretely typed. `atmosphere_config` and
`simulation_config` default to `nothing`; a river config hangs off the forcing config's `rivers`
field instead.

These are the setup-level drivers. Each builds the grid, checks that the previous step has run, calls
the generic pipeline, plots, and logs where the output went. **None is ever overloaded per source**,
and each maps one-to-one onto a CLI subcommand (`src/CLI.jl:18-26`), in this order:

| Driver | Subcommand | Defined at |
|---|---|---|
| `prepare_bathymetry(config)` | `prepare_bathymetry` | `src/Bathymetry/Bathymetry.jl:119` |
| `download_forcing(config)` | `download_forcing` | `src/Forcing/Forcing.jl:379` |
| `prepare_forcing(config)` | `prepare_forcing` | `src/Forcing/Forcing.jl:595` |
| `add_rivers(config)` | `add_rivers` | `src/Forcing/rivers.jl:188` |
| `download_atmosphere(config)` | `download_atmosphere` | `src/Atmospheres/Atmospheres.jl:169` |
| `prepare_atmosphere(config)` | `prepare_atmosphere` | `src/Atmospheres/Atmospheres.jl:487` |
| `run_simulation(config)` | `run_simulation` | `src/Simulations.jl:787` |

Two `FjordConfig` methods have no subcommand:

| Function | Defined at | What it is for |
|---|---|---|
| `build_simulation(config)` | `src/Simulations.jl:648` | returns the instrumented `Simulation` without running it — the REPL and debugger entry point |
| `simulation_forcing_path(config)` | `src/Simulations.jl:373` | which prepared forcing file the run reads, dispatching on `config.forcing_config.rivers` |

A step the setup opts out of returns `nothing` rather than raising. "You asked for a step this setup
does not configure" is user input, so it is reported by `CLI.main`, not by a pipeline.

---

## The shape of a new source

Concretely, adding a forcing dataset is one file:

```julia
Base.@kwdef mutable struct MyForcingConfig{R} <: AbstractForcingConfig
    data_root::String                      # the fields the supertype docstring lists
    output_directory::String
    output_file::String = "forcing.nc"
    plot_file::String = "forcing.png"
    relaxation_edge::Symbol = :south
    relaxation_cells::Int = 10
    relaxation_timescale::Float64 = 86400.0
    architecture::Symbol = :auto
    parameters::Vector{String}
    rivers::R = nothing                    # parametric so the field stays concrete
end

# The three required hooks, plus download_forcing since this one downloads.
Forcing.forcing_variable_names(config::MyForcingConfig) = Dict("thetao" => "T", "so" => "S")
Forcing.forcing_time_steps(config::MyForcingConfig) = ...
Forcing.forcing_source_grid(config::MyForcingConfig, filepath) = ...
Forcing.download_forcing(target_grid, config::MyForcingConfig) = ...
```

That is the whole extension. `prepare_forcing`, `add_rivers`, `forcing_from_file`, `plot_forcing`,
`forcing_path`, `forcing_directory`, `relaxation_lambda`, the `--config` handling and all seven
subcommands now work on it, because they dispatch on `AbstractForcingConfig` and never on the
concrete type.

Two conventions to keep (see the "Import Conventions" section of `CLAUDE.md`): extend a function as
`function Mod.foo(...)` or `Mod.foo(...) = ...`, never `import Mod: foo`; and keep every field
concretely typed or parameterized — no `Any`, no abstract field types.
