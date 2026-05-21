# FjordSim Public API Catalog

**Last Updated:** May 6, 2026  
**Purpose:** Maintain a precise record of all exported interfaces (functions, types, and their signatures) to ensure backward compatibility and catch breaking changes before they reach users.

---

## Module Exports by Location

### `FjordSim` (main module, [src/FjordSim.jl](../../src/FjordSim.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `ImmersedBoundaryGrid` | Constructor (Oceananigans overload) | `ImmersedBoundaryGrid(filepath::String, arch, halo)` | Create a bathymetry-based immersed boundary grid from NetCDF |
| `forcing_from_file` | Function | See [Forcing module](#forcing-module) | Re-export from Forcing |
| `top_bottom_boundary_conditions` | Function | See [BoundaryConditions module](#boundaryconditions-module) | Re-export from BoundaryConditions |
| `coupled_hydrostatic_simulation` | Function | `coupled_hydrostatic_simulation(grid, buoyancy, closure, tracer_advection, momentum_advection, tracers, initial_conditions, free_surface, coriolis, forcing, boundary_conditions, atmosphere, downwelling_radiation, sea_ice, biogeochemistry; results_dir=homedir()/FjordSim_results, stop_time=365days)` | Create and initialize a hydrostatic coupled simulation |
| `recursive_merge` | Function | See [Utils module](#utils-module) | Re-export from Utils |
| `progress` | Function | See [Utils module](#utils-module) | Re-export from Utils |
| `cell_advection_timescale_coupled_model` | Function | See [Utils module](#utils-module) | Re-export from Utils |
| `NORA3PrescribedAtmosphere` | Type/Struct | See [Atmospheres module](#atmospheres-module) | Re-export from Atmospheres |
| `NORA3PrescribedRadiation` | Function | See [Atmospheres module](#atmospheres-module) | Re-export from Atmospheres |
| `MultiYearNORA3` | Type/Struct | See [Atmospheres module](#atmospheres-module) | Re-export from Atmospheres |

---

### `BoundaryConditions` module ([src/BoundaryConditions.jl](../../src/BoundaryConditions.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `top_bottom_boundary_conditions` | Function | `top_bottom_boundary_conditions(; grid, bottom_drag_coefficient) -> NamedTuple` | Create boundary condition tuple with fields: `u(top, bottom)`, `v(top, bottom)`, `T(top)`, `S(top)` |

**Return Type Details:**
```julia
(
    u = (top = FluxBoundaryCondition(τx), bottom = FluxBoundaryCondition),
    v = (top = FluxBoundaryCondition(τy), bottom = FluxBoundaryCondition),
    T = (top = FluxBoundaryCondition(Jᵀ),),
    S = (top = FluxBoundaryCondition(Jˢ),)
)
```

---

### `Forcing` module ([src/Forcing.jl](../../src/Forcing.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `forcing_from_file` | Function | `forcing_from_file(...)` | *Signature not yet extracted; requires further inspection* |
| `ForcingFromFile` | Struct | `ForcingFromFile{FTS,V} where {FTS,V}` with fields: `fts_value::FTS`, `fts_λ::FTS`, `fieldname::V` | Internal callable struct for in-kernel forcing application |

---

### `FDatasets` module ([src/FDatasets.jl](../../src/FDatasets.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `DSForcing` | Struct | `DSForcing(metadata_filename::String, default_download_directory::String)` | Metadata wrapper for downscaled forcing datasets |
| `DSResults` | Struct | `DSResults(metadata_filename::String, default_download_directory::String; start_date_time::DateTime)` | Metadata wrapper for downscaled simulation results |
| `last_date` | Function | Re-export from NumericalEarth | Query last available date in dataset |

**DSForcing Fields:**
- `metadata_filename::String`
- `default_download_directory::String`
- `reversed_vertical_axis::Bool`
- `longitude_interfaces::Tuple`
- `latitude_interfaces::Tuple`
- `size::Tuple` (Nx, Ny, Nz)
- `all_dates::Vector`
- `first_date`
- `last_date`
- `z_interfaces::Vector`

**DSResults Fields:** Same as DSForcing (same struct definition reused for results)

---

### `Utils` module ([src/Utils.jl](../../src/Utils.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `compute_faces` | Function | `compute_faces(centers::Vector) -> Vector` | Convert grid centers to face locations |
| `progress` | Function | `progress(sim) -> nothing` | Print iteration progress with velocity and temperature info |
| `safe_execute` | Function | `safe_execute(callable) -> Function` | Return a safe wrapper that handles `nothing` inputs |
| `extract_z_faces` | Function | `extract_z_faces(grid) -> Vector` | Extract positive z-face coordinates from grid |
| `netcdf_to_jld2` | Function | `netcdf_to_jld2(netcdf_file::String, jld2_file::String) -> nothing` | Convert NetCDF file to JLD2 format |
| `save_fts` | Function | `save_fts(; jld2_filepath, fts_name, fts, grid, times, boundary_conditions) -> ...` | Save a FieldTimeSeries to disk |
| `recursive_merge` | Function | `recursive_merge(dict1, dict2) -> Dict` | *Signature not yet extracted; recursive merge of dicts* |
| `cell_advection_timescale_coupled_model` | Function | `cell_advection_timescale_coupled_model(coupled_model) -> Float` | Compute advection timescale from coupled model |

---

### `Grids` module ([src/Grids.jl](../../src/Grids.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `ImmersedBoundaryGrid` | Constructor | `ImmersedBoundaryGrid(filepath::String, arch, halo) -> ImmersedBoundaryGrid` | Create immersed boundary grid from bathymetry NetCDF |

---

### `Bathymetry` module ([src/Bathymetry.jl](../../src/Bathymetry.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `GeonorgeBathymetry` | Struct | `GeonorgeBathymetry(metadata_filename::String, default_download_directory::String, longitude_interfaces::NTuple{2,Float64}, latitude_interfaces::NTuple{2,Float64}, size::NTuple{3,Int})` | Metadata wrapper for generated raw Geonorge bathymetry |
| `prepare_geonorge_bathymetry` | Function | `prepare_geonorge_bathymetry(target_grid; output_path, geodatabase_path, raw_dir=download_bathymetry_cache, raw_resolution_factor=4, padding_cells=0, include_contours=true, contour_stride=1, geonorge_cache=true, regrid_cache=true, regrid_kw...) -> NamedTuple` | Build a regional raw bathymetry dataset from the local Geonorge FileGDB, regrid it, and write a FjordSim bathymetry NetCDF |
| `write_bathymetry_file` | Function | `write_bathymetry_file(filepath::String, target_grid, bottom_height) -> String` | Write a processed NetCDF bathymetry file compatible with `ImmersedBoundaryGrid` |

---

### `Atmospheres` module ([src/Atmospheres/Atmospheres.jl](../../src/Atmospheres/Atmospheres.jl))

| Symbol | Type | Signature | Purpose |
|--------|------|-----------|---------|
| `NORA3PrescribedAtmosphere` | Function | `NORA3PrescribedAtmosphere([architecture=CPU(), FT=Float32]; dataset, start_date, end_date, backend, time_indexing, surface_layer_height, other_kw...) -> PrescribedAtmosphere` | Build prescribed atmosphere from NORA3 reanalysis data |
| `NORA3PrescribedRadiation` | Function | `NORA3PrescribedRadiation([architecture=CPU(), FT=Float32]; dataset, start_date, end_date, backend, time_indexing, ocean_surface, sea_ice_surface, stefan_boltzmann_constant, other_kw...) -> PrescribedRadiation` | Build prescribed radiation (shortwave + longwave) from NORA3 data |
| `MultiYearNORA3` | Struct | `MultiYearNORA3(metadata_filename::String, default_download_directory::String)` | Metadata container for NORA3 atmosphere data |

**MultiYearNORA3 Fields:**
- `metadata_filename::String`
- `default_download_directory::String`
- `size::Tuple` (Nx, Ny) — 2D grid size
- `all_dates::Vector{<:Dates.Date}` — all available time steps

**NORA3 Variables (internal reference):**  
`(:freshwater_flux, :specific_humidity, :sea_level_pressure, :downwelling_longwave_radiation, :downwelling_shortwave_radiation, :temperature, :eastward_velocity, :northward_velocity)`

---

## Backward Compatibility Commitments

### Hard Constraints

1. **No removal of exported symbols** — once added to an `export` list, a symbol stays exported or becomes a deprecated alias.

2. **Function signatures must be preserved:**
   - Positional argument order, count, and types must not change.
   - Keyword argument names must not change; only new optional keywords (with defaults) may be added.
   - Return type/structure must remain identical in shape.

3. **Type/Struct fields must be preserved:**
   - Existing fields on exported structs must not be removed or renamed.
   - New fields may be added if they don't break construction (e.g., added with defaults or after construction).

4. **Overloaded methods (e.g., Oceananigans types) must preserve dispatch:**
   - Do not narrow the method signature (e.g., from `AbstractGrid` to `RectilinearGrid`).

### Process for Adding New Functions/Types

When adding a new exported entity:
1. Add the symbol to the `export` list in its module.
2. **Update this catalog immediately** with the full signature and return type.
3. Add test coverage if possible.
4. Document any non-obvious constraints in the docstring.

### Process for Changing Exported Interfaces

If you must change an existing exported interface:
1. **Verify with team** — breaking changes need explicit approval.
2. **Document the change here** in the changelog section below with:
   - Old and new signatures side-by-side
   - Migration path for users
   - Deprecation timeline (if applicable)
3. **Implement a deprecation wrapper** if removing old behavior; do not silently break.

---

## Maintenance Checklist

When preparing a release or PR that touches exported code:

- [ ] All new exported symbols added to this catalog with full signature
- [ ] All modified exported signatures updated here
- [ ] If removing an export, deprecation path documented
- [ ] Return types and argument order verified against catalog
- [ ] Tests added/updated for new exports or changed behavior
- [ ] Changelog entry added (if breaking)

---

## Related Files

- [Backward Compatibility Instructions](./backward-compatibility.instructions.md) — Principles and rules
- [Copilot Instructions](../copilot-instructions.md) — Project conventions and style
- [AGENTS.md](../../AGENTS.md) — Code simplicity and surgical changes guidelines
