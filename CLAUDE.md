# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run all tests
julia --project test/runtests.jl

# Run the Oslofjord example simulation (requires GPU + data files)
julia --project examples/oslofjord.jl

# Prepare bathymetry for a configured fjord (requires Geonorge GDB)
julia --project scripts/bathymetry_prepare.jl --config configs/drammensfjorden.toml
```

Activate the environment before running scripts interactively:
```julia
using Pkg; Pkg.activate(".")
```

## Architecture

FjordSim is a Julia package that wraps [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and [NumericalEarth.jl](https://github.com/CliMA/NumericalEarth.jl) to set up regional ocean simulations of Norwegian fjords.

A simulation is assembled from three components:

1. **Bathymetry file** (`src/Bathymetry.jl`) — processes Geonorge Sjøkart FileGDB data (EPSG:25833) into a NumericalEarth-compatible NetCDF using ArchGDAL. The pipeline: GDB → point/contour sampling → WGS84 raster → `NumericalEarth.regrid_bathymetry` → processed NetCDF. The struct `GeonorgeBathymetry <: AbstractStaticBathymetry` implements the NumericalEarth dataset interface.

2. **Grid** (`src/Grids.jl`) — adds a constructor `ImmersedBoundaryGrid(filepath, arch, halo)` that reads the processed bathymetry NetCDF and returns an `ImmersedBoundaryGrid` wrapping a `LatitudeLongitudeGrid` with `PartialCellBottom`.

3. **Forcing** (`src/Forcing.jl`) — loads river/relaxation forcing from NetCDF via `forcing_from_file`. The `ForcingFromFile` struct carries two `FieldTimeSeries` (values + lambdas) and dispatches to flux, advection, or relaxation terms based on the sign of lambda: `λ > 1` → x-flux, `λ < -1` → y-flux, `|λ| < 1` → relaxation. Uses a custom `NetCDFBackend` keeping 2 time indices in memory.

4. **Atmosphere** (`src/Atmospheres/NORA3.jl`) — `MultiYearNORA3` is a dataset wrapper for NORA3 reanalysis NetCDF. `NORA3PrescribedAtmosphere` and `NORA3PrescribedRadiation` construct NumericalEarth `PrescribedAtmosphere`/`PrescribedRadiation` objects backed by `NORA3FieldTimeSeries`. Default data path: `~/FjordSim_data/NORA3/NORA3.nc`.

5. **Boundary conditions** (`src/BoundaryConditions.jl`) — `top_bottom_boundary_conditions` creates wind/heat/salt flux fields at the top and quadratic bottom drag, returning a named tuple `(u, v, T, S)`.

6. **Top-level** (`src/FjordSim.jl`) — `coupled_hydrostatic_simulation` assembles a `HydrostaticFreeSurfaceModel` inside an `OceanSeaIceModel` (NumericalEarth) and returns a `Simulation`. Also patches `compute_bounding_indices` from NumericalEarth to prevent off-by-one errors with custom longitude/latitude grids.

7. **Dataset adapters** (`src/FDatasets.jl`) — `DSForcing` and `DSResults` are NumericalEarth dataset wrappers for local FjordSim NetCDF files (forcing inputs and simulation outputs), used for initial conditions and restart.

8. **Configs** (`configs/*.toml`) — TOML files describing grid dimensions, bathymetry prep parameters, and forcing download URLs for each fjord setup.

## Key conventions

- Bathymetry convention: `h < 0` = below sea level (bottom height), `h >= 0` = land. `ImmersedBoundaryGrid` loader detects and converts legacy files where `h` was stored as positive depth.
- `recursive_merge` (Utils) deep-merges named tuples; used to combine `top_bottom_boundary_conditions` output with open boundary conditions before `FieldBoundaryConditions`.
- `SetupConfig` module (referenced from scripts but currently deleted per git status) handled TOML config loading; its absence may affect `bathymetry_prepare.jl`.
- Data files default to `~/FjordSim_data/` and results to `~/FjordSim_results/`.
