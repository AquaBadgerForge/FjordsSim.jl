# FjordsSim Forcing Examples

This directory contains examples for using the new separate forcing module in FjordsSim.

## Files

### Julia Examples

- **`forcing_separate_example.jl`**: Complete Julia examples showing how to use the new forcing module
  - Combined boundary + river forcing
  - Boundary-only forcing
  - River-only forcing
  - Full working simulation example

### Python Examples

- **`create_forcing_files.py`**: Python script to generate NetCDF forcing files
  - Functions to create boundary condition files
  - Functions to create river forcing files
  - Example configurations for Isafjardardjup and custom domains

## Quick Start

### 1. Create Forcing Files (Python)

```bash
cd examples
python create_forcing_files.py
```

Edit the script to:
- Set your grid dimensions (Nx, Ny, Nz)
- Define boundary values for each tracer
- Specify river locations and properties
- Adjust relaxation timescales

### 2. Use in FjordsSim (Julia)

```julia
include("examples/forcing_separate_example.jl")

# Run test simulation
simulation = run_example()
```

Or integrate into your own setup:

```julia
using FjordsSim

setup = setup_region(
    # ... other parameters ...
    
    forcing_callable = forcing_from_separate_files,
    forcing_args = (
        grid_ref = grid_ref,
        boundary_filepath = "path/to/boundary_conditions.nc",
        river_filepath = "path/to/river_forcing.nc",
        tracers = (:T, :S),
        buffer_width = 10,
    ),
    
    # ... other parameters ...
)

simulation = coupled_hydrostatic_simulation(setup)
```

## Key Features

### Open Boundary Conditions
- Applied at domain edges (north, south, east, west)
- Includes buffer zone extending into domain (default: 10 cells)
- Uses relaxation to nudge model toward boundary values
- Prevents sharp discontinuities

### River Forcing
- Point sources at river mouth locations
- Volume flux with tracer properties
- Time-varying discharge and properties
- Can have multiple rivers

## File Formats

### Boundary Conditions File
```
Variables required for each tracer (e.g., T, S):
- T(time, Nz, Ny, Nx): boundary values
- T_lambda(time, Nz, Ny, Nx): relaxation timescales (1/s)
```

### River Forcing File
```
Variables required for each tracer:
- T_flux(time, Nz, Ny, Nx): volume flux (m³/s)
- T(time, Nz, Ny, Nx): tracer value in river water
```

## Documentation

For detailed documentation, see:
- **`../docs/forcing_separate_usage.md`**: Complete user guide with examples
  - Detailed explanations of concepts
  - File format specifications
  - Best practices and troubleshooting
  - Performance considerations

## Tips

1. **Grid dimensions**: Ensure NetCDF dimensions exactly match your FjordsSim grid
2. **Buffer width**: Choose based on grid resolution and physical scales (typically 5-20 cells)
3. **Relaxation timescales**: 
   - λ = 1e-5 → weak nudging (~1 day)
   - λ = 1e-4 → moderate nudging (~3 hours)
   - λ = 1e-3 → strong nudging (~10 minutes)
4. **River locations**: Verify river mouths are in water cells, not land
5. **Time variation**: Add seasonal or tidal variations for realism

## Troubleshooting

**Problem**: DimensionMismatch error
- Check that Nx, Ny, Nz in NetCDF files match grid exactly
- Verify dimension order: (time, Nz, Ny, Nx)

**Problem**: No forcing effect
- Check that λ or flux values are non-zero
- Verify data is not all fill values (-999)
- Confirm tracer names match between file and model

**Problem**: Instability near boundaries
- Reduce relaxation timescale (smaller λ)
- Increase buffer width
- Check for unrealistic boundary values

For more help, see the detailed documentation.
