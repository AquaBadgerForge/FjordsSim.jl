# New Forcing Module Summary

## Overview

A new forcing module has been created for FjordsSim that handles **open boundary conditions** and **river forcing** from **separate NetCDF files**. This provides more flexibility and clarity compared to combining both in a single file.

## Key Features

### 1. Open Boundary Conditions (OBC)
- Applied at domain boundaries (north, south, east, west)
- **Buffer zone**: Extends 10 grid points (configurable) into the domain
- **Smooth transition**: Linear decay from full forcing at boundary to zero in interior
- **Relaxation-based**: Nudges model state toward boundary values with adjustable timescale
- **Prevents discontinuities**: No sharp jumps at boundaries

### 2. River Forcing
- **Point sources**: Applied at specific river mouth locations
- **Volume flux**: Specifies discharge rate (m³/s) 
- **Tracer properties**: Temperature, salinity, etc. of river water
- **Time-varying**: Both discharge and properties can vary with time
- **Multiple rivers**: Support for any number of river sources

### 3. Flexible Configuration
Three usage modes:
1. **Combined**: Both boundaries and rivers (most common)
2. **Boundaries only**: For domains without rivers
3. **Rivers only**: When using alternative boundary methods

## Files Created

### Core Module
- **`src/forcing_separate.jl`** (708 lines)
  - `OpenBoundaryForcing`: Structure for boundary conditions
  - `RiverForcing`: Structure for river inputs
  - `BoundaryZoneMask`: Defines buffer zones
  - `create_boundary_mask()`: Generate buffer zone masks
  - `forcing_from_separate_files()`: Load both types
  - `forcing_from_boundaries_only()`: Load only boundaries
  - `forcing_from_rivers_only()`: Load only rivers
  - Helper functions for loading NetCDF data

### Documentation
- **`docs/forcing_separate_usage.md`** (664 lines)
  - Complete user guide
  - File format specifications
  - Detailed examples
  - Best practices
  - Troubleshooting guide
  - Python code examples for file creation

### Examples
- **`examples/forcing_separate_example.jl`** (189 lines)
  - Working Julia examples
  - Three configuration modes
  - Full test simulation
  - Ready-to-use setup functions

- **`examples/create_forcing_files.py`** (437 lines)
  - Python script to generate forcing files
  - Helper class `ForcingFileGenerator`
  - Example for Isafjardardjup
  - Template for custom domains
  - Handles boundary buffer zones automatically

- **`examples/README_forcing.md`** (97 lines)
  - Quick start guide
  - Links to detailed documentation
  - Common troubleshooting tips

### Tests
- **`test/test_forcing_separate.jl`** (149 lines)
  - Unit tests for boundary mask creation
  - Tests for different buffer widths
  - Edge case handling
  - Grid compatibility checks

### Integration
- **`src/FjordsSim.jl`** (modified)
  - Added `include("forcing_separate.jl")`
  - Module now available throughout FjordsSim

## File Format Specifications

### Boundary Conditions File

**Dimensions:**
```
Nx, Ny, Nz = grid dimensions
time = unlimited
```

**Variables (for each tracer, e.g., T, S):**
```
T(time, Nz, Ny, Nx)           - boundary values
T_lambda(time, Nz, Ny, Nx)    - relaxation timescale (1/s)
```

**Data structure:**
- Non-zero only at/near boundaries
- Interior: fill value (-999) or zero
- Lambda controls nudging strength

### River Forcing File

**Dimensions:**
```
Nx, Ny, Nz = grid dimensions
time = unlimited
```

**Variables (for each tracer):**
```
T_flux(time, Nz, Ny, Nx)      - volume flux (m³/s)
T(time, Nz, Ny, Nx)           - river water properties
```

**Data structure:**
- Non-zero only at river mouth locations
- Elsewhere: zero or fill value
- Positive flux for inflow

## Usage Example

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

## Key Concepts

### Buffer Zone

The buffer zone creates a smooth transition:

```
Domain edge → Buffer zone (10 cells) → Interior
     ↓             ↓                       ↓
  100% forcing  100% → 0% forcing      0% forcing
```

This prevents:
- Sharp discontinuities
- Numerical instabilities
- Unrealistic boundary effects

### Relaxation Timescale (λ)

Controls forcing strength:

```julia
λ = 1e-5  # Weak nudging (~1 day timescale)
λ = 1e-4  # Moderate (~3 hours) - typical choice
λ = 1e-3  # Strong (~10 minutes)
```

Formula: `tendency = -λ * (model_value - boundary_value)`

### River Volume Flux

Introduces water with different properties:

```julia
tendency = (flux / volume) * (river_value - model_value)
```

Accounts for:
- Discharge rate
- Grid cell geometry
- Property differences between river and ocean

## Implementation Details

### Kernel Functions

Both forcing types are implemented as callable structures that work in GPU/CPU kernels:

```julia
@inline function (p::OpenBoundaryForcing)(i, j, k, grid, clock, fields)
    # Get boundary value and mask weight
    # Apply masked relaxation
    # Return tendency
end

@inline function (p::RiverForcing)(i, j, k, grid, clock, fields)
    # Get flux and river properties
    # Apply volume flux forcing
    # Return tendency
end
```

### GPU Compatibility

All structures use `Adapt.jl` for seamless GPU transfer:
```julia
Adapt.adapt_structure(to, forcing) = ...
```

### Memory Efficiency

Uses `NetCDFBackend` to keep only 2 time steps in memory at once, loading data as needed during simulation.

## Python Helper Tools

The `create_forcing_files.py` script provides:

1. **`ForcingFileGenerator` class**:
   - Handles grid dimensions
   - Manages time coordinates
   - Creates NetCDF files with proper structure

2. **Automatic buffer zone creation**:
   - Calculates linear decay weights
   - Handles corners correctly
   - Works for any buffer width

3. **River point source setup**:
   - Place rivers at specific grid cells
   - Set discharge and properties
   - Support for multiple rivers

## Testing

Run tests with:

```julia
include("test/test_forcing_separate.jl")
```

Tests verify:
- Boundary mask creation
- Different buffer widths
- Grid compatibility
- Edge cases
- Structure definitions
- Function availability

## Integration with Existing Code

The new module:
- ✓ Works alongside existing `forcing.jl`
- ✓ Uses same `NetCDFBackend` and data loading infrastructure
- ✓ Compatible with all tracers (T, S, biogeochemical, etc.)
- ✓ Follows Oceananigans forcing conventions
- ✓ GPU/CPU compatible via Adapt.jl
- ✓ No changes needed to existing simulations

## Advantages Over Single-File Approach

1. **Clarity**: Separate concerns (boundaries vs. rivers)
2. **Flexibility**: Can use either or both independently
3. **Easier file creation**: Simpler data structure in each file
4. **Better organization**: Clear distinction between forcing types
5. **Debugging**: Easier to isolate issues
6. **Reusability**: Can reuse boundary file with different river scenarios

## Next Steps for Users

1. **Create forcing files**:
   ```bash
   python examples/create_forcing_files.py
   ```

2. **Adjust for your domain**:
   - Set correct grid dimensions
   - Define boundary values from observations/models
   - Identify river mouth locations
   - Estimate discharge rates

3. **Test with simulation**:
   ```julia
   include("examples/forcing_separate_example.jl")
   simulation = run_example()
   ```

4. **Validate results**:
   - Check forcing is applied correctly
   - Verify boundary transitions are smooth
   - Ensure rivers affect tracers as expected
   - Monitor for numerical stability

5. **Tune parameters**:
   - Adjust λ for boundary nudging strength
   - Modify buffer_width for smoother/sharper transitions
   - Calibrate river discharge and properties

## Documentation Structure

```
FjordsSim.jl/
├── src/
│   ├── forcing.jl                    # Original forcing module
│   ├── forcing_separate.jl           # NEW: Separate forcing module
│   └── FjordsSim.jl                  # Modified: includes new module
├── docs/
│   └── forcing_separate_usage.md     # NEW: Comprehensive user guide
├── examples/
│   ├── forcing_separate_example.jl   # NEW: Julia examples
│   ├── create_forcing_files.py       # NEW: Python file creator
│   └── README_forcing.md             # NEW: Quick start guide
└── test/
    └── test_forcing_separate.jl      # NEW: Unit tests
```

## Code Statistics

- **Total lines of code**: ~2,200 lines
- **Core module**: 708 lines
- **Documentation**: 664 lines
- **Examples**: 626 lines (Julia + Python)
- **Tests**: 149 lines

All code includes:
- Extensive comments
- Docstrings for all functions
- Type annotations
- Error handling
- Usage examples

## Comments and Explanations

Every function includes:
- **Purpose**: What it does
- **Parameters**: Detailed descriptions
- **Returns**: What it returns
- **Expected file format**: For I/O functions
- **Examples**: Usage patterns
- **Notes**: Important caveats or tips

## Support

For questions or issues:
1. Read `docs/forcing_separate_usage.md`
2. Check examples in `examples/`
3. Run tests to verify installation
4. Review troubleshooting section in docs

## Future Enhancements (Potential)

Ideas for extension:
- Tidal boundary conditions with harmonic constituents
- Seasonal river discharge patterns
- Data assimilation at boundaries
- Adaptive buffer zones based on flow conditions
- Boundary condition validation tools
- Visualization of forcing fields

## Conclusion

The new forcing module provides a robust, well-documented, and flexible system for applying open boundary conditions and river forcing in FjordsSim. It maintains compatibility with existing code while offering improved clarity and usability for complex forcing scenarios.
