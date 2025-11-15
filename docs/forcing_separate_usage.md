# Forcing from Separate Files: Usage Guide

This guide explains how to use the new `forcing_separate.jl` module to apply open boundary conditions and river forcing from separate NetCDF files.

## Overview

The module provides three main forcing configurations:

1. **Combined forcing**: Open boundaries + Rivers
2. **Boundary forcing only**: For domains without rivers
3. **River forcing only**: For domains with rivers but using other methods for boundaries

## Key Concepts

### Open Boundary Conditions (OBC)

Open boundary conditions are applied at the domain edges (north, south, east, west) with a buffer zone extending into the interior:

```
Domain boundary → Buffer zone (10 cells) → Interior
    ↓                    ↓                   ↓
 Full forcing       Gradual decay        No forcing
   (100%)           (100% → 0%)            (0%)
```

**How it works:**
- Uses relaxation to nudge model state towards boundary values
- Strength decreases linearly from boundary into interior
- Prevents sharp discontinuities at boundaries
- Default buffer width: 10 grid points (adjustable)

### River Forcing

Rivers are treated as volume flux sources with associated tracer properties:

```
River mouth location → Volume flux (m³/s) + Tracer values
```

**How it works:**
- Applied only at specified river mouth locations
- Introduces water with different temperature, salinity, etc.
- Flux rate determines how quickly river water mixes with ocean

## File Format Requirements

### Boundary Conditions File

**Required variables** (for each tracer, e.g., T, S):

```
dimensions:
    Nx = grid_size_x ;
    Ny = grid_size_y ;
    Nz = grid_size_z ;
    time = UNLIMITED ;

variables:
    // Boundary values
    float T(time, Nz, Ny, Nx) ;
        T:long_name = "Temperature at boundaries" ;
        T:units = "degrees_Celsius" ;
    
    // Relaxation timescale (1/seconds)
    float T_lambda(time, Nz, Ny, Nx) ;
        T_lambda:long_name = "Temperature relaxation timescale" ;
        T_lambda:units = "1/s" ;
    
    float S(time, Nz, Ny, Nx) ;
        S:long_name = "Salinity at boundaries" ;
        S:units = "psu" ;
    
    float S_lambda(time, Nz, Ny, Nx) ;
        S_lambda:long_name = "Salinity relaxation timescale" ;
        S_lambda:units = "1/s" ;
    
    double time(time) ;
        time:long_name = "time" ;
        time:units = "seconds since START_DATE" ;
```

**Data distribution:**
- Non-zero only at/near boundaries
- Interior cells should be zero or fill value (-999)
- Lambda typical range: 1e-5 to 1e-3 (1/s)
  - Smaller = weaker relaxation (longer timescale)
  - Larger = stronger relaxation (shorter timescale)

**Example values:**
```julia
# At boundary (i=1 or i=Nx, j=1 or j=Ny)
T[boundary_cells] = 10.0  # °C
T_lambda[boundary_cells] = 1e-4  # 1/s (relaxation time ~2.8 hours)

# In interior
T[interior_cells] = -999.0  # Fill value (ignored)
T_lambda[interior_cells] = 0.0  # No forcing
```

### River Forcing File

**Required variables** (for each tracer):

```
dimensions:
    Nx = grid_size_x ;
    Ny = grid_size_y ;
    Nz = grid_size_z ;
    time = UNLIMITED ;

variables:
    // Volume flux at river mouths
    float T_flux(time, Nz, Ny, Nx) ;
        T_flux:long_name = "Volume flux for temperature forcing" ;
        T_flux:units = "m3/s" ;
    
    // River water properties
    float T(time, Nz, Ny, Nx) ;
        T:long_name = "Temperature in river water" ;
        T:units = "degrees_Celsius" ;
    
    float S_flux(time, Nz, Ny, Nx) ;
        S_flux:long_name = "Volume flux for salinity forcing" ;
        S_flux:units = "m3/s" ;
    
    float S(time, Nz, Ny, Nx) ;
        S:long_name = "Salinity in river water" ;
        S:units = "psu" ;
    
    double time(time) ;
        time:long_name = "time" ;
        time:units = "seconds since START_DATE" ;
```

**Data distribution:**
- Non-zero only at river mouth locations
- All other cells should be zero or fill value
- Flux is positive for inflow

**Example values:**
```julia
# At river mouth (i=50, j=100, k=Nz surface layer)
T_flux[river_mouth] = 100.0  # m³/s
T[river_mouth] = 5.0  # °C (cold river water)
S_flux[river_mouth] = 100.0  # m³/s (same flux for all tracers)
S[river_mouth] = 0.1  # psu (fresh water)

# Elsewhere
T_flux[other_cells] = 0.0  # No flux
```

## Usage Examples

### Example 1: Combined Forcing (Boundaries + Rivers)

```julia
using FjordsSim

# Setup with both boundary and river forcing
setup = setup_region(
    # ... other parameters ...
    
    forcing_callable = forcing_from_separate_files,
    forcing_args = (
        grid_ref = grid_ref,
        boundary_filepath = joinpath(homedir(), "FjordsSim_data", "myfjord", "boundary_conditions.nc"),
        river_filepath = joinpath(homedir(), "FjordsSim_data", "myfjord", "river_forcing.nc"),
        tracers = tracers,
        buffer_width = 10,  # Optional: default is 10
    ),
    
    # ... other parameters ...
)
```

### Example 2: Boundary Forcing Only

```julia
setup = setup_region(
    # ... other parameters ...
    
    forcing_callable = forcing_from_boundaries_only,
    forcing_args = (
        grid_ref = grid_ref,
        boundary_filepath = joinpath(homedir(), "FjordsSim_data", "myfjord", "boundary_conditions.nc"),
        tracers = tracers,
        buffer_width = 15,  # Optional: wider buffer zone
    ),
    
    # ... other parameters ...
)
```

### Example 3: River Forcing Only

```julia
setup = setup_region(
    # ... other parameters ...
    
    forcing_callable = forcing_from_rivers_only,
    forcing_args = (
        grid_ref = grid_ref,
        river_filepath = joinpath(homedir(), "FjordsSim_data", "myfjord", "river_forcing.nc"),
        tracers = tracers,
    ),
    
    # ... other parameters ...
)
```

### Example 4: Adjusting Buffer Width

The buffer zone width can be adjusted based on your needs:

```julia
# Narrow buffer (5 cells) - sharper transition
buffer_width = 5

# Standard buffer (10 cells) - default
buffer_width = 10

# Wide buffer (20 cells) - smoother transition
buffer_width = 20

forcing_args = (
    grid_ref = grid_ref,
    boundary_filepath = "...",
    river_filepath = "...",
    tracers = tracers,
    buffer_width = buffer_width,
)
```

## Creating Input Files

### Python Example: Creating Boundary Conditions File

```python
import numpy as np
import netCDF4 as nc

# Grid dimensions
Nx, Ny, Nz = 100, 120, 30
Nt = 24  # 24 time steps

# Create file
ds = nc.Dataset('boundary_conditions.nc', 'w', format='NETCDF4')

# Dimensions
ds.createDimension('Nx', Nx)
ds.createDimension('Ny', Ny)
ds.createDimension('Nz', Nz)
ds.createDimension('time', None)  # Unlimited

# Time variable
time_var = ds.createVariable('time', 'f8', ('time',))
time_var.units = 'seconds since 2024-01-01'
time_var[:] = np.arange(0, 24*3600, 3600)  # Hourly

# Temperature at boundaries
T_var = ds.createVariable('T', 'f4', ('time', 'Nz', 'Ny', 'Nx'))
T_var.long_name = 'Temperature at boundaries'
T_var.units = 'degrees_Celsius'

# Temperature relaxation
T_lambda_var = ds.createVariable('T_lambda', 'f4', ('time', 'Nz', 'Ny', 'Nx'))
T_lambda_var.long_name = 'Temperature relaxation timescale'
T_lambda_var.units = '1/s'

# Initialize arrays
T_data = np.full((Nt, Nz, Ny, Nx), -999.0)  # Fill value
T_lambda_data = np.zeros((Nt, Nz, Ny, Nx))

# Set boundary values
buffer_width = 10
lambda_boundary = 1e-4  # 1/s

# West boundary
for i in range(buffer_width):
    weight = 1.0 - i / buffer_width
    T_data[:, :, :, i] = 10.0  # 10°C
    T_lambda_data[:, :, :, i] = lambda_boundary * weight

# East boundary
for i in range(Nx - buffer_width, Nx):
    weight = 1.0 - (Nx - 1 - i) / buffer_width
    T_data[:, :, :, i] = 12.0  # 12°C
    T_lambda_data[:, :, :, i] = lambda_boundary * weight

# Similar for north and south...

# Write data
T_var[:] = T_data
T_lambda_var[:] = T_lambda_data

# Repeat for salinity and other tracers...

ds.close()
```

### Python Example: Creating River Forcing File

```python
import numpy as np
import netCDF4 as nc

# Grid dimensions
Nx, Ny, Nz = 100, 120, 30
Nt = 24

# River location
river_i, river_j = 50, 100  # Grid indices
river_discharge = 150.0  # m³/s

# Create file
ds = nc.Dataset('river_forcing.nc', 'w', format='NETCDF4')

# Dimensions and time (same as above)
ds.createDimension('Nx', Nx)
ds.createDimension('Ny', Ny)
ds.createDimension('Nz', Nz)
ds.createDimension('time', None)

time_var = ds.createVariable('time', 'f8', ('time',))
time_var.units = 'seconds since 2024-01-01'
time_var[:] = np.arange(0, 24*3600, 3600)

# Temperature flux and value
T_flux_var = ds.createVariable('T_flux', 'f4', ('time', 'Nz', 'Ny', 'Nx'))
T_flux_var.units = 'm3/s'
T_var = ds.createVariable('T', 'f4', ('time', 'Nz', 'Ny', 'Nx'))
T_var.units = 'degrees_Celsius'

# Initialize
T_flux_data = np.zeros((Nt, Nz, Ny, Nx))
T_data = np.full((Nt, Nz, Ny, Nx), -999.0)

# Set river discharge at surface layer (k = Nz-1 for surface)
k_surface = Nz - 1
T_flux_data[:, k_surface, river_j, river_i] = river_discharge
T_data[:, k_surface, river_j, river_i] = 8.0  # 8°C river water

# Add time variation if needed
for t in range(Nt):
    # Seasonal variation
    seasonal_factor = 1.0 + 0.3 * np.sin(2 * np.pi * t / Nt)
    T_flux_data[t, k_surface, river_j, river_i] *= seasonal_factor

T_flux_var[:] = T_flux_data
T_var[:] = T_data

# Repeat for salinity and other tracers...

ds.close()
```

## Tips and Best Practices

### 1. Choosing Relaxation Timescales

The relaxation timescale λ controls how strongly the model is nudged toward boundary values:

```julia
# Weak nudging (long timescale ≈ 1 day)
λ = 1e-5  # 1/s

# Moderate nudging (timescale ≈ 3 hours)
λ = 1e-4  # 1/s

# Strong nudging (timescale ≈ 10 minutes)
λ = 1e-3  # 1/s
```

**Guidelines:**
- Use stronger nudging for well-constrained boundaries (e.g., connected to larger ocean)
- Use weaker nudging for uncertain boundary data
- Match nudging strength to expected boundary layer thickness

### 2. Buffer Zone Width

Choose buffer width based on:
- Grid resolution
- Physical scales (internal Rossby radius, boundary layer thickness)
- Numerical stability

```julia
# Fine grid (Δx ≈ 100m) → larger buffer
buffer_width = 20  # ~2 km buffer

# Coarse grid (Δx ≈ 1km) → smaller buffer
buffer_width = 5  # ~5 km buffer
```

### 3. River Flux Distribution

For wide rivers, distribute flux over multiple cells:

```python
# Single point source
T_flux_data[:, k, j, i] = total_discharge

# Distributed over 3x3 cells
discharge_per_cell = total_discharge / 9
for di in range(-1, 2):
    for dj in range(-1, 2):
        T_flux_data[:, k, j+dj, i+di] = discharge_per_cell
```

### 4. Handling Time Variation

Both boundary and river data can be time-varying:

```python
# Tidal variation at boundaries
for t in range(Nt):
    tidal_phase = 2 * np.pi * t / 12  # 12-hour period
    sea_level = 0.5 * np.sin(tidal_phase)  # ±0.5m
    T_data[t, :, :, boundary_indices] += sea_level * 0.1  # Slight T change

# Seasonal river flow
for t in range(Nt):
    day_of_year = t * 3600 / 86400  # Convert to days
    seasonal = 1.0 + 0.5 * np.sin(2 * np.pi * day_of_year / 365)
    T_flux_data[t, ...] *= seasonal
```

### 5. Debugging

Check your forcing files before running:

```julia
using NCDatasets

# Check boundary file
ds = Dataset("boundary_conditions.nc")
println("Variables: ", keys(ds))
println("T shape: ", size(ds["T"]))
println("T range: ", extrema(ds["T"][:]))
println("T_lambda range: ", extrema(ds["T_lambda"][:]))

# Find where boundary forcing is applied
T_lambda = ds["T_lambda"][:, :, :, 1]  # First time step
boundary_cells = findall(T_lambda .> 0)
println("Boundary cells: ", length(boundary_cells))

close(ds)
```

## Troubleshooting

### Problem: Model becomes unstable near boundaries

**Solutions:**
1. Reduce relaxation timescale (smaller λ)
2. Increase buffer zone width
3. Check boundary data for unrealistic values
4. Ensure smooth time variation in boundary data

### Problem: River water not mixing properly

**Solutions:**
1. Check river flux magnitude (not too large)
2. Verify river location is in water (not land)
3. Ensure tracer values are reasonable
4. Consider distributing flux over multiple cells

### Problem: "DimensionMismatch" error

**Solutions:**
1. Verify NetCDF dimensions match grid size
2. Check dimension order: (time, Nz, Ny, Nx)
3. Ensure all variables have correct dimensions

### Problem: No forcing effect observed

**Solutions:**
1. Check that λ or flux values are non-zero where expected
2. Verify tracer names match between file and model
3. Confirm file paths are correct
4. Check for fill values (-999) in forcing regions

## Performance Considerations

- **Memory**: Files use in-memory backend for 2 time steps (configurable)
- **I/O**: Data loaded on-demand as simulation progresses
- **GPU**: All forcing structures are GPU-compatible via Adapt.jl

## References

- Oceananigans.jl documentation: https://clima.github.io/OceananigansDocumentation/stable/
- NetCDF Climate and Forecast (CF) Conventions
- Ocean modeling best practices for boundary conditions
