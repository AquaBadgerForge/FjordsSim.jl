"""
Module for handling forcing from separate files for open boundary conditions and rivers.

This module provides functionality to:
1. Read open boundary conditions (OBC) from one NetCDF file
2. Read river forcing from a separate NetCDF file
3. Apply both forcings with proper spatial distribution

Open Boundary Conditions:
- Applied at domain boundaries (north, south, east, west)
- Extended into a buffer zone (default 10 grid points) inside the domain
- Use relaxation to nudge model state towards boundary values

River Forcing:
- Point sources or distributed sources (e.g., river mouths)
- Applied as volume flux with tracer properties
- Can be time-varying
"""

using Oceananigans.OutputReaders: FieldTimeSeries, Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.BoundaryConditions: fill_halo_regions!, FieldBoundaryConditions
using Oceananigans.Fields: interior
using Oceananigans.Forcings: Forcing
using Oceananigans.Grids: Center, Face, nodes
using Oceananigans.Units: hours
using Oceananigans.Operators: Ax, Ay, Az, volume
using Dates: DateTime, Year, Second
using Statistics: mean
using NCDatasets: Dataset
using Adapt

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend

# Note: This module uses NetCDFBackend, oceananigans_fieldname, DATA_LOCATION,
# and load_from_netcdf from forcing.jl which is included earlier in the module chain

#=============================================================================
# Open Boundary Condition Forcing
=============================================================================#

"""
    BoundaryZoneMask

Structure to define which cells are in the boundary relaxation zone.

# Fields
- `mask::Array{Float64, 3}`: 3D array where:
  - 0.0 = interior (no boundary forcing)
  - 1.0 = at boundary (strongest forcing)
  - 0.0-1.0 = buffer zone (linearly decreasing forcing strength)
"""
struct BoundaryZoneMask{A <: AbstractArray}
    mask::A
end

"""
    create_boundary_mask(grid, buffer_width::Int=10)

Create a mask for boundary relaxation zones.

The mask defines where open boundary conditions are applied:
- Boundaries: Full strength (mask = 1.0)
- Buffer zone: Linear decay from 1.0 to 0.0 over `buffer_width` points
- Interior: No forcing (mask = 0.0)

# Arguments
- `grid`: Oceananigans grid
- `buffer_width::Int`: Number of grid points for the buffer zone (default: 10)

# Returns
- `BoundaryZoneMask`: Mask structure containing 3D mask array

# Example
The mask transitions as:
boundary -> buffer zone -> interior
  1.0    ->  0.9, 0.8, ... 0.1 -> 0.0
"""
function create_boundary_mask(grid, buffer_width::Int=10)
    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz
    
    # Initialize mask to zero (interior cells)
    mask = zeros(Float64, Nx, Ny, Nz)
    
    # West boundary (i = 1) and buffer zone
    for i in 1:min(buffer_width, Nx)
        # Linear decay: 1.0 at boundary, decreasing inward
        weight = 1.0 - (i - 1) / buffer_width
        mask[i, :, :] .= max.(mask[i, :, :], weight)
    end
    
    # East boundary (i = Nx) and buffer zone
    for i in max(1, Nx - buffer_width + 1):Nx
        weight = 1.0 - (Nx - i) / buffer_width
        mask[i, :, :] .= max.(mask[i, :, :], weight)
    end
    
    # South boundary (j = 1) and buffer zone
    for j in 1:min(buffer_width, Ny)
        weight = 1.0 - (j - 1) / buffer_width
        mask[:, j, :] .= max.(mask[:, j, :], weight)
    end
    
    # North boundary (j = Ny) and buffer zone
    for j in max(1, Ny - buffer_width + 1):Ny
        weight = 1.0 - (Ny - j) / buffer_width
        mask[:, j, :] .= max.(mask[:, j, :], weight)
    end
    
    return BoundaryZoneMask(mask)
end

"""
    OpenBoundaryForcing{FTS, V, M}

Forcing structure for open boundary conditions with buffer zone.

# Fields
- `fts_value::FTS`: FieldTimeSeries containing boundary values
- `fts_λ::FTS`: FieldTimeSeries containing relaxation timescales
- `fieldname::V`: Type identifier for the field (Temperature, Salinity, etc.)
- `boundary_mask::M`: BoundaryZoneMask defining the spatial extent of forcing
"""
struct OpenBoundaryForcing{FTS, V, M}
    fts_value::FTS
    fts_λ::FTS
    fieldname::V
    boundary_mask::M
end

Adapt.adapt_structure(to, p::OpenBoundaryForcing) =
    OpenBoundaryForcing(
        Adapt.adapt(to, p.fts_value), 
        Adapt.adapt(to, p.fts_λ), 
        Adapt.adapt(to, p.fieldname),
        Adapt.adapt(to, p.boundary_mask)
    )

"""
    (p::OpenBoundaryForcing)(i, j, k, grid, clock, fields)

Kernel function for open boundary forcing.

Applies relaxation forcing at boundaries and buffer zones:
    tendency = -λ_effective * (model_value - boundary_value)

where λ_effective = λ * mask_weight

# Arguments
- `i, j, k`: Grid indices
- `grid`: Grid structure
- `clock`: Simulation clock
- `fields`: Named tuple of model fields

# Returns
- Forcing tendency (rate of change) for the field
"""
@inline function (p::OpenBoundaryForcing{FTS,V,M})(i, j, k, grid, clock, fields) where {FTS,V,M}
    # Get boundary value and base relaxation timescale from file
    value = @inbounds p.fts_value[i, j, k, Time(clock.time)]
    λ = @inbounds p.fts_λ[i, j, k, Time(clock.time)]
    
    # Get mask weight for this location (0.0 to 1.0)
    mask_weight = @inbounds p.boundary_mask.mask[i, j, k]
    
    FT = eltype(grid)
    
    # Only apply forcing if:
    # 1. We're in the boundary zone (mask_weight > 0)
    # 2. Valid data (value > -990 is a common fill value check)
    if mask_weight > 0.0 && value > -990
        # Get current model field value
        model_value = @inbounds fields[i, j, k, p.fieldname]
        
        # Apply masked relaxation: stronger near boundary, weaker in buffer
        λ_effective = convert(FT, λ * mask_weight)
        tendency = -λ_effective * (model_value - convert(FT, value))
        
        return tendency
    else
        return zero(FT)
    end
end

regularize_forcing(forcing::OpenBoundaryForcing, field, field_name, model_field_names) = forcing

#=============================================================================
# River Forcing
=============================================================================#

"""
    RiverForcing{FTS, V}

Forcing structure for river inputs (point or distributed sources).

Rivers are treated as volume flux sources with associated tracer properties.
The forcing is applied only at specific locations defined in the river file.

# Fields
- `fts_flux::FTS`: FieldTimeSeries containing volume flux (m³/s)
- `fts_value::FTS`: FieldTimeSeries containing tracer value (e.g., temperature, salinity)
- `fieldname::V`: Type identifier for the field
"""
struct RiverForcing{FTS, V}
    fts_flux::FTS
    fts_value::FTS
    fieldname::V
end

Adapt.adapt_structure(to, p::RiverForcing) =
    RiverForcing(
        Adapt.adapt(to, p.fts_flux), 
        Adapt.adapt(to, p.fts_value), 
        Adapt.adapt(to, p.fieldname)
    )

"""
    (p::RiverForcing)(i, j, k, grid, clock, fields)

Kernel function for river forcing.

Applies volume flux with associated tracer value:
    tendency = flux * (river_value - model_value) * area / volume

This formulation:
1. Introduces river water with its properties (river_value)
2. Scales by the flux rate
3. Accounts for grid cell geometry

# Arguments
- `i, j, k`: Grid indices
- `grid`: Grid structure
- `clock`: Simulation clock
- `fields`: Named tuple of model fields

# Returns
- Forcing tendency (rate of change) for the tracer
"""
@inline function (p::RiverForcing{FTS,V})(i, j, k, grid, clock, fields) where {FTS,V}
    # Get river flux (m³/s) and tracer value
    flux = @inbounds p.fts_flux[i, j, k, Time(clock.time)]
    value = @inbounds p.fts_value[i, j, k, Time(clock.time)]
    
    FT = eltype(grid)
    
    # Only apply forcing where flux > 0 and valid data
    if flux > 0.0 && value > -990
        # Get cell volume and characteristic area
        vol = volume(i, j, k, grid, Center(), Center(), Center())
        
        # Avoid division by zero in masked/land cells
        if vol == 0
            return zero(FT)
        end
        
        # Get current model field value
        model_value = @inbounds fields[i, j, k, p.fieldname]
        
        # River forcing: introduces water with different properties
        # Positive flux brings in river water, affecting the tracer
        area = Az(i, j, k, grid, Center(), Center(), Face())  # Horizontal area
        flux_per_volume = convert(FT, flux) * convert(FT, area) / convert(FT, vol)
        
        # Tendency: rate at which river properties replace model properties
        tendency = flux_per_volume * (convert(FT, value) - model_value)
        
        return tendency
    else
        return zero(FT)
    end
end

regularize_forcing(forcing::RiverForcing, field, field_name, model_field_names) = forcing

#=============================================================================
# File Loading Functions
=============================================================================#

"""
    load_open_boundary_forcing(
        filepath::String, 
        var_name::String, 
        grid, 
        time_indices_in_memory, 
        backend,
        buffer_width::Int=10
    )

Load open boundary condition forcing for a single variable.

# Arguments
- `filepath`: Path to NetCDF file containing boundary data
- `var_name`: Variable name (e.g., "T", "S", "u", "v")
- `grid`: Oceananigans grid
- `time_indices_in_memory`: Tuple of time indices to keep in memory
- `backend`: NetCDFBackend for time series management
- `buffer_width`: Width of buffer zone in grid points (default: 10)

# Returns
- Named tuple with OpenBoundaryForcing for the variable

# Expected NetCDF structure:
The file should contain:
- `var_name`: Variable data (Nx, Ny, Nz, time)
- `var_name * "_lambda"`: Relaxation timescale (Nx, Ny, Nz, time)
- `time`: Time coordinate

Note: Boundary data should be non-zero only at/near boundaries.
"""
function load_open_boundary_forcing(
    filepath::String, 
    var_name::String, 
    grid, 
    time_indices_in_memory, 
    backend,
    buffer_width::Int=10
)
    # Get field location and type
    field_name = oceananigans_fieldname[Symbol(var_name)]
    LX, LY, LZ = DATA_LOCATION[field_name]
    
    # Get grid dimensions
    grid_size_tupled = size.(nodes(grid, (LX(), LY(), LZ())))
    grid_size = Tuple(x[1] for x in grid_size_tupled)
    
    # Load boundary values and relaxation timescales from file
    data, times = load_from_netcdf(; 
        path = filepath, 
        var_name, 
        grid_size, 
        time_indices_in_memory
    )
    dataλ, timesλ = load_from_netcdf(; 
        path = filepath, 
        var_name = var_name * "_lambda", 
        grid_size, 
        time_indices_in_memory
    )
    
    # Create boundary conditions
    boundary_conditions = FieldBoundaryConditions(grid, (LX(), LY(), LZ()))
    
    # Create FieldTimeSeries for boundary values
    fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
        path = filepath,
        name = var_name,
    )
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)
    
    # Create FieldTimeSeries for relaxation timescales
    ftsλ = FieldTimeSeries{LX,LY,LZ}(
        grid,
        timesλ;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
        path = filepath,
        name = var_name * "_lambda",
    )
    copyto!(interior(ftsλ, :, :, :, :), dataλ)
    fill_halo_regions!(ftsλ)
    
    # Create boundary mask
    boundary_mask = create_boundary_mask(grid, buffer_width)
    
    # Create OpenBoundaryForcing structure
    _forcing = OpenBoundaryForcing(fts, ftsλ, field_name, boundary_mask)
    result = NamedTuple{(Symbol(var_name),)}((_forcing,))
    
    return result
end

"""
    load_river_forcing(
        filepath::String, 
        var_name::String, 
        grid, 
        time_indices_in_memory, 
        backend
    )

Load river forcing for a single variable.

# Arguments
- `filepath`: Path to NetCDF file containing river data
- `var_name`: Variable name (e.g., "T", "S")
- `grid`: Oceananigans grid
- `time_indices_in_memory`: Tuple of time indices to keep in memory
- `backend`: NetCDFBackend for time series management

# Returns
- Named tuple with RiverForcing for the variable

# Expected NetCDF structure:
The file should contain:
- `var_name * "_flux"`: Volume flux (m³/s) at river locations (Nx, Ny, Nz, time)
- `var_name`: Tracer value in river water (Nx, Ny, Nz, time)
- `time`: Time coordinate

Note: Flux should be non-zero only at river mouth locations.
"""
function load_river_forcing(
    filepath::String, 
    var_name::String, 
    grid, 
    time_indices_in_memory, 
    backend
)
    # Get field location and type
    field_name = oceananigans_fieldname[Symbol(var_name)]
    LX, LY, LZ = DATA_LOCATION[field_name]
    
    # Get grid dimensions
    grid_size_tupled = size.(nodes(grid, (LX(), LY(), LZ())))
    grid_size = Tuple(x[1] for x in grid_size_tupled)
    
    # Load river flux and tracer values from file
    data_flux, times_flux = load_from_netcdf(; 
        path = filepath, 
        var_name = var_name * "_flux", 
        grid_size, 
        time_indices_in_memory
    )
    data_value, times_value = load_from_netcdf(; 
        path = filepath, 
        var_name, 
        grid_size, 
        time_indices_in_memory
    )
    
    # Create boundary conditions
    boundary_conditions = FieldBoundaryConditions(grid, (LX(), LY(), LZ()))
    
    # Create FieldTimeSeries for river flux
    fts_flux = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times_flux;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
        path = filepath,
        name = var_name * "_flux",
    )
    copyto!(interior(fts_flux, :, :, :, :), data_flux)
    fill_halo_regions!(fts_flux)
    
    # Create FieldTimeSeries for river tracer values
    fts_value = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times_value;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
        path = filepath,
        name = var_name,
    )
    copyto!(interior(fts_value, :, :, :, :), data_value)
    fill_halo_regions!(fts_value)
    
    # Create RiverForcing structure
    _forcing = RiverForcing(fts_flux, fts_value, field_name)
    result = NamedTuple{(Symbol(var_name),)}((_forcing,))
    
    return result
end

#=============================================================================
# Combined Forcing Functions
=============================================================================#

"""
    forcing_from_separate_files(
        grid_ref, 
        boundary_filepath::String, 
        river_filepath::String, 
        tracers;
        buffer_width::Int=10
    )

Create forcing from separate open boundary and river files.

This function loads and combines:
1. Open boundary conditions with buffer zones
2. River forcing at specified locations

Both forcings are applied additively to each tracer.

# Arguments
- `grid_ref`: Reference to grid (typically a Ref{Grid})
- `boundary_filepath`: Path to NetCDF file with boundary conditions
- `river_filepath`: Path to NetCDF file with river forcing
- `tracers`: Tuple of tracer symbols (e.g., (:T, :S))
- `buffer_width`: Width of boundary buffer zone in grid points (default: 10)

# Returns
- Named tuple of combined forcings for each variable

# File Requirements

## Boundary File (boundary_filepath):
Should contain for each tracer (e.g., T, S):
- `T`: Temperature at boundaries (Nx, Ny, Nz, time)
- `T_lambda`: Relaxation timescale 1/s (Nx, Ny, Nz, time)
- Similar for S, and optionally u, v

Data should be non-zero only at/near boundaries.

## River File (river_filepath):
Should contain for each tracer:
- `T_flux`: Volume flux m³/s (Nx, Ny, Nz, time)
- `T`: River water temperature (Nx, Ny, Nz, time)
- Similar for S and other tracers

Flux should be non-zero only at river mouth locations.

# Example
```julia
tracers = (:T, :S)
forcing = forcing_from_separate_files(
    grid_ref,
    "boundary_conditions.nc",
    "river_forcing.nc",
    tracers;
    buffer_width = 10
)
```
"""
function forcing_from_separate_files(
    grid_ref, 
    boundary_filepath::String, 
    river_filepath::String, 
    tracers;
    buffer_width::Int=10
)
    grid = grid_ref[]
    
    # Check dimensions match for boundary file
    ds_boundary = Dataset(boundary_filepath)
    grid.underlying_grid.Nx == ds_boundary.dim["Nx"] &&
        grid.underlying_grid.Ny == ds_boundary.dim["Ny"] &&
        grid.underlying_grid.Nz == ds_boundary.dim["Nz"] ||
        throw(DimensionMismatch("boundary file dimensions not equal to grid dimensions"))
    
    # Get available variables in boundary file
    boundary_variables = map(String, tracers) ∩ keys(ds_boundary)
    close(ds_boundary)
    
    # Check dimensions match for river file
    ds_river = Dataset(river_filepath)
    grid.underlying_grid.Nx == ds_river.dim["Nx"] &&
        grid.underlying_grid.Ny == ds_river.dim["Ny"] &&
        grid.underlying_grid.Nz == ds_river.dim["Nz"] ||
        throw(DimensionMismatch("river file dimensions not equal to grid dimensions"))
    
    # Get available variables in river file (check for _flux suffix)
    river_flux_vars = filter(name -> endswith(name, "_flux"), keys(ds_river))
    river_variables = [replace(name, "_flux" => "") for name in river_flux_vars]
    river_variables = river_variables ∩ map(String, tracers)
    close(ds_river)
    
    # Setup backend for time series
    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    
    # Load open boundary forcings
    println("Loading open boundary conditions for: ", boundary_variables)
    boundary_forcings = if !isempty(boundary_variables)
        mapreduce(
            var_name -> load_open_boundary_forcing(
                boundary_filepath, 
                var_name, 
                grid, 
                time_indices_in_memory, 
                backend,
                buffer_width
            ),
            merge,
            boundary_variables,
        )
    else
        NamedTuple()
    end
    
    # Load river forcings
    println("Loading river forcing for: ", river_variables)
    river_forcings = if !isempty(river_variables)
        mapreduce(
            var_name -> load_river_forcing(
                river_filepath, 
                var_name, 
                grid, 
                time_indices_in_memory, 
                backend
            ),
            merge,
            river_variables,
        )
    else
        NamedTuple()
    end
    
    # Combine both forcings - they will be added together by Oceananigans
    # If a variable has both boundary and river forcing, both will be applied
    result = merge(boundary_forcings, river_forcings)
    
    return result
end

"""
    forcing_from_boundaries_only(
        grid_ref, 
        boundary_filepath::String, 
        tracers;
        buffer_width::Int=10
    )

Create forcing from open boundary conditions only (no rivers).

Convenience function for cases where only boundary forcing is needed.

# Arguments
- `grid_ref`: Reference to grid
- `boundary_filepath`: Path to NetCDF file with boundary conditions
- `tracers`: Tuple of tracer symbols
- `buffer_width`: Width of boundary buffer zone in grid points (default: 10)

# Returns
- Named tuple of forcings for each variable
"""
function forcing_from_boundaries_only(
    grid_ref, 
    boundary_filepath::String, 
    tracers;
    buffer_width::Int=10
)
    grid = grid_ref[]
    
    # Check dimensions
    ds = Dataset(boundary_filepath)
    grid.underlying_grid.Nx == ds.dim["Nx"] &&
        grid.underlying_grid.Ny == ds.dim["Ny"] &&
        grid.underlying_grid.Nz == ds.dim["Nz"] ||
        throw(DimensionMismatch("boundary file dimensions not equal to grid dimensions"))
    
    # Get available variables (include u, v if present)
    forcing_variables_names = (map(String, tracers) ∪ ("u", "v")) ∩ keys(ds)
    close(ds)
    
    # Setup backend
    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    
    # Load boundary forcings
    println("Loading open boundary conditions for: ", forcing_variables_names)
    result = mapreduce(
        var_name -> load_open_boundary_forcing(
            boundary_filepath, 
            var_name, 
            grid, 
            time_indices_in_memory, 
            backend,
            buffer_width
        ),
        merge,
        forcing_variables_names,
    )
    
    return result
end

"""
    forcing_from_rivers_only(
        grid_ref, 
        river_filepath::String, 
        tracers
    )

Create forcing from river inputs only (no boundary conditions).

Convenience function for cases where only river forcing is needed.

# Arguments
- `grid_ref`: Reference to grid
- `river_filepath`: Path to NetCDF file with river forcing
- `tracers`: Tuple of tracer symbols

# Returns
- Named tuple of forcings for each variable
"""
function forcing_from_rivers_only(
    grid_ref, 
    river_filepath::String, 
    tracers
)
    grid = grid_ref[]
    
    # Check dimensions
    ds = Dataset(river_filepath)
    grid.underlying_grid.Nx == ds.dim["Nx"] &&
        grid.underlying_grid.Ny == ds.dim["Ny"] &&
        grid.underlying_grid.Nz == ds.dim["Nz"] ||
        throw(DimensionMismatch("river file dimensions not equal to grid dimensions"))
    
    # Get available variables (check for _flux suffix)
    river_flux_vars = filter(name -> endswith(name, "_flux"), keys(ds))
    river_variables = [replace(name, "_flux" => "") for name in river_flux_vars]
    river_variables = river_variables ∩ map(String, tracers)
    close(ds)
    
    # Setup backend
    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    
    # Load river forcings
    println("Loading river forcing for: ", river_variables)
    result = mapreduce(
        var_name -> load_river_forcing(
            river_filepath, 
            var_name, 
            grid, 
            time_indices_in_memory, 
            backend
        ),
        merge,
        river_variables,
    )
    
    return result
end
