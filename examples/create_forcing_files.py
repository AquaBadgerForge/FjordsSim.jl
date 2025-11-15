"""
Python script to create NetCDF forcing files for FjordsSim

This script demonstrates how to create:
1. Boundary conditions file with buffer zones
2. River forcing file with point sources

Adjust the parameters and data sources to match your specific case.
"""

import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta


class ForcingFileGenerator:
    """Helper class to generate forcing files for FjordsSim"""
    
    def __init__(self, Nx, Ny, Nz, start_date='2024-01-01', time_step_hours=1, num_steps=24):
        """
        Initialize forcing file generator
        
        Parameters:
        -----------
        Nx, Ny, Nz : int
            Grid dimensions (must match FjordsSim grid)
        start_date : str
            Start date for time coordinate
        time_step_hours : float
            Time step in hours
        num_steps : int
            Number of time steps
        """
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.start_date = start_date
        self.time_step_hours = time_step_hours
        self.num_steps = num_steps
        
        # Generate time array (in seconds)
        self.times = np.arange(0, num_steps * time_step_hours * 3600, 
                              time_step_hours * 3600)
    
    def create_boundary_file(self, filepath, tracers_dict, buffer_width=10):
        """
        Create boundary conditions NetCDF file
        
        Parameters:
        -----------
        filepath : str
            Output file path
        tracers_dict : dict
            Dictionary with tracer configurations:
            {
                'T': {'west': 10.0, 'east': 12.0, 'north': 11.0, 'south': 10.5, 'lambda': 1e-4},
                'S': {'west': 34.0, 'east': 34.5, 'north': 34.2, 'south': 34.1, 'lambda': 1e-4},
            }
        buffer_width : int
            Width of buffer zone in grid cells
        """
        print(f"Creating boundary conditions file: {filepath}")
        
        # Create NetCDF file
        ds = nc.Dataset(filepath, 'w', format='NETCDF4')
        
        # Create dimensions
        ds.createDimension('Nx', self.Nx)
        ds.createDimension('Ny', self.Ny)
        ds.createDimension('Nz', self.Nz)
        ds.createDimension('time', None)  # Unlimited
        
        # Create time variable
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = f'seconds since {self.start_date}'
        time_var.long_name = 'time'
        time_var.calendar = 'gregorian'
        time_var[:] = self.times
        
        # Create variables for each tracer
        for tracer_name, tracer_config in tracers_dict.items():
            print(f"  Adding tracer: {tracer_name}")
            
            # Create value variable
            var = ds.createVariable(tracer_name, 'f4', 
                                   ('time', 'Nz', 'Ny', 'Nx'),
                                   fill_value=-999.0)
            var.long_name = f'{tracer_name} at boundaries'
            
            # Create lambda (relaxation) variable
            var_lambda = ds.createVariable(f'{tracer_name}_lambda', 'f4',
                                          ('time', 'Nz', 'Ny', 'Nx'),
                                          fill_value=0.0)
            var_lambda.long_name = f'{tracer_name} relaxation timescale'
            var_lambda.units = '1/s'
            
            # Initialize arrays with fill values
            data = np.full((self.num_steps, self.Nz, self.Ny, self.Nx), -999.0)
            data_lambda = np.zeros((self.num_steps, self.Nz, self.Ny, self.Nx))
            
            # Get boundary values and lambda
            west_val = tracer_config.get('west', -999.0)
            east_val = tracer_config.get('east', -999.0)
            north_val = tracer_config.get('north', -999.0)
            south_val = tracer_config.get('south', -999.0)
            lambda_val = tracer_config.get('lambda', 1e-4)
            
            # Apply to all time steps and depths
            for t in range(self.num_steps):
                for k in range(self.Nz):
                    
                    # West boundary (i = 0 to buffer_width-1)
                    if west_val > -990:
                        for i in range(min(buffer_width, self.Nx)):
                            weight = 1.0 - i / buffer_width
                            data[t, k, :, i] = west_val
                            data_lambda[t, k, :, i] = np.maximum(
                                data_lambda[t, k, :, i], 
                                lambda_val * weight
                            )
                    
                    # East boundary (i = Nx-buffer_width to Nx-1)
                    if east_val > -990:
                        for i in range(max(0, self.Nx - buffer_width), self.Nx):
                            weight = 1.0 - (self.Nx - 1 - i) / buffer_width
                            data[t, k, :, i] = east_val
                            data_lambda[t, k, :, i] = np.maximum(
                                data_lambda[t, k, :, i],
                                lambda_val * weight
                            )
                    
                    # South boundary (j = 0 to buffer_width-1)
                    if south_val > -990:
                        for j in range(min(buffer_width, self.Ny)):
                            weight = 1.0 - j / buffer_width
                            data[t, k, j, :] = south_val
                            data_lambda[t, k, j, :] = np.maximum(
                                data_lambda[t, k, j, :],
                                lambda_val * weight
                            )
                    
                    # North boundary (j = Ny-buffer_width to Ny-1)
                    if north_val > -990:
                        for j in range(max(0, self.Ny - buffer_width), self.Ny):
                            weight = 1.0 - (self.Ny - 1 - j) / buffer_width
                            data[t, k, j, :] = north_val
                            data_lambda[t, k, j, :] = np.maximum(
                                data_lambda[t, k, j, :],
                                lambda_val * weight
                            )
            
            # Write data
            var[:] = data
            var_lambda[:] = data_lambda
            
            # Print statistics
            boundary_cells = np.sum(data_lambda[0, :, :, :] > 0)
            print(f"    Boundary cells with forcing: {boundary_cells}")
            print(f"    Value range: {np.min(data[data > -990]):.2f} to {np.max(data[data > -990]):.2f}")
            print(f"    Lambda range: {np.min(data_lambda[data_lambda > 0]):.2e} to {np.max(data_lambda):.2e}")
        
        # Add global attributes
        ds.title = 'Open Boundary Conditions for FjordsSim'
        ds.created = datetime.now().isoformat()
        ds.buffer_width = buffer_width
        
        ds.close()
        print(f"✓ Boundary file created: {filepath}\n")
    
    def create_river_file(self, filepath, river_locations, river_properties):
        """
        Create river forcing NetCDF file
        
        Parameters:
        -----------
        filepath : str
            Output file path
        river_locations : list of tuples
            List of (i, j, k) grid indices for river mouths
        river_properties : dict
            Dictionary with river properties:
            {
                'discharge': [100.0, 50.0, ...],  # m³/s for each river
                'T': [8.0, 7.5, ...],  # Temperature for each river
                'S': [0.1, 0.2, ...],  # Salinity for each river
            }
        """
        print(f"Creating river forcing file: {filepath}")
        
        # Create NetCDF file
        ds = nc.Dataset(filepath, 'w', format='NETCDF4')
        
        # Create dimensions
        ds.createDimension('Nx', self.Nx)
        ds.createDimension('Ny', self.Ny)
        ds.createDimension('Nz', self.Nz)
        ds.createDimension('time', None)
        
        # Create time variable
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = f'seconds since {self.start_date}'
        time_var.long_name = 'time'
        time_var.calendar = 'gregorian'
        time_var[:] = self.times
        
        # Get list of tracers (exclude 'discharge')
        tracers = [key for key in river_properties.keys() if key != 'discharge']
        discharge_list = river_properties['discharge']
        
        # Validate inputs
        n_rivers = len(river_locations)
        assert len(discharge_list) == n_rivers, "Number of discharge values must match number of rivers"
        for tracer in tracers:
            assert len(river_properties[tracer]) == n_rivers, \
                f"Number of {tracer} values must match number of rivers"
        
        # Create variables for each tracer
        for tracer_name in tracers:
            print(f"  Adding tracer: {tracer_name}")
            
            # Create flux variable
            var_flux = ds.createVariable(f'{tracer_name}_flux', 'f4',
                                        ('time', 'Nz', 'Ny', 'Nx'),
                                        fill_value=0.0)
            var_flux.long_name = f'Volume flux for {tracer_name} forcing'
            var_flux.units = 'm3/s'
            
            # Create value variable
            var_value = ds.createVariable(tracer_name, 'f4',
                                         ('time', 'Nz', 'Ny', 'Nx'),
                                         fill_value=-999.0)
            var_value.long_name = f'{tracer_name} in river water'
            
            # Initialize arrays
            flux_data = np.zeros((self.num_steps, self.Nz, self.Ny, self.Nx))
            value_data = np.full((self.num_steps, self.Nz, self.Ny, self.Nx), -999.0)
            
            # Add river sources
            tracer_values = river_properties[tracer_name]
            for river_idx, (i, j, k) in enumerate(river_locations):
                discharge = discharge_list[river_idx]
                tracer_val = tracer_values[river_idx]
                
                # Apply constant discharge (can add time variation here)
                flux_data[:, k, j, i] = discharge
                value_data[:, k, j, i] = tracer_val
                
                print(f"    River {river_idx+1} at ({i}, {j}, {k}): "
                      f"Q={discharge:.1f} m³/s, {tracer_name}={tracer_val:.2f}")
            
            # Write data
            var_flux[:] = flux_data
            var_value[:] = value_data
            
            # Print statistics
            total_discharge = np.sum(flux_data[0, :, :, :])
            print(f"    Total discharge: {total_discharge:.1f} m³/s")
        
        # Add global attributes
        ds.title = 'River Forcing for FjordsSim'
        ds.created = datetime.now().isoformat()
        ds.number_of_rivers = len(river_locations)
        
        ds.close()
        print(f"✓ River file created: {filepath}\n")


def example_isafjardardjup():
    """Example: Create forcing files for Isafjardardjup"""
    
    # Grid dimensions (must match your FjordsSim grid!)
    Nx, Ny, Nz = 304, 320, 30
    
    # Time setup
    start_date = '2024-01-01'
    time_step_hours = 1.0
    num_steps = 168  # One week
    
    # Initialize generator
    generator = ForcingFileGenerator(Nx, Ny, Nz, start_date, time_step_hours, num_steps)
    
    # ========================================================================
    # 1. Create boundary conditions file
    # ========================================================================
    
    boundary_tracers = {
        'T': {
            'west': 8.0,    # °C - Atlantic water influence
            'east': 8.5,    # °C
            'north': 7.5,   # °C - Colder arctic influence
            'south': 9.0,   # °C - Warmer southern water
            'lambda': 1e-4, # 1/s - relaxation timescale ~2.8 hours
        },
        'S': {
            'west': 34.5,   # psu
            'east': 34.7,   # psu
            'north': 34.2,  # psu
            'south': 34.8,  # psu
            'lambda': 1e-4, # 1/s
        },
    }
    
    generator.create_boundary_file(
        'Isf_boundary_conditions.nc',
        boundary_tracers,
        buffer_width=10
    )
    
    # ========================================================================
    # 2. Create river forcing file
    # ========================================================================
    
    # Define river locations (i, j, k) - surface layer
    # NOTE: These are example locations, adjust to your actual river mouths!
    k_surface = Nz - 1  # Top layer
    river_locations = [
        (50, 100, k_surface),   # River 1
        (75, 150, k_surface),   # River 2
        (120, 200, k_surface),  # River 3
    ]
    
    # Define river properties
    river_properties = {
        'discharge': [100.0, 50.0, 30.0],  # m³/s
        'T': [5.0, 6.0, 4.5],               # °C - cold glacial melt
        'S': [0.1, 0.2, 0.05],              # psu - fresh water
    }
    
    generator.create_river_file(
        'Isf_river_forcing.nc',
        river_locations,
        river_properties
    )
    
    print("=" * 70)
    print("Forcing files created successfully!")
    print("=" * 70)
    print("\nNext steps:")
    print("1. Move the .nc files to your FjordsSim_data directory")
    print("2. Update the file paths in your simulation setup")
    print("3. Verify grid dimensions match your domain")
    print("4. Adjust boundary values and river locations as needed")


def example_custom_domain():
    """Template for creating forcing files for a custom domain"""
    
    # ========================================================================
    # CONFIGURE YOUR DOMAIN
    # ========================================================================
    
    # Grid dimensions - MUST MATCH YOUR FJORDSSIM GRID!
    Nx = 100  # Number of cells in x-direction
    Ny = 120  # Number of cells in y-direction
    Nz = 30   # Number of cells in z-direction (vertical)
    
    # Time configuration
    start_date = '2024-01-01 00:00:00'
    time_step_hours = 1.0  # Hourly data
    num_steps = 24         # 24 hours
    
    # ========================================================================
    # CREATE GENERATOR
    # ========================================================================
    
    generator = ForcingFileGenerator(Nx, Ny, Nz, start_date, time_step_hours, num_steps)
    
    # ========================================================================
    # DEFINE BOUNDARY CONDITIONS
    # ========================================================================
    
    # Adjust these values for your domain!
    boundary_tracers = {
        'T': {
            'west': 10.0,   # Temperature at west boundary (°C)
            'east': 11.0,   # Temperature at east boundary (°C)
            'north': 10.5,  # Temperature at north boundary (°C)
            'south': 10.2,  # Temperature at south boundary (°C)
            'lambda': 1e-4, # Relaxation timescale (1/s)
        },
        'S': {
            'west': 35.0,   # Salinity at west boundary (psu)
            'east': 35.2,   # Salinity at east boundary (psu)
            'north': 34.8,  # Salinity at north boundary (psu)
            'south': 35.1,  # Salinity at south boundary (psu)
            'lambda': 1e-4, # Relaxation timescale (1/s)
        },
    }
    
    # Create boundary file
    generator.create_boundary_file(
        'my_boundary_conditions.nc',
        boundary_tracers,
        buffer_width=10  # Adjust buffer width as needed
    )
    
    # ========================================================================
    # DEFINE RIVER FORCING
    # ========================================================================
    
    # Define river mouth locations (i, j, k grid indices)
    # k = Nz-1 is typically the surface layer
    k_surface = Nz - 1
    
    river_locations = [
        (30, 60, k_surface),   # River 1: example location
        (70, 90, k_surface),   # River 2: example location
        # Add more rivers as needed
    ]
    
    # Define river properties (one value per river)
    river_properties = {
        'discharge': [150.0, 80.0],  # Volume flux in m³/s
        'T': [8.0, 7.0],             # River water temperature (°C)
        'S': [0.5, 0.3],             # River water salinity (psu)
        # Add more tracers as needed (C, NUT, etc.)
    }
    
    # Create river file
    generator.create_river_file(
        'my_river_forcing.nc',
        river_locations,
        river_properties
    )
    
    print("\n" + "=" * 70)
    print("Custom forcing files created!")
    print("=" * 70)


if __name__ == '__main__':
    print("""
    =======================================================================
    FjordsSim Forcing File Generator
    =======================================================================
    
    This script creates NetCDF forcing files for FjordsSim with:
    - Open boundary conditions with buffer zones
    - River point source forcing
    
    Available examples:
    1. example_isafjardardjup() - Example for Isafjardardjup fjord
    2. example_custom_domain()  - Template for custom domains
    
    Usage:
    ------
    Uncomment one of the example functions below to generate files.
    Then adjust the parameters to match your specific case.
    
    =======================================================================
    """)
    
    # Uncomment to run:
    # example_isafjardardjup()
    # example_custom_domain()
