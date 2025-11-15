"""
Example: Using separate boundary and river forcing files

This example demonstrates how to set up a simulation with:
1. Open boundary conditions from one file
2. River forcing from another file
3. Adjustable buffer zone width
"""

using FjordsSim
using Oceananigans
using Oceananigans.Units

# Example for Isafjardardjup with separate forcing files
function example_separate_forcing_isafjardardjup()
    
    # Tracers to use
    tracers = (:T, :S, :e, :ϵ)
    
    # Setup model with separate forcing files
    setup = setup_region(
        # Grid setup
        grid_callable = grid_from_nc,
        grid_args = (
            arch = CPU(),
            halo = (7, 7, 7),
            filepath = joinpath(homedir(), "FjordsSim_data", "isafjardardjup", "Isafjardardjup_bathymetry304x320.nc"),
        ),
        
        # Buoyancy and physics
        buoyancy = SeawaterBuoyancy(
            equation_of_state = TEOS10EquationOfState(reference_density = 1020)
        ),
        
        # Turbulence closure
        closure = (
            TKEDissipationVerticalDiffusivity(minimum_tke = 7e-6),
            HorizontalScalarBiharmonicDiffusivity(ν = 15, κ = 10),
        ),
        
        # Advection schemes
        tracer_advection = (T = WENO(), S = WENO(), e = nothing, ϵ = nothing),
        momentum_advection = WENOVectorInvariant(),
        
        # Tracers and initial conditions
        tracers = tracers,
        initial_conditions = (T = 5.0, S = 33.0),
        
        # Free surface and Coriolis
        free_surface_callable = free_surface_default,
        free_surface_args = (grid_ref,),
        coriolis = HydrostaticSphericalCoriolis(),
        
        # ====================================================================
        # FORCING: Using separate files for boundaries and rivers
        # ====================================================================
        forcing_callable = forcing_from_separate_files,
        forcing_args = (
            grid_ref = grid_ref,
            # Boundary conditions file
            boundary_filepath = joinpath(
                homedir(), 
                "FjordsSim_data", 
                "isafjardardjup", 
                "Isf_boundary_conditions.nc"
            ),
            # River forcing file
            river_filepath = joinpath(
                homedir(), 
                "FjordsSim_data", 
                "isafjardardjup", 
                "Isf_river_forcing.nc"
            ),
            tracers = tracers,
            # Buffer zone: 10 cells (adjust as needed)
            buffer_width = 10,
        ),
        
        # Boundary conditions (surface fluxes, bottom drag)
        bc_callable = bc_ocean,
        bc_args = (grid_ref, 0.0),  # bottom_drag_coefficient = 0.0
        
        # Atmosphere forcing
        atmosphere = JRA55PrescribedAtmosphere(
            CPU(),
            latitude = (65.8, 66.41),
            longitude = (-23.58, -22.3)
        ),
        
        # Radiation
        radiation = ClimaOcean.Radiation(CPU(), ocean_emissivity = 0.96, ocean_albedo = 0.1),
        
        # Coupled model interfaces
        interfaces = ComponentInterfaces,
        interfaces_kwargs = (
            radiation = radiation,
            freshwater_density = 1000,
            atmosphere_ocean_fluxes = SimilarityTheoryFluxes(),
        ),
        
        # No biogeochemistry in this example
        biogeochemistry_callable = nothing,
        biogeochemistry_args = (nothing,),
        
        # Output directory
        results_dir = joinpath(homedir(), "FjordsSim_results", "isafjardardjup_separate_forcing"),
    )
    
    return setup
end

# Example: Boundary forcing only (no rivers)
function example_boundary_only()
    
    tracers = (:T, :S, :e, :ϵ)
    
    setup = setup_region(
        # ... (other parameters same as above) ...
        
        # Only boundary forcing, no rivers
        forcing_callable = forcing_from_boundaries_only,
        forcing_args = (
            grid_ref = grid_ref,
            boundary_filepath = joinpath(
                homedir(), 
                "FjordsSim_data", 
                "myfjord", 
                "boundary_conditions.nc"
            ),
            tracers = tracers,
            buffer_width = 15,  # Wider buffer zone
        ),
        
        # ... (rest of parameters) ...
    )
    
    return setup
end

# Example: River forcing only (no open boundaries)
function example_rivers_only()
    
    tracers = (:T, :S, :e, :ϵ)
    
    setup = setup_region(
        # ... (other parameters same as above) ...
        
        # Only river forcing, no boundaries
        forcing_callable = forcing_from_rivers_only,
        forcing_args = (
            grid_ref = grid_ref,
            river_filepath = joinpath(
                homedir(), 
                "FjordsSim_data", 
                "myfjord", 
                "river_forcing.nc"
            ),
            tracers = tracers,
        ),
        
        # ... (rest of parameters) ...
    )
    
    return setup
end

# Example: Creating a test simulation
function run_example()
    println("Setting up simulation with separate forcing files...")
    
    # Create setup
    setup = example_separate_forcing_isafjardardjup()
    
    # Build coupled simulation
    println("Building coupled simulation...")
    simulation = coupled_hydrostatic_simulation(setup)
    
    # Add some diagnostics
    model = simulation.model.ocean.model
    
    # Print forcing information
    println("\nForcing applied to fields:")
    for (name, forcing) in pairs(model.forcing)
        if !isnothing(forcing) && forcing != NamedTuple()
            println("  - $name")
            if forcing isa FjordsSim.OpenBoundaryForcing
                println("    Type: Open Boundary Condition with buffer zone")
            elseif forcing isa FjordsSim.RiverForcing
                println("    Type: River forcing")
            end
        end
    end
    
    # Set time step
    simulation.Δt = 20.0  # seconds
    
    # Run for a short test period
    simulation.stop_time = 1hour
    
    println("\nRunning test simulation for 1 hour...")
    run!(simulation)
    
    println("✓ Simulation completed successfully!")
    
    # Print final state
    T = model.tracers.T
    S = model.tracers.S
    println("\nFinal state:")
    println("  Temperature: ", extrema(interior(T)), " °C")
    println("  Salinity: ", extrema(interior(S)), " psu")
    
    return simulation
end

# Uncomment to run:
# simulation = run_example()

println("""
Example script loaded. Available functions:

1. example_separate_forcing_isafjardardjup() 
   - Full example with both boundary and river forcing

2. example_boundary_only()
   - Example with only boundary forcing

3. example_rivers_only()
   - Example with only river forcing

4. run_example()
   - Run a complete test simulation

To use, call the desired function:
  setup = example_separate_forcing_isafjardardjup()
  simulation = coupled_hydrostatic_simulation(setup)
  
Or run the test:
  simulation = run_example()
""")
