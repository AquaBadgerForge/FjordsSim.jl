module FjordSim

export
    # oceananigans methods
    ImmersedBoundaryGrid,
    LatitudeLongitudeGrid,
    # config supertypes and the setup container
    FjordConfig,
    AbstractGridConfig,
    AbstractBathymetryConfig,
    AbstractForcingConfig,
    AbstractRiverConfig,
    AbstractAtmosphereConfig,
    # generic entry points
    prepare_bathymetry,
    prepare_forcing,
    download_forcing,
    add_rivers,
    download_rivers,
    prepare_atmosphere,
    download_atmosphere,
    forcing_from_file,
    interpolation_architecture,
    plot_bathymetry,
    plot_forcing,
    plot_atmosphere,
    # path resolution, defined on the config supertypes
    bathymetry_path,
    forcing_path,
    forcing_directory,
    river_forcing_path,
    atmosphere_path,
    atmosphere_directory,
    plot_path,
    # extension hooks a new config subtype overloads
    bathymetry_dataset,
    regrid_options,
    forcing_time_steps,
    forcing_source_grid,
    forcing_variable_names,
    forcing_monthly_filename,
    river_locations,
    river_series,
    river_search_radius,
    atmosphere_time_steps,
    atmosphere_source_grid,
    atmosphere_variable_names,
    atmosphere_target_axes,
    ProjectedSourceGrid,
    ProjectedAtmosphereGrid,
    RiverLocation,
    AtmosphereRecord,
    # setups
    fjord_config,
    setup_names,
    oslofjorden,
    drammensfjorden,
    # built-in sources
    EvenGrid,
    DybdedataConfig,
    NorKystConfig,
    OF800RiversConfig,
    NORA3Config,
    geodatabase_path,
    # boundary conditions
    top_bottom_boundary_conditions,
    # simulations
    coupled_hydrostatic_simulation,
    # utils
    recursive_merge,
    progress,
    cell_advection_timescale_coupled_model,
    # atmosphere
    NORA3PrescribedAtmosphere,
    NORA3PrescribedRadiation,
    MultiYearNORA3

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Utils
using NumericalEarth
using NumericalEarth.DataWrangling.JRA55: compute_bounding_nodes, infer_longitudinal_topology
using NCDatasets
using Adapt

import Oceananigans: initialize!
import NumericalEarth.DataWrangling.JRA55: compute_bounding_indices

# Fix NumericalEarth for the custom longitude and latitude
# this is called from set! and uses grid to find the locations,
# which are 1 index more than necessary
function compute_bounding_indices(longitude::Nothing, latitude::Nothing, grid, LX, LY, λc, φc)
    λbounds = compute_bounding_nodes(longitude, grid, LX, λnodes)
    φbounds = compute_bounding_nodes(latitude, grid, LY, φnodes)

    i₁, i₂ = compute_bounding_indices(λbounds, λc)
    j₁, j₂ = compute_bounding_indices(φbounds, φc)
    TX = infer_longitudinal_topology(λbounds)

    # to prevent taking larger than grid areas
    i₁ = (i₂ - i₁ >= grid.Nx) ? (i₂ - grid.Nx + 1) : i₁
    j₁ = (j₂ - j₁ >= grid.Ny) ? (j₂ - grid.Ny + 1) : j₁

    return i₁, i₂, j₁, j₂, TX
end

include("Configs.jl")
include("Datasets.jl")
include("Utils.jl")
# Plotting comes before the pipelines that call it: their setup-level drivers plot as their last
# step, and Plotting itself only needs Configs.
include("Plotting.jl")
include("Bathymetry/Bathymetry.jl")
include("Atmospheres/Atmospheres.jl")
include("Forcing/Forcing.jl")
include("BoundaryConditions.jl")
include("Grids.jl")
# Setups builds every config type, so it comes after all of them — Grids' EvenGrid included.
include("Setups/Setups.jl")
# CLI names every driver and every setup, so it comes last.
include("CLI.jl")

using .Configs
using .Datasets
using .Utils
using .Plotting
using .Bathymetry
using .Atmospheres
using .Forcing
using .BoundaryConditions
using .Grids
using .Setups
using .CLI

function coupled_hydrostatic_simulation(
    grid,
    buoyancy,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    initial_conditions,
    free_surface,
    coriolis,
    forcing,
    boundary_conditions,
    atmosphere,
    downwelling_radiation,
    sea_ice,
    biogeochemistry;
    results_dir = joinpath(homedir(), "FjordSim_results"),
    stop_time = 365days,
)
    isdir(results_dir) || mkpath(results_dir)

    println("Start compiling HydrostaticFreeSurfaceModel")
    ocean_model = HydrostaticFreeSurfaceModel(
        grid;
        buoyancy,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        free_surface,
        coriolis,
        forcing,
        boundary_conditions,
        biogeochemistry,
    )
    println("Done compiling HydrostaticFreeSurfaceModel")
    set!(ocean_model; initial_conditions...)
    Δt = 1second
    ocean_sim = Simulation(ocean_model; Δt, stop_time)
    coupled_model = OceanSeaIceModel(sea_ice, ocean_sim; atmosphere, radiation = downwelling_radiation)
    println("Initialized coupled model")
    coupled_simulation = Simulation(coupled_model; Δt, stop_time)
    return coupled_simulation
end  # function coupled_hydrostatic_simulation

"""
    main(args)

Entry point for `julia --project -m FjordSim SUBCOMMAND --config SETUP`. Returns a process exit
code; see `FjordSim.CLI.USAGE` for the subcommands.

Deliberately *not* exported. Julia's startup runs `Main.main` after a script's body whenever that
binding resolves to an entry point, so exporting this would make every `using FjordSim` in a
script — `test/runtests.jl`, `examples/oslofjord.jl` — run the CLI on the way out.
"""
function main(args)
    return CLI.main(args)
end

# Bare `@main`, applied after the definition: `@main function main(args) ... end` expands to a
# *call*, which would run the CLI while the package precompiles.
@main

end  # module
