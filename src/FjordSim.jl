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
    AbstractSimulationConfig,
    # generic entry points
    prepare_bathymetry,
    prepare_forcing,
    download_forcing,
    add_rivers,
    download_rivers,
    prepare_atmosphere,
    download_atmosphere,
    build_simulation,
    run_simulation,
    forcing_from_file,
    interpolation_architecture,
    simulation_architecture,
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
    results_path,
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
    prescribed_atmosphere,
    prescribed_radiation,
    atmosphere_date_range,
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
    SimulationConfig,
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
using NCDatasets
using Adapt

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
# Simulations reads the grid back from the processed bathymetry, so it comes after Grids, and
# Setups builds a SimulationConfig, so it comes before Setups.
include("Simulations.jl")
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
using .Simulations
using .Setups
using .CLI

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
