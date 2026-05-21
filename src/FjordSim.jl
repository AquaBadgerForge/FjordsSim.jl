module FjordSim

export
    # oceananigans methods
    ImmersedBoundaryGrid,
    # forcings
    forcing_from_file,
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

include("FDatasets.jl")
include("Utils.jl")
include("Bathymetry.jl")
include("Atmospheres/Atmospheres.jl")
include("Forcing.jl")
include("BoundaryConditions.jl")
include("Grids.jl")

using .FDatasets
using .Utils
using .Bathymetry
using .Atmospheres
using .Forcing
using .BoundaryConditions
using .Grids

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

end  # module
