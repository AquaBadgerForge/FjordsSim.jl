__precompile__(false)  # Disable precompilation to allow method overwriting

module FjordsSim

using OpenSSL_jll
using Oceananigans.Models: HydrostaticFreeSurfaceModel, update_model_field_time_series!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: SplitExplicitFreeSurface
using Oceananigans.Simulations: Simulation
using Oceananigans.Fields: set!
using Oceananigans.TimeSteppers: tick!, Time
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
using Oceananigans.Forcings: DiscreteForcing
using Oceananigans.Units: second, seconds, day, days
using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel
using OceanBioME: LOBSTER

import Oceananigans.Advection: cell_advection_timescale

include("BGCModels/BGCModels.jl")

using .BGCModels

include("utils.jl")
include("grids.jl")
include("turbulence.jl")
include("initial_conditions.jl")
include("boundary_conditions.jl")
include("forcing.jl")
include("atmosphere.jl")
include("output.jl")

grid_ref = Ref{Any}(nothing)
biogeochemistry_ref = Ref{Any}(nothing)

"""A struct with callable and argument objects to construct a coupled hydrostatic simulation. """
mutable struct SetupModel
    grid_callable::Function
    grid_args::NamedTuple
    grid_ref::Ref
    buoyancy::Any
    closure::Any
    tracer_advection::Any
    momentum_advection::Any
    tracers::Tuple
    initial_conditions::NamedTuple
    free_surface_callable::Function
    free_surface_args::Tuple
    coriolis::Any
    forcing_callable::Any
    forcing_args::Any
    bc_callable::Any
    bc_args::Any
    atmosphere::Any
    radiation::Any
    interfaces::Any
    interfaces_kwargs::Any
    biogeochemistry_callable::Any
    biogeochemistry_args::Any
    biogeochemistry_ref::Ref
    results_dir::String
end

"""A constructor to create a folder with results. """
function SetupModel(
    grid_callable,
    grid_args,
    grid_ref,
    buoyance,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    initial_conditions,
    free_surface_callable,
    free_surface_args,
    coriolis,
    forcing_callable,
    forcing_args,
    bc_callable,
    bc_args,
    atmosphere,
    radiation,
    interfaces,
    interfaces_kwargs,
    biogeochemistry_callable,
    biogeochemistry_args,
    biogeochemistry_ref;
    results_dir = joinpath(homedir(), "FjordsSim_results"),
)

    !isdir(results_dir) && mkpath(results_dir)

    SetupModel(
        grid_callable,
        grid_args,
        grid_ref,
        buoyance,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        initial_conditions,
        free_surface_callable,
        free_surface_args,
        coriolis,
        forcing_callable,
        forcing_args,
        bc_callable,
        bc_args,
        atmosphere,
        radiation,
        interfaces,
        interfaces_kwargs,
        biogeochemistry_callable,
        biogeochemistry_args,
        biogeochemistry_ref,
        results_dir,
    )
end

"""
Creates a coupled (precomputed atmosphere + ocean) simulation following ClimaOcean with stubs for the not available components.
Stubs are `nothing`s and for some cases they don't work because some oceananigans callables expect empty tuples.
The hierarchy is as follows:
    ocean_model ->
    ocean_simulation ->
    coupled_model ->
    coupled_simulation
time_step! and update_state! methods are called recursively from coupled_simulation to ocean_model.
"""
function coupled_hydrostatic_simulation(sim_setup::SetupModel)
    grid = sim_setup.grid_callable(sim_setup.grid_args...)
    sim_setup.grid_ref[] = grid
    buoyancy = sim_setup.buoyancy
    closure = sim_setup.closure
    tracer_advection = sim_setup.tracer_advection
    momentum_advection = sim_setup.momentum_advection
    tracers = sim_setup.tracers
    free_surface = sim_setup.free_surface_callable(sim_setup.free_surface_args...)
    coriolis = sim_setup.coriolis
    forcing = sim_setup.forcing_callable(sim_setup.forcing_args...)
    biogeochemistry = safe_execute(sim_setup.biogeochemistry_callable)(sim_setup.biogeochemistry_args...)
    sim_setup.biogeochemistry_ref[] = biogeochemistry
    boundary_conditions = sim_setup.bc_callable(sim_setup.bc_args...)

    println("Start compiling HydrostaticFreeSurfaceModel")
    ## Model
    ocean_model = HydrostaticFreeSurfaceModel(;
        grid,
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

    # ## Simulation
    Δt = 1second
    ocean_sim = Simulation(ocean_model; Δt)
    println("Initialized simulation")

    ## Set initial conditions
    set!(ocean_model; sim_setup.initial_conditions...)
    # Ensure TKE-related tracers have small positive values if present to avoid zero-divide/NaN
    if hasproperty(ocean_model.tracers, :e)
        try
            set!(ocean_model.tracers.e, 1e-8)
        catch err
            @info "Could not initialize tracer e" exception = err
        end
    end
    if hasproperty(ocean_model.tracers, :ϵ)
        try
            set!(ocean_model.tracers.ϵ, 1e-9)
        catch err
            @info "Could not initialize tracer ϵ" exception = err
        end
    end
    println("Set initial conditions")
    try
        T = ocean_model.tracers.T
        S = ocean_model.tracers.S
        println("Initial extrema T: (", minimum(Oceananigans.Fields.interior(T)), ", ", maximum(Oceananigans.Fields.interior(T)), ") °C")
        println("Initial extrema S: (", minimum(Oceananigans.Fields.interior(S)), ", ", maximum(Oceananigans.Fields.interior(S)), ") psu")
    catch err
        @info "Could not print initial extrema for T/S" exception = err
    end

    ## Coupled model / simulation
    sea_ice = FreezingLimitedOceanTemperature(eltype(ocean_sim.model))
    interfaces = sim_setup.interfaces(sim_setup.atmosphere, ocean_sim, sea_ice; sim_setup.interfaces_kwargs...)
    coupled_model = OceanSeaIceModel(ocean_sim, sea_ice; sim_setup.atmosphere, sim_setup.radiation, interfaces)
    println("Initialized coupled model")
    coupled_simulation = Simulation(coupled_model; Δt)
    println("Initialized coupled simulation")

    return coupled_simulation
end

# to allow time step adjusting in OceanSeaIceModel
cell_advection_timescale(model::OceanSeaIceModel) = cell_advection_timescale(model.ocean.model)

free_surface_default(grid_ref) = SplitExplicitFreeSurface(grid_ref[]; cfl = 0.7)

biogeochemistry_LOBSTER(grid_ref) = LOBSTER(;
    grid = grid_ref[],
    surface_photosynthetically_active_radiation = PAR⁰,
    carbonates = true,
    variable_redfield = true,
    oxygen = true,
    scale_negatives = true,
)

biogeochemistry_OXYDEP(grid_ref, args_oxydep) = OXYDEP(;
    grid = grid_ref[],
    args_oxydep...,
    surface_photosynthetically_active_radiation = PAR⁰,
    TS_forced = false,
    Chemicals = false,
    scale_negatives = false,
)

end # module FjordsSim
