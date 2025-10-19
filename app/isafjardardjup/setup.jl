using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection
using Oceananigans.Architectures: GPU, CPU
using Oceananigans.BuoyancyFormulations: SeawaterBuoyancy
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis, BetaPlane
using Oceananigans.TurbulenceClosures: TKEDissipationVerticalDiffusivity, ScalarDiffusivity, HorizontalScalarBiharmonicDiffusivity
using ClimaOcean
using ClimaOcean.DataWrangling.JRA55
using ClimaOcean.OceanSeaIceModels.InterfaceComputations
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using FjordsSim:
    SetupModel,
    grid_from_bathymetry_file,
    grid_ref,
    forcing_from_file,
    regional_ocean_closure,
    bc_varna_bgh_oxydep,
    bgh_oxydep_boundary_conditions,
    bc_ocean,
    PAR⁰,
    free_surface_default,
    biogeochemistry_LOBSTER,
    biogeochemistry_OXYDEP,
    biogeochemistry_ref

const FT = Oceananigans.defaults.FloatType
const bottom_drag_coefficient = 0.0
const reference_density = 1020

function setup_region(;
    # Grid
    grid_callable = grid_from_bathymetry_file,
    grid_args = (
        arch = GPU(),
        halo = (7, 7, 7),
        filepath = joinpath(homedir(), "FjordsSim_data", "isafjardardjup", "Isf_topo299x320.jld2"),
        latitude = (65.76, 66.399),
        longitude = (-23.492, -22.3),
    ),
    # Buoyancy
    buoyancy = SeawaterBuoyancy(
        equation_of_state = TEOS10EquationOfState(
            reference_density = reference_density
            )
        ),
    # Closure
    closure = (
        TKEDissipationVerticalDiffusivity(minimum_tke = 7e-6),
        HorizontalScalarBiharmonicDiffusivity(ν = 15, κ = 10),
    ),
    # Tracer advection
    tracer_advection = (T = WENO(), S = WENO(), e = nothing, ϵ = nothing),
    # Momentum advection
    momentum_advection = WENOVectorInvariant(),
    # Tracers
    tracers = (:T, :S, :e, :ϵ),
    initial_conditions = (T = 5.0, S = 33.0),
    # Free surface
    free_surface_callable = free_surface_default,
    free_surface_args = (grid_ref,),
    # Coriolis
    coriolis = HydrostaticSphericalCoriolis(),
    # Forcing (disabled)
    forcing_callable = NamedTuple,
    forcing_args = (),
    # Boundary conditions
    bc_callable = bc_ocean,
    bc_args = (grid_ref, bottom_drag_coefficient),
    # Atmosphere
    atmosphere = JRA55PrescribedAtmosphere(
        grid_args.arch,
        latitude = (65.76, 66.399),
        longitude = (-23.492, -22.3)
    ),
    # Ocean emissivity from https://link.springer.com/article/10.1007/BF02233853
    # With suspended matter 0.96 https://www.sciencedirect.com/science/article/abs/pii/0034425787900095
    radiation = ClimaOcean.Radiation(grid_args.arch, ocean_emissivity = 0.96, ocean_albedo = 0.1),
    # coupled model different interfaces
    interfaces = ComponentInterfaces,
    interfaces_kwargs = (
        radiation = radiation,
        freshwater_density = 1000,
        atmosphere_ocean_fluxes = SimilarityTheoryFluxes(),
    ),

    # Biogeochemistry
    biogeochemistry_callable = nothing,
    biogeochemistry_args = (nothing,),

    # Output folder
    results_dir = joinpath(homedir(), "FjordsSim_results", "isafjardardjup"),
)

    return SetupModel(
        grid_callable,
        grid_args,
        grid_ref,
        buoyancy,
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

setup_region_3d() = setup_region()

args_oxydep = (
    initial_photosynthetic_slope = 0.1953 / day, # 1/(W/m²)/s
    Iopt = 80.0, # 50.0,     # (W/m2)
    alphaI = 1.8,   # [d-1/(W/m2)]
    betaI = 5.2e-4, # [d-1/(W/m2)]
    gammaD = 0.71,  # (-)
    Max_uptake = 1.7 / day,  # 1/d 2.0 4 5
    Knut = 1.5,            # (nd) 2.0
    r_phy_nut = 0.10 / day, # 1/d
    r_phy_pom = 0.15 / day, # 1/d
    r_phy_dom = 0.17 / day, # 1/d
    r_phy_het = 0.5 / day,  # 1/d 0.4 2.0
    Kphy = 0.1,             # (nd) 0.7
    r_pom_het = 0.7 / day,  # 1/d 0.7
    Kpom = 2.0,     # (nd)
    Uz = 0.6,       # (nd)
    Hz = 0.5,       # (nd)
    r_het_nut = 0.15 / day,      # 1/d 0.05
    r_het_pom = 0.15 / day,      # 1/d 0.02
    r_pom_nut_oxy = 0.006 / day, # 1/d
    r_pom_dom = 0.05 / day,      # 1/d
    r_dom_nut_oxy = 0.10 / day,  # 1/d
    O2_suboxic = 30.0,    # mmol/m3
    r_pom_nut_nut = 0.010 / day, # 1/d
    r_dom_nut_nut = 0.003 / day, # 1/d
    OtoN = 8.625, # (nd)
    CtoN = 6.625, # (nd)
    NtoN = 5.3,   # (nd)
    NtoB = 0.016, # (nd)
    sinking_speeds = (P = 0.15 / day, HET = 4.0 / day, POM = 10.0 / day),
)

setup_region_3d_OXYDEP() = setup_region(
    tracers = (:T, :S, :e, :ϵ, :C, :NUT, :P, :HET, :POM, :DOM, :O₂),
    initial_conditions = (
        T = 5.0,
        S = 33.0,
        C = 0.0,
        NUT = 10.0,
        P = 0.05,
        HET = 0.01,
        O₂ = 350.0,
        DOM = 1.0,
    ),
    biogeochemistry_callable = biogeochemistry_OXYDEP,
    biogeochemistry_args = (grid_ref, args_oxydep),
    bc_callable = bc_varna_bgh_oxydep,
    bc_args = (grid_ref, bottom_drag_coefficient, biogeochemistry_ref),
    tracer_advection = (
        T = WENO(),
        S = WENO(),
        C = WENO(),
        e = nothing,
        ϵ = nothing,
        NUT = WENO(),
        P = WENO(),
        HET = WENO(),
        POM = WENO(),
        DOM = WENO(),
        O₂ = WENO(),
    ),
)
