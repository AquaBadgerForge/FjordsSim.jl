using Oceananigans.BoundaryConditions: FluxBoundaryCondition, FieldBoundaryConditions, fill_halo_regions!
using Oceananigans.Fields: Field, set!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryCondition
using Oceananigans.Operators: ℑxyᶜᶠᵃ, ℑxyᶠᶜᵃ
using Oceananigans.Units: days
using Oceananigans.Architectures: GPU, CPU
using Oceananigans.Grids: Center, Face
using ClimaOcean.OceanSimulations: u_quadratic_bottom_drag, v_quadratic_bottom_drag
using OceanBioME: CarbonDioxideGasExchangeBoundaryCondition

function bc_ocean(grid_ref, bottom_drag_coefficient)
    grid = grid_ref[]
    # Set up boundary conditions using Field
    top_zonal_momentum_flux = τx = Field{Face,Center,Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center,Face,Nothing}(grid)
    top_ocean_heat_flux = Jᵀ = Field{Center,Center,Nothing}(grid)
    top_salt_flux = Jˢ = Field{Center,Center,Nothing}(grid)

    # Initialize flux fields to zero and fill halos to prevent garbage values
    set!(τx, 0.0); fill_halo_regions!(τx)
    set!(τy, 0.0); fill_halo_regions!(τy)
    set!(Jᵀ, 0.0); fill_halo_regions!(Jᵀ)
    set!(Jˢ, 0.0); fill_halo_regions!(Jˢ)

    u_bot_bc =
        FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    v_bot_bc =
        FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

    bc_ocean = (
        u = FieldBoundaryConditions(top = FluxBoundaryCondition(τx), bottom = u_bot_bc),
        v = FieldBoundaryConditions(top = FluxBoundaryCondition(τy), bottom = v_bot_bc),
        T = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ)),
        S = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ)),
    )
    return bc_ocean
end

function bc_varna_bgh_oxydep(grid_ref, bottom_drag_coefficient, biogeochemistry_ref)
    Nz = grid_ref[].Nz
    bgc_model = biogeochemistry_ref[]
    bc_varna_tuple = bc_ocean(grid_ref, bottom_drag_coefficient)
    bc_bgh_oxydep_tuple = bgh_oxydep_boundary_conditions(bgc_model, Nz)

    return merge(bc_varna_tuple, bc_bgh_oxydep_tuple)
end

function bc_lobster(grid_ref, bottom_drag_coefficient)
    bc_general = bc_ocean(grid_ref, bottom_drag_coefficient)
    CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition()
    bc_lobster = (DIC = FieldBoundaryConditions(top = CO₂_flux),)
    return merge(bc_general, bc_lobster)
end
