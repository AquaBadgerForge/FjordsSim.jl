using Oceananigans.Utils
using Oceananigans.TimeSteppers: implicit_step!
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    TKEDissipationDiffusivityFields,
    FlavorOfTD,
    compute_tke_dissipation_diffusivities!,
    compute_TKEDissipation_diffusivities!,
    substep_tke_dissipation!,
    get_time_step,
    CATKEVerticalDiffusivity,
    CATKEMixingLength,
    CATKEEquation

import Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    top_dissipation_flux, compute_diffusivities!, time_step_tke_dissipation_equations!

@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy)

# is this needed?
@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any,<:Any}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy) +
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy)

# is this needed?
@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any,<:Any,<:Any}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy) +
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy) +
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[3], buoyancy)

# TKE dissipation vertical diffusivity with another closure 
@inline top_dissipation_flux(
    i,
    j,
    grid,
    clock,
    fields,
    parameters,
    closure_tuple::Tuple{<:FlavorOfTD,<:Any},
    buoyancy,
) = top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy)

@inline top_dissipation_flux(
    i,
    j,
    grid,
    clock,
    fields,
    parameters,
    closure_tuple::Tuple{<:Any,<:FlavorOfTD},
    buoyancy,
) = top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy)

# overload to call the proper time_step_tke_dissipation_equations
function compute_diffusivities!(
    diffusivities::TKEDissipationDiffusivityFields,
    closure::FlavorOfTD,
    model;
    parameters = :xyz,
)

    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    clock = model.clock
    top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))

    if isfinite(model.clock.last_Δt) # Check that we have taken a valid time-step first.
        # Compute e at the current time:
        #   * update tendency Gⁿ using current and previous velocity field
        #   * use tridiagonal solve to take an implicit step
        time_step_tke_dissipation_equations!(model, diffusivities, closure)
    end

    # Update "previous velocities"
    u, v, w = model.velocities
    u⁻, v⁻ = diffusivities.previous_velocities
    parent(u⁻) .= parent(u)
    parent(v⁻) .= parent(v)

    launch!(
        arch,
        grid,
        parameters,
        compute_TKEDissipation_diffusivities!,
        diffusivities,
        grid,
        closure,
        velocities,
        tracers,
        buoyancy,
    )

    return nothing
end

# pass proper diffusivities and closure
function time_step_tke_dissipation_equations!(
    model,
    diffusivities::TKEDissipationDiffusivityFields,
    closure::FlavorOfTD,
)

    e = model.tracers.e
    ϵ = model.tracers.ϵ
    arch = model.architecture
    grid = model.grid
    Gⁿe = model.timestepper.Gⁿ.e
    G⁻e = model.timestepper.G⁻.e
    Gⁿϵ = model.timestepper.Gⁿ.ϵ
    G⁻ϵ = model.timestepper.G⁻.ϵ

    κe = diffusivities.κe
    κϵ = diffusivities.κϵ
    Le = diffusivities.Le
    Lϵ = diffusivities.Lϵ
    previous_velocities = diffusivities.previous_velocities
    e_index = findfirst(k -> k == :e, keys(model.tracers))
    ϵ_index = findfirst(k -> k == :ϵ, keys(model.tracers))
    implicit_solver = model.timestepper.implicit_solver

    Δt = model.clock.last_Δt
    Δτ = get_time_step(closure)

    if isnothing(Δτ)
        Δτ = Δt
        M = 1
    else
        M = ceil(Int, Δt / Δτ) # number of substeps
        Δτ = Δt / M
    end

    FT = eltype(grid)

    for m = 1:M # substep
        if m == 1 && M != 1
            χ = convert(FT, -0.5) # Euler step for the first substep
        else
            χ = model.timestepper.χ
        end

        launch!(
            arch,
            grid,
            :xyz,
            compute_tke_dissipation_diffusivities!,
            κe,
            κϵ,
            grid,
            closure,
            model.velocities,
            model.tracers,
            model.buoyancy,
        )

        # Compute the linear implicit component of the RHS (diffusivities, L)
        # and step forward
        launch!(
            arch,
            grid,
            :xyz,
            substep_tke_dissipation!,
            Le,
            Lϵ,
            grid,
            closure,
            model.velocities,
            previous_velocities, # try this soon: model.velocities, model.velocities,
            model.tracers,
            model.buoyancy,
            diffusivities,
            Δτ,
            χ,
            Gⁿe,
            G⁻e,
            Gⁿϵ,
            G⁻ϵ,
        )

    implicit_step!(e, implicit_solver, closure, diffusivities, Val(e_index), model.clock, model.tracers, Δτ)

    implicit_step!(ϵ, implicit_solver, closure, diffusivities, Val(ϵ_index), model.clock, model.tracers, Δτ)
    end

    return nothing
end

function regional_ocean_closure(FT = Oceananigans.defaults.FloatType)
    mixing_length = CATKEMixingLength(
        Cˢ = 1.131,  # Surface distance coefficient for shear length scale
        Cᵇ = 0.28,   # Bottom distance coefficient for shear length scale
        Cˢᵖ = 0.505,  # Sheared convective plume coefficient
        CRiᵟ = 1.02,   # Stability function width
        CRi⁰ = 0.254,  # Stability function lower Ri
        Cʰⁱu = 0.242,  # Shear mixing length coefficient for momentum at high Ri
        Cˡᵒu = 0.361,  # Shear mixing length coefficient for momentum at low Ri
        Cᵘⁿu = 0.370,  # Shear mixing length coefficient for momentum at negative Ri
        Cᶜu = 3.705,  # Convective mixing length coefficient for tracers
        Cᵉu = 0.0,    # Convective penetration mixing length coefficient for tracers
        Cʰⁱc = 0.098,  # Shear mixing length coefficient for tracers at high Ri
        Cˡᵒc = 0.369,  # Shear mixing length coefficient for tracers at low Ri
        Cᵘⁿc = 0.572,  # Shear mixing length coefficient for tracers at negative Ri
        Cᶜc = 4.793,  # Convective mixing length coefficient for tracers
        Cᵉc = 0.112,  # Convective penetration mixing length coefficient for tracers
        Cʰⁱe = 0.548,  # Shear mixing length coefficient for TKE at high Ri
        Cˡᵒe = 7.863,  # Shear mixing length coefficient for TKE at low Ri
        Cᵘⁿe = 1.447,  # Shear mixing length coefficient for TKE at negative Ri
        Cᶜe = 3.642,  # Convective mixing length coefficient for TKE
        Cᵉe = 0.0,    # Convective penetration mixing length coefficient for TKE
    )
    turbulent_kinetic_energy_equation = CATKEEquation(
        CʰⁱD = 0.579, # Dissipation length scale shear coefficient for high Ri
        CˡᵒD = 1.604, # Dissipation length scale shear coefficient for low Ri
        CᵘⁿD = 0.923, # Dissipation length scale shear coefficient for high Ri
        CᶜD = 3.254, # Dissipation length scale convecting layer coefficient
        CᵉD = 0.0,   # Dissipation length scale penetration layer coefficient
        Cᵂu★ = 3.179, # Surface shear-driven TKE flux coefficient
        CᵂwΔ = 0.383, # Surface convective TKE flux coefficient
        Cᵂϵ = 1.0,   # Dissipative near-bottom TKE flux coefficient
    )

    return CATKEVerticalDiffusivity(
        VerticallyImplicitTimeDiscretization(),
        FT;
        mixing_length,
        turbulent_kinetic_energy_equation,
    )
end
