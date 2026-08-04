module BoundaryConditions

export top_bottom_boundary_conditions,
    lateral_tracer_open_boundary_conditions,
    boundary_conditions,
    field_boundary_conditions,
    TopBottomFluxes,
    OpenLateralBoundary

using Oceananigans
using Oceananigans.BoundaryConditions:
    FluxBoundaryCondition,
    NormalFlowBoundaryCondition,
    PerturbationAdvection,
    ValueBoundaryCondition
using Oceananigans.Units: Time
using Oceananigans.Fields: ZeroField
using NumericalEarth.Oceans: u_quadratic_bottom_drag, v_quadratic_bottom_drag, build_tracer_top_bc
using ..Configs: AbstractBoundaryConditionConfig, AbstractForcingConfig
using ..Forcing: boundary_value_time_series
using ..Utils: recursive_merge

""" Return a named tuple with boundary conditions """
function top_bottom_boundary_conditions(; grid, bottom_drag_coefficient)
    top_zonal_momentum_flux = τx = Field{Face,Center,Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center,Face,Nothing}(grid)
    top_ocean_heat_flux = Jᵀ = Field{Center,Center,Nothing}(grid)
    top_salt_flux = Jˢ = Field{Center,Center,Nothing}(grid)

    # The tracer top conditions wrap their flux field in a `FreshwaterExchange` because
    # NumericalEarth's `net_fluxes` reads the exchange back out of the boundary condition to get the
    # freshwater volume flux and its heat content. A bare `FluxBoundaryCondition(Jᵀ)` has nothing to
    # read, so the coupled model raises a `MethodError` the first time it assembles fluxes.
    #
    # `Jʷ` stays zero here: NumericalEarth turns freshwater into a volume change through a
    # free-surface forcing it adds only on a mutable-z grid, and this grid's z faces are static.
    # So the surface T and S fluxes are still exactly `Jᵀ` and `Jˢ`.
    top_freshwater_volume_flux = Jʷ = Field{Center,Center,Nothing}(grid)
    freshwater_heat_content = Field{Center,Center,Nothing}(grid)
    freshwater_salt_content = ZeroField()

    u_bot_bc =
        FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    v_bot_bc =
        FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

    return (
        u = (top = FluxBoundaryCondition(τx), bottom = u_bot_bc),
        v = (top = FluxBoundaryCondition(τy), bottom = v_bot_bc),
        T = (top = build_tracer_top_bc(Jᵀ, Jʷ, freshwater_heat_content, nothing, :T),),
        S = (top = build_tracer_top_bc(Jˢ, Jʷ, freshwater_salt_content, nothing, :S),),
    )
end

const LATERAL_EDGES = (:south, :north, :west, :east)

"""
    x_boundary_tracer_value(j, k, grid, clock, model_fields, parameters)

Discrete-form boundary condition function for a west or east (x-normal) lateral boundary: the value
`parameters.fts` holds at `clock.time`, at the fixed x-index `parameters.boundary_index` and the
boundary's own tangential indices `j`, `k`. `parameters.fts` is the `FieldTimeSeries`
`boundary_value_time_series` returns; `parameters.boundary_index` is `1` for west or `grid.Nx` for
east — a Center-located tracer field has no `Nx + 1` index, unlike the velocity face row
`open_boundary_water!` unmasks.
"""
@inline function x_boundary_tracer_value(j, k, grid, clock, model_fields, parameters)
    return @inbounds parameters.fts[parameters.boundary_index, j, k, Time(clock.time)]
end

"""
    y_boundary_tracer_value(i, k, grid, clock, model_fields, parameters)

As `x_boundary_tracer_value`, for a south or north (y-normal) lateral boundary.
"""
@inline function y_boundary_tracer_value(i, k, grid, clock, model_fields, parameters)
    return @inbounds parameters.fts[i, parameters.boundary_index, k, Time(clock.time)]
end

"""
    tracer_open_boundary_condition(boundary_function, boundary_index, forcing, name, config)

One tracer's lateral open boundary condition at the exact domain edge: a `ValueBoundaryCondition`
wrapping `boundary_function` as a discrete-form condition, with a `PerturbationAdvection` scheme
relaxing towards `name`'s own exterior `FieldTimeSeries` — the same series the interior relaxation
band (`relaxation_lambda`) already reads, reused via `boundary_value_time_series` rather than
rebuilt.

The velocity component normal to this edge is a closed wall (`open_boundary_conditions`), so the
advecting velocity `PerturbationAdvection` reads at this boundary is always
exactly zero, which its own convention classifies as outflow — `inflow_timescale` therefore never
fires here. Both timescales are set to `config.relaxation_timescale` (the boundary-most, untapered
rate from `relaxation_lambda`) rather than exposing a separate, partly-inert knob.
"""
function tracer_open_boundary_condition(boundary_function, boundary_index, forcing, name, config::AbstractForcingConfig)
    haskey(forcing, name) || error(
        "relaxation_edge names a tracer open boundary but the forcing does not prepare `$name`; " *
        "add it to `config.parameters`.",
    )
    scheme = PerturbationAdvection(;
        inflow_timescale = config.relaxation_timescale,
        outflow_timescale = config.relaxation_timescale,
    )
    parameters = (fts = boundary_value_time_series(forcing, name), boundary_index)
    return ValueBoundaryCondition(boundary_function; discrete_form = true, parameters, scheme)
end

"""
    lateral_tracer_open_boundary_conditions(grid, forcing, config::AbstractForcingConfig, tracer_names)

An open `ValueBoundaryCondition` on `config.relaxation_edge` for every tracer in `tracer_names`
(e.g. `:T`, `:S`), relaxing towards the same exterior series the interior relaxation band already
nudges towards. Complements, rather than replaces, that interior forcing: this acts at the exact
boundary column, the band acts across the interior `relaxation_cells`. Errors if `forcing` does not
prepare one of `tracer_names` — naming a relaxation edge implies the forcing should supply every
simulated tracer there.

Velocity is untouched — the edge stays the closed wall `open_boundary_conditions` builds; see that
function's docstring for why a genuinely open velocity boundary is not attempted here.
"""
function lateral_tracer_open_boundary_conditions(grid, forcing, config::AbstractForcingConfig, tracer_names)
    edge = config.relaxation_edge
    edge in LATERAL_EDGES ||
        throw(ArgumentError("relaxation_edge must be one of $LATERAL_EDGES, got :$edge"))

    boundary_function, boundary_index = if edge === :west
        x_boundary_tracer_value, 1
    elseif edge === :east
        x_boundary_tracer_value, grid.Nx
    elseif edge === :south
        y_boundary_tracer_value, 1
    else
        y_boundary_tracer_value, grid.Ny
    end

    return NamedTuple(
        name => (; Symbol(edge) => tracer_open_boundary_condition(boundary_function, boundary_index, forcing, name, config))
        for name in tracer_names
    )
end

"""
    open_boundary_conditions(::Val{edge})

The open boundary condition for the forcing config's `relaxation_edge`, as a nested named tuple
`recursive_merge` can merge into the rest. The open component is the velocity normal to that edge.

Derived rather than configured: the relaxation edge is by construction the edge the regional domain
is open on, so a setup that named it twice could only ever disagree with itself.
"""
open_boundary_conditions(::Val{:south}) = (v = (south = NormalFlowBoundaryCondition(nothing),),)
open_boundary_conditions(::Val{:north}) = (v = (north = NormalFlowBoundaryCondition(nothing),),)
open_boundary_conditions(::Val{:west}) = (u = (west = NormalFlowBoundaryCondition(nothing),),)
open_boundary_conditions(::Val{:east}) = (u = (east = NormalFlowBoundaryCondition(nothing),),)

open_boundary_conditions(::Val{edge}) where {edge} = throw(
    ArgumentError("relaxation_edge must be one of $LATERAL_EDGES, got :$edge"),
)

"""
    TopBottomFluxes(; bottom_drag_coefficient)

The air-sea flux surface and the drag bottom: the boundary-condition piece every setup wants.

Contributes wind stress and heat/salt flux fields at the top and quadratic drag at the bottom, by
way of `top_bottom_boundary_conditions` — see that function for why the tracer top conditions carry
a `FreshwaterExchange` rather than being bare flux conditions.

Only `T` and `S` get a top condition, and that is deliberate rather than an oversight: NumericalEarth
assembles air-sea fluxes for heat and salt specifically, so a third tracer has no exchange to write
into and no `build_tracer_top_bc` method to build one with. A biogeochemical tracer's surface
condition belongs in its own `AbstractBoundaryConditionConfig`.
"""
Base.@kwdef struct TopBottomFluxes <: AbstractBoundaryConditionConfig
    bottom_drag_coefficient::Float64
end

boundary_conditions(config::TopBottomFluxes, grid, forcing, forcing_config, tracers) =
    top_bottom_boundary_conditions(;
        grid = grid,
        bottom_drag_coefficient = config.bottom_drag_coefficient,
    )

"""
    OpenLateralBoundary()

The domain's open edge: a closed wall on the velocity normal to the forcing config's
`relaxation_edge`, plus a relaxing `ValueBoundaryCondition` on every simulated tracer there.

Carries no fields, because everything it needs is already stated in the forcing config — which edge
(`relaxation_edge`) and how fast (`relaxation_timescale`) — and a second copy could only disagree
with the interior relaxation band that reads the same two. A setup with no open edge leaves this out
of its `boundary_conditions` tuple.
"""
struct OpenLateralBoundary <: AbstractBoundaryConditionConfig end

boundary_conditions(::OpenLateralBoundary, grid, forcing, forcing_config, tracers) =
    recursive_merge(
        open_boundary_conditions(Val(forcing_config.relaxation_edge)),
        lateral_tracer_open_boundary_conditions(grid, forcing, forcing_config, tracers),
    )

"""
    field_boundary_conditions(configs, grid, forcing, forcing_config, tracers)

Every boundary-condition config's contribution, merged and materialized into the `NamedTuple` of
`FieldBoundaryConditions` a `HydrostaticFreeSurfaceModel` consumes.

A separate name from `boundary_conditions` on purpose: one function returning a nested named tuple
for a config and a materialized one for a tuple of configs would be two shapes behind one name, and
a caller binding a local named `boundary_conditions` would shadow the hook it is trying to call.

Merging is last-wins, so the tuple's order is its precedence order. An empty tuple is a valid
statement — a model with default boundary conditions everywhere — and yields `(;)`.
"""
function field_boundary_conditions(configs::Tuple, grid, forcing, forcing_config, tracers)
    merged = mapreduce(
        config -> boundary_conditions(config, grid, forcing, forcing_config, tracers),
        recursive_merge,
        configs;
        init = (;),
    )

    return map(sides -> FieldBoundaryConditions(; sides...), merged)
end

end  # module BoundaryConditions
