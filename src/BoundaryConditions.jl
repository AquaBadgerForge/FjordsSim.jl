module BoundaryConditions

export air_sea_flux_boundary_conditions,
    quadratic_bottom_drag_boundary_conditions,
    top_bottom_boundary_conditions,
    lateral_tracer_open_boundary_conditions,
    boundary_condition_sides,
    field_boundary_conditions,
    AirSeaFluxes,
    QuadraticBottomDrag,
    OpenLateralBoundary,
    MergedBoundaryConditions

using Oceananigans
using Oceananigans.BoundaryConditions:
    FluxBoundaryCondition,
    NormalFlowBoundaryCondition,
    PerturbationAdvection,
    ValueBoundaryCondition
using Oceananigans.Units: Time
using Oceananigans.Fields: ZeroField
using NumericalEarth.Oceans: u_quadratic_bottom_drag, v_quadratic_bottom_drag, build_tracer_top_bc
using ..Configs:
    AbstractBoundaryConditionConfig, AbstractBoundaryConditionSetConfig, AbstractForcingConfig
using ..Forcing: boundary_value_time_series
using ..Utils: recursive_merge

"""
    air_sea_flux_boundary_conditions(grid)

The air-sea exchange surface: wind stress on `u` and `v`, and heat and salt flux on `T` and `S`,
each as a fresh flux field the coupled model writes into.

Only `T` and `S` get a top condition, deliberately: NumericalEarth assembles air-sea fluxes for heat
and salt specifically, so a third tracer has no exchange to write into and no `build_tracer_top_bc`
method to build one with.
"""
function air_sea_flux_boundary_conditions(grid)
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

    return (
        u = (top = FluxBoundaryCondition(τx),),
        v = (top = FluxBoundaryCondition(τy),),
        T = (top = build_tracer_top_bc(Jᵀ, Jʷ, freshwater_heat_content, nothing, :T),),
        S = (top = build_tracer_top_bc(Jˢ, Jʷ, freshwater_salt_content, nothing, :S),),
    )
end

"""
    quadratic_bottom_drag_boundary_conditions(coefficient)

Quadratic drag on `u` and `v` at the bottom, at the given dimensionless drag `coefficient`.

Needs neither the grid nor the forcing — `u_quadratic_bottom_drag` and `v_quadratic_bottom_drag` are
discrete-form functions taking the coefficient as their parameter — which is exactly why it is
separable from the surface half.
"""
quadratic_bottom_drag_boundary_conditions(coefficient) = (
    u = (
        bottom = FluxBoundaryCondition(
            u_quadratic_bottom_drag,
            discrete_form = true,
            parameters = coefficient,
        ),
    ),
    v = (
        bottom = FluxBoundaryCondition(
            v_quadratic_bottom_drag,
            discrete_form = true,
            parameters = coefficient,
        ),
    ),
)

"""
    top_bottom_boundary_conditions(; grid, bottom_drag_coefficient)

The air-sea surface and the drag bottom together, as one nested named tuple.

Kept as the merge of `air_sea_flux_boundary_conditions` and
`quadratic_bottom_drag_boundary_conditions` rather than as their common ancestor: the two halves
share no computation, and a setup usually wants them as two independent configs it can swap one at a
time. This remains the one place the whole `(u, v, T, S)` top-and-bottom shape is stated in a single
expression, which is what the tests assert against.
"""
top_bottom_boundary_conditions(; grid, bottom_drag_coefficient) = recursive_merge(
    air_sea_flux_boundary_conditions(grid),
    quadratic_bottom_drag_boundary_conditions(bottom_drag_coefficient),
)

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
    boundary_condition_sides(config, grid, forcing, forcing_config, tracers)

One `AbstractBoundaryConditionConfig`'s contribution, as a nested named tuple keyed by field and
then by side: `(u = (top = …, bottom = …), T = (top = …,))`.

The per-piece hook. Named for what it returns — the per-field *sides* that
`field_boundary_conditions` materializes — rather than `boundary_conditions`, which is what it used
to be called. That name shadowed `Oceananigans.Fields.boundary_conditions`, so it was the one
FjordSim hook that could not be re-exported and had to be reached as
`FjordSim.BoundaryConditions.boundary_conditions`; and any caller binding a local of that name
shadowed the hook it was trying to call.

Declared with no method, so a subtype that does not implement it fails as a `MethodError` naming the
hook.
"""
function boundary_condition_sides end

"""
    AirSeaFluxes()

The air-sea exchange surface: wind stress on `u` and `v`, heat flux on `T` and salt flux on `S`.

Carries no fields — the flux fields are allocated per grid, and everything that writes into them is
NumericalEarth's business. See `air_sea_flux_boundary_conditions` for why the tracer top conditions
carry a `FreshwaterExchange` rather than being bare flux conditions, and why a third tracer gets no
top condition here. A biogeochemical tracer's surface condition belongs in its own
`AbstractBoundaryConditionConfig`.

Split from the bottom drag it used to be fused with, because the two are independent scientific
statements sharing no computation: a setup can change its drag law, or drop drag entirely, without
restating the surface.
"""
struct AirSeaFluxes <: AbstractBoundaryConditionConfig end

boundary_condition_sides(::AirSeaFluxes, grid, forcing, forcing_config, tracers) =
    air_sea_flux_boundary_conditions(grid)

"""
    QuadraticBottomDrag(; coefficient)

Quadratic drag on `u` and `v` at the bottom.

`coefficient` is the dimensionless drag coefficient the discrete-form drag functions take as their
parameter.
"""
Base.@kwdef struct QuadraticBottomDrag <: AbstractBoundaryConditionConfig
    coefficient::Float64
end

QuadraticBottomDrag(coefficient::Real) = QuadraticBottomDrag(Float64(coefficient))

boundary_condition_sides(config::QuadraticBottomDrag, grid, forcing, forcing_config, tracers) =
    quadratic_bottom_drag_boundary_conditions(config.coefficient)

"""
    OpenLateralBoundary()

The domain's open edge: a closed wall on the velocity normal to the forcing config's
`relaxation_edge`, plus a relaxing `ValueBoundaryCondition` on every simulated tracer there.

Carries no fields, because everything it needs is already stated in the forcing config — which edge
(`relaxation_edge`) and how fast (`relaxation_timescale`) — and a second copy could only disagree
with the interior relaxation band that reads the same two. A setup with no open edge leaves this out
of its `MergedBoundaryConditions`.
"""
struct OpenLateralBoundary <: AbstractBoundaryConditionConfig end

boundary_condition_sides(::OpenLateralBoundary, grid, forcing, forcing_config, tracers) =
    recursive_merge(
        open_boundary_conditions(Val(forcing_config.relaxation_edge)),
        lateral_tracer_open_boundary_conditions(grid, forcing, forcing_config, tracers),
    )

"""
    field_boundary_conditions(config, grid, forcing, forcing_config, tracers)

The boundary conditions a whole `AbstractBoundaryConditionSetConfig` describes, materialized into
the `NamedTuple` a model constructor consumes.

The set-level hook, and a separate name from `boundary_condition_sides` on purpose: one function
returning a nested named tuple for one piece and a materialized one for a set of them would be two
shapes behind a single name.

Declared with no method, so a set config that does not implement it fails as a `MethodError` naming
the hook.
"""
function field_boundary_conditions end

"""
    MergedBoundaryConditions(pieces...)

The built-in `AbstractBoundaryConditionSetConfig`: merge every piece's contribution with
`recursive_merge`, then materialize each field's sides into `FieldBoundaryConditions`.

Merging is last-wins per leaf, so argument order is precedence order. Naming no pieces at all is a
valid statement — a model with default boundary conditions everywhere — and yields `(;)`.

A struct rather than the bare `Tuple` this used to dispatch on. A `Tuple` is not a type any subtype
can specialize, so the merge strategy and the `FieldBoundaryConditions` result shape were fixed for
every model FjordSim could ever assemble; and it said nothing about what the tuple *was*. A model
whose boundary objects are something else, or that wants precedence resolved another way, is now a
sibling of this type.
"""
struct MergedBoundaryConditions{P<:Tuple} <: AbstractBoundaryConditionSetConfig
    pieces::P
end

MergedBoundaryConditions(pieces::AbstractBoundaryConditionConfig...) =
    MergedBoundaryConditions(pieces)

function field_boundary_conditions(
    config::MergedBoundaryConditions,
    grid,
    forcing,
    forcing_config,
    tracers,
)
    merged = mapreduce(
        piece -> boundary_condition_sides(piece, grid, forcing, forcing_config, tracers),
        recursive_merge,
        config.pieces;
        init = (;),
    )

    return map(sides -> FieldBoundaryConditions(; sides...), merged)
end

end  # module BoundaryConditions
