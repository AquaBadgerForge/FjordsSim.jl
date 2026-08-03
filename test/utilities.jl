# Shared fixtures for the test suite. Included from runtests.jl at top level, before the first
# @testset, because a @testset body is a `let` scope and neither `struct` nor `const` may appear in
# one.
#
# The rule these fixtures exist to enforce: a test asserts *behavior*, never a value a setup chose.
# A setup file is a scientific statement about one fjord, and its author has to stay free to change
# every number in it without a test failing. So nothing here reads a setup, and the tests that used
# to borrow `oslofjorden().simulation_config` build their own config instead.

using Oceananigans:
    CATKEVerticalDiffusivity,
    HydrostaticSphericalCoriolis,
    SeawaterBuoyancy,
    WENO,
    WENOVectorInvariant
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: days, hour, minute, minutes
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature

# Alternative config types

# Used by the "Config extensibility" testset. They exist only to check that `FjordConfig` accepts any
# subtype of the abstract config supertypes, and that new behavior is added by overloading rather
# than by editing FjordSim.

struct SingleColumnGrid <: AbstractGridConfig
    depth::Float64
end

struct MinimalBathymetry <: AbstractBathymetryConfig
    data_root::String
    output_file::String
    plot_file::String
end

struct ConstantForcing <: AbstractForcingConfig
    data_root::String
    output_file::String
    plot_file::String
    temperature::Float64
    rivers::Nothing
end

struct MinimalRivers <: AbstractRiverConfig
    data_root::String
    output_file::String
    relaxation_timescale::Float64
    search_radius::Int
end

struct MinimalAtmosphere <: AbstractAtmosphereConfig
    data_root::String
    output_file::String
    plot_file::String
    output_directory::String
    resolution::Float64
    padding::Float64
end

struct MinimalSimulation <: AbstractSimulationConfig
    results_root::String
    output_file::String
    # `results_path` tags the output with `run_tag`, and `coverage_window` reports the interval the
    # prepare steps pad to, so these two are part of the supertype's field set as much as
    # `results_root` is — an alternative simulation config that omitted them would inherit neither.
    start_date::DateTime
    stop_time::Float64
end

# A stand-in for a coupled-model component that does have a clock, so `rewind_clock!`'s
# has-a-clock branch can be exercised without building a model.
struct ClockHolder{C}
    clock::C
end

# A river config the rivers round-trip drives `add_rivers` with, standing in for a real dataset:
# fixed outlets and a value that counts up in time, so a misplaced write is visible.
struct StubRivers <: AbstractRiverConfig
    data_root::String
    output_file::String
    relaxation_timescale::Float64
    search_radius::Int
    locations::Vector{FjordSim.Forcing.RiverLocation}
    series::Dict{String,Matrix{Float32}}
end

FjordSim.Forcing.river_locations(config::StubRivers) = config.locations
FjordSim.Forcing.river_series(config::StubRivers, times) = config.series

# New behavior for a new grid config, added without touching Grids.jl.
Oceananigans.LatitudeLongitudeGrid(architecture, config::SingleColumnGrid) = LatitudeLongitudeGrid(
    architecture;
    size = (1, 1, 2),
    halo = (1, 1, 1),
    longitude = (10.0, 11.0),
    latitude = (59.0, 60.0),
    z = [-config.depth, -config.depth / 2, 0.0],
)

# One forcing hook overloaded on the new forcing config, without touching Forcing.jl.
FjordSim.Forcing.forcing_variable_names(config::ConstantForcing) = Dict("temperature" => "T")

# A deeply parameterized type, so a frame mentioning it is long enough for `show_compact_error` to
# have something to abbreviate.
struct Nested{A,B} end

nested_failure(::Nested) = error("nested boom")

# The run tag

# `run_tag` is the instant this process started, so a test that wants an exact filename pins it
# rather than deriving one from the config, and puts it back afterwards.
const PINNED_RUN_TAG = "20260803T141530"

function with_run_tag(f, tag = PINNED_RUN_TAG)
    saved = FjordSim.Configs.LAUNCH_TAG[]
    FjordSim.Configs.LAUNCH_TAG[] = tag
    try
        return f()
    finally
        FjordSim.Configs.LAUNCH_TAG[] = saved
    end
end

# Config factories

"""
    test_grid_config(; kwargs...)
    test_bathymetry_config(; data_root = tempdir(), kwargs...)
    test_forcing_config(; data_root = tempdir(), kwargs...)
    test_atmosphere_config(; data_root = tempdir(), kwargs...)

The built-in configs with their required keywords filled in, so a test that cares about one field
names one field. A trailing keyword splat wins over an earlier explicit keyword, so `kwargs`
overrides anything named here.

These fill in *only* the keywords each config requires; every defaulted field stays defaulted, so a
test asserting a default asserts the config's own and not this file's.

A site that asserts an `UndefKeywordError` must keep spelling its keywords out — there, the omission
is the assertion, and a factory would fill in the very field being withheld.
"""
test_grid_config(;
    size = (2, 3, 2),
    halo = (1, 1, 1),
    longitude = (10.0, 12.0),
    latitude = (59.0, 62.0),
    z_faces = collect(range(-10.0 * size[3], 0.0, length = size[3] + 1)),
) = EvenGrid(; size, halo, longitude, latitude, z_faces)

test_bathymetry_config(; data_root = tempdir(), kwargs...) = DybdedataConfig(;
    data_root,
    output_file = "bathymetry.nc",
    plot_file = "bathymetry.png",
    kwargs...,
)

test_forcing_config(; data_root = tempdir(), kwargs...) = NorKystConfig(;
    data_root,
    output_directory = "norkyst",
    parameters = ["temperature", "salinity", "u_eastward", "v_northward"],
    years = [2020],
    kwargs...,
)

test_atmosphere_config(; data_root = tempdir(), kwargs...) = NORA3Config(;
    data_root,
    output_directory = "nora3",
    years = [2020],
    kwargs...,
)

# Simulation config

"""
Where a test simulation config's output resolves to.

Under `tempdir()` rather than the real `~/FjordSim_results`, so every `results_path` assertion is
exact and machine-independent and nothing a test builds can land beside a real run's snapshots.
Nothing creates it: a test that actually builds a simulation passes `results_root` from its own
`mktempdir`, because `coupled_hydrostatic_simulation` mkpaths whatever it is given.
"""
const TEST_RESULTS_ROOT = joinpath(tempdir(), "fjordsim_test_results")

"""
    test_simulation_fields(; kwargs...)

Every one of `SimulationConfig`'s fields, as a `Dict` the constructor can be splatted from.

A `Dict` rather than a `NamedTuple` because the "no field has a default" sweep needs to *remove* one
field at a time, which `delete!(copy(fields), name)` does and a `NamedTuple` cannot.

Three values are not free: `tracers` must name `T` and `S` (`top_bottom_boundary_conditions` and
`attach_writers!` reach for them by name); `buoyancy` must be a `SeawaterBuoyancy` over a
`TEOS10EquationOfState`, the only pair NumericalEarth's `reference_density`/`heat_capacity` have
methods for; and `sea_ice` must not be `nothing`. Everything else is a neutral value chosen here.
"""
function test_simulation_fields(; kwargs...)
    fields = Dict{Symbol,Any}(
        :results_root => TEST_RESULTS_ROOT,
        :output_file => "snapshots_test.nc",
        :architecture => :cpu,          # pinned, never :auto — no test may depend on a GPU
        :buoyancy => SeawaterBuoyancy(equation_of_state = TEOS10EquationOfState()),
        # CATKE plus a horizontal biharmonic operator at their own defaults. The shapes matter and
        # the numbers do not: CATKE is the case where the closure owns a tracer the config does not
        # name, and the pair is what sets the halo floor a coupled build needs.
        :closure => (CATKEVerticalDiffusivity(), HorizontalScalarBiharmonicDiffusivity()),
        :tracer_advection => (T = WENO(), S = WENO()),
        :momentum_advection => WENOVectorInvariant(),
        :tracers => (:T, :S),
        # A literal NamedTuple, so `resolve_initial_conditions` passes it through by identity and
        # nothing built from this fixture also depends on a forcing file having a record at
        # `start_date`. The other two shapes are exercised directly in "initial conditions".
        :initial_conditions => (T = 5.0, S = 33.0),
        :coriolis => HydrostaticSphericalCoriolis(),
        :sea_ice => FreezingLimitedOceanTemperature(),
        :biogeochemistry => nothing,
        :free_surface_cfl => 0.7,
        :bottom_drag_coefficient => 0.003,
        # Noon, matching what `write_prepared_forcing` puts on its time axis, so a config from this
        # fixture passes `validate_time_coverage` against a fixture forcing file.
        :start_date => DateTime(2020, 1, 1, 12),
        :stop_time => 1minute,
        :loops => 1,
        :output_interval => 1hour,
        :progress_interval => 1hour,
        :overwrite_existing => true,
        :checkpoint_interval => 30days,
        :pickup => false,
        :time_step_cfl => 0.1,
        :max_time_step => 3minutes,
        :max_time_step_change => 1.01,
    )

    for (name, value) in kwargs
        haskey(fields, name) || throw(
            ArgumentError(
                "`test_simulation_fields` does not name `$name`; `SimulationConfig`'s fields are " *
                "$(fieldnames(SimulationConfig))",
            ),
        )
        fields[name] = value
    end

    return fields
end

"""
    test_simulation_config(; kwargs...)

A `SimulationConfig` on the test suite's own values, with `kwargs` overriding named fields.

Every duration is written with an `Oceananigans.Units` constant, which is the style the field types
force: `Base.@kwdef` on a parametric struct dispatches rather than converting, so `stop_time = 3600`
is a `MethodError` where `1hour` is fine.

Overriding one of the nine parametric fields yields a different `SimulationConfig` type, and so a
fresh `build_simulation` specialization — a full recompile of the coupled model. Vary only the
`Float64`/`Int`/`Bool` fields inside a testset that builds one.
"""
test_simulation_config(; kwargs...) = SimulationConfig(; test_simulation_fields(; kwargs...)...)

# Grid and bathymetry fixtures

"""
    immersed_test_grid(filepath; size, halo = (1, 1, 1), longitude, latitude, z_faces,
                       bottom_height = ...)

The grid scaffold every pipeline test needs: an `EvenGrid`, its `LatitudeLongitudeGrid`, a bathymetry
file written for it, and the `ImmersedBoundaryGrid` read back out of that file.

Reading the file back rather than immersing the in-memory field is deliberate: it is that round-trip
through `write_bathymetry_file` and `PartialCellBottom` which every caller then exercises, and a mask
computed against a differently built grid would not be testing what the model sees.

`z_faces` defaults to 10 m cells and `bottom_height` to the deepest face everywhere. Everything stays
on the CPU. Returns `(; grid_config, underlying_grid, filepath, grid)`.
"""
function immersed_test_grid(
    filepath;
    size,
    halo = (1, 1, 1),
    longitude = (10.0, 11.0),
    latitude = (59.0, 60.0),
    z_faces = collect(range(-10.0 * size[3], 0.0, length = size[3] + 1)),
    bottom_height = fill(first(z_faces), size[1:2]),
)
    grid_config = test_grid_config(; size, halo, longitude, latitude, z_faces)
    underlying_grid = LatitudeLongitudeGrid(CPU(), grid_config)
    bottom = Field{Center, Center, Nothing}(underlying_grid)
    set!(bottom, bottom_height)
    write_bathymetry_file(filepath, underlying_grid, bottom)

    return (;
        grid_config,
        underlying_grid,
        filepath,
        grid = ImmersedBoundaryGrid(filepath, CPU(), halo),
    )
end

"""
    land_column_test_grid(filepath)

The 5 x 4 x 2 basin the river testsets share: a land column at `i = 2`, so columns 1 and 3 are
coastal while 4 and 5 are open water. The land is off the domain edge because an outlet has to sit
strictly inside the grid to be accepted at all.
"""
function land_column_test_grid(filepath)
    bottom_height = fill(-20.0, (5, 4))
    bottom_height[2, :] .= 0.0

    return immersed_test_grid(
        filepath;
        size = (5, 4, 2),
        longitude = (10.0, 10.5),
        latitude = (59.0, 59.4),
        bottom_height,
    )
end

# Prepared-file fixtures

"""
    write_prepared_forcing(filepath; size, names = ("T", "S", "u", "v"), dates = ...,
                           value = (name, index) -> 1.0f0, lambda = 0.0f0, land = nothing)

A NetCDF file in exactly the layout `prepare_forcing` writes: `Nx`/`Ny`/`Nz` centre dimensions with
`Nx_faces`/`Ny_faces` beside them, one `Float32` variable per name plus its `_lambda` twin, `NaN32`
as the fill value, and a CF-encoded `time`.

The staggering comes from `Forcing.forcing_dimension_names`, the same function `prepare_forcing`
writes with — spelling the convention out a second time is how a change to it earns a passing test
and a file the model cannot read.

- `value(name, index)`: what fills variable `name` at record `index`.
- `lambda`: what fills every `_lambda` twin. `forcing_from_file` dispatches on lambda's *sign*, so a
  test asserting a particular forcing term has to say which.
- `land`: an index tuple set to `NaN32`, the way `prepare_forcing` marks a dry cell.
"""
function write_prepared_forcing(
    filepath;
    size,
    names = ("T", "S", "u", "v"),
    dates = [DateTime(2020, 1, 1, 12), DateTime(2020, 1, 2, 12)],
    value = (name, index) -> 1.0f0,
    lambda = 0.0f0,
    land = nothing,
)
    Nx, Ny, Nz = size

    NCDataset(filepath, "c") do ds
        defDim(ds, "Nx", Nx)
        defDim(ds, "Ny", Ny)
        defDim(ds, "Nz", Nz)
        defDim(ds, "Nx_faces", Nx + 1)
        defDim(ds, "Ny_faces", Ny + 1)
        defDim(ds, "time", length(dates))
        defVar(ds, "time", dates, ("time",))

        for name in names
            dimensions = FjordSim.Forcing.forcing_dimension_names(name)
            shape = ntuple(index -> ds.dim[dimensions[index]], 3)
            values = defVar(ds, name, Float32, dimensions; attrib = ["_FillValue" => NaN32])
            lambdas =
                defVar(ds, name * "_lambda", Float32, dimensions; attrib = ["_FillValue" => NaN32])

            for index in eachindex(dates)
                slab = fill(Float32(value(name, index)), shape)
                isnothing(land) || (slab[land...] = NaN32)
                values[:, :, :, index] = slab
                lambdas[:, :, :, index] = fill(Float32(lambda), shape)
            end
        end
    end

    return filepath
end

"""
    write_nora3_stub(filepath; dates = [0.0, 3600.0], size = (2, 2), temperature = 273.15,
                     FT = Float64)

The smallest file `MultiYearNORA3` will open. The `lon`/`lat`/`time` dimension names are the prepared
file's contract, fixed by `NORA3FieldTimeSeries` reading its axes by them.

`dates` is written as given, which is the one thing callers differ on: a `Vector{Float64}` gives the
raw seconds axis a bare `MultiYearNORA3` accepts, a `Vector{DateTime}` the CF-encoded one a
`Metadatum` needs — `DataWrangling.Metadata` yields a `NORA3Metadatum` only when its date decodes, so
a raw axis would send `size` to the vector method and pass for the wrong reason.
"""
function write_nora3_stub(
    filepath;
    dates = [0.0, 3600.0],
    size = (2, 2),
    temperature = 273.15,
    FT = Float64,
)
    NCDataset(filepath, "c") do ds
        defDim(ds, "lon", size[1])
        defDim(ds, "lat", size[2])
        defDim(ds, "time", length(dates))
        defVar(ds, "time", dates, ("time",))
        defVar(ds, "air_temperature_2m", FT, ("lon", "lat", "time"))[:, :, :] .= temperature
    end

    return filepath
end
