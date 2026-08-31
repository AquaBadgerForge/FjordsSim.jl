# The `oslofjorden()` setup with the simplest OceanBioME.jl biogeochemistry bolted on: a
# four-compartment nitrogen NPZD model (nutrient, phytoplankton, zooplankton, detritus) in
# mmol N m⁻³.
#
# Nothing in `src/` changes. The point of this file is to show the three ways an out-of-tree config
# extends FjordSim, in increasing order of reach:
#
#   1. **Reuse.** `oslofjorden()` is called, and the pieces that need no change — the grid, the
#      bathymetry, the atmosphere — are passed straight through.
#   2. **A new config subtype with a new hook method.** `NPZDModel`, `NPZDRivers` and
#      `NPZDBoundaries` subtype an `Abstract*Config`, hold the built-in config they extend, and
#      delegate every hook but the one or two they mean to change. That is the documented way to add
#      a source (`docs/adding-a-source.md`) and it is all rivers and the model need.
#   3. **A new method on a *pipeline* generic.** `prepare_forcing` and `prepare_boundaries` are
#      generic functions dispatching on the config type, so `NPZDForcing` and `NPZDBoundaries` add a
#      method that runs the built-in pipeline and then *appends* the biogeochemical variables to the
#      file it wrote. Same copy-then-patch shape `add_rivers` already uses on the forcing file.
#
# (3) is what buys everything else for free. Once `forcing_npzd.nc` carries `N`/`N_lambda` and
# `boundaries_npzd.nc` carries `south_N`, three pieces of FjordSim that were written to be
# name-driven pick the new tracers up with no further code:
#
#   * `state_variables` intersects `model_tracers` with the file's variables, so `FromForcing()`
#     starts N, P, Z and D from the prepared field.
#   * `write_rivers` patches any series variable the file already carries, so river nitrate lands in
#     the same cells and on the same relaxation rates as river temperature and salinity.
#   * `open_tracer_boundary_conditions` radiates *every* model tracer towards its own prepared
#     series, so `OpenLateralBoundaryFromData` needs no change at all — it would in fact **error**
#     without `south_N`, which is why the boundary append is not optional.
#
# The prepared files are new names under the shared `data_root`, so the downloads are reused and
# `oslofjorden()`'s own `forcing.nc` and `boundaries.nc` are left untouched:
#
#   julia --project -m FjordSim prepare_bathymetry  --config oslofjorden          # shared, if not run
#   julia --project -m FjordSim download_forcing    --config oslofjorden          # shared
#   julia --project -m FjordSim download_boundaries --config oslofjorden          # shared
#   julia --project -m FjordSim download_atmosphere --config oslofjorden          # shared
#   julia --project -m FjordSim prepare_atmosphere  --config oslofjorden          # shared
#   julia --project -m FjordSim prepare_forcing     --config examples/oslofjorden_npzd.jl
#   julia --project -m FjordSim add_rivers          --config examples/oslofjorden_npzd.jl
#   julia --project -m FjordSim prepare_boundaries  --config examples/oslofjorden_npzd.jl
#   julia --project -m FjordSim run_simulation      --config examples/oslofjorden_npzd.jl
#
# `add_rivers` needs `NVE_API_KEY`, exactly as `oslofjorden()` does.

using FjordSim
using Oceananigans
using Oceananigans.Units
using OceanBioME
using NCDatasets
using Dates: dayofyear

# --- The literature ---------------------------------------------------------------------------
#
# Every number below is from NIVA report L.NR. 7626-2021, *Overvåking av Ytre Oslofjord 2019-2023.
# Tilførsler og undersøkelser i vannmassene i 2020* — the monitoring year this run simulates, so
# the profiles and the forcing describe the same twelve months. Where a value is an assumption
# rather than a measurement it says so.

# Station OF-1 (outer Oslofjord, 440 m) on 10 March 2020: the pre-bloom profile, the one visit of
# the year when nitrate is closest to conservative, so it is the shape to hang a seasonal cycle on.
# Reported as NO₃+NO₂ in µg N/L; divided by 14.007 to reach mmol N m⁻³.
#
# The surface node is not OF-1's — the report samples that station from 50 m down — but station
# LA-1 at 2 m, which read 94 and 95 µg N/L on 15 January and 12 February 2020.
const NITRATE_PROFILE_DEPTHS = [0.0, 50.0, 75.0, 100.0, 125.0, 150.0, 200.0, 250.0, 300.0, 400.0, 440.0]
const NITRATE_PROFILE_WINTER = [6.7, 6.4, 7.9, 8.6, 9.3, 8.6, 8.6, 8.6, 9.3, 12.8, 13.6]

# Summer surface nitrate, from the same LA-1 visits in June, July and August 2020, every one of
# which came back below the 1 µg N/L detection limit. Not zero, because a nutrient a biogeochemical
# model divides by should not be.
const NITRATE_SUMMER_SURFACE = 0.05

# The drawdown is a growing-season *window* rather than a peak, because that is what the data shows:
# LA-1 read 94 µg N/L at 2 m in January, then below detection on every visit from 18 June to
# 15 September, then 95 µg N/L again on 10 November. A Gaussian centred on midsummer cannot hold a
# four-month floor and still come back up by November — fitted to both ends it sat 15x above the
# observed August surface — so the shape is a product of two logistics: drawdown switching on
# through the spring bloom, and off again when autumn mixing re-supplies the surface.
#
# The onset is the Oslofjord spring bloom, which runs from late February into April; the end is
# mid-October. Both widths are the fortnight-scale of those transitions.
const BLOOM_ONSET_YEARDAY = 90.0
const BLOOM_ONSET_WIDTH = 18.0
const BLOOM_END_YEARDAY = 290.0
const BLOOM_END_WIDTH = 12.0

# How far down the drawdown reaches, as an e-folding scale. 50 m is what the OF-1 series asks for
# and it is deeper than a photic-zone argument would give: on 18 August 2020 the station read
# 49 µg N/L at 50 m against 115 at 100 m, and on 16 September 1.2 µg N/L at 50 m — so the summer
# deficit is still a third of the winter value two euphotic depths down.
const NITRATE_DEPLETION_SCALE = 50.0

# Summer chlorophyll a at 2 m, June-August 2020 (report Table 7): 1.0-4.7 µg/L across the sixteen
# stations, 2.5 µg/L being representative of the outer ones this domain is judged on. Converted at
# a carbon-to-chlorophyll ratio of 50 g C / g Chl and Redfield C:N = 106:16 by mole (5.68 g C / g N),
# so 2.5 µg Chl/L is 22 µg N/L, or 1.6 mmol N m⁻³.
const PHYTOPLANKTON_SUMMER = 1.6

# The report measures chlorophyll only from June to September, so the winter minimum is an
# assumption, not a measurement: 3 % of the summer value, which is the usual overwintering
# inoculum a spin-up needs to have something to grow from.
const PHYTOPLANKTON_WINTER = 0.05

# Both assumptions. No Oslofjord mesozooplankton or detrital nitrogen series was found for 2020, and
# these are the conventional NPZD spin-up fractions — a grazer stock well below its prey, and a
# detritus pool small enough that the first weeks of the run set it rather than this file.
const ZOOPLANKTON_FRACTION = 0.3
const DETRITUS_FRACTION = 0.2

# The euphotic e-folding scale the plankton profiles use. Shallower than the nitrate one because
# biomass tracks light, not the mixed layer.
const EUPHOTIC_SCALE = 25.0

# How well the nitrate profile above actually does, scored against all 23 OF-1 and LA-1 nitrate
# measurements the report publishes for 2020: mean bias +0.9, RMSE 2.2 mmol N m⁻³ on a field that
# spans 0-13. The two ends that matter are close — 6.65 against 6.71 at the surface in January,
# which is the instant `FromForcing()` actually reads, and 8.4 / 8.6 / 12.8 against 8.6 / 8.6 / 12.9
# at 100, 200 and 400 m in the March profile.
#
# Every large residual is subsurface and late-season: 4.1 against 0.1 at 50 m in September, 7.8
# against 2.2 at 100 m the same day. Those are not a seasonal cycle this shape is missing — OF-1's
# 50 m nitrate runs 90, 6.6, 11, 49, 1.2, 22 µg N/L through the year, which is Skagerrak water
# masses arriving and leaving, not a smooth drawdown. Reproducing them is the model's job, driven by
# the open boundary; an analytic profile that tried would be fitting advection with a curve.

# Flow-weighted mean total nitrogen, in mmol N m⁻³, for the rivers the report gauges (its Table 4,
# 2019): the published annual load divided by the published annual water volume, times the nitrate
# fraction below. Glomma is 16 007 t N on 68 724 × 10³ m³ d⁻¹, which is 638 µg TotN/L.
#
# Numedalslågen and Skienselva are outside this domain — both reach the sea west of 10.2°E — and are
# here only because they set `RIVER_NITRATE_DEFAULT`.
const RIVER_TOTAL_NITROGEN = Dict(
    "Glomma (Osterelva)" => 638.0,
    "Glomma (Vesterelva)" => 638.0,
    "Drammenselva" => 463.0,
)

# What fraction of a river's total nitrogen is dissolved inorganic, which is what NPZD's `N` is.
# Measured in the fjord rather than in the rivers, at the two Drammensfjord stations in 2020
# (report Table 6): NO₃+NO₂ over TOTN was 157/293 and 177/347 in summer, and 270/410 at D-2 on
# 14 January. 0.6 is the middle of that, and it is the weakest number on this page.
const RIVER_NITRATE_FRACTION = 0.6

# Every other river in the domain — the nineteen `NVERiversConfig` discovers and this file never
# names — gets the four gauged rivers' combined flow-weighted mean, 506 µg TotN/L: 24 849 t N on
# 134 516 × 10³ m³ d⁻¹. A domain-mean concentration is a much better guess for an ungauged stream
# than Glomma's, and it keeps the load scaling with each stream's own discharge, which is what the
# river forcing already does with volume.
const RIVER_TOTAL_NITROGEN_DEFAULT = 506.0

# µg N/L to mmol N m⁻³, at the atomic weight of nitrogen.
const NITROGEN_MOLAR_MASS = 14.007

# Surface irradiance. `PAR_FRACTION` is the standard photosynthetically active share of shortwave
# energy; `CLEARNESS_INDEX` is southern Norway's annual-mean ratio of surface to top-of-atmosphere
# insolation, from Oslo's ~950 kWh m⁻² yr⁻¹ global horizontal irradiation (108 W m⁻²) against an
# annual-mean top-of-atmosphere ~210 W m⁻² at this latitude.
#
# This is the one place the example fabricates a field FjordSim already has: `prescribed_radiation`
# reads NORA3's downwelling shortwave for the air-sea fluxes, and feeding *that* to the light model
# would be strictly better than a clear-sky formula with a constant cloud factor. It is not done
# here because it would mean reaching inside NumericalEarth's radiation object, which is a bigger
# change than this file is for.
const SOLAR_CONSTANT = 1361.0
const PAR_FRACTION = 0.43
const CLEARNESS_INDEX = 0.5

"""
    seasonal_drawdown(date)

How much of the winter-to-summer nitrate deficit is in place on `date`, from 0 in midwinter to 1
through the productive season.

The one seasonal shape this file has, shared by the nutrient and the plankton profiles so that
biomass appears exactly where nitrate disappears. Sharing it is what makes the plankton profiles a
summer plateau rather than a spring peak, which is as much as the report supports — it measures
chlorophyll only from June to September, and those are the values `PHYTOPLANKTON_SUMMER` is set by.
"""
function seasonal_drawdown(date)
    day = dayofyear(date)
    onset = 1 / (1 + exp(-(day - BLOOM_ONSET_YEARDAY) / BLOOM_ONSET_WIDTH))
    ending = 1 / (1 + exp((day - BLOOM_END_YEARDAY) / BLOOM_END_WIDTH))
    return onset * ending
end

"""
    winter_nitrate(z)

The March 2020 OF-1 nitrate profile at depth `z` (negative down, metres), linearly interpolated
between the reported levels and held constant below the deepest one.
"""
function winter_nitrate(z)
    depth = -z
    depth <= first(NITRATE_PROFILE_DEPTHS) && return first(NITRATE_PROFILE_WINTER)
    depth >= last(NITRATE_PROFILE_DEPTHS) && return last(NITRATE_PROFILE_WINTER)

    upper = searchsortedlast(NITRATE_PROFILE_DEPTHS, depth)
    span = NITRATE_PROFILE_DEPTHS[upper+1] - NITRATE_PROFILE_DEPTHS[upper]
    weight = (depth - NITRATE_PROFILE_DEPTHS[upper]) / span

    return (1 - weight) * NITRATE_PROFILE_WINTER[upper] + weight * NITRATE_PROFILE_WINTER[upper+1]
end

"""
    nitrate_profile(z, date)

Dissolved inorganic nitrogen at depth `z` on `date`, in mmol N m⁻³: the winter profile less a
surface-intensified seasonal deficit.

The deficit is scaled by the *surface* winter value rather than by the local one, so the same
absolute drawdown is removed at every depth the euphotic scale reaches — which is what OF-1 shows.
Floored at the summer surface value so the deep column can never be driven negative by a shallow
argument.
"""
function nitrate_profile(z, date)
    deficit = (first(NITRATE_PROFILE_WINTER) - NITRATE_SUMMER_SURFACE) * seasonal_drawdown(date)
    return max(NITRATE_SUMMER_SURFACE, winter_nitrate(z) - deficit * exp(z / NITRATE_DEPLETION_SCALE))
end

"""
    phytoplankton_profile(z, date)

Phytoplankton nitrogen at depth `z` on `date`, in mmol N m⁻³: the winter stock plus a summer bloom
confined to the euphotic layer, on `seasonal_drawdown`'s own calendar.
"""
phytoplankton_profile(z, date) =
    PHYTOPLANKTON_WINTER +
    (PHYTOPLANKTON_SUMMER - PHYTOPLANKTON_WINTER) * seasonal_drawdown(date) * exp(z / EUPHOTIC_SCALE)

zooplankton_profile(z, date) = ZOOPLANKTON_FRACTION * phytoplankton_profile(z, date)
detritus_profile(z, date) = DETRITUS_FRACTION * phytoplankton_profile(z, date)

# The four profiles, keyed by the tracer name they are written under. One `NamedTuple` because both
# the forcing append and the boundary append want exactly the same set, evaluated the same way.
const NPZD_PROFILES = (
    N = nitrate_profile,
    P = phytoplankton_profile,
    Z = zooplankton_profile,
    D = detritus_profile,
)

"""
    river_nitrate(name)

Dissolved inorganic nitrogen in the river whose outlet is called `name`, in mmol N m⁻³.

Discovered rivers the report does not gauge fall back to the domain-mean concentration, so an
ungauged stream still carries nitrogen in proportion to its own discharge — the volume, and hence
the load, comes from the river forcing's own relaxation rate.
"""
river_nitrate(name) =
    get(RIVER_TOTAL_NITROGEN, name, RIVER_TOTAL_NITROGEN_DEFAULT) * RIVER_NITRATE_FRACTION /
    NITROGEN_MOLAR_MASS

"""
    surface_par(λ, φ, t)

Daily-mean photosynthetically active radiation at the sea surface, W m⁻², at longitude `λ` and
latitude `φ` (degrees) and `t` seconds into the run.

Clear-sky top-of-atmosphere insolation — Cooper's (1969) declination series and the standard
daily-mean integral over the hour angle — scaled by a constant clearness index and the PAR fraction.
At 59.5°N it gives 110 W m⁻² at the solstice and 6 W m⁻² in January, either side of OceanBioME's
own 100 W m⁻² default.

This is a boundary condition on OceanBioME's `PAR` field, so it runs inside a kernel: `ifelse` and
`clamp` rather than branches, and no allocation. `t` is seconds from the run's `start_date`, which
this file sets to 1 January — so day-of-year is `t / 86400` with no offset, and moving `start_date`
means adding one here.
"""
@inline function surface_par(λ, φ, t)
    day = t / 86400
    declination = deg2rad(23.44) * sin(2π * (day + 284) / 365)
    latitude = deg2rad(φ)

    # Sunrise hour angle. Clamped because |tan φ · tan δ| exceeds one inside the polar circles, and
    # a kernel cannot afford `acos` to leave its domain even where this grid never goes.
    hour_angle = acos(clamp(-tan(latitude) * tan(declination), -1, 1))
    eccentricity = 1 + 0.033 * cos(2π * day / 365)

    insolation =
        SOLAR_CONSTANT / π * eccentricity * (
            hour_angle * sin(latitude) * sin(declination) +
            cos(latitude) * cos(declination) * sin(hour_angle)
        )

    return PAR_FRACTION * CLEARNESS_INDEX * max(0, insolation)
end

# --- Writing the profiles into the prepared files -----------------------------------------------

# Matches `FjordSim.Forcing`'s own writer settings, so the appended variables compress and chunk
# like the ones `prepare_forcing` wrote beside them.
const NPZD_DEFLATE_LEVEL = 5

"""
    append_forcing_variables!(filepath, profiles)

Add one `name`/`name_lambda` pair per entry of `profiles` to the prepared forcing file at
`filepath`, horizontally uniform and evaluated on the file's own depth and time axes.

Every relaxation rate is zero, exactly as `prepare_forcing` writes them for T and S: the prepared
forcing is a *state* to start from, not an interior nudging band. The only nonzero rates in the file
are the ones `add_rivers` puts at river mouths afterwards, and it can only put them where the
variable already exists — which is the whole reason this runs before it.

Land is copied from the file's own salinity mask rather than recomputed: `ForcingFromFile` reads a
non-finite value as its land sentinel, and a mask that disagreed with S's would force a cell the
rest of the file calls dry.
"""
function append_forcing_variables!(filepath, profiles)
    NCDataset(filepath, "a") do ds
        z = ds["Nz"][:]
        dates = ds["time"][:]
        shape = (ds.dim["Nx"], ds.dim["Ny"], ds.dim["Nz"])
        wet = map(value -> !ismissing(value) && isfinite(value), ds["S"][:, :, :, 1])
        lambda = zeros(Float32, shape)

        for (name, profile) in pairs(profiles)
            variable = defVar(
                ds, String(name), Float32, ("Nx", "Ny", "Nz", "time");
                chunksizes = [shape[1], shape[2], 1, 1],
                deflatelevel = NPZD_DEFLATE_LEVEL,
                attrib = ["_FillValue" => NaN32],
            )
            rates = defVar(
                ds, String(name) * "_lambda", Float32, ("Nx", "Ny", "Nz", "time");
                chunksizes = [shape[1], shape[2], 1, 1],
                deflatelevel = NPZD_DEFLATE_LEVEL,
                attrib = ["_FillValue" => NaN32],
            )

            @info "Writing $name into $(basename(filepath))"
            slab = Array{Float32}(undef, shape)
            for (index, date) in enumerate(dates)
                column = Float32[profile(depth, date) for depth in z]
                for k in axes(slab, 3)
                    slab[:, :, k] .= ifelse.(wet[:, :, k], column[k], NaN32)
                end
                variable[:, :, :, index] = slab
                rates[:, :, :, index] = lambda
            end
        end
    end

    return filepath
end

"""
    append_boundary_variables!(filepath, edges, profiles)

Add one `<edge>_<name>` variable per entry of `profiles` and per edge to the prepared boundary file
at `filepath`, on the layout `prepare_boundaries` writes: the along-edge axis, then depth, then time.

Uniform along the edge and finite everywhere, including under land. The boundary reader fills
non-finite cells from their nearest finite neighbour anyway (`fill_boundary_gaps!`), so a mask here
would be work with no effect — and the exterior state is a property of the water arriving, not of
which model cells happen to be wet.

`boundary_variable_name` rather than string interpolation, because the `<edge>_<name>` convention is
FjordSim's to state and `boundary_series` reads the file back through the same function.
"""
function append_boundary_variables!(filepath, edges, profiles)
    NCDataset(filepath, "a") do ds
        z = ds["Nz"][:]
        dates = ds["time"][:]
        depths = Float32[profile(depth, date) for profile in values(profiles), depth in z, date in dates]

        for edge in edges
            along = edge in (:south, :north) ? "Nx" : "Ny"
            n_along = ds.dim[along]

            for (row, name) in enumerate(keys(profiles))
                variable_name = boundary_variable_name(edge, String(name))
                @info "Writing $variable_name into $(basename(filepath))"

                variable = defVar(
                    ds, variable_name, Float32, (along, "Nz", "time");
                    deflatelevel = NPZD_DEFLATE_LEVEL,
                    attrib = ["_FillValue" => NaN32],
                )
                variable[:, :, :] = repeat(
                    reshape(depths[row, :, :], 1, length(z), length(dates)), n_along, 1, 1,
                )
            end
        end
    end

    return filepath
end

# --- NPZDModel: a new AbstractCoupledSimulationConfig -------------------------------------------

"""
    NPZDModel(base; tracers, surface_PAR, scale_negatives)

`base`'s model with OceanBioME's NPZD attached.

A config of its own rather than a `biogeochemistry` field on `base`, because `NPZD` needs the grid
and `CoupledHydrostaticSimulation` is built long before there is one — `build_simulation` reads the
grid out of the processed bathymetry. `free_surface` and `model_closure` already solve that by
taking a config and being called with the grid once it exists; `biogeochemistry` has no such hook,
so this type supplies one at the level above, by rebuilding `base` inside `coupled_simulation`.

`tracers` has to be spelled out in full. Oceananigans merges a biogeochemical model's required
tracers into the model's own only when *none* of them is named already
(`has_biogeochemical_tracers`), and NPZD requires `:T` — its growth rate carries a Q₁₀ — which every
FjordSim setup names. All-but-one present is an `ArgumentError`, so the list must be complete.
Naming it here is what `model_tracers` reports, and that is what decides the forcing terms, the open
boundary conditions and the initial conditions, all of which are resolved before the model exists.
"""
struct NPZDModel{M,T<:Tuple,P} <: AbstractCoupledSimulationConfig
    base::M
    tracers::T
    surface_PAR::P
    scale_negatives::Bool
end

NPZDModel(base; tracers, surface_PAR, scale_negatives) =
    NPZDModel(base, Tuple(Symbol.(tracers)), surface_PAR, Bool(scale_negatives))

FjordSim.model_tracers(model::NPZDModel) = model.tracers

"""
    negative_tracer_scaler(scale_negatives)

The `modifiers` an `NPZDModel` asks for: negative scaling over the four biogeochemical tracers, or
`nothing`.

Spelled out rather than left to `NPZD`'s `scale_negatives` keyword, which would include temperature
in the group — see `npzd_model`.

`invalid_fill_value = 0` rather than `ScaleNegativeTracers`' own `NaN` default: a cell whose N+P+Z+D
total goes negative is rare but not impossible even away from `T`, and the thinnest partial-bottom
cells this bathymetry has are numerically fragile enough to hit it within the first few hundred
seconds of a run. `NaN` there is not a local error — every neighbouring cell's tracers are `NaN` by
the next advection step, and the whole active domain follows within one output interval. `0` gives up
exact mass conservation at that one cell in exchange for not losing the run.
"""
negative_tracer_scaler(scale_negatives::Bool) =
    scale_negatives ? ScaleNegativeTracers(Tuple(keys(NPZD_PROFILES)); invalid_fill_value = 0) : nothing

"""
    npzd_model(model, grid)

`model.base` with the biogeochemistry and the tracer list filled in, as a plain
`CoupledHydrostaticSimulation` — so `coupled_simulation` below is a one-line delegation and every
other knob stays stated exactly once, in `base`.

`NPZD` returns an already-wrapped `Biogeochemistry`; passing it to `HydrostaticFreeSurfaceModel` as
it comes is correct, and wrapping it again is not.

`scale_negatives` is worth having on for a run like this one: WENO on a sharp detritus or
phytoplankton front will undershoot, and a negative nutrient is how these models blow up rather than
merely drift. It cannot be `NPZD`'s own `scale_negatives = true`, though, and that is not a matter of
taste. That keyword builds the modifier from `conserved_tracers(bgc)`, which groups the tracers by
element through `required_biogeochemical_tracers` — and `PhytoZoo`'s list is `(:P, :Z, :T)`, `:T`
being *temperature*, which it requires for its Q₁₀ and which is not a nitrogen pool at all. The
scaler would therefore rescale temperature along with the nutrients on any cell holding a negative
tracer, and set it to `NaN` wherever `N + P + Z + T + D` went negative. A `NaN` temperature is a
`NaN` buoyancy, hence a `NaN` velocity, hence a `NaN` `Δt` — which is exactly how this run used to
die, in `calculate_substeps`, before the first hour was out.

So the modifier is built here instead, over the four tracers that really do carry the nitrogen, and
`scale_negatives` stays at its `false` default. `NPZD_PROFILES` names them, so the set cannot drift
apart from the one the forcing and the boundary appends write.
"""
npzd_model(model::NPZDModel, grid) = CoupledHydrostaticSimulation(
    buoyancy = model.base.buoyancy,
    closure = model.base.closure,
    tracer_advection = model.base.tracer_advection,
    momentum_advection = model.base.momentum_advection,
    tracers = model.tracers,
    coriolis = model.base.coriolis,
    sea_ice = model.base.sea_ice,
    biogeochemistry = NPZD(
        grid;
        surface_PAR = model.surface_PAR,
        modifiers = negative_tracer_scaler(model.scale_negatives),
    ),
    free_surface = model.base.free_surface,
    extra_kwargs = model.base.extra_kwargs,
)

FjordSim.coupled_simulation(model::NPZDModel, grid; kwargs...) =
    coupled_simulation(npzd_model(model, grid), grid; kwargs...)

# --- NPZDRivers: a new AbstractRiverConfig ------------------------------------------------------

"""
    NPZDRivers(base; output_file, plot_file)

`base`'s rivers, carrying nitrate as well as temperature and salinity.

Every hook but `river_series` delegates, so the outlets are still the ones `NVERiversConfig`
discovers from NVE's ELVIS network, still placed at the same cells, still relaxed at the same
λ = Q̄/V. Only the values change, and only by one variable.

The alternative was `NVERiversConfig`'s `constants` field, which already writes one fixed value per
variable at every outlet. It is rejected here for what it cannot say: nitrogen in a Norwegian river
is a catchment property — Glomma drains the agricultural Østlandet at 638 µg TotN/L, Drammenselva
463 — and one number for all twenty-one outlets would throw that away.

`standalone` stays `false`: this patches a copy of the prepared forcing, so `FromForcing()` below
still has a full ocean state to start from.
"""
struct NPZDRivers{R} <: AbstractRiverConfig
    base::R
    data_root::String
    output_file::String
    plot_file::String
    relaxation_timescale::Float64
    search_radius::Int
    standalone::Bool
end

NPZDRivers(base; output_file, plot_file) = NPZDRivers(
    base,
    base.data_root,
    String(output_file),
    String(plot_file),
    base.relaxation_timescale,
    base.search_radius,
    base.standalone,
)

FjordSim.river_locations(rivers::NPZDRivers) = river_locations(rivers.base)
FjordSim.river_search_radius(rivers::NPZDRivers) = river_search_radius(rivers.base)
FjordSim.river_minimum_levels(rivers::NPZDRivers) = river_minimum_levels(rivers.base)
FjordSim.river_plume_depth(rivers::NPZDRivers, location::RiverLocation) =
    river_plume_depth(rivers.base, location)
FjordSim.river_lambdas(rivers::NPZDRivers, cells, target_grid) =
    river_lambdas(rivers.base, cells, target_grid)
FjordSim.download_rivers(target_grid, rivers::NPZDRivers) = download_rivers(target_grid, rivers.base)

"""
    river_series(rivers::NPZDRivers, times)

`base`'s series plus an `N` row per outlet, in `river_locations` order.

Constant in time, which is a deliberate limit rather than an oversight: NIVA publishes an annual
*load* and an annual water volume, and their ratio is a flow-weighted mean concentration with no
seasonal information in it. A seasonal cycle here would be invention, and the discharge these
concentrations multiply is already the observed HydAPI series — so the river nitrogen *flux* varies
through the year even though its concentration does not.

Keyed by outlet *name*, which is the only identifier a `RiverLocation` carries — `vassdragsnr`, what
`oslofjorden()` keys its own overrides by, does not survive into one. A name is a weaker key: the
gauged rivers are named only because that setup's overrides state a `name`, and renaming one there
would drop it to `RIVER_TOTAL_NITROGEN_DEFAULT` with nothing to notice. Hence the log line, which
says which outlets were matched and which took the default.
"""
function FjordSim.river_series(rivers::NPZDRivers, times)
    series = river_series(rivers.base, times)
    locations = river_locations(rivers)

    series["N"] = Float32[
        river_nitrate(location.name) for location in locations, _ in times
    ]

    gauged = [location.name for location in locations if haskey(RIVER_TOTAL_NITROGEN, location.name)]
    @info "River nitrate: $(length(gauged)) gauged outlet(s) — $(join(gauged, ", ")) — and " *
          "$(length(locations) - length(gauged)) at the domain mean " *
          "$(round(river_nitrate(""); digits = 1)) mmol N m⁻³"

    return series
end

# --- NPZDForcing: a new AbstractForcingConfig ---------------------------------------------------

"""
    NPZDForcing(base)

`base`'s interior forcing with the four NPZD tracers appended to whatever it prepares.

Only `prepare_forcing` is overloaded. `simulation_forcing` and `forcing_date_range` are the
supertype defaults, because the appended variables are in exactly the layout those readers already
expect — one `name`/`name_lambda` pair on `("Nx", "Ny", "Nz", "time")` — which is the point of
appending rather than writing a second file.

There is no hook inside `prepare_forcing` for a *synthesized* variable, and deliberately so: the
generic pipeline reads every variable it writes out of a downloaded source file, and validates it
against one. That is the right contract for a dataset adapter and the wrong one for a profile out of
a report, so this extends the pipeline at the level where it composes — a second method on
`prepare_forcing` itself — rather than pretending the values came from NorKyst.
"""
struct NPZDForcing{F,R} <: AbstractForcingConfig
    base::F
    data_root::String
    output_file::String
    plot_file::String
    rivers::R
end

NPZDForcing(base) =
    NPZDForcing(base, base.data_root, base.output_file, base.plot_file, base.rivers)

FjordSim.download_forcing(target_grid, config::NPZDForcing) = download_forcing(target_grid, config.base)

"""
    prepare_forcing(target_grid, config::NPZDForcing; coverage, edges)

Regrid `base`'s variables, then append the biogeochemical ones to the file it wrote.

`base` writes to the same path this config resolves — they share `data_root` and `output_file` — so
there is one prepared file, not two, and `add_rivers` copies it complete.
"""
function FjordSim.prepare_forcing(
    target_grid,
    config::NPZDForcing;
    coverage = nothing,
    edges = nothing,
)
    result = prepare_forcing(target_grid, config.base; coverage, edges)
    append_forcing_variables!(result.output_file, NPZD_PROFILES)

    return (; result..., variables = [result.variables; String.(collect(keys(NPZD_PROFILES)))])
end

# --- NPZDBoundaries: a new AbstractBoundaryDataConfig -------------------------------------------

"""
    NPZDBoundaries(base)

`base`'s open-boundary data with the four NPZD tracers appended to the file it prepares.

Not optional. `open_tracer_boundary_conditions` puts a radiating `Value` condition on *every* tracer
`model_tracers` names, and `boundary_series_value` raises rather than defaulting when the prepared
file has no series for one — naming an open edge is a statement that the boundary data supplies the
whole exterior state there. So `OpenLateralBoundaryFromData` below is `oslofjorden()`'s, unchanged,
and this is what makes it work.

`boundary_variable_names` is extended as well as `prepare_boundaries`: the map is what
`boundary_edge_series` reads the file back through, so a variable missing from it would be written
and then ignored. The four new entries are keyed by themselves rather than by a NorKyst variable —
nothing in the download corresponds to them, and nothing looks for them there, because
`prepare_boundaries` only ever considers the source names listed in `parameters`.
"""
struct NPZDBoundaries{B} <: AbstractBoundaryDataConfig
    base::B
    data_root::String
    output_file::String
    plot_file::String
    open_edges::Vector{Symbol}
end

NPZDBoundaries(base) =
    NPZDBoundaries(base, base.data_root, base.output_file, base.plot_file, base.open_edges)

FjordSim.boundary_variable_names(config::NPZDBoundaries) = merge(
    boundary_variable_names(config.base),
    Dict(String(name) => String(name) for name in keys(NPZD_PROFILES)),
)

FjordSim.download_boundaries(target_grid, config::NPZDBoundaries) =
    download_boundaries(target_grid, config.base)

"""
    prepare_boundaries(target_grid, config::NPZDBoundaries; coverage)

Regrid `base`'s exterior state, then append the biogeochemical variables to the file it wrote.

The same shape as `NPZDForcing`'s `prepare_forcing`, and for the same reason: `boundary_source_slab`
can *derive* a variable from others in the download, but `prepared_boundary_variable` still validates
the source name against the downloaded file, so a variable with no source at all has to be added
after the pipeline rather than inside it.
"""
function FjordSim.prepare_boundaries(target_grid, config::NPZDBoundaries; coverage = nothing)
    result = prepare_boundaries(target_grid, config.base; coverage)
    edges = config.open_edges
    append_boundary_variables!(result.output_file, edges, NPZD_PROFILES)

    appended = [
        boundary_variable_name(edge, String(name))
        for edge in edges for name in keys(NPZD_PROFILES)
    ]

    return (; result..., variables = [result.variables; appended])
end

# --- The config itself --------------------------------------------------------------------------

# `oslofjorden()` is the whole science of this run — grid, bathymetry, closure, rivers, boundary
# data, atmosphere, time stepping. Everything below either passes one of its pieces through
# untouched or wraps exactly one of them.
base = oslofjorden()

# The rivers first, since the forcing config holds them. `output_file` is a new name so
# `oslofjorden()`'s own `forcing_rivers_nve.nc` survives for comparison, and so does every other
# file this example writes.
rivers = NPZDRivers(
    base.forcing_config.rivers,
    output_file = "forcing_rivers_npzd.nc",
    plot_file = "forcing_rivers_npzd.png",
)

# A fresh `NorKystConfig` rather than a mutated copy of `base.forcing_config`, because `rivers` is a
# type parameter: swapping it is a new object either way. Same `output_directory`, so the downloaded
# NorKyst months are shared rather than fetched again — `download_forcing` is the one step in this
# pipeline that talks to MET's OPeNDAP server, and the only reason to re-run it is a new year.
forcing = NPZDForcing(
    NorKystConfig(
        data_root = base.forcing_config.data_root,
        output_directory = base.forcing_config.output_directory,
        output_file = "forcing_npzd.nc",
        plot_file = "forcing_npzd.png",
        architecture = base.forcing_config.architecture,
        parameters = base.forcing_config.parameters,
        years = base.forcing_config.years,
        rivers = rivers,
    ),
)

boundaries = NPZDBoundaries(
    NorKystBoundariesConfig(
        data_root = base.boundary_config.data_root,
        output_directory = base.boundary_config.output_directory,
        output_file = "boundaries_npzd.nc",
        plot_file = "boundaries_npzd.png",
        open_edges = base.boundary_config.open_edges,
        margin = base.boundary_config.margin,
        architecture = base.boundary_config.architecture,
        parameters = base.boundary_config.parameters,
        years = base.boundary_config.years,
    ),
)

simulation = base.simulation_config

FjordConfig(
    grid_config = base.grid_config,
    bathymetry_config = base.bathymetry_config,
    forcing_config = forcing,
    boundary_config = boundaries,
    atmosphere_config = base.atmosphere_config,
    simulation_config = SimulationConfig(
        # Its own root. Checkpoints are shared per `results_root`, so pointing this at
        # `oslofjorden()`'s would have the two runs pick each other's up.
        results_root = joinpath(homedir(), "FjordSim_results", "oslofjorden_npzd"),
        architecture = simulation.architecture,
        model = NPZDModel(
            simulation.model,
            # `:T` is in the list because NPZD's growth rate carries a Q₁₀ and therefore requires it;
            # see `NPZDModel`. `tracer_advection` is `oslofjorden()`'s scalar `WENO()`, which was
            # already chosen so that a tracer the setup does not enumerate — CATKE's `e` then, these
            # four now — does not silently fall back to an unbounded centered scheme.
            tracers = (:T, :S, :N, :P, :Z, :D),
            surface_PAR = surface_par,
            scale_negatives = true,
        ),
        # Unchanged. The air-sea flux piece gives a top boundary condition to u, v, T and S only, so
        # the four biogeochemical tracers get Oceananigans' default no-flux surface — right for a
        # first NPZD, which has no gas exchange and no atmospheric nitrogen deposition. The open
        # southern edge is what `NPZDBoundaries` exists to feed.
        boundary_conditions = simulation.boundary_conditions,
        writers = (
            SnapshotWriter(
                name = :ocean,
                output_file = "snapshots_ocean.nc",
                # The four new tracers beside `oslofjorden()`'s five. `PAR` is a biogeochemical
                # *auxiliary* field rather than a model field, so it is not nameable here.
                variables = (:T, :S, :u, :v, :e, :N, :P, :Z, :D),
                interval = 3hour,
                overwrite_existing = true,
            ),
            FieldSnapshotWriter(
                name = :surface,
                output_file = "snapshots_surface.jld2",
                variables = (:η,),
                interval = 3hour,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        callbacks = simulation.callbacks,
        time_stepping = simulation.time_stepping,
        # Unchanged, and it now reads nine fields rather than five: `state_variables` intersects
        # `model_tracers` with the prepared file's own variables, so appending N, P, Z and D to
        # `forcing_npzd.nc` is the entire initial-condition change.
        initial_conditions = FromForcing(),
        start_date = simulation.start_date,
        stop_time = simulation.stop_time,
        loops = simulation.loops,
        pickup = simulation.pickup,
    ),
)
