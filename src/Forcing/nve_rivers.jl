# NVE rivers adapter: the source-specific half of the river pipeline, backed by two NVE services.
# It owns a config subtype and the hooks `river_locations`, `river_series`, `download_rivers`,
# `river_plume_depth`, `river_lambdas` and `river_minimum_levels`.
#
# **Where the outlets come from** — NVE's map services, an unauthenticated ArcGIS REST endpoint at
# https://kart.nve.no/enterprise/rest/services. The ELVIS river network gives the river *mouths*
# and the REGINE catchment layer gives each one's size, so a setup states no coordinate at all: see
# `nve_outlets`. Hand-stated outlets remain possible and are what a config gets when it names no
# `minimum_discharge`.
#
# **Where the values come from** — HydAPI (https://hydapi.nve.no/UserDocumentation/), queried per
# station, parameter and year, so the window comes from `years` rather than from whatever a
# published artifact happens to cover. This is the half that needs an API key.
#
# Both are licensed under the Norwegian Licence for Open Government Data (NLOD), which requires
# attribution to NVE. Every HydAPI response carries its own `license` field.

using Downloads
using JSON
using Statistics: mean

const NVE_HYDAPI_URL = "https://hydapi.nve.no/api/v1"

# The environment variable holding the API key. Keys are free and self-service at
# https://hydapi.nve.no/Users — an email address and nothing else.
const NVE_API_KEY_VARIABLE = "NVE_API_KEY"

# HydAPI parameter codes, from `GET /Parameters` (47 of them). These four are the ones a river
# forcing could want; the full list is in the docstring of `NVERiversConfig`.
const NVE_WATER_TEMPERATURE = 1003     # Vanntemperatur, °C
const NVE_DISCHARGE = 1001             # Vannføring, m³/s

# The three resolutions HydAPI offers: 0 instantaneous, 60 hourly, 1440 daily. Daily is what the
# river pipeline wants — `river_series` matches records by calendar date, so a finer axis would
# only ask for the same day repeatedly.
const NVE_DAILY_RESOLUTION = 1440

# The rate limit is 5 requests per sliding second, reported in the `x-rate-limit-*` response
# headers. One pause per request keeps a whole setup's download comfortably inside it without
# having to read those headers back.
const NVE_REQUEST_INTERVAL = 0.25

# Seconds in a mean Julian year, for converting a normal annual runoff (millions of m³/yr) to m³/s.
const NVE_SECONDS_PER_YEAR = 365.25 * 24 * 3600

# NVE's ArcGIS REST endpoint, which needs no key. The `nve.geodataonline.no` host that older NVE
# documentation and Geonorge's own service register still name no longer resolves; this is the live
# one. Geonorge distributes the same ELVIS data, but only as whole-country SOSI/Shape/GML downloads
# plus a WMS — there is no queryable feature service there, which is why this reads NVE directly.
const NVE_MAP_URL = "https://kart.nve.no/enterprise/rest/services"

# ELVIS (Elvenettverk): the complete directed river network, every watercourse as a line with a flow
# direction. Segments share their endpoints exactly, which is what makes a river mouth findable.
const NVE_ELVENETT_LAYER = "Elvenett1/MapServer/2"
const NVE_ELVENETT_FIELDS = "elvenavn,vassdragsnr,vnrnfelt,elveordenstrahler"

# REGINE `Hovedfelt_Nedborfelt_til_hav`: the ~2000 catchments that drain to the sea (or out of the
# country), each with its normal annual runoff. That runoff is what sizes a river with no gauge.
const NVE_CATCHMENT_LAYER = "Nedborfelt2/MapServer/3"
const NVE_CATCHMENT_FIELDS = "vassdragsnr,nedborfeltnavn,areal_km2,qnormal_mm3aar"

# The layers' own `maxRecordCount`. A page carrying exactly this many features sets
# `exceededTransferLimit`, which is the signal to ask for the next `resultOffset`.
const NVE_MAP_PAGE_SIZE = 2000

"""
    NVERiver(; vassdragsnr = "", id = 0, name = "", longitude = NaN, latitude = NaN,
             temperature_station = "", discharge_station = "", mean_discharge = NaN,
             discharge_fraction = 1.0, plume_depth = 0.0)

One river of an `NVERiversConfig`: where its **outlet** is, which HydAPI stations describe it, and
how deep into the water column it is relaxed.

`vassdragsnr` decides which of the two things this is, and it is the only field that does:

- **Non-empty** — an *override* of the discovered outlet whose terminal ELVIS segment carries that
  watercourse number (`002.A21` is Glomma's eastern mouth). Every other field is optional and
  replaces what `nve_outlets` derived: `name = ""` keeps NVE's catchment name, `NaN` coordinates
  keep NVE's river mouth, `id = 0` keeps the discharge rank. An override matching no discovered
  outlet is an error rather than a warning — a typo would otherwise silently drop a river's gauges
  and leave it on the catchment normal, which is a plausible-looking wrong answer.
- **Empty** — a *manually stated* outlet, which must carry `id`, `name`, `longitude` and
  `latitude`. This is what a config with no `minimum_discharge` uses for all its rivers.

`longitude` and `latitude` are the **river mouth**. They are not the gauge's coordinates and cannot
be derived from them: a gauging station sits tens of kilometres upstream at 4 to 100 m elevation,
and HydAPI's own catchment-outlet fields are unusable — `utmEastOutlet`/`utmNorthOutlet` are `null`
for every station checked, and `utmEastInlet` is a placeholder repeated verbatim across unrelated
stations. The station's coordinates are read only to log a sanity check. Getting the mouth from the
river network instead is what `nve_outlets` is for.

Two stations rather than one, because discharge and temperature are frequently not co-located.
Glomma is the case in point: discharge at Solbergfoss (`2.605.0`, 40 464 km², 97 % of the outlet
catchment) and temperature 40 km downstream at Sarpfossen (`2.1087.0`), with no station carrying
both. Either may be `""`:

- `temperature_station = ""` writes `NaN` for that river, which `ForcingFromFile` reads as the
  `-999.0` sentinel and ignores — so the river still freshens its cells, it just does not force
  their temperature. This is the escape for a station whose sensor is not trustworthy, and one is
  needed: Akerselva's 2020 water temperature averages 21.2 °C and peaks at 31.1 °C.
- `discharge_station = ""` skips the discharge download. `river_lambdas` then falls back to the
  catchment normal or to the config's `relaxation_timescale`.

`mean_discharge`, in m³/s, overrides both, for a river whose size is known from somewhere other
than HydAPI. `NaN` means "derive it"; see `river_lambdas`.

`discharge_fraction` multiplies whatever `Q̄` that precedence resolves to, and exists for a river
that reaches the sea through more than one mouth. Glomma is the case: it divides around Kråkerøy
into Østerelva and Vesterelva, ELVIS gives both mouths, and both take Solbergfoss's series scaled
by a fraction. NVE publishes no split, so the fraction is the setup's stated assumption.

`plume_depth` is metres below the surface, read by `river_plume_depth`. `0.0` relaxes the surface
cell alone; `Inf` relaxes the whole wet column, which is the right reading where the model cell is
inside a river bed — Glomma and Drammenselva both are. On an *override*, `0.0` means "keep the
config's `default_plume_depth`", since `0.0` is a real depth rather than a sentinel.

`catchment_discharge` is derived, not stated: `nve_outlets` fills it with what REGINE says the
catchment's normal annual runoff is, and `nve_mean_discharge` falls back to it. A setup that wants
to state a size uses `mean_discharge`, which outranks it.
"""
struct NVERiver
    vassdragsnr::String
    id::Int
    name::String
    longitude::Float64
    latitude::Float64
    temperature_station::String
    discharge_station::String
    mean_discharge::Float64
    discharge_fraction::Float64
    plume_depth::Float64
    catchment_discharge::Float64
end

function NVERiver(;
    vassdragsnr = "",
    id = 0,
    name = "",
    longitude = NaN,
    latitude = NaN,
    temperature_station = "",
    discharge_station = "",
    mean_discharge = NaN,
    discharge_fraction = 1.0,
    plume_depth = 0.0,
    catchment_discharge = NaN,
)
    river = NVERiver(
        String(vassdragsnr),
        Int(id),
        String(name),
        Float64(longitude),
        Float64(latitude),
        String(temperature_station),
        String(discharge_station),
        Float64(mean_discharge),
        Float64(discharge_fraction),
        Float64(plume_depth),
        Float64(catchment_discharge),
    )

    isempty(river.vassdragsnr) &&
        !(isfinite(river.longitude) && isfinite(river.latitude) && !isempty(river.name)) &&
        throw(
            ArgumentError(
                "An NVE river naming no `vassdragsnr` is a manually stated outlet and needs a " *
                "`name`, a `longitude` and a `latitude`. Give it a `vassdragsnr` to make it an " *
                "override of a discovered outlet instead.",
            ),
        )

    return river
end

"""
    NVERiversConfig(; data_root, outlets, years, kwargs...)

River forcing from NVE: river mouths and catchment sizes from the ELVIS/REGINE map services, and
observed daily water temperature and discharge from HydAPI, turned into surface — or full-column —
relaxation towards river properties.

`years` is the calendar years to download and is required, along with `data_root` and `outlets`;
everything else has a default.

# Outlet discovery
`minimum_discharge` decides where the outlets come from, and it is the field to reach for first.

- **`Inf`, the default** — no discovery. `outlets` is the complete list of rivers and every entry
  states its own mouth coordinates.
- **A discharge in m³/s** — every river mouth NVE has in the model domain at or above that size is
  discovered from the map services, and `outlets` becomes a list of *overrides* keyed by
  `vassdragsnr`. No coordinate is stated anywhere. `oslofjorden()` uses `0.5`, which gives 21
  rivers covering 99.6 % of the domain's 1184 m³/s of freshwater.

Discovery is what `nve_outlets` does, and the reason it exists is that hand-copied outlets are
wrong in ways nothing checks. Measured against NVE's own river mouths, the Oslofjord coordinates
this config first carried — lifted from the OF800 dataset's outlet table — put Mosseleva 7 km and
Gjersjøelva 6 km from their rivers, dropped Drammenselva 602 m from the nearest water into a column
at the 2 m `minimum_depth` floor, missed Glomma's western arm entirely, and put Glomma itself on a
shelf beside both of its beds rather than in either.

`margin`, in degrees, is how far past the domain the river network is fetched. It has to be
non-zero: a spatial query returns whole features, but a segment whose *successor* lies wholly
outside the box would then look like a river mouth, so the network is read wider than it is used
and the mouths are filtered back to the domain.

`default_plume_depth` is the `river_plume_depth` a discovered outlet gets when nothing overrides it.

# Requesting data
An API key is required for every **HydAPI** endpoint (no key answers HTTP 401 with an empty body).
Keys are free and self-service at <https://hydapi.nve.no/Users>. `api_key` may state one directly;
left `""`, it is read from the `$NVE_API_KEY_VARIABLE` environment variable. The map services at
`map_url` need no key at all, so a config that discovers its outlets and names no station downloads
its rivers without one.

# Parameters
`GET /Parameters` lists 47. Four could matter to a river forcing, and only two are used:

| code | name | unit | used |
|---|---|---|---|
| 1001 | Vannføring / Discharge | m³/s | yes, to size `river_lambdas` |
| 1003 | Vanntemperatur / Water temperature | °C | yes, as the `T` river value |
| 1002 | Vannhastighet / Water speed | m/s | no — **9 stations nationwide**, and it is the gauge cross-section mean velocity tens of km upstream, not a river-mouth inflow. Rivers enter here as relaxation, not as a momentum or volume flux |
| 1006 | Ledningsevne / Conductivity | µS/cm | no — 16 stations nationwide, and river values (20–200 µS/cm) are fresh anyway |

There is **no salinity parameter**, which is why `constants` writes `S = 0` rather than reading it.

# Finding stations for a new fjord
`GET /Stations?Polygon=POLYGON((lat lon, lat lon, …))` — **latitude first**, unusually for WKT —
returns every station in a box with its full parameter/period/resolution inventory in one call. Two
traps: `Active` is an inverted enum (`0` or omitted is active-only, `1` is all), and a discontinued
station is absent from `/Stations` while `/Series?StationId=…` still has it. Neither endpoint
accepts a comma-separated station list; only `/Observations` does.

# Fields
- `data_root`, `output_file`, `plot_file`: resolved by `river_forcing_path` and `plot_path`. Both
  defaults differ from `OF800RiversConfig`'s so the two sources can be prepared under one
  `data_root` without either overwriting the other.
- `download_directory`: where cached responses go, under `data_root`.
- `outlets`, `years`: the rivers — stated or overridden — and the calendar years to fetch.
- `minimum_discharge`, `margin`, `default_plume_depth`: outlet discovery, above.
- `discovered`: where `nve_outlets` memoises what it derived. Never written by a setup.
- `base_url`, `map_url`, `api_key`, `resolution_time`: the requests.
- `temperature_name`: the FjordSim variable the temperature series is written into.
- `constants`: variable => constant value, written at every plume level of every river.
- `relaxation_timescale`: the fallback λ⁻¹ for a river whose size is unknown.
- `minimum_relaxation_timescale`: the floor on λ⁻¹ derived from discharge. See `river_lambdas`.
- `search_radius`, `minimum_levels`, `standalone`: as `AbstractRiverConfig` documents.

`minimum_levels` defaults to `0`, unlike the value `oslofjorden()` gives its `OF800RiversConfig`,
because a plume is the better answer to a shallow outlet than relocation: see `river_plume_depth`.
"""
Base.@kwdef mutable struct NVERiversConfig <: AbstractRiverConfig
    data_root::String
    output_file::String = "forcing_rivers_nve.nc"
    plot_file::String = "forcing_rivers_nve.png"
    download_directory::String = "nve"
    outlets::Vector{NVERiver}
    years::Vector{Int}
    base_url::String = NVE_HYDAPI_URL
    map_url::String = NVE_MAP_URL
    api_key::String = ""
    resolution_time::Int = NVE_DAILY_RESOLUTION
    minimum_discharge::Float64 = Inf
    margin::Float64 = 0.1
    default_plume_depth::Float64 = 0.0
    discovered::Vector{NVERiver} = NVERiver[]
    temperature_name::String = "T"
    constants::Dict{String,Float64} = Dict("S" => 0.0)
    relaxation_timescale::Float64 = 3600.0
    minimum_relaxation_timescale::Float64 = 600.0
    search_radius::Int = 10
    minimum_levels::Int = 0
    standalone::Bool = false
end

"""
    river_minimum_levels(config::NVERiversConfig)

The wet-level floor for outlet placement, from `config.minimum_levels`. See the generic
`river_minimum_levels` for why the supertype fallback is a plain `0` rather than this field read.
"""
river_minimum_levels(config::NVERiversConfig) = config.minimum_levels

"""
    river_plume_depth(config::NVERiversConfig, location)

How deep `location`'s relaxation reaches, from the matching outlet's `plume_depth`. See the generic
`river_plume_depth` for what the number means and why it is a depth rather than a level count.
"""
function river_plume_depth(config::NVERiversConfig, location::RiverLocation)
    outlet = nve_outlet(config, location.id)
    return outlet.plume_depth
end

"""
    nve_outlet(config, id)

The `NVERiver` with the given `id`. Ids are the `RiverLocation` ids `river_locations` hands out, so
this is how a hook taking a `RiverLocation` gets back to the outlet that produced it.
"""
function nve_outlet(config::NVERiversConfig, id)
    outlets = nve_outlets(config)
    index = findfirst(outlet -> outlet.id == id, outlets)
    isnothing(index) && error("No NVE river with id $id in this config.")
    return outlets[index]
end

"""
    nve_outlets(config)

Every river this config forces, as `NVERiver`s with a mouth, a name and an id.

With no `minimum_discharge` this is just `config.outlets`, and nothing is read. Otherwise it is the
**discovered** set, derived once from the cached map-service responses and memoised into
`config.discovered` — the config is mutable for the same reason `DybdedataConfig` is, so that a
derivation shared by several hooks happens once per run rather than once per river cell.

The derivation, in order:

1. Every ELVIS segment is directed downstream and shares its endpoints exactly with its neighbours,
   so a **river mouth is an end-vertex that is not the start-vertex of any segment**. On the
   Oslofjord domain, 16 365 segments give 463 such terminals.
2. Keep the terminals whose catchment (`vnrnfelt`) is one REGINE calls sea-draining — 57 of the 463
   — since the rest are network endpoints at lakes and at the edge of the mapped area.
3. Drop the ones outside the domain itself — 49 survive on Oslofjorden. `margin` had the network
   fetched wider than the domain precisely so this step has something to cut, and cutting it here is
   also what removes the inland terminals that a large catchment has upstream of the fjord.
4. Drop `Q̄ < minimum_discharge`, where `Q̄` is the catchment's normal annual runoff in m³/s — 21
   survive at Oslofjorden's 0.5 m³/s.
5. Keep **one mouth per catchment**, the one whose segment has the highest Strahler order, ties
   broken by `vassdragsnr` so the choice is deterministic. Three Oslofjord catchments offer two
   terminals and in each case one of them is an inland artifact tens of kilometres away. A
   catchment offering more than one *inside* the domain is warned about rather than silently
   halved, because that is a genuine distributary and the setup may want both.
6. Sort by descending `Q̄`, `vassdragsnr` behind it, and assign ids by that rank — so an id also
   states how large a river is, and two rivers of equal size keep the same ids run after run.
7. Apply `config.outlets` as overrides, matched on the terminal segment's `vassdragsnr`.

Note that a real distributary survives step 5 whenever NVE files it under its own catchment, which
is what happens to Glomma's western arm: ELVIS puts Vesterelva in the small Seutelva catchment, so
it is discovered as a river in its own right and the setup only has to correct its name and its
share of Glomma's discharge.

Nothing here tests the model's land mask. It does not need to: an outlet that survives all seven
steps but is nonetheless inland is dropped by `river_cells`, which already refuses to place a river
with no coastal water cell within `river_search_radius` of it and says so.
"""
function nve_outlets(config::NVERiversConfig)
    isfinite(config.minimum_discharge) || return config.outlets
    isempty(config.discovered) || return config.discovered

    catchments, domain = nve_read_catchments(config)
    segments, _ = nve_read_network(config)

    outlets = nve_discover_outlets(segments, catchments, domain, config)
    isempty(outlets) && error(
        "No NVE river mouth in the domain reaches the config's `minimum_discharge` of " *
        "$(config.minimum_discharge) m³/s. The domain holds $(length(catchments)) sea-draining " *
        "catchments; lower the threshold, or state the outlets directly.",
    )

    config.discovered = nve_apply_overrides(outlets, config)

    return config.discovered
end

"""
    nve_network_path(config)
    nve_catchment_path(config)

Where the two cached map-service responses live, under `nve_download_directory(config)`.

Neither is a verbatim mirror. The network file keeps each segment's **two endpoints** and its
attributes and throws the intermediate vertices away, which is a quarter of the bytes and all of
the information a river mouth needs; the catchment file keeps no geometry at all. Both record the
domain they were fetched for, so a grid change re-downloads instead of silently reusing another
domain's network — see `nve_cache_covers`.
"""
nve_network_path(config::NVERiversConfig) =
    joinpath(nve_download_directory(config), "elvis_elvenett.json")

nve_catchment_path(config::NVERiversConfig) =
    joinpath(nve_download_directory(config), "regine_hovedfelt.json")

"""
    nve_domain(target_grid)
    nve_grown_domain(domain, margin)

The model domain as ArcGIS wants an envelope — `(west, south, east, north)` — and that box grown by
`margin` degrees on every side.

The grown box is what the network is *queried* with and the plain one is what its terminals are
*filtered* to, and the two must differ. A spatial query returns whole features, so a segment
crossing the box edge comes back intact — but its downstream neighbour, lying wholly outside, does
not, which leaves the crossing segment looking like a river mouth. Reading wider than the domain
and cutting back to it is what makes those false mouths land outside the filter.

Measured on Oslofjorden: querying the domain alone returns 14 132 segments with 481 terminals,
grown by 0.1° it returns 16 365 with 463. Some of those 481 were segments the box had cut off from
their own downstream neighbour, and the margin gives them one back; the margin ring then introduces
false terminals of its own at *its* edge, which is what step 3 of `nve_outlets` cuts.
"""
function nve_domain(target_grid)
    west, east = x_domain(target_grid)
    south, north = y_domain(target_grid)
    return (west, south, east, north)
end

nve_grown_domain(domain, margin) =
    (domain[1] - margin, domain[2] - margin, domain[3] + margin, domain[4] + margin)

"""
    nve_map_request(config, url, description)

GET one ArcGIS REST query and parse it. No API key: NVE's map services are open, unlike HydAPI.

The service reports a rejected query as **HTTP 200 carrying an `error` object**, not as a status
code, so the body is checked as well as the transport — the same shape of trap as HydAPI's
zero-observation 200, and it would otherwise surface much later as an empty river list.
"""
function nve_map_request(config::NVERiversConfig, url, description)
    @info "Downloading $description from NVE's map services"

    filepath = try
        Downloads.download(url)
    catch exception
        exception isa Downloads.RequestError || rethrow()
        status = exception.response.status
        error(
            "NVE map service request for $description failed" *
            (iszero(status) ? "" : " with HTTP $status") *
            ": $(exception.message). Request was $url",
        )
    end

    response = nothing
    try
        response = JSON.parsefile(filepath)
    finally
        rm(filepath; force = true)
    end

    failure = get(response, "error", nothing)
    isnothing(failure) || error(
        "NVE's map service rejected the request for $description: " *
        "$(get(failure, "message", failure)). Request was $url",
    )

    return response
end

"""
    nve_map_features(config, layer, fields, domain, geometry, description)

Every feature of one layer intersecting `domain`, paging until the service stops asking for more.

A page carrying `maxRecordCount` features sets `exceededTransferLimit`; the next page is the same
query at the next `resultOffset`. The Oslofjord river network is 14 132 segments, so eight pages.
`geometryPrecision=6` is about 0.1 m and halves the payload without moving a node — the mouths it
finds are the same set the full-precision geometry gives, and the same set the native EPSG:25833
geometry gives at centimetre precision.

The Oslofjord river network is 16 365 segments over the domain grown by `margin`, so nine pages.
"""
function nve_map_features(config::NVERiversConfig, layer, fields, domain, geometry, description)
    features = Any[]
    offset = 0

    while true
        url = string(
            config.map_url, "/", layer, "/query",
            "?where=1%3D1",
            "&geometry=", join(domain, ","),
            "&geometryType=esriGeometryEnvelope",
            "&inSR=4326",
            "&spatialRel=esriSpatialRelIntersects",
            "&outFields=", fields,
            "&returnGeometry=", geometry,
            "&outSR=4326",
            "&geometryPrecision=6",
            "&f=geojson",
            "&resultOffset=", offset,
            "&resultRecordCount=", NVE_MAP_PAGE_SIZE,
        )

        page = nve_map_request(config, url, "$description from record $offset")
        append!(features, get(page, "features", ()))

        get(page, "exceededTransferLimit", false) === true || break
        offset += NVE_MAP_PAGE_SIZE
    end

    return features
end

"""
    nve_segment_endpoints(feature)

The first and last vertex of one ELVIS segment, as `(longitude, latitude)` pairs, or `nothing` for
a feature carrying no usable geometry.

ELVIS lines are directed downstream, so the last vertex is the downstream end — which is the whole
basis of `nve_discover_outlets`. A `MultiLineString` is flattened across its parts rather than
rejected: the parts of one feature are consecutive, so its first and last vertex are still its two
ends.
"""
function nve_segment_endpoints(feature)
    geometry = get(feature, "geometry", nothing)
    isnothing(geometry) && return nothing

    coordinates = get(geometry, "coordinates", nothing)
    isnothing(coordinates) && return nothing

    points =
        get(geometry, "type", "") == "MultiLineString" ?
        collect(Iterators.flatten(coordinates)) : coordinates
    isempty(points) && return nothing

    return nve_node(first(points)), nve_node(last(points))
end

"""
    nve_node(point)

One river-network node as a `(longitude, latitude)` tuple.

Nodes are compared **exactly**. That is safe rather than fragile: ELVIS is a topological network
whose segments genuinely share their endpoints, and the service reprojects a shared node to the
same output coordinate every time, so two segments meeting at a node agree bit for bit. Measured on
the Oslofjord network, exact matching finds the same 481 mouths in EPSG:4326 as in the native
EPSG:25833.
"""
nve_node(point) = (Float64(point[1]), Float64(point[2]))

"""
    nve_cache_covers(filepath, domain)

Whether a cached map-service response is present and was fetched for exactly this domain.
"""
function nve_cache_covers(filepath, domain)
    isfile(filepath) || return false

    cached = get(JSON.parsefile(filepath), "domain", nothing)
    (cached isa AbstractVector && length(cached) == 4) || return false

    return all(isapprox(Float64(cached[k]), domain[k]; atol = 1e-9) for k = 1:4)
end

"""
    nve_write_cache(filepath, domain, key, value)
    nve_read_cache(filepath, key, description)

Write and read one reduced map-service cache: the domain it was fetched for, plus `value` under
`key`.
"""
function nve_write_cache(filepath, domain, key, value)
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        write(io, JSON.json(Dict("domain" => collect(domain), key => value)))
    end
    return filepath
end

function nve_read_cache(filepath, key, description)
    isfile(filepath) || error(
        "Cached $description $filepath does not exist. Run " *
        "`julia --project -m FjordSim add_rivers` for this setup, which downloads it.",
    )

    cache = JSON.parsefile(filepath)

    return get(cache, key, []), Tuple(Float64.(get(cache, "domain", zeros(4))))
end

nve_read_catchments(config::NVERiversConfig) =
    nve_read_cache(nve_catchment_path(config), "catchments", "REGINE catchment list")

nve_read_network(config::NVERiversConfig) =
    nve_read_cache(nve_network_path(config), "segments", "ELVIS river network")

"""
    nve_download_network(config, target_grid)

Fetch and cache the two map-service layers the outlet discovery reads: REGINE's sea-draining
catchments over the domain, and the ELVIS river network over the domain grown by `config.margin`.

Both are skipped when a cache for the same domain is already there, so this is cheap to re-run;
changing the grid changes the domain and re-downloads both.
"""
function nve_download_network(config::NVERiversConfig, target_grid)
    domain = nve_domain(target_grid)
    query = nve_grown_domain(domain, config.margin)

    catchment_file = nve_catchment_path(config)
    if !nve_cache_covers(catchment_file, domain)
        catchments = Dict{String,Any}()

        for feature in nve_map_features(
            config, NVE_CATCHMENT_LAYER, NVE_CATCHMENT_FIELDS, query, false,
            "REGINE catchments draining to the sea",
        )
            properties = get(feature, "properties", nothing)
            isnothing(properties) && continue
            number = get(properties, "vassdragsnr", nothing)
            number isa AbstractString || continue

            name = get(properties, "nedborfeltnavn", nothing)
            catchments[number] = Dict{String,Any}(
                "name" => name isa AbstractString ? name : number,
                "area_km2" => get(properties, "areal_km2", nothing),
                "qnormal_mm3aar" => get(properties, "qnormal_mm3aar", nothing),
            )
        end

        nve_write_cache(catchment_file, domain, "catchments", catchments)
        @info "Cached $(length(catchments)) sea-draining REGINE catchments in $catchment_file"
    end

    network_file = nve_network_path(config)
    if !nve_cache_covers(network_file, domain)
        segments = Any[]

        for feature in nve_map_features(
            config, NVE_ELVENETT_LAYER, NVE_ELVENETT_FIELDS, query, true,
            "the ELVIS river network",
        )
            endpoints = nve_segment_endpoints(feature)
            isnothing(endpoints) && continue
            properties = get(feature, "properties", Dict{String,Any}())

            number = get(properties, "vassdragsnr", "")
            catchment = get(properties, "vnrnfelt", "")
            strahler = get(properties, "elveordenstrahler", 0)

            push!(
                segments,
                Dict{String,Any}(
                    "start" => collect(endpoints[1]),
                    "stop" => collect(endpoints[2]),
                    "vassdragsnr" => number isa AbstractString ? number : "",
                    "vnrnfelt" => catchment isa AbstractString ? catchment : "",
                    "strahler" => strahler isa Number ? Int(strahler) : 0,
                ),
            )
        end

        nve_write_cache(network_file, domain, "segments", segments)
        @info "Cached $(length(segments)) ELVIS river segments in $network_file"
    end

    return nve_download_directory(config)
end

"""
    nve_catchment_discharge(catchment)

A catchment's normal annual runoff as a mean discharge in m³/s, or `nothing` when it has none.

`qnormal_mm3aar` is the 1961–90 normal in millions of m³ per year. It is the size signal that makes
a river with no gauge still scale by how large it is, and it checks out against observation where
both exist: NVE gives Glomma 775.0 m³/s where Solbergfoss — 97 % of the catchment — measured a 2020
mean of about 846 in a wetter-than-normal year.
"""
function nve_catchment_discharge(catchment)
    runoff = get(catchment, "qnormal_mm3aar", nothing)
    (runoff isa Number && runoff > 0) || return nothing
    return Float64(runoff) * 1e6 / NVE_SECONDS_PER_YEAR
end

"""
    nve_terminal_rank(segment)

How strongly one terminal segment claims to be its catchment's river mouth: Strahler order first,
`vassdragsnr` as a deterministic tie-break. See step 5 of `nve_outlets`.
"""
nve_terminal_rank(segment) = (get(segment, "strahler", 0), get(segment, "vassdragsnr", ""))

"""
    nve_discover_outlets(segments, catchments, domain, config)

The river mouths in `domain` at or above `config.minimum_discharge`, as `NVERiver`s with ids
assigned by descending discharge. Steps 1 to 6 of `nve_outlets`, which documents why each is there.
"""
function nve_discover_outlets(segments, catchments, domain, config::NVERiversConfig)
    starts = Set(nve_node(segment["start"]) for segment in segments)
    west, south, east, north = domain

    terminals = Dict{String,Vector{Any}}()

    for segment in segments
        node = nve_node(segment["stop"])
        node in starts && continue

        catchment_number = get(segment, "vnrnfelt", "")
        haskey(catchments, catchment_number) || continue
        (west <= node[1] <= east && south <= node[2] <= north) || continue

        push!(get!(terminals, catchment_number, Any[]), segment)
    end

    found = Tuple{Float64,NVERiver}[]

    for (catchment_number, candidates) in terminals
        catchment = catchments[catchment_number]
        discharge = nve_catchment_discharge(catchment)
        isnothing(discharge) && continue
        discharge < config.minimum_discharge && continue

        sort!(candidates; by = nve_terminal_rank, rev = true)
        length(candidates) > 1 && @warn(
            "Catchment $catchment_number ($(catchment["name"])) reaches the sea at " *
            "$(length(candidates)) places inside the domain; using the one at " *
            "$(nve_node(first(candidates)["stop"])). If this river genuinely has more than one " *
            "mouth, give each a `vassdragsnr` override with its own `discharge_fraction`."
        )

        segment = first(candidates)
        node = nve_node(segment["stop"])

        push!(
            found,
            (
                discharge,
                NVERiver(
                    vassdragsnr = segment["vassdragsnr"],
                    name = catchment["name"],
                    longitude = node[1],
                    latitude = node[2],
                    catchment_discharge = discharge,
                    plume_depth = config.default_plume_depth,
                ),
            ),
        )
    end

    # Descending discharge, `vassdragsnr` ascending behind it. The tie-break is load-bearing rather
    # than tidy: `terminals` is a `Dict`, so its iteration order is not something to rely on, and
    # equal runoffs are not rare — Borreelva and Selvikelva are both 0.67 m³/s on this domain. With
    # discharge alone the two would swap ids between runs, and an id is what a log line, a plot
    # legend and `nve_outlet` all identify a river by.
    sort!(found; by = entry -> (-first(entry), last(entry).vassdragsnr))

    return [nve_river_with_id(river, id) for (id, (_, river)) in enumerate(found)]
end

"""
    nve_river_with_id(river, id)

`river` with a different id, everything else untouched.
"""
nve_river_with_id(river::NVERiver, id) = NVERiver(
    river.vassdragsnr,
    Int(id),
    river.name,
    river.longitude,
    river.latitude,
    river.temperature_station,
    river.discharge_station,
    river.mean_discharge,
    river.discharge_fraction,
    river.plume_depth,
    river.catchment_discharge,
)

"""
    nve_apply_overrides(outlets, config)

Step 7 of `nve_outlets`: fold `config.outlets` into the discovered mouths, matching on
`vassdragsnr`.

Three things are errors rather than warnings, because each would otherwise be a plausible-looking
wrong answer that no later stage can notice: an override naming no `vassdragsnr` (which is a
manually stated outlet, and this config does not use those), an override matching no discovered
mouth (a typo, or a river below the threshold), and two overrides naming the same mouth. So is an
id collision after merging, since two rivers sharing an id would make `nve_outlet` return the wrong
one for every hook that looks a `RiverLocation` back up.
"""
function nve_apply_overrides(outlets, config::NVERiversConfig)
    numbers = [outlet.vassdragsnr for outlet in outlets]
    seen = String[]

    for override in config.outlets
        isempty(override.vassdragsnr) && error(
            "River \"$(override.name)\" states no `vassdragsnr`, so it is a manually stated " *
            "outlet — but this config discovers its outlets (`minimum_discharge` is " *
            "$(config.minimum_discharge) m³/s). Give it the `vassdragsnr` of the mouth it means " *
            "to override, or drop `minimum_discharge` to state every outlet by hand.",
        )

        override.vassdragsnr in seen && error(
            "Two overrides both name `vassdragsnr` \"$(override.vassdragsnr)\". One mouth takes " *
            "one override.",
        )
        push!(seen, override.vassdragsnr)

        override.vassdragsnr in numbers || error(
            "No discovered river mouth carries `vassdragsnr` \"$(override.vassdragsnr)\" " *
            "(override \"$(override.name)\"). Either the number is wrong, or that river is " *
            "below the config's `minimum_discharge` of $(config.minimum_discharge) m³/s. The " *
            "mouths that were discovered, largest first, are: " *
            join(("$(o.vassdragsnr) ($(o.name))" for o in outlets), ", "),
        )
    end

    merged = [nve_merge_override(outlet, config) for outlet in outlets]

    ids = [river.id for river in merged]
    allunique(ids) || error(
        "Two rivers share an id after applying overrides: $(sort(ids)). Ids are the discharge " *
        "rank unless an override states one, so an override naming an `id` can collide with a " *
        "rank — leave `id` unset to keep the rank.",
    )

    return merged
end

"""
    nve_merge_override(outlet, config)

One discovered mouth with its override folded in, or the mouth unchanged when nothing overrides it.

An unset field keeps what discovery derived: `id = 0` keeps the discharge rank, `name = ""` keeps
NVE's catchment name, a `NaN` coordinate keeps NVE's mouth, and `plume_depth = 0.0` keeps
`config.default_plume_depth`. That last one costs an override the ability to ask for the surface
level alone on a config whose default is deeper, which is the price of `0.0` being a real depth
rather than a sentinel; nothing needs it yet.

`catchment_discharge` is never overridden — it is what NVE says the catchment is, and an override
that wants a different size states `mean_discharge`, which outranks it.
"""
function nve_merge_override(outlet::NVERiver, config::NVERiversConfig)
    index = findfirst(river -> river.vassdragsnr == outlet.vassdragsnr, config.outlets)
    isnothing(index) && return outlet

    override = config.outlets[index]

    return NVERiver(
        outlet.vassdragsnr,
        iszero(override.id) ? outlet.id : override.id,
        isempty(override.name) ? outlet.name : override.name,
        isfinite(override.longitude) ? override.longitude : outlet.longitude,
        isfinite(override.latitude) ? override.latitude : outlet.latitude,
        override.temperature_station,
        override.discharge_station,
        override.mean_discharge,
        override.discharge_fraction,
        iszero(override.plume_depth) ? outlet.plume_depth : override.plume_depth,
        outlet.catchment_discharge,
    )
end

"""
    nve_download_directory(config)
    nve_observations_path(config, station, parameter, year)
    nve_station_path(config, station)

Where cached HydAPI responses live. One file per (station, parameter, year) for observations, so
each station's download is independently cacheable and re-running `download_rivers` is a no-op.

That granularity costs nothing. A year of daily data for one station and one parameter is 366
records in a single request, which is far inside both documented limits — 10 series and 150 000
observations per request — so there is no batching to be gained, only the ability to skip work
already done.
"""
nve_download_directory(config::NVERiversConfig) =
    joinpath(config.data_root, config.download_directory)

nve_observations_path(config::NVERiversConfig, station, parameter, year) = joinpath(
    nve_download_directory(config),
    "observations_$(station)_$(parameter)_$(year).json",
)

nve_station_path(config::NVERiversConfig, station) =
    joinpath(nve_download_directory(config), "station_$(station).json")

"""
    nve_api_key(config)

The HydAPI key: `config.api_key` if it states one, else the `$NVE_API_KEY_VARIABLE` environment
variable. Raises naming both and where to get a key, rather than letting the request come back as
an HTTP 401 with an empty body.
"""
function nve_api_key(config::NVERiversConfig)
    isempty(config.api_key) || return config.api_key

    key = get(ENV, NVE_API_KEY_VARIABLE, "")
    isempty(key) && error(
        "No NVE HydAPI key. Set the $NVE_API_KEY_VARIABLE environment variable, or give the " *
        "river config an `api_key`. Keys are free at https://hydapi.nve.no/Users — an email " *
        "address and nothing else.",
    )

    return key
end

"""
    nve_observations_url(config, station, parameter, year)

The `/Observations` request for one station, parameter and calendar year.

`ReferenceTime` is an ISO-8601 interval and is **inclusive at both ends**, so asking for
`year-01-01/year+1-01-01` returns the following 1 January as well; `river_series` matches by
calendar date and simply never asks for it. Omitting `ReferenceTime` altogether returns only the
single latest observation, which is why it is never optional here.
"""
nve_observations_url(config::NVERiversConfig, station, parameter, year) = string(
    config.base_url,
    "/Observations?StationId=", station,
    "&Parameter=", parameter,
    "&ResolutionTime=", config.resolution_time,
    "&ReferenceTime=", year, "-01-01/", year + 1, "-01-01",
)

"""
    nve_request(config, url, filepath, description)

Fetch `url` into `filepath` with the API key header, skipping the request when the file is already
there. Returns `filepath`, or `nothing` when the server answered 404.

HydAPI's status codes each mean something specific and none of them is a generic failure, so they
are translated rather than passed through as a `RequestError` naming only a number:

- **401** is a missing or rejected key, and its body is empty, so there is nothing to report but
  the key itself.
- **404** is *the series does not exist* — a station HydAPI has never heard of, or one that carries
  no such parameter at all. The caller decides whether that is fatal.
- **400** is a validation failure carrying a parseable `errors` map, worth surfacing verbatim.
- **429** is the 5-requests-per-second limit.

**A window with no records is not a 404.** A series that exists but holds nothing in the requested
interval comes back as **HTTP 200 with `observationCount: 0`** and an empty `observations` array —
Sarpsfoss (`2.31.0`) has water temperature for 2002–2009 and answers exactly that for 2020. So a
successful request is not the same as data, which is why `download_rivers` checks the count rather
than trusting the status; see `nve_observation_count`.

`validate_river_download` then catches the other way this can fail quietly: a response that is an
HTML page rather than the JSON body, which starts `<` instead of `{`.
"""
function nve_request(config::NVERiversConfig, url, filepath, description)
    isfile(filepath) && return filepath

    mkpath(dirname(filepath))
    @info "Downloading $description from HydAPI"
    headers = ["X-API-Key" => nve_api_key(config), "Accept" => "application/json"]

    try
        Downloads.download(url, filepath; headers)
    catch exception
        rm(filepath; force = true)
        exception isa Downloads.RequestError || rethrow()
        status = exception.response.status

        status == 404 && return nothing
        status == 401 && error(
            "HydAPI rejected the API key (HTTP 401, empty body) while fetching $description. " *
            "Check the $NVE_API_KEY_VARIABLE environment variable, or the config's `api_key`.",
        )
        status == 429 && error(
            "HydAPI rate limit reached (HTTP 429) while fetching $description. The limit is 5 " *
            "requests per second; NVE_REQUEST_INTERVAL is meant to keep a download inside it.",
        )
        status == 400 && error(
            "HydAPI rejected the request for $description (HTTP 400): $(exception.message). " *
            "Request was $url",
        )
        # `status` is 0 for a failure below HTTP — DNS, TLS, a dropped connection — so the curl
        # message is the only thing that says what went wrong there.
        error(
            "HydAPI request for $description failed" *
            (iszero(status) ? "" : " with HTTP $status") *
            ": $(exception.message). Request was $url",
        )
    end

    sleep(NVE_REQUEST_INTERVAL)

    return validate_river_download(filepath, url)
end

"""
    download_rivers(target_grid, config::NVERiversConfig)

Fetch everything this config's rivers need into `nve_download_directory(config)`, in two halves.

First the **map services**, which need no key: REGINE's sea-draining catchments and the ELVIS river
network over `target_grid`'s domain, from which `nve_outlets` derives the river mouths. This is
what `target_grid` is for, and it is why this hook takes one — the same argument
`download_forcing` and `download_boundaries` already make, that a dataset needs the domain bounds
and `x_domain`/`y_domain` provide them for any grid. A config naming no `minimum_discharge` states
its own outlets, so the download is skipped and the grid goes unused.

Then **HydAPI**: every series the outlets name — water temperature and discharge, one calendar year
at a time — plus each station's metadata. Files already present are left alone, so this is cheap to
re-run and a partial download resumes.

A configured station with no data for a requested year is an **error**, not a `NaN` row: the setup
named that station and that year, so an empty answer is a statement about the setup. A river that
genuinely has no temperature says so with `temperature_station = ""` instead. Station *metadata*
is the exception — it is only ever a sanity log and a fallback size, so a station missing from
`/Stations` (which hides discontinued stations) is a warning.

"No data" has **two** forms and both are checked here, because only one of them fails the request:
a series HydAPI does not have at all is a 404, while a series with nothing in the requested window
is a perfectly successful 200 carrying zero observations. The empty cache file is removed so a
later run re-requests rather than reading the emptiness back as a gap-filled row of `NaN`.
"""
function download_rivers(target_grid, config::NVERiversConfig)
    isempty(config.years) && error("The NVE river config names no years.")
    !isfinite(config.minimum_discharge) && isempty(config.outlets) && error(
        "The NVE river config names neither outlets nor a `minimum_discharge`, so it has no " *
        "rivers at all. State the outlets, or give a discharge threshold to discover them.",
    )

    isfinite(config.minimum_discharge) && nve_download_network(config, target_grid)

    outlets = nve_outlets(config)
    isfinite(config.minimum_discharge) &&
        @info "Discovered $(length(outlets)) NVE river mouths at or above " *
              "$(config.minimum_discharge) m³/s: " *
              join(("$(o.name) [$(o.vassdragsnr)] " *
                    "$(round(o.catchment_discharge, digits = 2)) m³/s" for o in outlets), ", ")

    for outlet in outlets, (parameter, station) in nve_outlet_series(outlet)
        for year in config.years
            filepath = nve_observations_path(config, station, parameter, year)
            description = "$(outlet.name) parameter $parameter at station $station for $year"
            result = nve_request(
                config, nve_observations_url(config, station, parameter, year),
                filepath, description,
            )

            if !isnothing(result) && iszero(nve_observation_count(filepath))
                rm(filepath; force = true)
                result = nothing
            end

            isnothing(result) && error(
                "HydAPI has no parameter $parameter at resolution $(config.resolution_time) for " *
                "station $station in $year (requested for $(outlet.name)). Check the station's " *
                "coverage with $(config.base_url)/Series?StationId=$station, or leave the " *
                "station name empty if this river has no such series.",
            )
        end
    end

    for station in nve_stations(config)
        nve_download_station(config, station)
    end

    return nve_download_directory(config)
end

"""
    nve_outlet_series(outlet)

The `(parameter, station)` pairs `outlet` needs downloading, skipping the ones it leaves unnamed.
"""
function nve_outlet_series(outlet::NVERiver)
    series = Tuple{Int,String}[]
    isempty(outlet.temperature_station) ||
        push!(series, (NVE_WATER_TEMPERATURE, outlet.temperature_station))
    isempty(outlet.discharge_station) ||
        push!(series, (NVE_DISCHARGE, outlet.discharge_station))
    return series
end

"""
    nve_stations(config)

Every distinct station id the config names, in a stable order.
"""
function nve_stations(config::NVERiversConfig)
    stations = String[]
    for outlet in nve_outlets(config), (_, station) in nve_outlet_series(outlet)
        station in stations || push!(stations, station)
    end
    return stations
end

"""
    nve_download_station(config, station)

Cache one station's metadata, and log what HydAPI says the gauge is.

`/Stations` needs `Active=1` — an inverted enum where `0` or omitted means active-only — and even
then it does not carry stations that have been decommissioned, while `/Series` does. So a 404 falls
back to `/Series?StationId=…`, and a station absent from both is a warning: metadata is a sanity
log and a fallback river size, never something the pipeline needs.
"""
function nve_download_station(config::NVERiversConfig, station)
    filepath = nve_station_path(config, station)
    description = "metadata for station $station"

    result = nve_request(
        config, "$(config.base_url)/Stations?StationId=$station&Active=1", filepath, description,
    )
    isnothing(result) && (result = nve_request(
        config, "$(config.base_url)/Series?StationId=$station", filepath, description,
    ))

    if isnothing(result)
        @warn "HydAPI has no metadata for station $station; its catchment size is unavailable"
        return nothing
    end

    entry = nve_first_entry(filepath)
    isnothing(entry) && return filepath

    name = get(entry, "stationName", "?")
    river = get(entry, "riverName", "")
    masl = get(entry, "masl", nothing)
    area = get(entry, "drainageBasinArea", nothing)
    @info "Station $station: $name" *
          (isempty(river) ? "" : " on $river") *
          (isnothing(masl) ? "" : ", $masl m a.s.l.") *
          (isnothing(area) ? "" : ", catchment $(round(area, digits = 1)) km²")

    return filepath
end

"""
    nve_first_entry(filepath)

The first element of a cached response's `data` array, or `nothing` when the file holds no data.
Every HydAPI response is the same envelope — `currentLink`, `apiVersion`, `license`, `createdAt`,
`queryTime`, `itemCount`, `data` — so one accessor serves observations, stations and series alike.
"""
function nve_first_entry(filepath)
    isfile(filepath) || return nothing
    data = get(JSON.parsefile(filepath), "data", nothing)
    (isnothing(data) || isempty(data)) && return nothing
    return first(data)
end

"""
    nve_observation_count(filepath)

How many observations a cached observations response actually carries.

This is the check a status code cannot make. A series HydAPI does not have answers 404, but a series
that exists and simply holds nothing in the requested window answers **200** with
`observationCount: 0` and an empty `observations` array — so a request that succeeded is not
evidence of data, and without this a configured station-year outside its record would be cached as
an empty file and read back as a row of `NaN` that nothing ever complained about.
"""
function nve_observation_count(filepath)
    entry = nve_first_entry(filepath)
    isnothing(entry) && return 0
    return length(get(entry, "observations", ()))
end

"""
    nve_daily_values(filepath)

`Date => value` for one cached observations response.

`value` may be `null` on an otherwise gap-free axis, and such a record is dropped here rather than
carried as a missing marker — `river_series` fills the gap from the nearest date it does have.
Daily records are timestamped **11:00 UTC**, not midnight (HydAPI's convention for a daily mean),
which is exactly why the key is a `Date` and matching is by calendar date.
"""
function nve_daily_values(filepath)
    records = Dict{Date,Float64}()
    entry = nve_first_entry(filepath)
    isnothing(entry) && return records

    for observation in get(entry, "observations", ())
        value = get(observation, "value", nothing)
        value isa Number || continue
        time = get(observation, "time", nothing)
        time isa AbstractString || continue
        records[Date(DateTime(first(time, 19)))] = Float64(value)
    end

    return records
end

"""
    nve_station_series(config, station, parameter)

Every cached daily record of one station and parameter, merged across the config's `years`.
"""
function nve_station_series(config::NVERiversConfig, station, parameter)
    records = Dict{Date,Float64}()
    for year in config.years
        filepath = nve_observations_path(config, station, parameter, year)
        isfile(filepath) || error(
            "Cached HydAPI response $filepath does not exist. Run `download_rivers` first.",
        )
        merge!(records, nve_daily_values(filepath))
    end
    return records
end

"""
    river_locations(config::NVERiversConfig)

The outlets, as `RiverLocation`s, from `nve_outlets`.

The coordinates are river mouths — either discovered from the ELVIS river network or stated by the
setup, never a gauging station's, for the reasons `NVERiver` gives. When they are discovered this
reads the cached map-service responses, so `download_rivers` has to have run; `add_rivers` calls
that first, and the error names the step when it has not.
"""
river_locations(config::NVERiversConfig) = [
    RiverLocation(outlet.id, outlet.name, outlet.longitude, outlet.latitude) for
    outlet in nve_outlets(config)
]

"""
    river_series(config::NVERiversConfig, times)

River values on `times`: the temperature series read from the cached HydAPI responses, plus the
config's `constants` — `S = 0` by default, since HydAPI has no salinity parameter and a river is
fresh.

Rows are in `river_locations` order, matched to `times` by **calendar date**, because HydAPI stamps
a daily mean at 11:00 UTC while a prepared forcing file's axis sits wherever its own source put it.

A river with no temperature station, and a date no response covers, both come out as `NaN`, which
is the file's `_FillValue` and reads back as the `-999.0` sentinel every `ForcingFromFile` branch
gates on — so those cells are simply not temperature-forced. `nve_fill_gaps!` first fills the
short gaps that a single missing observation leaves, since a river blinking out of its temperature
forcing for one day of a year is noise rather than information.
"""
function river_series(config::NVERiversConfig, times)
    outlets = nve_outlets(config)
    temperatures = fill(NaN32, length(outlets), length(times))
    dates = Date.(times)

    for (row, outlet) in enumerate(outlets)
        isempty(outlet.temperature_station) && continue
        records = nve_station_series(config, outlet.temperature_station, NVE_WATER_TEMPERATURE)

        for (column, date) in enumerate(dates)
            temperatures[row, column] = Float32(get(records, date, NaN32))
        end

        nve_fill_gaps!(view(temperatures, row, :), "$(outlet.name) temperature")
    end

    series = Dict{String,Matrix{Float32}}(config.temperature_name => temperatures)
    for (name, value) in config.constants
        series[name] = fill(Float32(value), length(outlets), length(times))
    end

    return series
end

"""
    nve_fill_gaps!(row, label)

Replace each non-finite entry of `row` with the nearest finite one in time, in place. A row with
nothing finite in it at all is left alone and stays inert.

Nearest-in-time is the choice `SourceFill` and `fill_boundary_gaps!` already make elsewhere in this
module, for the same reason: the alternative to a neighbouring value is no forcing at that step,
which is a larger and less physical discontinuity than a repeated one.
"""
function nve_fill_gaps!(row, label)
    finite = findall(isfinite, row)
    isempty(finite) && return row

    gap_count = length(row) - length(finite)
    iszero(gap_count) && return row

    for index in eachindex(row)
        isfinite(row[index]) && continue
        nearest = finite[argmin(abs.(finite .- index))]
        row[index] = row[nearest]
    end

    @warn "Filled $gap_count of $(length(row)) $label records from the nearest date"

    return row
end

"""
    river_lambdas(config::NVERiversConfig, cells, target_grid)

The relaxation coefficient at each river cell, `λ = Q̄ / V`, capped at
`1 / config.minimum_relaxation_timescale`.

Relaxing a plume of volume `V` towards `S = 0` at rate `λ` removes salt at `λ V S`, while a
discharge `Q` dilutes it at `Q S` — so `Q̄ / V` is the rate a river of that size actually implies,
and it is the reason this hook exists. One shared timescale cannot express it: across Oslofjord's
rivers `Q̄` spans Akerselva's ≈3 m³/s to Glomma's ≈846 m³/s, a factor of 300, which becomes a
factor of 50 in `λ` once each river's plume volume is accounted for.

`V` is summed over the cell's own `levels`, so a river given a deeper plume by `river_plume_depth`
is nudged proportionally more gently — the freshwater is spread through more water, which is what a
plume depth means.

`Q̄` resolves in a stated order, and only the first that is available is used, with the outlet's
`discharge_fraction` applied to whichever answered:

1. the outlet's `mean_discharge`, if it states one;
2. the mean of the downloaded daily discharge series;
3. the **REGINE catchment normal** the outlet was discovered with — the 1961–90 normal annual
   runoff of the catchment that actually drains through this mouth — which is what makes a river
   with no gauge at all still scale by its size, and is why most of a discovered river list needs
   no station;
4. `annualRunoff` from the cached station metadata, the same quantity for the *gauge's* catchment
   rather than the outlet's, which is all a manually stated outlet can reach;
5. nothing, in which case `1 / config.relaxation_timescale` is used and the outlet is named in a
   warning.

# The cap is not optional
`ForcingFromFile` reads `λ > 1` as an x-flux, and the relaxation term is explicit, so `λ Δt` must
stay well below 1. Neither bound is comfortable at the top of the range: Drammenselva's peak
1012 m³/s into a single surface cell of Drammensfjord's 99 m grid gives `λ = 0.102 s⁻¹`, i.e.
`λ Δt = 1.02` at the 10 s steps that run actually takes — an unstable relaxation, from the
physically correct coefficient. `minimum_relaxation_timescale` is what keeps `λ Δt` bounded
regardless of how large a river is or how small its cell.
"""
function river_lambdas(config::NVERiversConfig, cells, target_grid)
    maximum_lambda = 1 / config.minimum_relaxation_timescale
    fallback = 1 / config.relaxation_timescale

    # A timescale of a second or less puts λ outside the relaxation regime entirely, where
    # `ForcingFromFile` reads it as an x-flux — a silent change of forcing term, not a strong nudge.
    # Checked here rather than left to the cap because the cap is what would otherwise enforce it.
    for (name, timescale) in
        (:minimum_relaxation_timescale => config.minimum_relaxation_timescale,
         :relaxation_timescale => config.relaxation_timescale)
        timescale > 1 || error(
            "The river config's `$name` is $timescale s, which makes λ = $(1 / timescale) ≥ 1. " *
            "`ForcingFromFile` reads λ > 1 as an x-flux, so a relaxation coefficient has to stay " *
            "inside (-1, 1) — give it a timescale of more than one second.",
        )
    end

    lambdas = Vector{Float32}(undef, length(cells))

    for (n, cell) in enumerate(cells)
        outlet = nve_outlet(config, cell.location.id)
        discharge, source = nve_mean_discharge(config, outlet)

        if isnothing(discharge)
            @warn "No discharge for river $(outlet.id) ($(outlet.name)): using the config's " *
                  "relaxation_timescale of $(config.relaxation_timescale) s. Give the outlet a " *
                  "`discharge_station` or a `mean_discharge` to scale it by river size."
            lambdas[n] = Float32(fallback)
            continue
        end

        plume_volume = sum(
            volume(cell.i, cell.j, level, target_grid, Center(), Center(), Center()) for
            level in cell.levels
        )
        lambda = min(discharge / plume_volume, maximum_lambda)
        capped = lambda == maximum_lambda ? " (capped)" : ""
        @info "River $(outlet.id) ($(outlet.name)): Q̄ = $(round(discharge, digits = 2)) m³/s " *
              "from $source, plume $(round(plume_volume, digits = 0)) m³, " *
              "λ = $lambda s⁻¹ (τ = $(round(Int, 1 / lambda)) s)$capped"
        lambdas[n] = Float32(lambda)
    end

    return lambdas
end

"""
    nve_mean_discharge(config, outlet)

`(discharge, source)` in m³/s for one outlet, or `(nothing, source)` when its size is unknown.
`source` names which rule supplied it, so a run's log states where the number came from. The
precedence is documented on `river_lambdas`.

`discharge_fraction` is applied last, to whatever that precedence resolved, so a river reaching the
sea through two mouths can share one gauge between them.
"""
function nve_mean_discharge(config::NVERiversConfig, outlet::NVERiver)
    discharge, source = nve_resolve_discharge(config, outlet)

    isnothing(discharge) && return nothing, source
    isone(outlet.discharge_fraction) && return discharge, source

    return discharge * outlet.discharge_fraction,
    "$source, scaled by $(round(outlet.discharge_fraction, digits = 3))"
end

"""
    nve_resolve_discharge(config, outlet)

`nve_mean_discharge` before `discharge_fraction` is applied: the four rules in order, and the name
of the one that answered.
"""
function nve_resolve_discharge(config::NVERiversConfig, outlet::NVERiver)
    isfinite(outlet.mean_discharge) && return outlet.mean_discharge, "the config"

    if !isempty(outlet.discharge_station)
        records = nve_station_series(config, outlet.discharge_station, NVE_DISCHARGE)
        isempty(records) || return mean(values(records)), "station $(outlet.discharge_station)"
    end

    isfinite(outlet.catchment_discharge) &&
        return outlet.catchment_discharge, "the REGINE catchment normal"

    for station in (outlet.discharge_station, outlet.temperature_station)
        isempty(station) && continue
        entry = nve_first_entry(nve_station_path(config, station))
        isnothing(entry) && continue
        runoff = get(entry, "annualRunoff", nothing)
        runoff isa Number && runoff > 0 &&
            return runoff * 1e6 / NVE_SECONDS_PER_YEAR, "the catchment normal at station $station"
    end

    return nothing, "nowhere"
end
