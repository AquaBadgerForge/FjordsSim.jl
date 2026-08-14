# Open-boundary data: the generic half of the boundary pipeline. A dataset adapter subtypes
# `AbstractBoundaryDataConfig` and implements `boundary_time_steps`, `boundary_source_grid` and
# `boundary_variable_names` (plus `download_boundaries` if it fetches data); everything below is
# shared. `src/Forcing/norkyst_boundaries.jl` is the template to copy.
#
# What this exists for: an open lateral boundary needs the exterior state *at* the boundary —
# elevation and barotropic transport for the Flather condition, velocity and tracers for the
# Orlanski one — at a cadence that resolves the tide. The interior forcing is daily means, which
# have the tide averaged out of them. Only the boundary row needs hourly data, so it is a separate
# file of a few hundred MB rather than an hourly interior forcing of hundreds of GB.
#
# Every variable in the file is named for its side (`south_T`, `south_eta`, ...), so a second open
# edge is a set of new variables in the same file rather than a new file or a format change. The
# edge is a `Val` everywhere it decides anything.

const BOUNDARY_DEFLATE_LEVEL = 5
# Longest run of missing hours interpolated without a separate warning.
const BOUNDARY_MAX_GAP = 6

"""
Prepared boundary variables with no depth axis: the free-surface elevation and the two barotropic
velocity components.

`eta` is the exterior elevation the Flather condition (`GravityWaveRadiation`) compares the model's
own boundary `η` against; `ubar`/`vbar` become the exterior barotropic *transport* once multiplied
by the model's column depth. The rest of the prepared variables are full-depth.
"""
const SURFACE_BOUNDARY_VARIABLES = ("eta", "ubar", "vbar")

is_surface_boundary_variable(name) = String(name) in SURFACE_BOUNDARY_VARIABLES

"""
    boundary_full_location(name)

The model location a boundary variable has *before* the edge dimension is reduced away — the
staggering that decides which node row it is sampled at and which NetCDF dimension it is written on.

Separate from `data_location`, which serves the interior forcing file: that one has no surface
variables and no notion of `ubar`/`vbar`, and widening it would change what `forcing_from_file`
reads.
"""
function boundary_full_location(name)
    variable = Symbol(name)
    variable === :u && return (Face, Center, Center)
    variable === :v && return (Center, Face, Center)
    variable === :ubar && return (Face, Center, Nothing)
    variable === :vbar && return (Center, Face, Nothing)
    variable === :eta && return (Center, Center, Nothing)
    return (Center, Center, Center)
end

"""
    boundary_variable_name(edge, name)

The prepared variable's name in the boundary file: its side, an underscore, and the FjordSim
variable name — `boundary_variable_name(:south, "T") == "south_T"`.

The side is in the name rather than in the filename so one file holds however many boundaries a
setup opens.
"""
boundary_variable_name(edge::Symbol, name) = string(validate_open_edge(edge), "_", name)

# --- Edge dispatch ---
#
# One method per edge and a catch-all that raises, rather than a four-branch `if`. `first`/`last`
# select the boundary node along the cross-edge axis: south and west are the low end of their axis,
# north and east the high end, for a Center row and a Face row alike (the Face node list is one
# longer, and its boundary node is still its first or last).

"""
    boundary_node_selector(::Val{edge})

`first` for the edge at the low end of its axis, `last` for the high end. Applied to a node vector
or to a mask axis to pick the boundary row.
"""
boundary_node_selector(::Val{:south}) = first
boundary_node_selector(::Val{:north}) = last
boundary_node_selector(::Val{:west}) = first
boundary_node_selector(::Val{:east}) = last
boundary_node_selector(::Val{edge}) where {edge} =
    throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))

"""
    boundary_along_axis(::Val{edge})

Which horizontal axis the boundary *runs along*: `1` (x) for a south or north edge, `2` (y) for a
west or east one. The other horizontal axis is the one reduced away.
"""
boundary_along_axis(::Val{:south}) = 1
boundary_along_axis(::Val{:north}) = 1
boundary_along_axis(::Val{:west}) = 2
boundary_along_axis(::Val{:east}) = 2
boundary_along_axis(::Val{edge}) where {edge} =
    throw(ArgumentError("open_edge must be one of $LATERAL_EDGES, got :$edge"))

"""
    boundary_location(::Val{edge}, name)

The reduced Oceananigans location of a boundary variable's `FieldTimeSeries`: the cross-edge
horizontal location becomes `Nothing`, and a surface variable's vertical location is `Nothing` too.

These are exactly the reduced flavours `Oceananigans.BoundaryConditions.getbc` accepts as a boundary
condition — `XZFTS` for a south or north edge, `YZFTS` for west or east — which is what lets the
boundary condition pass a series straight through with no discrete-form wrapper.
"""
function boundary_location(::Val{edge}, name) where {edge}
    LX, LY, LZ = boundary_full_location(name)
    return boundary_along_axis(Val(edge)) == 1 ? (LX, Nothing, LZ) : (Nothing, LY, LZ)
end

"""
    boundary_dimension_names(::Val{edge}, name)

NetCDF dimension names for a prepared boundary variable, in Julia (fastest-first) order: the axis
the boundary runs along, then `Nz` unless the variable is a surface one, then `time`.

The along-boundary dimension is staggered with the variable — `south_u` and `south_ubar` are on
`Nx_faces`, everything else on a south edge is on `Nx`. All six spatial dimensions are defined in
every boundary file (`define_forcing_dimensions!`), so adding a second edge writes new variables
into an unchanged dimension set.
"""
function boundary_dimension_names(::Val{edge}, name) where {edge}
    LX, LY, _ = boundary_full_location(name)
    along = if boundary_along_axis(Val(edge)) == 1
        LX === Face ? "Nx_faces" : "Nx"
    else
        LY === Face ? "Ny_faces" : "Ny"
    end

    return is_surface_boundary_variable(name) ? (along, "time") : (along, "Nz", "time")
end

"""
    boundary_mask_slice(::Val{edge}, mask, surface)

The boundary row of a full `water_mask`, as a `(nx, 1, nz)` or `(1, ny, nz)` array — one cell thick
across the edge, so every three-dimensional structure in the forcing core (`PreparedVariable`, the
interpolation kernel's index space) applies unchanged.

`surface` additionally takes the topmost level, giving a `(nx, 1, 1)` mask for `eta`, `ubar` and
`vbar`. The topmost level is the surface because `znodes` increases upwards.
"""
function boundary_mask_slice(::Val{edge}, mask, surface) where {edge}
    select = boundary_node_selector(Val(edge))
    levels = surface ? (size(mask, 3):size(mask, 3)) : (1:size(mask, 3))

    return if boundary_along_axis(Val(edge)) == 1
        row = select(axes(mask, 2))
        mask[:, row:row, levels]
    else
        column = select(axes(mask, 1))
        mask[column:column, :, levels]
    end
end

"""
    boundary_target_nodes(::Val{edge}, target_grid, LX, LY)

Longitudes and latitudes of the boundary row at the location `(LX, LY)`: the full node list along
the edge, and the single cross-edge node the boundary sits at.

A Face-located variable across the edge (`v` on a south edge, `u` on a west one) lands exactly on the
domain boundary, which is the face the open condition writes. A Center-located one is sampled at the
first interior cell centre rather than at the halo centre just outside it — half a cell, ~200 m on
the Oslofjord grid, against an 800 m source: below the resolution of the data being interpolated.
"""
function boundary_target_nodes(::Val{edge}, target_grid, LX, LY) where {edge}
    select = boundary_node_selector(Val(edge))
    longitude = Array(λnodes(target_grid, LX()))
    latitude = Array(φnodes(target_grid, LY()))

    return if boundary_along_axis(Val(edge)) == 1
        longitude, [select(latitude)]
    else
        [select(longitude)], latitude
    end
end

"""
    boundary_domain(edges, target_grid, margin)
    boundary_domain(::Val{edge}, target_grid, margin)

The longitude/latitude box a download needs to cover the boundary rows: for one edge, the full extent
along it and a `margin`-degree band around the edge itself across it; for several, the smallest box
containing each of their bands.

A thin band rather than the whole domain box, because the boundary file is hourly: the full box at
24 records a day is more than twenty times the interior download, while the band is a tenth of the
box. `margin` must still leave room for the source cells surrounding the boundary row, since the
interpolation is bilinear.

**That saving is a single-edge one.** Two opposite edges, or a domain in the open ocean naming all
four, have bands whose bounding box is the whole domain plus the margin — so a multi-edge download is
the full box at hourly cadence, and the interior of it is downloaded and never read. Keeping it one
box is what lets `boundary_time_steps` and `boundary_source_grid` stay single-file hooks and
`prepare_boundaries` interpolate every edge from one source. Per-edge downloads would be the
optimization, at the price of both hooks becoming per-edge.
"""
function boundary_domain(edges, target_grid, margin)
    boxes = [boundary_domain(Val(edge), target_grid, margin) for edge in lateral_edges(edges)]
    isempty(boxes) && throw(ArgumentError("boundary_domain needs at least one open edge, got none."))

    longitude = (minimum(box -> box[1][1], boxes), maximum(box -> box[1][2], boxes))
    latitude = (minimum(box -> box[2][1], boxes), maximum(box -> box[2][2], boxes))

    return longitude, latitude
end

function boundary_domain(::Val{edge}, target_grid, margin) where {edge}
    longitude = x_domain(target_grid)
    latitude = y_domain(target_grid)
    select = boundary_node_selector(Val(edge))

    return if boundary_along_axis(Val(edge)) == 1
        edge_latitude = select(latitude)
        (longitude[1] - margin, longitude[2] + margin), (edge_latitude - margin, edge_latitude + margin)
    else
        edge_longitude = select(longitude)
        (edge_longitude - margin, edge_longitude + margin), (latitude[1] - margin, latitude[2] + margin)
    end
end

# --- Extension hooks ---

"""
    boundary_time_steps(config)

Every source time record available to `config`, as `SourceRecord`s sorted by date with duplicates
dropped. `prepare_boundaries` completes them to a gap-free hourly axis.

A new dataset implements this on its `AbstractBoundaryDataConfig` subtype.
"""
function boundary_time_steps end

"""
    boundary_source_grid(config, filepath)

Geometry of the downloaded boundary data in `filepath`. Return a `ProjectedSourceGrid` for a source
on a regular grid in projected meters, which reuses `source_field_grid` and
`projected_target_nodes` unchanged.
"""
function boundary_source_grid end

"""
    boundary_variable_names(config)

Mapping from source variable name to the FjordSim boundary name it becomes, e.g. `"zeta" => "eta"`.
Only the intersection with `config.parameters` is prepared, so this also declares which variables
the dataset can supply.
"""
function boundary_variable_names end

"""
    boundary_source_slab(config, reader, step, source_name)

One source slab of `source_name` for `step`, as `write_boundaries_file` interpolates from.

Defaults to `blended_slab`, which reads the variable straight out of the downloaded file — the right
answer for any variable that already means what its name says. A source overrides this when a
variable needs deriving from more than one of its own, which a vector component does: the pipeline
prepares one variable at a time, so a rotation from grid-relative to geographic axes has nowhere else
to happen. Rotating after interpolation is not an option, since on a south edge `ubar` lands on
`Nx_faces` and `vbar` on `Nx`, different node counts under different masks.

Both slabs of a pair are available here, from the same `reader` at the same `step`, so a component
derived this way stays consistent with its partner.
"""
boundary_source_slab(config::AbstractBoundaryDataConfig, reader, step, source_name) =
    blended_slab(reader, step, source_name)

"""
    download_boundaries(config::FjordConfig)
    download_boundaries(target_grid, config::AbstractBoundaryDataConfig)

Fetch and subset the hourly source data a later `prepare_boundaries` call reads, into
`boundary_data_directory(config)`.

The `FjordConfig` form is the generic driver, the same shape as `download_forcing`: it builds the
setup's grid on the CPU and dispatches on the setup's boundary config, which is also where the open
edges are stated — so a dataset reads them with `open_edges(config)` rather than being handed them. It
says
nothing about forcing: a setup can name boundary data and no interior forcing, or the other way
round.
"""
function download_boundaries(config::FjordConfig)
    boundaries = config.boundary_config
    isnothing(boundaries) && return nothing

    return download_boundaries(domain_grid(config.grid_config, CPU()), boundaries)
end

download_boundaries(target_grid, ::Nothing) = nothing

"""
    boundary_date_range(config)

First and last date the prepared boundary file covers, as a `(first, last)` tuple, or `nothing` for
a source that cannot report its dates.

The mirror of `forcing_date_range` and `atmosphere_date_range`: `validate_time_coverage` checks the
run's window against it, which is what keeps the reader's `Cyclical()` time indexing from quietly
wrapping a boundary series back to the far end of the year.
"""
boundary_date_range(config::AbstractBoundaryDataConfig) = NCDataset(boundary_data_path(config)) do ds
    dates = ds["time"][:]
    (first(dates), last(dates))
end

boundary_date_range(::Nothing) = nothing

# --- Preparation ---

"""
    hourly_time_steps(records; max_gap = BOUNDARY_MAX_GAP)

`uniform_time_steps` at a one-hour spacing: the cadence of the boundary data, and the reason it is a
separate file from the daily interior forcing.
"""
hourly_time_steps(records; max_gap = BOUNDARY_MAX_GAP) =
    uniform_time_steps(records, Hour(1); max_gap, unit = "hour")

"""
    surface_source_field_grid(source::ProjectedSourceGrid, architecture)

The source subset as a single-level `RectilinearGrid`, for a source variable with no depth axis.

`solve_vertical_faces` cannot build this: for one level it leaves the lower bound at `-Inf` and puts
the first face there. The single cell is centred on `SURFACE_TARGET_DEPTH` instead, which is also the
vertical node `prepared_boundary_variable` gives a surface variable — so the trilinear interpolation
lands exactly on the cell centre and the vertical direction contributes nothing.
"""
function surface_source_field_grid(source::ProjectedSourceGrid, architecture = CPU())
    Δx = source.x[2] - source.x[1]
    Δy = source.y[2] - source.y[1]

    return RectilinearGrid(
        architecture;
        size = (length(source.x), length(source.y), 1),
        x = collect(range(source.x[1] - Δx / 2, step = Δx, length = length(source.x) + 1)),
        y = collect(range(source.y[1] - Δy / 2, step = Δy, length = length(source.y) + 1)),
        z = [SURFACE_TARGET_DEPTH - 1, SURFACE_TARGET_DEPTH + 1],
        topology = (Bounded, Bounded, Bounded),
    )
end

"""
The vertical node a surface boundary variable is interpolated at, and the centre of
`surface_source_field_grid`'s single cell. Any value works as long as the two agree; zero is the
free surface.
"""
const SURFACE_TARGET_DEPTH = 0.0

"""
    prepared_boundary_variable(source_name, edge, target_grid, source, filepath, config)

Build the boundary-row mask, projected target nodes and source fill for one boundary variable.

Reuses `PreparedVariable`, whose three-dimensional mask and lambda fields hold a one-cell-thick
slab perfectly well. `lambda` is a zero array that is never written — the boundary file carries no
relaxation rates, the nudging timescales being a property of the boundary *condition* — and
`dimensions` is the three- or two-name tuple this file uses rather than the forcing file's four.
"""
function prepared_boundary_variable(
    source_name,
    edge,
    target_grid,
    source,
    filepath,
    config::AbstractBoundaryDataConfig,
)
    name = boundary_variable_names(config)[source_name]
    LX, LY, _ = boundary_full_location(name)
    surface = is_surface_boundary_variable(name)
    @info "Preparing boundary target nodes and source fill for $source_name -> $name"

    mask = boundary_mask_slice(Val(edge), water_mask(target_grid, LX, LY, edge), surface)
    longitude, latitude = boundary_target_nodes(Val(edge), target_grid, LX, LY)
    size(mask)[1:2] == (length(longitude), length(latitude)) || error(
        "Boundary node count $(length(longitude))x$(length(latitude)) does not match the mask " *
        "$(size(mask)[1:2]) for $name; the grid topology must be bounded in both directions.",
    )

    x, y = projected_target_nodes(longitude, latitude, source)
    shape = (length(longitude), length(latitude))
    z = surface ? [SURFACE_TARGET_DEPTH] : Array(znodes(target_grid, Center()))

    return PreparedVariable(
        source_name,
        boundary_variable_name(edge, name),
        boundary_dimension_names(Val(edge), name),
        mask,
        reshape(x, shape),
        reshape(y, shape),
        z,
        source_fill(source_validity(filepath, source_name)),
        zeros(Float32, size(mask)),
    )
end

"""
    prepare_boundaries(target_grid, config::AbstractBoundaryDataConfig; coverage = nothing)

Regrid the hourly source files already downloaded into `boundary_data_directory(config)` onto the
open edge rows of `target_grid`, and write the boundary NetCDF at `boundary_data_path(config)` that
`boundary_series` reads. The edges are `open_edges(config)`, so they are read from the one config that
states them rather than passed alongside that config, where the two could disagree.

Every named edge is prepared into the **same file**, which the layout already allowed for: variables
are named for their side (`south_T`, `west_T`), all six spatial dimensions are defined in every
boundary file, and one time axis serves all of them. A domain in the open ocean naming all four edges
therefore writes one file with four sides in it, not four files.

The dataset enters only through `boundary_variable_names`, `boundary_time_steps` and
`boundary_source_grid`, exactly as `prepare_forcing`'s does; everything after them is shared, and
most of it *is* `prepare_forcing`'s own machinery applied to a one-cell-thick target slab —
`water_mask`, `SourceFill`, `projected_target_nodes` and the trilinear interpolation kernel are all
reused unchanged.

`coverage` is the `(first, last)` interval the run needs, from `coverage_window`, padded to by
`pad_time_steps` at most one record spacing at each end. `nothing` prepares exactly the downloaded
range.

# Returns
A named tuple with `output_file`, `times` and `variables` (the written variable names).
"""
function prepare_boundaries(target_grid, config::AbstractBoundaryDataConfig; coverage = nothing)
    architecture = interpolation_architecture(config)
    edges = open_edges(config)

    variable_names = boundary_variable_names(config)
    source_names = [name for name in config.parameters if haskey(variable_names, name)]
    isempty(source_names) && error(
        "None of the configured parameters $(config.parameters) map to a boundary variable. " *
        "Known $(nameof(typeof(config))) variables: $(sort(collect(keys(variable_names)))).",
    )

    steps = pad_time_steps(hourly_time_steps(boundary_time_steps(config)), coverage)
    reference_file = first(steps).lower.filepath
    source = boundary_source_grid(config, reference_file)
    @info "Preparing the $(join(edges, ", ")) boundary from " *
          "$(length(unique(step -> step.lower.filepath, steps))) file(s), " *
          "$(length(steps)) time steps: $(first(steps).date) to $(last(steps).date)"

    # Every edge's variables go into one vector and one file: they are named for their side, so the
    # writer needs no notion of which edge it is on, and one time axis serves all of them.
    variables = [
        prepared_boundary_variable(name, edge, target_grid, source, reference_file, config)
        for edge in edges for name in source_names
    ]

    output_file = boundary_data_path(config)
    @info "Writing boundary file to $output_file, interpolating on $(summary(architecture))"
    write_boundaries_file(
        output_file, edges, target_grid, variables, steps, source, architecture, config,
    )
    @info "Finished preparing the $(join(edges, ", ")) boundary"

    return (;
        output_file,
        times = [step.date for step in steps],
        variables = [variable.name for variable in variables],
    )
end

prepare_boundaries(target_grid, ::Nothing; coverage = nothing) = nothing

"""
    prepare_boundaries(config::FjordConfig)

Regrid the open-boundary data a whole setup names onto its simulation grid, and write the diagnostic
plot. Returns `nothing` when the setup names no boundary data.

The setup-level driver, the same shape as `prepare_forcing(config::FjordConfig)`: the grid comes
from the processed bathymetry so the boundary row's land mask matches the model exactly, which means
`prepare_bathymetry` must have run first, and the coverage window comes from the simulation config.

# Returns
The `prepare_boundaries(target_grid, config)` named tuple with `plot_file` added.
"""
function prepare_boundaries(config::FjordConfig)
    boundaries = config.boundary_config
    isnothing(boundaries) && return nothing

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run `julia --project -m FjordSim prepare_bathymetry` for this setup first.",
    )

    grid = simulation_grid(config.grid_config, bathymetry_file, CPU())
    result = prepare_boundaries(grid, boundaries; coverage = coverage_window(config.simulation_config))
    plot_file = plot_boundaries(boundaries)

    @info "Prepared boundary variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Boundary file saved to $(result.output_file)"
    @info "Boundary plot saved to $plot_file"

    return (; result..., plot_file)
end

"""
    write_boundaries_file(filepath, edges, target_grid, variables, steps, source, architecture, config)

Write the boundary NetCDF, streaming one time step at a time so peak memory stays at a single
boundary slab. Interpolation runs on `architecture`.

Two source fields rather than `write_forcing_file`'s one: the full-depth variables share the source's
own depth axis, and the surface variables need the single-level grid `surface_source_field_grid`
builds.

`config` is carried this far only to reach `boundary_source_slab`, which is the one step of the loop a
source can override — everything else here is dataset-agnostic.
"""
function write_boundaries_file(
    filepath,
    edges,
    target_grid,
    variables,
    steps,
    source,
    architecture,
    config,
)
    isdir(dirname(filepath)) || mkpath(dirname(filepath))
    isfile(filepath) && rm(filepath; force = true)

    ds = NCDataset(filepath, "c")
    try
        define_forcing_dimensions!(ds, target_grid, [step.date for step in steps])
        ds.attrib[BOUNDARY_EDGE_ATTRIBUTE] = join(edges, ",")

        for variable in variables
            # One boundary section of one time step per chunk, matching how both this writer and
            # `load_boundary_data` touch the file: one time index at a time.
            chunksizes = [count_along_dimensions(variable)..., 1]
            defVar(ds, variable.name, Float32, variable.dimensions; chunksizes,
                deflatelevel = BOUNDARY_DEFLATE_LEVEL, attrib = ["_FillValue" => NaN32])
        end

        full_depth_field = Field{Center,Center,Center}(source_field_grid(source, architecture))
        surface_field = Field{Center,Center,Center}(surface_source_field_grid(source, architecture))

        devices = [(
            x = on_architecture(architecture, variable.x),
            y = on_architecture(architecture, variable.y),
            z = on_architecture(architecture, variable.z),
            mask = on_architecture(architecture, variable.mask),
            output = on_architecture(architecture, zeros(Float32, size(variable.mask))),
            host = Array{Float32}(undef, size(variable.mask)),
            # A surface variable is the one with no depth dimension, which is the same test as
            # `is_surface_boundary_variable` on its bare name but needs no notion of which edge this
            # variable belongs to — and this file may hold several.
            surface = length(variable.dimensions) == 2,
        ) for variable in variables]

        reader = SourceReader(first(steps).lower.filepath)
        announced = ""

        try
            for (index, step) in enumerate(steps)
                if step.lower.filepath != announced
                    announced = step.lower.filepath
                    @info "Regridding $(basename(announced))"
                end

                for (variable, device) in zip(variables, devices)
                    source_field = device.surface ? surface_field : full_depth_field
                    slab = boundary_source_slab(config, reader, step, variable.source_name)
                    set_source_field!(source_field, slab, variable.mask_fill)
                    interpolate_to_target!(
                        device.output, source_field, device.x, device.y, device.z, device.mask, architecture,
                    )
                    copyto!(device.host, device.output)
                    ds[variable.name][ntuple(_ -> Colon(), length(variable.dimensions) - 1)..., index] =
                        reshape(device.host, count_along_dimensions(variable))
                end
            end
        finally
            close(reader)
        end
    finally
        close(ds)
    end

    return filepath
end

"""
    count_along_dimensions(variable::PreparedVariable)

The variable's shape in the boundary file: its `(nx, ny, nz)` mask with the one-cell cross-edge
dimension, and a surface variable's one-level vertical dimension, dropped.
"""
count_along_dimensions(variable::PreparedVariable) =
    Tuple(extent for extent in size(variable.mask) if extent > 1)

"""
Global attribute naming the edge a boundary file was written for. The filename carries no side, so
this is where a file states its own — a reader that finds an unexpected edge has read the wrong file
rather than silently interpolating the wrong boundary.
"""
const BOUNDARY_EDGE_ATTRIBUTE = "open_edge"

"""
    unprefixed_boundary_name(name, edge)

`name` with its `"<edge>_"` prefix removed, the inverse of `boundary_variable_name`.
"""
unprefixed_boundary_name(name, edge) = replace(String(name), string(edge, "_") => ""; count = 1)

# --- Reading ---

"""
    BoundaryBackend

`FieldTimeSeries` backend for the prepared boundary file: two time indices in memory, streamed from
NetCDF, exactly like `NetCDFBackend`.

A separate type rather than a reuse of `NetCDFBackend` because the two read different ranks. A
boundary variable is stored with its cross-edge dimension dropped, so `set!` has to reshape on the
way in; making `load_from_netcdf` do that for every forcing variable too would put a rank branch in
the middle of the interior forcing's hot read path for no benefit there.
"""
struct BoundaryBackend <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
end

BoundaryBackend(length) = BoundaryBackend(1, length)
new_backend(::BoundaryBackend, start, length) = BoundaryBackend(start, length)

Base.length(backend::BoundaryBackend) = backend.length
Base.summary(backend::BoundaryBackend) = string("BoundaryBackend(", backend.start, ", ", backend.length, ")")

const BoundaryFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:BoundaryBackend}

"""
    load_boundary_data(; path, variable_name, series_size, time_indices_in_memory, reference_date)

One page of a prepared boundary variable, reshaped from the file's `(along, Nz, time)` or
`(along, time)` layout to the `series_size` of the reduced `FieldTimeSeries` that holds it, plus the
time axis in seconds from `reference_date`.

Land is `NaN` in the file and is filled here by `fill_boundary_gaps!`, so the series a boundary
condition reads is finite everywhere. That is not cosmetic — see that function.
"""
function load_boundary_data(;
    path::String,
    variable_name::String,
    series_size::Tuple,
    time_indices_in_memory,
    reference_date = nothing,
)
    ds = NCDataset(path)
    try
        variable = ds[variable_name]
        native_times = ds["time"]

        data = zeros(Float64, series_size..., length(time_indices_in_memory))
        page = 1
        for index in time_indices_in_memory
            slice = coalesce.(read_last_dimension(variable, index), NaN)
            @views data[:, :, :, page] .= reshape(slice, series_size)
            page += 1
        end

        fill_boundary_gaps!(data)

        start_time = isnothing(reference_date) ? native_times[1] : reference_date
        times = convert.(Int64, native_times_to_seconds(native_times, start_time))

        return data, times
    finally
        close(ds)
    end
end

"""
    fill_boundary_gaps!(data)

Replace every non-finite value in `data` with the nearest finite one along the boundary, per level
and per time step; zero where a whole line is dry. Mutates and returns `data`.

The prepared file marks a dry boundary cell `NaN`, as every FjordSim prepared file does. Unlike the
interior forcing, though, a boundary series is read as an *exterior value* by schemes that do
arithmetic on it directly, and only some of the nodes they write are ones the grid calls closed.

The tangential velocity is the case that forces this. On a south edge it is `u`, whose faces at
`i = 1` and `i = Nx + 1` are peripheral — the halo column outside them is inactive — so `water_mask`
marks them dry and the prepared file has a hole there. But `immersed_peripheral_node`, which is what
`NormalRadiation` zeroes on, is about the *immersed* boundary and is false at a wet domain-edge face,
so the hole would be read as a real exterior velocity and nudged into the halo at both bottom corners
of the domain. A sentinel value would be worse than a neighbouring one, and a neighbouring one is
what every other gap-filling step in this package already chooses (`SourceFill`, `daily_time_steps`).
"""
function fill_boundary_gaps!(data)
    # One of the first two axes is the reduced cross-edge one and has length 1 — `vec` flattens the
    # pair into the line along the boundary whichever it is, so a west edge needs no separate case.
    for time in axes(data, 4), level in axes(data, 3)
        line = vec(view(data, :, :, level, time))

        if !any(isfinite, line)
            line .= 0
            continue
        end

        # Two sweeps, each recording how far away the value it carries came from, so a gap in the
        # middle of the line takes whichever side is genuinely nearer rather than always the left —
        # a land segment several cells wide is otherwise filled entirely from one end of it.
        count = length(line)
        behind = fill(NaN, count)
        behind_distance = fill(count + 1, count)
        ahead = fill(NaN, count)
        ahead_distance = fill(count + 1, count)

        carried = NaN
        distance = count + 1
        for index = 1:count
            if isfinite(line[index])
                carried = line[index]
                distance = 0
            else
                distance += 1
            end
            behind[index] = carried
            behind_distance[index] = distance
        end

        carried = NaN
        distance = count + 1
        for index = count:-1:1
            if isfinite(line[index])
                carried = line[index]
                distance = 0
            else
                distance += 1
            end
            ahead[index] = carried
            ahead_distance[index] = distance
        end

        for index = 1:count
            isfinite(line[index]) && continue
            line[index] = behind_distance[index] <= ahead_distance[index] ? behind[index] : ahead[index]
        end
    end

    return data
end

"""
    set!(fts::BoundaryFTS, path, name)

Page a boundary series' in-memory window in from the file. Called by
`Oceananigans.OutputReaders.update_field_time_series!` as the run advances, which reaches these
series because `extract_field_time_series` walks the `parameters` of a discrete boundary function and
the boundary conditions of every field in `Oceananigans.fields(model)` — barotropic `U` and `V`
included.
"""
function set!(fts::BoundaryFTS, path::String = fts.path, name::String = fts.name)
    data, _ = load_boundary_data(;
        path,
        variable_name = name,
        series_size = size(fts)[1:end-1],
        time_indices_in_memory = time_indices(fts),
    )

    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    return nothing
end

"""
    boundary_series(config, grid, reference_date)

The prepared exterior state along every edge `open_edges(config)` names, as a `NamedTuple` keyed by
*edge* whose entries are `NamedTuple`s of reduced `FieldTimeSeries` keyed by the bare FjordSim names —
`(; south = (; T, S, u, v, eta, ubar, vbar), west = …)`, each restricted to what the file carries.

Keyed by edge because a domain may be open on several, and keyed on the bare name inside because the
side is already the outer key — so no boundary condition spells a side into a variable name, and the
`(; T, S, …)` group a scheme consumes is the same shape whichever edge it came from.

`reference_date` is the instant the time axes are zeroed at — the simulation config's `start_date`,
the same instant the interior forcing and the atmosphere are zeroed at, which is the only thing
keeping the three in phase. Indexing is `Cyclical()` like every other FjordSim reader, and
`validate_time_coverage` is what keeps the wrap unreachable.
"""
function boundary_series(config::AbstractBoundaryDataConfig, grid, reference_date = nothing)
    edges = open_edges(config)
    filepath = boundary_data_path(config)
    isfile(filepath) || error(
        "Prepared boundary file $filepath does not exist. " *
        "Run `julia --project -m FjordSim prepare_boundaries` for this setup first.",
    )

    validate_boundary_file_edges(filepath, edges)

    return NamedTuple(edge => boundary_edge_series(config, grid, edge, filepath, reference_date) for edge in edges)
end

"""
    validate_boundary_file_edges(filepath, edges)

Check the prepared file states every edge the setup opens.

The file records its own sides in the `open_edge` global attribute, so a setup that has since opened
another edge is caught here — before the missing variables surface as a bare `haskey` failure — and
told to re-run the step that would write them.
"""
function validate_boundary_file_edges(filepath, edges)
    stated = NCDataset(filepath) do ds
        get(ds.attrib, BOUNDARY_EDGE_ATTRIBUTE, nothing)
    end
    isnothing(stated) && return nothing

    written = Symbol.(split(stated, ","))
    missing_edges = setdiff(edges, written)
    isempty(missing_edges) || error(
        "Boundary file $filepath was written for the $(join(written, ", ")) edge(s) but this setup " *
        "opens $(join(edges, ", ")). Re-run `julia --project -m FjordSim prepare_boundaries` for it.",
    )

    return nothing
end

"""
    boundary_edge_series(config, grid, edge, filepath, reference_date)

One edge's group: the reduced `FieldTimeSeries` the file carries for that side, keyed by bare name.
"""
function boundary_edge_series(config, grid, edge, filepath, reference_date)
    names = sort(collect(values(boundary_variable_names(config))))
    present = NCDataset(filepath) do ds
        [name for name in names if haskey(ds, boundary_variable_name(edge, name))]
    end
    isempty(present) &&
        error("Boundary file $filepath carries no :$edge variables. Re-run `prepare_boundaries`.")

    backend = BoundaryBackend(2)
    time_indices_in_memory = (1, length(backend))

    series = map(present) do name
        LX, LY, LZ = boundary_location(Val(edge), name)
        variable_name = boundary_variable_name(edge, name)
        series_size = reduced_series_size(grid, (LX, LY, LZ))

        data, times = load_boundary_data(;
            path = filepath,
            variable_name,
            series_size,
            time_indices_in_memory,
            reference_date,
        )

        fts = FieldTimeSeries{LX,LY,LZ}(
            grid,
            times;
            backend,
            time_indexing = Cyclical(),
            path = filepath,
            name = variable_name,
        )
        copyto!(interior(fts, :, :, :, :), data)
        fill_halo_regions!(fts)

        Symbol(name) => fts
    end

    return NamedTuple(series)
end

boundary_series(::Nothing, grid, reference_date = nothing) = nothing

"""
    reduced_series_size(grid, locations)

Spatial shape of a reduced `FieldTimeSeries` at `locations` on `grid`: the node count per direction,
`1` where the location is `Nothing`. The same shape `interior(fts, :, :, :, :)` has, which is what
`load_boundary_data` reshapes the file's slice into.
"""
function reduced_series_size(grid, locations)
    Nx, Ny, Nz = size(grid)

    return ntuple(3) do direction
        location = locations[direction]
        interior_count = (Nx, Ny, Nz)[direction]
        location === Nothing ? 1 : (location === Face ? interior_count + 1 : interior_count)
    end
end
