#!/usr/bin/env julia

using Oceananigans
using CairoMakie
using CUDA  # needed for GPU() to be usable
using NCDatasets
using FjordSim
using Oceananigans.Grids: x_domain, y_domain

# Logical figure pixels per grid cell and fixed margin, sized from the grid so each cell is a
# distinct block, following scripts/bathymetry_prepare.jl.
const FORCING_PLOT_PIXELS_PER_CELL = 3
const FORCING_PLOT_MARGIN = (200, 150)

# Surface fields worth eyeballing after preparation, with the colormap suiting each.
const FORCING_PLOT_PANELS = (
    (name = "T", label = "Temperature (ᵒC)", colormap = :thermal),
    (name = "S", label = "Salinity (g kg⁻¹)", colormap = :haline),
    (name = "u", label = "u (m s⁻¹)", colormap = :balance),
    (name = "v", label = "v (m s⁻¹)", colormap = :balance),
    (name = "T_lambda", label = "T relaxation rate (s⁻¹)", colormap = :viridis),
)

function default_forcing_figure_size(grid, panel_count)
    Nx, Ny, _ = size(grid)
    return (
        panel_count * (Nx * FORCING_PLOT_PIXELS_PER_CELL + FORCING_PLOT_MARGIN[1]),
        Ny * FORCING_PLOT_PIXELS_PER_CELL + FORCING_PLOT_MARGIN[2],
    )
end

"""
    surface_slice(ds, name)

The top vertical level of `name` at the first time step, as a `Float32` array with missing
values turned into `NaN` so Makie leaves land blank.
"""
function surface_slice(ds, name)
    Nz = ds.dim["Nz"]
    return Float32.(coalesce.(ds[name][:, :, Nz, 1], NaN32))
end

function plot_forcing(grid, config)
    forcing_file = forcing_path(config)
    plot_file = forcing_plot_path(config)
    isdir(dirname(plot_file)) || mkpath(dirname(plot_file))

    NCDataset(forcing_file) do ds
        panels = [panel for panel in FORCING_PLOT_PANELS if haskey(ds, panel.name)]
        figure = Figure(size = default_forcing_figure_size(grid, length(panels)))
        date = ds["time"][1]

        for (column, panel) in enumerate(panels)
            data = surface_slice(ds, panel.name)
            Nx, Ny = size(data)
            longitude = collect(range(x_domain(grid)[1], x_domain(grid)[2], length = Nx))
            latitude = collect(range(y_domain(grid)[1], y_domain(grid)[2], length = Ny))

            axis = Axis(
                figure[1, column];
                xlabel = "Longitude",
                ylabel = column == 1 ? "Latitude" : "",
                title = panel.name,
            )
            finite = filter(isfinite, data)
            colorrange = isempty(finite) || allequal(finite) ? (0, 1) : extrema(finite)
            plot = heatmap!(
                axis,
                longitude,
                latitude,
                data;
                colormap = panel.colormap,
                colorrange,
                nan_color = :ivory,
            )
            Colorbar(figure[2, column], plot; label = panel.label, vertical = false)
        end

        Label(figure[0, :], "Prepared NorKyst forcing, surface level at $date", fontsize = 20)
        save(plot_file, figure)
    end

    return plot_file
end

function parse_args(args = ARGS)
    config_path = nothing
    force_cpu = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i == length(args) && error("--config requires a value")
            config_path = args[i + 1]
            i += 2
        elseif startswith(arg, "--config=")
            config_path = split(arg, "=", limit = 2)[2]
            i += 1
        elseif arg == "--cpu"
            force_cpu = true
            i += 1
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isnothing(config_path) && error("--config PATH is required")

    return (; config_path, force_cpu)
end

"""
    interpolation_architecture(force_cpu)

Where to run the interpolation kernel: the GPU when one is usable, unless `--cpu` was passed.
Only the kernel moves; the grid and masks stay on the CPU.
"""
function interpolation_architecture(force_cpu)
    if force_cpu
        @info "Interpolating on the CPU (--cpu)"
        return CPU()
    elseif CUDA.functional()
        @info "Interpolating on the GPU ($(CUDA.name(CUDA.device())))"
        return GPU()
    else
        @info "No usable GPU found; interpolating on the CPU"
        return CPU()
    end
end

function print_usage()
    println("""
    Regrid downloaded NorKyst-800m forcing onto a configured FjordSim setup's grid.

    Reads the monthly files written by scripts/forcing_download_norkyst.jl and the processed
    bathymetry written by scripts/bathymetry_prepare.jl, so run both first.

    Usage:
      julia --project scripts/forcing_prepare.jl --config PATH

    Options:
      --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
      --cpu           Interpolate on the CPU even when a GPU is available
    """)
end

function main()
    args = parse_args()
    config = include(abspath(args.config_path))

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run scripts/bathymetry_prepare.jl --config $(args.config_path) first.",
    )

    # The grid stays on the CPU: building the land masks walks `peripheral_node` cell by cell.
    grid = ImmersedBoundaryGrid(bathymetry_file, CPU(), config.grid_config.halo)
    architecture = interpolation_architecture(args.force_cpu)
    result = prepare_norkyst_forcing(grid, config.forcing_config; architecture)
    plot_file = plot_forcing(grid, config.forcing_config)

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing file saved to $(result.output_file)"
    @info "Forcing plot saved to $plot_file"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
