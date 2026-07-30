using Oceananigans
using Oceananigans.Units
using CairoMakie
using Printf
using Statistics
using FjordSim
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain

const DEFAULT_CONFIG_PATH = joinpath(@__DIR__, "..", "configs", "drammensfjorden.jl")
# Logical figure pixels per grid cell, and fixed margin for title/labels/colorbar.
# Sized from the grid so every cell is rendered as a distinct, unambiguous block
# rather than being blurred together by too few output pixels per cell.
const BATHYMETRY_PLOT_PIXELS_PER_CELL = 6
const BATHYMETRY_PLOT_MARGIN = (300, 150)

function default_bathymetry_figure_size(grid)
    Nx, Ny, _ = size(grid)
    return (
        Nx * BATHYMETRY_PLOT_PIXELS_PER_CELL + BATHYMETRY_PLOT_MARGIN[1],
        Ny * BATHYMETRY_PLOT_PIXELS_PER_CELL + BATHYMETRY_PLOT_MARGIN[2],
    )
end

function plot_bathymetry(grid, bathymetry; plot_path, title = "Bathymetry", figure_size = default_bathymetry_figure_size(grid))
    isdir(dirname(plot_path)) || mkpath(dirname(plot_path))

    cpu_bathymetry = on_architecture(CPU(), bathymetry)
    bathymetry_data = Array(interior(cpu_bathymetry, :, :, 1))

    Nx, Ny, _ = size(grid)
    longitude = collect(range(x_domain(grid)[1], x_domain(grid)[2], length = Nx))
    latitude = collect(range(y_domain(grid)[1], y_domain(grid)[2], length = Ny))

    figure = Figure(size = figure_size)
    axis = Axis(figure[1, 1]; xlabel = "Longitude", ylabel = "Latitude", title)

    plot = heatmap!(axis, longitude, latitude, bathymetry_data; colormap = :deep, colorrange = extrema(bathymetry_data))
    land_mask = ifelse.(bathymetry_data .>= 0, 1.0f0, NaN32)
    heatmap!(
        axis,
        longitude,
        latitude,
        land_mask;
        colormap = [:ivory, :ivory],
        colorrange = (0, 1),
        nan_color = RGBAf(0, 0, 0, 0),
    )
    Colorbar(figure[1, 2], plot; label = "Bottom height (m)")
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = -collect(25:25:300),
        color = (:white, 0.35),
        linewidth = 1,
    )
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = [0.0],
        color = :black,
        linewidth = 4,
    )
    save(plot_path, figure)

    return plot_path
end

function parse_args(args = ARGS)
    config_path = DEFAULT_CONFIG_PATH

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
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    return (; config_path)
end

function print_usage()
    println("""
    Prepare Geonorge bathymetry for a configured FjordSim setup.

    Usage:
      julia --project scripts/bathymetry_prepare.jl [--config PATH]

    Options:
      --config PATH   Setup config (Julia file). Default: configs/drammensfjorden.jl
    """)
end

function print_grid_extents(grid)
    dx = xspacings(grid)
    dy = yspacings(grid)
    dz = zspacings(grid)

    println("Grid cell side extents (m):")
    @printf("  N = (%d, %d, %d)\n", grid.Nx, grid.Ny, grid.Nz)
    @printf("  Δx min/max/mean = %.3f / %.3f / %.3f\n", minimum(dx), maximum(dx), mean(dx))
    @printf("  Δy min/max/mean = %.3f / %.3f / %.3f\n", minimum(dy), maximum(dy), mean(dy))
    @printf("  Δz min/max/mean = %.3f / %.3f / %.3f\n", minimum(dz), maximum(dz), mean(dz))
end

function main()
    args = parse_args()
    config = include(abspath(args.config_path))

    grid = LatitudeLongitudeGrid(CPU(), config.grid)
    mkpath(dirname(config.bathymetry.output_path))
    print_grid_extents(grid)

    result = prepare_geonorge_bathymetry(grid, config.bathymetry)

    plot_bathymetry(grid, result.bottom_height; plot_path = config.bathymetry.plot_path)

    @info "Raw Geonorge bathymetry saved to $(result.raw_path)"
    @info "Processed FjordSim bathymetry saved to $(result.output_path)"
    @info "Bathymetry plot saved to $(config.bathymetry.plot_path)"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
