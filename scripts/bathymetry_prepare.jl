#!/usr/bin/env julia

using Oceananigans
using Printf
using Statistics
using FjordSim

include("cli.jl")

const USAGE = """
Prepare bathymetry for a configured FjordSim setup.

Usage:
  julia --project scripts/bathymetry_prepare.jl --config PATH

Options:
  --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
"""

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
    args = parse_config_args(USAGE)
    config = include(abspath(args.config_path))

    grid = LatitudeLongitudeGrid(CPU(), config.grid_config)
    mkpath(dirname(bathymetry_path(config.bathymetry_config)))
    print_grid_extents(grid)

    result = prepare_bathymetry(grid, config.bathymetry_config)

    plot_file = plot_bathymetry(grid, result.bottom_height, config.bathymetry_config)

    @info "Raw bathymetry saved to $(result.raw_file)"
    @info "Processed FjordSim bathymetry saved to $(result.output_file)"
    @info "Bathymetry plot saved to $plot_file"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
