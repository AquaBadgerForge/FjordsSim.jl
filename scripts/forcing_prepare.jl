#!/usr/bin/env julia

using Oceananigans
using FjordSim

include("cli.jl")

const USAGE = """
Regrid downloaded forcing onto a configured FjordSim setup's grid.

Reads the source files written by scripts/forcing_download.jl and the processed bathymetry
written by scripts/bathymetry_prepare.jl, so run both first.

Where the interpolation runs is the setup's `architecture` field (`:auto`, `:cpu` or `:gpu`),
not a command-line option.

Usage:
  julia --project scripts/forcing_prepare.jl --config PATH

Options:
  --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
"""

function main()
    args = parse_config_args(USAGE)
    config = include(abspath(args.config_path))

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run scripts/bathymetry_prepare.jl --config $(args.config_path) first.",
    )

    # The grid stays on the CPU: building the land masks walks `peripheral_node` cell by cell.
    grid = ImmersedBoundaryGrid(bathymetry_file, CPU(), config.grid_config.halo)
    result = prepare_forcing(grid, config.forcing_config)
    plot_file = plot_forcing(grid, config.forcing_config)

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing file saved to $(result.output_file)"
    @info "Forcing plot saved to $plot_file"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
