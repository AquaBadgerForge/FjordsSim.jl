#!/usr/bin/env julia

using Oceananigans
using CUDA  # needed for GPU() to be usable
using FjordSim

include("cli.jl")

const USAGE = """
Regrid downloaded forcing onto a configured FjordSim setup's grid.

Reads the source files written by scripts/forcing_download.jl and the processed bathymetry
written by scripts/bathymetry_prepare.jl, so run both first.

Usage:
  julia --project scripts/forcing_prepare.jl --config PATH

Options:
  --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
  --cpu           Interpolate on the CPU even when a GPU is available
"""

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

function main()
    args = parse_config_args(USAGE; flags = ("--cpu",))
    config = include(abspath(args.config_path))

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run scripts/bathymetry_prepare.jl --config $(args.config_path) first.",
    )

    # The grid stays on the CPU: building the land masks walks `peripheral_node` cell by cell.
    grid = ImmersedBoundaryGrid(bathymetry_file, CPU(), config.grid_config.halo)
    architecture = interpolation_architecture("--cpu" in args.flags)
    result = prepare_forcing(grid, config.forcing_config; architecture)
    plot_file = plot_forcing(grid, config.forcing_config)

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing file saved to $(result.output_file)"
    @info "Forcing plot saved to $plot_file"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
