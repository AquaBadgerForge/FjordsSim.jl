#!/usr/bin/env julia

using Oceananigans
using FjordSim

include("cli.jl")

const USAGE = """
Write river relaxation on top of a configured FjordSim setup's prepared forcing.

Reads the forcing file written by scripts/forcing_prepare.jl and the processed bathymetry
written by scripts/bathymetry_prepare.jl, so run both first. The result is a copy of the
forcing file with the river cells patched, leaving the original untouched.

Only a setup whose forcing config has a `rivers` field is supported; `rivers = nothing` means
the setup has no rivers.

Usage:
  julia --project scripts/forcing_add_rivers.jl --config PATH

Options:
  --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
"""

function main()
    args = parse_config_args(USAGE)
    config = include(abspath(args.config_path))

    rivers = config.forcing_config.rivers
    isnothing(rivers) && error(
        "No rivers are configured for $(args.config_path). " *
        "Set the forcing config's `rivers` field to an AbstractRiverConfig to add them.",
    )

    bathymetry_file = bathymetry_path(config.bathymetry_config)
    isfile(bathymetry_file) || error(
        "Processed bathymetry $bathymetry_file does not exist. " *
        "Run scripts/bathymetry_prepare.jl --config $(args.config_path) first.",
    )

    forcing_file = forcing_path(config.forcing_config)
    isfile(forcing_file) || error(
        "Prepared forcing $forcing_file does not exist. " *
        "Run scripts/forcing_prepare.jl --config $(args.config_path) first.",
    )

    download_rivers(rivers)

    # The grid stays on the CPU: building the land masks walks `peripheral_node` cell by cell.
    grid = ImmersedBoundaryGrid(bathymetry_file, CPU(), config.grid_config.halo)
    result = add_rivers(grid, config.forcing_config)

    @info "Placed $(length(result.cells)) of $(length(river_locations(rivers))) rivers"
    @info "Patched variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Forcing with rivers saved to $(result.output_file)"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
