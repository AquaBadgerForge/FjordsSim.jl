#!/usr/bin/env julia

using Oceananigans
using FjordSim

include("cli.jl")

const USAGE = """
Regrid downloaded atmosphere data onto a regular longitude/latitude grid covering a configured
FjordSim setup, writing the single file the simulation reads.

Reads the files written by scripts/atmosphere_download.jl, so run that first. It needs no
bathymetry: the atmosphere grid is derived from the setup's longitude/latitude domain, not from the
ocean grid's land mask, so unlike forcing_prepare.jl this does not depend on
scripts/bathymetry_prepare.jl.

This step is offline and takes minutes, so it is cheap to re-run after changing `resolution` or
`padding` — no re-download is needed.

Usage:
  julia --project scripts/atmosphere_prepare.jl --config PATH

Options:
  --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
"""

function main()
    args = parse_config_args(USAGE)
    config = include(abspath(args.config_path))

    atmosphere_config = config.atmosphere_config
    isnothing(atmosphere_config) &&
        error("$(abspath(args.config_path)) names no atmosphere_config, so there is nothing to prepare.")

    source_directory = atmosphere_directory(atmosphere_config)
    isdir(source_directory) || error(
        "Atmosphere source directory $source_directory does not exist. " *
        "Run scripts/atmosphere_download.jl --config $(args.config_path) first.",
    )

    grid = LatitudeLongitudeGrid(CPU(), config.grid_config)
    result = prepare_atmosphere(grid, atmosphere_config)
    plot_file = plot_atmosphere(atmosphere_config)

    @info "Prepared variables: $(join(result.variables, ", "))"
    @info "Time range: $(first(result.times)) to $(last(result.times)) ($(length(result.times)) steps)"
    @info "Atmosphere file saved to $(result.output_file)"
    @info "Atmosphere plot saved to $plot_file"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
