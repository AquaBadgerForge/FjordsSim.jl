#!/usr/bin/env julia

using FjordSim

include("cli.jl")

const USAGE = """
Download and subset the forcing dataset a configured FjordSim setup names.

Usage:
  julia --project scripts/forcing_download.jl --config PATH

Options:
  --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
"""

function main()
    args = parse_config_args(USAGE)
    config = include(abspath(args.config_path))

    @info "Config: $(abspath(args.config_path))"
    output_directory = download_forcing(config)

    @info "Forcing source files saved to $output_directory"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
