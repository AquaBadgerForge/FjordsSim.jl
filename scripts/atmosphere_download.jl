#!/usr/bin/env julia

using FjordSim

include("cli.jl")

const USAGE = """
Download and subset the atmosphere dataset a configured FjordSim setup names.

This is the network-bound half of the atmosphere pipeline and by far the slower one: NORA3 is
served one file per forecast lead hour, so a year is close to 10000 OPeNDAP reads. A month whose
file already exists is skipped, so an interrupted run resumes where it stopped.

Usage:
  julia --project scripts/atmosphere_download.jl --config PATH

Options:
  --config PATH   Setup config (Julia file). Required, e.g. configs/oslofjorden.jl
"""

function main()
    args = parse_config_args(USAGE)
    config = include(abspath(args.config_path))

    isnothing(config.atmosphere_config) &&
        error("$(abspath(args.config_path)) names no atmosphere_config, so there is nothing to download.")

    @info "Config: $(abspath(args.config_path))"
    output_directory = download_atmosphere(config)

    @info "Atmosphere source files saved to $output_directory"
end

if abspath(PROGRAM_FILE) == @__FILE__() || get(ENV, "FJORDSIM_RUN_MAIN", "") == "1"
    main()
end
