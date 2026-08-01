module CLI

export parse_arguments

using ..Configs: FjordConfig
using ..Setups: fjord_config, setup_names
using ..Bathymetry: prepare_bathymetry
using ..Atmospheres: prepare_atmosphere, download_atmosphere
using ..Forcing: prepare_forcing, download_forcing, add_rivers
using ..Simulations: run_simulation

"""
Every subcommand, named exactly like the function it calls, so `-m FjordSim prepare_forcing` and
`prepare_forcing(config)` in the REPL are the same thing spelled the same way.

Ordered by the sequence a setup is prepared and run in, which is the order `USAGE` lists them in.
"""
const SUBCOMMANDS = [
    "prepare_bathymetry" => prepare_bathymetry,
    "download_forcing" => download_forcing,
    "prepare_forcing" => prepare_forcing,
    "add_rivers" => add_rivers,
    "download_atmosphere" => download_atmosphere,
    "prepare_atmosphere" => prepare_atmosphere,
    "run_simulation" => run_simulation,
]

"""
    subcommand_names()

The subcommand names, in pipeline order.
"""
subcommand_names() = [name for (name, _) in SUBCOMMANDS]

"""
    subcommand_driver(name)

The driver function `name` runs, or an error listing the known subcommands.
"""
function subcommand_driver(name)
    index = findfirst(pair -> pair.first == name, SUBCOMMANDS)
    isnothing(index) && throw(
        ArgumentError("Unknown subcommand \"$name\". Available: $(join(subcommand_names(), ", "))."),
    )

    return SUBCOMMANDS[index].second
end

const USAGE = """
Prepare the input data for a FjordSim setup, and run it.

Usage:
  julia --project -m FjordSim SUBCOMMAND --config SETUP

Subcommands, in the order a setup is prepared and run:
  prepare_bathymetry    Regrid the bathymetry source onto the setup's grid. Downloads the Geonorge
                        Sjøkart FileGDB on first use (~2.3 GB).
  download_forcing      Download and subset the forcing dataset over the setup's region and years.
  prepare_forcing       Regrid the download onto the simulation grid. Needs prepare_bathymetry.
                        Where the interpolation runs is the forcing config's `architecture` field,
                        not a command-line option.
  add_rivers            Write river relaxation into a copy of the prepared forcing, leaving the
                        original untouched. A setup with no rivers does nothing.
  download_atmosphere   Download and subset the atmosphere dataset. By far the slowest step: NORA3
                        is served one file per forecast lead hour, so a year is close to 10000
                        OPeNDAP reads. A month already downloaded is skipped, so an interrupted run
                        resumes. A setup with no atmosphere does nothing.
  prepare_atmosphere    Regrid the download onto a regular longitude/latitude grid. Needs
                        download_atmosphere but not prepare_bathymetry, and is cheap to re-run
                        after changing `resolution` or `padding`.
  run_simulation        Build and run the coupled simulation, writing NetCDF snapshots into the
                        setup's results directory. Needs every step above that the setup
                        configures: prepare_bathymetry, prepare_forcing, add_rivers if it names
                        rivers, and prepare_atmosphere if it names an atmosphere. Where it runs
                        is the simulation config's `architecture` field. A setup with no
                        simulation config does nothing.

Options:
  --config SETUP   Which setup to prepare. Required. One of: $(join(setup_names(), ", ")).
                   A path ending in .jl is loaded as an out-of-tree config file instead.
  -h, --help       Show this message.

Everything other than which setup to prepare is stated in the setup itself, so this is the whole
command-line surface. The same steps are available from the REPL:

  using FjordSim
  config = $(first(setup_names()))()
  prepare_bathymetry(config)
"""

"""
    parse_arguments(args)

Parse one positional subcommand plus `--config SETUP` (or `--config=SETUP`) and `-h`/`--help`.

Returns `(; subcommand, config, help)`. `--help` sets `help` and leaves the rest `nothing` rather
than printing and exiting, so printing stays in `main` and this function is testable. Invalid
input throws `ArgumentError`.
"""
function parse_arguments(args)
    subcommand = nothing
    config = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            return (; subcommand = nothing, config = nothing, help = true)
        elseif arg == "--config"
            i == length(args) && throw(ArgumentError("--config requires a value"))
            config = args[i + 1]
            i += 2
        elseif startswith(arg, "--config=")
            config = split(arg, "=", limit = 2)[2]
            i += 1
        elseif startswith(arg, "-")
            throw(ArgumentError("Unknown option: $arg"))
        else
            isnothing(subcommand) ||
                throw(ArgumentError("Expected one subcommand, got \"$subcommand\" and \"$arg\""))
            subcommand = arg
            i += 1
        end
    end

    isnothing(subcommand) &&
        throw(ArgumentError("A subcommand is required. Available: $(join(subcommand_names(), ", "))."))
    (isnothing(config) || isempty(config)) && throw(ArgumentError("--config SETUP is required"))

    return (; subcommand, config, help = false)
end

"""
    main(args)

Run one subcommand and return a process exit code: 0 on success, 2 for bad arguments.

A step a setup opts out of — `add_rivers` on a setup with no rivers, the atmosphere steps on a
setup with no atmosphere — returns `nothing` from the driver and is reported here rather than
raising, because that is a property of the setup and not a failure.
"""
function main(args)
    arguments = try
        parse_arguments(args)
    catch exception
        exception isa ArgumentError || rethrow()
        println(stderr, "fjordsim: ", exception.msg)
        println(stderr, "Try `julia --project -m FjordSim --help`.")
        return 2
    end

    if arguments.help
        print(USAGE)
        return 0
    end

    driver, config = try
        subcommand_driver(arguments.subcommand), fjord_config(arguments.config)
    catch exception
        exception isa ArgumentError || rethrow()
        println(stderr, "fjordsim: ", exception.msg)
        return 2
    end

    @info "Setup: $(arguments.config), step: $(arguments.subcommand)"
    isnothing(driver(config)) &&
        @info "$(arguments.subcommand) is a no-op for $(arguments.config): the setup does not configure it"

    return 0
end

end  # module CLI
