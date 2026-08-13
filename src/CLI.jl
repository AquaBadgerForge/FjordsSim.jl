module CLI

export parse_arguments

using ..Configs: AbstractSimulationConfig, run_tag
using ..Setups: fjord_config, setup_names
using ..Bathymetry: prepare_bathymetry
using ..Atmospheres: prepare_atmosphere, download_atmosphere
using ..Forcing:
    prepare_forcing, download_forcing, add_rivers, download_boundaries, prepare_boundaries
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
    "download_boundaries" => download_boundaries,
    "prepare_boundaries" => prepare_boundaries,
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
  download_boundaries   Download and subset hourly source data along the open lateral boundary. A
                        thin band rather than the whole domain, since the exterior state an open
                        boundary needs has to resolve the tide and daily means do not. A setup
                        whose forcing config names no `boundaries` does nothing.
  prepare_boundaries    Regrid the download onto the open edge of the simulation grid. Needs
                        prepare_bathymetry and download_boundaries.
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
                        rivers, prepare_boundaries if it names open-boundary data, and
                        prepare_atmosphere if it names an atmosphere. Where it runs
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
    log_path(config)

Where to write this run's transcript: `fjordsim_<run_tag>.log` under the setup's `results_root`, so a
log lands beside the output it describes rather than wherever the command happened to be run from,
and so each launch keeps its own transcript.

Takes the *simulation* config, because `results_root` and the run tag are its fields and only
`run_simulation` is logged: that is the one step whose failure is a stacktrace through the whole
coupled model, and the one step a setup without a simulation config does not run at all.
"""
log_path(config::AbstractSimulationConfig) =
    joinpath(config.results_root, string("fjordsim_", run_tag(config), ".log"))

"""
The width, in columns, the error reporter tells `showerror` it has.

`Base.type_depth_limit` floors this at 120, so a smaller number changes nothing. Stated explicitly
rather than left to `displaysize` so the log does not depend on how wide the terminal was — under the
tee `stdout` is a `Pipe`, which reports the 24x80 default anyway.
"""
const STACKTRACE_WIDTH = 120

"""
    show_compact_error(io, exception, backtrace)

Report `exception` and its backtrace with the type parameters in each frame abbreviated to `{…}`.

Julia only abbreviates them when the IO carries `:stacktrace_types_limited`, which is what
`REPL.repl_display_error` sets and a script's bare `showerror` does not — which is why a stacktrace
through the coupled model spells out every parameter of every `HydrostaticFreeSurfaceModel` and
`ImmersedBoundaryGrid` in it, and stays unreadable even in a file.

It shortens *frames*, not messages: a `GPUCompiler` `InvalidIRError` body is as long as it ever was.
"""
function show_compact_error(io, exception, backtrace)
    abbreviated = Ref(false)
    context = IOContext(
        io,
        :limit => true,
        :displaysize => (first(displaysize(io)), STACKTRACE_WIDTH),
        :stacktrace_types_limited => abbreviated,
    )

    showerror(context, exception, backtrace)
    println(io)
    abbreviated[] && println(io, "fjordsim: type parameters above were abbreviated to `{…}`.")

    return nothing
end

"""
    tee_output(f, log_file)

Run `f` with `stdout` and `stderr` redirected into `log_file`, still printing everything live to the
real `stdout`, and return whatever `f` returns.

`redirect_stdio` only accepts fd-backed streams, so the tee is a `Pipe` plus a task copying each
chunk to both destinations. Both streams share one pipe, so the log interleaves them in the order
they were written.

`log_file`'s directory is created if absent, since a setup's results directory need not exist
before its first run.
"""
function tee_output(f, log_file)
    original_stdout = stdout
    directory = dirname(log_file)
    isempty(directory) || mkpath(directory)

    return open(log_file, "w") do log_stream
        pipe = Pipe()
        Base.link_pipe!(pipe)

        reader = Threads.@spawn while !eof(pipe)
            chunk = readavailable(pipe)
            write(original_stdout, chunk)
            write(log_stream, chunk)
            flush(log_stream)
        end

        try
            return redirect_stdio(f; stdout = pipe, stderr = pipe)
        finally
            close(pipe.in)
            wait(reader)
            close(pipe)
        end
    end
end

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
    run_step(driver, config, subcommand, setup)

Run one driver and return its exit code: 0 on success, 1 if it threw.

A step a setup opts out of — `add_rivers` on a setup with no rivers, the atmosphere steps on a setup
with no atmosphere — returns `nothing` from the driver and is reported here rather than raising,
because that is a property of the setup and not a failure.

A failure is caught rather than propagated: under `tee_output` an exception left to `Base._start`
would be printed after the redirect had already been torn down, so the error — the one thing worth
having in the log — would be the only thing missing from it.

`driver` runs through `Base.invokelatest` because `config` may carry types and methods that did not
exist until `fjord_config` loaded them: an out-of-tree `--config path.jl` file is `Base.include`d
from within this same `main` call, and any new struct or method it defines is invisible to code
already running in the world active when `main` started — a plain `driver(config)` would throw
"method too new" the first time a config file adds a new dispatch, e.g. a new
`AbstractBoundaryConditionConfig` subtype, rather than only new field values.
"""
function run_step(driver, config, subcommand, setup)
    @info "Setup: $setup, step: $subcommand"

    try
        isnothing(Base.invokelatest(driver, config)) &&
            @info "$subcommand is a no-op for $setup: the setup does not configure it"
        return 0
    catch exception
        println(stderr, "fjordsim: $subcommand failed on $setup:")
        show_compact_error(stderr, exception, catch_backtrace())
        return 1
    end
end

"""
    main(args)

Run one subcommand and return a process exit code: 0 on success, 1 if the step failed, 2 for bad
arguments.

`run_simulation` runs inside `tee_output`, writing to `log_path`; every other step prints to the
terminal only. It is the one step long enough and loud enough to be worth a transcript — a stacktrace
through the coupled model scrolls its own error message away — the one whose output is not
reproducible by re-running it, and the only one whose setup is guaranteed to name a `results_root` to
put the transcript in.

Parsing, `--help` and config resolution stay outside the tee, so a usage error leaves no log file
behind — which also means the log path is only ever resolved from a config that loaded.
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

    simulation_config = config.simulation_config
    if driver === run_simulation && !isnothing(simulation_config)
        log_file = log_path(simulation_config)
        @info "Logging to $(abspath(log_file))"

        return tee_output(log_file) do
            run_step(driver, config, arguments.subcommand, arguments.config)
        end
    end

    return run_step(driver, config, arguments.subcommand, arguments.config)
end

end  # module CLI
