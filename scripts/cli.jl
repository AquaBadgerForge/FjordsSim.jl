# Shared command-line handling for the scripts in this directory. Each script owns a usage
# string and its own `main`; everything else about the CLI is here.
#
# The `abspath(PROGRAM_FILE) == @__FILE__()` run guard cannot live here — `@__FILE__` would
# expand to this file's path — so each script keeps its own three-line guard.

"""
    parse_config_args(usage; args = ARGS)

Parse `--config PATH` (or `--config=PATH`) and `-h`/`--help`.

Returns `(; config_path)`. Prints `usage` and exits on `--help`; errors on an unknown argument
or a missing `--config`. Every knob other than which setup to run belongs in the setup config,
so this is the whole command-line surface.
"""
function parse_config_args(usage; args = ARGS)
    config_path = nothing

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
            println(usage)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isnothing(config_path) && error("--config PATH is required")

    return (; config_path)
end
