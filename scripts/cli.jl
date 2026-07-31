# Shared command-line handling for the scripts in this directory. Each script owns a usage
# string and its own `main`; everything else about the CLI is here.
#
# The `abspath(PROGRAM_FILE) == @__FILE__()` run guard cannot live here — `@__FILE__` would
# expand to this file's path — so each script keeps its own three-line guard.

"""
    parse_config_args(usage; args = ARGS, flags = ())

Parse `--config PATH` (or `--config=PATH`), `-h`/`--help` and any option named in `flags`.

Returns `(; config_path, flags)`, where `flags` is the set of flag names present. Prints
`usage` and exits on `--help`; errors on an unknown argument or a missing `--config`.
"""
function parse_config_args(usage; args = ARGS, flags = ())
    config_path = nothing
    present = Set{String}()

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
        elseif arg in flags
            push!(present, arg)
            i += 1
        elseif arg in ("-h", "--help")
            println(usage)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isnothing(config_path) && error("--config PATH is required")

    return (; config_path, flags = present)
end
