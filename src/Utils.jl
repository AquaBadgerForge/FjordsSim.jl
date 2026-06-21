module Utils

export compute_faces,
    progress,
    progress_callback,
    safe_execute,
    extract_z_faces,
    netcdf_to_jld2,
    save_fts,
    recursive_merge,
    cell_advection_timescale_coupled_model

using Oceananigans
using Oceananigans.Fields: interior
using Oceananigans.OutputReaders: FieldTimeSeries, OnDisk
using Oceananigans.Utils: prettytime
using Oceananigans.Advection: cell_advection_timescale
using JLD2: @save
using Printf: @sprintf

cell_advection_timescale_coupled_model(coupled_model) = cell_advection_timescale(coupled_model.ocean.model)

function compute_faces(centers)
    spacing = diff(centers)[1]  # Assuming uniform spacing
    faces = vcat([centers[1] - spacing / 2], (centers[1:end-1] .+ centers[2:end]) / 2, [centers[end] + spacing / 2])
    return faces
end

const wall_time = Ref(time_ns())

_progress_tracers(tracers::Symbol) = (tracers,)

function _progress_tracers(tracers)
    tracer_names = Tuple(tracers)
    all(tracer -> tracer isa Symbol, tracer_names) ||
        throw(ArgumentError("progress tracers must be Symbols."))
    return tracer_names
end

function _tracer_extrema_message(model_tracers, tracer_name)
    if !hasproperty(model_tracers, tracer_name)
        available_tracers = join([":" * String(name) for name in propertynames(model_tracers)], ", ")
        throw(ArgumentError("Tracer :$tracer_name is not defined. Available tracers: $available_tracers."))
    end

    tracer = getproperty(model_tracers, tracer_name)
    tracer_max = maximum(interior(tracer))
    tracer_min = minimum(interior(tracer))
    units = tracer_name === :T ? " ᵒC" : ""

    return @sprintf("extrema(%s): (%.2f, %.2f)%s", String(tracer_name), tracer_max, tracer_min, units)
end

function progress(sim; tracers = (:T,), wall_time_reference = wall_time)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities

    umax = (maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time_reference[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹", umax...)

    tracer_names = _progress_tracers(tracers)
    for tracer_name in tracer_names
        msg *= ", " * _tracer_extrema_message(ocean.model.tracers, tracer_name)
    end

    msg *= @sprintf(", wall time: %s", prettytime(step_time))

    @info msg

    wall_time_reference[] = time_ns()
end

function progress_callback(; tracers = (:T,))
    wall_time_reference = Ref(time_ns())
    return sim -> progress(sim; tracers, wall_time_reference)
end

function safe_execute(callable)
    return function (args...)
        if callable === nothing || args === nothing
            return nothing
        elseif isa(callable, Function)
            return callable(args...)
        else
            return nothing
        end
    end
end

function extract_z_faces(grid)
    bar = grid["zᵃᵃᶜ"]
    zero_index = findfirst(x -> x > 0.0, bar)
    n = grid["Nz"] + 1
    if zero_index > 1
        start_index = max(1, zero_index - n)
        z = bar[start_index:zero_index-1]
    else
        z = Float64[]
    end
    return z
end

function netcdf_to_jld2(netcdf_file::String, jld2_file::String)
    ds = NCDataset(netcdf_file, "r")
    data_dict = Dict()
    for varname in keys(ds)
        arr = convert(Array, ds[varname])
        data_dict[varname] = arr
        print(size(arr))
    end

    @save jld2_file data_dict
    close(ds)
    println("Conversion completed: NetCDF to JLD2")
end

function save_fts(; jld2_filepath, fts_name, fts, grid, times, boundary_conditions)
    isfile(jld2_filepath) && rm(jld2_filepath)
    on_disk_fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        boundary_conditions,
        backend = OnDisk(),
        path = jld2_filepath,
        name = fts_name,
    )
    for i = 1:size(fts)[end]
        set!(on_disk_fts, fts[i], i, times[i])
    end
end

function recursive_merge(nt1::NamedTuple, nt2::NamedTuple)
    # Get all unique keys from both NamedTuples
    all_keys = union(keys(nt1), keys(nt2))

    # Initialize an empty NamedTuple for the result
    result_pairs = Pair{Symbol,Any}[]

    for key in all_keys
        val1 = get(nt1, key, nothing)
        val2 = get(nt2, key, nothing)

        if val1 isa NamedTuple && val2 isa NamedTuple
            # If both values are NamedTuples, recursively merge them
            push!(result_pairs, key => recursive_merge(val1, val2))
        elseif val2 !== nothing
            # If only val2 exists or is not a NamedTuple, use val2
            push!(result_pairs, key => val2)
        else
            # Otherwise, use val1 (if it exists)
            push!(result_pairs, key => val1)
        end
    end
    return (; result_pairs...)
end

end # module
