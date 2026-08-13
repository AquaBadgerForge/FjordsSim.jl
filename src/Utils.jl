module Utils

export compute_faces,
    progress,
    safe_execute,
    extract_z_faces,
    netcdf_to_jld2,
    save_fts,
    recursive_merge,
    cell_advection_timescale_coupled_model

using Oceananigans
using Oceananigans.Fields: interior, location
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

const WALL_TIME = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))

    umax = (maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - WALL_TIME[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(
        ", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
        umax...,
        Tmax,
        Tmin,
        prettytime(step_time)
    )

    @info msg

    WALL_TIME[] = time_ns()
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
    for variable_name in keys(ds)
        array = convert(Array, ds[variable_name])
        data_dict[variable_name] = array
        print(size(array))
    end

    @save jld2_file data_dict
    close(ds)
    println("Conversion completed: NetCDF to JLD2")
end

function save_fts(; jld2_filepath, fts_name, fts, grid, times, boundary_conditions)
    isfile(jld2_filepath) && rm(jld2_filepath)
    LX, LY, LZ = location(fts)
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

function recursive_merge(a::NamedTuple, b::NamedTuple)
    # Get all unique keys from both NamedTuples
    all_keys = union(keys(a), keys(b))

    # Initialize an empty NamedTuple for the result
    result_pairs = Pair{Symbol,Any}[]

    for key in all_keys
        a_value = get(a, key, nothing)
        b_value = get(b, key, nothing)

        if a_value isa NamedTuple && b_value isa NamedTuple
            # If both values are NamedTuples, recursively merge them
            push!(result_pairs, key => recursive_merge(a_value, b_value))
        elseif b_value !== nothing
            # If only b_value exists or is not a NamedTuple, use b_value
            push!(result_pairs, key => b_value)
        else
            # Otherwise, use a_value (if it exists)
            push!(result_pairs, key => a_value)
        end
    end
    return (; result_pairs...)
end

"""
    recursive_merge(a, b, rest...)

Fold `recursive_merge` left to right over any number of named tuples, so later arguments win per
leaf exactly as they do for two. A boundary config contributing several independent groups — an open
edge contributes one per part of the state — merges them in one call rather than nesting pairs.
"""
recursive_merge(a::NamedTuple, b::NamedTuple, rest::NamedTuple...) =
    foldl(recursive_merge, rest; init = recursive_merge(a, b))

end # module