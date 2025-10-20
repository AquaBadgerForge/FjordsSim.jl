using Oceananigans.Fields: interior
using Oceananigans.OutputReaders: FieldTimeSeries, OnDisk
using Oceananigans.Utils: prettytime, pretty_filesize
using NCDatasets: Dataset
using JLD2: @save
using Printf: @sprintf
using Logging: @warn

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))

    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    @info msg

    wall_time[] = time_ns()
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
        z = []
    end
    return z
end

function netcdf_to_jld2(netcdf_file::String, jld2_file::String)
    ds = Dataset(netcdf_file, "r")
    data_dict = Dict()
    for varname in keys(ds)
        data_dict[varname] = convert(Array, ds[varname])
        print(size(convert(Array, ds[varname])))
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

"""
Check ocean state for NaN/Inf and abort early with a helpful message.

Checks u, v, w and all tracers on the ocean model inside a coupled simulation.
Attach as a Simulation callback, e.g.
    coupled_simulation.callbacks[:nan_guard] = Callback(check_for_nans, IterationInterval(1))
"""
function check_for_nans(sim)
    ocean_sim = sim.model.ocean
    ocean_model = ocean_sim.model

    # Check velocities
    for (nm, fld) in ((:u, ocean_model.velocities.u), (:v, ocean_model.velocities.v), (:w, ocean_model.velocities.w))
        A = interior(fld)
        if !all(isfinite, A)
            lin = findfirst(x -> !isfinite(x), A)
            val = isnothing(lin) ? NaN : A[lin]
            I = isnothing(lin) ? (missing, missing, missing) : ind2sub(size(A), lin)
            @warn "Detected non-finite value in field" field = nm value = val index = I iteration = iteration(sim) time = ocean_model.clock.time
            error("NaN/Inf detected in $(nm); aborting simulation to prevent corrupt outputs.")
        end
    end

    # Check tracers
    for nm in propertynames(ocean_model.tracers)
        fld = getproperty(ocean_model.tracers, nm)
        A = interior(fld)
        if !all(isfinite, A)
            lin = findfirst(x -> !isfinite(x), A)
            val = isnothing(lin) ? NaN : A[lin]
            I = isnothing(lin) ? (missing, missing, missing) : ind2sub(size(A), lin)
            @warn "Detected non-finite value in tracer" tracer = nm value = val index = I iteration = iteration(sim) time = ocean_model.clock.time
            error("NaN/Inf detected in tracer $(nm); aborting simulation to prevent corrupt outputs.")
        end
    end

    return nothing
end

