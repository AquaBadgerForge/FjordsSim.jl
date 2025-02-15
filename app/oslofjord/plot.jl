using Oceananigans
using JLD2
using Oceananigans.Units
using CairoMakie: Axis, Figure, Colorbar, Observable, Reverse, record, heatmap!, @lift

function record_variable(
    variable,
    var_name,
    Nz,
    times,
    folder,
    figsize;
    colorrange = (0, 0.5),
    colormap = :deep,
    framerate = 12,
)
    Nt = length(times)
    iter = Observable(Nt)

    si = @lift begin
        s = variable[$iter]
        s = interior(s, :, :, Nz)
        s[s.==0] .= NaN
        s
    end

    fig = Figure(size = figsize)
    title = @lift "$(var_name) at " * prettytime(times[$iter])
    ax = Axis(
        fig[1, 1];
        title = title,
        xlabel = "Grid points, eastward direction",
        ylabel = "Grid points, northward direction",
    )
    hm = heatmap!(ax, si, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(var_name) (ms⁻¹)")

    record(fig, joinpath(folder, "$(var_name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

folder = joinpath(homedir(), "FjordsSim_results", "oslofjord")
filename = joinpath(folder, "snapshots")
T = FieldTimeSeries("$filename.jld2", "T")
S = FieldTimeSeries("$filename.jld2", "S")
u = FieldTimeSeries("$filename.jld2", "u")
v = FieldTimeSeries("$filename.jld2", "v")

grid = jldopen("$filename.jld2")["grid"]
Nz = grid["underlying_grid"]["Nz"]

record_variable(T, "temperature surface", Nz, T.times, folder, (300, 700); colorrange = (-1, 20))
record_variable(S, "salinity surface", Nz, S.times, folder, (300, 700); colorrange = (-1, 40))
record_variable(u, "u velocity surface", Nz, u.times, folder, (300, 700); colorrange = (-5, 5))
record_variable(v, "v velocity surface", Nz, v.times, folder, (300, 700); colorrange = (-5, 5))

record_variable(T, "temperature bottom", 1, T.times, folder, (300, 700); colorrange = (-1, 20))
record_variable(S, "salinity bottom", 1, S.times, folder, (300, 700); colorrange = (-1, 40))
record_variable(u, "u velocity bottom", 1, u.times, folder, (300, 700); colorrange = (-5, 5))
record_variable(v, "v velocity bottom", 1, v.times, folder, (300, 700); colorrange = (-5, 5))
