using Oceananigans
using JLD2
using Oceananigans.Units
using CairoMakie
using FjordsSim:
    plot_1d_phys,
    extract_z_faces,
    record_vertical_tracer,
    record_surface_speed,
    record_horizontal_tracer,
    plot_ztime,
    plot_bottom_tracer,
    record_bottom_tracer,
    record_vertical_diff

folder = joinpath(homedir(), "FjordsSim_results", "varna")
filename = joinpath(folder, "snapshots")
T = FieldTimeSeries("$filename.jld2", "T")
S = FieldTimeSeries("$filename.jld2", "S")
u = FieldTimeSeries("$filename.jld2", "u")
v = FieldTimeSeries("$filename.jld2", "v")
O₂ = FieldTimeSeries("$filename.jld2", "O₂")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
PHY = FieldTimeSeries("$filename.jld2", "P")
HET = FieldTimeSeries("$filename.jld2", "HET")
DOM = FieldTimeSeries("$filename.jld2", "DOM")
POM = FieldTimeSeries("$filename.jld2", "POM")
C = FieldTimeSeries("$filename.jld2", "C")
times = T.times

# filename10 = joinpath(folder, "varna_snapshotsDOM10")
# filename20 = joinpath(folder, "varna_snapshots2024-11-19")

# O₂10 =  FieldTimeSeries("$filename10.jld2", "O₂")[:,:,:,1:2:end]
# O₂20 =  FieldTimeSeries("$filename20.jld2", "O₂")


# O₂diff = O₂10 .- O₂20

grid = jldopen("$filename.jld2")["grid"]
Nz = grid["underlying_grid"]["Nz"]


# bottom_z evaluation
bottom_z = ones(Int, size(O₂, 1), size(O₂, 2))
for i = 1:size(O₂, 1)
    for j = 1:size(O₂, 2)
        for k = size(O₂, 3):-1:1  # Loop backwards to find the latest non-zero
            if O₂[i, j, k, 1] == 0
                bottom_z[i, j] = k
                if k != Nz
                    bottom_z[i, j] = k + 1
                end
                break
            end
        end
    end
end


println(keys(grid["underlying_grid"]))
println(grid["underlying_grid"]["Δyᶠᶜᵃ"])

# stupid, but I cannot find a right way with znodes
# znodes(grid["underlying_grid"], with_halos=false)
z = grid["underlying_grid"]["zᵃᵃᶜ"][8:Nz+7]

# z = extract_z_faces(grid)

plot_ztime(PHY, HET, POM, DOM, NUT, O₂, T, S, 84, 14, times, z, folder)

# HORIZONTAL
# plot_1d_phys(T, S, z, times, folder)

# plot_bottom_oxygen(O₂, bottom_z, 1000 , folder)

record_surface_speed(u, v, Nz, times, folder)

# record_bottom_oxygen(O₂, times, bottom_z, folder)

record_horizontal_tracer(
    C,
    times,
    folder,
    "Contsurf",
    "Contaminant (% of max. concentration)",
    colorrange = (0, 100),
    colormap = :matter,
    iz = Nz,
)

record_horizontal_tracer(
    T,
    times,
    folder,
    "Tsurf",
    "Temperature (°C)",
    colorrange = (5, 40),
    colormap = Reverse(:RdYlBu),
    iz = Nz,
)

record_horizontal_tracer(
    S,
    times,
    folder,
    "Ssurf",
    "Salinity (PSU)",
    iz = Nz,
    colorrange = (0, 17),
    colormap = :viridis,
)

record_horizontal_tracer(
    O₂,
    times,
    folder,
    "O2bottom",
    "Dissolved oxygen (μM)",
    colorrange = (100, 350),
    colormap = :turbo,
    iz = 3,
)

record_horizontal_tracer(
    NUT,
    times,
    folder,
    "NUTsurf",
    "Nutrients (μM N)",
    colorrange = (0, 20),
    colormap = Reverse(:cherry),
    iz = Nz,
)

record_horizontal_tracer(
    PHY,
    times,
    folder,
    "PHYsurf",
    "Phytoplankton (μM N)",
    colorrange = (0, 2),
    colormap = Reverse(:cubehelix),
    iz = Nz,
)

record_horizontal_tracer(
    DOM,
    times,
    folder,
    "DOMsurf10",
    "DOM (μM N)",
    colorrange = (0, 50),
    colormap = Reverse(:CMRmap),
    iz = Nz,
)


# VERTICAL
record_vertical_tracer(
    T,
    z,
    18,
    times,
    folder,
    "Tprofile",
    "Temperature (°C)",
    colorrange = (5, 40),
    colormap = Reverse(:RdYlBu),
)

record_vertical_tracer(
    u,
    z,
    18,
    times,
    folder,
    "Uprofile",
    "u velocity component (ms⁻¹)",
    colorrange = (-0.5, 0.5),
    colormap = :deep,
)

record_vertical_tracer(S, z, 18, times, folder, "Sprofile", "Salinity (PSU)", colorrange = (0, 17), colormap = :viridis)

record_vertical_tracer(
    PHY,
    z,
    18,
    times,
    folder,
    "PHYprofile",
    "Phytoplankton (μM N)",
    colorrange = (0, 2),
    colormap = Reverse(:cubehelix),
)

record_vertical_tracer(
    DOM,
    z,
    18,
    times,
    folder,
    "DOMprofile",
    "DOM (μM N)",
    colorrange = (0, 20),
    colormap = Reverse(:CMRmap),
)

record_vertical_tracer(
    POM,
    z,
    18,
    times,
    folder,
    "POMprofile",
    "POM (μM N)",
    colorrange = (0, 100),
    colormap = Reverse(:CMRmap),
)

record_vertical_tracer(
    NUT,
    z,
    18,
    times,
    folder,
    "NUTprofile",
    "Nutrients (μM N)",
    colorrange = (0, 10),
    colormap = Reverse(:cherry),
)

record_vertical_tracer(
    O₂,
    z,
    18,
    times,
    folder,
    "O2profile",
    "Dissolved Oxygen (μM)",
    colorrange = (100, 350),
    colormap = :turbo,
)

# record_vertical_diff(
#     O₂diff, z, 18, times, folder, "O2profile_diff", "O₂ (DOM=10) - O₂ (DOM=20)",
#     colorrange=(-200, 200), colormap=:balance,
#         )
