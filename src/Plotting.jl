module Plotting

export plot_bathymetry, plot_forcing, plot_atmosphere

using CairoMakie
using NCDatasets
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interior
using Oceananigans.Grids: x_domain, y_domain

using ..Configs:
    AbstractBathymetryConfig,
    AbstractForcingConfig,
    AbstractAtmosphereConfig,
    forcing_path,
    atmosphere_path,
    plot_path

# Logical figure pixels per grid cell, and fixed margin for title/labels/colorbar. Sized from
# the grid so every cell is rendered as a distinct, unambiguous block rather than being blurred
# together by too few output pixels per cell.
const BATHYMETRY_PLOT_PIXELS_PER_CELL = 6
const BATHYMETRY_PLOT_MARGIN = (300, 150)
const FORCING_PLOT_PIXELS_PER_CELL = 3
const FORCING_PLOT_MARGIN = (200, 150)
const ATMOSPHERE_PLOT_PIXELS_PER_CELL = 6
const ATMOSPHERE_PLOT_MARGIN = (200, 200)
const ATMOSPHERE_PLOT_COLUMNS = 4

# Surface fields worth eyeballing after preparation, with the colormap suiting each.
const FORCING_PLOT_PANELS = (
    (name = "T", label = "Temperature (ᵒC)", colormap = :thermal),
    (name = "S", label = "Salinity (g kg⁻¹)", colormap = :haline),
    (name = "u", label = "u (m s⁻¹)", colormap = :balance),
    (name = "v", label = "v (m s⁻¹)", colormap = :balance),
    (name = "T_lambda", label = "T relaxation rate (s⁻¹)", colormap = :viridis),
)

# Every variable of a prepared atmosphere file, with the colormap suiting each.
const ATMOSPHERE_PLOT_PANELS = (
    (name = "air_temperature_2m", label = "Air temperature (K)", colormap = :thermal),
    (name = "specific_humidity_2m", label = "Specific humidity (kg kg⁻¹)", colormap = :viridis),
    (name = "air_pressure_at_sea_level", label = "Sea level pressure (Pa)", colormap = :viridis),
    (name = "precipitation", label = "Precipitation (kg m⁻² s⁻¹)", colormap = :dense),
    (name = "u_wind_10m", label = "u wind (m s⁻¹)", colormap = :balance),
    (name = "v_wind_10m", label = "v wind (m s⁻¹)", colormap = :balance),
    (name = "swrad", label = "Downwelling shortwave (W m⁻²)", colormap = :solar),
    (name = "lwrad", label = "Downwelling longwave (W m⁻²)", colormap = :solar),
)

"""
    default_figure_size(grid, pixels_per_cell, margin; columns = 1)

Figure size in logical pixels for a `columns`-panel plot of `grid`, one block per grid cell.
"""
function default_figure_size(grid, pixels_per_cell, margin; columns = 1)
    Nx, Ny, _ = size(grid)
    return (
        columns * (Nx * pixels_per_cell + margin[1]),
        Ny * pixels_per_cell + margin[2],
    )
end

"""
    plot_axes(grid, Nx, Ny)

Longitude and latitude cell coordinates spanning the domain of `grid`, for data of shape
`(Nx, Ny)`. The shape is passed in rather than taken from `grid` because a staggered forcing
variable carries one extra face row.
"""
function plot_axes(grid, Nx, Ny)
    longitude = collect(range(x_domain(grid)[1], x_domain(grid)[2], length = Nx))
    latitude = collect(range(y_domain(grid)[1], y_domain(grid)[2], length = Ny))
    return longitude, latitude
end

prepare_plot_file(config) = let file = plot_path(config)
    isdir(dirname(file)) || mkpath(dirname(file))
    file
end

"""
    dataset_label(config)

The dataset a config names, for plot titles: its type name without the `Config` suffix, so
`NorKystConfig` reads as "NorKyst".
"""
dataset_label(config) = replace(string(nameof(typeof(config))), r"Config$" => "")

"""
    plot_bathymetry(grid, bottom_height, config::AbstractBathymetryConfig; title, figure_size)

Write a diagnostic bathymetry plot to `plot_path(config)`: bottom height with land blanked out,
depth contours every 25 m and the coastline at `h = 0`. Returns the plot path.
"""
function plot_bathymetry(
    grid,
    bottom_height,
    config::AbstractBathymetryConfig;
    title = "Bathymetry",
    figure_size = default_figure_size(grid, BATHYMETRY_PLOT_PIXELS_PER_CELL, BATHYMETRY_PLOT_MARGIN),
)
    plot_file = prepare_plot_file(config)

    cpu_bottom_height = on_architecture(CPU(), bottom_height)
    bathymetry_data = Array(interior(cpu_bottom_height, :, :, 1))

    Nx, Ny, _ = size(grid)
    longitude, latitude = plot_axes(grid, Nx, Ny)

    figure = Figure(size = figure_size)
    axis = Axis(figure[1, 1]; xlabel = "Longitude", ylabel = "Latitude", title)

    plot = heatmap!(axis, longitude, latitude, bathymetry_data; colormap = :deep, colorrange = extrema(bathymetry_data))
    land_mask = ifelse.(bathymetry_data .>= 0, 1.0f0, NaN32)
    heatmap!(
        axis,
        longitude,
        latitude,
        land_mask;
        colormap = [:ivory, :ivory],
        colorrange = (0, 1),
        nan_color = RGBAf(0, 0, 0, 0),
    )
    Colorbar(figure[1, 2], plot; label = "Bottom height (m)")
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = -collect(25:25:300),
        color = (:white, 0.35),
        linewidth = 1,
    )
    contour!(
        axis,
        longitude,
        latitude,
        bathymetry_data;
        levels = [0.0],
        color = :black,
        linewidth = 4,
    )
    save(plot_file, figure)

    return plot_file
end

"""
    surface_slice(ds, name)

The top vertical level of `name` at the first time step, as a `Float32` array with missing
values turned into `NaN` so Makie leaves land blank.
"""
function surface_slice(ds, name)
    Nz = ds.dim["Nz"]
    return Float32.(coalesce.(ds[name][:, :, Nz, 1], NaN32))
end

"""
    plot_forcing(grid, config::AbstractForcingConfig)

Write a diagnostic plot of the prepared forcing file at `forcing_path(config)` to
`plot_path(config)`: one panel per `FORCING_PLOT_PANELS` entry present in the file, at the
surface level of the first time step. Returns the plot path.
"""
function plot_forcing(grid, config::AbstractForcingConfig)
    forcing_file = forcing_path(config)
    plot_file = prepare_plot_file(config)

    NCDataset(forcing_file) do ds
        panels = [panel for panel in FORCING_PLOT_PANELS if haskey(ds, panel.name)]
        figure = Figure(
            size = default_figure_size(
                grid,
                FORCING_PLOT_PIXELS_PER_CELL,
                FORCING_PLOT_MARGIN;
                columns = length(panels),
            ),
        )
        date = ds["time"][1]

        for (column, panel) in enumerate(panels)
            data = surface_slice(ds, panel.name)
            Nx, Ny = size(data)
            longitude, latitude = plot_axes(grid, Nx, Ny)

            axis = Axis(
                figure[1, column];
                xlabel = "Longitude",
                ylabel = column == 1 ? "Latitude" : "",
                title = panel.name,
            )
            finite = filter(isfinite, data)
            colorrange = isempty(finite) || allequal(finite) ? (0, 1) : extrema(finite)
            plot = heatmap!(
                axis,
                longitude,
                latitude,
                data;
                colormap = panel.colormap,
                colorrange,
                nan_color = :ivory,
            )
            Colorbar(figure[2, column], plot; label = panel.label, vertical = false)
        end

        Label(figure[0, :], "Prepared $(dataset_label(config)) forcing, surface level at $date", fontsize = 20)
        save(plot_file, figure)
    end

    return plot_file
end

"""
    plot_atmosphere(config::AbstractAtmosphereConfig)

Write a diagnostic plot of the prepared atmosphere file at `atmosphere_path(config)` to
`plot_path(config)`: one panel per `ATMOSPHERE_PLOT_PANELS` entry present in the file, at the first
time step. Returns the plot path.

Unlike `plot_forcing` this takes no grid, because the prepared file carries its own `lon` and `lat`
axes — the atmosphere grid is independent of the ocean grid.
"""
function plot_atmosphere(config::AbstractAtmosphereConfig)
    atmosphere_file = atmosphere_path(config)
    plot_file = prepare_plot_file(config)

    NCDataset(atmosphere_file) do ds
        panels = [panel for panel in ATMOSPHERE_PLOT_PANELS if haskey(ds, panel.name)]
        isempty(panels) && error("$atmosphere_file holds none of the expected atmosphere variables")

        longitude = Array(ds["lon"][:])
        latitude = Array(ds["lat"][:])
        columns = min(length(panels), ATMOSPHERE_PLOT_COLUMNS)
        rows = cld(length(panels), columns)
        figure = Figure(
            size = (
                columns * (length(longitude) * ATMOSPHERE_PLOT_PIXELS_PER_CELL + ATMOSPHERE_PLOT_MARGIN[1]),
                rows * (length(latitude) * ATMOSPHERE_PLOT_PIXELS_PER_CELL + ATMOSPHERE_PLOT_MARGIN[2]),
            ),
        )
        date = ds["time"][1]

        for (position, panel) in enumerate(panels)
            # Two figure rows per panel row: the heatmap and its horizontal colorbar.
            row = 2 * cld(position, columns) - 1
            column = mod1(position, columns)
            data = Float32.(coalesce.(ds[panel.name][:, :, 1], NaN32))

            axis = Axis(
                figure[row, column];
                xlabel = "Longitude",
                ylabel = column == 1 ? "Latitude" : "",
                title = panel.name,
            )
            finite = filter(isfinite, data)
            colorrange = isempty(finite) || allequal(finite) ? (0, 1) : extrema(finite)
            plot = heatmap!(
                axis,
                longitude,
                latitude,
                data;
                colormap = panel.colormap,
                colorrange,
                nan_color = :ivory,
            )
            Colorbar(figure[row+1, column], plot; label = panel.label, vertical = false)
        end

        Label(figure[0, :], "Prepared $(dataset_label(config)) atmosphere at $date", fontsize = 20)
        save(plot_file, figure)
    end

    return plot_file
end

end  # module Plotting
