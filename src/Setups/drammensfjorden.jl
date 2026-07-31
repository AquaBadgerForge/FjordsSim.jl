"""
    drammensfjorden()

The Drammensfjord setup: Geonorge Sjøkart Dybdedata bathymetry and NorKyst-800m forcing on a
150 x 200 x 11 grid covering 10.20-10.45°E, 59.58-59.75°N.

It names neither rivers nor an atmosphere, so both default to `nothing` and `add_rivers`,
`download_atmosphere` and `prepare_atmosphere` are no-ops for it.
"""
function drammensfjorden()
    data_root = joinpath(homedir(), "FjordSim_data", "drammensfjorden")

    return FjordConfig(
        grid_config = EvenGrid(
            size      = (150, 200, 11),
            halo      = (7, 7, 7),
            longitude = (10.20, 10.45),
            latitude  = (59.58, 59.75),
            z_faces   = [-100.0, -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0],
        ),
        bathymetry_config = DybdedataConfig(
            data_root             = data_root,
            output_file           = "bathymetry.nc",
            plot_file             = "bathymetry.png",
            raw_resolution_factor = 2,
            padding_cells         = 2,
            include_contours      = false,
            contour_stride        = 10,
            interpolation_passes  = 1,
            major_basins          = 1,
            geonorge_cache        = false,
            regrid_cache          = false,
        ),
        forcing_config = NorKystConfig(
            data_root            = data_root,
            output_directory     = "norkyst",
            output_file          = "forcing.nc",
            plot_file            = "forcing.png",
            relaxation_edge      = :south,
            relaxation_cells     = 10,
            relaxation_timescale = 86400.0,
            architecture         = :auto,
            parameters           = ["temperature", "salinity", "u_eastward", "v_northward"],
            years                = [2020],
        ),
    )
end
