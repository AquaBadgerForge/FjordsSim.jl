using FjordSim

_data_root = expanduser("~/FjordSim_data/oslofjorden")

FjordConfig(
    grid_config = EvenGrid(
        size      = (240, 520, 18),
        halo      = (7, 7, 7),
        longitude = (10.2, 11.02),
        latitude  = (59.0, 59.93),
        z_faces   = [
            -450.0, -400.0, -350.0, -300.0, -250.0, -200.0, -150.0, -100.0,
            -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0,
        ],
    ),
    bathymetry_config = DybdedataConfig(
        data_root             = _data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        raw_resolution_factor = 2,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        geonorge_cache        = true,
        regrid_cache          = false,
    ),
    forcing_config = NorKystConfig(
        data_root = _data_root,
        years     = [2020],
    ),
)
