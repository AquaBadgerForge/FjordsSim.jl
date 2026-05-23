using FjordSim

_data_root = expanduser("~/FjordSim_data/drammensfjorden")

(
    grid = GridConfig(
        size      = (150, 200, 11),
        halo      = (7, 7, 7),
        longitude = (10.20, 10.45),
        latitude  = (59.58, 59.75),
        z_faces   = [-100.0, -75.0, -50.0, -25.0, -15.0, -10.0, -7.5, -5.0, -3.0, -2.0, -1.0, 0.0],
    ),
    bathymetry = BathymetryConfig(
        output_path           = joinpath(_data_root, "bathymetry.nc"),
        plot_path             = joinpath(_data_root, "bathymetry.png"),
        geodatabase_path      = joinpath(expanduser("~/FjordSim_data"), "Basisdata_0000_Norge_25833_Dybdedata_FGDB.gdb"),
        raw_resolution_factor = 2,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        geonorge_cache        = true,
        regrid_cache          = false,
    ),
    norkyst = NorKystConfig(
        output_dir  = joinpath(_data_root, "norkyst"),
        catalog_url = "https://thredds.met.no/thredds/catalog/fou-hi/norkyst800m/catalog.xml",
        opendap_url = "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m/",
        parameters  = ["temperature", "salinity", "u_eastward", "v_northward"],
        years       = [2020],
    ),
)
