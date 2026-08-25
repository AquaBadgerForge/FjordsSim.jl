using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FjordSim

config = oslofjorden()
prepare_bathymetry(config)
# download_forcing(config)
# prepare_forcing(config)
# add_rivers(config)
# download_boundaries(config)
# prepare_boundaries(config)
# download_atmosphere(config)
# prepare_atmosphere(config)

simulation = build_simulation(config)
