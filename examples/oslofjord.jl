# End-to-end Oslofjord simulation.
#
# Everything this run needs is in the setup: the grid, the prepared bathymetry, forcing and
# atmosphere files, and the `SimulationConfig` naming the closure, the run length and the
# outputs. To change any of it, edit `src/Setups/oslofjorden.jl` rather than this file.
#
# The preparation steps have to have run first:
#
#   julia --project -m FjordSim prepare_bathymetry  --config oslofjorden
#   julia --project -m FjordSim download_forcing    --config oslofjorden
#   julia --project -m FjordSim prepare_forcing     --config oslofjorden
#   julia --project -m FjordSim add_rivers          --config oslofjorden
#   julia --project -m FjordSim download_atmosphere --config oslofjorden
#   julia --project -m FjordSim prepare_atmosphere  --config oslofjorden
#
# This script is the same thing as `julia --project -m FjordSim run_simulation --config
# oslofjorden`. To step through the assembly instead of running it, build the simulation without
# starting it:
#
#   simulation = build_simulation(oslofjorden())
#   run!(simulation)

using FjordSim
using CUDA  # it should be here to make GPU() not throw an error

run_simulation(oslofjorden())
