# Copyright 2024 The FjordsSim Authors.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
using CUDA
CUDA.set_runtime_version!(v"12.2")

using Oceananigans.Units: second, seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Utils: TimeInterval, IterationInterval
using Oceananigans.Simulations: Callback, conjure_time_step_wizard!, run!
using Oceananigans.OutputWriters: NetCDFWriter
using Oceanostics
using FjordsSim: coupled_hydrostatic_simulation, progress
using Printf
using Logging
using Dates

# Set up real-time logging to file
log_file = open("simulation.log", "w")

# Custom logging function that writes to both console and file
function log_message(msg)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    full_msg = "[$timestamp] $msg"
    println(stdout, full_msg)
    println(log_file, full_msg)
    flush(stdout)
    flush(log_file)
end

log_message("Starting simulation setup")

include("setup.jl")

## Model Setup
sim_setup = setup_region_3d()

# Ensure the output directory exists
mkpath(sim_setup.results_dir)
log_message("Output directory created: $(sim_setup.results_dir)")

coupled_simulation = coupled_hydrostatic_simulation(sim_setup)
coupled_simulation.callbacks[:progress] = Callback(progress, TimeInterval(3hours))

## Set up output writers
ocean_sim = coupled_simulation.model.ocean
ocean_model = ocean_sim.model

prefix = joinpath(sim_setup.results_dir, "snapshots")
log_message("Output file prefix: $prefix")

ocean_sim.output_writers[:nc_writer] = NetCDFWriter(
    ocean_model, merge(ocean_model.tracers, ocean_model.velocities);
    schedule = TimeInterval(1hours),
    filename = "$prefix.nc",
    overwrite_existing = true,
)

## Spinning up the simulation
# We use an adaptive time step that maintains the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) equal to 0.1.
log_message("Starting spin-up phase (10 days)")
ocean_sim.stop_time = 10days
coupled_simulation.stop_time = 10days

conjure_time_step_wizard!(ocean_sim; cfl=0.1, max_Δt=1.5minutes, max_change=1.01)
log_message("Running spin-up simulation...")
run!(coupled_simulation)
log_message("Spin-up phase completed successfully")

## Running the simulation
# This time, we set the CFL in the time_step_wizard to be 0.25 as this is the maximum recommended CFL to be
# used in conjunction with Oceananigans' hydrostatic time-stepping algorithm ([two step Adams-Bashfort](https://en.wikipedia.org/wiki/Linear_multistep_method))
log_message("Starting main simulation phase (355 days)")
ocean_sim.stop_time = 355days
coupled_simulation.stop_time = 355days

conjure_time_step_wizard!(ocean_sim; cfl=0.25, max_Δt=10minutes, max_change=1.01)
log_message("Running main simulation...")
run!(coupled_simulation)
log_message("Simulation completed successfully at $(now())")
