# FjordSim.jl

A framework for ocean simulations built on top of [Oceananigans](https://github.com/CliMA/Oceananigans.jl) and [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl).

FjordSim's main contribution is a streamlined way to set up regional simulations.
A simulation is assembled from the following components:

1. **Grid**
   Domain bounds, horizontal size, and vertical faces.

2. **Bathymetry file**
   Contains the domain coordinates (longitude, latitude, depth) along with the bathymetric data.
   It can be generated from the public [Geonorge](https://www.geonorge.no/) Sjøkart Dybdedata dataset, see below.

3. **Forcing file**
   Includes information about sinks and sources (e.g., rivers), boundary conditions, and custom forcings.
   It can also be used to load initial conditions.

4. **Atmospheric data**
   Supports JRA55 from [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl) and [NORA3](https://thredds.met.no/thredds/projects/nora3.html).

## Installation

There are several options:

- Install from github:
1. Clone the git repository `git clone https://github.com/NIVANorge/FjordSim.jl.git`.
2. Move to the downloaded folder `cd FjordSim.jl`.
3. Run Julia REPL and activate the FjordSim environment `julia --project`.
4. Enter the Pkg REPL by pressing `]` from Julia REPL.
5. Type `instantiate` to 'resolve' a `Manifest.toml` from a `Project.toml` to install and precompile dependency packages.

- Add the latest FjordSim to your Julia project: `add https://github.com/NIVANorge/FjordSim.jl.git`.

- Add from the Julia registry: `add FjordSim`.

## Setups (work in progress)

Each fjord is described by a `FjordConfig` in a Julia file under `configs/`
(`oslofjorden.jl`, `drammensfjorden.jl`), which groups a grid, a bathymetry, and a forcing
configuration:

```julia
FjordConfig(
    grid_config       = EvenGrid(size = (150, 200, 11), longitude = (10.20, 10.45), ...),
    bathymetry_config = DybdedataConfig(data_root = _data_root, output_file = "bathymetry.nc", ...),
    forcing_config    = NorKystConfig(data_root = _data_root, years = [2020]),
)
```

By default, input data goes to `$HOME/FjordSim_data/<fjord>/` and results to
`$HOME/FjordSim_results/<fjord>/`. The config fields naming individual files are relative to the
setup's `data_root`; setting one to an absolute path relocates just that file, which is how a
single copy of the Geonorge database is shared between fjords.

To add a fjord, copy an existing config and adjust it.

To add a new kind of grid, bathymetry source, or forcing dataset, define a struct subtyping
`AbstractGridConfig`, `AbstractBathymetryConfig`, or `AbstractForcingConfig` and overload that
supertype's hooks — `FjordConfig`, the generic pipelines and the scripts stay unchanged. Each
pipeline is one generic function plus a handful of dispatch points:

| Pipeline | Generic entry point | Hooks a new source implements |
|---|---|---|
| Grid | — | `LatitudeLongitudeGrid(arch, config)` |
| Bathymetry | `prepare_bathymetry(target_grid, config)` | `bathymetry_dataset(target_grid, config)`; optionally `regrid_options(config)` |
| Forcing | `prepare_forcing(target_grid, config)`, `download_forcing(config)` | `forcing_time_steps`, `forcing_source_grid`, `forcing_variable_names`; `download_forcing(target_grid, config)` if it downloads |

Path resolution (`bathymetry_path`, `forcing_path`, `forcing_directory`, `plot_path`) and the
diagnostic plots (`plot_bathymetry`, `plot_forcing`) are defined on the supertypes, so a new
source inherits them. `src/Bathymetry/Geonorge.jl` and `src/Forcing/NorKyst.jl` are the built-in
adapters and the templates to copy; each supertype's docstring in `src/Configs.jl` lists the
fields and hooks it expects.

## Prepare input data

Every script takes the setup to prepare via `--config`, which is required:

```bash
# Bathymetry: downloads the Geonorge Sjøkart FileGDB on first use (~2.3 GB), regrids it onto the
# configured grid, and writes the bathymetry NetCDF plus a diagnostic plot
julia --project scripts/bathymetry_prepare.jl --config configs/oslofjorden.jl

# Forcing: downloads and subsets the configured dataset (NorKyst-800m) over the setup's region
# and years
julia --project scripts/forcing_download.jl --config configs/oslofjorden.jl

# Forcing: regrids the download onto the simulation grid, writing the forcing NetCDF and a plot
julia --project scripts/forcing_prepare.jl --config configs/oslofjorden.jl
```

The scripts are thin CLI wrappers: the work is `prepare_bathymetry`, `download_forcing` and
`prepare_forcing`, each a generic function on the matching config supertype.

## Run an example Oslofjord simulation

1. Download the [grid, forcing, atmospheric forcing](https://www.dropbox.com/scl/fo/gc3yc155b5eohi7998wgh/AGN2Yt3HyQ0LlZGImpcca6o?rlkey=x6okc3uxe2avud6sbxgd00l14&st=093llyqp&dl=0).
2. In `FjordSim.jl/examples/oslofjord.jl` it is possible to specify the location of the input data files.
By default, the files should be in `$HOME/FjordSim_data/oslofjord/` and `$HOME/FjordSim_data/JRA55/` or `$HOME/FjordSim_data/NORA3/`.
Also, it is possible to specify the results folder destination.
By default, the result will go to `$HOME/FjordSim_results/oslofjord/`.
3. Run `julia --project examples/oslofjord.jl`.
This will generate a netcdf results file.

The example does not go through `FjordConfig`: it builds the grid straight from a bathymetry
NetCDF and wires the components together by hand.

Further preparation scripts for the Oslofjord are available in the following repository:
[https://github.com/NIVANorge/oslofjord-sim](https://github.com/NIVANorge/oslofjord-sim)

![example_result](./artifacts/phytoplankton_multi.png)
