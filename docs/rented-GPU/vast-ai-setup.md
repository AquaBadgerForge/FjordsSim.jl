# Run FjordsSim on a Vast.ai GPU instance

## What is a Vast.ai instance and why use it?

Vast.ai is a marketplace for renting GPU compute from many providers. An "instance" is a remote machine (VM or container)
with one or more NVIDIA GPUs, CPU, RAM, and disk, that you access via SSH. You pay per-hour and can attach storage volumes
to persist input/output data between runs.

Why use Vast.ai for FjordsSim:
- Cost-effective access to high-end GPUs (e.g., A100/H100) without owning hardware.
- Flexibility to pick GPU type, VRAM size, CPU/RAM, and attach storage volumes sized for your datasets.
- Fast spin-up: start a prebuilt image like `fisa47/fjordssim` (Julia 1.11 + CUDA 12.2 deps) and run immediately.
- Reproducible environments via Docker images; easy to share setup across collaborators.
- Good for one-off long runs or frequent experiments because of hourly billing.

This guide sets up a fresh Vast.ai GPU rental to run the Isafjardardjup simulation with GPU acceleration.

The steps assume a Linux Ubuntu base (22.04 recommended). You can either:
- Use the prebuilt Docker image `fisa47/fjordssim` (Julia 1.11 preinstalled, CUDA 12.2 deps included), or
- Use an NVIDIA CUDA base image and install Julia yourself (1.10.x recommended by this repo).

---

## 0) Requirements for an instance

- GPU: NVIDIA RTX 3060 or better (compute capability >= 7.0).
- VRAM: Minimum 8 GB (16+ GB recommended for larger domains).
- Driver/runtime: CUDA 12.x support (12.2 preferred; matches project runtime pin by default).
- Runtime: Docker image with CUDA 12.x userland, or rely on CUDA.jl’s binary runtime.
- Julia: 1.10.x recommended (repo compat), 1.11 works in practice; see notes below.
- Disk: Minimum 50 GB (100–200 GB recommended). Prefer attaching a storage volume on Vast.ai.

---

## 1) Rent and start an instance on Vast.ai

1. On vast.ai, click "Rent GPU" and choose an instance with:
  - GPU: NVIDIA 3060+
  - VRAM: >= 8 GB (>= 16 GB recommended)
  - Disk: Attach a storage volume (>= 50 GB; 100+ GB recommended) for inputs and results
  - Host driver: supports CUDA 12.x
  - Suggested Docker image options:
    - Recommended: `fisa47/fjordssim` (preinstalled Julia 1.11 and CUDA 12.2 dependencies)
    - Alternative: `nvidia/cuda:12.2.0-devel-ubuntu22.04` (or `runtime`) + manual Julia install

>You can now work in the shell like on any remote server — run commands, edit files, and launch Julia from the terminal. **Note**: Vast.ai instances often run the container/VM as root, so your home directory will be `/root` and your effective user is root. Be cautious with commands that assume a non-root user and with file permissions on attached volumes.

Practical tips:
- Expect warnings from tools such as Jupyter Notebook/Lab about running as root; if you start Jupyter add `--allow-root` (e.g., `jupyter lab --allow-root`).
- If you prefer a non-root workflow, create a user and adjust ownership of your data/results directories (`sudo useradd ...`; `sudo chown -R user:group /path`).


---


## 2) Install Julia 1.10.x (not needed for `fisa47/fjordssim` image)

If you used the `fisa47/fjordssim` image, Julia 1.11 is already preinstalled. You can usually use it as-is.  
If you used an NVIDIA CUDA base image, install Julia yourself

```bash
# Install Julia (1.10.x)
JVER=1.10.5
cd /opt
sudo mkdir -p julia && cd julia
sudo wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-$JVER-linux-x86_64.tar.gz
sudo tar -xf julia-$JVER-linux-x86_64.tar.gz
sudo ln -s /opt/julia/julia-$JVER/bin/julia /usr/local/bin/julia

# Verify
julia -e 'println(VERSION)'
```
---

## 3) Clone the project repos

```bash
mkdir -p ~/workspace && cd ~/workspace
# If you already have the repo in this path, skip clone
git clone https://github.com/AquaBadgerForge/FjordsSim.jl.git
# Optional: clone ClimaOcean.jl to watch for recent deps
git clone https://github.com/CliMA/ClimaOcean.jl.git
```

---

## 4) Instantiate Julia environment and GPU toolchain

Enter the project and instantiate packages.

```bash
cd ~/workspace/FjordsSim.jl
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

GPU toolchain sanity checks:

```bash
julia --project -e '
    using CUDA;
    CUDA.versioninfo();
    @show CUDA.functional();
'
```

**Notes:**
- This project pins CUDA.jl (compat `4, 5`) and sets the CUDA runtime in `app/isafjardardjup/simulation.jl` via `CUDA.set_runtime_version!(v"12.2")`.
- If your host/driver supports a different runtime (e.g., 12.4), change that line or set `CUDA.set_runtime_version!(nothing)` to let CUDA.jl auto-pick.

---

## 5) Prepare input data and directories

This app expects bathymetry and writes results under your home directory.

- Input bathymetry (JLD2):
  - Path expected by `app/isafjardardjup/setup.jl`:
    `~/FjordsSim_data/isafjardardjup/Isf_topo299x320.jld2`

Create data and results directories:

```bash
mkdir -p ~/FjordsSim_data/isafjardardjup
mkdir -p ~/FjordsSim_results/isafjardardjup
```

Copy/upload your `Isf_topo299x320.jld2` into the data path. To verify it contains `depth` and `z_faces`:

```bash
julia -e '
using JLD2; using Printf
f = joinpath(homedir(),"FjordsSim_data","isafjardardjup","Isf_topo299x320.jld2")
@info "Loading" f
JLD2.@load f depth z_faces
@info "depth size" size(depth)
@info "z_faces length" length(z_faces)
'
```
Also `fjordssim-notebooks` can be used to check the input data.

### IMPORTANT: z_faces must be positive and strictly increasing

Oceananigans (recent versions) require `z_faces` to be positive and strictly increasing. If your file has negative depths (older convention).

The project also sanitizes and sorts `z_faces` in `src/grids.jl`, but your source bathymetry file should be correct to avoid surprises.

### Tip: Upload input data from cloud

On Vast.ai, it’s convenient to pull inputs from cloud services:
- Google Drive/Dropbox/OneDrive
- Alternatively, mount an additional Vast.ai storage volume with your data.

---

## 6) Run the simulation (Isafjardardjup)

```bash
cd ~/workspace/FjordsSim.jl
julia --project app/isafjardardjup/simulation.jl
```

What you should see:
- Logs about model compilation and initialization, some warnings and a warning showing an error with IntervalArithmetics are fine
- Progress lines like `Iter: N, time: X, Δt: Y, max|u|: (...), extrema(T): (...)
- NetCDF snapshots under `/root/FjordsSim_results/isafjardardjup/snapshots.nc`

The simulation script is set to:
- Spin-up: 10 days with conservative CFL and small max Δt
- Main run: additional 355 days with CFL 0.25 and larger max Δt

---

## 7) Troubleshooting

- CUDA runtime mismatch:
  - If you see errors about CUDA runtime selection, edit `CUDA.set_runtime_version!` in `app/isafjardardjup/simulation.jl` to match the host driver’s supported CUDA version.

- NaN/Inf during run:
  - This repo includes a NaN guard callback (`nan_guard`) that aborts early and prints the offending field and index.
  - Check bathymetry: ensure non-zero wet cells and reasonable depths.
  - Ensure `e` and `ϵ` are initialized (done by the code) and check time step wizard settings (CFL).

- Missing data / downloads fail:
  - ClimaOcean will download JRA55 and other datasets via DataDeps on first use. Ensure outbound internet is allowed.

- Permissions/Paths:
  - Results are written to `~/FjordsSim_results/isafjardardjup`

