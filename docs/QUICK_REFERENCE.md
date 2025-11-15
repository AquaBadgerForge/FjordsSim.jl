# Forcing Separate Module - Quick Reference Card

## 🎯 Purpose
Handle open boundary conditions and river forcing from separate NetCDF files.

## 📦 Installation
Already included in FjordsSim! Just use:
```julia
using FjordsSim
```

## 🚀 Quick Start

### 1. Create Forcing Files (Python)
```python
python examples/create_forcing_files.py
```

### 2. Use in Julia
```julia
using FjordsSim

setup = setup_region(
    forcing_callable = forcing_from_separate_files,
    forcing_args = (
        grid_ref = grid_ref,
        boundary_filepath = "path/to/boundary.nc",
        river_filepath = "path/to/river.nc",
        tracers = (:T, :S),
        buffer_width = 10,  # optional
    ),
    # ... other args ...
)
```

## 📊 File Formats

### Boundary File
```
Dimensions: Nx, Ny, Nz, time
Variables:
  T(time, Nz, Ny, Nx)         # boundary values
  T_lambda(time, Nz, Ny, Nx)  # relaxation rate (1/s)
  S(time, Nz, Ny, Nx)
  S_lambda(time, Nz, Ny, Nx)
```

### River File
```
Dimensions: Nx, Ny, Nz, time
Variables:
  T_flux(time, Nz, Ny, Nx)    # volume flux (m³/s)
  T(time, Nz, Ny, Nx)         # river water temperature
  S_flux(time, Nz, Ny, Nx)
  S(time, Nz, Ny, Nx)
```

## 🔧 Three Usage Modes

### Combined (Boundaries + Rivers)
```julia
forcing_callable = forcing_from_separate_files
forcing_args = (
    grid_ref = grid_ref,
    boundary_filepath = "boundary.nc",
    river_filepath = "river.nc",
    tracers = (:T, :S),
    buffer_width = 10,
)
```

### Boundaries Only
```julia
forcing_callable = forcing_from_boundaries_only
forcing_args = (
    grid_ref = grid_ref,
    boundary_filepath = "boundary.nc",
    tracers = (:T, :S),
    buffer_width = 10,
)
```

### Rivers Only
```julia
forcing_callable = forcing_from_rivers_only
forcing_args = (
    grid_ref = grid_ref,
    river_filepath = "river.nc",
    tracers = (:T, :S),
)
```

## ⚙️ Key Parameters

### Buffer Width
- Controls boundary transition zone
- Typical: 5-20 grid cells
- Smaller = sharper transition
- Larger = smoother transition

```julia
buffer_width = 5   # Sharp (narrow transition)
buffer_width = 10  # Default (balanced)
buffer_width = 20  # Smooth (wide transition)
```

### Relaxation Timescale (λ)
- Controls nudging strength at boundaries
- Units: 1/seconds
- In your NetCDF file as `T_lambda`, `S_lambda`, etc.

```julia
λ = 1e-5  # Weak   (~1 day relaxation)
λ = 1e-4  # Normal (~3 hours) ← typical
λ = 1e-3  # Strong (~10 minutes)
```

## 🎨 Boundary Buffer Zone Visualization

```
Domain Edge → Buffer Zone → Interior
    |            |            |
   100%      75%-25%-0%       0%
    |            |            |
  [0][1][2][3][4][5][6]...[Nx]
   ^^^        ^^^         ^^^
  Full     Gradual        No
forcing    decay       forcing
```

## 📝 Python File Creation Template

```python
from examples.create_forcing_files import ForcingFileGenerator

gen = ForcingFileGenerator(Nx=100, Ny=120, Nz=30)

# Boundary file
gen.create_boundary_file(
    'boundary.nc',
    {
        'T': {'west': 10.0, 'east': 11.0, 'lambda': 1e-4},
        'S': {'west': 34.0, 'east': 35.0, 'lambda': 1e-4},
    },
    buffer_width=10
)

# River file
gen.create_river_file(
    'river.nc',
    river_locations=[(50, 100, 29)],  # (i, j, k)
    river_properties={
        'discharge': [100.0],  # m³/s
        'T': [8.0],           # °C
        'S': [0.5],           # psu
    }
)
```

## 🧪 Testing

```julia
# Load and run tests
include("test/test_forcing_separate.jl")

# Try examples
include("examples/forcing_separate_example.jl")
simulation = run_example()
```

## 📖 Full Documentation
- **Complete guide**: `docs/forcing_separate_usage.md`
- **Examples**: `examples/forcing_separate_example.jl`
- **Python tools**: `examples/create_forcing_files.py`
- **Summary**: `docs/FORCING_MODULE_SUMMARY.md`

## 🐛 Common Issues

### DimensionMismatch
**Problem**: Grid sizes don't match
**Solution**: Check Nx, Ny, Nz in NetCDF == grid dimensions

### No forcing effect
**Problem**: Values all zero or fill value
**Solution**: Check λ > 0 at boundaries, flux > 0 at rivers

### Instability
**Problem**: Model crashes near boundaries
**Solution**: 
- Decrease λ (weaker nudging)
- Increase buffer_width
- Check for unrealistic boundary values

## 📞 Getting Help

1. Read `docs/forcing_separate_usage.md` (comprehensive)
2. Check examples in `examples/`
3. Review troubleshooting in docs
4. Verify with test suite

## ✅ Checklist for First Use

- [ ] Grid dimensions correct (Nx, Ny, Nz)?
- [ ] Boundary values reasonable?
- [ ] Lambda values set (typically 1e-4)?
- [ ] Buffer width chosen (5-20)?
- [ ] River locations in water cells?
- [ ] River discharge positive?
- [ ] Tracer names match model?
- [ ] File paths correct?
- [ ] Time coordinates correct?

## 🎓 Key Concepts

**Open Boundary Conditions**: Nudge model toward specified values at domain edges

**Buffer Zone**: Smooth transition region preventing sharp discontinuities

**Relaxation**: `tendency = -λ * (model - boundary)`, λ controls strength

**River Forcing**: `tendency = (flux/volume) * (river - model)`

## 💡 Tips

1. **Start conservative**: Use moderate λ (~1e-4) and standard buffer (10)
2. **Test incrementally**: Try boundaries first, add rivers second
3. **Visualize**: Plot forcing fields before running simulation
4. **Monitor**: Check temperature/salinity at boundaries during run
5. **Calibrate**: Adjust λ and buffer_width based on results

## 🔗 Related Functions

From existing `forcing.jl`:
- `forcing_from_file()` - Original combined forcing
- `ForcingFromFile` - Original forcing structure
- `NetCDFBackend` - Time series backend (reused)

## 📊 Performance

- **Memory**: 2 time steps in memory (configurable)
- **GPU**: Fully compatible via Adapt.jl
- **I/O**: Loads data on-demand
- **Overhead**: Minimal, same as original forcing

---

**Version**: 1.0  
**Created**: November 2025  
**Module**: `src/forcing_separate.jl`
