# FjordsSim.jl Installation Notes

## Overview
This document provides detailed installation instructions and troubleshooting notes for FjordsSim.jl, including dependency management and compatibility requirements.

## System Requirements

### Julia Version
- **Supported versions**: Julia 1.10 or 1.11
- **Note**: Oceananigans (a core dependency) is officially tested on Julia 1.10, but works on Julia 1.11 with minor warnings

### Operating System
- Linux (tested)
- macOS (should work)
- Windows (should work, may require WSL for optimal performance)

## Installation Methods

### Method 1: Standard Installation (Registry Packages)

If all packages are available in the Julia registry:

```julia
using Pkg
Pkg.activate("path/to/FjordsSim.jl")
Pkg.instantiate()
```

### Method 2: Installing with ClimaOcean from GitHub

**Important**: ClimaOcean v0.8.7 may not be available in the Julia registry. If you encounter errors, install it directly from GitHub:

```julia
using Pkg
Pkg.activate("path/to/FjordsSim.jl")

# Install ClimaOcean from GitHub first
Pkg.add(PackageSpec(url="https://github.com/CliMA/ClimaOcean.jl.git", rev="v0.8.7"))

# Then resolve and instantiate other dependencies
Pkg.instantiate()
```

## Dependency Compatibility Notes

### Critical Package Versions

The following packages require specific version ranges for compatibility:

1. **ClimaOcean**: v0.8.x
   - Must be installed from GitHub if v0.8.7 is not in registry
   - URL: `https://github.com/CliMA/ClimaOcean.jl.git`

2. **GPUArrays**: v9, v10, or v11
   - ClimaOcean v0.8.7 requires v11+
   - Updated from v9-only to support newer versions

3. **OceanBioME**: v0.1.x or v0.10.x - v0.14.x
   - Must support both legacy (0.1) and modern (0.10+) versions
   - Required for compatibility with Adapt v4

4. **Oceanostics**: v0.1.x, v0.15.x, or v0.16.x
   - Must support v0.16.x for compatibility with Oceananigans v0.100+

5. **Oceananigans**: v0.96.x or v0.100.x
   - Core simulation engine
   - Note: Will show warnings on Julia 1.11 (expected behavior)

6. **Interpolations**: v0.14.x or v0.15.x
   - Allow v0.15 for modern dependency resolution

7. **CairoMakie**: v0.10.x - v0.13.x
   - Visualization package
   - Broader range for compatibility with newer plotting features

## Installation Order

When installing from scratch or resolving conflicts, follow this order:

1. **Update Project.toml compatibility constraints** (if needed)
2. **Install ClimaOcean from GitHub** (if v0.8.7 not in registry)
3. **Run Pkg.resolve()** to update manifest
4. **Run Pkg.instantiate()** to install remaining packages
5. **Run Pkg.precompile()** to precompile all packages

## Common Installation Issues

### Issue 1: ClimaOcean Version Mismatch

**Error**: `empty intersection between ClimaOcean@0.6.11 and project compatibility 0.8`

**Solution**:
```julia
using Pkg
Pkg.rm("ClimaOcean")
Pkg.add(PackageSpec(url="https://github.com/CliMA/ClimaOcean.jl.git", rev="v0.8.7"))
```

### Issue 2: GPUArrays Version Conflict

**Error**: `Unsatisfiable requirements detected for package GPUArrays`

**Solution**: Update `Project.toml` compat section:
```toml
GPUArrays = "9, 10, 11"
```

### Issue 3: OceanBioME/Oceanostics Compatibility

**Error**: Restrictions from Adapt or Oceananigans preventing installation

**Solution**: Update `Project.toml` compat section:
```toml
OceanBioME = "0.1, 0.10, 0.11, 0.12, 0.13, 0.14"
Oceanostics = "0.1, 0.15, 0.16"
```

### Issue 4: Interpolations Version Conflict

**Error**: `empty intersection between Interpolations@0.15.1 and project compatibility 0.14`

**Solution**: Update `Project.toml` compat section:
```toml
Interpolations = "0.14, 0.15"
```

### Issue 5: Precompilation Failures

**Error**: `✗ Oceananigans → OceananigansMakieExt` or `✗ ClimaOcean`

**Possible Causes**:
- Missing ClimaOcean in manifest
- Version conflicts
- Julia 1.11 compatibility issues (warnings only, not critical)

**Solution**:
1. Ensure ClimaOcean is properly installed (see Issue 1)
2. Run `Pkg.resolve()` to update manifest
3. Run `Pkg.precompile()` again

## Complete Installation Script

For a fresh installation with all compatibility fixes:

```julia
using Pkg

# Navigate to project
cd("path/to/FjordsSim.jl")
Pkg.activate(".")

# Install ClimaOcean from GitHub
println("Installing ClimaOcean v0.8.7 from GitHub...")
Pkg.add(PackageSpec(url="https://github.com/CliMA/ClimaOcean.jl.git", rev="v0.8.7"))

# Resolve dependencies
println("Resolving dependencies...")
Pkg.resolve()

# Install remaining packages
println("Installing packages...")
Pkg.instantiate()

# Precompile
println("Precompiling packages...")
Pkg.precompile()

println("Installation complete!")

# Verify installation
using ClimaOcean
println("ClimaOcean v", pkgversion(ClimaOcean), " loaded successfully!")
```

## Testing Installation

To verify your installation is working correctly:

```julia
using Pkg
Pkg.activate(".")

# Load core packages
using ClimaOcean
using Oceananigans
using OceanBioME
using Oceanostics

println("All core packages loaded successfully!")

# Check versions
println("ClimaOcean: v", pkgversion(ClimaOcean))
println("Oceananigans: v", pkgversion(Oceananigans))
println("OceanBioME: v", pkgversion(OceanBioME))
println("Oceanostics: v", pkgversion(Oceanostics))
```

## Expected Warnings

### Julia Version Warning
```
Warning: You are using Julia v1.11 or later!
Oceananigans is currently tested on Julia v1.10.
```
**Status**: This is expected and non-critical. Oceananigans works on Julia 1.11 despite this warning.

### Sediments Warning
```
WARNING: could not import TimeSteppers.store_tendencies! into Sediments
```
**Status**: This warning appears during OceanBioME precompilation and does not affect functionality.

## Troubleshooting Tips

1. **Always check `Pkg.status()` first** to see current package versions
2. **Use `Pkg.status("PackageName")` for specific packages** to check their installation status
3. **Check compatibility with `Pkg.status(mode=PKGMODE_MANIFEST)`** to see resolved versions
4. **Use `Pkg.status(; outdated=true)`** to see available updates
5. **Clear package cache if persistent issues**: Remove `~/.julia/packages/` and reinstall
6. **Check manifest**: Ensure `Manifest.toml` matches `Project.toml` constraints

## Updating Packages

To update packages while respecting compatibility constraints:

```julia
using Pkg
Pkg.update()
```

To update specific packages:

```julia
Pkg.update("PackageName")
```

## Support

If you encounter issues not covered here:

1. Check the [Oceananigans.jl issues](https://github.com/CliMA/Oceananigans.jl/issues)
2. Check the [ClimaOcean.jl issues](https://github.com/CliMA/ClimaOcean.jl/issues)
3. Open an issue in the FjordsSim.jl repository

## Last Updated

- **Date**: October 17, 2025
- **Julia Version**: 1.11.5
- **ClimaOcean Version**: 0.8.7 (from GitHub)
- **Oceananigans Version**: 0.100.6
