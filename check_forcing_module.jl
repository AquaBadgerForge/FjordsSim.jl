using FjordsSim

println("="^70)
println("Testing FjordsSim Forcing Separate Module")
println("="^70)

println("\nChecking if module loads...")
println("✓ FjordsSim loaded successfully")

println("\nChecking new functions:")
functions_to_check = [
    :forcing_from_separate_files,
    :forcing_from_boundaries_only,
    :forcing_from_rivers_only,
    :create_boundary_mask,
    :OpenBoundaryForcing,
    :RiverForcing,
    :BoundaryZoneMask,
]

all_found = true
for func in functions_to_check
    found = isdefined(FjordsSim, func)
    status = found ? "✓" : "✗"
    println("  $status $func")
    if !found
        all_found = false
    end
end

println("\n" * "="^70)
if all_found
    println("SUCCESS: All new functions are available!")
    println("="^70)
    println("\nNext steps:")
    println("1. Create forcing files: python examples/create_forcing_files.py")
    println("2. Run tests: include(\"test/test_forcing_separate.jl\")")
    println("3. Try examples: include(\"examples/forcing_separate_example.jl\")")
    println("4. Read docs: docs/forcing_separate_usage.md")
else
    println("ERROR: Some functions are missing!")
    println("="^70)
end
