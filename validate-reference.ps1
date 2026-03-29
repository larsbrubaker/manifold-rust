param(
    [string]$Phase = "all",
    [switch]$SkipCppBuild,
    [switch]$SkipRustTests
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$cppRoot = Join-Path $repoRoot "cpp-reference/manifold"
$cppBuild = Join-Path $cppRoot "build"

if (-not (Test-Path $cppRoot)) {
    throw "Missing C++ reference at '$cppRoot'. Run: git submodule update --init --recursive"
}

$phaseMap = @{
    "all" = @{
        CppFilter = "*"
        RustFilters = @()
    }
    "phase3" = @{
        CppFilter = "Polygon*"
        RustFilters = @("polygon::")
    }
    "phase6" = @{
        CppFilter = "*Constructor*:*Warp*"
        RustFilters = @("constructors::")
    }
    "phase7" = @{
        CppFilter = "*Normal*:*Coplanar*:*Degenerate*"
        RustFilters = @("face_op::", "edge_op::")
    }
    "phase8" = @{
        CppFilter = "*Property*:*Volume*:*Area*:*Convex*"
        RustFilters = @("properties::")
    }
    "phase9" = @{
        CppFilter = "*Smooth*:*Tangent*"
        RustFilters = @("smoothing::", "svd::")
    }
    "phase11" = @{
        CppFilter = "*Boolean*"
        RustFilters = @("boolean3::", "boolean_result::")
    }
    "phase13" = @{
        CppFilter = "*Manifold*"
        RustFilters = @("manifold::")
    }
    "phase14" = @{
        CppFilter = "*CrossSection*"
        RustFilters = @("cross_section::")
    }
    "phase16" = @{
        CppFilter = "*SDF*:*LevelSet*"
        RustFilters = @("sdf::")
    }
    "phase17" = @{
        CppFilter = "*Hull*:*Minkowski*"
        RustFilters = @("quickhull::", "minkowski::")
    }
}

if (-not $phaseMap.ContainsKey($Phase)) {
    $valid = ($phaseMap.Keys | Sort-Object) -join ", "
    throw "Unknown phase '$Phase'. Valid values: $valid"
}

if (-not $SkipCppBuild) {
    cmake -S $cppRoot -B $cppBuild -DMANIFOLD_TEST=ON -DMANIFOLD_PAR=OFF -DMANIFOLD_DEBUG=OFF
    cmake --build $cppBuild --config Release
}

$cppFilter = $phaseMap[$Phase].CppFilter
$testExeCandidates = @(
    (Join-Path $cppBuild "bin/Release/manifold_test.exe"),
    (Join-Path $cppBuild "bin/manifold_test.exe"),
    (Join-Path $cppBuild "test/Release/manifold_test.exe"),
    (Join-Path $cppBuild "test/manifold_test.exe")
)
$testExe = $testExeCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1

if (-not $testExe) {
    throw "Could not find built C++ test executable under '$cppBuild'."
}

Write-Host "Running C++ reference tests for phase '$Phase'..."
& $testExe "--gtest_filter=$cppFilter"

if (-not $SkipRustTests) {
    $rustFilters = $phaseMap[$Phase].RustFilters
    Push-Location $repoRoot
    try {
        Write-Host "Running Rust tests for phase '$Phase'..."
        if ($rustFilters.Count -eq 0) {
            cargo test --lib
        } else {
            foreach ($rustFilter in $rustFilters) {
                cargo test --lib $rustFilter
            }
        }
    } finally {
        Pop-Location
    }
}

Write-Host "Validation finished for phase '$Phase'."
