param([int]$Parallel = 4)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
& (Join-Path $repoRoot 'tests/gfn1_fast_2_2_6/run_windows_build.ps1') `
    -Parallel $Parallel -BuildName 'build-gfn1-fast-2.3.2-windows-ifx'
if ($LASTEXITCODE -ne 0) { throw 'Native build regression failed' }
$buildDir = Join-Path $repoRoot 'build-gfn1-fast-2.3.2-windows-ifx-release'
$mathDir = Join-Path $buildDir 'math-test'
New-Item -ItemType Directory -Force $mathDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
Push-Location $mathDir
try {
    # Compile the new module itself with runtime checking, not just the test
    # caller. Link the existing native wrappers from the Release archive.
    $command = "call `"$setvars`" intel64 >nul && ifx /nologo /libs:dll /threads /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/gfn1_subspace.f90`" " +
        "`"$PSScriptRoot/test_subspace.f90`" `"$buildDir/xtb.lib`" /Qmkl:sequential " +
        '/exe:test_subspace.exe && test_subspace.exe'
    & cmd /d /c $command
    if ($LASTEXITCODE -ne 0) { throw 'Strict subspace math regression failed' }
} finally { Pop-Location }
$command = "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/test_runtime.py`" " +
    "--exe `"$buildDir/fast-gfn1-xtb.exe`" --output `"$buildDir/reuse-regression`""
& cmd /d /c $command
if ($LASTEXITCODE -ne 0) { throw 'Subspace runtime regression failed' }
