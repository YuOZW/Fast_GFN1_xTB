param([ValidateSet('Release','Debug')][string]$BuildType = 'Release', [switch]$Benchmark)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$arguments = "--exe `"$buildDir/fast-gfn1-xtb.exe`" --output `"$buildDir/solvent-cli-test`""
if ($Benchmark) {
    if ($BuildType -ne 'Release') { throw 'Benchmarks require Release builds' }
    $arguments += " --reference `"$repoRoot/build-reference-xtb-6.7.1-windows-ifx-release/xtb.exe`""
}
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/test_runtime.py`" $arguments"
if ($LASTEXITCODE -ne 0) { throw 'Solvent Hessian CLI regression failed' }
