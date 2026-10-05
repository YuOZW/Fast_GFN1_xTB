param([ValidateSet('Release','Debug')][string]$BuildType = 'Release', [switch]$Benchmark)
$ErrorActionPreference = 'Stop'
if ($Benchmark -and $BuildType -ne 'Release') { throw 'Benchmarks require Release builds' }
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$arguments = "--exe `"$buildDir/fast-gfn1-xtb.exe`" --output `"$buildDir/halogen-cli-test`""
if ($Benchmark) { $arguments += " --reference `"$repoRoot/build-reference-xtb-6.7.1-windows-ifx-release/xtb.exe`"" }
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/test_runtime.py`" $arguments"
if ($LASTEXITCODE -ne 0) { throw 'Halogen Hessian CLI regression failed' }
