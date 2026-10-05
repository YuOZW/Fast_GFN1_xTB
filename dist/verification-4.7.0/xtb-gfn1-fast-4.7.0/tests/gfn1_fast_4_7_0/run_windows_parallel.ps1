param([ValidateSet('Release','Debug')][string]$BuildType = 'Release')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/test_parallel_runtime.py`" --exe `"$buildDir/xtb.exe`" --output `"$buildDir/parallel-hessian-test`""
if ($LASTEXITCODE -ne 0) { throw 'Actual analytic OpenMP regression failed' }
