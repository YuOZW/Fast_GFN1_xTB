param([ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$OutputName = 'parallel-large-probe-4.7.0-final')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$buildDir = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release'
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/probe_large_parallel.py`" --exe `"$buildDir/xtb.exe`" --output `"$buildDir/$OutputName`""
if ($LASTEXITCODE -ne 0) { throw 'Large OpenMP response diagnostic failed' }
