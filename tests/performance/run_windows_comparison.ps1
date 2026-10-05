param([int]$Repeats = 7)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$fast = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release/fast-gfn1-xtb.exe'
$reference = Join-Path $repoRoot 'build-reference-xtb-6.7.1-windows-ifx-release/xtb.exe'
$output = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release/stock-comparison'
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/compare_stages.py`" --exe `"$fast`" --reference `"$reference`" --output `"$output`" --repeats $Repeats"
if ($LASTEXITCODE -ne 0) { throw 'Stock/stage comparison failed' }
