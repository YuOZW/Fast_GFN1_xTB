param([ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$OutputName = 'large-hessian-benchmark-4.4.0',
      [string]$Executable = '')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$buildDir = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release'
if (!$Executable) { $Executable = Join-Path $buildDir 'fast-gfn1-xtb.exe' }
$exePath = (Resolve-Path -LiteralPath $Executable).Path
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$arguments = "--exe `"$exePath`" --reference `"$repoRoot/build-reference-xtb-6.7.1-windows-ifx-release/xtb.exe`" --output `"$buildDir/$OutputName`""
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/compare_large_hessian.py`" $arguments"
if ($LASTEXITCODE -ne 0) { throw 'Large molecular Hessian comparison failed' }
