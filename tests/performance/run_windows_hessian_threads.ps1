param([string]$Executable = '', [ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$OutputName = 'hessian-thread-crossover-4.5.1')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$buildDir = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release'
if (!$Executable) { $Executable = Join-Path $buildDir 'xtb.exe' }
$exePath = (Resolve-Path -LiteralPath $Executable).Path
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$arguments = "--exe `"$exePath`" --output `"$buildDir/$OutputName`""
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/compare_hessian_threads.py`" $arguments"
if ($LASTEXITCODE -ne 0) { throw 'Hessian thread crossover measurements failed' }
