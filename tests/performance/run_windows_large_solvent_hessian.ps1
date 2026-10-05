param([string]$Executable = '', [switch]$AuditOnly,
      [ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$OutputName = 'large-solvent-hessian-benchmark-4.5.1')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$buildDir = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release'
if (!$Executable) { $Executable = Join-Path $buildDir 'fast-gfn1-xtb.exe' }
$exePath = (Resolve-Path -LiteralPath $Executable).Path
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$arguments = "--exe `"$exePath`" --output `"$buildDir/$OutputName`""
if ($AuditOnly) { $arguments += ' --audit-only' }
else { $arguments += " --reference `"$repoRoot/build-reference-xtb-6.7.1-windows-ifx-release/xtb.exe`"" }
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/compare_large_solvent_hessian.py`" $arguments"
if ($LASTEXITCODE -ne 0) { throw 'Large solvent Hessian audit/comparison failed' }
