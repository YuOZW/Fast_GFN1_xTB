param(
    [ValidateSet('Release','Debug')][string]$BuildType = 'Release',
    [int]$Parallel = 4
)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
& (Join-Path $repoRoot 'tests/gfn1_fast_2_2_6/run_windows_build.ps1') `
    -BuildType $BuildType -Parallel $Parallel -BuildName 'build-gfn1-fast-3.0.0-windows-ifx'
if ($LASTEXITCODE -ne 0) { throw 'Native build regression failed' }
$buildDir = Join-Path $repoRoot "build-gfn1-fast-3.0.0-windows-ifx-$($BuildType.ToLowerInvariant())"
$mathDir = Join-Path $buildDir 'foe-math-test'
New-Item -ItemType Directory -Force $mathDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
& $PY (Join-Path $PSScriptRoot 'test_source_guards.py')
if ($LASTEXITCODE -ne 0) { throw 'FOE architecture guards failed' }
Push-Location $mathDir
try {
    $runtimeFlags = '/libs:dll /threads'
    if ($BuildType -eq 'Debug') { $runtimeFlags += ' /dbglibs' }
    $command = "call `"$setvars`" intel64 >nul && ifx /nologo $runtimeFlags /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/gfn1_fermi_operator.f90`" " +
        "`"$PSScriptRoot/test_fermi_operator.f90`" `"$buildDir/xtb.lib`" /Qmkl:sequential " +
        '/exe:test_fermi_operator.exe && test_fermi_operator.exe'
    & cmd /d /c $command
    if ($LASTEXITCODE -ne 0) { throw 'Strict Fermi operator regression failed' }
} finally { Pop-Location }
$command = "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/test_runtime.py`" " +
    "--exe `"$buildDir/fast-gfn1-xtb.exe`" --output `"$buildDir/foe-regression`""
if ($BuildType -eq 'Debug') { $command += ' --small-only' }
& cmd /d /c $command
if ($LASTEXITCODE -ne 0) { throw 'FOE runtime regression failed' }
