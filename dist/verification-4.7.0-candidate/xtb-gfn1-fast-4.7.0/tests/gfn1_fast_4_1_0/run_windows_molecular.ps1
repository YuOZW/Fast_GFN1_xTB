param(
    [ValidateSet('Release','Debug')][string]$LibraryBuildType = 'Release',
    [ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$BuildName = 'build-gfn1-fast-current-windows-ifx',
    [ValidateRange(1,28)][int]$Threads = 1
)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "$BuildName-$($LibraryBuildType.ToLowerInvariant())"
$mathDir = Join-Path $buildDir 'molecular-response-test'
if (!(Test-Path (Join-Path $buildDir 'xtb.lib'))) { throw 'Run the native build first' }
New-Item -ItemType Directory -Force $mathDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
Push-Location $mathDir
try {
    $runtimeFlags = '/libs:dll /threads'
    if ($LibraryBuildType -eq 'Debug') { $runtimeFlags += ' /dbglibs' }
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=$Threads&& set MKL_NUM_THREADS=1&& ifx /nologo $runtimeFlags /fpp /Qopenmp /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/xtb/repulsion.F90`" `"$repoRoot/src/disp/dftd3.f90`" " +
        "`"$repoRoot/src/xtb/response.f90`" " +
        "`"$repoRoot/src/xtb/hessian_response.f90`" `"$PSScriptRoot/test_molecular_response.f90`" " +
        "`"$buildDir/xtb.lib`" `"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" " +
        '/Qmkl:sequential /Qoption,link,/STACK:67108864 /exe:test_molecular_response.exe && test_molecular_response.exe'
    & cmd /d /c $command *> molecular-response-test.log
    Get-Content molecular-response-test.log
    if ($LASTEXITCODE -ne 0) { throw 'Strict actual-molecule response test failed' }
} finally { Pop-Location }
