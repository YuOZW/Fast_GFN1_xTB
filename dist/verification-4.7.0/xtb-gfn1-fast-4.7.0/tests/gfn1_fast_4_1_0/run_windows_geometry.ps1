param(
    [ValidateSet('Release','Debug')][string]$LibraryBuildType = 'Release',
    [ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$BuildName = 'build-gfn1-fast-current-windows-ifx'
)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "$BuildName-$($LibraryBuildType.ToLowerInvariant())"
$mathDir = Join-Path $buildDir 'geometry-math-test'
if (!(Test-Path (Join-Path $buildDir 'xtb.lib'))) { throw 'Run the native build first' }
New-Item -ItemType Directory -Force $mathDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
Push-Location $mathDir
try {
    $runtimeFlags = '/libs:dll /threads'
    if ($LibraryBuildType -eq 'Debug') { $runtimeFlags += ' /dbglibs' }
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& ifx /nologo $runtimeFlags /Qopenmp /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/intgrad.f90`" " +
        "`"$repoRoot/src/disp/coordinationnumber.f90`" `"$PSScriptRoot/test_geometry.f90`" " +
        "`"$buildDir/xtb.lib`" `"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" " +
        '/Qmkl:sequential /Qoption,link,/STACK:67108864 /exe:test_geometry.exe && test_geometry.exe'
    & cmd /d /c $command *> geometry-test.log
    Get-Content geometry-test.log
    if ($LASTEXITCODE -ne 0) { throw 'Strict geometry-Hessian test failed' }
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& ifx /nologo $runtimeFlags /Qopenmp /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/intgrad.f90`" " +
        "`"$repoRoot/src/disp/coordinationnumber.f90`" `"$repoRoot/src/xtb/hessian_integrals.f90`" " +
        "`"$PSScriptRoot/test_integral_response.f90`" " +
        "`"$buildDir/xtb.lib`" `"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" " +
        '/Qmkl:sequential /Qoption,link,/STACK:67108864 /exe:test_integral_response.exe && test_integral_response.exe'
    & cmd /d /c $command *> integral-response-test.log
    Get-Content integral-response-test.log
    if ($LASTEXITCODE -ne 0) { throw 'Strict molecular integral-response test failed' }
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& ifx /nologo $runtimeFlags /Qopenmp /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/coulomb/klopmanohno.f90`" " +
        "`"$PSScriptRoot/test_coulomb_geometry.f90`" " +
        "`"$buildDir/xtb.lib`" `"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" " +
        '/Qmkl:sequential /Qoption,link,/STACK:67108864 /exe:test_coulomb_geometry.exe && test_coulomb_geometry.exe'
    & cmd /d /c $command *> coulomb-geometry-test.log
    Get-Content coulomb-geometry-test.log
    if ($LASTEXITCODE -ne 0) { throw 'Strict Coulomb geometry-response test failed' }
} finally { Pop-Location }
