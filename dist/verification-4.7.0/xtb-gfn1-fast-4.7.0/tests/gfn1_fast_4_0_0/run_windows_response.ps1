param([ValidateSet('Release','Debug')][string]$LibraryBuildType = 'Release',
      [ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$BuildName = 'build-gfn1-fast-current-windows-ifx')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "$BuildName-$($LibraryBuildType.ToLowerInvariant())"
$mathDir = Join-Path $buildDir 'response-math-test'
if (!(Test-Path (Join-Path $buildDir 'xtb.lib'))) { throw 'Run the requested native build first' }
New-Item -ItemType Directory -Force $mathDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
Push-Location $mathDir
try {
    $runtimeFlags = '/libs:dll /threads'
    if ($LibraryBuildType -eq 'Debug') { $runtimeFlags += ' /dbglibs' }
    # Compile response itself with runtime checks, independently of archive.
    $command = "call `"$setvars`" intel64 >nul && ifx /nologo $runtimeFlags /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/xtb/response.f90`" " +
        "`"$PSScriptRoot/test_response.f90`" `"$buildDir/xtb.lib`" /Qmkl:sequential " +
        '/exe:test_response.exe && test_response.exe'
    & cmd /d /c $command
    if ($LASTEXITCODE -ne 0) { throw 'Strict electronic-response test failed' }
    $command = "call `"$setvars`" intel64 >nul && ifx /nologo $runtimeFlags /Od /check:all /fpe:0 " +
        "/I:`"$buildDir/include`" `"$repoRoot/src/xtb/coulomb.f90`" " +
        "`"$PSScriptRoot/test_coulomb_response.f90`" `"$buildDir/xtb.lib`" /Qmkl:sequential " +
        '/exe:test_coulomb_response.exe && test_coulomb_response.exe'
    & cmd /d /c $command
    if ($LASTEXITCODE -ne 0) { throw 'Strict Coulomb-response test failed' }
} finally { Pop-Location }
