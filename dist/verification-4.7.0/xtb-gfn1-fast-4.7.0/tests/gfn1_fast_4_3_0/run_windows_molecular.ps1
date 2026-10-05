param([ValidateSet('Release','Debug')][string]$BuildType = 'Release')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'molecular-halogen-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$runtimeFlags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $runtimeFlags += ' /dbglibs' }
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& " +
        "ifx /nologo $runtimeFlags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$repoRoot/src/xtb/halogen.f90`" `"$repoRoot/src/xtb/hessian_response.f90`" " +
        "`"$PSScriptRoot/test_molecular_halogen.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        '/Qoption,link,/STACK:67108864 /exe:test_molecular_halogen.exe && test_molecular_halogen.exe'
    & cmd /d /c $command *> molecular-halogen-test.log
    Get-Content molecular-halogen-test.log
    if ($LASTEXITCODE -ne 0) { throw 'Strict molecular halogen Hessian test failed' }
} finally { Pop-Location }
