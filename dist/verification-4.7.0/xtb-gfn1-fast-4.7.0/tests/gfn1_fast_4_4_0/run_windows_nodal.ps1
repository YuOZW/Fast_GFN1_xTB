param([ValidateSet('Release','Debug')][string]$BuildType = 'Release')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'nodal-gradient-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$flags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $flags += ' /dbglibs' }
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& " +
        "set XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT=0&& " +
        "ifx /nologo $flags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/probe_nodal_gradient.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        '/Qoption,link,/STACK:67108864 /exe:test_nodal_gradient.exe && test_nodal_gradient.exe'
    & cmd /d /c $command *> nodal-gradient-test.log
    Get-Content nodal-gradient-test.log
    if ($LASTEXITCODE -ne 0) { throw 'Strict nodal Gradient/Hessian test failed' }
    $legacy = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& " +
        'set XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT=1&& test_nodal_gradient.exe legacy'
    & cmd /d /c $legacy *> legacy-nodal-test.log
    Get-Content legacy-nodal-test.log
    if ($LASTEXITCODE -ne 0) { throw 'Legacy diagnostic defect was not reproduced' }
} finally { Pop-Location }
