param([ValidateSet('Release','Debug')][string]$BuildType = 'Release')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'solvent-response-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$flags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $flags += ' /dbglibs' }
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& set XTBPATH=$repoRoot&& " +
        'set XTB_GFN1_FAST_ENABLE_ANALYTIC_HESSIAN=1&& set XTB_GFN1_FAST_PROFILE=0&& ' +
        "ifx /nologo $flags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/test_molecular_solvent.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        '/Qoption,link,/STACK:67108864 /exe:test_molecular_solvent.exe && test_molecular_solvent.exe'
    & cmd /d /c $command *> molecular-solvent-response-test.log
    if ($LASTEXITCODE -ne 0) { Get-Content molecular-solvent-response-test.log -Tail 40; throw 'Molecular solvent response regression failed' }
    Get-Content molecular-solvent-response-test.log -Tail 10
} finally { Pop-Location }
