param([ValidateSet('Release','Debug')][string]$BuildType = 'Release', [switch]$Diagnose)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'screened-response-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$flags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $flags += ' /dbglibs' }
$mode = ''
if ($Diagnose) { $mode = ' diagnose' }
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& " +
        "set XTBPATH=$repoRoot&& set XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT=0&& " +
        "ifx /nologo $flags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/probe_screened_response.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        "/Qoption,link,/STACK:67108864 /exe:probe_screened_response.exe && probe_screened_response.exe$mode"
    & cmd /d /c $command *> screened-response-test.log
    if ($LASTEXITCODE -ne 0) { Get-Content screened-response-test.log -Tail 30; throw 'Screened response validation failed' }
    Get-Content screened-response-test.log -Tail 15
} finally { Pop-Location }
