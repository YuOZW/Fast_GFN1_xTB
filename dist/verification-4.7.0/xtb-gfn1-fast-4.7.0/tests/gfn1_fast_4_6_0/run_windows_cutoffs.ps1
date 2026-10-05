param([ValidateSet('Release','Debug')][string]$BuildType = 'Release', [switch]$Diagnose)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'cutoff-response-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$flags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $flags += ' /dbglibs' }
$mode = ''
if ($Diagnose) { $mode = ' diagnose' }
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& " +
        "set XTBPATH=$repoRoot&& ifx /nologo $flags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/probe_cutoff_response.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        "/Qoption,link,/STACK:67108864 /exe:probe_cutoff_response.exe && probe_cutoff_response.exe$mode"
    & cmd /d /c $command *> cutoff-response-test.log
    if ($LASTEXITCODE -ne 0) { Get-Content cutoff-response-test.log -Tail 30; throw 'Cutoff response validation failed' }
    Get-Content cutoff-response-test.log -Tail 15
} finally { Pop-Location }
