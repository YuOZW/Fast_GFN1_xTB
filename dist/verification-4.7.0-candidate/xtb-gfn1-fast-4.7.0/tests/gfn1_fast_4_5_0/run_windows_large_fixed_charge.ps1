param([ValidateSet('Release','Debug')][string]$BuildType = 'Release')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'large-solvent-response-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$flags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $flags += ' /dbglibs' }
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& set XTBPATH=$repoRoot&& " +
        "ifx /nologo $flags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/probe_large_fixed_charge.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        "/Qoption,link,/STACK:67108864 /exe:probe_large_fixed_charge.exe && probe_large_fixed_charge.exe `"$repoRoot/assets/inputs/xyz/taxol.xyz`""
    & cmd /d /c $command *> large-fixed-charge-probe.log
    if ($LASTEXITCODE -ne 0) { Get-Content large-fixed-charge-probe.log -Tail 40; throw 'Large fixed-charge probe failed' }
    Get-Content large-fixed-charge-probe.log -Tail 20
} finally { Pop-Location }
