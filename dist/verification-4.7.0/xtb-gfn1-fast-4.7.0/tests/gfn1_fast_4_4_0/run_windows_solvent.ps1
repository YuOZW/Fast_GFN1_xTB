param([ValidateSet('Release','Debug')][string]$BuildType = 'Release')
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'solvent-gradient-diagnostic'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$flags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $flags += ' /dbglibs' }
Push-Location $testDir
try {
    $setup = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& "
    $command = $setup + "ifx /nologo $flags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/probe_solvent_gradient.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        '/Qoption,link,/STACK:67108864 /exe:probe_solvent_gradient.exe && probe_solvent_gradient.exe'
    & cmd /d /c $command *> screened-boundary.log
    if ($LASTEXITCODE -ne 0) { Get-Content screened-boundary.log -Tail 40; throw 'Surface boundary probe failed' }
    Get-Content screened-boundary.log -Tail 3
    & cmd /d /c ($setup + 'probe_solvent_gradient.exe stable') *> smooth-nodal-gbsa.log
    if ($LASTEXITCODE -ne 0) { Get-Content smooth-nodal-gbsa.log -Tail 40; throw 'Strict nodal GBSA test failed' }
    Get-Content smooth-nodal-gbsa.log -Tail 3
} finally { Pop-Location }
