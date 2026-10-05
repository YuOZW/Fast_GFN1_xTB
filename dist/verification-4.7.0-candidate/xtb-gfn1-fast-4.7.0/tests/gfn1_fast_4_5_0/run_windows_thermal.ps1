param([switch]$Reference)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$name = 'build-gfn1-fast-current-windows-ifx-release'
if ($Reference) { $name = 'build-reference-xtb-6.7.1-windows-ifx-release' }
$buildDir = Join-Path $repoRoot $name
$testDir = Join-Path $buildDir 'thermal-solvent-diagnostic'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=1&& set MKL_NUM_THREADS=1&& set XTBPATH=$repoRoot&& " +
        "ifx /nologo /libs:dll /threads /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/probe_thermal_energy.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        '/Qoption,link,/STACK:67108864 /exe:probe_thermal_energy.exe && probe_thermal_energy.exe'
    & cmd /d /c $command *> thermal-solvent-diagnostic.log
    if ($LASTEXITCODE -ne 0) { Get-Content thermal-solvent-diagnostic.log -Tail 40; throw 'Thermal solvent diagnostic failed' }
    Get-Content thermal-solvent-diagnostic.log -Tail 11
} finally { Pop-Location }
