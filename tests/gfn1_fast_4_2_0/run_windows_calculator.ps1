param([ValidateSet('Release','Debug')][string]$BuildType = 'Release',
      [ValidateRange(1,28)][int]$Threads = 1)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "build-gfn1-fast-current-windows-ifx-$($BuildType.ToLowerInvariant())"
$testDir = Join-Path $buildDir 'hessian-calculator-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$runtimeFlags = '/libs:dll /threads'
if ($BuildType -eq 'Debug') { $runtimeFlags += ' /dbglibs' }
Push-Location $testDir
try {
    $command = "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=$Threads&& set MKL_NUM_THREADS=1&& " +
        'set XTB_GFN1_FAST_ENABLE_ANALYTIC_HESSIAN=1&& set XTB_GFN1_FAST_PROFILE=1&& ' +
        "ifx /nologo $runtimeFlags /Od /check:all /fpe:0 /Qopenmp /I:`"$buildDir/include`" " +
        "`"$PSScriptRoot/test_calculator.f90`" `"$buildDir/xtb.lib`" " +
        "`"$buildDir/subprojects/mctc-lib/mctc-lib.lib`" /Qmkl:sequential " +
        '/Qoption,link,/STACK:67108864 /exe:test_calculator.exe && test_calculator.exe'
    & cmd /d /c $command *> calculator-test.log
    if ($LASTEXITCODE -ne 0) { Get-Content calculator-test.log -Tail 40; throw 'Calculator Hessian API test failed' }
    $used = @(Select-String -Path calculator-test.log -Pattern '4\.\d+\.\d+ analytic Hessian: used=([TF])')
    if ($used.Count -ne 2 -or $used[0].Matches[0].Groups[1].Value -ne 'T' -or $used[1].Matches[0].Groups[1].Value -ne 'F') {
        throw 'Actual analytic/polarizability fallback dispatch not verified'
    }
    Get-Content calculator-test.log
} finally { Pop-Location }
