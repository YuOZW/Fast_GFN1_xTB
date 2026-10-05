$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$testDir = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release/resource-policy-test'
New-Item -ItemType Directory -Force $testDir | Out-Null
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
Push-Location $testDir
try {
    foreach ($openmp in @($true,$false)) {
        $flag = if ($openmp) { '/Qopenmp' } else { '' }
        $name = if ($openmp) { 'omp' } else { 'serial' }
        $compile = "call `"$setvars`" intel64 >nul && ifx /nologo /libs:dll /threads /Od /check:all /fpe:0 $flag `"$repoRoot/src/gfn1_fast_policy.f90`" `"$PSScriptRoot/test_resources.f90`" /exe:resources-$name.exe"
        & cmd /d /c $compile
        if ($LASTEXITCODE -ne 0) { throw 'Resource test compilation failed' }
        $cases = @(@{ expected=8; setting='' },
                   @{ expected=2; setting='set XTB_GFN1_FAST_HESSIAN_OMP_MAX_THREADS=2&& ' },
                   @{ expected=1; setting='set XTB_GFN1_FAST_DISABLE_HESSIAN_OPENMP=1&& ' },
                   @{ expected=1; setting='set XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO=999&& ' },
                   @{ expected=3; setting='set XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB=900&& ' },
                   @{ expected=1; setting='set XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB=1&& ' })
        foreach ($case in $cases) {
            $expected = if ($openmp) { $case.expected } else { 1 }
            $clean = 'set XTB_GFN1_FAST_HESSIAN_OMP_MAX_THREADS=8&& set XTB_GFN1_FAST_DISABLE_HESSIAN_OPENMP=0&& set XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO=64&& set XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB=1024&& '
            & cmd /d /c "call `"$setvars`" intel64 >nul && set OMP_NUM_THREADS=8&& set OMP_DYNAMIC=FALSE&& $clean$($case.setting)resources-$name.exe $expected"
            if ($LASTEXITCODE -ne 0) { throw 'Resource policy behavior failed' }
        }
    }
} finally { Pop-Location }
