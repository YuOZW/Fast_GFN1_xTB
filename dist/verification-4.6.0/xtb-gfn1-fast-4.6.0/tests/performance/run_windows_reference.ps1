param([int]$Parallel = 4)
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$referenceRoot = Join-Path $repoRoot '_reference/xtb-6.7.1'
$buildDir = Join-Path $repoRoot 'build-reference-xtb-6.7.1-windows-ifx-release'
$commit = '26b28010e805f7d1aeeef39813feb473e69cc4be'
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
if (!(Test-Path (Join-Path $referenceRoot 'CMakeLists.txt'))) {
    throw 'Fetch the official grimme-lab/xtb v6.7.1 into _reference/xtb-6.7.1 first.'
}
$actual = & git -c "safe.directory=$($referenceRoot.Replace('\','/'))" -C $referenceRoot rev-parse HEAD
if ($LASTEXITCODE -ne 0 -or $actual.Trim() -ne $commit) { throw 'Unexpected reference source commit' }

# Upstream 6.7.1 unconditionally fetches/builds its unit-test dependency.
# Exclude test targets for this executable-only performance configuration.
$mainCmakePath = Join-Path $referenceRoot 'CMakeLists.txt'
$mainCmakeSource = [IO.File]::ReadAllText($mainCmakePath)
$mainCmakeSource = $mainCmakeSource.Replace('if(NOT TARGET "test-drive::test-drive")',
    'if(WITH_TESTS AND NOT TARGET "test-drive::test-drive")')
if ($mainCmakeSource -notmatch 'if\(WITH_TESTS\)\s+add_subdirectory\("test"\)') {
    $mainCmakeSource = $mainCmakeSource.Replace('add_subdirectory("test")',
        "if(WITH_TESTS)`n  add_subdirectory(`"test`")`nendif()")
}
[IO.File]::WriteAllText($mainCmakePath,$mainCmakeSource)

# Windows compatibility only; no fast GFN1 algorithms/property removal.
$compilerPath = Join-Path $referenceRoot 'cmake/CMakeLists.txt'
$compilerSource = [IO.File]::ReadAllText($compilerPath)
$compilerSource = $compilerSource.Replace('-axAVX2 -r8 -traceback','/QaxAVX2 /real-size:64 /traceback')
$compilerSource = $compilerSource.Replace('-check all -fpe0','/check:all /fpe:0')
[IO.File]::WriteAllText($compilerPath,$compilerSource)
$environmentPath = Join-Path $referenceRoot 'src/type/environment.f90'
$environmentSource = [IO.File]::ReadAllText($environmentPath)
if (!$environmentSource.Contains("rdvar('USERPROFILE'")) {
    $environmentSource = $environmentSource.Replace("   call rdvar('HOME', self%home, err)",
        "   call rdvar('HOME', self%home, err)`n" +
        "   if (.not.allocated(self%home)) call rdvar('USERPROFILE', self%home, err)`n" +
        "   if (.not.allocated(self%home)) self%home = ''")
    [IO.File]::WriteAllText($environmentPath,$environmentSource)
}
# Stock's nested allocatable descriptors crash at OpenMP entry with this
# Windows ifx runtime, including one-thread runs. Suppress only this routine's
# OpenMP directives for the serial benchmark. Numerical displacements, SCC
# calculations and contractions are unchanged; other parallelism is untouched.
$calculatorPath = Join-Path $referenceRoot 'src/type/calculator.f90'
$calculatorSource = [IO.File]::ReadAllText($calculatorPath)
if (!$calculatorSource.Contains('Windows reference: serial numerical Hessian')) {
    $start = $calculatorSource.IndexOf('   !$omp parallel do if(self%threadsafe)')
    $finish = $calculatorSource.IndexOf('end subroutine hessian', $start)
    if ($start -lt 0 -or $finish -lt 0) { throw 'Unexpected stock numerical Hessian OpenMP layout' }
    $body = $calculatorSource.Substring($start,$finish-$start).Replace('!$omp','! reference-omp')
    $calculatorSource = $calculatorSource.Substring(0,$start) +
        "   ! Windows reference: serial numerical Hessian (ifx OpenMP descriptor workaround).`n" +
        $body + $calculatorSource.Substring($finish)
    [IO.File]::WriteAllText($calculatorPath,$calculatorSource)
}
$dependencyDir = Join-Path $referenceRoot 'subprojects/mctc-lib'
if (!(Test-Path (Join-Path $dependencyDir 'CMakeLists.txt'))) {
    New-Item -ItemType Directory -Force $dependencyDir | Out-Null
    Get-ChildItem -LiteralPath (Join-Path $repoRoot 'subprojects/mctc-lib') -Force |
        Where-Object { $_.Name -ne '.git' } | Copy-Item -Destination $dependencyDir -Recurse -Force
}
$previousEnvironment = @{}
try {
    $lines = & cmd /d /c "call `"$setvars`" intel64 >nul && set"
    if ($LASTEXITCODE -ne 0) { throw 'oneAPI initialization failed' }
    foreach ($line in $lines) {
        if ($line -match '^([^=]+)=(.*)$') {
            $name=$Matches[1]
            if ($previousEnvironment.ContainsKey($name)) { continue }
            $previousEnvironment[$name]=[Environment]::GetEnvironmentVariable($name,'Process')
            [Environment]::SetEnvironmentVariable($name,$Matches[2],'Process')
        }
    }
    New-Item -ItemType Directory -Force $buildDir | Out-Null
    $log = Join-Path $buildDir 'build.log'
    & cmake -S $referenceRoot -B $buildDir -G Ninja -DCMAKE_BUILD_TYPE=Release `
        -DCMAKE_C_COMPILER=cl -DCMAKE_Fortran_COMPILER=ifx -DMCTCLIB_FIND_METHOD=subproject `
        -DWITH_TBLITE=OFF -DWITH_CPCMX=OFF -DWITH_TESTS=OFF -DWITH_JSON=OFF -DWITH_OBJECT=OFF `
        -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF `
        '-DCMAKE_Fortran_FLAGS_RELEASE=/O3 /DNDEBUG' `
        '-DCMAKE_EXE_LINKER_FLAGS=/Qoption,link,/STACK:67108864' *> $log
    if ($LASTEXITCODE -ne 0) { Get-Content $log -Tail 40; throw 'Reference configure failed' }
    & cmake --build $buildDir --target xtb-exe --parallel $Parallel *>> $log
    if ($LASTEXITCODE -ne 0) { Get-Content $log -Tail 40; throw 'Reference build failed' }
    Write-Output "PASS native original xTB 6.7.1 build ($commit): $buildDir/xtb.exe"
} finally {
    foreach ($name in $previousEnvironment.Keys) {
        [Environment]::SetEnvironmentVariable($name,$previousEnvironment[$name],'Process')
    }
}
