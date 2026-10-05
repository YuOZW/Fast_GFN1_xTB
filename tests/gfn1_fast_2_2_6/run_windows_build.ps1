param(
    [ValidateSet('Release', 'Debug')][string]$BuildType = 'Release',
    [int]$Parallel = 4,
    [ValidatePattern('^[a-zA-Z0-9.-]+$')][string]$BuildName = 'build-gfn1-fast-current-windows-ifx',
    [string]$OneApiSetvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $repoRoot "$BuildName-$($BuildType.ToLowerInvariant())"
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$dependencyCommit = '77f65c6f2cf6330d05d0757ca173da097096780e'

function Invoke-Checked {
    param([string]$Program, [string[]]$Arguments)
    & $Program @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Program failed with exit code $LASTEXITCODE" }
}

if (!(Test-Path $PY)) { throw "Required Python not found: $PY" }
if (!(Test-Path $OneApiSetvars)) { throw "Intel oneAPI environment script not found: $OneApiSetvars" }
if (!(Test-Path (Join-Path $repoRoot 'subprojects/mctc-lib/CMakeLists.txt'))) {
    throw "Prepare subprojects/mctc-lib at commit $dependencyCommit first (see WINDOWS_TEST_REPORT.md)."
}
$dependencyDir = Join-Path $repoRoot 'subprojects/mctc-lib'
if (Test-Path -LiteralPath (Join-Path $dependencyDir '.git')) {
    $actualCommit = & git -c "safe.directory=$($dependencyDir.Replace('\', '/'))" -C $dependencyDir rev-parse HEAD
    if ($LASTEXITCODE -ne 0 -or $actualCommit.Trim() -ne $dependencyCommit) {
        throw "mctc-lib must be pinned at $dependencyCommit"
    }
} else {
    # A source ZIP intentionally contains no Git database. Do not let git
    # walk up and mistake an enclosing project's HEAD for the dependency.
    $packageManifest = Join-Path $repoRoot 'SOURCE_SHA256.json'
    $packageVerifier = Join-Path $repoRoot 'tests/release/verify_source_package.py'
    if (!(Test-Path -LiteralPath $packageManifest) -or !(Test-Path -LiteralPath $packageVerifier)) {
        throw "Use the pinned mctc-lib Git checkout or a hash-verified source package ($dependencyCommit)."
    }
    Invoke-Checked $PY @($packageVerifier, '--root', $repoRoot, '--dependencies-only')
}

# Import into this PowerShell process only; do not change persistent user settings.
$previousEnvironment = @{}
try {
    $environmentLines = & cmd /d /c "call `"$OneApiSetvars`" intel64 >nul && set"
    if ($LASTEXITCODE -ne 0) { throw 'Intel oneAPI environment setup failed' }
    foreach ($line in $environmentLines) {
        if ($line -match '^([^=]+)=(.*)$') {
            $name = $Matches[1]
            # Some launchers supply both PATH and Path. cmd's updated PATH
            # comes first; do not overwrite it with the inherited duplicate.
            if ($previousEnvironment.ContainsKey($name)) { continue }
            $previousEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
            [Environment]::SetEnvironmentVariable($name, $Matches[2], 'Process')
        }
    }

    New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
    $logPath = Join-Path $buildDir 'build.log'
    $configureArguments = @(
        '-S', $repoRoot, '-B', $buildDir, '-G', 'Ninja',
        "-DCMAKE_BUILD_TYPE=$BuildType", '-DCMAKE_C_COMPILER=cl', '-DCMAKE_Fortran_COMPILER=ifx',
        '-DMCTCLIB_FIND_METHOD=subproject', '-DWITH_TBLITE=OFF', '-DWITH_CPCMX=OFF',
        '-DWITH_TESTS=OFF', '-DWITH_JSON=OFF', '-DWITH_OBJECT=OFF',
        '-DBLA_VENDOR=Intel10_64lp_seq', '-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF',
        '-DCMAKE_EXE_LINKER_FLAGS=/Qoption,link,/STACK:67108864'
    )
    if ($BuildType -eq 'Release') {
        $configureArguments += '-DCMAKE_Fortran_FLAGS_RELEASE=/O3 /DNDEBUG'
    }
    try {
        Invoke-Checked 'cmake' $configureArguments *> $logPath
        Invoke-Checked 'cmake' @('--build', $buildDir, '--target', 'fast-gfn1-xtb', '--parallel', "$Parallel") *>> $logPath
    } catch {
        Get-Content $logPath -Tail 60
        throw
    }
    Write-Host "PASS Windows $BuildType build: $(Join-Path $buildDir 'fast-gfn1-xtb.exe')"
    Invoke-Checked $PY @('--version')
    Invoke-Checked $PY @('-c', 'import sys; print(sys.executable)')
    Invoke-Checked $PY @((Join-Path $PSScriptRoot 'test_source_guards.py'))
    Invoke-Checked $PY @((Join-Path $PSScriptRoot 'test_compact_density_math.py'))

    $policyDir = Join-Path $buildDir 'policy-test'
    New-Item -ItemType Directory -Force -Path $policyDir | Out-Null
    Push-Location $policyDir
    try {
        Invoke-Checked 'ifx' @('/nologo', '/check:all', '/fpe:0',
            (Join-Path $repoRoot 'src/gfn1_fast_policy.f90'),
            (Join-Path $PSScriptRoot 'test_partial_dispatch_policy.f90'),
            '/exe:test_partial_dispatch_policy.exe')
        Invoke-Checked (Join-Path $policyDir 'test_partial_dispatch_policy.exe') @()
    } finally { Pop-Location }

    $smokeArguments = @((Join-Path $PSScriptRoot 'test_windows_smoke.py'),
        '--exe', (Join-Path $buildDir 'fast-gfn1-xtb.exe'), '--output', (Join-Path $buildDir 'smoke'))
    # Full runtime checking is expensive for the 350-AO gradient. Release
    # exercises the large-system parity cases; Debug checks startup and ALPB.
    if ($BuildType -eq 'Debug') { $smokeArguments += '--small-only' }
    Invoke-Checked $PY $smokeArguments
} finally {
    foreach ($name in $previousEnvironment.Keys) {
        [Environment]::SetEnvironmentVariable($name, $previousEnvironment[$name], 'Process')
    }
}
