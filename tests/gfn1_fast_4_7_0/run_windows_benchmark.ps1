$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$exe = Join-Path $repoRoot '_reference/gfn1-fast-4.7.0-source-snapshot/binaries/release/xtb.exe'
$reference = Join-Path $repoRoot '_reference/gfn1-fast-4.6.0-source-snapshot/binaries/release/xtb.exe'
$output = Join-Path $repoRoot 'build-gfn1-fast-current-windows-ifx-release/parallel-frozen-benchmark-4.7.0'
& cmd /d /c "call `"$setvars`" intel64 >nul && `"$PY`" `"$PSScriptRoot/benchmark_parallel.py`" --exe `"$exe`" --reference `"$reference`" --output `"$output`" --repeats 5"
if ($LASTEXITCODE -ne 0) { throw 'Frozen analytic Hessian parallel benchmark failed' }
