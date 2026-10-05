# Native GFN1 Hessian performance (2026-10-05)

Windows ifx 2025.2 `/O3`, oneMKL LP64 sequential, one thread, the same GFN1
parameters/input and accuracy, no fast-math/LTO. Actual CLI process wall
includes startup/output, Hessian projection and frequency postprocessing.
One warmup per policy, then seven randomized interleaved trials; profiling
is off while timing. Numerical Hessian step is 5e-4 bohr and SCC accuracy 1e-7.

| Molecule | Standard numerical, s | Fast numerical, s | Analytic, s | Speed vs standard | Speed vs fast numerical |
|---|---:|---:|---:|---:|---:|
| Disilane (8 atoms / 30 AO) | 0.199770 | 0.150448 | 0.044638 | 4.48× | 3.37× |
| Benzene (12 atoms / 30 AO) | 0.406882 | 0.295154 | 0.051674 | 7.87× | 5.71× |

Every trial checks the actual Hessian against the standard numerical output:
fast numerical error at most 1e-10, analytic error at most 1.67e-7 Eh/bohr².
Separate profiled calls confirm the timed analytic policy used the analytic
implementation rather than silently falling back. This comparison covers
small gas-phase molecules; no large-molecule or solvent analytic speed claim.
The fast executable retains the existing Energy/Gradient-only property policy;
full stock dipole/bond-order/IR property compatibility is outside project scope.

## Reviewable standard reference

Official grimme-lab/xtb v6.7.1 commit
`26b28010e805f7d1aeeef39813feb473e69cc4be`, built separately with identical
compiler/BLAS and a 64 MiB executable stack. Compatibility changes are Windows
flag spelling, HOME/USERPROFILE initialization, gating unused test targets,
and the following numerical-Hessian workaround.

The official OpenMP Hessian also crashes at entry with this ifx runtime while
privatizing nested allocatable descriptors. Its OpenMP directives **only in
the Hessian routine** are suppressed for this one-thread reference. The
original +/- displacement sequence, molecular/restart copying, full SCC
Gradient and derivative arithmetic are unchanged. Parallelism elsewhere
remains enabled. This is a serial numerical reference; it makes no statement
about stock parallel Hessian performance. The executable's Energy/Gradient
algorithms and property calculations retain the official implementation.

Measured executable SHA-256:

- Fast: `0f36ab0a4de856726a604353f9973ccb8b04c80bbca3784d2f251fda45b24bf2`
- Stock: `4ee8e9e26109625e51398dcf1e6552cebbcc66ccc62fba667f5674d2af8a292a`

Raw seven samples, per-trial matrices/logs and hashes are in the current
Release build's `hessian-cli-test/benchmark_summary.json` and `bench_*` directories.
`tests/performance/run_windows_reference.ps1` reproduces every compatibility
change. Re-run with `tests/gfn1_fast_4_2_0/run_windows_regression.ps1 -Benchmark`.

`BENCHMARK_RESULTS_20261005.json` preserves the dated raw samples independently
of mutable build directories. The corresponding source hashes and full stock
compatibility diff are saved in `../performance/SOURCE_MANIFEST_20261005.json`
and `../performance/reference_compatibility.patch`.
