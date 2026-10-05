# GFN1-fast 4.3.0 halogen Hessian performance (2026-10-05)

Windows native ifx 2025.2 `/O3`, oneMKL LP64 sequential, one thread. Actual CLI
process wall time includes startup, output, Hessian projection and frequencies.
Identical input, GFN1 parameters, SCC accuracy 1e-7 and numerical step 4e-4 bohr.
One warmup per policy, seven randomly interleaved trials, profiling off.

| Complex (8 atoms) | Standard numerical / s | Fast numerical / s | Analytic / s | Speed vs standard | Speed vs fast numerical |
|---|---:|---:|---:|---:|---:|
| CH3Br + water | 0.176210 | 0.137331 | 0.044180 | 3.99× | 3.11× |
| CH3I + water | 0.174350 | 0.134881 | 0.043821 | 3.98× | 3.08× |

Every timed trial compares the actual matrices with the standard numerical
implementation: fast numerical matrices agree exactly at written precision;
analytic maximum difference is `6.95e-8 Eh/bohr²`. Reported Energy/Gradient
norms agree within 1e-10. Separate profiled calls verify analytic dispatch and
one active halogen triplet. Independent unprojected full-Gradient derivative
tests are described in `RESULTS.md`.

The standard reference is official xTB v6.7.1 commit
`26b28010e805f7d1aeeef39813feb473e69cc4be`, built with the same compiler and
BLAS. Its Hessian routine's OpenMP directives are suppressed for this
**one-thread comparison** to avoid Windows ifx's nested-allocatable descriptor
crash. The original displacement, SCC, Gradient and finite-difference
arithmetic are retained. This makes no claim about stock parallel performance.
The exact compatibility patch and reproduction script are in
`../performance/reference_compatibility.patch` and
`../performance/run_windows_reference.ps1`.

The correction extends the existing analytic strategy with measured benefit
on these complexes. It remains available by explicit enablement while solvent
response and broader molecule coverage are developed. Reuse and dense FOE
remain disabled by default because their measured whole-calculation cost is
higher than SYEVD.

Raw samples and hashes are preserved in `BENCHMARK_RESULTS_20261005.json`;
corresponding source hashes in `SOURCE_MANIFEST_20261005.json`. Re-run:

```powershell
.\tests\gfn1_fast_4_3_0\run_windows_regression.ps1 -BuildType Release -Benchmark
```

Per-trial matrices/logs are in the current Release build's `halogen-cli-test/`
directories. Future reruns can update those mutable directories; the dated
JSON is a snapshot of the measurements reported above.
