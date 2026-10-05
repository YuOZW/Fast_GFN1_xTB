# 4.5.0 solvent Hessian performance

Windows native ifx 2025.2, i7-14700F, oneMKL LP64 sequential. Each fixture/model
uses one warmup and seven randomized interleaved trials of stock numerical,
fast numerical and fast analytic Hessians. One thread, profile disabled while
timing; a separate profile audit confirms complete solvent analytic dispatch.
All timings are process wall time, including startup and the whole CLI run.
No builds or other runtime tests ran during these measurements.

| Fixture | Solvent | Stock numerical | Fast numerical | Fast analytic | Speedup over stock |
|---|---|---:|---:|---:|---:|
| Disilane, 8 atoms / 30 AO | GBSA water | 0.200495 s | 0.158725 s | 0.046011 s | 4.36× |
| Disilane, 8 atoms / 30 AO | ALPB water | 0.202878 s | 0.154238 s | 0.046431 s | 4.37× |
| Smooth CH3Br-water, 8 atoms | GBSA water | 0.168987 s | 0.131437 s | 0.043943 s | 3.85× |
| Smooth CH3Br-water, 8 atoms | ALPB water | 0.172349 s | 0.132788 s | 0.045953 s | 3.75× |

Analytic speedup over the existing fast numerical path is respectively 3.45×,
3.32×, 2.99× and 2.89×. All trials compare the entire 24×24 output Hessian:
analytic/stock maximum difference is at most 6.7e-9 Eh/bohr²; numerical/stock
at most 1e-10. Printed Energy and Gradient norm agree within 1e-10 each trial.

**Stock compatibility scope:** official xTB 6.7.1 commit
`26b28010e805f7d1aeeef39813feb473e69cc4be`, built with matching Windows flags,
parameters, oneMKL and stack settings. Its numerical Hessian crashes under
ifx at entry to the original OpenMP region with nested allocatable descriptors.
Only Hessian OpenMP directives are suppressed in the comparison build, retaining
original displacement/SCC/Gradient arithmetic and postprocessing. This is a
one-thread comparison; it does not establish stock multithread performance.

The common Hessian SCC accuracy is 1e-7 and numerical step/analytic branch
window is 1e-4 bohr. The CLI's ordinary `--acc` clamp to 1e-4 is distinct from
the `$hess sccacc` control. GBSA uses Still and ALPB the CLI default P16.
The smooth halogen coordinates were fixed using original SASA branch margins,
not timed-path or Hessian-error selection. Tests also retain the original
halogen/water surface boundary and the disilane ALPB 2e-4 conservative fallback.
These measurements support the opt-in solvent analytic path for the tested
smooth settings; they do not prove speed at default steps or large N.

Run:

```powershell
.\tests\gfn1_fast_4_5_0\run_windows_regression.ps1 -BuildType Release -Benchmark
```

Samples, randomized order, per-trial Energy/Gradient observations and fixed
executable hashes are preserved in `BENCHMARK_RESULTS_20261005.json`; the
complete log is `benchmark_run_20261005.log`. The measured Release SHA-256 is
`193cae00930589ac20108c184d65d6377a2366f77ba8527ee945e31fac7dab78`.
