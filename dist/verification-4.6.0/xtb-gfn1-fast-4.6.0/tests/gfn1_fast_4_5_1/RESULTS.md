# GFN1-fast 4.5.1: guarded solvent response and resource order

The deeply buried surface-point prepass proves an entire point contribution
remains zero across the requested displacement window before checking branches
of its other neighbours. This removes irrelevant branch rejections and dense
point-derivative work. Active surface points retain the original screening,
branch guards and tolerances. The original SASA model and SCC algorithm are
unchanged.

Solvent geometry temporaries now finish before the large AO/CN response arrays
are allocated. The workspace bound uses the larger of these two stages, plus
persistent solvent fields and caches. The default cap remains 1024 MiB.

## Exact source and binary scope

The current working tree was clean-built on Windows native with ifx 2025.2.0,
MSVC 19.44.35217, sequential LP64 MKL, CMake 4.1.1 and Ninja 1.12.1.
Python is the prescribed Windows Conda `fastxtb` 3.13.15 executable.

- Release SHA-256: `bcb61b376e0612704901e51f44319db458b1656b4e077685ec29e4d40080048e`.
- Debug SHA-256: `c767c06d8f65b729189420cbfe096d49e13eba5102200a0f499c9d95544e853f`.
- `SOURCE_MANIFEST_20261005.json` verifies the compiled sources and frozen
  executables in `_reference/gfn1-fast-4.5.1-source-snapshot`. This is a recovery
  snapshot, not a standalone distributable checkout.
- `PRE_REORDER_TAXOL_20261005.json` and the `pre_reorder_*` files describe an
  earlier allocation order. Its fixed-charge probe log is from the original
  surface implementation; it is diagnostic historical evidence.

## Correctness after the allocation change

Release Energy/Gradient smoke tests pass all seven cases; Debug water gas/ALPB
passes with `/check:all /fpe:0`. Source, weighted-density math and compiled
dispatcher tests pass. Both builds pass all 15 actual solvent CLI Hessian
conditions: 10 analytic and five fallback cases. CLI matrix differences are
at most `1.94e-8 Eh/bohrﾂｲ`. The detailed native-recheck evidence is referenced
by `tests/WINDOWS_CURRENT_TEST_REPORT.md`.

The final Release library passes independent converged nonlinear SCC finite
differences for every raw Hessian column of taxol (113 atoms, 350 AO). These
checks do not symmetrize or project the analytic result. Reciprocity,
translation and original Gradient parity checks also pass.

| Native solvent model | Columns | Step / bohr | Maximum raw H difference / Eh/bohrﾂｲ | Original Gradient difference / Eh/bohr |
|---|---:|---:|---:|---:|
| GBSA Still | 339 | 1e-5 | 2.85377e-8 | 1.35525e-15 |
| ALPB P16 | 339 | 1e-5 | 1.59194e-8 | 1.01221e-15 |

`LARGE_NATIVE_RESULTS_20261005.json` records exact library/driver/source hashes
and all results; `large-molecular-solvent-20261005.log` preserves the output.
Raw Energy finite-difference discrepancies are also recorded (about 5e-8
Eh/bohr); they are not asserted as satisfying a tighter Energy-derivative
bound. The original finite-temperature occupation and screened surface
behavior remains unchanged. CPU times in this native correctness driver are
diagnostics, not repeated CLI performance measurements.

## Small-system thread crossover

`THREAD_CROSSOVER_20261005.json` contains seven randomized interleaved trials,
one warmup, profile-off wall timings and separate successful profile audits
at 1, 4 and 8 OpenMP threads. MKL uses one thread. Compare the fast numerical
and analytic algorithms at identical thread counts; no stock multithread
performance claim is made.

| Case | Analytic speedup, 1 thread | 4 threads | 8 threads |
|---|---:|---:|---:|
| Water, gas, strict SCC | 1.205 | 1.258 | 1.224 |
| Water, ordinary CLI defaults | 1.145 | 1.235 | 1.279 |
| Disilane, gas, strict SCC | 3.403 | 1.604 | 1.357 |
| Disilane, ordinary CLI defaults | 2.630 | 1.494 | 1.216 |
| Benzene, gas, strict SCC | 5.790 | 2.145 | 1.651 |
| Disilane, ALPB | 3.261 | 1.554 | 1.316 |

Numerical 1/4/8-thread matrices agree exactly in these measurements. Strict
gas steps are `5e-4 bohr`, ALPB `1e-4 bohr`; maximum differences are respectively
`1.67e-7` and `6.00e-9 Eh/bohrﾂｲ`. Ordinary defaults use the existing `0.005 bohr`
numerical step and normal SCC accuracy; their maximum difference is
`2.75135e-5 Eh/bohrﾂｲ`. The analytic response avoids that finite-step truncation
error. These tests do not establish universal speedups across all molecules.

## Large CLI resource and original-code comparison

The new `tests/performance/run_windows_large_solvent_hessian.ps1` first audits
the real CLI dispatch under the unchanged 1 GiB cap. It then measures five
randomized interleaved stock/analytic trials after one warmup for each of GBSA
and ALPB, using the same taxol input, `1e-5 bohr` step and `1e-7` Hessian SCC
accuracy. Each trial checks the entire 339-by-339 written matrix and displayed
Energy/Gradient norm. Windows process counters record peak working set,
peak pagefile usage and sampled private bytes. Timings use frozen executables.

The stock comparison preserves numerical Hessian formulas and displacements;
only Hessian OpenMP directives are removed to avoid the documented Windows
ifx allocation-descriptor crash. This is a one-thread comparison.

The completed frozen-binary CLI audits are in `LARGE_CLI_AUDIT_20261005.json`:

| Model | Analytic used | Estimated workspace / bytes | Peak working set / bytes | Peak pagefile usage / bytes |
|---|---|---:|---:|---:|
| GBSA | true | 880428680 | 837627904 | 857718784 |
| ALPB | true | 880428680 | 837058560 | 857808896 |

Both audited processes use the unchanged default cap of `1073741824 bytes`.
These process counters and the conservative workspace estimate measure
different quantities; neither is a universal bound on every input.

The five-repeat original-code comparison has completed for both models.
`LARGE_BENCHMARK_20261005.json` preserves all trials, full matrix errors,
input/binary hashes and process resource counters; `large-benchmark-20261005.log`
preserves the complete run log.

| Model | Stock median / s | Analytic median / s | Speedup | Maximum complete-matrix difference / Eh/bohr² |
|---|---:|---:|---:|---:|
| GBSA water, taxol | 229.775900 | 35.167327 | 6.533789 | 7.20e-9 |
| ALPB water, taxol | 229.952206 | 33.692745 | 6.824977 | 1.00e-8 |

These are frozen 4.5.1, one-thread, profile-off results at `1e-5 bohr` and
`1e-7` Hessian SCC accuracy. Resource audits are separate from timing results.
Analytic Hessian remained opt-in in this measured binary. Subsequent production
policy and screening/cutoff rejection changes require their own validation.
Salt, unsupported models, branch-window failures and insufficient memory
retain numerical fallback. Distribution validation remains outstanding.
