# GFN1-fast 4.7.0 OpenMP implementation and validation

Windows native, Intel ifx 2025.2 / oneMKL sequential, i7-14700F,
Python 3.13.15 in the prescribed fastxtb environment. WSL/GPU are not used.

## Implementation

- Parallel shell susceptibility construction and coordinate response columns;
  initialized electronic/SCC caches are immutable, with worker-local scratch.
- Parallel atomic AO integral derivatives, per-worker geometric/CN reductions,
  and disjoint per-atom SASA derivatives. No OpenMP copies of nested allocatable
  wavefunction or solvent descriptors.
- Four-coordinate BLAS trace contractions reuse the existing derivative
  tensors, rather than allocating a complete AO response tensor.
- Worker-local moving-overlap transforms are reused between the two coupled
  electronic applications. Zero-overlap susceptibility columns skip that
  transform. Canonical finite-temperature P/W response remains unchanged.
- Born pair and SASA screening factor Hessians scatter six-coordinate blocks,
  keeping the complete product Hessian and all original branch-window guards.
- One resource selector serves the native API and CLI. Nested OpenMP and
  systems below 64 AO use one response worker; otherwise the maximum is eight,
  limited by available runtime threads, shells/coordinates and the 1 GiB
  analytic workspace estimate. A smaller budget reduces workers before falling
  back to numerical differentiation if the serial response does not fit.
- Concurrent policy initialization/readers use a common lock. There are no
  changes to the final SCF density, eigensolver default or Energy/Gradient
  formulas. Existing response-conditioning and fallback tolerances are kept.

## Correctness and guards

Release and Debug pass 19 actual parallel CLI comparisons, including 2/4-worker
teams, eight-worker disilane, a two-worker policy cap, a two-worker OpenMP runtime
limit, an explicit response-serial switch, zero/high temperature and open shell.
Four parallel fallback cases (memory, SASA branch, screened overlap, constraints)
match their numerical baselines. Printed full matrices agree within 1e-10.

The existing 21 default dispatch cases, screening and D3 cutoff/window checks
also pass in both builds. All 38 native solvent molecular conditions pass with
forced parallel response teams; these independently differentiate converged
nonlinear SCC gradients for every Cartesian column. The gas native test checks
the optional density and shell-charge response tensors against SCC differences.
A real nested analytic API call and partial additive calculator columns pass.
Resource tests cover concurrent first calls, 900 MiB worker reduction, 1 MiB
rejection, nested operation, and a standalone compile without OpenMP.

For taxol (113 atoms, 350 AO), eight-thread native analytic Hessians pass all
339 raw, unprojected Gradient-FD columns for both solvent models, with no Hessian
symmetrization in this validation:

| Model | Step (bohr) | Maximum raw H error (Eh/bohr^2) |
|---|---:|---:|
| GBSA | 1e-5 | 9.4174784e-9 |
| ALPB | 1e-5 | 1.5953890e-8 |

Original production/response Gradient differences are about 2e-15 Eh/bohr.
The acceptance tolerance remains 2e-6 Eh/bohr^2; none was relaxed.

## Repeated frozen-binary performance

Timing uses complete taxol CLI `--hess`, the same input/parameters, 1e-5 bohr
and SCC accuracy 1e-7. Each configuration has one warm-up and five randomized
measured repetitions. Profiling is off during timing; separate profiled runs
verify actual team sizes and resources. MKL uses one thread. Every measured
full matrix, original Energy and Gradient norm is checked against frozen
4.6.0. This compares versions of Fast GFN1 and does not claim a multithreaded
stock-xTB comparison.

| Model | 4.6.0, 1 thread (s) | 4.7.0, 1 thread (s) | 4.7.0, 4 threads (s) | 4.7.0, 8 threads (s) |
|---|---:|---:|---:|---:|
| Gas | 30.83 | 15.60 | 5.54 | 4.27 |
| GBSA | 35.11 | 16.21 | 5.74 | 4.51 |
| ALPB | 33.86 | 16.00 | 5.81 | 4.55 |

Eight-thread speedups over 4.6.0 are 7.22x / 7.78x / 7.44x (gas/GBSA/ALPB);
speedups over the optimized serial path are 3.65x / 3.59x / 3.51x. The optimized
serial path alone is 1.98x / 2.17x / 2.12x faster than 4.6.0. Eight workers are
faster than four in all three five-repeat medians, so the eight-worker maximum
is adopted. Small systems retain the 64-AO crossover.

All printed full matrices have zero difference from 4.6.0 in every trial;
Energy/Gradient norm agree within 1e-10. The largest measured eight-worker
Windows peak pagefile-usage counter is 1,039,687,680 bytes (about 0.968 GiB),
below the 1 GiB value on this input. This is a process commit metric; the
policy cap is an analytic workspace estimate, not a universal whole-process
limit. `BENCHMARK_20261005.json` includes every trial and resource record.
`VALIDATION_20261005.json` and the per-build CLI records bind correctness
evidence to the frozen source and executable hashes.

Distribution: `dist/xtb-gfn1-fast-4.7.0.zip`, with pinned mctc-lib source.
The external `dist/REEXTRACTION_4.7.0.json` records source/dependency checksums,
fresh native Release/Debug builds and actual parallel/native FD calculations.
The ZIP is created once; re-extraction evidence is written outside it.

Frozen source/test-code manifest: `SOURCE_MANIFEST_20261005.json` (603 files).
Release executable SHA-256:
`ebf4e04d6d85dc22d1d07f379389ae164e0d2b6b56e13a54901421fcbf349f37`.
Debug executable SHA-256:
`eccf25cc2d63323ef9a42d43f930f26cf80b63ce68e0fbe9e0d314a1efa2390f`.
Earlier source snapshots and distribution packages are preserved.
