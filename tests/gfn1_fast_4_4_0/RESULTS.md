# GFN1-fast 4.4.0 / 2.4.1: Gradient at overlap nodes

Verified on Windows native on 2026-10-05 using ifx 2025.2.0, oneMKL sequential
and the prescribed Conda Python 3.13.15. Both complete Release and Debug builds
were rebuilt from scratch (782 steps each). This is a development snapshot.

## Correction

The previous Gradient reconstructed the Hamiltonian prefactor from H0/S.
When overlap S is exactly zero, H0 is also zero but its Cartesian derivative
can be nonzero. Setting this prefactor to zero discards a physical derivative.
An actual converged asymmetric water SCC reproduced a Gradient error of
`0.046801854207 Eh/bohr`; independent Energy differences confirmed it.

The production GFN1 kernel now constructs the prefactor directly from shell
energies, shell scaling, Slater exponents and the distance polynomial. Polynomial
and CN contractions also avoid division by S, the polynomial or average shell
energy. `XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT=1` restores the diagnostic
legacy behavior. The generic Gradient override also retains the legacy path.

Small nonzero overlaps screened out by the existing SCF cutoff remain distinct
from exact nodes: raw packed H0 classifies this case and preserves the prior
screening behavior. Ordinary water/taxol Energy and Gradient are preserved to
the existing strict tolerances. Analytic Hessian dispatch rejects legacy
Gradient policies and uses numerical fallback.

## Evidence

- Strict native nodal tests, both builds: actual SCC accuracy `1e-7`, all nine
  Cartesian Energy/Gradient differences, raw Hessian and rotation covariance.
  Corrected production/analytic Gradient difference is below `3e-16 Eh/bohr`.
  Energy derivative error is `9.61e-9 Eh/bohr` at step `2e-4 bohr`.
  Hessian errors decrease from `1.35e-7` to `3.37e-8 Eh/bohr²` when halving
  the step. The separate legacy process reproduces the original defect.
- Actual nodal CLI Hessian, both builds: analytic/numerical difference
  `1.176e-7 Eh/bohr²`; legacy and generic policies produce identical fallback
  matrices. CLI all-coordinate Energy derivative errors: gas `1.08e-8`,
  ALPB `9.55e-9`, smooth nodal GBSA `9.27e-9 Eh/bohr`.
- CLI `--acc` is clamped to a minimum of `1e-4`. The Energy derivative CLI
  tests now explicitly request this actual value. Native API tests use `1e-7`.
  Hessian CLI tests set their separate `$hess sccacc=1e-7` control.
- Strict GBSA API tests, both builds, four steps down to `5e-5 bohr`:
  smooth nodal geometry has maximum Energy derivative error below
  `6.2e-10 Eh/bohr`; fixed-charge solvent error below `3.6e-12 Eh/bohr`.
- Rebuilt ordinary Release smoke: seven cases; Debug smoke: gas/ALPB water.
  Energy differences are zero at output precision; largest path Gradient
  difference `2.005e-15 Eh/bohr`. Existing source/math/dispatcher checks pass.
- Debug gas CLI thermal/open-shell/default/fallback cases and seven halogen
  molecular/CLI cases pass with the corrected Gradient. Smaller-step halogen
  raw Hessian error is at most `1.99e-8 Eh/bohr²`. Partial additive calculator
  API, restart preservation and polarizability fallback pass in both builds.

## Existing GBSA surface discontinuity

The original nodal fixture is close to a screened surface-quadrature branch.
`compute_numsa` omits points whose surface weight is at most `1e-6`.
Crossing that threshold changes the SASA and hydrogen-bond energies; a finite
difference crossing the branch need not equal a derivative within one branch.

The isolated fixed-charge solvent reproduces the full SCC error, excluding
electronic or charge-response errors. Component tests isolate SASA/H-bond
terms; the Born derivative remains within `1e-9 Eh/bohr`. At step `5e-5`,
the full/fixed errors are respectively `4.252e-7` / `4.252e-7 Eh/bohr`.
The official stock executable independently produces exactly the same
Energy values and derivative errors in the tested z component at all four
steps. The stock executable's different full Gradient norm is due to its
separate overlap-node defect.

This boundary fixture remains an explicit reproduction test. The smooth
GBSA fixture changes only the second H y-coordinate from 0.65 to 0.67 Å;
the first H remains at y=0, retaining the exact O_py/H_s overlap node.
No tolerance was relaxed and no solvent theory/threshold was changed.
Future solvent analytic Hessians must reject non-smooth screening branches.

## Re-run and artifacts

```powershell
.\tests\gfn1_fast_4_4_0\run_windows_nodal.ps1 -BuildType Release
.\tests\gfn1_fast_4_4_0\run_windows_solvent.ps1 -BuildType Release
.\tests\gfn1_fast_4_4_0\run_windows_regression.ps1 -BuildType Release
# Repeat with -BuildType Debug.
.\tests\gfn1_fast_4_4_0\run_windows_regression.ps1 -BuildType Release -Benchmark
```

Dated copies of native logs and CLI summaries are in this directory. Source
and executable hashes are in `SOURCE_MANIFEST_20261005.json`. Compiled source
and executable snapshots are preserved under the ignored
`_reference/gfn1-fast-4.4.0-source-snapshot`; this is a recovery snapshot, not
a complete distributable checkout. The repeated taxol full-Hessian comparison
has completed: five profile-off trials, stock/analytic median times
250.238602 / 30.070810 seconds, **8.32×** speedup. All full matrices agree to
`2.991e-7 Eh/bohr²`; the separate profile audit confirms actual analytic
dispatch. See `BENCHMARK.md` and `LARGE_BENCHMARK_20261005.json` for scope,
reference compatibility and raw evidence.

Analytic Hessians remain opt-in. Solvent analytic response, broader screening
consistency and performance-based default adoption remain project work.
