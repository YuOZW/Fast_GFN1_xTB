# GFN1-fast 4.2.0: molecular analytic Hessian CLI dispatch

Development snapshot, 2026-10-05; the complete requested project remains open.

`TxTBCalculator%hessian` now dispatches a supported gas-phase GFN1 reference to
the analytic molecular Hessian assembled from canonical electronic P/W response,
coupled shell charge response, overlap/H0/CN geometry derivatives, Coulomb,
repulsion and pairwise D3. It adds only the requested atom columns to the
Hessian, preserves unrequested columns and the supplied restart reference,
and assigns the requested dipole columns to zero, matching the existing
Energy/Gradient-only SCF's deliberately zero dipole response.

Enable explicitly with `XTB_GFN1_FAST_ENABLE_ANALYTIC_HESSIAN=1`.
`XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN=1` takes precedence. The default is
disabled while broader model/size validation is completed. No SCC mixing,
occupation, physical response denominators or Energy/Gradient tolerances change.

Solvation, halogen corrections (Cl/Br/I/At), PBC, point charges/electric fields,
constraints/fixed/frozen atoms, wall/metadynamic potentials and requested
polarizability derivatives dispatch to the original numerical formula.
Invalid references, failed response gates, allocation failure, asymmetric or
non-translational Hessians also fall back. A conservative workspace estimate
is limited by `XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB` (default 1024 MiB).
Reciprocity/translation checks inspect the raw matrix before any CLI projection.
No post-hoc symmetrization is used to pass these gates.

Profile (`XTB_GFN1_FAST_PROFILE=1`) reports actual analytic use, fallback reason,
estimated workspace and analytic attempt wall time. The CLI's existing
`Numerical Hessian` heading and final timing label predate this dispatcher;
use the explicit `4.2.0 analytic Hessian: used=T/F` diagnostic to identify it.

## Actual CLI and API checks

Release and Debug current-tree builds and their ordinary Energy/Gradient tests
pass. Actual `--hess` commands are used, without combining `--grad` and `--hess`
(the run-type selection would otherwise retain Gradient only).

All six Release analytic/numerical CLI comparisons passed with Hessian step
5e-4 bohr and reference/displacement SCC accuracy 1e-7:

| Input | Largest Cartesian Hessian difference, Eh/bohr² |
|---|---:|
| Water, 300 K | 1.259e-7 |
| Water, 0 K | 1.259e-7 |
| Water, 30000 K | 9.420e-8 |
| Water +1, open shell, 300 K | 1.298e-7 |
| Disilane, 300 K | 1.487e-7 |
| Benzene, 300 K | 1.670e-7 |

These are the CLI's projected, symmetrized Cartesian output matrices. The
independent 4.1.0 molecular tests additionally compare the **unprojected raw**
analytic matrix to every coordinate of the converged Gradient finite
difference at two step sizes; errors decrease quadratically. Those strict
Release/Debug tests cover water at 0/30000 K, open-shell water at 3000 K and
disilane at 300 K, with smaller-step full-Hessian error at most 4.17e-8 Eh/bohr².

Ordinary CLI defaults are also checked separately for water and disilane:
no `--acc` override and no Hessian step/SCC settings. Their explicit analytic
request is accepted and agrees with the numerical Hessian at the default
0.005 bohr step within 2e-4 Eh/bohr² (finite-difference truncation error is
larger at that step). These checks supplement the tighter validation above.

ALPB water, HCl, a distance constraint, explicit disable and a 1 MiB memory
budget each exercise numerical fallback and reproduce the numerical
reference exactly at printed Hessian precision. Disilane's 4-thread numerical
fallback reports `outer OpenMP=T` and matches the 1-thread matrix exactly.
Debug repeats the water thermal/open-shell and ALPB/HCl/constraint/disable
cases under bounds/allocation checks and floating-point traps.

Strict Fortran calculator API tests pass in both Release/Debug: partial atom
lists, additive Hessian columns, unrequested Hessian/dipole columns and
restart density/charges are preserved. A requested polarizability array forces
numerical fallback and preserves its unrequested columns. The runner checks
actual profile output is T for the analytic call and F for this fallback.

## Windows numerical fallback repair

The existing numerical Hessian crashed during Windows ifx OpenMP entry while
privatizing derived scratch with nested allocatable descriptors, even with a
serialized region. Each worker now constructs its scratch in a separate
internal routine, retains the existing displacement work sharing and computes
the same +/- Gradient formula. Optional polarizability output uses a valid
shared buffer, copied back outside the region. One- and four-thread numerical
results above verify the repaired path. The active wall count is initialized
and reset separately from allocated wall capacity, so unused capacity does not
incorrectly disable the analytic path.

## Reproduce

From the repository root (Windows native, specified fastxtb Python):

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Release
.\tests\gfn1_fast_4_2_0\run_windows_regression.ps1 -Benchmark
.\tests\gfn1_fast_4_2_0\run_windows_calculator.ps1
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Debug
.\tests\gfn1_fast_4_2_0\run_windows_regression.ps1 -BuildType Debug
.\tests\gfn1_fast_4_2_0\run_windows_calculator.ps1 -BuildType Debug
```

The benchmark requires the standard reference build from
`tests/performance/run_windows_reference.ps1`. CLI logs/matrices and summary
are in each current build's `hessian-cli-test`; benchmark samples are saved
separately in `benchmark_summary.json` so a regression rerun preserves them.
API logs are in
`hessian-calculator-test/calculator-test.log`. See `BENCHMARK.md` for exact
reference compatibility changes and measured speed, without extrapolating to
large molecules or solvents.

## Remaining work

Complete halogen and solvent/Born response, extend molecular/size coverage,
resolve screened/nodal-overlap derivatives consistently with the Gradient,
reduce dense derivative storage and establish a production crossover. Finish
versioned source artifacts, ZIP/re-extraction checks and SHA deliverables.
The feature is usable on validated gas systems but does not prove the whole
4.x/project scope complete.
