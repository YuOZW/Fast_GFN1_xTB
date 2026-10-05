# GFN1-fast 4.3.0: gas-phase halogen Hessian

Verified on Windows native, 2026-10-05, Intel ifx 2025.2 / oneMKL sequential,
with the prescribed Conda Python 3.13.15. Release and Debug libraries/executables
were rebuilt. Strict component programs additionally use `/Od /check:all /fpe:0`.

## Change

`xbMolecularHessian` differentiates the original GFN1 B-X...A correction
analytically, including its damped Lennard-Jones radial function and sixth-power
angular factor. It retains SCF's donor/acceptor groups, 20-bohr cutoff, nearest
atom selection and exponent 12. The original `xbpot` Energy/Gradient is unchanged.
The generic kernel also supports custom positive exponents and donor parameters.

The complete gas-phase response now adds this correction to electronic/SCC,
Coulomb, repulsion and pairwise D3 Hessians. The calculator's blanket halogen
rejection is removed; `XTB_GFN1_FAST_ENABLE_ANALYTIC_HESSIAN=1` enables the path.
The profile reports `GFN1-fast 4.3.0 analytic Hessian: used=T/F` and the number
of nonzero-strength, distinct-neighbour halogen triplets assembled.

Equal nearest-neighbour distances, coincident active atoms and cutoff boundaries
reject the analytic kernel. The calculator then uses the original numerical
derivative. If the acceptor is itself the unique nearest atom, its angular
contribution is identically zero locally. Default GFN1 Cl strength is zero;
Br, I and At have nonzero default corrections in the loaded parameter set.

## Evidence

- **96 component cases**, both Release and Debug: Cl/Br/I/At donors, N/O/P/S
  acceptors, three geometries, default and custom parameters/exponent. Compared
  directly with the original `xbpot`: maximum Energy difference `4.34e-18 Eh`,
  Gradient difference `1.05e-17 Eh/bohr`. All-coordinate original Gradient
  finite differences at 1e-3 and 5e-4 bohr: largest error `1.41e-8 Eh/bohr²`,
  with quadratic convergence. Raw symmetry, translation and branch guards pass.
- **Seven molecular cases**, both configurations: CH3Cl/Br/I/At + water at
  300 K, Br at 0/30000 K, I cation open shell at 3000 K. Real GFN1 SCC, all
  24 Cartesian Gradient columns, steps 4e-4/2e-4 bohr. Nonzero production XB
  Energy and analytic triplet counts explicitly checked for Br/I/At. Maximum
  complete Gradient difference `1.39e-14 Eh/bohr`; smaller-step raw Hessian
  error `1.99e-8 Eh/bohr²`. Halving the step reduces error by about four.
  No post hoc symmetrization/projection in this molecular test.
- **Actual CLI**, both configurations: the same seven cases use analytic
  Hessians; maximum difference from numerical Hessian `6.95e-8 Eh/bohr²`.
  Nearest-distance tie and exact cutoff fixtures fall back and produce the
  identical numerical matrix. This CLI comparison includes normal projection.
- **Existing CLI/API regressions**, rebuilt libraries: water thermal/open-shell,
  disilane, benzene, HCl, default settings; solvent, constraints, disable,
  memory-budget and polarizability fallback; partial-column additive API,
  restart preservation, numerical 1/4-thread agreement. Release/Debug pass
  their respective established scopes. HCl now exercises supported gas-phase
  dispatch. Energy/Gradient smoke and existing source/math/policy guards pass.

Re-run from the repository root, after the current native build:

```powershell
.\tests\gfn1_fast_4_3_0\run_windows_halogen.ps1 -BuildType Release
.\tests\gfn1_fast_4_3_0\run_windows_molecular.ps1 -BuildType Release
.\tests\gfn1_fast_4_3_0\run_windows_regression.ps1 -BuildType Release
# Repeat these with -BuildType Debug for runtime checks.
```

Build logs remain under the current build directories. `RESULTS_20261005.json`
and `strict_*_20261005.log` preserve this stage's dated evidence;
`SOURCE_MANIFEST_20261005.json` records source and executable hashes. Mutable
build folders must not be treated as fixed versioned releases.

## Large-molecule preliminary check

Taxol (113 atoms / 350 AO) completed analytic dispatch with a predicted
workspace of 872,726,600 bytes. Its complete projected 339×339 matrix agrees
with the official serial numerical Hessian to `2.991e-7 Eh/bohr²` at step
5e-4 bohr and SCC accuracy 1e-7. Reported Energy differs by about 1e-12 Eh;
Gradient norm agrees at printed precision. Evidence and matrix/executable
hashes are in `LARGE_PROBE_20261005.json`.

The initial stock/analytic process times were 250.94 / 31.43 seconds, but
profiling differed and these are single runs. They are **not** a repeated
performance result. `../performance/run_windows_large_hessian.ps1` performs
one warmup and five randomized profile-off trials with full-matrix comparisons
and a separate dispatch audit; its timings must be inspected before adopting
a policy for large molecules.

## Remaining scope

Analytic Hessians remain opt-in. Solvation, external fields/charges, constraints,
polarizability derivatives and periodic models retain numerical fallback.
Broader large-molecule response/memory tests, screened/nodal integral consistency
and solvent analytic response remain project work. The performance evidence in
`BENCHMARK.md` concerns small gas-phase complexes; it proves no solvent or
large-molecule analytic speedup.
