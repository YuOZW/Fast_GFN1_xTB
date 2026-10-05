# 4.5.0 solvent response implementation status

The full Born-solvent response passes independent native and actual CLI tests
in Release and Debug (Windows ifx 2025.2 / oneMKL sequential, 2026-10-05).
It is connected to the opt-in calculator Hessian for GBSA/ALPB Still/P16
without salt, including all solvent-induced SCC and CM5 response terms.
The frozen 4.4 five-repeat taxol benchmark has completed. These
component builds/tests started after its terminal success, avoiding competing
compilation or runtime test work during the performance measurements.

Added components:

- `getBornRadiusResponse`: analytic Cartesian first/second Born-radius
  derivatives through every directed descreening branch and the original
  OBC transformation. Retains the near-equal reduced-radius special case,
  actual neighbor cutoff, and cached-value/first-derivative parity checks.
- `getSurfaceResponse`: exact screened angular quadrature derivatives,
  with the original surface constants/weights/neighbor order. Optional
  displacement-scale boundary checks reject soft-sphere or point-screening
  references too close to a branch change.
- `calc_cm5` optional response: streamed Hessian of a weighted sum of CM5
  corrections. Ordinary calls preserve their operations; the new contraction
  avoids storing every atom's complete CM5 Hessian simultaneously.
- Native Born/CM5 and surface finite-difference programs, with independent
  original-gradient comparisons and explicit boundary/branch cases.
- Small scalar derivative algebra for Still/P16 and ALPB inertia, complete
  fixed-bare-charge functional derivatives and mapped shell-charge coupling.
- Full molecular solvent response, conservative memory/branch dispatch,
  selected additive calculator columns, numerical fallback and repeated
  stock/fast numerical/analytic CLI measurements.

Remaining work before broader default adoption / completion of the full goal:

1. Larger solvent systems: response, peak resource behavior, conservative
   branch-window acceptance and repeated performance.
2. Numerical/analytic thread crossover, ordinary default-step and default
   policy adoption based on measurements rather than small-input extrapolation.
3. General screening consistency and validated source/distribution archive.
4. Salt and other unsupported solvent models retain numerical fallback; their
   responses require separate independent verification before any extension.

The complete solvent functional uses Q=q+CM5(R),
E=0.5 Q^T A(R) Q + E_SASA(R) + constant. Thus the shell-charge response kernel
must include the mapped A, and the fixed-q nuclear derivative of the potential
is A_x Q + A CM5_x. Geometric Hessians need both mixed A_x/CM5_y terms,
CM5_x^T A CM5_y, and the potential-weighted CM5 second derivative.
Gas-only response at a solvated reference is insufficient.

## Evidence

`COMPONENT_RESULTS_20261005.json` and `SOURCE_MANIFEST_20261005.json` preserve
the earlier component-only evidence and its executable hashes. The frozen
component snapshot retains those exact sources. Full-response current evidence
is in `SOLVENT_RESULTS_20261005.json`, `SOURCE_SOLVENT_MANIFEST_20261005.json`,
`RESULTS.md` and `BENCHMARK.md`. Dated native logs are saved separately.

- Born radii: 32 model/pair/geometries, all six Cartesian columns, four step
  sizes down to `5e-5 bohr`; smaller-step maximum difference `5.39e-9` in native
  radius-derivative units. Directed branch counts: 16 fully buried, 16 partial,
  32 nonoverlapping. The near-equal reduced-radius branch also uses deliberately
  distinct radii separated by `5e-9`, and agrees with the original derivative.
- Weighted CM5: 16 distinct pair/geometries, repeated in both solvent fixtures,
  four steps down to `2.5e-5 bohr`; maximum difference `2.70e-8` in native
  charge-derivative units. Raw symmetry, total-charge Hessian conservation,
  incomplete optional arguments and invalid output shape are checked.
- Angular surfaces: actual water GBSA/ALPB parameter fixtures, all nine columns,
  steps `1e-4`/`5e-5 bohr`; maximum smaller-step difference `1.50e-8` in native
  surface-derivative units, raw reciprocity and translation invariance pass.
  The known screened-surface fixture rejects a `2e-4 bohr` displacement window.
  A fixed quadrature is not exactly rotation invariant; no false continuum
  rotational identity is asserted for its discretized surface.
- Initial coarse finite differences at 0.4-bohr atom separation exceeded the
  absolute target. Four-step sweeps show quadratic reduction and satisfy the
  original `2e-7` final bound; no analytic formula or accuracy tolerance was
  altered to accommodate the discrepancy.
- Rebuilt usual Release smoke: seven cases; Debug: water gas/ALPB. Existing
  source/math/dispatcher checks pass. All-coordinate nodal gas/ALPB/GBSA Energy
  derivatives and actual CLI analytic/fallback Hessians pass in both builds.

The list above describes the initial component-only checks. The later full
response stage adds 38 converged molecular SCC cases, all-coordinate raw
Hessians, CLI/fallback checks and seven-repeat solvent performance comparisons.
Small tested solvent Hessians are 3.75–4.37 times faster than stock numerical;
this does not establish large-N or universal default-step speedups.

Re-run after rebuilding the current libraries:

```powershell
.\tests\gfn1_fast_4_5_0\run_windows_born.ps1 -BuildType Release
.\tests\gfn1_fast_4_5_0\run_windows_surface.ps1 -BuildType Release
# Repeat with -BuildType Debug.
```
