# 4.1.0 molecular geometry derivative development

Verified on Windows native, 2026-10-04, with ifx 2025.2 and oneMKL sequential.
This is an uncommitted development stage, not a complete analytic molecular
Hessian or a released package. The production CLI still uses its existing
numerical Hessian; the new geometry APIs are currently exercised by tests.

## Implemented geometry terms

- `get_hess_overlap`: contracted Gaussian overlap value, first and second
  center derivatives for s/p/d Cartesian shells. Spherical transformations
  use the existing `dtrf2`. Primitive and distance screening are preserved.
- `getCoordinationNumberHessian`: molecular second coordinate derivatives
  for exp/erf/cov/gfn CN counting functions. GFN1 uses the exp variant.
- `buildGFN1IntegralResponse`: actual GFN1 H0/overlap first derivatives for
  all coordinates, including linear CN-dependent shell self energies, and
  the fixed-density geometric Hessian `Tr(P H0_xy) - Tr(W S_xy)`.
  H0 inputs are native eV; derivative outputs and W use Eh.
- `TKlopmanOhno%getMolecularHessianResponse`: fixed-shell-charge `d(J q)/dR`,
  Coulomb energy Gradient and its geometric Hessian. It uses the existing
  shell hardnesses and averaging function with the actual generalized
  Klopman-Ohno exponent. Geometry-independent third-order terms belong to
  the separate charge-response kernel implemented in 4.0.0.

Second integral/Coulomb derivatives are directly contracted into 3N x 3N
Hessians. No NAO/NSH squared times (3N) squared second-derivative tensor is
stored. First H0/overlap derivative tensors are still dense; this is not
yet a low-memory production Hessian implementation.

Unsupported periodic cases return `ok=.false.`. Nonfinite geometry/charges,
inconsistent shapes, nonsymmetric test densities and coincident Coulomb
centers are rejected by tested guards. Callers must ignore outputs on failure.
The future molecular assembler must use the numerical fallback on these paths.

## Verification

The new modules and test programs were independently compiled using
`/Od /check:all /fpe:0`, against both freshly rebuilt Release and Debug
libraries. Debug dependencies retain their runtime checks. All final runs pass.

- 27 overlap cases: every s/p/d pair at coincident, near and separated centers.
  Existing values/Gradients agree to below 2e-13; second derivatives agree
  with Gradient finite differences to at most 2.01e-10. Both-center signs,
  spherical transforms and distance screening are checked.
- All four CN counting functions, all coordinates of an O/H/H/Si molecule:
  second derivatives agree with existing Gradient differences to at most
  1.87e-9. Symmetry, translational identities and rejection paths pass.
- Actual GFN1 water (8 AO) and disilane (30 AO), in a general orientation:
  every coordinate's dH0/dS agrees with finite differences of the existing
  H0/overlap builder. Independent arbitrary symmetric P/W matrices are held
  fixed to isolate geometric terms. Their Gradient contractions agree with
  the existing production Gradient to below 1.4e-17.
- All entries of their fixed-density Hessians are checked against finite
  differences of the existing Gradient, at 4e-4 and 2e-4 bohr steps.
  The smaller-step maximum errors are about 2.08e-10 (water) and
  9.04e-11 (disilane), and decrease by approximately four when the step halves.
- A planar water configuration exercises strictly zero overlap entries.
  H0/overlap first derivatives are checked against the original builder.
  The fixed-density Hessian is independently checked by scalar
  `Tr(P H0) - Tr(W S)` mixed energy differences, at 1e-3 and 5e-4 bohr steps.
  Smaller-step errors are 2.85e-9 in Release and 3.21e-9 in Debug.
- Actual GFN1 shell-hardness Coulomb models for water (6 shells) and Si2H2
  (10 shells), with the GFN1 exponent 2 and additional exponents 1.3/3:
  Gradient agreement is below 1.8e-18. Every-coordinate `d(Jq)` is compared
  with original matrix differences; Hessians are compared with original
  Coulomb Gradient differences. Smaller-step maxima are 1.99e-10 for
  potential derivatives and 6.14e-11 for Hessians. Quadratic convergence,
  symmetry and three translational null directions pass.

No post-hoc symmetrization is applied to analytic test Hessians. Full molecular
rotational identities require electronic response and all energy terms and
have not been established by these fixed-density/fixed-charge tests.

## Nodal overlap limitation found in the existing Gradient

The existing packed Gradient reconstructs an H0 prefactor from `H0/S` and
sets it to zero when the supplied overlap is exactly zero. An arbitrary
test density with odd AO contributions then gives a fixed-density Gradient
difference of about 7.08e-3 Eh/bohr at planar water. This diagnostic is not
an observed error in the self-consistent water calculation; such a density
need not respect the molecule's symmetry.

The new H0 coordinate API forms the prefactor directly from shell parameters
and handles nodal overlaps without division. The scalar-energy difference
test above provides an independent check in this case. The nodal test does
not claim agreement with the legacy H0/S Gradient. Before integrating the
complete Hessian, its Gradient reference and screened-overlap treatment must
be resolved consistently. Ordinary CLI screening remains unchanged.

## Native regression and reproduction

The full current tree also passes fresh Release/Debug compilation and the
ordinary source/math/dispatcher tests. Release's seven water/taxol gas/ALPB
Energy/Gradient cases and Debug's two water cases pass. Taxol comparisons
against legacy density/generic Gradient retain zero Energy difference at
printed precision and at most 2.01e-15 Eh/bohr Gradient difference.

From the repository root:

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Release
.\tests\gfn1_fast_4_1_0\run_windows_geometry.ps1
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Debug
.\tests\gfn1_fast_4_1_0\run_windows_geometry.ps1 -LibraryBuildType Debug
```

Final logs are in each build's `geometry-math-test` directory:
`geometry-test.log`, `integral-response-test.log`, `coulomb-geometry-test.log`.
An earlier Debug integral run was stopped because temporary-array warning
backtraces dominated execution. The original H0 builder's constant origin
argument is now a named array constant with identical values; the final
Debug run completed. The stopped run is not counted as a pass.

## Remaining work

The subsequent molecular assembly now connects converged MO/eigenvalue/
occupation and SCC charge references to the electronic/coupled and geometry
responses. It includes the overlap-dependent SCC potential, charge response,
repulsion and pairwise D3. Strict Release/Debug molecular tests compare the
raw unprojected matrix against all coordinates of the converged Gradient at
steps 5e-4 and 2.5e-4 bohr for water at 0/30000 K, open-shell water at 3000 K,
and disilane at 300 K. Smaller-step maximum full-Hessian errors are 4.11e-8,
3.92e-8, 4.17e-8 and 3.72e-8 Eh/bohr², decreasing approximately quadratically.
Raw reciprocity, translation and production Gradient agreement also pass.
Reproduce with `run_windows_molecular.ps1` (`-LibraryBuildType Debug` for Debug).

4.2.0 connects that gas-phase API to an opt-in CLI dispatcher, verifies actual
Hessian commands, repairs Windows numerical fallback and compares speed with
the official serial numerical reference; see `../gfn1_fast_4_2_0/RESULTS.md`.
Halogen and solvent/Born response, screened/nodal-overlap consistency, broad
model/size validation and versioned release artifacts remain open. Neither
the component tests nor the small-molecule speedups complete the full scope.
