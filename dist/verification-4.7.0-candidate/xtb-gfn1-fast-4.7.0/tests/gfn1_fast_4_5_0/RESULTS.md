# 4.5.0 complete Born-solvent SCC Hessian

The Windows ifx Release and strict Debug builds pass the complete solvent
response tests (2026-10-05). The opt-in calculator/CLI analytic Hessian now
supports GFN1 GBSA/ALPB with Still/P16 kernels and no salt. Energy/Gradient
calculations retain their existing production paths. Analytic Hessians remain
disabled by default pending wider performance, resource and adoption checks.

## Implementation and independent verification

`src/solv/derivative.f90` supplies small scalar first/second derivative algebra
for Born kernels and the six independent inertia entries. The complete solvent
functional is Q=q+CM5(R), E=Q^T A(R) Q/2 + SASA + constant. The response returns
the fixed-bare-q Gradient/Hessian, charge potential A Q, charge kernel A and
fixed-q coordinate derivative A_x Q + A CM5_x. It includes self/pair Born,
H-bond surface terms, screened angular SASA and the ALPB shape coefficient.
The assembly retains original cached-value/first-derivative/matrix parity
guards. It streams pair coefficients without a molecular-size A_xx tensor.

Mapped solvent potentials, kernels and coordinate derivatives enter the
existing shell-charge linear response before its Hamiltonian residual and
coupled solve. The full Hessian includes the existing repulsion, D3 and
halogen terms. Both the geometric CM5 Hessian and all mixed charge/coordinate
terms are retained; a gas-only response at a solvated reference rejects.

The native test programs are compiled with `/Od /check:all /fpe:0`; they link
either the optimized Release or strict Debug library. Independent references
use the original TBorn Energy/Gradient and nonlinear converged SCC calculations.
No Hessian symmetrization/projection is used in native comparisons.

- Eight fixed-bare-charge solvent cases: GBSA/ALPB, Still/P16, neutral and
  charged. Original Energy, Gradient, potential and Born matrix agree. All nine
  Cartesian Hessian/potential columns use three step sizes down to 2.5e-5 bohr.
  Charge finite differences independently check A, potential and the mixed
  Maxwell identity. Raw symmetry and translation pass.
- 38 real molecular SCC cases: 32 asymmetric-water conditions covering all
  model/kernel combinations, 0/300/1000/30000 K and neutral/open-shell cation;
  four disilane cases with d AOs; two CH3Br-water cases with a nonzero halogen
  correction. Every Cartesian column is compared at 2e-4 and 1e-4 bohr.
  Smaller-step maximum raw Hessian error is 1.172e-8 Eh/bohr² across both builds.
  Raw symmetry, translation, production Gradient and quadratic step convergence
  pass. Partial additive solvent calculator Hessian/dipole columns also pass.
- Actual CLI: ten analytic cases plus five numerical-fallback cases in both
  configurations. Full output matrices agree with numerical Hessians; larger
  fixtures at step 1e-4 differ by at most 6.7e-9 Eh/bohr². Boundary, conservative
  displacement window, memory, legacy Gradient and salt fallback outputs are
  exactly the corresponding numerical outputs.
- Ordinary water/taxol Energy/Gradient, gas CLI Hessians, nonzero halogen
  molecular response and nodal Energy/Gradient/fallback checks were rerun.
  Standard water/taxol values and path parity are retained.

## Diagnosed pre-existing behavior

The charged GBSA/Still water at 1000 K has a raw Energy finite-difference
discrepancy about 3.15e-8 Eh/bohr. A separate program using only established
calculator APIs reproduces it with stock xTB. At 2.5e-5 bohr, stock/current
coordinate 3 give discrepancies about 3.17e-8/3.19e-8 and dN/dR about 4.10e-8.
`fermismear` stops at its original 1e-9 electron-number tolerance; the computed
Energy derivative includes mu*dN/dR. The canonical Energy derivative after
accounting for independently measured alpha/beta dN is within 1.194e-8 over
all molecular test steps and both builds. This is a test-side diagnostic;
production occupations, tolerances and Energy arithmetic were not changed.
Stock's separate known nodal Gradient defect appears in coordinates 2/5;
it is not misidentified as the thermal discrepancy.

The original CH3Br-water surface also crosses an existing screened SASA branch
over a 2e-4-bohr displacement window. The new guard rejects it. A separate
bounded geometry probe selects trial 5 solely using original surface branch
margins, before any SCC Hessian comparison; that frozen geometry is used for
the smooth halogen test and benchmark. The original boundary remains tested.
The disilane ALPB guard is conservative: it rejects step 2e-4 despite the raw
finite-difference check passing, and accepts 1e-4. These branch decisions and
fallback costs are retained, not overridden to force a faster measurement.

An initial fused shape-moment loop produced incorrect ALPB moments under
ifx `/O3`, while Debug was correct. Separating explicit center/moment loops
from derivative construction fixes the Release discrepancy. Original ALPB
factor/matrix parity and independent original Gradient comparisons pass in
both builds, with temporary diagnostics removed.

## Resource and dispatch limits

The calculator includes an additional conservative bound
8*(27*N³ + 200*N² + 4000*N) bytes for solvent tensors/caches. The existing
1024-MiB default cap remains. This is an estimate, not a measured peak.
Solvent tensors are still dense and scale cubically; large solvent Hessian
performance/memory is not established by these small-molecule tests.
COSMO/CPCM-X, salt, unsupported kernels, external potentials, constraints,
polarizability requests, failed response/branch/resource checks continue to
use numerical fallback. The existing gas restrictions are preserved.

The supplied Hessian finite-difference step defines the conservative surface
branch window. Ordinary CLI default steps can therefore fall back where the
smaller validated step uses analytic response. No universal solvent speedup
or default adoption is claimed.

## Reproduce

From the repository root, after the standard Windows builds:

```powershell
.\tests\gfn1_fast_4_5_0\run_windows_solvent.ps1 -BuildType Release
.\tests\gfn1_fast_4_5_0\run_windows_molecular.ps1 -BuildType Release
.\tests\gfn1_fast_4_5_0\run_windows_regression.ps1 -BuildType Release
# Repeat with -BuildType Debug.
.\tests\gfn1_fast_4_5_0\run_windows_thermal.ps1
.\tests\gfn1_fast_4_5_0\run_windows_thermal.ps1 -Reference
.\tests\gfn1_fast_4_5_0\run_windows_regression.ps1 -BuildType Release -Benchmark
```

The native molecular runner explicitly enables analytic Hessians for the
partial-column calculator test. Native/CLI logs, measured hashes and copied
JSON reports are recorded in `SOLVENT_RESULTS_20261005.json` and
`SOURCE_SOLVENT_MANIFEST_20261005.json`. The earlier component manifest and
component snapshot are historical evidence, preserved separately.

Remaining full-goal work includes larger solvent response/performance and
memory checks, thread crossover/default adoption, general screening consistency
and a validated distributable source/archive workflow.
