# GFN1-fast 4.0.0 electronic Hessian response building blocks

Implemented and tested on Windows native, 2026-10-04. This stage supplies
electronic response and a coupled shell-charge solve. It does **not** yet
replace the molecular numerical Hessian.

## Implemented

`src/xtb/response.f90` now provides:

- `TGFN1ElectronicResponse`: cached first-order P/W response for a generalized
  eigenproblem, separate fixed alpha/beta particle counts, finite-temperature
  occupation response and moving AO overlap. Energies and dH use Eh; T uses K.
- Stable finite-temperature divided differences at/near degeneracy; no
  denominator shift. Zero-temperature occupied/virtual degeneracy rejects.
- `TGFN1CoupledResponse`: builds/factors `I-chi` once and reuses the shell-charge
  response factorization for coordinates. Includes the isotropic GFN1
  Hamiltonian shift `-S_ij*(v_i+v_j)/2` and validates the solved charge response.
  Singular/ill-conditioned SCC response returns failure, without regularization.
- Existing closed-shell response API remains compatible.

`TxTBCoulomb%getResponseKernel` supplies the fixed-geometry shell-potential
Jacobian from the actual lower-triangle Coulomb matrix plus atomic/shell
third-order terms. Solvent and geometry-dependent potential derivatives must
be added by the future molecular assembler.

The canonical chemical-potential response follows the particle conservation
principle used in [canonical density matrix perturbation theory](https://arxiv.org/abs/1503.07037).
This implementation evaluates the response in the full reference MO basis;
it does not claim the paper's sparse scaling or reuse its recursive algorithm.

## Verification

The new response code itself is compiled with `/Od /check:all /fpe:0`, linked
against native wrappers. Independent displaced generalized eigensolves and
canonical occupation bisection supply finite-difference references.

- Closed/open shell: 0/300/3000/30000 K, non-diagonal SPD overlap.
- AO electron conservation and differentiated HP=SW identity.
- Closed-shell zero-temperature parity against the previous response routine.
- Two displacement sizes: quadratic finite-difference convergence where
  truncation error is resolvable. At the smaller step, maximum dP error is
  approximately 2.85e-9 and dW error 7e-12 across these cases.
- Shell Mulliken charge response finite difference and charge conservation.
- Finite-temperature exact/near degeneracy, including orbital gauge invariance.
- Nonfinite/nonsymmetric perturbations, incompatible shapes/occupations,
  missing reference, failed-cache invalidation and zero-temperature degeneracy.
- Coupled response vs an independently converged nonlinear SCC model, closed
  and open shell at 0/3000 K. Maximum dP error approximately 3.46e-9,
  dW error 8.3e-12. These are synthetic SCC models, not molecular benchmarks.
- Electronic free-energy gradient and Hessian from gradient differentiation,
  plus mixed Hessian finite differences and unsymmetrized reciprocity.
- A deliberately critical attractive potential rejects the singular SCC
  linearization.
- Actual Coulomb `addShift` finite differences, including third-order terms,
  agree with the new kernel to less than 9e-14.

## Reproduce

First run the native 3.0.0 build. Then:

```powershell
.\tests\gfn1_fast_4_0_0\run_windows_response.ps1
.\tests\gfn1_fast_4_0_0\run_windows_response.ps1 -LibraryBuildType Debug
```

The runner builds the response module and Coulomb module independently with
strict checks, using the selected Release/Debug archive for dependencies.
Executables are in that native build directory's `response-math-test` folder.

## Remaining molecular Hessian work

Connect converged molecular references and first/second geometry derivatives;
include H0/overlap integral, CN, Coulomb/Born, repulsion, dispersion and
halogen terms, and assemble the complete Cartesian Hessian by differentiating
the analytic gradient. Preserve the numerical Hessian fallback and validate
against molecular gradient differences, symmetry and translational/rotational
identities. A passing electronic kernel does not establish these requirements.
