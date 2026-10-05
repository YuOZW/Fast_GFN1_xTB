# GFN1-fast 4.6.0: default analytic Hessian and boundary fallback

The analytic GFN1 Hessian is now the default for supported gas and Born-solvated
models. Explicit `ENABLE_ANALYTIC_HESSIAN=0` disables it; `DISABLE_ANALYTIC_HESSIAN=1`
takes precedence. Existing model, memory, solvent branch, reciprocity and
translation guards retain the numerical fallback. Subspace reuse and the
fermi-operator SCC experiments remain disabled by default.

## Reproduced defects and safeguards

The old 4.5.1 response accepted a small nonzero SCF-screened overlap in asymmetric
water with a production/response Gradient difference of 0.046801854207 Eh/bohr.
The exact node and ordinary nonzero overlap matched to about 1e-15. The new
response compares its electronic Gradient with the actual production derivative
using the raw packed H0 and screened S, including CN chain derivatives. It rejects
at differences above 1e-10 Eh/bohr; it does not change the parent screening model.
`run_windows_screening.ps1` passes in Release and Debug: two smooth references
accepted and the two screened nonvariational references rejected.

The original Ag-pair D3 Gradient jumps at its hard 60-bohr cutoff. Previously
accepted analytic Hessians differ from all-coordinate Gradient finite differences
by about 9.77e-7 Eh/bohr² there. `run_windows_cutoffs.ps1` checks the smooth sides,
the exact boundary and points whose displacement window straddles it. Release
and Debug pass. D3, CN and repulsion Hessians now reject cutoff windows; halogen
Hessians guard both cutoff and nearest-neighbour changes in the requested window.
The CLI forwards its actual Hessian displacement to gas and solvent response.

## Native and actual CLI validation

`run_windows_default.ps1` passes 21 conditions in both configurations: 11 valid
analytic defaults, nine explicit/model/boundary fallbacks and analytic thread
parity. Numerical baselines explicitly set DISABLE, while empty policy exercises
the real default. The CLI applies a small Euler rotation before Hessian evaluation;
the screened-overlap fixture undoes it in the input. The D3 CLI fixture uses two
water fragments because distant isolated Ag2 fails the parent SCC eigenproblem;
the independent component reproduction continues to use real Ag parameters.

The 38 actual solvated SCC native conditions and the existing gas/solvent CLI
regressions pass with the new source. Component response and partial additive
calculator tests also pass. Full hash-bound results and performance/distribution
validation are recorded separately as they complete. Frozen 4.5.1 benchmark
results are historical evidence and are not relabelled as 4.6 timings.

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildType Release -BuildName build-gfn1-fast-current-windows-ifx
.\tests\gfn1_fast_4_6_0\run_windows_screening.ps1
.\tests\gfn1_fast_4_6_0\run_windows_cutoffs.ps1
.\tests\gfn1_fast_4_6_0\run_windows_default.ps1
```

Use `-BuildType Debug` for the checked build. Source package verification uses
the pinned mctc-lib source manifest if the extracted package has no Git database.
