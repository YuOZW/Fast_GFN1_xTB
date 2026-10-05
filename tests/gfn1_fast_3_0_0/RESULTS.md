# GFN1-fast 3.0.0 experimental diagonalization-free SCC

Implemented direct P/W construction for intermediate SCC iterations using
zero-temperature SP2 purification and a finite-temperature recursive rational
Fermi operator. The final verification and orbitals use full diagonalization.
Electronic number, HP=SW, matrix refinement and periodic same-H full audits
gate acceptance. Charge mixing and convergence conditions are unchanged.

The implementation is dense and **disabled by default**. This is a verified
algorithm prototype, not a production speed improvement.

## Regressions completed on Windows, 2026-10-04

Release and Debug pass independent known-spectrum P/W, entropy and electron
count tests, including open shell, degeneracy and bounded-work rejection.
Runtime comparisons use energy tolerance 1e-11 Eh and maximum gradient
component tolerance 1e-10 Eh/bohr, and assert actual accepted/audited paths.

- Disilane at 0/300/1000/5000 K: 2 density-only iterations, 1 full audit,
  no failures; full solves 9 -> 7.
- Disilane ALPB and open shell: density-only steps and successful audits.
- Taxol at 0/300 K and ALPB: 3 density-only iterations, 2 full audits,
  no failures; full solves 11 -> 8.
- Taxol 1000/5000 K and open shell: audit rejects candidate; full fallback
  preserves final Energy/Gradient. These cases do not establish accepted FOE.
- Maximum-work, minimum-size, explicit disable precedence and warmup gates pass.
- Autotune execution preserves final Energy/Gradient; it can stop further FOE
  work after a slow trial. Its initial trial still has measurable cost.
- Profile-on and profile-off repeated paired runs also check Energy/Gradient.

Debug runtime regression covers disilane only; large-system Debug execution
is not included. Release covers the taxol cases.

## High-temperature audit diagnosis

| Taxol case | FOE electron count error | Full-reference electron count error |
|---|---:|---:|
| 1000 K | 2.27e-13 | 1.21e-9 |
| 5000 K | 6.25e-13 | 5.02e-10 |

The full reference uses fermismear's 1e-9 stopping tolerance per spin, and
returns occupations evaluated before its last chemical-potential correction.
FOE enforces a tighter canonical count. These observed number differences
explain why a strictly canonical candidate may differ from the legacy
reference even when HP=SW holds. The measurements establish a source of
disagreement, not a decomposition proving that it accounts for all P/W/energy
error. Audit gates were not relaxed, and fermismear was not changed.

## Performance: five paired runs per profile mode

End-to-end process wall times, one OpenMP and one MKL thread, Windows ifx/MKL.
Medians on the 350-AO taxol input; startup and I/O are included.

| Profile | Full / s | Forced FOE / s | FOE / full | FOE with autotune / s |
|---|---:|---:|---:|---:|
| Off | 0.220003 | 0.753812 | 3.426 | 0.343082 |
| On | 0.220311 | 0.760821 | 3.453 | 0.346812 |

Dense FOE is slower for this case. Sparse/local algorithms or a different
temperature/size range would need new measurements before changing the default.

## Reproduce

```powershell
.\tests\gfn1_fast_3_0_0\run_windows_regression.ps1
.\tests\gfn1_fast_3_0_0\run_windows_regression.ps1 -BuildType Debug
```

For explicit runtime experimentation set `XTB_GFN1_FAST_ENABLE_FERMI_OPERATOR=1`.
`XTB_GFN1_FAST_DISABLE_FERMI_OPERATOR=1` takes precedence. Other controls are
`FERMI_OPERATOR_MIN_NAO` (128), `FERMI_OPERATOR_WARMUP` (4),
`FERMI_OPERATOR_MAX_STEPS` (100), `FERMI_OPERATOR_AUDIT_PERIOD` (3), and
`FERMI_OPERATOR_MAX_CHARGE_RMS` (1e-2), all prefixed with `XTB_GFN1_FAST_`.
The test-only forced mode sets min NAO=1 and disables autotune.

Authoritative result files are `foe-regression/summary.json` and the per-case
logs in the respective native build directory. `probe_runtime.py` remains a
development diagnostic; `test_runtime.py` is the regression entry point.
The earlier 594-AO benchmark input is absent. Release ZIP/SHA/re-extraction
checks and the full requested project, including a molecular analytic Hessian,
are still outstanding.
