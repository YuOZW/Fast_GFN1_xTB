# 2.3.2 subspace reuse: implemented, production speed target not met

The source includes exact GFN1 shell-shift projection into the preceding full
MO basis, an invariant-subspace graph correction, strict generalized AO
residual validation, same-H full-spectrum audits, and final full verification.
The independent standard symmetric solver also contains block Rayleigh--Ritz
corrections and thick restart. The GFN1 SCC dispatcher uses the graph path;
failed graph iterations return directly to the full factorized solver.

The first eligible candidate is audited; subsequently every third eligible
step is audited. Audits compare eigenvalues, density, shell populations and
band free energy at the same H and always continue with the full solution.
SCC mixing and convergence thresholds are unchanged. Reference caches are
local to an SCC invocation and only created near convergence.

This is an opt-in experimental feature. Normal runs keep SYEVD because the
350-AO taxol benchmark does not benefit. The original 594-AO input is absent,
and the intended reduction to four to six full solves is **not achieved**.
This result does not establish performance on larger systems.

## Native checks on 2026-10-04

Intel ifx 2025.2 / LP64 sequential oneMKL, one thread, Windows x64:

- Release build and all seven baseline Energy/Gradient smoke cases pass.
- The new module itself is compiled with `/Od /check:all /fpe:0` for math
  tests. Standard symmetric and non-diagonal SPD generalized problems agree
  with full SYEVD; AO residual and occupied density are independently checked.
- Missing reference, mismatched shapes, nonfinite shift and zero iteration
  limit are rejected. Bounded block correction/thick restart is exercised.
  Gershgorin block bounds reject a disconnected virtual intruder, which an
  eigenvector residual alone cannot detect.
- Full/reuse comparisons pass at 0 K, 300 K, 1000 K, with ALPB water, and for
  the +1/one-unpaired-electron taxol case. Energy differences are zero at
  printed precision; max gradient component differences are below 1e-10
  Eh/bohr. At 0/300 K, two SCC steps use reuse and one is fully audited;
  the other cases safely fall back. Tests do not count fallback as accepted
  reuse.
- Iteration-limit fallback and minimum-size gating are exercised.
- Five paired process wall measurements per mode: with profile off, full
  median 0.219624 s versus forced reuse 0.266783 s (+21.5%); with profile on,
  full 0.222320 s versus reuse 0.267806 s (+20.5%). Times include process
  startup and file output. These measurements precede the final cache-timing
  instrumentation; the runner writes current results to JSON.

## Reproduce

```powershell
.\tests\gfn1_fast_2_3_2\run_windows_regression.ps1
```

The native setup and pinned dependency are documented in
`tests/gfn1_fast_2_2_6/WINDOWS_TEST_REPORT.md`. Outputs are in
`build-gfn1-fast-2.3.2-windows-ifx-release`, with machine-readable
`reuse-regression/summary.json` and per-case logs. Build directories contain
the most recently built source, not immutable version releases.

## Controls

Enable with `XTB_GFN1_FAST_ENABLE_SUBSPACE_REUSE=1`. Disable always takes
precedence: `XTB_GFN1_FAST_DISABLE_SUBSPACE_REUSE=1`.

The measured cost gate is enabled by default. Set
`XTB_GFN1_FAST_DISABLE_REUSE_AUTOTUNE=1` only for forced comparisons.
`XTB_GFN1_FAST_REUSE_MIN_NAO`, `REUSE_GUARD`, `REUSE_MAX_ITER`,
`REUSE_AUDIT_PERIOD` and `REUSE_MAX_CHARGE_RMS` (each with the same
`XTB_GFN1_FAST_` prefix) adjust the probe policy, not physical precision.
Defaults are 128, 16, 8, 3 and 1e-3. Residual and parity tolerances cannot be
relaxed by these controls.
