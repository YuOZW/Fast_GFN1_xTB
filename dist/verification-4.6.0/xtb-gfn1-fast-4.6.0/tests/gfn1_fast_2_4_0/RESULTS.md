# 2.4.0 Gradient implementation and native verification

The molecular GFN1 overlap-only gradient now optionally uses the existing
packed H0 directly. H0/S is evaluated only for nonzero screened overlaps,
matching matlist2; screened pairs still retain their density and
energy-weighted-density derivative terms. No dense H0/S matrix is allocated
on this path. Generic gradients and periodic systems retain their existing
paths.

For electronic derivative contractions, the spherical coefficient matrix is
mapped back to Cartesian once per shell pair. Primitive derivatives are
contracted directly with these coefficients. This replaces three derivative
transformations and a stored derivative-matrix accumulation. No tolerance
or integral cutoff is relaxed. Both improvements are enabled by default.

## Checks on 2026-10-04

Windows ifx 2025.2, sequential LP64 oneMKL, one thread unless specified:

- The seven baseline Release Energy/Gradient smoke checks pass.
- The integral module itself is compiled with `/Od /check:all /fpe:0` for
  the independent math tests. All nine s/p/d shell-pair combinations match
  the old transform-then-contract path within 1e-12. Their contracted
  derivatives also match finite differences of overlap integrals.
- Disilane (30 AO) and taxol (350 AO) pass old/new, packed-only,
  contraction-only, generic-gradient and four-OpenMP-thread comparisons.
- ALPB-water and +1/one-unpaired-electron cases pass for both molecules.
  All total energies agree at printed precision; max gradient differences
  are below 1e-10 Eh/bohr (largest observed 7.9e-14).
- Full Energy finite differences for two disilane coordinates agree with
  the reported Gradient to 3.85e-9 and 1.19e-8 Eh/bohr.
- The full Debug build passes water gas/ALPB, disilane old/new/fallback,
  four-thread, ALPB/open-shell and full Energy finite-difference checks.
  The Debug math executable uses the matching debug runtime libraries.
- Five paired process wall measurements per mode on taxol: profile off,
  old median 0.223360 s / new 0.220483 s (-1.3%); profile on, old 0.221700 s /
  new 0.220050 s (-0.7%). The overlap-gradient kernel median is 0.036 s /
  0.035 s. These small differences require broader benchmarks before
  claiming a general speed improvement; times include startup/file output.

The original 594-AO TMS input is unavailable; no comparison to its WSL
timings is claimed. This is a development increment, not a packaged release.

## Reproduce / fallback

```powershell
.\tests\gfn1_fast_2_4_0\run_windows_regression.ps1
.\tests\gfn1_fast_2_4_0\run_windows_regression.ps1 -BuildType Debug
```

Debug tests use water and disilane; the much more expensive taxol bounds-
checked gradient is covered by Release comparisons instead. Current output
is in `build-gfn1-fast-2.4.0-windows-ifx-release/gradient-regression/summary.json`
and per-case logs. Debug outputs use the corresponding `-debug` directory.

Set `XTB_GFN1_FAST_DISABLE_PACKED_GRADIENT=1` to retain the dense H0/S
construction, or `XTB_GFN1_FAST_DISABLE_GRADIENT_CONTRACTION=1` to retain
derivative matrices and the three spherical transformations. Set both to
recover the 2.2.6 overlap-only gradient. The original
`XTB_GFN1_FAST_DISABLE_GRADIENT_KERNEL=1` still selects the generic gradient.
