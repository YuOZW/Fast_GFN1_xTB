# Current working-tree Windows verification

See `RESULTS.md` for the subsequent formal Release/Debug FOE regressions,
paired five-run benchmarks and high-temperature electron-count diagnostics.

Verified 2026-10-04, Windows native; WSL was not used. This report covers the
current uncommitted working tree, including the 2.3.2, 2.4.0 and experimental
3.0.0 changes. Git HEAD remains the 2.2.6 baseline. Build folder names identify
test configurations, not immutable release packages.

## Toolchain

Intel ifx 2025.2.0, Visual Studio 2022 MSVC, oneMKL LP64 sequential,
CMake/Ninja, and the AGENTS.md Python executable:
`C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe` (3.13.15).
The pinned mctc-lib dependency is already present.

Windows compatibility settings and initial baseline verification are documented
in `../gfn1_fast_2_2_6/WINDOWS_TEST_REPORT.md`.

## Build and ordinary execution

Both Release and Debug compile and link successfully. Release passes seven
Energy/Gradient cases: water gas/ALPB, taxol production/legacy density/generic
gradient, and taxol ALPB production/legacy density. Debug passes water gas and
ALPB with bounds/allocation checking and floating-point traps enabled.

Both builds also pass the existing source guards, compact-density mathematical
test and compiled dispatcher test.

| Case | NAO | Energy / Eh | Gradient norm / Eh per bohr |
|---|---:|---:|---:|
| Water gas | 8 | -5.764209595203 | 0.094423990378 |
| Water ALPB water | 8 | -5.787247944321 | 0.075935279788 |
| Taxol gas | 350 | -195.911950900123 | 0.149903893916 |
| Taxol ALPB water | 350 | -195.958435141829 | 0.136167479182 |

Taxol converges in 11 SCC iterations, using compact density for 10 iterations
in the production configuration. Energy differences against the legacy density
and generic gradient runs are zero at printed precision; the maximum gradient
component difference is 2.01e-15 Eh/bohr. Acceptance limits are 1e-11 Eh and
1e-10 Eh/bohr.

## Experimental density-only SCC

The independent Fermi-operator module was compiled with `/Od /check:all /fpe:0`
and tested against known-spectrum density and energy-weighted density references.
Zero temperature, 50/300/1000 K, closed/open shell, clipped entropy, electron
count, HP=SW, bounded-work rejection, degeneracy, nonfinite Hamiltonian and
indefinite overlap tests pass. Maximum P/W errors are approximately 2e-15.

Release integration probes pass Energy/Gradient comparisons for disilane and
taxol at 0/300/1000/5000 K. With experimental FOE forced and its performance
autotuner disabled:

| Case | Density-only accepted iterations | Full audits | Audit failures |
|---|---:|---:|---:|
| Disilane, each tested temperature | 2 | 1 | 0 |
| Taxol, 0/300 K | 3 | 2 | 0 |
| Taxol, 1000/5000 K | 0 | 1 | 1 |

Taxol's two high-temperature candidates fail strict full-solve audits and are
discarded. Their final results match the full-solve reference. This is fallback
coverage, not evidence that FOE is accurate enough for those cases. The largest
final gradient difference across the probes is 9.01e-14 Eh/bohr; printed energies
agree. Dense FOE is slower on taxol and remains disabled by default.

Debug disilane integration also passes at 0 and 5000 K, with two density-only
iterations and one successful full audit in each case. Maximum gradient
differences against full solve are 9.99e-16 and 1.01e-14 Eh/bohr respectively.
Bounds/allocation checks and floating-point traps remain enabled. Results are recorded in
`build-gfn1-fast-3.0.0-windows-ifx-debug/foe-debug/summary.json`.

## Reproduce

From the repository root in PowerShell:

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-3.0.0-windows-ifx -BuildType Release
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-3.0.0-windows-ifx -BuildType Debug
```

The scripts initialize oneAPI for their own process, build and execute the
tests, then restore the calling environment. Executables are in the respective
`build-gfn1-fast-3.0.0-windows-ifx-release` and `-debug` directories; logs and
machine-readable results are in `build.log` and `smoke/summary.json`.

For a direct Release run, use a oneAPI environment and a suitable input:

```powershell
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_STACKSIZE = '64M'
cmd /d /c 'call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 >nul && "C:\GitHub_repository\Fast_GFN1_xTB\build-gfn1-fast-3.0.0-windows-ifx-release\xtb.exe" input.xyz --gfn 1 --grad --norestart'
```

The old 594-AO TMS coordinate input is absent. Cross-WSL parity and performance
have not been verified. These checks establish native build/run compatibility
and the tested numerical comparisons; they are not release packaging checks.
