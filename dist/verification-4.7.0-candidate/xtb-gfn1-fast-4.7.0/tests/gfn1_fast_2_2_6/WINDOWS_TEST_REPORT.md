# Windows native build and execution check

Verified on 2026-10-04 (JST), starting from GFN1-fast 2.2.6, commit
`c2ccf660ffe2ddc5f96971b77fdb325c82a8fbeb`, plus the Windows compatibility
changes described below. WSL was not used.

## Toolchain and configuration

- Windows x64 / PowerShell
- Intel Fortran `ifx` 2025.2.0 (20250605)
- Visual Studio 2022 MSVC 19.44.35217.0
- Intel oneMKL 2025.2, LP64 sequential BLAS/LAPACK
- CMake 4.1.1 / Ninja 1.12.1
- Python 3.13.15: `C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe`
- mctc-lib commit `77f65c6f2cf6330d05d0757ca173da097096780e`
- Native GFN1: tblite, CPCM-X, upstream unit-test targets and JSON disabled
- Release: `/O3 /QaxAVX2 /real-size:64 /traceback`, without fast-math or LTO
- Debug: `/Od /debug:full /check:all /fpe:0`, plus the same dialect flags
- PE x64 executable, 64 MiB main-thread stack reserve
- Runtime checks use one BLAS/OpenMP thread and 64 MiB OpenMP stacks

## Compatibility changes

1. `cmake/CMakeLists.txt` now selects Windows Intel option spelling. Previously
   `-axAVX2` and `-r8` were silently ignored by Windows ifx.
2. The PowerShell build uses `WITH_OBJECT=OFF` and builds `xtb-exe`. The default
   object/shared build produced two rules for `xtb.lib` on Windows, preventing
   Ninja from starting.
3. The executable linker reserves 64 MiB of stack. Windows' default stack caused
   taxol's Coulomb-gradient OpenMP reduction to stop with error 170.
4. `src/type/environment.f90` uses `USERPROFILE` when `HOME` is unavailable and
   initializes an empty fallback if neither exists. Debug previously stopped
   before CLI startup with an unallocated `HOME` error. The integration tests
   explicitly remove `HOME` and `XTBHOME` from the child process environment.
5. The scripts handle launchers containing both `PATH` and `Path`, and supply
   the Intel/MKL DLL directories to Python's child executable environment.

No GFN1 model, SCC mixing, convergence threshold or solver policy was changed.

## Results

Both Release and Debug compilation: **PASS**.

Both configurations pass the seven existing source guard groups, the
compact-density mathematical test and the compiled Fortran dispatcher test.

Release passes seven native Energy/Gradient runs:

| Case | NAO | Energy / Eh | Gradient norm / Eh per bohr |
|---|---:|---:|---:|
| Water, gas | 8 | -5.764209595203 | 0.094423990378 |
| Water, ALPB water | 8 | -5.787247944321 | 0.075935279788 |
| Taxol, production | 350 | -195.911950900123 | 0.149903893916 |
| Taxol, legacy density | 350 | -195.911950900123 | 0.149903893916 |
| Taxol, generic gradient | 350 | -195.911950900123 | 0.149903893916 |
| Taxol, ALPB water | 350 | -195.958435141829 | 0.136167479182 |
| Taxol, ALPB water, legacy density | 350 | -195.958435141829 | 0.136167479182 |

Taxol converges in 11 SCC iterations; production uses compact density for 10
iterations. Disabling compact density gives zero compact iterations. The
generic-gradient run confirms that the specialized kernel is disabled.

| Release comparison | Absolute energy difference / Eh | Max gradient component difference / Eh per bohr |
|---|---:|---:|
| Production vs legacy density | 0 at printed precision | 2.00534e-15 |
| Production vs generic gradient | 0 at printed precision | 9.99201e-16 |
| ALPB production vs legacy density | 0 at printed precision | 1.09981e-15 |

Acceptance: energy difference <= 1e-11 Eh and gradient difference <= 1e-10
Eh/bohr. Every run checks SCC convergence, finite gradient components and
agreement between the written gradient and reported norm.

Debug passes water gas and water ALPB Energy/Gradient with allocation/bounds
checks and floating-point traps enabled. Its results agree with the Release
values shown above. The larger taxol Debug run was interrupted because of
runtime-checking cost; it is not counted as a completed check. The default
Debug script therefore runs the two small-system checks only.

The former 594-AO TMS benchmark input is not present here. Cross-WSL numerical
parity and performance comparisons were not performed. These are functional
checks, not repeated performance benchmarks.

## Reproduce

From `C:\GitHub_repository\Fast_GFN1_xTB`:

```powershell
& .\tests\gfn1_fast_2_2_6\run_windows_build.ps1
& .\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildType Debug
```

The pinned dependency has already been prepared in `subprojects/mctc-lib`.
For a new checkout where it is absent:

```powershell
git init subprojects/mctc-lib
git -C subprojects/mctc-lib fetch --depth 1 https://github.com/grimme-lab/mctc-lib.git 77f65c6f2cf6330d05d0757ca173da097096780e
git -C subprojects/mctc-lib checkout --detach FETCH_HEAD
```

The build script imports oneAPI into its own process and restores the caller's
environment. The executable dynamically uses Intel Fortran/OpenMP/MKL DLLs,
so direct runs should also initialize oneAPI, for example:

```powershell
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_STACKSIZE = '64M'
cmd /d /c 'call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 >nul && "C:\GitHub_repository\Fast_GFN1_xTB\build-gfn1-fast-2.2.6-windows-ifx-release\xtb.exe" input.xyz --gfn 1 --grad --norestart'
```

Outputs are in `build-gfn1-fast-2.2.6-windows-ifx-release` and
`build-gfn1-fast-2.2.6-windows-ifx-debug`: `xtb.exe`, `build.log`, the compiled
policy test, individual `smoke/<case>/run.log` and gradient files, and
`smoke/summary.json`. Full Debug smoke tests can be requested separately by
omitting `--small-only` from `test_windows_smoke.py`, with oneAPI initialized.
