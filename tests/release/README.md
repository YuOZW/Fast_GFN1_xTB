# Verified source package workflow

The 4.6.0 source ZIP has been created, SHA-256 verified and extracted into a
new directory. All 921 payload hashes and all 124 pinned mctc-lib source files
passed. Both native Release/Debug builds completed 786 steps with no Git
metadata, then passed smoke, 21 actual default/fallback CLI conditions, the
screening and D3 cutoff-window tests and 38 solvated SCC Hessian conditions.
The authoritative immutable-package receipt is `dist/REEXTRACTION_4.6.0.json`.
The ZIP preserves its creation-time documentation; completion results are
recorded outside it so its verified hash remains stable.

`create_source_package.py` requires matching final source/binary hash evidence,
a source version matching the profile marker, and the exact clean pinned
mctc-lib commit. It writes per-file internal hashes and an external ZIP SHA.
It does not overwrite an existing output, publish, modify Git or bundle Intel
runtime/compiler binaries. `verify_source_package.py` checks extracted file
hashes and safe paths; `--dependencies-only` verifies the bundled dependency
inventory for the gitless Windows build runner.

Use the prescribed Conda Python directly:

```powershell
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
& $PY tests/release/verify_source_package.py --root <extracted-source-root>
```

Build with the extracted `tests/gfn1_fast_2_2_6/run_windows_build.ps1`, setting
`-BuildName build-gfn1-fast-current-windows-ifx` and `-BuildType Release` or
`Debug`. The example source manifest and ZIP name are version-specific; for
later releases use newly validated evidence and a fresh output path.

## Windows runtime ZIP (oneAPI is needed only on the packaging PC)

`create_windows_package.py` produces a ZIP containing the Release executable,
app-local Intel Fortran/OpenMP, sequential MKL and Microsoft VC CRT DLLs,
parameters, notices, the matching source ZIP and a per-file SHA-256 manifest.
MKL CPU dispatch DLLs are included for different CPUs, including the SSE4.2
path, rather than copying only the DLLs loaded by the packaging PC.
The PE dependency closure includes delayed imports. Compiler DLLs must appear
on the installed Intel Fortran redistribution list; VC CRT DLLs are taken
from Visual Studio's x64 REDIST directory, never from System32.

In an Intel oneAPI developer PowerShell, after the README Release build:

```powershell
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
& $PY tests/release/create_windows_package.py `
    --build-dir build/windows-release `
    --compiler-root 'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.2' `
    --mkl-root 'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.2' `
    --vc-crt-dir 'C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Redist\MSVC\14.44.35112\x64\Microsoft.VC143.CRT' `
    --output dist/fast-gfn1-xtb-4.7.0-windows-x64-portable.zip
& $PY tests/release/verify_windows_package.py `
    --zip dist/fast-gfn1-xtb-4.7.0-windows-x64-portable.zip `
    --output dist/verification-windows-portable
```

Use the matching toolchain/runtime directories for another installation.
Output ZIPs, staging directories and verification directories must be new;
earlier evidence is never overwritten. The binary is rebuilt before packaging.
Review `RESULTS.json` in the verification directory before distributing the ZIP.
This workflow creates local artifacts; it does not upload or publish them.

The verifier extracts the ZIP, verifies its payload/source hashes, runs the
PowerShell installer into a fresh path containing Japanese characters and
spaces, and invokes the installed launcher. It removes all development paths
and oneAPI environment variables, verifies that a bare executable fails with
STATUS_DLL_NOT_FOUND, then checks Energy/Gradient, gas/ALPB analytic Hessian,
actual two-worker OpenMP, the SSE4.2 dispatch path and loaded DLL origins.
Only Windows system libraries may be loaded outside the installed package.
The verification installer omits `-AddToPath` and checks that the real user
PATH is unchanged. The public `install.cmd` uses `-AddToPath` as documented.

The runtime redistribution notices and the consumer instructions are maintained
under `packaging/windows/`. End users extract the entire ZIP and run
`install.cmd` (per-user PATH installation), or use `fast-gfn1-xtb.cmd` directly.
No compiler, Python, administrator rights or network access is needed on the
receiving Windows PC.
