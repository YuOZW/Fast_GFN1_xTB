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
