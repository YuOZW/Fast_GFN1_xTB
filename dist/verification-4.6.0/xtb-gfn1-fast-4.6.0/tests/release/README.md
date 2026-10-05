# Source package preparation

These helpers are prepared for final source/archive verification. They have
not yet been run as a completed release workflow. Do not count their presence
as proof of ZIP, SHA or re-extraction validation.

`create_source_package.py` writes an internal per-file SHA-256 manifest,
includes the pinned mctc-lib source without Git metadata, and writes an
external ZIP SHA-256 file. It checks the requested version against the source
profile marker, checks the pinned dependency commit/clean tracked source, and
checks the source/binary hashes of the supplied final evidence manifest.
Existing outputs are never overwritten. Project source is not modified.

`verify_source_package.py` verifies an extracted package's recorded file
hashes and path containment. Its `--dependencies-only` mode can validate the
bundled dependency when a packaged checkout has no `.git` database.

The Windows build runner now has a manifest-based dependency fallback for
checkouts without dependency Git metadata; prepared Git checkouts retain the
existing pinned-commit check. This new fallback has not yet been executed.
After the benchmark ends, verify that fallback, perform final build/regressions, generate matching evidence,
create the ZIP, verify its external hash, safely extract into a new workspace
directory, verify internal source hashes, and rebuild/re-run regressions
from that extracted tree. Preserve those logs and new binary hashes.

Use only the AGENTS.md Python executable:

```powershell
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
& $PY tests/release/create_source_package.py --version 4.6.0 --evidence tests/gfn1_fast_4_6_0/SOURCE_MANIFEST_20261005.json --output dist/xtb-gfn1-fast-4.6.0-source.zip
& $PY tests/release/verify_source_package.py --root <extracted-source-root>
```

The paths above are prospective until final 4.6 evidence exists. Packaging is
local; these helpers do not publish, deploy, alter Git state, install runtime
dependencies or bundle Intel compiler/runtime binaries.
