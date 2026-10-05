"""Save a new source/executable snapshot, refusing to overwrite earlier evidence."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
snapshot = ROOT/'_reference/gfn1-fast-4.7.0-source-snapshot'
snapshot.mkdir()
previous = json.loads((ROOT/'tests/gfn1_fast_4_6_0/SOURCE_MANIFEST_20261005.json').read_text(encoding='utf-8'))
sources = set(previous['source_sha256'])
sources.update(p.relative_to(ROOT).as_posix() for p in (ROOT/'tests/gfn1_fast_4_7_0').rglob('*')
               if p.is_file() and p.suffix.lower() in ('.py','.ps1','.f90','.sh'))
evidence = {'stage':'4.7.0 analytic Hessian OpenMP and exact contractions',
            'created_at':datetime.now().astimezone().isoformat(),
            'scope':'Compiled source and test code with frozen Windows executables; timings recorded separately',
            'source_sha256':{},'binaries':{}}
for name in sorted(sources):
    source = ROOT/name
    target = snapshot/name
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,target)
    evidence['source_sha256'][name] = hashlib.sha256(target.read_bytes()).hexdigest()
for configuration in ('release','debug'):
    source = ROOT/f'build-gfn1-fast-current-windows-ifx-{configuration}/xtb.exe'
    target = snapshot/f'binaries/{configuration}/xtb.exe'
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,target)
    evidence['binaries'][configuration] = {'path':str(target),'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
data = json.dumps(evidence,indent=2)+'\n'
(snapshot/'SOURCE_MANIFEST_20261005.json').write_text(data,encoding='utf-8')
(ROOT/'tests/gfn1_fast_4_7_0/SOURCE_MANIFEST_20261005.json').write_text(data,encoding='utf-8')
print('PASS new immutable candidate snapshot',len(sources),evidence['binaries'])
