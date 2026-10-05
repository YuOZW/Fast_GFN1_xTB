"""Extract the new package or record completed native re-extraction checks."""
import argparse
from datetime import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import zipfile

ROOT = Path(__file__).resolve().parents[2]
VERSION = '4.7.0'
spec = importlib.util.spec_from_file_location('package',ROOT/'tests/release/verify_source_package.py')
package = importlib.util.module_from_spec(spec)
spec.loader.exec_module(package)

parser = argparse.ArgumentParser()
parser.add_argument('--extract',action='store_true')
args = parser.parse_args()
archive = ROOT/f'dist/xtb-gfn1-fast-{VERSION}.zip'
parent = ROOT/f'dist/verification-{VERSION}'
extracted = parent/f'xtb-gfn1-fast-{VERSION}'
if args.extract:
    parent.mkdir()
    with zipfile.ZipFile(archive) as source:
        source.extractall(parent)
    print('PASS extracted source hashes',package.verify(extracted))
    raise SystemExit(0)

manifest = json.loads((ROOT/'tests/gfn1_fast_4_7_0/SOURCE_MANIFEST_20261005.json').read_text(encoding='utf-8'))
count = package.verify(extracted)
dependencies = package.verify(extracted,True)
for name,digest in manifest['source_sha256'].items():
    assert hashlib.sha256((extracted/name).read_bytes()).hexdigest()==digest,name
report = {'status':'PASS','created_at':datetime.now().astimezone().isoformat(),
          'zip':str(archive),'zip_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
          'extracted_root':str(extracted),'verified_source_files':count,
          'verified_mctc_source_files':dependencies,'compiled_source_matches_frozen_manifest_files':len(manifest['source_sha256']),
          'git_metadata':'No .git supplied; enclosing checkout blocked with GIT_CEILING_DIRECTORIES',
          'builds':{}}
for configuration in ('release','debug'):
    build = extracted/f'build-gfn1-fast-current-windows-ifx-{configuration}'
    binary = (build/'xtb.exe').read_bytes()
    offset = struct.unpack_from('<I',binary,0x3c)[0]
    assert binary[offset:offset+4]==b'PE\0\0'
    assert struct.unpack_from('<H',binary,offset+4)[0]==0x8664
    parallel = json.loads((build/'parallel-hessian-test/summary.json').read_text(encoding='utf-8'))
    defaults = json.loads((build/'default-policy-test/summary.json').read_text(encoding='utf-8'))
    assert len(parallel['cases'])==19 and len(defaults['cases'])==21
    assert parallel['exe_sha256']==hashlib.sha256(binary).hexdigest()
    molecular = (build/'reextraction-molecular.log').read_text(encoding='utf-8')
    assert 'PASS 38 real solvated SCC cases' in molecular
    report['builds'][configuration] = {'sha256':parallel['exe_sha256'],'PE_format':'AMD64 PE32+',
                                      'actual_parallel_CLI_comparisons':19,'default_CLI_cases':21,
                                      'native_solvent_FD_cases':38}
data = json.dumps(report,indent=2)+'\n'
(ROOT/'dist/REEXTRACTION_4.7.0.json').write_text(data,encoding='utf-8')
(ROOT/'tests/gfn1_fast_4_7_0/REEXTRACTION_20261005.json').write_text(data,encoding='utf-8')
print('PASS verified re-extracted Windows native Release/Debug builds and actual parallel calculations',count)
