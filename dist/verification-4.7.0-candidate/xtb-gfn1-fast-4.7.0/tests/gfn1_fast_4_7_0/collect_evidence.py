"""Collect completed hash-bound validation and performance records."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
manifest = json.loads((HERE/'SOURCE_MANIFEST_20261005.json').read_text(encoding='utf-8'))
for name,digest in manifest['source_sha256'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
for binary in manifest['binaries'].values():
    assert hashlib.sha256(Path(binary['path']).read_bytes()).hexdigest()==binary['sha256']
release = ROOT/'build-gfn1-fast-current-windows-ifx-release'
benchmark = json.loads((release/'parallel-frozen-benchmark-4.7.0/summary.json').read_text(encoding='utf-8'))
assert benchmark['exe_sha256']['current']==manifest['binaries']['release']['sha256']
assert set(benchmark['models'])=={'gas','gbsa','alpb'}
for record in benchmark['models'].values():
    assert len(record['trials'])==5
    assert record['median_seconds']['current_8']<record['median_seconds']['current_4']
    assert record['median_seconds']['current_1']<record['median_seconds']['baseline_4_6_0']
shutil.copy2(release/'parallel-frozen-benchmark-4.7.0/summary.json',HERE/'BENCHMARK_20261005.json')
report = {'status':'PASS','created_at':datetime.now().astimezone().isoformat(),
          'source_manifest':'SOURCE_MANIFEST_20261005.json','builds':{},
          'production_defaults':{'minimum_AO':64,'maximum_workers':8,'workspace_cap_MiB':1024},
          'large_raw_Hessian_FD':{},'evidence_log_sha256':{}}
for configuration in ('release','debug'):
    build = ROOT/f'build-gfn1-fast-current-windows-ifx-{configuration}'
    assert hashlib.sha256((build/'xtb.exe').read_bytes()).hexdigest()==manifest['binaries'][configuration]['sha256']
    actual = json.loads((build/'parallel-hessian-test/summary.json').read_text(encoding='utf-8'))
    defaults = json.loads((build/'default-policy-test/summary.json').read_text(encoding='utf-8'))
    assert len(actual['cases'])==19
    assert len(defaults['cases'])==21
    assert actual['exe_sha256']==manifest['binaries'][configuration]['sha256']
    shutil.copy2(build/'parallel-hessian-test/summary.json',HERE/f'PARALLEL_CLI_{configuration.upper()}_20261005.json')
    shutil.copy2(build/'default-policy-test/summary.json',HERE/f'DEFAULT_CLI_{configuration.upper()}_20261005.json')
    molecular = build/'parallel-final-molecular-4.7.0.log'
    assert 'PASS 38 real solvated SCC cases' in molecular.read_text(encoding='utf-8')
    logs = [molecular,build/'parallel-complete-cli-4.7.0.log',build/'parallel-default-4.7.0.log',
            build/'parallel-screening-4.7.0.log',build/'parallel-cutoffs-4.7.0.log']
    for path in logs:
        target = HERE/(configuration+'-'+path.name)
        shutil.copy2(path,target)
        report['evidence_log_sha256'][target.name]=hashlib.sha256(target.read_bytes()).hexdigest()
    report['builds'][configuration]={'exe_sha256':actual['exe_sha256'],'actual_parallel_CLI_comparisons':19,
                                    'default_CLI_cases':21,'native_solvent_FD_cases':38,
                                    'screening':'PASS','D3_cutoff_windows':'PASS'}
large_log = release/'parallel-large-full-fd-4.7.0.log'
text = large_log.read_text(encoding='utf-8')
assert 'PASS two large solvated SCC Hessians, all 339 raw columns' in text
numbers = re.findall(r'large full SCC model/columns/step/H error/E-G error/FD CPU:\s+(\d+)\s+(\d+)\s+([\d.E+-]+)\s+([\d.E+-]+)',text)
assert len(numbers)==2
for model,columns,step,error in numbers:
    assert int(columns)==339 and float(error)<2e-6
    report['large_raw_Hessian_FD']['gbsa' if model=='1' else 'alpb'] = {
        'requested_threads':8,'columns':int(columns),'step_bohr':float(step),'maximum_error_Eh_bohr2':float(error)}
for path in [large_log,release/'parallel-resource-tests-4.7.0.log',release/'parallel-nested-api-4.7.0.log',
             release/'parallel-optional-response-4.7.0.log',release/'parallel-sparse-sasa-response-4.7.0.log',
             release/'parallel-sparse-sasa-component-4.7.0.log']:
    target=HERE/path.name
    shutil.copy2(path,target)
    report['evidence_log_sha256'][target.name]=hashlib.sha256(target.read_bytes()).hexdigest()
(HERE/'VALIDATION_20261005.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print('PASS final source/binary, actual teams, fallback, native FD and repeated performance evidence')
