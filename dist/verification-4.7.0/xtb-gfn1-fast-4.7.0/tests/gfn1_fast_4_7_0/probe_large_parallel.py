"""Single-run stage/resource diagnostics; not a repeated benchmark."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('large', ROOT/'tests/performance/compare_large_solvent_hessian.py')
large = importlib.util.module_from_spec(spec)
spec.loader.exec_module(large)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', nargs='+', type=int, default=[1, 4, 8])
    args = parser.parse_args()
    os.environ['OMP_DYNAMIC'] = 'FALSE'
    exe, output = args.exe.resolve(), args.output.resolve()
    xyz = (ROOT/'assets/inputs/xyz/taxol.xyz').read_text(encoding='utf-8')
    digest = hashlib.sha256(exe.read_bytes()).hexdigest()
    report = {'exe_sha256': digest, 'scope': 'Single profiled stage/resource diagnostic, not repeated performance', 'cases': []}
    for model in ('gbsa', 'alpb'):
        serial = None
        for threads in args.threads:
            result = large.run(exe, output, f'{model}_threads{threads}', xyz, model, 1e-5, True, True, threads=threads)
            assert result[0]['analytic_used'] is True
            assert re.search(r'analytic response threads =\s+'+str(threads)+r'\b', result[2])
            if serial is None:
                serial = result
            error = large.hess.compare(serial, result, model+' workers '+str(threads), 1e-10)
            record = result[0]
            record['maximum_H_difference_from_serial'] = error
            stage = re.search(r'analytic stage wall .*? = (.*)', result[2])
            record['stage_wall_seconds'] = [float(v) for v in stage[1].split()]
            report['cases'].append(record)
            (output/'progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
            assert hashlib.sha256(exe.read_bytes()).hexdigest() == digest
    (output/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS large worker resource/parity diagnostics', flush=True)


if __name__ == '__main__':
    main()
