"""Frozen executable comparison: warm-up plus randomized complete Hessians.

Every trial checks the full Cartesian matrix, energy and Gradient norm.
Profiled audits are separate from timing trials. All BLAS calls use MKL=1.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('large', ROOT/'tests/performance/compare_large_solvent_hessian.py')
large = importlib.util.module_from_spec(spec)
spec.loader.exec_module(large)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 5:
        parser.error('At least five measured repetitions are required')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    executables = {'baseline_4_6_0': args.reference.resolve(), 'current': args.exe.resolve()}
    hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name,path in executables.items()}
    os.environ['OMP_DYNAMIC'] = 'FALSE'
    os.environ.pop('OMP_THREAD_LIMIT', None)
    os.environ['OMP_MAX_ACTIVE_LEVELS'] = '1'
    xyz = (ROOT/'assets/inputs/xyz/taxol.xyz').read_text(encoding='utf-8')
    report = {'scope': 'Frozen 4.6.0 versus optimized 4.7.0, actual complete CLI Hessians',
              'exe_sha256': hashes, 'input_sha256': hashlib.sha256(xyz.encode()).hexdigest(),
              'threads': {'baseline_4_6_0':1,'current_1':1,'current_4':4,'current_8':8},
              'warmups_per_configuration': 1, 'measured_repeats':args.repeats,
              'MKL_NUM_THREADS':1,'OMP_DYNAMIC':False,'profile_during_timing':False,'models':{}}
    rng = random.Random(4700)
    for model in ('gas','gbsa','alpb'):
        audits = {}
        for name,threads in report['threads'].items():
            exe = executables['baseline_4_6_0' if name=='baseline_4_6_0' else 'current']
            result = large.run(exe,output,model+'_audit_'+name,xyz,model,1e-5,True,True,threads=threads)
            assert result[0]['analytic_used'] is True
            if name != 'baseline_4_6_0':
                assert re.search(r'analytic response threads =\s+'+str(threads)+r'\b',result[2])
            audits[name] = result
        for name,result in audits.items():
            large.hess.compare(audits['baseline_4_6_0'],result,model+' audit '+name,1e-10)
        record = {'audits':{k:v[0] for k,v in audits.items()},'trials':[]}
        report['models'][model] = record
        samples = {name:[] for name in report['threads']}
        for trial in range(-1,args.repeats):
            order = list(samples);rng.shuffle(order)
            results = {}
            for name in order:
                exe = executables['baseline_4_6_0' if name=='baseline_4_6_0' else 'current']
                result = large.run(exe,output,f'{model}_trial{trial}_{name}',xyz,model,1e-5,True,False,
                                   threads=report['threads'][name])
                error = large.hess.compare(audits['baseline_4_6_0'],result,model+' '+name,1e-10)
                for observable in ('energy_Eh','gradient_norm_Eh_bohr'):
                    assert abs(result[0][observable]-audits['baseline_4_6_0'][0][observable])<=1e-10
                result[0]['maximum_H_difference_from_4_6_0'] = error
                results[name] = result[0]
                if trial>=0:
                    samples[name].append(result[0]['process_wall_seconds'])
                for key,path in executables.items():
                    assert hashlib.sha256(path.read_bytes()).hexdigest()==hashes[key]
            if trial>=0:
                record['trials'].append({'order':order,'results':results})
            (output/'progress.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        record['samples_seconds'] = samples
        medians = {k:statistics.median(v) for k,v in samples.items()}
        record['median_seconds'] = medians
        record['speedup_over_4_6_0'] = {k:medians['baseline_4_6_0']/v for k,v in medians.items()}
        record['parallel_speedup_over_current_1'] = {k:medians['current_1']/v for k,v in medians.items()}
        print('BENCH frozen parallel Hessian',model,medians,flush=True)
    (output/'summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PASS frozen randomized 1/4/8-thread comparison, all matrices and observables',flush=True)


if __name__ == '__main__':
    main()
