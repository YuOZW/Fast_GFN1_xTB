"""Full taxol Hessian parity and interleaved stock/analytic process timings."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('hess',ROOT/'tests/gfn1_fast_4_2_0/test_runtime.py')
hess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hess)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=5)
    args = parser.parse_args()
    if args.repeats<5:
        parser.error('At least five repeats required')
    exe,stock,output = args.exe.resolve(),args.reference.resolve(),args.output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    xyz = (ROOT/'assets/inputs/xyz/taxol.xyz').read_text()
    hashes = {name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in [('fast',exe),('stock',stock)]}
    report = {'molecule':'taxol','atoms':113,'nao':350,'n3':339,'exe_sha256':hashes,
              'input_sha256':hashlib.sha256(xyz.encode()).hexdigest(),'threads':1,'profile':False,
              'warmups':1,'scc_accuracy':1e-7,'step_bohr':.0005,'trials':[],
              'stock_hessian_compatibility':'Original numerical formula; only Hessian OpenMP directives suppressed for Windows ifx; one thread'}
    times = {'stock':[],'analytic':[]}
    rng = random.Random(20261005)
    for repeat in range(-1,args.repeats):
        order = list(times); rng.shuffle(order); trial = {}
        for name in order:
            trial[name] = hess.run(stock if name=='stock' else exe,output,f'taxol_{repeat}_{name}',xyz,
                policy=hess.ENABLE if name=='analytic' else {},profile=False,timeout=900)
            if repeat>=0:
                times[name].append(trial[name][0]['process_wall_seconds'])
        error = hess.compare(trial['stock'],trial['analytic'],f'taxol trial {repeat}',2e-6)
        for pattern in ('TOTAL ENERGY','GRADIENT NORM'):
            vals = [float(re.search(pattern+r'\s+('+hess.NUMBER+')',trial[name][2])[1]) for name in times]
            if abs(vals[0]-vals[1])>1e-10:
                raise AssertionError('Taxol stock Energy/Gradient norm parity failed')
        if hashlib.sha256(exe.read_bytes()).hexdigest()!=hashes['fast']:
            raise AssertionError('Fast executable changed during measurement')
        if repeat>=0:
            report['trials'].append({'order':order,'results':{k:v[0] for k,v in trial.items()},
                                    'max_abs_H_difference_Eh_bohr2':error})
        (output/'progress.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    audit = hess.run(exe,output,'taxol_profile',xyz,policy=hess.ENABLE,timeout=900)
    if audit[0]['analytic_used'] is not True:
        raise AssertionError('Taxol analytic benchmark silently fell back')
    report['profile_audit'] = audit[0]
    report['samples_seconds'] = times
    report['median_seconds'] = {k:statistics.median(v) for k,v in times.items()}
    report['speedup_over_stock'] = report['median_seconds']['stock']/report['median_seconds']['analytic']
    (output/'benchmark_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('BENCH taxol Hessian',report['median_seconds'],report['speedup_over_stock'],flush=True)
    print('PASS large-molecule full Hessian parity and timing',flush=True)


if __name__=='__main__':
    main()
