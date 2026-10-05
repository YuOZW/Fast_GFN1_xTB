"""Halogen-containing CLI Hessians, branch fallback and stock comparison."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('hess', ROOT/'tests/gfn1_fast_4_2_0/test_runtime.py')
hess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hess)


def complex_xyz(element='Br', cutoff=False):
    bond = {'Cl': 1.78, 'Br': 1.94, 'I': 2.15, 'At': 2.25}[element]
    carbon = [-bond, .04, -.05]
    oxygen = [20*.529177210903, 0, 0] if cutoff else [3.2, .17, -.13]
    coords = [[0, 0, 0], carbon]
    coords += [[a+b for a,b in zip(carbon,delta)] for delta in
               [[-.36,1.026,.03],[-.37,-.51,.89],[-.35,-.52,-.88]]]
    coords += [oxygen]+[[a+b for a,b in zip(oxygen,delta)] for delta in [[.57,.74,.12],[.58,-.75,-.09]]]
    symbols = [element,'C','H','H','H','O','H','H']
    return '8\nCH3X-water\n'+'\n'.join(f'{s} '+ ' '.join(f'{v:.15f}' for v in xyz)
                                        for s,xyz in zip(symbols,coords))+'\n'


def active_terms(result, expected):
    match = re.search(r'analytic halogen terms = (\d+)',result[2])
    if not result[0]['analytic_used'] or not match or int(match[1]) != expected:
        raise AssertionError('Expected analytic halogen contribution was not exercised')


def observations(result):
    log = result[2]
    return {name: float(re.search(pattern+r'\s+('+hess.NUMBER+')',log)[1].replace('D','E'))
            for name,pattern in [('energy_Eh','TOTAL ENERGY'),('gradient_norm_Eh_bohr','GRADIENT NORM')]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--repeats', type=int, default=7)
    args = parser.parse_args()
    exe,output = args.exe.resolve(),args.output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    report = {'exe':str(exe),'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'cases':[],'benchmarks':[]}
    cases = [(element.lower(),complex_xyz(element),(),0 if element=='Cl' else 1)
             for element in ('Cl','Br','I','At')]
    cases += [('br_zero',complex_xyz(),('--etemp','0'),1),
              ('br_highT',complex_xyz(),('--etemp','30000'),1),
              ('i_open',complex_xyz('I'),('--etemp','3000','--chrg','1','--uhf','1'),1)]
    for name,xyz,options,terms in cases:
        numerical = hess.run(exe,output,name+'_numerical',xyz,options,step=.0004)
        analytic = hess.run(exe,output,name+'_analytic',xyz,options,hess.ENABLE,step=.0004)
        active_terms(analytic,terms)
        error = hess.compare(numerical,analytic,name,2e-6)
        report['cases'].append({'case':name,'halogen_terms':terms,'max_abs_H_difference_Eh_bohr2':error,
                                **observations(analytic)})
    # Equal nearest distances and the SCF 20-bohr pair cutoff are nonsmooth.
    # CODATA2018 Bohr radius above matches the pinned mctc-lib XYZ reader.
    tie = '3\nnearest-neighbour tie\nBr 0 0 0\nC -1.587531632709 0 0\nN 1.587531632709 0 0\n'
    for name,xyz in [('tie',tie),('cutoff',complex_xyz(cutoff=True))]:
        candidate = hess.run(exe,output,name+'_fallback',xyz,policy=hess.ENABLE)
        reference = hess.run(exe,output,name+'_numerical',xyz)
        if candidate[0]['analytic_used'] is not False:
            raise AssertionError(f'{name}: nonsmooth halogen branch accepted')
        error = hess.compare(reference,candidate,name+' fallback',1e-10)
        report['cases'].append({'case':name+'_fallback','dispatch_reason':candidate[0]['reason'],
                                'max_abs_H_difference_Eh_bohr2':error})
    if args.reference:
        if args.repeats<5:
            parser.error('At least five repeats for performance comparisons')
        stock = args.reference.resolve()
        report['stock_sha256'] = hashlib.sha256(stock.read_bytes()).hexdigest()
        report['stock_hessian_compatibility'] = ('Stock numerical formula; OpenMP directives suppressed '
            'for Windows ifx descriptor crash; one thread only')
        rng = random.Random(20261005)
        for element in ('Br','I'):
            xyz = complex_xyz(element)
            times = {k:[] for k in ('stock','fast_numerical','fast_analytic')}
            max_h = 0.
            for repeat in range(-1,args.repeats):
                order = list(times); rng.shuffle(order)
                trial = {}
                for name in order:
                    trial[name] = hess.run(stock if name=='stock' else exe,output,
                        f'bench_{element}_{repeat}_{name}',xyz,
                        policy=hess.ENABLE if name=='fast_analytic' else {},step=.0004,profile=False)
                    if repeat>=0:
                        times[name].append(trial[name][0]['process_wall_seconds'])
                hess.compare(trial['stock'],trial['fast_numerical'],element+' stock numerical',1e-9)
                max_h = max(max_h,hess.compare(trial['stock'],trial['fast_analytic'],element+' stock analytic',2e-6))
                obs = {k:observations(v) for k,v in trial.items()}
                for name in ('fast_numerical','fast_analytic'):
                    for key in obs[name]:
                        if abs(obs[name][key]-obs['stock'][key])>1e-10:
                            raise AssertionError('Stock Energy/Gradient norm parity failed')
            active_terms(hess.run(exe,output,f'bench_{element}_profile',xyz,policy=hess.ENABLE),1)
            medians = {k:statistics.median(v) for k,v in times.items()}
            report['benchmarks'].append({'molecule':'CH3'+element+'-water','profile':False,'warmups':1,'threads':1,
                'samples_seconds':times,'median_seconds':medians,'max_abs_H_difference_Eh_bohr2':max_h,
                'speedup_over_stock':medians['stock']/medians['fast_analytic'],
                'speedup_over_fast_numerical':medians['fast_numerical']/medians['fast_analytic']})
            print('BENCH halogen Hessian',element,medians,flush=True)
    filename = 'benchmark_summary.json' if args.reference else 'summary.json'
    (output/filename).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PASS halogen CLI Hessian regression',flush=True)


if __name__=='__main__':
    main()
