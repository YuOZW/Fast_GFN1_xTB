"""Complete CLI solvent Hessians: supported dispatch, parity and fallback."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('hess', ROOT/'tests/gfn1_fast_4_2_0/test_runtime.py')
hess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hess)
SMOOTH = '3\nsmooth nodal water\nO 0 0 0\nH .758602 0 .504284\nH -.3 .67 .62\n'
BOUNDARY = SMOOTH.replace('-.3 .67 .62', '-.3 .65 .62')


def smooth_halogen():
    carbon = [-1.94, .04, -.05]
    oxygen = [3.2, .17, -.13]
    coordinates = [[0, 0, 0], carbon]
    coordinates += [[a+b for a, b in zip(carbon, delta)] for delta in
                    [[-.36, 1.026, .03], [-.37, -.51, .89], [-.35, -.52, -.88]]]
    coordinates += [oxygen]+[[a+b for a, b in zip(oxygen, delta)] for delta in [[.57, .74, .12], [.58, -.75, -.09]]]
    # Frozen trial 5 selected by the independent original surface margins.
    coordinates = [[value+.1*math.sin(95+7*(i+1)+3*(a+1)) for a, value in enumerate(row)]
                   for i, row in enumerate(coordinates)]
    symbols = ['Br', 'C', 'H', 'H', 'H', 'O', 'H', 'H']
    return '8\nsmooth CH3Br-water\n'+'\n'.join(symbol+' '+' '.join(f'{v:.16f}' for v in row)
                                                for symbol, row in zip(symbols, coordinates))+'\n'


def observed(result):
    return {name: float(re.search(pattern+r'\s+('+hess.NUMBER+')', result[2])[1].replace('D', 'E'))
            for name, pattern in [('energy_Eh', 'TOTAL ENERGY'), ('gradient_norm_Eh_bohr', 'GRADIENT NORM')]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    exe, output = args.exe.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = {'exe': str(exe), 'exe_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(), 'cases': [], 'benchmarks': []}
    for model in ('gbsa', 'alpb'):
        options = ('--'+model, 'water')
        for name, extra in [('neutral', ()), ('open', ('--chrg', '1', '--uhf', '1')),
                            ('highT', ('--etemp', '30000'))]:
            label = model+'_'+name
            numerical = hess.run(exe, output, label+'_numerical', SMOOTH, (*options, *extra), step=.0002)
            analytic = hess.run(exe, output, label+'_analytic', SMOOTH, (*options, *extra), hess.ENABLE, step=.0002)
            if analytic[0]['analytic_used'] is not True or analytic[0]['reason'] != 'analytic solvated GFN1':
                raise AssertionError('Complete solvent analytic path was not used: '+str(analytic[0]))
            error = hess.compare(numerical, analytic, label, 2e-6)
            report['cases'].append({'case': label, 'max_abs_H_difference_Eh_bohr2': error})
    for molecule, xyz in [('disilane', hess.grad.DISILANE), ('smooth_halogen', smooth_halogen())]:
        for model in ('gbsa', 'alpb'):
            label = molecule+'_'+model
            options = ('--'+model, 'water')
            numerical = hess.run(exe, output, label+'_numerical', xyz, options, step=.0001)
            analytic = hess.run(exe, output, label+'_analytic', xyz, options, hess.ENABLE, step=.0001)
            if analytic[0]['analytic_used'] is not True:
                raise AssertionError('Broader solvent response rejected: '+str(analytic[0]))
            error = hess.compare(numerical, analytic, label, 2e-6)
            report['cases'].append({'case': label, 'max_abs_H_difference_Eh_bohr2': error})
    for label, xyz, options, policy, control in [
        ('screened_surface', BOUNDARY, ('--gbsa', 'water'), hess.ENABLE, ''),
        ('conservative_disilane_window', hess.grad.DISILANE, ('--alpb', 'water'), hess.ENABLE, ''),
        ('solvent_memory', hess.grad.DISILANE, ('--alpb', 'water'),
         {**hess.ENABLE, 'XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB': '1'}, ''),
        ('solvent_legacy', SMOOTH, ('--alpb', 'water'),
         {**hess.ENABLE, 'XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT': '1'}, ''),
        ('solvent_salt', SMOOTH, ('--gbsa', 'water'), hess.ENABLE, '$gbsa\nion_st=0.1\n$end\n'),
    ]:
        numerical_policy = {**policy, **hess.DISABLE}
        numerical = hess.run(exe, output, label+'_numerical', xyz, options, numerical_policy, step=.0002, control=control)
        candidate = hess.run(exe, output, label+'_fallback', xyz, options, policy, step=.0002, control=control)
        if candidate[0]['analytic_used'] is not False:
            raise AssertionError('Unsupported solvent request did not fall back: '+str(candidate[0]))
        error = hess.compare(numerical, candidate, label, 1e-10)
        report['cases'].append({'case': label, 'max_abs_H_difference_Eh_bohr2': error,
                                'reason': candidate[0]['reason']})
    if args.reference:
        stock = args.reference.resolve()
        report['stock_sha256'] = hashlib.sha256(stock.read_bytes()).hexdigest()
        report['stock_hessian_compatibility'] = 'Stock numerical arithmetic; Hessian OpenMP directives suppressed for Windows ifx descriptor crash; one thread only'
        rng = random.Random(20261005)
        for molecule, xyz in [('disilane', hess.grad.DISILANE), ('smooth_halogen', smooth_halogen())]:
            for model in ('gbsa', 'alpb'):
                options = ('--'+model, 'water')
                label = molecule+'_'+model
                samples = {mode: [] for mode in ('stock', 'fast_numerical', 'fast_analytic')}
                trials = []
                for repeat in range(-1, 7):
                    order = list(samples);rng.shuffle(order)
                    trial = {}
                    for mode in order:
                        trial[mode] = hess.run(stock if mode == 'stock' else exe, output,
                            f'bench_{label}_{repeat}_{mode}', xyz, options,
                            hess.ENABLE if mode == 'fast_analytic' else hess.DISABLE, step=.0001, profile=False)
                        if repeat >= 0:
                            samples[mode].append(trial[mode][0]['process_wall_seconds'])
                    hess.compare(trial['stock'], trial['fast_numerical'], 'stock numerical '+label, 1e-9)
                    hess.compare(trial['stock'], trial['fast_analytic'], 'stock analytic '+label, 2e-6)
                    observations = {mode: observed(value) for mode, value in trial.items()}
                    for mode in ('fast_numerical', 'fast_analytic'):
                        if any(abs(observations[mode][key]-observations['stock'][key]) > 1e-10
                               for key in observations['stock']):
                            raise AssertionError('Stock Energy/Gradient parity failed: '+label)
                    if repeat >= 0:
                        trials.append({'order': order, 'observations': observations})
                audit = hess.run(exe, output, label+'_profile_audit', xyz, options, hess.ENABLE, step=.0001)
                if audit[0]['analytic_used'] is not True:
                    raise AssertionError('Timed analytic policy fell back: '+label)
                medians = {mode: statistics.median(values) for mode, values in samples.items()}
                report['benchmarks'].append({'case': label, 'threads': 1, 'profile': False, 'warmups': 1, 'step_bohr': .0001,
                    'samples_seconds': samples, 'median_seconds': medians, 'trials': trials,
                    'speedup_over_stock': medians['stock']/medians['fast_analytic'],
                    'speedup_over_fast_numerical': medians['fast_numerical']/medians['fast_analytic']})
                print('BENCH solvent Hessian', label, medians, flush=True)
        if hashlib.sha256(exe.read_bytes()).hexdigest() != report['exe_sha256']:
            raise AssertionError('Executable changed during benchmark')
        if hashlib.sha256(stock.read_bytes()).hexdigest() != report['stock_sha256']:
            raise AssertionError('Stock executable changed during benchmark')
    filename = 'benchmark_summary.json' if args.reference else 'summary.json'
    (output/filename).write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS actual solvent CLI analytic Hessian/fallback regression', flush=True)


if __name__ == '__main__':
    main()
