"""Interleaved numerical/analytic crossover at 1/4/8 CPU threads.

Stock Hessian comparison remains one-thread only because of the documented
Windows ifx OpenMP descriptor crash. Here compare validated fast algorithms.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import statistics

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('hess', ROOT/'tests/gfn1_fast_4_2_0/test_runtime.py')
hess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hess)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=7)
    args = parser.parse_args()
    if args.repeats < 5:
        parser.error('At least five repeats required for crossover measurements')
    exe, output = args.exe.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    executable_hash = hashlib.sha256(exe.read_bytes()).hexdigest()
    report = {'exe_sha256': executable_hash, 'exe': str(exe), 'profile': False, 'warmups': 1,
              'mkl_threads': 1, 'cases': [], 'stock_scope': 'No stock multithread timing claim; stock one-thread comparisons are preserved separately'}
    rng = random.Random(20261005)
    cases = [('water_gas', hess.WATER, (), .0005, False),
             ('water_gas_defaults', hess.WATER, (), .005, True),
             ('disilane_gas', hess.grad.DISILANE, (), .0005, False),
             ('disilane_gas_defaults', hess.grad.DISILANE, (), .005, True),
             ('benzene_gas', hess.benzene(), (), .0005, False),
             ('disilane_alpb', hess.grad.DISILANE, ('--alpb', 'water'), .0001, False)]
    modes = [(method, threads) for method in ('numerical', 'analytic') for threads in (1, 4, 8)]
    for label, xyz, options, step, defaults in cases:
        times = {f'{method}_{threads}': [] for method, threads in modes}
        trials = []
        maximum_error = 0.0
        for repeat in range(-1, args.repeats):
            order = list(modes);rng.shuffle(order)
            trial = {}
            for method, threads in order:
                name = f'{method}_{threads}'
                policy = hess.ENABLE if method == 'analytic' else {'XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN': '1'}
                trial[name] = hess.run(exe, output, f'{label}_{repeat}_{name}', xyz, options,
                    policy, step=step, threads=threads, profile=False, defaults=defaults)
                if repeat >= 0:
                    times[name].append(trial[name][0]['process_wall_seconds'])
            for method, threads in modes:
                name = f'{method}_{threads}'
                tolerance = 2e-4 if defaults else 2e-6
                error = hess.compare(trial['numerical_1'], trial[name], label+' '+name, tolerance)
                if method == 'numerical' and error > 1e-9:
                    raise AssertionError('Threaded numerical arithmetic changed: '+label+' '+name)
                maximum_error = max(maximum_error, error)
            if hashlib.sha256(exe.read_bytes()).hexdigest() != executable_hash:
                raise AssertionError('Executable changed during measurement')
            if repeat >= 0:
                trials.append({'order': [f'{m}_{t}' for m, t in order],
                               'process_wall_seconds': {k: v[0]['process_wall_seconds'] for k, v in trial.items()}})
            (output/'progress.json').write_text(json.dumps({'report': report, 'active_case': label,
                'completed_trials': len(trials), 'trials': trials}, indent=2), encoding='utf-8')
        for threads in (1, 4, 8):
            audit = hess.run(exe, output, f'{label}_audit_{threads}', xyz, options, hess.ENABLE,
                step=step, threads=threads, defaults=defaults)
            if audit[0]['analytic_used'] is not True:
                raise AssertionError('Timed analytic policy fell back: '+label)
        medians = {k: statistics.median(v) for k, v in times.items()}
        case = {'case': label, 'step_bohr': step, 'default_accuracy': defaults, 'samples_seconds': times,
                'median_seconds': medians, 'trials': trials, 'maximum_H_difference_Eh_bohr2': maximum_error,
                'speedup_analytic_over_numerical_same_threads': {
                    str(t): medians[f'numerical_{t}']/medians[f'analytic_{t}'] for t in (1, 4, 8)}}
        report['cases'].append(case)
        (output/'progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print('BENCH Hessian thread crossover', label, medians, flush=True)
    (output/'benchmark_summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS Hessian crossover measurements', flush=True)


if __name__ == '__main__':
    main()
