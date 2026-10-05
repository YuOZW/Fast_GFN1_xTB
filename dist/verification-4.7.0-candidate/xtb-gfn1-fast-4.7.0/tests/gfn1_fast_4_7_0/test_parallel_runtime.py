"""Exercise actual analytic workers, numerical parity and serial switch.

Small molecules force the worker threshold so correctness does not depend on
the production crossover threshold. Timings here are diagnostic, not benchmarks.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('hess', ROOT/'tests/gfn1_fast_4_2_0/test_runtime.py')
hess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hess)
spec = importlib.util.spec_from_file_location('defaults', ROOT/'tests/gfn1_fast_4_6_0/test_default_runtime.py')
defaults = importlib.util.module_from_spec(spec)
spec.loader.exec_module(defaults)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    os.environ['OMP_DYNAMIC'] = 'FALSE'
    exe, output = args.exe.resolve(), args.output.resolve()
    report = {'exe_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(), 'cases': []}
    policy = {**hess.ENABLE, 'XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO': '1'}
    cases = [('water', hess.WATER, ()),
             ('water_zero', hess.WATER, ('--etemp', '0')),
             ('water_open', hess.WATER, ('--chrg', '1', '--uhf', '1')),
             ('water_hot', hess.WATER, ('--etemp', '30000')),
             ('disilane', hess.grad.DISILANE, ()),
             ('benzene', hess.benzene(), ())]
    for name, xyz, options in cases:
        serial = hess.run(exe, output, name+'_serial', xyz, options,
                          {**policy, 'XTB_GFN1_FAST_DISABLE_HESSIAN_OPENMP': '1'}, threads=4)
        assert serial[0]['analytic_used'] is True
        assert re.search(r'analytic response threads =\s+1\b', serial[2])
        for threads in (2, 4):
            parallel = hess.run(exe, output, name+f'_threads{threads}', xyz, options, policy, threads=threads)
            assert parallel[0]['analytic_used'] is True
            assert re.search(r'analytic response threads =\s+'+str(threads)+r'\b', parallel[2])
            error = hess.compare(serial, parallel, name+f' actual workers {threads}', 1e-10)
            report['cases'].append({'case': name, 'threads': threads, 'maximum_H_difference': error})
    reference = hess.run(exe, output, 'disilane_extra_serial', hess.grad.DISILANE, policy=policy)
    for name, overrides, expected in (
            ('eight_workers', {}, 8),
            ('worker_limit', {'XTB_GFN1_FAST_HESSIAN_OMP_MAX_THREADS': '2'}, 2),
            ('runtime_team_limit', {'OMP_THREAD_LIMIT': '2'}, 2)):
        result = hess.run(exe, output, name, hess.grad.DISILANE, policy={**policy, **overrides}, threads=8)
        assert result[0]['analytic_used'] is True
        assert re.search(r'analytic response threads =\s+'+str(expected)+r'\b', result[2])
        error = hess.compare(reference, result, name, 1e-10)
        report['cases'].append({'case': name, 'actual_threads': expected, 'maximum_H_difference': error})
    for name, xyz, options, overrides, control in (
            ('parallel_memory_fallback', hess.grad.DISILANE, (),
             {'XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB': '1'}, ''),
            ('parallel_surface_fallback', defaults.solvent.BOUNDARY, ('--gbsa','water'), {}, ''),
            ('parallel_screening_fallback', defaults.TINY, (), {}, ''),
            ('parallel_constraint_fallback', hess.WATER, (), {}, '$constrain\ndistance: 1,2,1.0\n$end\n')):
        numerical = hess.run(exe, output, name+'_numerical', xyz, options,
                             {**policy, **overrides, **hess.DISABLE}, control=control, threads=4)
        candidate = hess.run(exe, output, name, xyz, options, {**policy, **overrides}, control=control, threads=4)
        assert candidate[0]['analytic_used'] is False
        error = hess.compare(numerical, candidate, name, 1e-10)
        report['cases'].append({'case': name, 'fallback': True, 'maximum_H_difference': error})
    assert hashlib.sha256(exe.read_bytes()).hexdigest() == report['exe_sha256']
    (output/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS actual analytic coordinate/susceptibility workers and serial switch', flush=True)


if __name__ == '__main__':
    main()
