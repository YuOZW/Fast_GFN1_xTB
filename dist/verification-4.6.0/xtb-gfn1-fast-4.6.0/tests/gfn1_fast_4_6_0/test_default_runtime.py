"""Prepared gate for real default analytic dispatch and explicit fallback.

Run after production policy adoption. Empty policy={} deliberately exercises
the production default without ENABLE; DISABLE selects the numerical baseline.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('solvent', ROOT/'tests/gfn1_fast_4_5_0/test_runtime.py')
solvent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solvent)
hess = solvent.hess
DISABLE = {'XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN': '1'}
def screened_fixture():
    # The CLI rotates every Hessian geometry by small fixed Euler angles.
    # Undo that rotation in the input so the classified native overlap node
    # remains a small nonzero node at the actual Hessian reference.
    points = [[0., 0., 0.], [.758602, 1e-16, .504284], [-.3, .65, .62]]
    for point in points:
        for i, j, angle in [(0, 1, -.0003), (0, 2, -.0002), (1, 2, -.0001)]:
            c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
            x, y = point[i], point[j]
            point[i], point[j] = x*c+y*s, -x*s+y*c
    return '3\nnonzero screened overlap after CLI rotation\n'+''.join(
        element+' '+' '.join(f'{value:.17g}' for value in point)+'\n'
        for element, point in zip(('O', 'H', 'H'), points))


TINY = screened_fixture()
# Two well-conditioned water fragments, with the oxygen pair at 60 bohr.
# Isolated Ag2 at this separation fails the parent generalized eigenproblem;
# Ag parameters remain independently exercised in the component cutoff probe.
D3_CUTOFF = ('6\nwater fragments at 60-bohr D3 cutoff\n'
             'O 0 0 0\nH .758602 .13 .504284\nH -.758602 -.18 .62\n'
             'O 31.750632654180 0 0\nH 32.509234654180 .13 .504284\n'
             'H 30.992030654180 -.18 .62\n')


def compare_observables(reference, candidate):
    for pattern in ('TOTAL ENERGY', 'GRADIENT NORM'):
        values = [float(re.search(pattern+r'\s+('+hess.NUMBER+')', item[2])[1])
                  for item in (reference, candidate)]
        if abs(values[0]-values[1]) > 1e-10:
            raise AssertionError('Default dispatch changed original '+pattern)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    exe, output = args.exe.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = {'exe': str(exe), 'exe_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
              'baseline_policy': DISABLE, 'default_policy': {}, 'cases': []}
    cases = [('water', hess.WATER, (), .0005, False),
             ('disilane', hess.grad.DISILANE, (), .0005, False),
             ('benzene', hess.benzene(), (), .0005, False),
             ('water_defaults', hess.WATER, (), .005, True),
             ('disilane_defaults', hess.grad.DISILANE, (), .005, True),
             ('water_zero', hess.WATER, ('--etemp', '0'), .0005, False),
             ('water_open', hess.WATER, ('--chrg', '1', '--uhf', '1'), .0005, False),
             ('water_highT', hess.WATER, ('--etemp', '30000'), .0005, False),
             ('disilane_alpb', hess.grad.DISILANE, ('--alpb', 'water'), .0001, False),
             ('smooth_halogen_gas', solvent.smooth_halogen(), (), .0002, False),
             ('smooth_halogen_gbsa', solvent.smooth_halogen(), ('--gbsa', 'water'), .0001, False)]
    for name, xyz, options, step, defaults in cases:
        reference = hess.run(exe, output, name+'_numerical', xyz, options, DISABLE,
                             step=step, defaults=defaults)
        candidate = hess.run(exe, output, name+'_bare_default', xyz, options, {},
                             step=step, defaults=defaults)
        if reference[0]['analytic_used'] is not False or candidate[0]['analytic_used'] is not True:
            raise AssertionError('Real default/numerical dispatch was not exercised: '+name)
        error = hess.compare(reference, candidate, name, 2e-4 if defaults else 2e-6)
        compare_observables(reference, candidate)
        report['cases'].append({'case': name, 'default_analytic_used': True,
                                'maximum_H_difference_Eh_bohr2': error})
    fallback = [('explicit_disable', hess.WATER, (), {**hess.ENABLE, **DISABLE}, ''),
                ('explicit_enable_zero', hess.WATER, (), {'XTB_GFN1_FAST_ENABLE_ANALYTIC_HESSIAN': '0'}, ''),
                ('constraint', hess.WATER, (), {}, '$constrain\ndistance: 1,2,1.0\n$end\n'),
                ('memory', hess.grad.DISILANE, (), {'XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB': '1'}, ''),
                ('legacy', hess.WATER, (), {'XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT': '1'}, ''),
                ('salt', solvent.SMOOTH, ('--gbsa', 'water'), {}, '$gbsa\nion_st=0.1\n$end\n'),
                ('surface_branch', solvent.BOUNDARY, ('--gbsa', 'water'), {}, ''),
                ('d3_cutoff', D3_CUTOFF, (), {}, ''),
                ('screened_overlap', TINY, (), {}, '')]
    for name, xyz, options, policy, control in fallback:
        reference = hess.run(exe, output, name+'_numerical', xyz, options, {**policy, **DISABLE},
                             step=.0002, control=control)
        candidate = hess.run(exe, output, name+'_bare_default', xyz, options, policy,
                             step=.0002, control=control)
        if candidate[0]['analytic_used'] is not False:
            raise AssertionError('Actual default fallback failed: '+name+' '+str(candidate[0]))
        error = hess.compare(reference, candidate, name, 1e-10)
        compare_observables(reference, candidate)
        report['cases'].append({'case': name, 'default_analytic_used': False,
                                'reason': candidate[0]['reason'], 'maximum_H_difference_Eh_bohr2': error})
    parallel = hess.run(exe, output, 'disilane_default_threads4', hess.grad.DISILANE, policy={}, threads=4)
    reference = hess.run(exe, output, 'disilane_default_threads1', hess.grad.DISILANE, policy={})
    if parallel[0]['analytic_used'] is not True or reference[0]['analytic_used'] is not True:
        raise AssertionError('Default analytic thread dispatch failed')
    error = hess.compare(reference, parallel, 'default analytic thread parity', 1e-10)
    report['cases'].append({'case': 'default_analytic_threads', 'default_analytic_used': True,
                            'maximum_H_difference_Eh_bohr2': error})
    if hashlib.sha256(exe.read_bytes()).hexdigest() != report['exe_sha256']:
        raise AssertionError('Executable changed during validation')
    (output/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS actual default analytic/fallback CLI policy', flush=True)


if __name__ == '__main__':
    main()
