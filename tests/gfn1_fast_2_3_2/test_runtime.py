"""Independent same-input full/reuse Energy and Gradient regressions."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('smoke', ROOT / 'tests/gfn1_fast_2_2_6/test_windows_smoke.py')
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    exe, output = args.exe.resolve(), args.output.resolve()
    xyz = (ROOT / 'assets/inputs/xyz/taxol.xyz').read_text()
    force = {'XTB_GFN1_FAST_ENABLE_SUBSPACE_REUSE': '1',
             'XTB_GFN1_FAST_DISABLE_REUSE_AUTOTUNE': '1',
             'XTB_GFN1_FAST_REUSE_MAX_ITER': '24'}
    summary = {'cases': [], 'benchmarks': []}
    cases = [('thermal300', ()), ('zeroT', ('--etemp', '0')),
             ('thermal1000', ('--etemp', '1000')),
             ('alpb', ('--alpb', 'water')),
             ('open_shell', ('--chrg', '1', '--uhf', '1'))]
    for name, options in cases:
        baseline, bg, _ = smoke.run_case(exe, output, name + '_full', xyz, options,
                                         {'XTB_GFN1_FAST_DISABLE_SUBSPACE_REUSE': '1'})
        result, rg, log = smoke.run_case(exe, output, name + '_reuse', xyz, options, force)
        de = abs(result['energy_Eh'] - baseline['energy_Eh'])
        dg = max(abs(a - b) for a, b in zip(bg, rg))
        if de > 1e-11 or dg > 1e-10:
            raise AssertionError(f'{name}: dE={de}, dG={dg}')
        stats = re.search(r'\n\s+attempts=(\d+), accepted=(\d+), audits=(\d+), failures=(\d+)', log)
        attempted, accepted, audits, failures = map(int, stats.groups())
        if name in ('thermal300', 'zeroT') and (accepted < 1 or audits < 1):
            raise AssertionError(f'{name}: reuse and same-H full audit were not exercised')
        if failures and name in ('thermal300', 'zeroT'):
            raise AssertionError(f'{name}: unexpected solver/audit failure')
        summary['cases'].append({'case': name, 'abs_dE_Eh': de, 'max_dG_Eh_bohr': dg,
                                 'attempts': attempted, 'accepted': accepted, 'audits': audits,
                                 'solver_or_audit_failures': failures})
        print(f'PASS reuse parity {name}: dE={de:.3e}, max dG={dg:.3e}, accepted={accepted}', flush=True)
        if name == 'thermal300':
            reference, reference_gradient = baseline, bg

    limited_result, limited_gradient, log = smoke.run_case(exe, output, 'iteration_limit', xyz, overrides={
        **force, 'XTB_GFN1_FAST_REUSE_MAX_ITER': '1'})
    if 'failures=1' not in log:
        raise AssertionError('iteration-limit fallback not exercised')
    if abs(limited_result['energy_Eh']-reference['energy_Eh'])>1e-11 or max(
        abs(a-b) for a,b in zip(limited_gradient,reference_gradient))>1e-10:
        raise AssertionError('iteration-limit fallback changed Energy/Gradient')
    limited = smoke.run_case(exe, output, 'size_gate', xyz, overrides={
        **force, 'XTB_GFN1_FAST_REUSE_MIN_NAO': '100000'})[2]
    if 'enabled=F' not in limited:
        raise AssertionError('minimum-size gate not effective')

    # Compare repeated end-to-end runs with profiles both on and off. A
    # performance loss is recorded, never silently labelled a speedup.
    for profile in ('0', '1'):
        times = {'full': [], 'reuse': []}
        for repeat in range(args.repeats):
            for name, policy in [('full', {'XTB_GFN1_FAST_DISABLE_SUBSPACE_REUSE': '1'}), ('reuse', force)]:
                result, _, _ = smoke.run_case(exe, output, f'bench_{profile}_{repeat}_{name}', xyz,
                    overrides={**policy, 'XTB_GFN1_FAST_PROFILE': profile})
                times[name].append(result['wall_seconds'])
        medians = {name: statistics.median(values) for name, values in times.items()}
        summary['benchmarks'].append({'profile': profile, 'samples_seconds': times,
            'median_seconds': medians, 'reuse_over_full': medians['reuse'] / medians['full']})
        print(f'BENCH profile={profile}: {medians}, reuse/full={medians["reuse"] / medians["full"]:.3f}', flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
