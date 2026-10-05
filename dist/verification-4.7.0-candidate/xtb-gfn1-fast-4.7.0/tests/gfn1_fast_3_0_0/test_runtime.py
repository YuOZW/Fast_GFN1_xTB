"""FOE integration: actual accepted steps, audits, fallback and paired timings."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[2]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


smoke = load('smoke', 'tests/gfn1_fast_2_2_6/test_windows_smoke.py')
gradient = load('gradient', 'tests/gfn1_fast_2_4_0/test_runtime.py')
FORCE = {'XTB_GFN1_FAST_ENABLE_FERMI_OPERATOR': '1',
         'XTB_GFN1_FAST_DISABLE_FERMI_OPERATOR_AUTOTUNE': '1',
         'XTB_GFN1_FAST_FERMI_OPERATOR_MIN_NAO': '1'}


def profile(log):
    block = log.split('GFN1-fast 3.0.0 Fermi operator profile')[1]
    counts = re.search(r'density-only attempts=(\d+), accepted=(\d+), audits=(\d+), failures=(\d+)', block)
    result = dict(zip(('attempts', 'accepted', 'audits', 'failures'), map(int, counts.groups())))
    result['disabled'] = 'disabled=T' in block.splitlines()[1]
    result['enabled'] = 'enabled=T' in block.splitlines()[1]
    result['full_solves'] = int(re.search(r'full solves=(\d+)', log).group(1))
    number = re.search(r'audit electron count error: FOE=\s*(\S+), full=\s*(\S+)', block)
    result['audit_electron_error'] = float(number[1])
    result['reference_electron_error'] = float(number[2])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--small-only', action='store_true')
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 5:
        parser.error('paired benchmarks require at least five repeats')
    exe, output = args.exe.resolve(), args.output.resolve()
    summary = {'cases': [], 'benchmarks': []}
    references = {}
    molecules = [('disilane', gradient.DISILANE)]
    if not args.small_only:
        molecules.append(('taxol', (ROOT / 'assets/inputs/xyz/taxol.xyz').read_text()))
    for molecule, xyz in molecules:
        cases = [(str(t), ('--etemp', str(t))) for t in (0, 300, 1000, 5000)]
        cases += [('alpb', ('--alpb', 'water')), ('open', ('--chrg', '1', '--uhf', '1'))]
        for suffix, options in cases:
            name = molecule + '_' + suffix
            reference = smoke.run_case(exe, output, name + '_full', xyz, options)
            candidate = smoke.run_case(exe, output, name + '_foe', xyz, options, FORCE)
            parity = gradient.compare(reference, candidate, name)
            ref_stats, stats = profile(reference[2]), profile(candidate[2])
            if ref_stats['enabled'] or ref_stats['attempts']:
                raise AssertionError('experimental FOE must be disabled by default')
            if suffix in ('0', '300', '1000', '5000'):
                if molecule == 'taxol' and suffix in ('1000', '5000'):
                    if not (stats['audits'] >= 1 and stats['failures'] >= 1 and stats['disabled']):
                        raise AssertionError(f'{name}: high-temperature audit rejection not exercised')
                elif not (stats['accepted'] >= 1 and stats['audits'] >= 1 and stats['failures'] == 0):
                    raise AssertionError(f'{name}: density-only steps and successful audit required')
            if stats['accepted'] and stats['full_solves'] >= ref_stats['full_solves']:
                raise AssertionError(f'{name}: accepted FOE did not reduce full solves')
            summary['cases'].append({'case': name, **parity, **stats,
                                     'reference_full_solves': ref_stats['full_solves']})
            if suffix == '300':
                references[molecule] = reference

    xyz = gradient.DISILANE
    gates = [
        ('bounded', {'XTB_GFN1_FAST_FERMI_OPERATOR_MAX_STEPS': '1'}),
        ('size', {'XTB_GFN1_FAST_FERMI_OPERATOR_MIN_NAO': '100000'}),
        ('disable', {'XTB_GFN1_FAST_DISABLE_FERMI_OPERATOR': '1'}),
        ('warmup', {'XTB_GFN1_FAST_FERMI_OPERATOR_WARMUP': '100000'}),
        ('autotune', {'XTB_GFN1_FAST_DISABLE_FERMI_OPERATOR_AUTOTUNE': '0'}),
    ]
    for name, overrides in gates:
        candidate = smoke.run_case(exe, output, 'gate_' + name, xyz, overrides={**FORCE, **overrides})
        parity = gradient.compare(references['disilane'], candidate, 'gate_' + name)
        stats = profile(candidate[2])
        if name == 'bounded' and not (stats['failures'] == 1 and stats['accepted'] == 0 and stats['disabled']):
            raise AssertionError('bounded expansion must fall back')
        if name in ('size', 'disable', 'warmup') and stats['attempts']:
            raise AssertionError(f'{name} gate must prevent FOE attempts')
        summary['cases'].append({'case': 'gate_' + name, **parity, **stats})

    if not args.small_only:
        xyz = (ROOT / 'assets/inputs/xyz/taxol.xyz').read_text()
        for enabled in ('0', '1'):
            samples = {'full': [], 'foe': [], 'autotune': []}
            for repeat in range(args.repeats):
                for mode, policy in [('full', {}), ('foe', FORCE),
                                     ('autotune', {**FORCE, 'XTB_GFN1_FAST_DISABLE_FERMI_OPERATOR_AUTOTUNE': '0'})]:
                    result = smoke.run_case(exe, output, f'bench_{enabled}_{repeat}_{mode}', xyz,
                        overrides={**policy, 'XTB_GFN1_FAST_PROFILE': enabled})
                    gradient.compare(references['taxol'], result, f'bench_{enabled}_{repeat}_{mode}')
                    samples[mode].append(result[0]['wall_seconds'])
            medians = {name: statistics.median(values) for name, values in samples.items()}
            summary['benchmarks'].append({'profile': enabled, 'samples_seconds': samples,
                'median_seconds': medians, 'foe_over_full': medians['foe']/medians['full'],
                'autotune_over_full': medians['autotune']/medians['full']})
            print(f'BENCH profile={enabled}: {medians}', flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
