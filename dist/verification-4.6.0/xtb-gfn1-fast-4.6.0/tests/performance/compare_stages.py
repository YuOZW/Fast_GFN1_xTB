"""Same-host stock/reference and stage-policy comparison, with numerical parity."""
import argparse
import hashlib
import importlib.util
import json
import random
import re
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


smoke = load('smoke', 'tests/gfn1_fast_2_2_6/test_windows_smoke.py')
gradient = load('gradient', 'tests/gfn1_fast_2_4_0/test_runtime.py')
EARLIER_GRADIENT = {
    'XTB_GFN1_FAST_DISABLE_DIRECT_H0_GRADIENT': '1',
    'XTB_GFN1_FAST_DISABLE_PACKED_GRADIENT': '1',
    'XTB_GFN1_FAST_DISABLE_GRADIENT_CONTRACTION': '1',
}
REUSE = {'XTB_GFN1_FAST_ENABLE_SUBSPACE_REUSE': '1',
         'XTB_GFN1_FAST_REUSE_MAX_ITER': '24'}
FOE = {'XTB_GFN1_FAST_ENABLE_FERMI_OPERATOR': '1',
       'XTB_GFN1_FAST_FERMI_OPERATOR_MIN_NAO': '1'}
POLICIES = {
    'stock_6_7_1': {},
    'baseline_2_2_6_policy': EARLIER_GRADIENT,
    'stage_2_3_2_autotune': {**EARLIER_GRADIENT, **REUSE},
    'stage_2_3_2_forced': {**EARLIER_GRADIENT, **REUSE,
                           'XTB_GFN1_FAST_DISABLE_REUSE_AUTOTUNE': '1'},
    'stage_2_4_gradient': {},
    'stage_3_0_0_autotune': FOE,
    'stage_3_0_0_forced': {**FOE, 'XTB_GFN1_FAST_DISABLE_FERMI_OPERATOR_AUTOTUNE': '1'},
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=7)
    args = parser.parse_args()
    if args.repeats < 5:
        parser.error('At least five interleaved repeats are required')
    exe, reference, output = args.exe.resolve(), args.reference.resolve(), args.output.resolve()
    official = ROOT / '_reference/xtb-6.7.1'
    expected = '26b28010e805f7d1aeeef39813feb473e69cc4be'
    commit = subprocess.check_output(['git', '-c', f'safe.directory={official.as_posix()}',
                                     '-C', str(official), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != expected:
        raise AssertionError('Unexpected official reference commit')
    parameter = ROOT / 'param_gfn1-xtb.txt'
    if parameter.read_text() != (official / parameter.name).read_text():
        raise AssertionError('Stock and fast parameter files differ')
    output.mkdir(parents=True, exist_ok=True)
    report = {'reference_commit': commit, 'threads': 1, 'profile_during_timing': False,
              'warmups_per_policy': 1, 'repeats': args.repeats, 'seed': 20261005,
              'timing': 'process wall, including startup and normal CLI output',
              'stage_scope': 'Current executable with explicit stage policies; not archived release binaries',
              'executable_sha256': {name: hashlib.sha256(path.read_bytes()).hexdigest()
                                    for name, path in [('fast', exe), ('stock', reference)]},
              'parameter_sha256': hashlib.sha256(parameter.read_bytes()).hexdigest(),
              'parameter_content_identical_after_newline_normalization': True, 'cases': []}
    # Capture the compatibility-only reference diff so the comparator is reviewable.
    diff = subprocess.check_output(['git', '-c', f'safe.directory={official.as_posix()}',
                                    '-C', str(official), 'diff', '--'], text=True)
    (output / 'reference_compatibility.patch').write_text(diff, encoding='utf-8')
    taxol = (ROOT / 'assets/inputs/xyz/taxol.xyz').read_text()
    cases = [('disilane_gas', gradient.DISILANE, ()),
             ('taxol_gas', taxol, ()), ('taxol_alpb', taxol, ('--alpb', 'water'))]
    rng = random.Random(report['seed'])
    for case, xyz, options in cases:
        runs, times, max_de, max_dg = {}, {name: [] for name in POLICIES}, 0.0, 0.0
        for repeat in range(-1, args.repeats):
            order = list(POLICIES)
            rng.shuffle(order)
            for name in order:
                executable = reference if name == 'stock_6_7_1' else exe
                result = smoke.run_case(executable, output, f'{case}_{repeat}_{name}', xyz, options,
                                        {**POLICIES[name], 'XTB_GFN1_FAST_PROFILE': '0'})
                runs.setdefault(name, []).append(result)
                if repeat >= 0:
                    times[name].append(result[0]['wall_seconds'])
        stock_g = runs['stock_6_7_1'][0][1]
        stock_e = runs['stock_6_7_1'][0][0]['energy_Eh']
        for name, results in runs.items():
            for result, g, log in results:
                de, dg = abs(result['energy_Eh']-stock_e), max(abs(a-b) for a, b in zip(g, stock_g))
                if de > 1e-11 or dg > 1e-10:
                    raise AssertionError(f'{case}/{name} stock parity: dE={de}, dG={dg}')
                max_de, max_dg = max(max_de, de), max(max_dg, dg)
        medians = {name: statistics.median(values) for name, values in times.items()}
        # Profile only after timed trials: confirm experimental paths really ran.
        diagnostics = {}
        for name in POLICIES:
            if name == 'stock_6_7_1':
                continue
            _, _, log = smoke.run_case(exe, output, f'{case}_profile_{name}', xyz, options,
                                       {**POLICIES[name], 'XTB_GFN1_FAST_PROFILE': '1'})
            full = re.search(r'full solves=(\d+)', log)
            diagnostics[name] = {'full_solves': int(full[1]) if full else None}
            for title, label in [('GFN1-fast 2.3.2 subspace reuse profile', 'reuse'),
                                 ('GFN1-fast 3.0.0 Fermi operator profile', 'foe')]:
                if title in log:
                    block = log.split(title, 1)[1]
                    counts = re.search(r'attempts=(\d+), accepted=(\d+), audits=(\d+), failures=(\d+)', block)
                    if counts:
                        diagnostics[name][label] = dict(zip(('attempts','accepted','audits','failures'),
                                                           map(int, counts.groups())))
        item = {'case': case, 'input_sha256': hashlib.sha256(xyz.encode()).hexdigest(),
                'options': options, 'nao': runs['stock_6_7_1'][0][0]['nao'],
                'max_abs_energy_difference_Eh': max_de, 'max_abs_gradient_difference_Eh_bohr': max_dg,
                'samples_seconds': times, 'median_seconds': medians,
                'speedup_over_stock': {name: medians['stock_6_7_1']/t for name,t in medians.items()},
                'relative_to_previous': {
                    'reuse_autotune_over_baseline': medians['stage_2_3_2_autotune']/medians['baseline_2_2_6_policy'],
                    'reuse_forced_over_baseline': medians['stage_2_3_2_forced']/medians['baseline_2_2_6_policy'],
                    'gradient_over_baseline': medians['stage_2_4_gradient']/medians['baseline_2_2_6_policy'],
                    'foe_autotune_over_gradient': medians['stage_3_0_0_autotune']/medians['stage_2_4_gradient'],
                    'foe_forced_over_gradient': medians['stage_3_0_0_forced']/medians['stage_2_4_gradient']},
                'diagnostics': diagnostics}
        report['cases'].append(item)
        (output/'summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(f'BENCH {case}: {medians}',flush=True)
    print(f'PASS same-host stock/stage benchmark: {output / "summary.json"}',flush=True)


if __name__ == '__main__':
    main()
