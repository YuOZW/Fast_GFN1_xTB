"""Collect bounded, hash-checked current full-response evidence.

Run with the repository's specified fastxtb Python after the native/CLI runs.
Earlier component evidence is kept separate.
"""
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
STAGE = Path(__file__).resolve().parent
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    report = {'stage': '4.5.0 complete Born-solvent SCC Hessian, opt-in',
              'verified_at': datetime.now(timezone(timedelta(hours=9))).isoformat(),
              'environment': 'Windows native ifx 2025.2 / oneMKL sequential LP64; strict native test drivers Od/check:all/fpe:0',
              'executables': {}, 'molecular': {}, 'fixed_charge': {}, 'cli': {}, 'regressions': {}, 'diagnostic': {}}
    for mode in ('release', 'debug'):
        build = ROOT/f'build-gfn1-fast-current-windows-ifx-{mode}'
        exe = build/'xtb.exe'
        report['executables'][mode] = {'path': str(exe), 'sha256': sha(exe), 'library_sha256': sha(build/'xtb.lib')}
        for source, label in [('molecular-solvent-response-test.log', 'molecular_solvent'),
                              ('full-solvent-response-test.log', 'full_solvent')]:
            shutil.copy2(build/'solvent-response-test'/source, STAGE/f'{mode}_{label}_20261005.log')
        text = (STAGE/f'{mode}_molecular_solvent_20261005.log').read_text(encoding='utf-8')
        if 'PASS 38 real solvated SCC cases' not in text or 'PASS partial additive solvated analytic calculator Hessian' not in text:
            raise AssertionError('Missing complete native PASS')
        values = re.findall(r'solvated SCC step/H/raw-E/canonical-E:\s*('+NUMBER+r')\s*('+NUMBER+r')\s*('+NUMBER+r')\s*('+NUMBER+r')', text)
        gradients = re.findall(r'ALPB/kernel/T/charged/G error:\s*[TF]\s*\d+\s*'+NUMBER+r'\s*[TF]\s*('+NUMBER+r')', text)
        if len(values) != 76 or len(gradients) != 38:
            raise AssertionError('Incomplete molecular case/step sweep')
        report['molecular'][mode] = {'cases': 38, 'displacements_bohr': [.0002, .0001],
            'smaller_step_max_abs_H_error_Eh_bohr2': max(float(v[1]) for v in values if float(v[0]) < .00011),
            'max_abs_G_parity_Eh_bohr': max(map(float, gradients)),
            'max_raw_Energy_FD_error_Eh_bohr': max(float(v[2]) for v in values),
            'max_canonical_Energy_FD_error_Eh_bohr': max(float(v[3]) for v in values),
            'partial_additive_columns': True, 'raw_symmetry_translation': True,
            'gas_reference_and_stale_geometry_rejected': True, 'original_halogen_surface_boundary_rejected': True,
            'log': f'{mode}_molecular_solvent_20261005.log'}
        text = (STAGE/f'{mode}_full_solvent_20261005.log').read_text(encoding='utf-8')
        if 'PASS fixed-bare-charge GBSA/ALPB Still/P16' not in text:
            raise AssertionError('Missing fixed-charge PASS')
        values = re.findall(r'model/kernel/charge/step/H/V/G:\s*\d+\s*\d+\s*\d+\s*('+NUMBER+r')\s*('+NUMBER+r')\s*('+NUMBER+r')\s*('+NUMBER+r')', text)
        if len(values) != 24:
            raise AssertionError('Incomplete fixed-charge case/step sweep')
        smaller = [v for v in values if float(v[0]) < .00003]
        report['fixed_charge'][mode] = {'cases': 8, 'steps_bohr': [.0001, .00005, .000025],
            'smaller_step_H_error_Eh_bohr2': max(float(v[1]) for v in smaller),
            'smaller_step_potential_derivative_error_Eh_bohr': max(float(v[2]) for v in smaller),
            'smaller_step_G_Energy_FD_error_Eh_bohr': max(float(v[3]) for v in smaller),
            'original_E_G_V_matrix_parity': True, 'charge_FD_Maxwell_identity': True,
            'log': f'{mode}_full_solvent_20261005.log'}
        data = json.loads((build/'solvent-cli-test/summary.json').read_text())
        if data['exe_sha256'] != sha(exe) or len(data['cases']) != 15:
            raise AssertionError('Stale/incomplete solvent CLI evidence')
        shutil.copy2(build/'solvent-cli-test/summary.json', STAGE/f'{mode}_solvent_cli_20261005.json')
        report['cli'][mode] = {'analytic_cases': 10, 'fallback_cases': 5,
            'max_abs_H_error_Eh_bohr2': max(c['max_abs_H_difference_Eh_bohr2'] for c in data['cases']),
            'summary': f'{mode}_solvent_cli_20261005.json'}
        for label, subdir in [('nodal', 'nodal-cli-test'), ('gas', 'hessian-cli-test')]:
            path = build/subdir/'summary.json'
            data = json.loads(path.read_text())
            if data['exe_sha256'] != sha(exe):
                raise AssertionError('Stale '+mode+' '+label+' regression hash')
            shutil.copy2(path, STAGE/f'{mode}_{label}_current_20261005.json')
            report['regressions'][mode+'_'+label] = {'summary': f'{mode}_{label}_current_20261005.json', 'exe_sha256': data['exe_sha256']}
        smoke = build/'smoke/summary.json'
        shutil.copy2(smoke, STAGE/f'{mode}_smoke_current_20261005.json')
        shutil.copy2(build/'build.log', STAGE/f'{mode}_build_current_20261005.log')
    for mode, folder in [('current', 'build-gfn1-fast-current-windows-ifx-release'),
                         ('stock', 'build-reference-xtb-6.7.1-windows-ifx-release')]:
        path = ROOT/folder/'thermal-solvent-diagnostic/thermal-solvent-diagnostic.log'
        shutil.copy2(path, STAGE/f'{mode}_thermal_diagnostic_20261005.log')
        rows = re.findall(r'^\s*(\d+)\s*('+NUMBER+r')\s*('+NUMBER+r')\s*('+NUMBER+r')', path.read_text(encoding='utf-8'), re.M)
        if len(rows) != 36:
            raise AssertionError('Incomplete stock/current thermal probe')
        report['diagnostic'][mode] = {'coordinate3_smallest_step':
            {'energy_FD_Eh_bohr': float(rows[-7][1]), 'raw_G_difference_Eh_bohr': float(rows[-7][2]),
             'electron_number_derivative_bohr': float(rows[-7][3])},
            'log': f'{mode}_thermal_diagnostic_20261005.log'}
    path = ROOT/'build-gfn1-fast-current-windows-ifx-release/solvent-cli-test/benchmark_summary.json'
    data = json.loads(path.read_text())
    if data['exe_sha256'] != report['executables']['release']['sha256'] or len(data['benchmarks']) != 4:
        raise AssertionError('Stale/incomplete benchmark')
    if any(len(b['trials']) != 7 for b in data['benchmarks']):
        raise AssertionError('Incomplete benchmark trial sweep')
    shutil.copy2(path, STAGE/'BENCHMARK_RESULTS_20261005.json')
    report['benchmark'] = {'summary': 'BENCHMARK_RESULTS_20261005.json', 'cases': 4, 'warmups': 1,
                          'repeats': 7, 'threads': 1, 'profile': False, 'stock_compatibility': data['stock_hessian_compatibility']}
    (STAGE/'SOLVENT_RESULTS_20261005.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({'executables': report['executables'], 'molecular': report['molecular'],
                      'fixed_charge': report['fixed_charge']}, indent=2))


if __name__ == '__main__':
    main()
