"""Investigate --parallel precedence, actual Hessian teams and timings.

Uses the portable executable with no oneAPI environment. Timings are process
wall time, randomized, one warm-up plus five measured repeats; profiling is
enabled only in separate audits. Does not change production code or packages.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'
WATER = '3\nwater\nO 0 0 0\nH .758602 .13 .504284\nH -.758602 -.18 .62\n'
DISILANE = '8\ndisilane\nSi -1.15 0 0\nSi 1.15 .1 .05\nH -1.65 1.15 .1\nH -1.65 -.6 1\nH -1.7 -.5 -1.05\nH 1.7 1.2 .15\nH 1.65 -.5 1.1\nH 1.7 -.4 -1\n'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    exe, output = args.exe.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    windows = Path(os.environ['SystemRoot'])
    env = {'SystemRoot': str(windows), 'WINDIR': str(windows), 'PATH': str(windows/'System32'),
           'TEMP': os.environ['TEMP'], 'TMP': os.environ['TEMP'], 'USERPROFILE': os.environ['USERPROFILE'],
           'XTBPATH': str(exe.parent.parent/'share/xtb'), 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
           'OMP_DYNAMIC': 'FALSE', 'OMP_STACKSIZE': '64M', 'KMP_STACKSIZE': '64M'}
    report = {'exe': str(exe), 'exe_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
              'environment': env, 'audits': [], 'benchmark': {}, 'scope': 'CLI flag, portable Windows Release'}

    def save():
        (output/'RESULTS.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')

    def run(name, xyz, options, override=None, profile=True, allow_failure=False):
        work = output/name
        work.mkdir()
        (work/'input.xyz').write_text(xyz, encoding='ascii')
        (work/'control.inp').write_text('$hess\nsccacc=0.0000001\nstep=0.00001\n$end\n', encoding='ascii')
        runtime = dict(env, XTB_GFN1_FAST_PROFILE=str(int(profile)))
        runtime.update(override or {})
        command = [str(exe), 'input.xyz', '--gfn', '1', '--norestart', '--acc', '0.0000001',
                   '--input', 'control.inp', *options]
        start = time.perf_counter()
        process = subprocess.run(command, cwd=work, env=runtime, capture_output=True, text=True,
                                 encoding='utf-8', errors='replace', timeout=180)
        elapsed = time.perf_counter()-start
        log = process.stdout+'\n'+process.stderr
        (work/'run.log').write_text(log, encoding='utf-8')
        record = {'case': name, 'command': options, 'environment_overrides': override or {},
                  'returncode': process.returncode, 'process_seconds': elapsed}
        for field, pattern in (
            ('startup_team_threads', r'omp threads\s*:\s*(\d+)'),
            ('actual_response_threads', r'analytic response threads =\s*(\d+)'),
            ('requested_response_workers', r'analytic requested workers =\s*(\d+)'),
            ('energy_Eh', r'TOTAL ENERGY\s+('+NUMBER+')'),
            ('gradient_norm', r'GRADIENT NORM\s+('+NUMBER+')'),
        ):
            match = re.search(pattern, log)
            if match:
                record[field] = float(match[1]) if field in ('energy_Eh', 'gradient_norm') else int(match[1])
        match = re.search(r'analytic Hessian: used=([TF])', log)
        if match:
            record['analytic_used'] = match[1] == 'T'
        match = re.search(r'outer OpenMP\s*=\s*([TF]), max threads=\s*(\d+)', log)
        if match:
            record['numerical_outer_parallel'] = match[1] == 'T'
            record['numerical_max_threads'] = int(match[2])
        matrix = None
        if (work/'hessian').exists():
            matrix = [float(v.replace('D', 'E')) for line in (work/'hessian').read_text().splitlines()
                      if not line.startswith('$') for v in re.findall(NUMBER, line)]
            assert len(matrix) == (3*int(xyz.splitlines()[0]))**2 and all(math.isfinite(v) for v in matrix)
        if not allow_failure:
            assert process.returncode == 0 and 'energy_Eh' in record, (record, log[-1200:])
        return record, matrix, log

    audits = [
        ('water_default_p8', WATER, ['--hess', '--parallel', '8'], {}, 8, 1),
        ('cli_overrides_env1', WATER, ['--hess', '--parallel', '4'], {'XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO': '1'}, 4, 4),
        ('cli_overrides_env8', WATER, ['--hess', '--parallel', '1'], {'OMP_NUM_THREADS': '8', 'XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO': '1'}, 1, 1),
        ('short_P2', WATER, ['--hess', '-P', '2'], {'XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO': '1'}, 2, 2),
        ('last_option_wins', WATER, ['--hess', '-P', '2', '--parallel', '4'], {'XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO': '1'}, 4, 4),
        ('thread_limit2', DISILANE, ['--hess', '--parallel', '8'], {'OMP_THREAD_LIMIT': '2', 'XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO': '1'}, 2, 2),
        ('alpb_p4', DISILANE, ['--hess', '--alpb', 'water', '--parallel', '4'], {'XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO': '1'}, 4, 4),
    ]
    for name, xyz, options, overrides, maximum, actual in audits:
        result, _, _ = run(name, xyz, options, overrides)
        assert result['startup_team_threads'] == maximum and result['actual_response_threads'] == actual, result
        report['audits'].append(result)
        save()
        print('PASS', name, 'OMP', maximum, 'actual Hessian workers', actual, flush=True)
    for name, options in [('zero', ['--grad', '--parallel', '0']),
                          ('negative', ['--grad', '--parallel', '-1']),
                          ('invalid', ['--grad', '--parallel', 'abc']),
                          ('missing', ['--grad', '--parallel'])]:
        result, _, log = run(name, WATER, options, {'OMP_NUM_THREADS': '8'}, allow_failure=True)
        result['diagnostics'] = [line for line in log.splitlines() if any(word in line.lower() for word in ('warning', 'error', 'invalid', 'thread', 'integer'))]
        report['audits'].append(result)
        save()
        print('OBSERVE', name, result['returncode'], result.get('startup_team_threads'), flush=True)
    xyz = (ROOT/'assets/inputs/xyz/taxol.xyz').read_text(encoding='utf-8')
    serial = None
    for threads in (1, 4, 8, 16):
        result, matrix, _ = run(f'taxol_hess_audit_p{threads}', xyz, ['--hess', '--parallel', str(threads)])
        assert result['analytic_used'] and result['actual_response_threads'] == min(threads, 8), result
        if serial is None:
            serial = (result, matrix)
        result['maximum_H_difference'] = max(abs(a-b) for a, b in zip(serial[1], matrix))
        assert result['maximum_H_difference'] < 1e-10 and abs(result['energy_Eh']-serial[0]['energy_Eh']) < 1e-10
        report['audits'].append(result)
        save()
        print('PASS', result['case'], 'actual workers', result['actual_response_threads'], flush=True)
    rng = random.Random(20261006)
    for task in ('grad', 'hess'):
        samples = {str(t): [] for t in (1, 4, 8)}
        # Hessian audit above warms this executable; warm each timing configuration too.
        for repeat in range(-1, 5):
            order = [1, 4, 8]
            rng.shuffle(order)
            for threads in order:
                result, matrix, _ = run(f'taxol_{task}_repeat{repeat}_p{threads}', xyz,
                                       ['--'+task, '--parallel', str(threads)], profile=False)
                assert abs(result['energy_Eh']-serial[0]['energy_Eh']) < 1e-9
                if matrix is not None:
                    assert max(abs(a-b) for a, b in zip(serial[1], matrix)) < 1e-10
                if repeat >= 0:
                    samples[str(threads)].append(result['process_seconds'])
            report['benchmark'][task] = {'samples': samples, 'repeats_completed': repeat+1}
            save()
            print('PROGRESS', task, 'repeat', repeat, flush=True)
        medians = {t: statistics.median(times) for t, times in samples.items()}
        report['benchmark'][task].update(median_seconds=medians,
                                         speedup={t: medians['1']/seconds for t, seconds in medians.items()})
        save()
    report['status'] = 'PASS'
    save()
    print(json.dumps(report['benchmark'], indent=2), flush=True)


if __name__ == '__main__':
    main()
