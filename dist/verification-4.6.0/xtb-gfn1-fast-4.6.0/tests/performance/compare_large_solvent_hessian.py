"""Guarded large solvent CLI Hessians, resource audit and stock comparison.

Run against a frozen executable. The stock build preserves numerical Hessian
arithmetic with only its Windows ifx-incompatible OpenMP directives removed.
All comparisons use the same input, SCC accuracy, displacement and one thread.
"""
import argparse
import ctypes
from ctypes import wintypes
import hashlib
import importlib.util
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
spec = importlib.util.spec_from_file_location('hess', ROOT/'tests/gfn1_fast_4_2_0/test_runtime.py')
hess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hess)


class ProcessMemory(ctypes.Structure):
    _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD)] + [
        (name, ctypes.c_size_t) for name in (
            'PeakWorkingSetSize', 'WorkingSetSize', 'QuotaPeakPagedPoolUsage',
            'QuotaPagedPoolUsage', 'QuotaPeakNonPagedPoolUsage',
            'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage', 'PrivateUsage')]


def run(exe, output, name, xyz, model, step, analytic, profile):
    work = output/name
    work.mkdir(parents=True, exist_ok=True)
    (work/'input.xyz').write_text(xyz, encoding='ascii')
    (work/'control.inp').write_text(f'$hess\nsccacc=0.0000001\nstep={step}\n$end\n', encoding='ascii')
    env = {k: v for k, v in os.environ.items() if not k.startswith('XTB_GFN1_FAST_')}
    runtime = [str(Path(env[k])/'bin') for k in ('CMPLR_ROOT', 'MKLROOT') if k in env]
    env['PATH'] = os.pathsep.join([*runtime, env.get('PATH', '')])
    env.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_STACKSIZE='64M',
               KMP_STACKSIZE='64M', XTBPATH=str(ROOT), XTB_GFN1_FAST_PROFILE=str(int(profile)))
    env.update(hess.ENABLE if analytic else {'XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN': '1'})
    command = [str(exe), 'input.xyz', '--gfn', '1', '--hess', '--norestart',
               '--acc', '0.0000001', '--input', 'control.inp', '--'+model, 'water']
    query = ctypes.WinDLL('psapi', use_last_error=True).GetProcessMemoryInfo
    query.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessMemory), wintypes.DWORD]
    query.restype = wintypes.BOOL
    resources = {'peak_working_set_bytes': 0, 'peak_pagefile_usage_bytes': 0,
                 'maximum_sampled_private_bytes': 0, 'samples': 0, 'sample_interval_seconds': .05}
    def sample(process):
        info = ProcessMemory(); info.cb = ctypes.sizeof(info)
        if not query(wintypes.HANDLE(int(process._handle)), ctypes.byref(info), info.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        for key, value in [('peak_working_set_bytes', info.PeakWorkingSetSize),
                           ('peak_pagefile_usage_bytes', info.PeakPagefileUsage),
                           ('maximum_sampled_private_bytes', info.PrivateUsage)]:
            resources[key] = max(resources[key], value)
        resources['samples'] += 1
    started = time.perf_counter()
    with (work/'run.log').open('w', encoding='utf-8') as log_file:
        with subprocess.Popen(command, cwd=work, env=env, stdout=log_file,
                              stderr=subprocess.STDOUT) as process:
            try:
                while True:
                    sample(process)
                    try:
                        status = process.wait(timeout=.05)
                        sample(process)
                        break
                    except subprocess.TimeoutExpired:
                        if time.perf_counter()-started > 1200:
                            raise TimeoutError('Large solvent Hessian exceeded 1200 seconds')
            except BaseException:
                process.kill(); process.wait()
                raise
    elapsed = time.perf_counter()-started
    log = (work/'run.log').read_text(encoding='utf-8')
    if status or not (work/'hessian').is_file():
        raise AssertionError(f'{name}: exit {status}; see {work / "run.log"}')
    n3 = 3*int(xyz.splitlines()[0])
    matrix = [float(v.replace('D', 'E').replace('d', 'e'))
              for line in (work/'hessian').read_text().splitlines() if not line.startswith('$')
              for v in re.findall(hess.NUMBER, line)]
    if len(matrix) != n3*n3 or not all(math.isfinite(v) for v in matrix):
        raise AssertionError('Invalid complete Hessian: '+name)
    match = re.search(r'GFN1-fast 4\.\d+\.\d+ analytic Hessian: used=([TF])', log)
    reason = re.search(r'dispatch reason = (.*)', log)
    budget = re.search(r'estimated workspace bytes = (\d+)', log)
    result = {'name': name, 'model': model, 'process_wall_seconds': elapsed, 'n3': n3,
              'step_bohr': step, 'resources': resources, 'log': str(work/'run.log'),
              'analytic_used': match[1] == 'T' if match else None,
              'reason': reason[1].strip() if reason else None,
              'estimated_workspace_bytes': int(budget[1]) if budget else None}
    for key, pattern in [('energy_Eh', 'TOTAL ENERGY'), ('gradient_norm_Eh_bohr', 'GRADIENT NORM')]:
        result[key] = float(re.search(pattern+r'\s+('+hess.NUMBER+')', log)[1])
    print('PASS large solvent CLI', name, f'{elapsed:.4f}s', resources, flush=True)
    return result, matrix, log


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--models', nargs='+', choices=['gbsa', 'alpb'], default=['gbsa', 'alpb'])
    parser.add_argument('--step', type=float, default=1e-5)
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if os.name != 'nt':
        parser.error('This resource audit requires Windows')
    if args.step <= 0 or not math.isfinite(args.step):
        parser.error('A positive finite displacement is required')
    if not args.audit_only and (args.reference is None or args.repeats < 5):
        parser.error('Comparison requires --reference and at least five repeats')
    exe, output = args.exe.resolve(), args.output.resolve()
    stock = args.reference.resolve() if args.reference else None
    paths = {'analytic': exe, **({'stock': stock} if stock else {})}
    hashes = {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in paths.items()}
    xyz = (ROOT/'assets/inputs/xyz/taxol.xyz').read_text()
    output.mkdir(parents=True, exist_ok=True)
    report = {'exe_sha256': hashes, 'exe': {k: str(v) for k, v in paths.items()},
              'input_sha256': hashlib.sha256(xyz.encode()).hexdigest(), 'molecule': 'taxol',
              'atoms': 113, 'nao': 350, 'n3': 339, 'step_bohr': args.step, 'scc_accuracy': 1e-7,
              'threads': 1, 'profile_for_timings': False, 'warmups': 1, 'models': [],
              'stock_scope': 'Original numerical formula; only Hessian OpenMP directives removed for Windows ifx; one thread'}
    rng = random.Random(20261005)
    for model in args.models:
        audit = run(exe, output, model+'_profile_audit', xyz, model, args.step, True, True)
        if audit[0]['analytic_used'] is not True:
            raise AssertionError('Guarded large solvent analytic path rejected: '+str(audit[0]))
        if audit[0]['estimated_workspace_bytes'] > 1024*1048576:
            raise AssertionError('Default workspace cap exceeded')
        record = {'model': model, 'profile_audit': audit[0], 'trials': []}
        report['models'].append(record)
        (output/'progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        if args.audit_only:
            continue
        times = {k: [] for k in paths}
        for repeat in range(-1, args.repeats):
            order = list(paths); rng.shuffle(order); trial = {}
            for name in order:
                (output/'live_state.json').write_text(json.dumps({'model': model, 'repeat': repeat,
                    'method': name, 'status': 'running', 'exe_sha256': hashes}, indent=2), encoding='utf-8')
                trial[name] = run(paths[name], output, f'{model}_{repeat}_{name}', xyz,
                                  model, args.step, name == 'analytic', False)
                if repeat >= 0:
                    times[name].append(trial[name][0]['process_wall_seconds'])
                record['completed_process'] = trial[name][0]
                (output/'progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
                for key, path in paths.items():
                    if hashlib.sha256(path.read_bytes()).hexdigest() != hashes[key]:
                        raise AssertionError('Executable changed during comparison: '+key)
            error = hess.compare(trial['stock'], trial['analytic'], model+f' trial {repeat}', 2e-6)
            for key in ('energy_Eh', 'gradient_norm_Eh_bohr'):
                if abs(trial['stock'][0][key]-trial['analytic'][0][key]) > 1e-10:
                    raise AssertionError('Original Energy/Gradient norm differs: '+model)
            if repeat >= 0:
                record['trials'].append({'order': order, 'results': {k: v[0] for k, v in trial.items()},
                                        'maximum_H_difference_Eh_bohr2': error})
            (output/'progress.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        record['samples_seconds'] = times
        record['median_seconds'] = {k: statistics.median(v) for k, v in times.items()}
        record['speedup_over_stock'] = record['median_seconds']['stock']/record['median_seconds']['analytic']
        print('BENCH large solvent Hessian', model, record['median_seconds'], record['speedup_over_stock'], flush=True)
    filename = 'audit_summary.json' if args.audit_only else 'benchmark_summary.json'
    (output/filename).write_text(json.dumps(report, indent=2), encoding='utf-8')
    (output/'live_state.json').write_text(json.dumps({'status': 'completed'}), encoding='utf-8')
    print('PASS complete large solvent Hessian '+('resource audit' if args.audit_only else 'comparison'), flush=True)


if __name__ == '__main__':
    main()
