"""Extract, install and run the ZIP with a clean Windows-only environment.

Checks runtime module origins, MKL CPU dispatch, OpenMP, source/PE closure,
checksums and missing-DLL failure. Tests never modify the real user PATH.
"""
import argparse
import ctypes
from ctypes import wintypes
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import zipfile

from windows_pe import imports, is_system

NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def modules(pid):
    kernel, psapi = ctypes.WinDLL('kernel32', use_last_error=True), ctypes.WinDLL('psapi', use_last_error=True)
    kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    psapi.EnumProcessModulesEx.argtypes = [wintypes.HANDLE, ctypes.POINTER(ctypes.c_void_p), wintypes.DWORD, ctypes.POINTER(wintypes.DWORD), wintypes.DWORD]
    psapi.GetModuleFileNameExW.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.LPWSTR, wintypes.DWORD]
    handle = kernel.OpenProcess(0x410, False, pid)
    if not handle:
        return set()
    result = set()
    try:
        handles = (ctypes.c_void_p*1024)()
        needed = wintypes.DWORD()
        if psapi.EnumProcessModulesEx(handle, handles, ctypes.sizeof(handles), ctypes.byref(needed), 3):
            for module in handles[:min(1024, needed.value//ctypes.sizeof(ctypes.c_void_p))]:
                name = ctypes.create_unicode_buffer(32768)
                if psapi.GetModuleFileNameExW(handle, module, name, len(name)):
                    result.add(str(Path(name.value).resolve()))
    finally:
        kernel.CloseHandle(handle)
    return result


def run(command, cwd, env, log, audit=False, timeout=180):
    loaded, started = set(), time.monotonic()
    with log.open('w', encoding='utf-8') as stream:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT)
        try:
            while process.poll() is None:
                if time.monotonic()-started > timeout:
                    raise TimeoutError(f'Process timed out: {log}')
                if audit:
                    loaded.update(modules(process.pid))
                time.sleep(.01)
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()
    return process.returncode, log.read_text(encoding='utf-8', errors='replace'), sorted(loaded)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--zip', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    archive, output = args.zip.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError('Use a fresh validation directory')
    output.mkdir(parents=True)
    expected = archive.with_suffix('.zip.sha256').read_text(encoding='ascii').split()[0]
    assert digest(archive) == expected
    extracted = output/'extracted'
    with zipfile.ZipFile(archive) as source:
        for name in source.namelist():
            assert (extracted/name).resolve().is_relative_to(extracted.resolve())
        source.extractall(extracted)
    package = extracted/archive.stem
    manifest = json.loads((package/'PACKAGE_MANIFEST.json').read_text(encoding='utf-8'))
    for name, expected in manifest['files'].items():
        path = (package/name).resolve()
        assert path.is_relative_to(package) and digest(path) == expected, name
    dlls = {p.name.lower() for p in (package/'bin').glob('*.dll')}
    for binary in (package/'bin').iterdir():
        if binary.suffix.lower() in ('.dll', '.exe'):
            assert all(is_system(name) or name in dlls for name in imports(binary)), binary
    with zipfile.ZipFile(package/'source/source.zip') as source:
        marker = next(n for n in source.namelist() if n.endswith('/SOURCE_SHA256.json'))
        prefix = marker.removesuffix('SOURCE_SHA256.json')
        source_manifest = json.loads(source.read(marker))
        for name, expected in source_manifest['files'].items():
            assert hashlib.sha256(source.read(prefix+name)).hexdigest() == expected, name
        assert any(n.startswith('subprojects/mctc-lib/src/') for n in source_manifest['files'])
    # No compiler/runtime/Conda/WSL paths or environment variables survive.
    windows = Path(os.environ['SystemRoot']).resolve()
    powershell = windows/'System32/WindowsPowerShell/v1.0/powershell.exe'
    profile, temp = output/'empty profile', output/'temp'
    profile.mkdir()
    temp.mkdir()
    env = {'SystemRoot': str(windows), 'WINDIR': str(windows), 'COMSPEC': str(windows/'System32/cmd.exe'),
           'PATH': str(windows/'System32'), 'PATHEXT': '.COM;.EXE;.BAT;.CMD',
           'TEMP': str(temp), 'TMP': str(temp), 'USERPROFILE': str(profile),
           'LOCALAPPDATA': str(profile/'AppData/Local')}
    installed = output/'日本語 installation with spaces'
    before_path = os.environ.get('PATH')
    user_path = subprocess.check_output([str(powershell), '-NoProfile', '-Command',
                                        "[Environment]::GetEnvironmentVariable('Path','User')"], env=env)
    code, log, _ = run([str(powershell), '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File',
                        str(package/'install.ps1'), '-Destination', str(installed)], output, env, output/'install.log')
    assert code == 0 and (installed/'bin/fast-gfn1-xtb.exe').exists(), log
    assert os.environ.get('PATH') == before_path
    assert subprocess.check_output([str(powershell), '-NoProfile', '-Command',
                                    "[Environment]::GetEnvironmentVariable('Path','User')"], env=env) == user_path
    for name, expected in manifest['files'].items():
        assert digest(installed/name) == expected
    # Confirm the clean environment really fails without Intel runtimes.
    ctypes.WinDLL('kernel32').SetErrorMode(0x8003)
    bare = output/'missing runtimes'
    bare.mkdir()
    shutil.copy2(installed/'bin/fast-gfn1-xtb.exe', bare/'fast-gfn1-xtb.exe')
    code, _, _ = run([str(bare/'fast-gfn1-xtb.exe'), '--version'], bare, env, bare/'run.log')
    assert code == 0xc0000135, f'Expected STATUS_DLL_NOT_FOUND, got {code:#x}'
    launcher = installed/'fast-gfn1-xtb.cmd'

    def command(options):
        return f'"{env["COMSPEC"]}" /d /s /c ""{launcher}" '+subprocess.list2cmdline(options)+'"'

    code, log, _ = run(command(['--version']), output, env, output/'version.log')
    assert code == 0 and 'normal termination of fast-gfn1-xtb' in log, log
    cases, module_paths = {}, set()
    fixture = (installed/'examples/water.xyz').read_text(encoding='ascii')
    # Avoid a solvent surface branch in the symmetric demonstration geometry.
    hessian_fixture = '3\nasymmetric water\nO 0 0 0\nH .758602 .13 .504284\nH -.758602 -.18 .62\n'
    for name, threads, options, isa in (
        ('gas_gradient', 1, ['--grad'], None),
        ('gas_hessian_serial', 1, ['--hess'], None),
        ('gas_hessian_threads2', 2, ['--hess'], None),
        ('alpb_hessian_threads2', 2, ['--hess', '--alpb', 'water'], None),
        ('sse42_gradient', 1, ['--grad'], 'SSE4_2'),
    ):
        work = output/name
        work.mkdir()
        (work/'input.xyz').write_text(hessian_fixture if '--hess' in options else fixture, encoding='ascii')
        accuracy = []
        if '--hess' in options:
            (work/'control.inp').write_text('$hess\nsccacc=0.0000001\nstep=0.0005\n$end\n', encoding='ascii')
            accuracy = ['--acc', '0.0000001', '--input', 'control.inp']
        runtime = dict(env, OMP_NUM_THREADS=str(threads), XTB_GFN1_FAST_PROFILE='1',
                       XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO='1')
        if isa:
            runtime['MKL_ENABLE_INSTRUCTIONS'] = isa
        code, log, loaded = run(command(['input.xyz', '--gfn', '1', '--norestart', *accuracy, *options]), work,
                                runtime, work/'run.log')
        assert code == 0 and 'convergence criteria satisfied' in log, log[-2000:]
        energy = float(re.search(r'TOTAL ENERGY\s+('+NUMBER+')', log)[1])
        assert math.isfinite(energy)
        cases[name] = {'energy_Eh': energy, 'threads': threads}
        if '--hess' in options:
            assert 'analytic Hessian: used=T' in log
            assert re.search(r'analytic response threads =\s+'+str(threads)+r'\b', log)
            matrix = [float(v.replace('D', 'E')) for line in (work/'hessian').read_text().splitlines()
                      if not line.startswith('$') for v in re.findall(NUMBER, line)]
            assert len(matrix) == 81 and all(math.isfinite(v) for v in matrix)
            cases[name]['hessian'] = matrix
        print('PASS isolated package', name, flush=True)
    assert abs(cases['gas_gradient']['energy_Eh']-cases['sse42_gradient']['energy_Eh']) < 1e-9
    assert max(abs(a-b) for a, b in zip(cases['gas_hessian_serial']['hessian'],
                                      cases['gas_hessian_threads2']['hessian'])) < 1e-10
    # A longer SCC run gives time to inspect actual DLL origins, including dispatch kernels.
    work = output/'module-audit'
    work.mkdir()
    source_manifest_taxol = 'assets/inputs/xyz/taxol.xyz'
    with zipfile.ZipFile(installed/'source/source.zip') as source:
        (work/'input.xyz').write_bytes(source.read(prefix+source_manifest_taxol))
    for isa in (None, 'SSE4_2'):
        runtime = dict(env, XTB_GFN1_FAST_PROFILE='1')
        if isa:
            runtime['MKL_ENABLE_INSTRUCTIONS'] = isa
        code, log, loaded = run([str(installed/'bin/fast-gfn1-xtb.exe'), 'input.xyz', '--gfn', '1',
                                '--grad', '--norestart'], work, runtime, work/(str(isa)+'.log'), audit=True)
        assert code == 0, log[-2000:]
        module_paths.update(loaded)
    assert module_paths, 'No live modules were captured'
    for path in map(Path, module_paths):
        assert path.is_relative_to(installed/'bin') or path.is_relative_to(windows/'System32') or path.is_relative_to(windows/'WinSxS'), path
    for name in dlls:
        matching = [Path(p) for p in module_paths if Path(p).name.lower() == name]
        assert all(p.is_relative_to(installed/'bin') for p in matching), matching
    for required in ('mkl_core.2.dll', 'mkl_sequential.2.dll', 'libifcoremd.dll', 'libiomp5md.dll', 'vcruntime140.dll'):
        assert any(Path(p).name.lower() == required for p in module_paths), required
    assert not any('oneapi' in p.lower() or 'mambaforge' in p.lower() for p in module_paths)
    for case in cases.values():
        case.pop('hessian', None)
    report = {'zip': str(archive), 'zip_sha256': digest(archive), 'environment': env,
              'files_verified': len(manifest['files']), 'source_files_verified': len(source_manifest['files']),
              'missing_runtime_negative_test': 'STATUS_DLL_NOT_FOUND', 'cases': cases,
              'loaded_modules': sorted(module_paths), 'installed': str(installed), 'status': 'PASS'}
    (output/'RESULTS.json').write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print('PASS extracted ZIP / no-oneAPI installation / runtime module origins', flush=True)


if __name__ == '__main__':
    main()
