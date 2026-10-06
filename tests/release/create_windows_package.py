"""Build a portable Windows ZIP with app-local CPU runtimes and source.

Run with the repository's fastxtb Python. Requires a previously validated
Release build and Intel/VS redistribution files only on the packaging PC.
Never installs software, modifies PATH or overwrites an existing artifact.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

from create_source_package import ROOT, tracked_payload
from windows_pe import imports, is_system


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--build-dir', type=Path, default=ROOT/'build/windows-release')
    parser.add_argument('--compiler-root', type=Path, required=True)
    parser.add_argument('--mkl-root', type=Path, required=True)
    parser.add_argument('--vc-crt-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build = args.build_dir.resolve()
    exe = build/'fast-gfn1-xtb.exe'
    cache = (build/'CMakeCache.txt').read_text(encoding='utf-8')
    for setting in ('CMAKE_BUILD_TYPE:STRING=Release', 'BLA_VENDOR:UNINITIALIZED=Intel10_64lp_seq'):
        # CMake may normalize the cache variable type after configuration.
        key, value = setting.split('=', 1)
        if not re.search(r'^'+key.split(':')[0]+r':[^=]+='+re.escape(value)+r'$', cache, re.M):
            raise ValueError(f'Packaging requires {setting}')
    version = re.search(r'GFN1-fast ([0-9.]+) analytic Hessian:',
                        (ROOT/'src/xtb/calculator.f90').read_text(encoding='utf-8'))[1]
    # Rebuild using the configured Release toolchain to prevent stale binaries.
    subprocess.run(['cmake', '--build', str(build), '--target', 'fast-gfn1-xtb', '--parallel', '4'], check=True)
    compiler, mkl, crt = args.compiler_root.resolve(), args.mkl_root.resolve(), args.vc_crt_dir.resolve()
    fredist = compiler/'share/doc/compiler/fredist.txt'
    allowed_compiler = fredist.read_text(encoding='utf-8').lower()
    roots = {'intel-fortran': compiler/'bin', 'intel-mkl': mkl/'bin', 'microsoft-crt': crt}
    available = {p.name.lower(): (p, component) for component, folder in roots.items()
                 for p in folder.glob('*.dll')}
    if 'mkl_sequential.2.dll' not in imports(exe):
        raise ValueError('Expected the validated sequential MKL build')
    # MKL uses LoadLibrary for CPU dispatch, outside the PE import table.
    dispatch = ['mkl_core.2.dll', 'mkl_def.2.dll', 'mkl_mc3.2.dll', 'mkl_avx2.2.dll', 'mkl_avx512.2.dll']
    dispatch += ['mkl_vml_'+isa+'.2.dll' for isa in ('def', 'cmpt', 'mc3', 'avx2', 'avx512')]
    queue, selected, graph = [exe]+[available[name][0] for name in dispatch], {}, {}
    while queue:
        path = queue.pop()
        name = path.name.lower()
        if name in graph:
            continue
        graph[name] = imports(path)
        if path != exe:
            origin, component = available[name]
            if component == 'intel-fortran' and '/bin/'+name not in allowed_compiler:
                raise ValueError(f'Not on Intel Fortran redistribution list: {name}')
            selected[name] = (origin, component)
        for dependency in graph[name]:
            if not is_system(dependency):
                if dependency not in available:
                    raise FileNotFoundError(f'Unresolved dependency: {name} -> {dependency}')
                queue.append(available[dependency][0])
    output = args.output.resolve()
    staging = output.parent/(output.stem+'-staging')
    source_zip = output.parent/(output.stem+'-source.zip')
    for path in (output, staging, source_zip, output.with_suffix('.zip.sha256')):
        if path.exists():
            raise FileExistsError(f'Refuse to overwrite {path}')
    output.parent.mkdir(parents=True, exist_ok=True)
    staging.mkdir()
    template = ROOT/'packaging/windows'
    shutil.copytree(template, staging, dirs_exist_ok=True)
    # Windows PowerShell 5.1 requires a BOM for UTF-8 scripts with non-ASCII paths/text.
    for script in staging.glob('*.ps1'):
        script.write_text(script.read_text(encoding='utf-8-sig'), encoding='utf-8-sig')
    (staging/'bin').mkdir()
    shutil.copy2(exe, staging/'bin'/exe.name)
    for name, (origin, component) in selected.items():
        shutil.copy2(origin, staging/'bin'/name)
    parameters = staging/'share/xtb'
    parameters.mkdir(parents=True)
    for path in ROOT.glob('param_*.txt'):
        shutil.copy2(path, parameters/path.name)
    shutil.copy2(ROOT/'COPYING', staging/'COPYING')
    for name in ('fortran', 'c', 'openmp'):
        shutil.copytree(compiler/'share/doc/compiler/licensing'/name, staging/'licenses'/('intel-'+name))
    shutil.copy2(fredist, staging/'licenses/intel-fortran/fredist.txt')
    shutil.copytree(mkl/'share/doc/mkl/licensing', staging/'licenses/intel-mkl')
    # Include the exact project source and pinned dependency rather than a URL-only offer.
    evidence = {'source_sha256': {p.relative_to(ROOT).as_posix(): sha(p)
                                 for p in tracked_payload(ROOT) if p.is_file()
                                 and p.relative_to(ROOT).parts[0] not in ('dist', '_reference', '.codex', '.agents')},
                'binaries': {'windows-release': {'path': str(exe), 'sha256': sha(exe)}}}
    evidence_path = build/'windows-package-evidence.json'
    evidence_path.write_text(json.dumps(evidence, indent=2)+'\n', encoding='utf-8')
    subprocess.run([sys.executable, str(ROOT/'tests/release/create_source_package.py'),
                    '--version', version, '--output', str(source_zip), '--evidence', str(evidence_path)], check=True)
    (staging/'source').mkdir()
    shutil.copy2(source_zip, staging/'source/source.zip')
    shutil.copy2(source_zip.with_suffix('.zip.sha256'), staging/'source/source.zip.sha256')
    # Hash-file entry names must match the filename inside the bundle.
    (staging/'source/source.zip.sha256').write_text(sha(source_zip)+'  source.zip\n', encoding='ascii')
    manifest = {'version': version, 'platform': 'windows-x64', 'blas': 'MKL sequential',
                'runtime_origins': {name: component for name, (_, component) in selected.items()},
                'pe_imports': graph, 'files': {p.relative_to(staging).as_posix(): sha(p)
                                             for p in sorted(staging.rglob('*')) if p.is_file()}}
    (staging/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    with zipfile.ZipFile(output, 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(staging.rglob('*')):
            if path.is_file():
                archive.write(path, output.stem+'/'+path.relative_to(staging).as_posix())
    output.with_suffix('.zip.sha256').write_text(sha(output)+'  '+output.name+'\n', encoding='ascii')
    print(f'Created {output} ({output.stat().st_size/1048576:.1f} MiB), {len(selected)} runtime DLLs')


if __name__ == '__main__':
    main()
