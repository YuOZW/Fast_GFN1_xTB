"""Package a verified working tree with pinned source-only dependencies.

Does not build, publish, delete or modify project source. Run only after the
final source/binary/regression evidence is ready. The ZIP hash is external;
the internal manifest covers every payload file, excluding the manifest itself.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]
MCTC_COMMIT = '77f65c6f2cf6330d05d0757ca173da097096780e'


def tracked_payload(repository, include_untracked=True):
    options = ['--cached', '--others', '--exclude-standard'] if include_untracked else ['--cached']
    process = subprocess.run(['git', '-c', 'safe.directory='+str(repository), '-C', str(repository),
                              'ls-files', *options, '-z'],
                             check=True, capture_output=True)
    return [repository/entry.decode('utf-8') for entry in process.stdout.split(b'\0') if entry]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, required=True,
                        help='Final source/binary manifest already validated by the regression workflow')
    args = parser.parse_args()
    if not re.fullmatch(r'\d+\.\d+\.\d+', args.version):
        parser.error('Expected a numeric version such as 4.6.0')
    marker = re.search(r'GFN1-fast ([0-9.]+) analytic Hessian:',
                       (ROOT/'src/xtb/calculator.f90').read_text(encoding='utf-8'))
    if marker is None or marker[1] != args.version:
        parser.error('Package version does not match the compiled-source profile marker')
    dependency = ROOT/'subprojects/mctc-lib'
    actual = subprocess.check_output(['git', '-c', 'safe.directory='+str(dependency), '-C',
                                      str(dependency), 'rev-parse', 'HEAD'], text=True).strip()
    if actual != MCTC_COMMIT:
        raise AssertionError('Dependency is not at the pinned commit')
    subprocess.run(['git', '-c', 'safe.directory='+str(dependency), '-C', str(dependency),
                    'diff', '--quiet', 'HEAD', '--'], check=True)
    evidence = json.loads(args.evidence.read_text(encoding='utf-8'))
    if not evidence.get('source_sha256') or not evidence.get('binaries'):
        raise AssertionError('Final source/binary hash evidence is absent')
    for name, expected in evidence['source_sha256'].items():
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT):
            raise AssertionError('Evidence source path escapes project root')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise AssertionError('Source has changed since validated evidence: '+name)
    for binary in evidence['binaries'].values():
        if hashlib.sha256(Path(binary['path']).read_bytes()).hexdigest() != binary['sha256']:
            raise AssertionError('Validated binary has changed: '+binary['path'])
    files = []
    for path in tracked_payload(ROOT):
        rel = path.relative_to(ROOT)
        if rel.parts[0] in ('.git', '.agents', '.codex', '_reference', 'dist'):
            continue
        if rel.name == 'geometry-test.log' or '__pycache__' in rel.parts:
            continue
        if not path.is_file():
            continue
        files.append(path)
    files.extend(p for p in tracked_payload(dependency, include_untracked=False) if p.is_file())
    payload = {p.relative_to(ROOT).as_posix(): p for p in sorted(set(files))}
    if not any(name.startswith('subprojects/mctc-lib/src/') for name in payload):
        raise AssertionError('Pinned dependency source is missing from payload')
    # Re-extraction tests must use the prescribed Python and the real bundled
    # dependency, without carrying Git databases or compiler build directories.
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() or output.with_suffix(output.suffix+'.sha256').exists():
        raise FileExistsError('Refuse to overwrite an existing package or hash')
    prefix = 'xtb-gfn1-fast-'+args.version
    manifest = {'version': args.version, 'dependency_mctc_lib_commit': MCTC_COMMIT,
                'scope': 'Project source and pinned mctc-lib source; Intel runtime/compiler not bundled',
                'files': {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in payload.items()}}
    with zipfile.ZipFile(output, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, path in payload.items():
            data = path.read_bytes()
            if hashlib.sha256(data).hexdigest() != manifest['files'][name]:
                raise AssertionError('Source changed during packaging: '+name)
            archive.writestr(prefix+'/'+name, data)
        archive.writestr(prefix+'/SOURCE_SHA256.json', json.dumps(manifest, indent=2)+'\n')
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(output.suffix+'.sha256').write_text(digest+'  '+output.name+'\n', encoding='ascii')
    print('PASS source package created:', output, len(payload), 'files', digest)


if __name__ == '__main__':
    main()
