"""Verify an extracted source tree or its bundled dependency before building."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath

MCTC_COMMIT = '77f65c6f2cf6330d05d0757ca173da097096780e'


def verify(root, dependencies_only=False):
    root = root.resolve()
    manifest = json.loads((root/'SOURCE_SHA256.json').read_text(encoding='utf-8'))
    if manifest['dependency_mctc_lib_commit'] != MCTC_COMMIT:
        raise AssertionError('Unexpected dependency commit in package manifest')
    files = manifest['files']
    if dependencies_only:
        files = {name: digest for name, digest in files.items()
                 if name.startswith('subprojects/mctc-lib/')}
        if not files:
            raise AssertionError('Bundled dependency hashes are absent')
        actual = {p.relative_to(root).as_posix() for p in (root/'subprojects/mctc-lib').rglob('*') if p.is_file()}
        if actual != set(files):
            raise AssertionError('Bundled dependency file inventory differs from package manifest')
    for name, expected in files.items():
        rel = PurePosixPath(name)
        if rel.is_absolute() or '..' in rel.parts or '\\' in name or ':' in name:
            raise AssertionError('Invalid manifest path: '+name)
        path = (root/rel).resolve()
        if not path.is_relative_to(root):
            raise AssertionError('Manifest path escapes package root')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise AssertionError('Source checksum differs: '+name)
    return len(files)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--dependencies-only', action='store_true')
    args = parser.parse_args()
    count = verify(args.root, args.dependencies_only)
    print('PASS extracted '+('dependency' if args.dependencies_only else 'source')+' hashes:', count)


if __name__ == '__main__':
    main()
