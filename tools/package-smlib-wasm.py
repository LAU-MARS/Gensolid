#!/usr/bin/env python3
"""Build a deployable Rust + SMLib browser package (no native binaries required at runtime)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=os.environ.get('SMLIB_SOURCE_DIR'))
    parser.add_argument('--output', type=Path, default=ROOT/'target/smlib-web')
    parser.add_argument('--validation', action='store_true')
    parser.add_argument('--sanitize', action='store_true')
    args = parser.parse_args()
    if not args.source:
        parser.error('--source or SMLIB_SOURCE_DIR is required')
    if args.sanitize and not args.validation:
        parser.error('--sanitize requires --validation')
    out = args.output.resolve()
    cache = ROOT/'target'/('smlib-wasm-asan-cache' if args.sanitize else 'smlib-wasm-cache')
    command = [sys.executable, str(ROOT/'tools/build-smlib.py'), '--source', str(args.source.resolve()), '--output', str(cache), '--wasm']
    if args.sanitize:
        command += ['--sanitize']
    subprocess.run(command, cwd=ROOT, check=True)
    subprocess.run(['cargo', 'build', '--locked', '-p', 'gensolid-wasm', '--target', 'wasm32-unknown-unknown',
                    '--release', '--features', 'smlib-validation' if args.validation else 'smlib'], cwd=ROOT, check=True)
    target = Path(os.environ.get('CARGO_TARGET_DIR', ROOT/'target'))
    if not target.is_absolute():
        target = ROOT/target
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='smlib-package-', dir=out.parent) as work:
        staging = Path(work)
        subprocess.run(['wasm-bindgen', '--target', 'web', '--out-dir', str(staging), '--out-name', 'gensolid',
                        str(target/'wasm32-unknown-unknown/release/gensolid_wasm.wasm')], check=True)
        for name in ('smlib.mjs', 'smlib.wasm', 'LICENSE', 'THIRD_PARTY_NOTICES.md'):
            shutil.copy2(cache/'wasm32'/name, staging/name)
        shutil.copy2(ROOT/'NOTICE', staging/'GENSOLID_NOTICE')
        shutil.copy2(ROOT/'LICENSE', staging/'GENSOLID_LICENSE')
        shutil.copy2(ROOT/'tools/wasm/index.mjs', staging/'index.mjs')
        if args.validation:
            for name in ('validation.html', 'validation-worker.mjs'):
                shutil.copy2(ROOT/'tools/wasm'/name, staging/name)
        files = {str(p.relative_to(staging)): hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in staging.rglob('*') if p.is_file()}
        (staging/'manifest.json').write_text(json.dumps({
            'smlibRevision': 'c9c5979f04cea8c38213004d4f04df0dcfb40d05',
            'validation': args.validation, 'sanitizer': args.sanitize,
            'emscripten': subprocess.check_output(['em++', '--version'], text=True).splitlines()[0],
            'rust': subprocess.check_output(['rustc', '--version'], text=True).strip(),
            'files': files,
        }, indent=2)+'\n')
        # Keep staging disposable and overwrite only generated package paths.
        shutil.copytree(staging, out, dirs_exist_ok=True)
        if not args.validation:
            for name in ('validation.html', 'validation-worker.mjs'):
                (out/name).unlink(missing_ok=True)
    print(out/'index.mjs')

if __name__ == '__main__':
    main()
