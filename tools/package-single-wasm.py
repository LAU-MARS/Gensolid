#!/usr/bin/env python3
"""Link Rust and SMLib into one browser Wasm using a unified Emscripten runtime."""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
TOOLCHAIN = 'nightly-2026-10-01'
TARGET = 'wasm32-unknown-emscripten'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=os.environ.get('SMLIB_SOURCE_DIR'))
    parser.add_argument('--output', type=Path, default=ROOT/'target/smlib-single-web')
    parser.add_argument('--validation', action='store_true')
    parser.add_argument('--sanitize', action='store_true', help='Instrument the C++ side with AddressSanitizer')
    args = parser.parse_args()
    if args.sanitize and not args.validation:
        parser.error('--sanitize requires --validation')
    if not args.source:
        parser.error('--source or SMLIB_SOURCE_DIR is required')
    version = subprocess.check_output(['em++', '--version'], text=True).splitlines()[0]
    if not re.search(r'\b6\.0\.12\b', version):
        parser.error('Activate Emscripten 6.0.12; the Rust std and C++ runtime must use matching settings')
    for key in ('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'EMCC_CFLAGS'):
        if os.environ.get(key):
            parser.error(f'Unset {key}; this package fixes the Rust/C++ ABI flags')
    env = os.environ.copy()
    env['SMLIB_SOURCE_DIR'] = str(args.source.resolve())
    env['SMLIB_BUILD_DIR'] = str(ROOT/'target'/('smlib-wasm-asan-cache' if args.sanitize else 'smlib-wasm-cache'))
    env['SMLIB_SANITIZE'] = '1' if args.sanitize else '0'
    # Separate Cargo output prevents accidental reuse of std with different ABI flags.
    target = ROOT/'target'/('single-wasm-rust-asan' if args.sanitize else 'single-wasm-rust')
    env['CARGO_TARGET_DIR'] = str(target)
    env['CARGO_TARGET_WASM32_UNKNOWN_EMSCRIPTEN_RUSTFLAGS'] = '-Cpanic=abort'
    subprocess.run(['cargo', '+'+TOOLCHAIN, 'rustc', '-Zbuild-std=std,panic_abort', '--locked',
                    '-p', 'gensolid-ffi', '--lib', '--crate-type', 'staticlib', '--target', TARGET,
                    '--release', '--features', 'smlib-validation' if args.validation else 'smlib'],
                   cwd=ROOT, env=env, check=True)
    exports = ['_malloc','_free','_gs_web_call','_gs_web_code','_gs_web_values',
               '_gs_web_values_len','_gs_web_message','_gs_web_message_len','_gs_web_free']
    if args.validation:
        exports += ['_gs_sm_create','_gs_sm_destroy','_gs_sm_call','_gs_sm_result_free']
    out = args.output.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='single-wasm-', dir=out.parent) as work:
        staging = Path(work)
        subprocess.run(['em++', str(target/TARGET/'release/libgensolid.a'), *(['-fsanitize=address'] if args.sanitize else []), '-O2', '-fexceptions',
                        '--no-entry', '--profiling-funcs', '-sMODULARIZE=1', '-sEXPORT_ES6=1',
                        '-sENVIRONMENT=web,worker,node', '-sALLOW_MEMORY_GROWTH=1',
                        '-sSTACK_SIZE=8388608', '-sSTACK_OVERFLOW_CHECK=2', '-sFILESYSTEM=0',
                        '-sDISABLE_EXCEPTION_CATCHING=0', '-sEXPORTED_FUNCTIONS='+json.dumps(exports),
                        '-sEXPORTED_RUNTIME_METHODS=["HEAPU8","HEAPU32","HEAPF64"]',
                        '-o', str(staging/'gensolid.mjs')], cwd=ROOT, check=True)
        loader = (ROOT/'tools/wasm/index-single.mjs').read_text()
        loader = loader.replace('/* VALIDATION_EXPORT */',
                'export function validateSmlib() { return call(100)[0]; }' if args.validation else '')
        (staging/'index.mjs').write_text(loader)
        types = (ROOT/'tools/wasm/index-single.d.ts').read_text()
        if args.validation:
            types += '\nexport function validateSmlib(): number;\n'
        (staging/'index.d.ts').write_text(types)
        for filename in ('LICENSE','THIRD_PARTY_NOTICES.md'):
            shutil.copy2(args.source/filename, staging/filename)
        shutil.copy2(ROOT/'NOTICE', staging/'GENSOLID_NOTICE')
        shutil.copy2(ROOT/'LICENSE', staging/'GENSOLID_LICENSE')
        if args.validation:
            shutil.copy2(ROOT/'tools/wasm/validation-single.html', staging/'validation.html')
            shutil.copy2(ROOT/'tools/wasm/validation-worker.mjs', staging/'validation-worker.mjs')
        files = {p.relative_to(staging).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in staging.rglob('*') if p.is_file()}
        (staging/'manifest.json').write_text(json.dumps({'mode':'single','validation':args.validation,'cppSanitizer':args.sanitize,
            'smlibRevision':'c9c5979f04cea8c38213004d4f04df0dcfb40d05','emscripten':version,
            'rustToolchain':TOOLCHAIN,'rustStdRebuilt':True,'panic':'abort','files':files},indent=2)+'\n')
        if out.exists() and any(out.iterdir()):
            old = out/'manifest.json'
            if not old.exists() or json.loads(old.read_text()).get('mode') != 'single':
                parser.error('Use a new output directory; this directory is not a single-Wasm package')
        shutil.copytree(staging,out,dirs_exist_ok=True)
        if not args.validation:
            for filename in ('validation.html','validation-worker.mjs'):
                (out/filename).unlink(missing_ok=True)
        if len(list(out.glob('*.wasm'))) != 1:
            raise RuntimeError('Single-Wasm package must contain exactly one .wasm file')
    print(out/'index.mjs')

if __name__ == '__main__':
    main()
