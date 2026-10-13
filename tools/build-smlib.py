#!/usr/bin/env python3
"""Build the pinned, unmodified SMLib core and Gensolid's private C ABI (no network)."""
import argparse
import concurrent.futures
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

REVISION = 'c9c5979f04cea8c38213004d4f04df0dcfb40d05'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--jobs', type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument('--wasm', action='store_true', help='Build a browser companion module with Emscripten')
    parser.add_argument('--sanitize', action='store_true', help='Enable Wasm AddressSanitizer')
    parser.add_argument('--archive-only', action='store_true', help='Skip companion-module linking')
    args = parser.parse_args()
    if args.sanitize and not args.wasm:
        parser.error('--sanitize requires --wasm')
    root = args.source.resolve()
    out = args.output.resolve()
    if args.wasm:
        out = out / 'wasm32'
    if sys.platform not in ('darwin', 'linux') and not (sys.platform == 'win32' and args.wasm):
        parser.error('This native build supports macOS and Linux; Windows is not yet validated.')
    revision = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != REVISION:
        parser.error(f'Expected SMLib revision {REVISION}, got {revision}')
    if subprocess.check_output(['git', '-C', str(root), 'status', '--porcelain', '--', 'source/SMLib', 'source/SM_API'], text=True).strip():
        parser.error('SMLib and SM_API source must be unmodified')
    out.mkdir(parents=True, exist_ok=True)
    build_lock = (out/'build.lock').open('a+b')
    if sys.platform == 'win32':
        import msvcrt
        build_lock.seek(0)
        build_lock.write(b'0')
        build_lock.flush()
        while True:
            build_lock.seek(0)
            try:
                msvcrt.locking(build_lock.fileno(), msvcrt.LK_NBLCK, 1)
                break
            except OSError:
                time.sleep(1)
    else:
        import fcntl
        fcntl.flock(build_lock, fcntl.LOCK_EX)
    native = Path(__file__).resolve().parents[1] / 'crates/gensolid-ffi/native'
    compiler = 'em++' if args.wasm else os.environ.get('CXX', 'clang++' if sys.platform == 'darwin' else 'c++')
    flags = ['-std=c++17', '-O1', '-fPIC', '-DNDEBUG', '-DSM_NO_DEBUG_CODE', '-DSM_NO_GFX_OUTPUT_CODE', '-w',
             '-I' + str(root/'source/SMLib/inc'), '-I' + str(root/'source/SM_API/inc')]
    if args.wasm:
        flags += ['-fexceptions']
        if args.sanitize:
            flags += ['-fsanitize=address', '-g']
    sources = sorted((root/'source/SMLib/src').glob('*.cpp')) + [root/'source/SM_API/src'/name for name in ('SmApiGeneral.cpp', 'SmApiCurves.cpp', 'SmApiSurfaces.cpp')] + [native/'smlib.cpp']
    fingerprint = hashlib.sha256((revision + str(root) + repr(flags) + subprocess.check_output([compiler, '--version'], text=True)).encode()).hexdigest()
    def compile_one(source):
        obj = out / (source.parent.parent.name + '_' + source.stem + '.o')
        key = hashlib.sha256((fingerprint + str(source)).encode() + source.read_bytes())
        if source.parent == native:
            key.update((native/'smlib.h').read_bytes())
        key = key.hexdigest()
        stamp = obj.with_suffix('.fingerprint')
        changed = not obj.exists() or not stamp.exists() or stamp.read_text() != key
        if changed:
            command = [compiler, *flags, '-c', str(source), '-o', str(obj)]
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                obj.unlink(missing_ok=True)
                raise RuntimeError(f'{source.name}:\n{result.stderr}')
            stamp.write_text(key)
        return obj
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        objects = list(pool.map(compile_one, sources))
    library = out/'libgensolid_smlib.a'
    staging = out/f'libgensolid_smlib.{os.getpid()}.a'
    staging.unlink(missing_ok=True)
    subprocess.run(['emar' if args.wasm else 'ar', 'rcs', str(staging), *map(str, objects)], check=True)
    os.replace(staging, library)
    if args.wasm and not args.archive_only:
        subprocess.run([compiler, str(library), *(['-fsanitize=address', '-g'] if args.sanitize else []), '-O2', '-fexceptions', '--no-entry', '--profiling-funcs',
            '-sMODULARIZE=1', '-sEXPORT_ES6=1', '-sENVIRONMENT=web,worker,node',
            '-sALLOW_MEMORY_GROWTH=1', '-sSTACK_SIZE=8388608', '-sSTACK_OVERFLOW_CHECK=2', '-sFILESYSTEM=0', '-sDISABLE_EXCEPTION_CATCHING=0',
            '-sEXPORTED_FUNCTIONS=["_malloc","_free","_gs_sm_create","_gs_sm_destroy","_gs_sm_call","_gs_sm_result_free"]',
            '-sEXPORTED_RUNTIME_METHODS=["HEAPU8","HEAPU32","HEAPF64"]',
            '-o', str(out/'smlib.mjs')], check=True)
    for name in ('LICENSE', 'THIRD_PARTY_NOTICES.md'):
        shutil.copy2(root/name, out/name)
    print(library)

if __name__ == '__main__':
    try:
        main()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
