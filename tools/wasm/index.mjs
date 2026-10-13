// Ship this whole directory, including snippets/ and both .wasm files.
import initRust, { initializeSmlib } from './gensolid.js';
import createSmlib from './smlib.mjs';
export * from './gensolid.js';
let pending;
/** Load once per JS realm. Optional byte buffers also allow offline/Node loading. */
export default function load({rustWasm, smlibWasm} = {}) {
    if (!pending) pending = (async () => {
        await initRust(rustWasm ? {module_or_path: rustWasm} : undefined);
        const native = await createSmlib(smlibWasm ? {wasmBinary: smlibWasm} : {});
        initializeSmlib(native);
    })().catch(error => { pending = undefined; throw error; });
    return pending;
}
