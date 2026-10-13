// Private wasm32 C ABI transport. No native pointer is dereferenced in Rust.
let native;
const sessions = new Set();
export function initialize(module) {
    if (sessions.size) throw new Error('SMLib has live sessions');
    for (const name of ['_malloc', '_free', '_gs_sm_create', '_gs_sm_call', '_gs_sm_destroy', '_gs_sm_result_free']) {
        if (typeof module?.[name] !== 'function') throw new Error(`Invalid SMLib module: ${name}`);
    }
    native = module;
}
function ready() {
    if (!native) throw new Error('Initialize SMLib before creating a backend');
    return native;
}
function allocate(m, bytes) {
    const p = m._malloc(bytes);
    if (!p) throw new Error('SMLib allocation failed');
    return p;
}
export function create(linear, angular, size) {
    const m = ready(), out = allocate(m, 4);
    try {
        m.HEAPU32[out / 4] = 0;
        const code = m._gs_sm_create(linear, angular, size, out);
        const session = m.HEAPU32[out / 4];
        if (!code) sessions.add(session);
        return [code, session];
    } finally { m._free(out); }
}
export function destroy(session) {
    if (sessions.delete(session)) native._gs_sm_destroy(session);
}
export function call(session, op, ids, values, sizes) {
    const m = ready();
    if (!sessions.has(session)) throw new Error('Invalid SMLib session');
    const allocations = [];
    let out = 0;
    function copy(data, width) {
        if (!data.length) return 0;
        const p = allocate(m, data.length * width);
        allocations.push(p);
        // malloc can grow memory: fetch the current view after allocation.
        (width === 8 ? m.HEAPF64 : m.HEAPU32).set(data, p / width);
        return p;
    }
    try {
        const i = copy(ids, 4), v = copy(values, 8), s = copy(sizes, 4);
        out = allocate(m, 32); // GsSmResult wasm32 layout; i64 aligned at byte 24.
        m.HEAPU8.fill(0, out, out + 32);
        const code = m._gs_sm_call(session, op, i, ids.length, v, values.length, s, sizes.length, out);
        const view = new DataView(m.HEAPU8.buffer);
        const vp = view.getUint32(out + 4, true), vn = view.getUint32(out + 8, true);
        const ip = view.getUint32(out + 12, true), ni = view.getUint32(out + 16, true);
        return [code, view.getBigInt64(out + 24, true).toString(),
            m.HEAPF64.slice(vp / 8, vp / 8 + vn), m.HEAPU32.slice(ip / 4, ip / 4 + ni)];
    } finally {
        if (out) { m._gs_sm_result_free(m.HEAPU32[out / 4]); m._free(out); }
        for (const p of allocations) m._free(p);
    }
}
