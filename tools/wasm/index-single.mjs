import createKernel from './gensolid.mjs';
let runtime, pending;
/** Load the single Rust+C++ Wasm module once per JS realm. */
export default function load(options = {}) {
    if (!pending) pending = createKernel(options).then(module => {
        runtime = module;
        return module;
    }).catch(error => { pending = undefined; throw error; });
    return pending;
}
function call(op, id = 0, values = []) {
    if (!runtime) throw new Error('Load Gensolid before creating a session');
    const m = runtime;
    let input = 0, response = 0;
    try {
        if (values.length) {
            input = m._malloc(values.length * 8);
            if (!input) throw new Error('Gensolid input allocation failed');
            m.HEAPF64.set(values, input / 8);
        }
        response = m._gs_web_call(op, id, input, values.length);
        if (!response) throw new Error('Gensolid response allocation failed');
        const code = m._gs_web_code(response);
        if (code) {
            const ptr = m._gs_web_message(response), size = m._gs_web_message_len(response);
            const error = new Error(new TextDecoder().decode(m.HEAPU8.subarray(ptr, ptr + size)));
            error.code = code;
            throw error;
        }
        const ptr = m._gs_web_values(response), size = m._gs_web_values_len(response);
        return m.HEAPF64.slice(ptr / 8, ptr / 8 + size);
    } finally {
        if (response) m._gs_web_free(response);
        if (input) m._free(input);
    }
}
const finalize = new FinalizationRegistry(id => { if (runtime) call(0, id); });
const geometryKey = Symbol();
class Handle {
    constructor(id) { this.#id = id; finalize.register(this, id, this); }
    #id;
    free() {
        if (this.#id) { call(0, this.#id); finalize.unregister(this); this.#id = 0; }
    }
    run(op, values) {
        if (!this.#id) throw new Error('Gensolid handle has been freed');
        return call(op, this.#id, values);
    }
}
function triple(values) {
    if (!values || values.length !== 3) throw new Error('Expected three coordinates');
    return Array.from(values);
}
export class Session extends Handle {
    constructor() { super(call(1)[0]); }
    bodyCreateSolidBlock(x,y,z) { return this.run(10,[x,y,z])[0]; }
    bodyCreateSolidBlockAt(x,y,z,location,axis,refDirection) {
        return this.run(11,[x,y,z,...triple(location),...triple(axis),...triple(refDirection)])[0];
    }
    bodyAskFaces(body) { return Int32Array.from(this.run(12,[body])); }
    bodyAskEdges(body) { return Int32Array.from(this.run(13,[body])); }
    bodyAskVertices(body) { return Int32Array.from(this.run(14,[body])); }
    bodyCheck(body) { return this.run(15,[body])[0]; }
    bodyAskVolume(body) { return this.run(16,[body])[0]; }
    bodyAskBoundingBox(body) { return this.run(17,[body]); }
    vertexAskPoint(vertex) { return this.run(18,[vertex]); }
}
export class SmlibBackend extends Handle {
    constructor() { super(call(2)[0]); }
    circle(center,radius) { return new SmlibCurve(this.run(20,[...triple(center),radius])[0], geometryKey); }
    lineSegment(start,end) { return new SmlibCurve(this.run(21,[...triple(start),...triple(end)])[0], geometryKey); }
    cornerSurface(corners) {
        if (!corners || corners.length !== 12) throw new Error('Expected four XYZ corners');
        return new SmlibSurface(this.run(22,corners)[0], geometryKey);
    }
}
export class SmlibCurve extends Handle {
    constructor(id, key) { if (key !== geometryKey) throw new Error('Use a backend to create curves'); super(id); }
    domain() { return this.run(30); }
    evaluate(t) { return this.run(31,[t]); }
    reversed() { return new SmlibCurve(this.run(32)[0], geometryKey); }
}
export class SmlibSurface extends Handle {
    constructor(id, key) { if (key !== geometryKey) throw new Error('Use a backend to create surfaces'); super(id); }
    domain() { return this.run(40); }
    evaluate(u,v) { return this.run(41,[u,v]); }
}
if (Symbol.dispose) Handle.prototype[Symbol.dispose] = Handle.prototype.free;
/* VALIDATION_EXPORT */
