/** Loads one Rust+C++ module. Optional Emscripten options include wasmBinary and locateFile. */
export default function load(options?: {
    wasmBinary?: Uint8Array | ArrayBuffer;
    locateFile?: (path: string, prefix: string) => string;
}): Promise<unknown>;
export class Session {
    constructor();
    free(): void;
    bodyCreateSolidBlock(x: number, y: number, z: number): number;
    bodyCreateSolidBlockAt(x: number, y: number, z: number, location: ArrayLike<number>, axis: ArrayLike<number>, refDirection: ArrayLike<number>): number;
    bodyAskFaces(body: number): Int32Array;
    bodyAskEdges(body: number): Int32Array;
    bodyAskVertices(body: number): Int32Array;
    bodyCheck(body: number): number;
    bodyAskVolume(body: number): number;
    bodyAskBoundingBox(body: number): Float64Array;
    vertexAskPoint(vertex: number): Float64Array;
}
export class SmlibBackend {
    constructor();
    free(): void;
    circle(center: ArrayLike<number>, radius: number): SmlibCurve;
    lineSegment(start: ArrayLike<number>, end: ArrayLike<number>): SmlibCurve;
    cornerSurface(corners: ArrayLike<number>): SmlibSurface;
}
export class SmlibCurve {
    private constructor();
    free(): void;
    domain(): Float64Array;
    evaluate(t: number): Float64Array;
    reversed(): SmlibCurve;
}
export class SmlibSurface {
    private constructor();
    free(): void;
    domain(): Float64Array;
    evaluate(u: number, v: number): Float64Array;
}
