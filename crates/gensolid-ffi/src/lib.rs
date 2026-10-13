//! PK-style C ABI.
//!
//! Conventions follow Parasolid's PK interface so existing integrations port easily:
//!
//! - One global session, started with [`GS_SESSION_start`] and stopped with
//!   [`GS_SESSION_stop`].
//! - Entities are integer tags; `0` is the null tag.
//! - Every function returns a [`GS_ERROR_code_t`]; `GS_ERROR_no_errors` (0) on success.
//! - Options structs carry an `o_t_version` field and are initialised with their
//!   `*_o_init` function.
//! - Arrays returned by the kernel are released with [`GS_MEMORY_free`].
//!
//! The header `include/gensolid.h` is generated from this file with cbindgen.

#![allow(
    unsafe_code,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals
)]

/// Optional native SMLib curve and surface backend.
///
/// cbindgen:ignore
#[cfg(feature = "smlib")]
pub mod smlib;

use std::alloc::{Layout, alloc, dealloc};
use std::ffi::c_void;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::ptr;
use std::sync::Mutex;

use gensolid_api::{Axis2, Error, Precision, Session, SessionOptions, Tag};

/// Error code returned by every API function.
pub type GS_ERROR_code_t = i32;
/// Tag of any entity.
pub type GS_ENTITY_t = i32;
/// Tag of a body.
pub type GS_BODY_t = i32;
/// Tag of a face.
pub type GS_FACE_t = i32;
/// Tag of an edge.
pub type GS_EDGE_t = i32;
/// Tag of a vertex.
pub type GS_VERTEX_t = i32;

/// Success.
pub const GS_ERROR_no_errors: GS_ERROR_code_t = 0;
// Literal values so cbindgen can emit them; `codes_match_kernel_enums` keeps them in sync.
/// An argument is outside its documented domain.
pub const GS_ERROR_invalid_argument: GS_ERROR_code_t = 1;
/// A tag does not refer to an existing entity.
pub const GS_ERROR_not_found: GS_ERROR_code_t = 2;
/// A tag refers to an entity of the wrong class.
pub const GS_ERROR_wrong_entity_class: GS_ERROR_code_t = 3;
/// Geometry is degenerate at session precision.
pub const GS_ERROR_degenerate: GS_ERROR_code_t = 4;
/// A parameter lies outside the valid range.
pub const GS_ERROR_out_of_range: GS_ERROR_code_t = 5;
/// Valid input outside what the kernel supports by design.
pub const GS_ERROR_unsupported: GS_ERROR_code_t = 6;
/// Valid, in-scope input whose algorithm is not written yet.
pub const GS_ERROR_unimplemented: GS_ERROR_code_t = 7;
/// An internal limit was reached.
pub const GS_ERROR_capacity_exceeded: GS_ERROR_code_t = 8;
/// Input could not be parsed.
pub const GS_ERROR_parse: GS_ERROR_code_t = 9;
/// Internal kernel error.
pub const GS_ERROR_internal: GS_ERROR_code_t = 10;
/// No session is running.
pub const GS_ERROR_session_not_started: GS_ERROR_code_t = 100;
/// A session is already running.
pub const GS_ERROR_session_already_started: GS_ERROR_code_t = 101;
/// A required pointer argument was null.
pub const GS_ERROR_null_argument: GS_ERROR_code_t = 102;
/// The options struct has an unknown `o_t_version`.
pub const GS_ERROR_bad_option_version: GS_ERROR_code_t = 103;
/// The kernel panicked; the session should be stopped.
pub const GS_ERROR_panic: GS_ERROR_code_t = 104;
/// Memory allocation failed.
pub const GS_ERROR_out_of_memory: GS_ERROR_code_t = 105;

/// A referenced entity or geometry does not exist.
pub const GS_BODY_state_missing_entity: i32 = 1;
/// A child and its parent do not reference each other.
pub const GS_BODY_state_back_pointer: i32 = 2;
/// The region list is malformed.
pub const GS_BODY_state_region_layout: i32 = 3;
/// A loop has no fins.
pub const GS_BODY_state_empty_loop: i32 = 4;
/// Consecutive fins of a loop do not share a vertex.
pub const GS_BODY_state_open_loop: i32 = 5;
/// A solid edge does not have exactly two fins.
pub const GS_BODY_state_non_manifold_edge: i32 = 6;
/// Both fins of a manifold edge run the same way.
pub const GS_BODY_state_fin_senses: i32 = 7;
/// A face has no surface or an edge has no curve.
pub const GS_BODY_state_missing_geometry: i32 = 8;
/// A vertex is off its edge curve.
pub const GS_BODY_state_vertex_off_curve: i32 = 9;
/// A vertex is off its face surface.
pub const GS_BODY_state_vertex_off_surface: i32 = 10;
/// An edge is off its face surface.
pub const GS_BODY_state_edge_off_surface: i32 = 11;

/// Current version of every options struct.
pub const GS_O_T_VERSION: i32 = 1;

/// Options for [`GS_SESSION_start`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct GS_SESSION_start_o_t {
    /// Must be `GS_O_T_VERSION`.
    pub o_t_version: i32,
    /// Linear precision in model units.
    pub linear_precision: f64,
    /// Angular precision in radians.
    pub angular_precision: f64,
    /// Edge length of the size box in model units.
    pub size_box: f64,
}

/// Coordinate system: location, principal axis, reference direction.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct GS_AXIS2_sf_t {
    /// Origin.
    pub location: [f64; 3],
    /// Principal (z) axis.
    pub axis: [f64; 3],
    /// Reference (x) direction.
    pub ref_direction: [f64; 3],
}

/// One fault found by [`GS_BODY_check`]; `state` is a `GS_BODY_state_*` value.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct GS_check_fault_t {
    /// Fault category.
    pub state: i32,
    /// Entity at fault.
    pub entity: GS_ENTITY_t,
}

static SESSION: Mutex<Option<Session>> = Mutex::new(None);

/// Header in front of every kernel allocation, holding its total size.
const HEADER: usize = 16;

fn code(error: &Error) -> GS_ERROR_code_t {
    error.code() as i32
}

fn guard(f: impl FnOnce() -> Result<(), GS_ERROR_code_t>) -> GS_ERROR_code_t {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(Ok(())) => GS_ERROR_no_errors,
        Ok(Err(code)) => code,
        Err(_) => GS_ERROR_panic,
    }
}

fn with_session<R>(f: impl FnOnce(&mut Session) -> Result<R, Error>) -> Result<R, GS_ERROR_code_t> {
    let mut lock = SESSION.lock().map_err(|_| GS_ERROR_panic)?;
    let session = lock.as_mut().ok_or(GS_ERROR_session_not_started)?;
    f(session).map_err(|e| code(&e))
}

/// Copies `items` into a kernel allocation releasable with [`GS_MEMORY_free`].
fn alloc_array<T: Copy>(items: &[T]) -> Result<*mut T, GS_ERROR_code_t> {
    if items.is_empty() {
        return Ok(ptr::null_mut());
    }
    const { assert!(align_of::<T>() <= HEADER) };
    let bytes = size_of_val(items)
        .checked_add(HEADER)
        .ok_or(GS_ERROR_out_of_memory)?;
    let layout = Layout::from_size_align(bytes, HEADER).map_err(|_| GS_ERROR_out_of_memory)?;
    // SAFETY: `layout` has non-zero size. The size header is written at the start of the
    // block and the items after it; both fit inside the `bytes` allocated.
    unsafe {
        let base = alloc(layout);
        if base.is_null() {
            return Err(GS_ERROR_out_of_memory);
        }
        base.cast::<usize>().write(bytes);
        let data = base.add(HEADER).cast::<T>();
        ptr::copy_nonoverlapping(items.as_ptr(), data, items.len());
        Ok(data)
    }
}

/// Writes `value` through `out`, or fails with `GS_ERROR_null_argument`.
///
/// # Safety
/// `out` must be null or valid for writes of `T`.
unsafe fn write_out<T>(out: *mut T, value: T) -> Result<(), GS_ERROR_code_t> {
    if out.is_null() {
        return Err(GS_ERROR_null_argument);
    }
    // SAFETY: non-null and valid for writes by the caller's contract.
    unsafe { out.write(value) };
    Ok(())
}

/// Writes a tag array result.
///
/// # Safety
/// `n_out` and `array_out` must be null or valid for writes.
unsafe fn write_tags(
    tags: &[Tag],
    n_out: *mut i32,
    array_out: *mut *mut Tag,
) -> Result<(), GS_ERROR_code_t> {
    if n_out.is_null() || array_out.is_null() {
        return Err(GS_ERROR_null_argument);
    }
    let n = i32::try_from(tags.len()).map_err(|_| GS_ERROR_capacity_exceeded)?;
    let array = alloc_array(tags)?;
    // SAFETY: both pointers are non-null and valid for writes by the caller's contract.
    unsafe {
        n_out.write(n);
        array_out.write(array);
    }
    Ok(())
}

/// Initialises session options to the defaults.
///
/// # Safety
/// `options` must be null or point to writable memory for one `GS_SESSION_start_o_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_SESSION_start_o_init(
    options: *mut GS_SESSION_start_o_t,
) -> GS_ERROR_code_t {
    let p = Precision::DEFAULT;
    let defaults = GS_SESSION_start_o_t {
        o_t_version: GS_O_T_VERSION,
        linear_precision: p.linear(),
        angular_precision: p.angular(),
        size_box: p.size(),
    };
    // SAFETY: forwarded caller contract.
    guard(|| unsafe { write_out(options, defaults) })
}

/// Starts the session. `options` may be null for defaults.
///
/// # Safety
/// `options` must be null or point to a valid `GS_SESSION_start_o_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_SESSION_start(options: *const GS_SESSION_start_o_t) -> GS_ERROR_code_t {
    guard(|| {
        let precision = if options.is_null() {
            Precision::DEFAULT
        } else {
            // SAFETY: non-null and valid by the caller's contract.
            let o = unsafe { *options };
            if o.o_t_version != GS_O_T_VERSION {
                return Err(GS_ERROR_bad_option_version);
            }
            Precision::new(o.linear_precision, o.angular_precision, o.size_box)
                .map_err(|e| code(&e))?
        };
        let mut lock = SESSION.lock().map_err(|_| GS_ERROR_panic)?;
        if lock.is_some() {
            return Err(GS_ERROR_session_already_started);
        }
        *lock = Some(Session::start(SessionOptions { precision }));
        Ok(())
    })
}

/// Stops the session and deletes all entities. Kernel allocations stay valid until freed.
#[unsafe(no_mangle)]
pub extern "C" fn GS_SESSION_stop() -> GS_ERROR_code_t {
    guard(|| {
        let mut lock = SESSION
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        lock.take().map(|_| ()).ok_or(GS_ERROR_session_not_started)
    })
}

/// Creates a solid block of size `x × y × z` centred on the basis location (world origin if
/// `basis_set` is null).
///
/// # Safety
/// `basis_set` must be null or valid; `body` must be valid for writes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_BODY_create_solid_block(
    x: f64,
    y: f64,
    z: f64,
    basis_set: *const GS_AXIS2_sf_t,
    body: *mut GS_BODY_t,
) -> GS_ERROR_code_t {
    guard(|| {
        let basis = (!basis_set.is_null()).then(|| {
            // SAFETY: non-null and valid by the caller's contract.
            let b = unsafe { *basis_set };
            Axis2 {
                location: b.location,
                axis: b.axis,
                ref_direction: b.ref_direction,
            }
        });
        let tag = with_session(|s| s.body_create_solid_block(x, y, z, basis.as_ref()))?;
        // SAFETY: forwarded caller contract.
        unsafe { write_out(body, tag) }
    })
}

/// Faces of a body.
///
/// # Safety
/// The output pointers must be valid for writes. Free `*faces` with `GS_MEMORY_free`; it is
/// null when `*n_faces` is zero.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_BODY_ask_faces(
    body: GS_BODY_t,
    n_faces: *mut i32,
    faces: *mut *mut GS_FACE_t,
) -> GS_ERROR_code_t {
    guard(|| {
        let tags = with_session(|s| s.body_ask_faces(body))?;
        // SAFETY: forwarded caller contract.
        unsafe { write_tags(&tags, n_faces, faces) }
    })
}

/// Edges of a body.
///
/// # Safety
/// The output pointers must be valid for writes. Free `*edges` with `GS_MEMORY_free`; it is
/// null when `*n_edges` is zero.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_BODY_ask_edges(
    body: GS_BODY_t,
    n_edges: *mut i32,
    edges: *mut *mut GS_EDGE_t,
) -> GS_ERROR_code_t {
    guard(|| {
        let tags = with_session(|s| s.body_ask_edges(body))?;
        // SAFETY: forwarded caller contract.
        unsafe { write_tags(&tags, n_edges, edges) }
    })
}

/// Vertices of a body.
///
/// # Safety
/// The output pointers must be valid for writes. Free `*vertices` with `GS_MEMORY_free`; it
/// is null when `*n_vertices` is zero.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_BODY_ask_vertices(
    body: GS_BODY_t,
    n_vertices: *mut i32,
    vertices: *mut *mut GS_VERTEX_t,
) -> GS_ERROR_code_t {
    guard(|| {
        let tags = with_session(|s| s.body_ask_vertices(body))?;
        // SAFETY: forwarded caller contract.
        unsafe { write_tags(&tags, n_vertices, vertices) }
    })
}

/// Checks a body. `*n_faults == 0` means the body is valid.
///
/// # Safety
/// The output pointers must be valid for writes. Free `*faults` with `GS_MEMORY_free`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_BODY_check(
    body: GS_BODY_t,
    n_faults: *mut i32,
    faults: *mut *mut GS_check_fault_t,
) -> GS_ERROR_code_t {
    guard(|| {
        if n_faults.is_null() || faults.is_null() {
            return Err(GS_ERROR_null_argument);
        }
        let found = with_session(|s| s.body_check(body))?;
        let items: Vec<_> = found
            .iter()
            .map(|f| GS_check_fault_t {
                state: f.state as i32,
                entity: f.entity,
            })
            .collect();
        let n = i32::try_from(items.len()).map_err(|_| GS_ERROR_capacity_exceeded)?;
        let array = alloc_array(&items)?;
        // SAFETY: both pointers are non-null and valid for writes by the caller's contract.
        unsafe {
            n_faults.write(n);
            faults.write(array);
        }
        Ok(())
    })
}

/// Volume of the solid regions of a body.
///
/// # Safety
/// `volume` must be valid for writes of one double.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_BODY_ask_volume(body: GS_BODY_t, volume: *mut f64) -> GS_ERROR_code_t {
    guard(|| {
        let value = with_session(|s| s.body_ask_volume(body))?;
        // SAFETY: forwarded caller contract.
        unsafe { write_out(volume, value) }
    })
}

/// Axis-aligned bounding box of a body.
///
/// # Safety
/// `min_corner` and `max_corner` must be valid for writes of three doubles each.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_BODY_ask_bounding_box(
    body: GS_BODY_t,
    min_corner: *mut [f64; 3],
    max_corner: *mut [f64; 3],
) -> GS_ERROR_code_t {
    guard(|| {
        let (min, max) = with_session(|s| s.body_ask_bounding_box(body))?;
        // SAFETY: forwarded caller contract.
        unsafe {
            write_out(min_corner, min)?;
            write_out(max_corner, max)
        }
    })
}

/// Position of a vertex, written to `point[0..3]`.
///
/// # Safety
/// `point` must be valid for writes of three doubles.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_VERTEX_ask_point(
    vertex: GS_VERTEX_t,
    point: *mut [f64; 3],
) -> GS_ERROR_code_t {
    guard(|| {
        let p = with_session(|s| s.vertex_ask_point(vertex))?;
        // SAFETY: forwarded caller contract.
        unsafe { write_out(point, p) }
    })
}

/// Releases memory returned by the kernel. Null is ignored.
///
/// # Safety
/// `memory` must be null or a pointer returned by this library and not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn GS_MEMORY_free(memory: *mut c_void) -> GS_ERROR_code_t {
    if memory.is_null() {
        return GS_ERROR_no_errors;
    }
    // SAFETY: `memory` came from `alloc_array`, so the size header sits `HEADER` bytes
    // before it and describes the original layout.
    unsafe {
        let base = memory.cast::<u8>().sub(HEADER);
        let bytes = base.cast::<usize>().read();
        dealloc(base, Layout::from_size_align_unchecked(bytes, HEADER));
    }
    GS_ERROR_no_errors
}

#[cfg(test)]
mod tests {
    use gensolid_api::{ErrorCode, FaultKind};

    use super::*;

    #[test]
    fn codes_match_kernel_enums() {
        let errors = [
            (GS_ERROR_invalid_argument, ErrorCode::InvalidArgument),
            (GS_ERROR_not_found, ErrorCode::NotFound),
            (GS_ERROR_wrong_entity_class, ErrorCode::WrongEntityClass),
            (GS_ERROR_degenerate, ErrorCode::Degenerate),
            (GS_ERROR_out_of_range, ErrorCode::OutOfRange),
            (GS_ERROR_unsupported, ErrorCode::Unsupported),
            (GS_ERROR_unimplemented, ErrorCode::Unimplemented),
            (GS_ERROR_capacity_exceeded, ErrorCode::CapacityExceeded),
            (GS_ERROR_parse, ErrorCode::Parse),
            (GS_ERROR_internal, ErrorCode::Internal),
        ];
        for (c, e) in errors {
            assert_eq!(c, e as i32, "{e:?}");
        }
        let states = [
            (GS_BODY_state_missing_entity, FaultKind::MissingEntity),
            (GS_BODY_state_back_pointer, FaultKind::BackPointer),
            (GS_BODY_state_region_layout, FaultKind::RegionLayout),
            (GS_BODY_state_empty_loop, FaultKind::EmptyLoop),
            (GS_BODY_state_open_loop, FaultKind::OpenLoop),
            (GS_BODY_state_non_manifold_edge, FaultKind::NonManifoldEdge),
            (GS_BODY_state_fin_senses, FaultKind::FinSenses),
            (GS_BODY_state_missing_geometry, FaultKind::MissingGeometry),
            (GS_BODY_state_vertex_off_curve, FaultKind::VertexOffCurve),
            (
                GS_BODY_state_vertex_off_surface,
                FaultKind::VertexOffSurface,
            ),
            (GS_BODY_state_edge_off_surface, FaultKind::EdgeOffSurface),
        ];
        for (c, k) in states {
            assert_eq!(c, k as i32, "{k:?}");
        }
    }

    // One test: the session is global.
    #[test]
    #[allow(clippy::float_cmp)] // Half-sizes of the block are exactly representable.
    fn c_api_session_flow() {
        unsafe {
            let mut body = 0;
            assert_eq!(
                GS_BODY_create_solid_block(1.0, 1.0, 1.0, ptr::null(), &mut body),
                GS_ERROR_session_not_started
            );

            let mut options = std::mem::zeroed::<GS_SESSION_start_o_t>();
            assert_eq!(GS_SESSION_start_o_init(&mut options), GS_ERROR_no_errors);
            assert_eq!(GS_SESSION_start(&options), GS_ERROR_no_errors);
            assert_eq!(
                GS_SESSION_start(ptr::null()),
                GS_ERROR_session_already_started
            );

            assert_eq!(
                GS_BODY_create_solid_block(1.0, 2.0, 3.0, ptr::null(), &mut body),
                GS_ERROR_no_errors
            );
            assert_ne!(body, 0);
            assert_eq!(
                GS_BODY_create_solid_block(-1.0, 2.0, 3.0, ptr::null(), &mut body),
                GS_ERROR_invalid_argument
            );

            let (mut n, mut faces) = (0, ptr::null_mut());
            assert_eq!(
                GS_BODY_ask_faces(body, &mut n, &mut faces),
                GS_ERROR_no_errors
            );
            assert_eq!(n, 6);
            let face_tags = std::slice::from_raw_parts(faces, n as usize).to_vec();
            GS_MEMORY_free(faces.cast());
            assert_eq!(
                GS_BODY_ask_faces(face_tags[0], &mut n, &mut faces),
                GS_ERROR_wrong_entity_class
            );

            let (mut n_faults, mut faults) = (-1, ptr::null_mut());
            assert_eq!(
                GS_BODY_check(body, &mut n_faults, &mut faults),
                GS_ERROR_no_errors
            );
            assert_eq!(n_faults, 0);
            assert!(faults.is_null());

            let (mut nv, mut vertices) = (0, ptr::null_mut());
            assert_eq!(
                GS_BODY_ask_vertices(body, &mut nv, &mut vertices),
                GS_ERROR_no_errors
            );
            let mut point = [0.0; 3];
            assert_eq!(
                GS_VERTEX_ask_point(*vertices, &mut point),
                GS_ERROR_no_errors
            );
            assert_eq!(point, [-0.5, -1.0, -1.5]);
            let mut volume = 0.0;
            assert_eq!(GS_BODY_ask_volume(body, &mut volume), GS_ERROR_no_errors);
            assert_eq!(volume, 6.0);
            let (mut min_corner, mut max_corner) = ([0.0; 3], [0.0; 3]);
            assert_eq!(
                GS_BODY_ask_bounding_box(body, &mut min_corner, &mut max_corner),
                GS_ERROR_no_errors
            );
            assert_eq!(min_corner, [-0.5, -1.0, -1.5]);
            assert_eq!(max_corner, [0.5, 1.0, 1.5]);
            GS_MEMORY_free(vertices.cast());

            assert_eq!(
                GS_BODY_ask_faces(body, ptr::null_mut(), &mut faces),
                GS_ERROR_null_argument
            );
            assert_eq!(GS_SESSION_stop(), GS_ERROR_no_errors);
            assert_eq!(GS_SESSION_stop(), GS_ERROR_session_not_started);
        }
    }
}

/// Single-module browser ABI; isolated from the public native C API.
///
/// cbindgen:ignore
#[cfg(all(target_os = "emscripten", feature = "smlib"))]
mod single_wasm;

#[cfg(all(target_os = "emscripten", feature = "smlib-validation"))]
extern crate self as gensolid;
/// Shared geometry assertions compiled only into diagnostic Wasm packages.
///
/// cbindgen:ignore
#[cfg(all(target_os = "emscripten", feature = "smlib-validation"))]
#[path = "../tests/smlib.rs"]
mod single_validation;
