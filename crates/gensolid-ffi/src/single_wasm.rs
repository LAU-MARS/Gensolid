//! Private single-Wasm browser ABI. Rust calls SMLib directly in the same memory.
use crate::smlib::{Backend, Curve, Surface};
use gensolid_api::{Axis2, Session, SessionOptions};
use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_math::Point3;
use std::{cell::RefCell, collections::BTreeMap};

enum Entity {
    Session(Box<Session>),
    Backend(Backend),
    Curve(Curve),
    Surface(Surface),
}
#[derive(Default)]
struct Arena {
    next: usize,
    entities: BTreeMap<usize, Entity>,
}
thread_local! { static ARENA: RefCell<Arena> = RefCell::new(Arena::default()); }
impl Arena {
    fn insert(&mut self, entity: Entity) -> Result<Vec<f64>> {
        self.next = self
            .next
            .checked_add(1)
            .ok_or_else(|| Error::new(ErrorCode::CapacityExceeded, "browser handles exhausted"))?;
        self.entities.insert(self.next, entity);
        Ok(vec![self.next as f64])
    }
    fn get(&self, id: usize) -> Result<&Entity> {
        self.entities.get(&id).ok_or_else(|| {
            Error::new(
                ErrorCode::NotFound,
                "browser handle was released or does not exist",
            )
        })
    }
    fn backend(&self, id: usize) -> Result<Backend> {
        match self.get(id)? {
            Entity::Backend(v) => Ok(v.clone()),
            _ => Err(wrong()),
        }
    }
    fn curve(&self, id: usize) -> Result<Curve> {
        match self.get(id)? {
            Entity::Curve(v) => Ok(v.clone()),
            _ => Err(wrong()),
        }
    }
    fn surface(&self, id: usize) -> Result<Surface> {
        match self.get(id)? {
            Entity::Surface(v) => Ok(v.clone()),
            _ => Err(wrong()),
        }
    }
    fn session(&mut self, id: usize) -> Result<&mut Session> {
        match self.entities.get_mut(&id) {
            Some(Entity::Session(v)) => Ok(v),
            Some(_) => Err(wrong()),
            None => Err(Error::new(
                ErrorCode::NotFound,
                "browser session does not exist",
            )),
        }
    }
}
fn wrong() -> Error {
    Error::new(ErrorCode::WrongEntityClass, "wrong browser handle class")
}
fn number(values: &[f64], i: usize) -> Result<f64> {
    values
        .get(i)
        .copied()
        .ok_or_else(|| Error::invalid_argument("missing browser argument"))
}
fn point(values: &[f64], start: usize) -> Result<Point3> {
    Ok(Point3::new(
        number(values, start)?,
        number(values, start + 1)?,
        number(values, start + 2)?,
    ))
}
fn triple(values: &[f64], start: usize) -> Result<[f64; 3]> {
    let p = point(values, start)?;
    Ok([p.x, p.y, p.z])
}
fn tag(values: &[f64]) -> Result<i32> {
    let v = number(values, 0)?;
    if v < f64::from(i32::MIN) || v > f64::from(i32::MAX) || v.fract().abs() > 0. {
        return Err(Error::invalid_argument(
            "entity tag must be a signed 32-bit integer",
        ));
    }
    Ok(v as i32)
}
fn dispatch(arena: &mut Arena, op: u32, id: usize, v: &[f64]) -> Result<Vec<f64>> {
    if v.iter().any(|x| !x.is_finite()) {
        return Err(Error::invalid_argument("browser inputs must be finite"));
    }
    // All geometry calls below enter the existing safe Rust API, then the C ABI.
    match op {
        0 => {
            arena.entities.remove(&id).ok_or_else(|| {
                Error::new(ErrorCode::NotFound, "browser handle already released")
            })?;
            Ok(vec![])
        }
        1 => arena.insert(Entity::Session(Box::new(Session::start(
            SessionOptions::default(),
        )))),
        2 => arena.insert(Entity::Backend(Backend::new(Precision::DEFAULT)?)),
        10 => Ok(vec![f64::from(
            arena.session(id)?.body_create_solid_block(
                number(v, 0)?,
                number(v, 1)?,
                number(v, 2)?,
                None,
            )?,
        )]),
        11 => {
            let frame = Axis2 {
                location: triple(v, 3)?,
                axis: triple(v, 6)?,
                ref_direction: triple(v, 9)?,
            };
            Ok(vec![f64::from(
                arena.session(id)?.body_create_solid_block(
                    number(v, 0)?,
                    number(v, 1)?,
                    number(v, 2)?,
                    Some(&frame),
                )?,
            )])
        }
        12 => Ok(arena
            .session(id)?
            .body_ask_faces(tag(v)?)?
            .into_iter()
            .map(f64::from)
            .collect()),
        13 => Ok(arena
            .session(id)?
            .body_ask_edges(tag(v)?)?
            .into_iter()
            .map(f64::from)
            .collect()),
        14 => Ok(arena
            .session(id)?
            .body_ask_vertices(tag(v)?)?
            .into_iter()
            .map(f64::from)
            .collect()),
        15 => Ok(vec![arena.session(id)?.body_check(tag(v)?)?.len() as f64]),
        16 => Ok(vec![arena.session(id)?.body_ask_volume(tag(v)?)?]),
        17 => {
            let (a, b) = arena.session(id)?.body_ask_bounding_box(tag(v)?)?;
            Ok(a.into_iter().chain(b).collect())
        }
        18 => Ok(arena.session(id)?.vertex_ask_point(tag(v)?)?.to_vec()),
        20 => {
            let c = arena.backend(id)?.circle(point(v, 0)?, number(v, 3)?)?;
            arena.insert(Entity::Curve(c))
        }
        21 => {
            let c = arena
                .backend(id)?
                .line_segment(point(v, 0)?, point(v, 3)?)?;
            arena.insert(Entity::Curve(c))
        }
        22 => {
            let s = arena.backend(id)?.corner_surface([
                point(v, 0)?,
                point(v, 3)?,
                point(v, 6)?,
                point(v, 9)?,
            ])?;
            arena.insert(Entity::Surface(s))
        }
        30 => Ok(arena.curve(id)?.domain()?.to_vec()),
        31 => {
            let e = arena.curve(id)?.evaluate(number(v, 0)?)?;
            Ok(vec![
                e.point.x, e.point.y, e.point.z, e.d1.x, e.d1.y, e.d1.z, e.d2.x, e.d2.y, e.d2.z,
            ])
        }
        32 => {
            let c = arena.curve(id)?.reversed()?;
            arena.insert(Entity::Curve(c))
        }
        40 => Ok(arena.surface(id)?.domain()?.to_vec()),
        41 => {
            let e = arena
                .surface(id)?
                .evaluate([number(v, 0)?, number(v, 1)?])?;
            Ok(vec![
                e.point.x, e.point.y, e.point.z, e.du.x, e.du.y, e.du.z, e.dv.x, e.dv.y, e.dv.z,
                e.duu.x, e.duu.y, e.duu.z, e.duv.x, e.duv.y, e.duv.z, e.dvv.x, e.dvv.y, e.dvv.z,
            ])
        }
        #[cfg(feature = "smlib-validation")]
        100 => Ok(vec![
            crate::single_validation::run_smlib_validation()? as f64
        ]),
        _ => Err(Error::new(
            ErrorCode::Unimplemented,
            "unknown single-Wasm operation",
        )),
    }
}
/// Owned reply for the private browser ABI.
pub struct Response {
    code: i32,
    values: Vec<f64>,
    message: String,
}

/// Calls the private browser dispatcher. Returned response must be freed once.
/// # Safety
/// `values` must point to `len` readable, aligned doubles, or be null when len=0.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gs_web_call(
    op: u32,
    id: usize,
    values: *const f64,
    len: usize,
) -> *mut Response {
    let result = if len > isize::MAX as usize / size_of::<f64>()
        || (len != 0 && (values.is_null() || !values.is_aligned()))
    {
        Err(Error::invalid_argument("invalid browser input span"))
    } else {
        // SAFETY: the private JS caller owns this allocation for the synchronous call.
        let input = if len == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(values, len) }
        };
        ARENA.with(|a| match a.try_borrow_mut() {
            Ok(mut arena) => dispatch(&mut arena, op, id, input),
            Err(_) => Err(Error::internal("reentrant browser operation")),
        })
    };
    let response = match result {
        Ok(values) => Response {
            code: 0,
            values,
            message: String::new(),
        },
        Err(e) => Response {
            code: e.code() as i32,
            values: vec![],
            message: e.to_string(),
        },
    };
    Box::into_raw(Box::new(response))
}
/// # Safety
/// `p` is a live response returned by `gs_web_call`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gs_web_code(p: *const Response) -> i32 {
    // SAFETY: the caller retains the response until after these accessor calls.
    unsafe { (*p).code }
}
/// # Safety
/// `p` is a live response returned by `gs_web_call`; pointer expires when freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gs_web_values(p: *const Response) -> *const f64 {
    // SAFETY: the caller retains the response until after copying the values.
    unsafe { (*p).values.as_ptr() }
}
/// # Safety
/// `p` is a live response returned by `gs_web_call`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gs_web_values_len(p: *const Response) -> usize {
    // SAFETY: the caller retains the response until after copying the values.
    unsafe { (*p).values.len() }
}
/// # Safety
/// `p` is a live response returned by `gs_web_call`; pointer expires when freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gs_web_message(p: *const Response) -> *const u8 {
    // SAFETY: the caller retains the response until after copying the message.
    unsafe { (*p).message.as_ptr() }
}
/// # Safety
/// `p` is a live response returned by `gs_web_call`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gs_web_message_len(p: *const Response) -> usize {
    // SAFETY: the caller retains the response until after copying the message.
    let response = unsafe { &*p };
    response.message.len()
}
/// # Safety
/// `p` must be null or an unfreed pointer returned by `gs_web_call`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gs_web_free(p: *mut Response) {
    if !p.is_null() {
        // SAFETY: each response is returned exactly once to its allocating Rust module.
        unsafe {
            drop(Box::from_raw(p));
        }
    }
}
