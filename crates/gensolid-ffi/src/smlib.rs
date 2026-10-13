//! Optional, serialized SMLib curve/surface backend (enable `smlib`).
//!
//! Native objects stay in a reference-counted session arena until its last handle
//! is dropped. Handles are deliberately neither Send nor Sync. Operations that
//! edit geometry return copies, so failure cannot partially mutate Rust inputs.
//! Upstream numerical results are not Gensolid intersection certificates or bodies.
//! See `docs/smlib.md` for the coverage and upstream tolerance limitations.

use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_geom::{BCurve, BSurface};
use gensolid_math::{Point3, Vec3};
use std::{ffi::c_void, ptr::NonNull, rc::Rc};

/// Pinned upstream revision used by the native build.
pub const UPSTREAM_REVISION: &str = "c9c5979f04cea8c38213004d4f04df0dcfb40d05";

#[cfg(not(target_arch = "wasm32"))]
#[repr(C)]
struct RawResult {
    owner: *mut c_void,
    values: *const f64,
    values_len: usize,
    indices: *const usize,
    indices_len: usize,
    upstream_status: i64,
}
#[cfg(not(target_arch = "wasm32"))]
unsafe extern "C" {
    fn gs_sm_create(linear: f64, angular: f64, size: f64, session: *mut *mut c_void) -> i32;
    fn gs_sm_destroy(session: *mut c_void);
    fn gs_sm_call(
        session: *mut c_void,
        op: u32,
        ids: *const usize,
        ni: usize,
        values: *const f64,
        nv: usize,
        sizes: *const usize,
        ns: usize,
        result: *mut RawResult,
    ) -> i32;
    fn gs_sm_result_free(owner: *mut c_void);
}
#[cfg(not(target_arch = "wasm32"))]
impl Drop for RawResult {
    fn drop(&mut self) {
        // SAFETY: owner is null or the unique response returned by gs_sm_call.
        unsafe { gs_sm_result_free(self.owner) };
    }
}
struct Inner {
    ptr: NonNull<c_void>,
    precision: Precision,
}
impl Drop for Inner {
    fn drop(&mut self) {
        // SAFETY: the last Rc owns the live native session; all handles are gone.
        #[cfg(not(target_arch = "wasm32"))]
        unsafe {
            gs_sm_destroy(self.ptr.as_ptr())
        };
        #[cfg(target_arch = "wasm32")]
        browser::destroy(self.ptr.as_ptr() as usize);
    }
}
struct Output {
    values: Vec<f64>,
    indices: Vec<usize>,
}
fn native_error(code: i32, upstream: i64, operation: u32) -> Error {
    let category = match code {
        1 => ErrorCode::InvalidArgument,
        3 => ErrorCode::WrongEntityClass,
        4 => ErrorCode::Degenerate,
        6 => ErrorCode::Unsupported,
        7 => ErrorCode::Unimplemented,
        8 => ErrorCode::CapacityExceeded,
        _ => ErrorCode::Internal,
    };
    Error::new(
        category,
        format!("SMLib operation {operation}: bridge status {code}, upstream status {upstream}"),
    )
}
impl Inner {
    #[cfg(not(target_arch = "wasm32"))]
    fn call(&self, op: u32, ids: &[usize], values: &[f64], sizes: &[usize]) -> Result<Output> {
        if values.iter().any(|x| !x.is_finite()) {
            return Err(Error::invalid_argument("SMLib inputs must be finite"));
        }
        let mut raw = RawResult {
            owner: std::ptr::null_mut(),
            values: std::ptr::null(),
            values_len: 0,
            indices: std::ptr::null(),
            indices_len: 0,
            upstream_status: 0,
        };
        // SAFETY: self owns a live session; spans remain borrowed through this
        // synchronous call. The native adapter catches exceptions and serializes access.
        let code = unsafe {
            gs_sm_call(
                self.ptr.as_ptr(),
                op,
                ids.as_ptr(),
                ids.len(),
                values.as_ptr(),
                values.len(),
                sizes.as_ptr(),
                sizes.len(),
                &mut raw,
            )
        };
        if code != 0 {
            return Err(native_error(code, raw.upstream_status, op));
        }
        // SAFETY: the private ABI supplies initialized spans owned by raw.owner;
        // raw remains alive until both spans are copied. Empty spans may be null.
        let values = if raw.values_len == 0 {
            Vec::new()
        } else {
            unsafe { std::slice::from_raw_parts(raw.values, raw.values_len) }.to_vec()
        };
        // SAFETY: same ownership and span contract as values above.
        let indices = if raw.indices_len == 0 {
            Vec::new()
        } else {
            unsafe { std::slice::from_raw_parts(raw.indices, raw.indices_len) }.to_vec()
        };
        if values.iter().any(|v| !v.is_finite()) {
            return Err(Error::internal("SMLib returned non-finite output"));
        }
        Ok(Output { values, indices })
    }
}
#[cfg(target_arch = "wasm32")]
#[path = "smlib_browser.rs"]
mod browser;
#[cfg(target_arch = "wasm32")]
pub use browser::initialize;
#[cfg(target_arch = "wasm32")]
impl Inner {
    fn call(&self, op: u32, ids: &[usize], values: &[f64], sizes: &[usize]) -> Result<Output> {
        if values.iter().any(|x| !x.is_finite()) {
            return Err(Error::invalid_argument("SMLib inputs must be finite"));
        }
        browser::call(self.ptr.as_ptr() as usize, op, ids, values, sizes)
    }
}

fn fixed<const N: usize>(v: &[f64]) -> Result<[f64; N]> {
    v.try_into()
        .map_err(|_| Error::internal("unexpected SMLib output length"))
}
fn first(o: &Output) -> Result<usize> {
    o.indices
        .first()
        .copied()
        .ok_or_else(|| Error::internal("SMLib returned no geometry"))
}
fn xyz(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn vector(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn points(p: &[Point3]) -> Vec<f64> {
    p.iter().flat_map(|&p| xyz(p)).collect()
}
fn positive(value: f64, precision: Precision) -> Result<()> {
    if !value.is_finite() || value <= precision.degenerate_length() {
        Err(Error::degenerate(
            "dimension must exceed session linear precision",
        ))
    } else {
        Ok(())
    }
}

/// An independent arena of native geometry; cloning shares ownership.
#[derive(Clone)]
pub struct Backend {
    inner: Rc<Inner>,
}
/// Owned handle to a curve in a native session. Edits produce new handles.
#[derive(Clone)]
pub struct Curve {
    inner: Rc<Inner>,
    id: usize,
}
/// Owned handle to a surface in a native session. Edits produce new handles.
#[derive(Clone)]
pub struct Surface {
    inner: Rc<Inner>,
    id: usize,
}
/// Curve position and first/second derivatives returned by SMLib.
#[derive(Clone, Copy, Debug)]
pub struct CurveEvaluation {
    /// Position.
    pub point: Point3,
    /// First derivative.
    pub d1: Vec3,
    /// Second derivative.
    pub d2: Vec3,
}
/// Surface position and first/second partial derivatives returned by SMLib.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceEvaluation {
    /// Position.
    pub point: Point3,
    /// First U derivative.
    pub du: Vec3,
    /// First V derivative.
    pub dv: Vec3,
    /// Second U derivative.
    pub duu: Vec3,
    /// Mixed U/V derivative.
    pub duv: Vec3,
    /// Second V derivative.
    pub dvv: Vec3,
}
/// Parameterization used by the interpolating curve fitter.
#[derive(Clone, Copy, Debug)]
#[repr(usize)]
pub enum Parameterization {
    /// Equally spaced parameters.
    Uniform = 0,
    /// Parameters proportional to chord lengths.
    ChordLength = 1,
    /// Parameters proportional to square roots of chord lengths.
    Centripetal = 2,
}
/// Parameter direction of a surface.
#[derive(Clone, Copy, Debug)]
pub enum Direction {
    /// First surface parameter.
    U,
    /// Second surface parameter.
    V,
}
/// A native solver result, including intervals for overlaps.
/// No global completeness certificate is implied by this type.
#[derive(Clone, Debug)]
pub struct Solution<const D: usize> {
    /// Upstream SmSolutionType value, retained without losing unknown variants.
    pub native_kind: usize,
    /// Start parameters (curve/curve: t1,t2; curve/surface: t,u,v).
    pub start: [f64; D],
    /// End parameters, meaningful for range solutions.
    pub end: [f64; D],
}
fn solutions<const D: usize>(o: Output) -> Result<Vec<Solution<D>>> {
    if o.values.len() != o.indices.len() * D * 2 {
        return Err(Error::internal("invalid solution layout"));
    }
    o.values
        .chunks_exact(D * 2)
        .zip(o.indices)
        .map(|(v, kind)| {
            Ok(Solution {
                native_kind: kind,
                start: fixed(&v[..D])?,
                end: fixed(&v[D..])?,
            })
        })
        .collect()
}
/// One numerical surface intersection branch, with paired parameter curves.
pub struct IntersectionBranch {
    /// Three-dimensional curve.
    pub curve: Curve,
    /// Parameter curve on the first surface (x=u, y=v).
    pub first_uv: Curve,
    /// Parameter curve on the second surface (x=u, y=v).
    pub second_uv: Curve,
    /// Upstream SmTsectCurveType value.
    pub native_kind: usize,
    /// Upstream-reported deviation; not a certified error bound.
    pub reported_deviation: f64,
}
/// Curves resulting from projection/drop onto a surface.
pub struct ProjectedCurves {
    /// Curves in the surface parameter space.
    pub uv: Vec<Curve>,
    /// Curves in world coordinates.
    pub spatial: Vec<Curve>,
}
impl Backend {
    /// Creates a session. All native entry points are serialized process-wide.
    pub fn new(precision: Precision) -> Result<Self> {
        #[cfg(target_arch = "wasm32")]
        let ptr = browser::create(precision)? as *mut c_void;
        #[cfg(not(target_arch = "wasm32"))]
        let ptr = {
            let mut ptr = std::ptr::null_mut();
            // SAFETY: ptr is a writable output; precision is validated by its constructor.
            let code = unsafe {
                gs_sm_create(
                    precision.linear(),
                    precision.angular(),
                    precision.size(),
                    &mut ptr,
                )
            };
            if code != 0 {
                return Err(native_error(code, 0, 0));
            }
            ptr
        };
        let ptr =
            NonNull::new(ptr).ok_or_else(|| Error::internal("SMLib returned a null session"))?;
        Ok(Self {
            inner: Rc::new(Inner { ptr, precision }),
        })
    }
    /// Session precision passed to low-level intersections and point projections.
    pub fn precision(&self) -> Precision {
        self.inner.precision
    }
    fn curve(&self, op: u32, ids: &[usize], d: &[f64], n: &[usize]) -> Result<Curve> {
        let o = self.inner.call(op, ids, d, n)?;
        Ok(Curve {
            inner: self.inner.clone(),
            id: first(&o)?,
        })
    }
    fn surface(&self, op: u32, ids: &[usize], d: &[f64], n: &[usize]) -> Result<Surface> {
        let o = self.inner.call(op, ids, d, n)?;
        Ok(Surface {
            inner: self.inner.clone(),
            id: first(&o)?,
        })
    }
    fn curves(&self, op: u32, ids: &[usize], d: &[f64], n: &[usize]) -> Result<Vec<Curve>> {
        Ok(self
            .inner
            .call(op, ids, d, n)?
            .indices
            .into_iter()
            .map(|id| Curve {
                inner: self.inner.clone(),
                id,
            })
            .collect())
    }
    fn same(&self, other: &Rc<Inner>) -> Result<()> {
        if Rc::ptr_eq(&self.inner, other) {
            Ok(())
        } else {
            Err(Error::invalid_argument(
                "geometry belongs to another SMLib session",
            ))
        }
    }
    fn curve_ids(&self, curves: &[Curve]) -> Result<Vec<usize>> {
        curves
            .iter()
            .map(|c| {
                self.same(&c.inner)?;
                Ok(c.id)
            })
            .collect()
    }
    /// Imports validated Gensolid NURBS data, preserving rational control points.
    pub fn import_curve(&self, curve: &BCurve) -> Result<Curve> {
        let mut d = curve.knots().to_vec();
        for (i, p) in curve.points().iter().enumerate() {
            let w = curve.weights().map_or(1.0, |w| w[i]);
            d.extend([p.x, p.y, p.z, w]);
        }
        self.curve(
            1,
            &[],
            &d,
            &[curve.degree(), curve.points().len(), curve.knots().len()],
        )
    }
    /// Imports a row-major (u then v) Gensolid NURBS surface.
    pub fn import_surface(&self, surface: &BSurface) -> Result<Surface> {
        let (u, v) = surface.knots();
        let (pu, pv) = surface.degrees();
        let (nu, nv) = surface.counts();
        let mut d = u.to_vec();
        d.extend(v);
        for (i, p) in surface.points().iter().enumerate() {
            let w = surface.weights().map_or(1.0, |w| w[i]);
            d.extend([p.x, p.y, p.z, w]);
        }
        self.surface(2, &[], &d, &[pu, pv, nu, nv, u.len(), v.len()])
    }
    /// Finite analytic line segment.
    pub fn line_segment(&self, start: Point3, end: Point3) -> Result<Curve> {
        positive(start.distance(end), self.precision())?;
        self.curve(20, &[], &[xyz(start), xyz(end)].concat(), &[])
    }
    /// Rational circle in the XY plane.
    pub fn circle(&self, center: Point3, radius: f64) -> Result<Curve> {
        positive(radius, self.precision())?;
        let mut d = xyz(center).to_vec();
        d.push(radius);
        self.curve(21, &[], &d, &[])
    }
    /// Analytic ellipse in the XY plane.
    pub fn ellipse(&self, center: Point3, radii: [f64; 2]) -> Result<Curve> {
        for r in radii {
            positive(r, self.precision())?;
        }
        let mut d = xyz(center).to_vec();
        d.extend(radii);
        self.curve(22, &[], &d, &[])
    }
    /// Rational circular arc in the XY plane; angles are in radians.
    pub fn arc(&self, center: Point3, radius: f64, angles: [f64; 2]) -> Result<Curve> {
        positive(radius, self.precision())?;
        let mut d = xyz(center).to_vec();
        d.extend([radius, angles[0].to_degrees(), angles[1].to_degrees()]);
        self.curve(23, &[], &d, &[])
    }
    /// Helix along +Z, fitted using session linear precision. Upstream may reject
    /// requested tolerances outside its supported fitting range.
    pub fn helix(
        &self,
        origin: Point3,
        height: f64,
        radii: [f64; 2],
        turns: f64,
        right_handed: bool,
    ) -> Result<Curve> {
        for x in [height, radii[0], radii[1]] {
            positive(x, self.precision())?;
        }
        let mut d = xyz(origin).to_vec();
        d.extend([height, radii[0], radii[1], turns]);
        self.curve(24, &[], &d, &[usize::from(right_handed)])
    }
    /// Cubic interpolation through at least four points.
    pub fn interpolate_curve(
        &self,
        input: &[Point3],
        parameterization: Parameterization,
    ) -> Result<Curve> {
        at_least(input.len(), 4)?;
        self.curve(25, &[], &points(input), &[parameterization as usize])
    }
    /// Upstream curve approximation; uses upstream fitting policy, not a
    /// certified session-precision error bound.
    pub fn approximate_curve(&self, input: &[Point3]) -> Result<Curve> {
        at_least(input.len(), 4)?;
        self.curve(26, &[], &points(input), &[])
    }
    /// Cubic curve from at least four control points and upstream-generated knots.
    pub fn control_curve(&self, input: &[Point3]) -> Result<Curve> {
        at_least(input.len(), 4)?;
        self.curve(27, &[], &points(input), &[])
    }
    /// Four NURBS sides of an XY rectangle centered at the supplied point.
    pub fn rectangle(&self, center: Point3, dimensions: [f64; 2]) -> Result<Vec<Curve>> {
        for x in dimensions {
            positive(x, self.precision())?;
        }
        let mut d = xyz(center).to_vec();
        d.extend(dimensions);
        self.curves(28, &[], &d, &[])
    }
    /// NURBS sides of a regular XY polygon.
    pub fn regular_polygon(&self, center: Point3, sides: usize, radius: f64) -> Result<Vec<Curve>> {
        at_least(sides, 3)?;
        if sides > i32::MAX as usize {
            return Err(Error::new(
                ErrorCode::CapacityExceeded,
                "polygon side count exceeds native indexing",
            ));
        }
        positive(radius, self.precision())?;
        let mut d = xyz(center).to_vec();
        d.push(radius);
        self.curves(29, &[], &d, &[sides])
    }
    /// Bilinear patch: corners ordered (umin,vmin), (umax,vmin),
    /// (umin,vmax), (umax,vmax).
    pub fn corner_surface(&self, corners: [Point3; 4]) -> Result<Surface> {
        self.surface(30, &[], &points(&corners), &[])
    }
    /// Bicubic surface from a control grid, at least four rows and columns.
    pub fn control_surface(
        &self,
        input: &[Point3],
        rows: usize,
        columns: usize,
    ) -> Result<Surface> {
        grid(input, rows, columns)?;
        self.surface(31, &[], &points(input), &[rows, columns])
    }
    /// Interpolating surface from an ordered point grid.
    pub fn interpolate_surface(
        &self,
        input: &[Point3],
        rows: usize,
        columns: usize,
    ) -> Result<Surface> {
        grid(input, rows, columns)?;
        self.surface(32, &[], &points(input), &[rows, columns])
    }
    /// Fits a surface to an unordered point set using upstream fitting policy.
    pub fn approximate_surface(&self, input: &[Point3]) -> Result<Surface> {
        at_least(input.len(), 4)?;
        self.surface(33, &[], &points(input), &[])
    }
    /// Extrudes a NURBS profile by a vector, returning geometry only.
    pub fn extrude(&self, profile: &Curve, direction: Vec3) -> Result<Surface> {
        self.same(&profile.inner)?;
        positive(direction.length(), self.precision())?;
        self.surface(34, &[profile.id], &vector(direction), &[])
    }
    /// Sweeps a NURBS profile along a NURBS path.
    pub fn sweep(&self, profile: &Curve, path: &Curve) -> Result<Surface> {
        self.same(&profile.inner)?;
        self.same(&path.inner)?;
        self.surface(35, &[profile.id, path.id], &[], &[])
    }
    /// Ruled surface between two NURBS curves.
    pub fn ruled(&self, first: &Curve, second: &Curve) -> Result<Surface> {
        self.same(&first.inner)?;
        self.same(&second.inner)?;
        self.surface(36, &[first.id, second.id], &[], &[])
    }
    /// Revolves a NURBS profile; angle is in radians.
    pub fn revolve(
        &self,
        profile: &Curve,
        origin: Point3,
        axis: Vec3,
        angle: f64,
    ) -> Result<Surface> {
        self.same(&profile.inner)?;
        let mut d = xyz(origin).to_vec();
        let axis = axis
            .normalized(self.precision().angular())
            .ok_or_else(|| Error::degenerate("revolution axis is zero"))?;
        d.extend(vector(axis));
        d.push(angle.to_degrees());
        self.surface(37, &[profile.id], &d, &[])
    }
    /// Surface skin through at least two NURBS section curves.
    pub fn loft(&self, sections: &[Curve]) -> Result<Surface> {
        at_least(sections.len(), 2)?;
        self.surface(38, &self.curve_ids(sections)?, &[], &[])
    }
    /// Makes copies with compatible NURBS degrees and knot vectors.
    pub fn compatible_curves(&self, curves: &[Curve]) -> Result<Vec<Curve>> {
        at_least(curves.len(), 2)?;
        self.curves(44, &self.curve_ids(curves)?, &[], &[])
    }
    /// Orders copies using upstream loop ordering. Curves are not reversed and
    /// this does not certify that a closed, consistently oriented loop exists.
    pub fn order_curves(&self, curves: &[Curve]) -> Result<Vec<Curve>> {
        at_least(curves.len(), 1)?;
        self.curves(45, &self.curve_ids(curves)?, &[], &[])
    }
}
fn at_least(n: usize, min: usize) -> Result<()> {
    if n < min {
        Err(Error::invalid_argument(format!(
            "expected at least {min} items"
        )))
    } else {
        Ok(())
    }
}
fn grid(points: &[Point3], rows: usize, cols: usize) -> Result<()> {
    at_least(rows, 4)?;
    at_least(cols, 4)?;
    if rows.checked_mul(cols) != Some(points.len()) {
        return Err(Error::invalid_argument(
            "grid dimensions do not match point count",
        ));
    }
    Ok(())
}
impl Curve {
    fn backend(&self) -> Backend {
        Backend {
            inner: self.inner.clone(),
        }
    }
    /// Natural parameter interval.
    pub fn domain(&self) -> Result<[f64; 2]> {
        fixed(&self.inner.call(7, &[self.id], &[], &[])?.values)
    }
    /// Position, first derivative and second derivative without tangent replacement.
    pub fn evaluate(&self, t: f64) -> Result<CurveEvaluation> {
        let a = fixed::<9>(&self.inner.call(5, &[self.id], &[t], &[])?.values)?;
        Ok(CurveEvaluation {
            point: Point3::new(a[0], a[1], a[2]),
            d1: Vec3::new(a[3], a[4], a[5]),
            d2: Vec3::new(a[6], a[7], a[8]),
        })
    }
    /// Exports native NURBS data. Analytic curve handles require conversion first.
    pub fn nurbs(&self) -> Result<BCurve> {
        let o = self.inner.call(3, &[self.id], &[], &[])?;
        let [degree, count, knots]: [usize; 3] = o
            .indices
            .as_slice()
            .try_into()
            .map_err(|_| Error::internal("invalid curve dimensions"))?;
        if o.values.len() != knots + 4 * count {
            return Err(Error::internal("invalid curve data"));
        }
        let (p, w) = unpack_points(&o.values[knots..]);
        BCurve::new(degree, o.values[..knots].to_vec(), p, Some(w))
    }
    /// Numerical closest-point parameter solutions, including range solutions.
    pub fn project_point(&self, point: Point3) -> Result<Vec<Solution<1>>> {
        solutions(self.inner.call(9, &[self.id], &xyz(point), &[])?)
    }
    /// Numerical curve intersections, preserving overlap parameter intervals.
    pub fn intersect_curve(&self, other: &Curve) -> Result<Vec<Solution<2>>> {
        self.backend().same(&other.inner)?;
        solutions(self.inner.call(11, &[self.id, other.id], &[], &[])?)
    }
    /// Numerical curve/surface intersections; parameters are t,u,v.
    pub fn intersect_surface(&self, other: &Surface) -> Result<Vec<Solution<3>>> {
        self.backend().same(&other.inner)?;
        solutions(self.inner.call(12, &[self.id, other.id], &[], &[])?)
    }
    /// Offsets a planar NURBS curve according to upstream orientation conventions.
    pub fn offset(&self, distance: f64) -> Result<Curve> {
        self.backend().curve(40, &[self.id], &[distance], &[])
    }
    /// Offsets in an explicit plane; useful when a straight curve has no unique plane.
    pub fn offset_in_plane(&self, distance: f64, normal: Vec3) -> Result<Curve> {
        let n = normal
            .normalized(self.inner.precision.angular())
            .ok_or_else(|| Error::degenerate("offset plane normal is zero"))?;
        self.backend()
            .curve(40, &[self.id], &[distance, n.x, n.y, n.z], &[])
    }
    fn onto(&self, op: u32, surface: &Surface) -> Result<ProjectedCurves> {
        self.backend().same(&surface.inner)?;
        let o = self.inner.call(op, &[self.id, surface.id], &[], &[])?;
        let [count] = fixed::<1>(&o.values)?;
        let count = count as usize;
        if count > o.indices.len() {
            return Err(Error::internal("invalid projection layout"));
        }
        let mut curves: Vec<_> = o
            .indices
            .into_iter()
            .map(|id| Curve {
                inner: self.inner.clone(),
                id,
            })
            .collect();
        let spatial = curves.split_off(count);
        Ok(ProjectedCurves {
            uv: curves,
            spatial,
        })
    }
    /// Lifts a UV NURBS curve through a surface. UV coordinates are curve x,y.
    pub fn lift(&self, surface: &Surface) -> Result<Vec<Curve>> {
        Ok(self.onto(41, surface)?.spatial)
    }
    /// Drops a NURBS curve to a surface, returning UV and spatial pieces.
    pub fn drop_to(&self, surface: &Surface) -> Result<ProjectedCurves> {
        self.onto(42, surface)
    }
    /// Projects a NURBS curve to a surface using upstream projection semantics.
    pub fn project_to(&self, surface: &Surface) -> Result<ProjectedCurves> {
        self.onto(43, surface)
    }
    /// Removes extra knots from a copy using session linear tolerance.
    pub fn simplify(&self) -> Result<Curve> {
        self.backend().curve(46, &[self.id], &[], &[])
    }
    /// Converts a NURBS curve to line/arc pieces using session linear precision.
    pub fn lines_and_arcs(&self) -> Result<Vec<Curve>> {
        self.backend().curves(47, &[self.id], &[], &[])
    }
    /// Native SmContinuityType at a parameter.
    pub fn continuity(&self, t: f64) -> Result<usize> {
        first(&self.inner.call(48, &[self.id], &[t], &[])?)
    }
    /// Trims a NURBS copy; upstream may snap boundaries to existing knots.
    pub fn trim(&self, interval: [f64; 2]) -> Result<Curve> {
        validate_interval(interval, self.domain()?)?;
        self.backend().curve(49, &[self.id], &interval, &[])
    }
    /// Reverses a NURBS copy while preserving the natural parameter interval.
    pub fn reversed(&self) -> Result<Curve> {
        self.backend().curve(51, &[self.id], &[], &[])
    }
    /// Degree-elevates a NURBS copy.
    pub fn elevate_degree(&self, degree: usize) -> Result<Curve> {
        let original = self.nurbs()?.degree();
        if degree < original || degree > i16::MAX as usize {
            return Err(Error::invalid_argument("invalid elevated degree"));
        }
        self.backend().curve(53, &[self.id], &[], &[degree])
    }
    /// Inserts knots into a NURBS copy. Values must be within the natural domain.
    pub fn insert_knots(&self, knots: &[f64]) -> Result<Curve> {
        let spline = self.nurbs()?;
        let sorted = validated_insertions(spline.knots(), spline.degree(), knots)?;
        self.backend().curve(54, &[self.id], &sorted, &[])
    }
}
impl Surface {
    fn backend(&self) -> Backend {
        Backend {
            inner: self.inner.clone(),
        }
    }
    /// Natural domain as [umin,umax,vmin,vmax].
    pub fn domain(&self) -> Result<[f64; 4]> {
        fixed(&self.inner.call(8, &[self.id], &[], &[])?.values)
    }
    /// Position and all first/second partial derivatives.
    pub fn evaluate(&self, uv: [f64; 2]) -> Result<SurfaceEvaluation> {
        let a = fixed::<18>(&self.inner.call(6, &[self.id], &uv, &[])?.values)?;
        Ok(SurfaceEvaluation {
            point: Point3::new(a[0], a[1], a[2]),
            du: Vec3::new(a[3], a[4], a[5]),
            dv: Vec3::new(a[6], a[7], a[8]),
            duu: Vec3::new(a[9], a[10], a[11]),
            duv: Vec3::new(a[12], a[13], a[14]),
            dvv: Vec3::new(a[15], a[16], a[17]),
        })
    }
    /// Exports native NURBS surface data in Gensolid row-major order.
    pub fn nurbs(&self) -> Result<BSurface> {
        let o = self.inner.call(4, &[self.id], &[], &[])?;
        let [pu, pv, nu, nv, ku, kv]: [usize; 6] = o
            .indices
            .as_slice()
            .try_into()
            .map_err(|_| Error::internal("invalid surface dimensions"))?;
        if o.values.len() != ku + kv + 4 * nu * nv {
            return Err(Error::internal("invalid surface data"));
        }
        let (p, w) = unpack_points(&o.values[ku + kv..]);
        BSurface::new(
            pu,
            pv,
            o.values[..ku].to_vec(),
            o.values[ku..ku + kv].to_vec(),
            nu,
            nv,
            p,
            Some(w),
        )
    }
    /// Numerical closest-point UV solutions.
    pub fn project_point(&self, point: Point3) -> Result<Vec<Solution<2>>> {
        solutions(self.inner.call(10, &[self.id], &xyz(point), &[])?)
    }
    /// Numerical surface intersection branches. Empty success means that the
    /// upstream solver returned no curves, not a Gensolid certified empty set.
    pub fn intersect(&self, other: &Surface) -> Result<Vec<IntersectionBranch>> {
        self.backend().same(&other.inner)?;
        let o = self.inner.call(13, &[self.id, other.id], &[], &[])?;
        let n = o.values.len();
        if o.indices.len() != 4 * n {
            return Err(Error::internal("invalid surface intersection layout"));
        }
        Ok((0..n)
            .map(|i| IntersectionBranch {
                curve: Curve {
                    inner: self.inner.clone(),
                    id: o.indices[i],
                },
                first_uv: Curve {
                    inner: self.inner.clone(),
                    id: o.indices[n + i],
                },
                second_uv: Curve {
                    inner: self.inner.clone(),
                    id: o.indices[2 * n + i],
                },
                native_kind: o.indices[3 * n + i],
                reported_deviation: o.values[i],
            })
            .collect())
    }
    /// Approximate offset surface, requesting session linear precision.
    pub fn offset(&self, distance: f64) -> Result<Surface> {
        self.backend().surface(39, &[self.id], &[distance], &[])
    }
    /// Isoparametric curve with the selected coordinate held constant.
    pub fn iso_curve(&self, direction: Direction, parameter: f64) -> Result<Curve> {
        let domain = self.domain()?;
        let offset = 2 * direction.index();
        if !parameter.is_finite() || parameter < domain[offset] || parameter > domain[offset + 1] {
            return Err(Error::invalid_argument(
                "iso parameter lies outside the domain",
            ));
        }
        self.backend()
            .curve(55, &[self.id], &[parameter], &[direction.index()])
    }
    /// Degree-elevates a copy along one surface direction.
    pub fn elevate_degree(&self, direction: Direction, degree: usize) -> Result<Surface> {
        let (u, v) = self.nurbs()?.degrees();
        let original = match direction {
            Direction::U => u,
            Direction::V => v,
        };
        if degree < original || degree > i16::MAX as usize {
            return Err(Error::invalid_argument("invalid elevated degree"));
        }
        self.backend()
            .surface(56, &[self.id], &[], &[direction.index(), degree])
    }
    /// Inserts knots along one direction into a copy.
    pub fn insert_knots(&self, direction: Direction, knots: &[f64]) -> Result<Surface> {
        let spline = self.nurbs()?;
        let (u, v) = spline.knots();
        let (pu, pv) = spline.degrees();
        let (original, degree) = match direction {
            Direction::U => (u, pu),
            Direction::V => (v, pv),
        };
        let sorted = validated_insertions(original, degree, knots)?;
        self.backend()
            .surface(57, &[self.id], &sorted, &[direction.index()])
    }
    /// Splits a NURBS surface at an interior U or V parameter.
    pub fn split(&self, direction: Direction, parameter: f64) -> Result<(Surface, Surface)> {
        let domain = self.domain()?;
        let offset = 2 * direction.index();
        if !parameter.is_finite() || parameter <= domain[offset] || parameter >= domain[offset + 1]
        {
            return Err(Error::invalid_argument("split must be inside the domain"));
        }
        let o = self
            .inner
            .call(58, &[self.id], &[parameter], &[direction.index()])?;
        let [left, right]: [usize; 2] = o
            .indices
            .as_slice()
            .try_into()
            .map_err(|_| Error::internal("invalid split output"))?;
        Ok((
            Surface {
                inner: self.inner.clone(),
                id: left,
            },
            Surface {
                inner: self.inner.clone(),
                id: right,
            },
        ))
    }
    /// Rectangular parameter-domain trim of a NURBS copy; not a trimmed B-rep face.
    pub fn trim(&self, domain: [f64; 4]) -> Result<Surface> {
        let [u0, u1, v0, v1] = self.domain()?;
        validate_interval([domain[0], domain[1]], [u0, u1])?;
        validate_interval([domain[2], domain[3]], [v0, v1])?;
        self.backend().surface(50, &[self.id], &domain, &[])
    }
    /// Reverses one parameter direction of a NURBS copy.
    pub fn reversed(&self, direction: Direction) -> Result<Surface> {
        self.backend().surface(
            52,
            &[self.id],
            &[],
            &[match direction {
                Direction::U => 0,
                Direction::V => 1,
            }],
        )
    }
}
fn validate_interval([lo, hi]: [f64; 2], [min, max]: [f64; 2]) -> Result<()> {
    if !lo.is_finite() || !hi.is_finite() || lo >= hi || lo < min || hi > max {
        Err(Error::invalid_argument("invalid trim domain"))
    } else {
        Ok(())
    }
}
fn unpack_points(data: &[f64]) -> (Vec<Point3>, Vec<f64>) {
    data.as_chunks::<4>()
        .0
        .iter()
        .map(|a| (Point3::new(a[0], a[1], a[2]), a[3]))
        .unzip()
}

impl Direction {
    fn index(self) -> usize {
        match self {
            Self::U => 0,
            Self::V => 1,
        }
    }
}
// Knot equality here compares the NURBS representation, not spatial topology.
fn validated_insertions(original: &[f64], degree: usize, new: &[f64]) -> Result<Vec<f64>> {
    let lo = original[degree];
    let hi = original[original.len() - degree - 1];
    if new.iter().any(|&x| !x.is_finite() || x <= lo || x >= hi) {
        return Err(Error::invalid_argument("insertions must be interior knots"));
    }
    let mut sorted = new.to_vec();
    sorted.sort_by(f64::total_cmp);
    for &value in &sorted {
        let count = original
            .iter()
            .chain(&sorted)
            .filter(|&&x| x.total_cmp(&value).is_eq())
            .count();
        if count > degree {
            return Err(Error::invalid_argument("knot multiplicity exceeds degree"));
        }
    }
    Ok(sorted)
}
