//! Pure-Rust curve and surface backend powered by `solidcore` (Apache-2.0).
//!
//! Source provenance: adapts the published `solidcore` 1.0.1 crate
//! (<https://github.com/tangency-labs/solidcore>, commit `7267671f`), a safe-Rust port
//! of NVIDIA's SMLib/SM_API (<https://github.com/NVIDIA-Omniverse/usd-brep>, Apache-2.0).
//! See `NOTICE` for copyright and license notices.
//!
//! No `unsafe` code, C++ toolchain, companion WebAssembly module or process-global
//! mutex is used. Numerical results from `solidcore` remain candidate geometry until
//! validated by Gensolid's interval/predicate certificates and `gensolid-check`.

use std::cell::RefCell;
use std::rc::Rc;

use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_math::{Point3, Vec3};
use solidcore::curve::curve::{Curve as ScCurve, CurveKind as ScCurveKind};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::point_set::PointSet3d;
use solidcore::geom::solution_array::{ObjectKind, SolutionArray, SolutionObject, SolutionType};
use solidcore::geom::vector2d::Point2 as ScPoint2;
use solidcore::kernel::{Kernel, KernelOptions};
use solidcore::math::{Point3 as ScPoint3, Vec3 as ScVec3};
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::b_spline_curve_n_lib::TrimMode;
use solidcore::nurbs::b_spline_surface::{BSplineSurface, SurfParamType};
use solidcore::nurbs::crv_fit::p1::ParamMethod;
use solidcore::nurbs::math::{CPoint, NOW, NOZ};
use solidcore::nurbs::srf_fit::lst_sq_approx;
use solidcore::primitive::api_curves::{self, HelixParams};
use solidcore::primitive::{api_curve_ops, api_surface_ops};
use solidcore::surface::surface::{Surface as ScSurface, SurfaceKind as ScSurfaceKind};
use solidcore::sweep::api_surfaces;
use solidcore::{NlibCode, SolidError};

use crate::backend::{GeometryBackend, SurfaceConstruction, SurfaceIntersectionCandidate};
use crate::{BCurve, BSurface};

/// Upstream NVIDIA `usd-brep` revision matched by `solidcore` 1.0.1.
pub const UPSTREAM_REVISION: &str = "c9c5979f04cea8c38213004d4f04df0dcfb40d05";

/// Pinned `solidcore` crate version backing this module.
pub const SOLIDCORE_VERSION: &str = "1.0.1";

/// Alias of [`Backend`] for callers importing alongside [`crate::rust_surface::RustGeometryBackend`].
pub type SolidcoreBackend = Backend;

struct Inner {
    kernel: RefCell<Kernel>,
    precision: Precision,
}

/// An independent, reference-counted session of `solidcore` curve/surface geometry.
///
/// Unlike the legacy C++ bridge, separate `Backend` instances hold isolated [`Kernel`]
/// contexts without any process-global lock.
#[derive(Clone)]
pub struct Backend {
    inner: Rc<Inner>,
}

/// Owned handle to a curve in a `solidcore` session. Editing operations return new handles.
#[derive(Clone)]
pub struct Curve {
    inner: Rc<Inner>,
    curve: Rc<ScCurve>,
}

/// Owned handle to a surface in a `solidcore` session. Editing operations return new handles.
#[derive(Clone)]
pub struct Surface {
    inner: Rc<Inner>,
    surface: Rc<ScSurface>,
}

/// Curve position and first/second derivatives.
#[derive(Clone, Copy, Debug)]
pub struct CurveEvaluation {
    /// Position on the curve.
    pub point: Point3,
    /// First derivative with respect to parameter.
    pub d1: Vec3,
    /// Second derivative with respect to parameter.
    pub d2: Vec3,
}

/// Surface position and first/second partial derivatives.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceEvaluation {
    /// Position on the surface.
    pub point: Point3,
    /// First partial derivative in U.
    pub du: Vec3,
    /// First partial derivative in V.
    pub dv: Vec3,
    /// Second partial derivative in U.
    pub duu: Vec3,
    /// Mixed partial derivative in U and V.
    pub duv: Vec3,
    /// Second partial derivative in V.
    pub dvv: Vec3,
}

/// Parameterization method used by the interpolating curve fitter.
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

/// A numerical solver result, including parameter intervals for overlaps.
///
/// No global completeness certificate is implied by this type.
#[derive(Clone, Debug)]
pub struct Solution<const D: usize> {
    /// Upstream `SmSolutionType` discriminant (`1` = single value, `2` = range).
    pub native_kind: usize,
    /// Start parameters (`[t]` for point/curve, `[t1, t2]` for curve/curve, `[t, u, v]` for curve/surface).
    pub start: [f64; D],
    /// End parameters, meaningful when `native_kind == 2` (`RangeOfValues`).
    pub end: [f64; D],
}

/// One numerical surface-surface intersection branch with paired UV parameter curves.
pub struct IntersectionBranch {
    /// Three-dimensional world-space intersection curve.
    pub curve: Curve,
    /// Parameter-space curve on the first surface (`x = u`, `y = v`).
    pub first_uv: Curve,
    /// Parameter-space curve on the second surface (`x = u`, `y = v`).
    pub second_uv: Curve,
    /// Upstream `SmTsectCurveType` discriminant.
    pub native_kind: usize,
    /// Solver-reported deviation diagnostic; never treated as a certified bound.
    pub reported_deviation: f64,
}

/// Curves produced by projecting or dropping a 3D curve onto a surface.
pub struct ProjectedCurves {
    /// Curves in the surface UV parameter space.
    pub uv: Vec<Curve>,
    /// Lifted curves in 3D world coordinates.
    pub spatial: Vec<Curve>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct UnnamedObject;

impl SolutionObject for UnnamedObject {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Other
    }
}

#[derive(Clone, Copy)]
enum CurveSurfaceOperand<'a> {
    Curve(&'a ScCurve),
    Surface(&'a ScSurface),
}

impl PartialEq for CurveSurfaceOperand<'_> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Curve(a), Self::Curve(b)) => core::ptr::eq(*a, *b),
            (Self::Surface(a), Self::Surface(b)) => core::ptr::eq(*a, *b),
            _ => false,
        }
    }
}

impl SolutionObject for CurveSurfaceOperand<'_> {
    fn kind(&self) -> ObjectKind {
        match self {
            Self::Curve(_) => ObjectKind::Curve,
            Self::Surface(_) => ObjectKind::Surface,
        }
    }

    fn natural_interval(&self) -> Option<Extent1d> {
        match self {
            Self::Curve(c) => Some(c.natural_interval()),
            Self::Surface(_) => None,
        }
    }

    fn natural_uv_domain(&self) -> Option<Extent2d> {
        match self {
            Self::Curve(_) => None,
            Self::Surface(s) => Some(s.natural_uv_domain_if_built()),
        }
    }
}

fn map_error(err: SolidError) -> Error {
    let message = format!("solidcore error: {err}");
    let code = match err {
        SolidError::InvalidInput { .. } => ErrorCode::InvalidArgument,
        SolidError::NotFound { .. } => ErrorCode::NotFound,
        SolidError::Degenerate { .. } => ErrorCode::Degenerate,
        SolidError::Unsupported { .. } => ErrorCode::Unsupported,
        SolidError::Unported { .. } => ErrorCode::Unimplemented,
        SolidError::Nlib { code, .. } => match code {
            NlibCode::Memory => ErrorCode::CapacityExceeded,
            NlibCode::Input => ErrorCode::InvalidArgument,
            _ => ErrorCode::Internal,
        },
        SolidError::Upstream { code, .. } => {
            if code == api_curves::SM_ERR_NOT_WITHIN_TOLERANCE {
                ErrorCode::Unsupported
            } else {
                ErrorCode::Internal
            }
        }
        _ => ErrorCode::Internal,
    };
    Error::new(code, message)
}

fn sc_point(p: Point3) -> Result<ScPoint3> {
    if !p.is_finite() {
        return Err(Error::invalid_argument("point coordinates must be finite"));
    }
    Ok(ScPoint3::new(p.x, p.y, p.z))
}

fn sc_vec(v: Vec3) -> Result<ScVec3> {
    if !v.is_finite() {
        return Err(Error::invalid_argument("vector components must be finite"));
    }
    Ok(ScVec3::new(v.x, v.y, v.z))
}

fn gs_point(p: ScPoint3) -> Result<Point3> {
    if p.is_undef() || !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
        return Err(Error::internal("solidcore returned non-finite point"));
    }
    Ok(Point3::new(p.x, p.y, p.z))
}

fn gs_vec(v: ScVec3) -> Result<Vec3> {
    if v.is_undef() || !v.x.is_finite() || !v.y.is_finite() || !v.z.is_finite() {
        return Err(Error::internal("solidcore returned non-finite derivative"));
    }
    Ok(Vec3::new(v.x, v.y, v.z))
}

fn sc_points(points: &[Point3]) -> Result<Vec<ScPoint3>> {
    points.iter().copied().map(sc_point).collect()
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

fn require_finite(value: f64, what: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::invalid_argument(format!("{what} must be finite")))
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

fn validate_interval([lo, hi]: [f64; 2], [min, max]: [f64; 2]) -> Result<()> {
    if !lo.is_finite() || !hi.is_finite() || lo >= hi || lo < min || hi > max {
        Err(Error::invalid_argument("invalid trim domain"))
    } else {
        Ok(())
    }
}

fn knots_match(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.total_cmp(y).is_eq())
}

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

fn decode_cpoint(cp: CPoint) -> Result<(Point3, f64)> {
    let rational = !cp.w.total_cmp(&NOW).is_eq();
    let w = if rational { cp.w } else { 1.0 };
    if !w.is_finite() || w <= 0.0 {
        return Err(Error::invalid_argument(
            "NURBS weight must be positive and finite",
        ));
    }
    let z = if cp.z.total_cmp(&NOZ).is_eq() {
        0.0
    } else {
        cp.z / w
    };
    let x = cp.x / w;
    let y = cp.y / w;
    if !x.is_finite() || !y.is_finite() || !z.is_finite() {
        return Err(Error::internal(
            "solidcore returned non-finite control point",
        ));
    }
    Ok((Point3::new(x, y, z), w))
}

fn bspline_of_curve(k: &Kernel, curve: &ScCurve) -> Result<BSplineCurve> {
    let mut bs = match curve.kind() {
        ScCurveKind::BSpline(c) => c.clone(),
        ScCurveKind::Line(l) => {
            let mut copy = l.clone();
            copy.make_nurb().map_err(map_error)?;
            copy.base().clone()
        }
        ScCurveKind::Ellipse(e) => {
            let mut copy = e.clone();
            copy.make_nurb(k).map_err(map_error)?;
            copy.base().clone()
        }
        ScCurveKind::Hyperbola(h) => {
            let mut copy = h.clone();
            copy.make_nurb(None).map_err(map_error)?;
            copy.base().clone()
        }
        ScCurveKind::Parabola(p) => {
            let mut copy = p.clone();
            copy.make_nurb(None).map_err(map_error)?;
            copy.base().clone()
        }
        _ => {
            return Err(Error::new(
                ErrorCode::WrongEntityClass,
                "curve handle is not a B-spline-backed curve",
            ));
        }
    };
    bs.get_or_create_gw_nurb().map_err(map_error)?;
    Ok(bs)
}

fn bspline_of_surface(k: &Kernel, surface: &ScSurface) -> Result<BSplineSurface> {
    let mut bs = match surface.kind() {
        ScSurfaceKind::BSpline(s) => s.clone(),
        ScSurfaceKind::Plane(p) => p.base_with_nurb().map_err(map_error)?.into_owned(),
        ScSurfaceKind::SurfOfExtrusion(e) => {
            let mut copy = (**e).clone();
            copy.make_nurb(k).map_err(map_error)?;
            copy.base().clone()
        }
        ScSurfaceKind::Cone(c) => {
            let mut copy = (**c).clone();
            copy.make_nurb(k).map_err(map_error)?;
            copy.rev().base().clone()
        }
        ScSurfaceKind::Sphere(s) => {
            let mut copy = (**s).clone();
            copy.make_nurb(k).map_err(map_error)?;
            copy.base().base().clone()
        }
        ScSurfaceKind::Torus(t) => {
            let mut copy = (**t).clone();
            copy.make_nurb(k).map_err(map_error)?;
            copy.base().base().clone()
        }
        ScSurfaceKind::SurfOfRevolution(r) => {
            let mut copy = (**r).clone();
            copy.make_nurb(k).map_err(map_error)?;
            copy.base().clone()
        }
        _ => {
            return Err(Error::new(
                ErrorCode::WrongEntityClass,
                "surface handle is procedural and requires NURBS approximation",
            ));
        }
    };
    bs.get_or_create_gw_nurb().map_err(map_error)?;
    Ok(bs)
}

fn extract_solutions<const D: usize, O: SolutionObject>(
    arr: &SolutionArray<O>,
) -> Result<Vec<Solution<D>>> {
    let mut out = Vec::with_capacity(arr.as_slice().len());
    for sol in arr.as_slice() {
        let mut start = [0.0; D];
        let mut end = [0.0; D];
        let is_range = sol.solution_type == SolutionType::RangeOfValues;
        for i in 0..D {
            let s = sol.start.param(i).map_err(map_error)?;
            let e = if is_range {
                sol.end.param(i).map_err(map_error)?
            } else {
                s
            };
            if !s.is_finite() || !e.is_finite() {
                return Err(Error::internal(
                    "solidcore returned non-finite solution parameter",
                ));
            }
            start[i] = s;
            end[i] = e;
        }
        out.push(Solution {
            native_kind: sol.solution_type as usize,
            start,
            end,
        });
    }
    Ok(out)
}

impl Direction {
    fn to_sc(self) -> SurfParamType {
        match self {
            Self::U => SurfParamType::U,
            Self::V => SurfParamType::V,
        }
    }

    fn index(self) -> usize {
        match self {
            Self::U => 0,
            Self::V => 1,
        }
    }
}

impl Parameterization {
    fn to_sc(self) -> ParamMethod {
        match self {
            Self::Uniform => ParamMethod::Uniform,
            Self::ChordLength => ParamMethod::ChordLength,
            Self::Centripetal => ParamMethod::Centripetal,
        }
    }
}

impl Backend {
    /// Creates an isolated `solidcore` geometry session with tolerances derived from `precision`.
    pub fn new(precision: Precision) -> Result<Self> {
        let mut kernel = Kernel::new(&KernelOptions::default());
        solidcore::geom::tol::set_context_model_size_estimate(
            &mut kernel.context,
            precision.size(),
        );
        let zone = kernel.context.this_zone_tol_3d();
        kernel.tolerances.zone_tol_3d = zone;
        kernel.tolerances.xsect_tol_3d = 2.0 * zone;
        kernel.tolerances.approx_tol_3d = 0.5 * zone;
        kernel.tolerances.model_size_estimate = precision.size();
        Ok(Self {
            inner: Rc::new(Inner {
                kernel: RefCell::new(kernel),
                precision,
            }),
        })
    }

    /// Session precision passed to intersections, projections and approximations.
    pub fn precision(&self) -> Precision {
        self.inner.precision
    }

    fn same(&self, other: &Rc<Inner>) -> Result<()> {
        if Rc::ptr_eq(&self.inner, other) {
            Ok(())
        } else {
            Err(Error::invalid_argument(
                "geometry belongs to another solidcore session",
            ))
        }
    }

    fn wrap_curve(&self, curve: ScCurve) -> Curve {
        Curve {
            inner: self.inner.clone(),
            curve: Rc::new(curve),
        }
    }

    fn wrap_surface(&self, surface: ScSurface) -> Surface {
        Surface {
            inner: self.inner.clone(),
            surface: Rc::new(surface),
        }
    }

    /// Imports validated Gensolid NURBS curve data, preserving rational control points and knots.
    pub fn import_curve(&self, curve: &BCurve) -> Result<Curve> {
        let degree = curve.degree();
        let count = curve.points().len();
        let knots = curve.knots();
        if degree == 0
            || degree > i16::MAX as usize
            || count <= degree
            || knots.len() != count + degree + 1
        {
            return Err(Error::invalid_argument("invalid B-curve dimensions"));
        }
        if knots.iter().any(|k| !k.is_finite()) {
            return Err(Error::invalid_argument("curve knots must be finite"));
        }
        let mut pw = Vec::with_capacity(count);
        for (i, p) in curve.points().iter().enumerate() {
            let w = curve.weights().map_or(1.0, |weights| weights[i]);
            if !p.is_finite()
                || !w.is_finite()
                || w <= 0.0
                || !(p.x * w).is_finite()
                || !(p.y * w).is_finite()
                || !(p.z * w).is_finite()
            {
                return Err(Error::invalid_argument(
                    "curve control points and weights must be finite and positive",
                ));
            }
            pw.push(CPoint {
                x: p.x * w,
                y: p.y * w,
                z: p.z * w,
                w,
            });
        }
        let bs = BSplineCurve::from_data(3, count - 1, knots.len() - 1, degree, &pw, knots)
            .map_err(map_error)?;
        let stored = bs
            .gw_nurb()
            .ok_or_else(|| Error::internal("constructed B-spline curve lacks NURBS data"))?;
        if !knots_match(knots, stored.knots()) {
            return Err(Error::unsupported(
                "solidcore normalized input curve knots; representation was not preserved",
            ));
        }
        Ok(self.wrap_curve(ScCurve::from(bs)))
    }

    /// Imports a row-major (`u` then `v`) Gensolid NURBS surface, preserving knots and weights.
    pub fn import_surface(&self, surface: &BSurface) -> Result<Surface> {
        let (u, v) = surface.knots();
        let (pu, pv) = surface.degrees();
        let (nu, nv) = surface.counts();
        if pu == 0
            || pv == 0
            || pu > i16::MAX as usize
            || pv > i16::MAX as usize
            || nu <= pu
            || nv <= pv
            || u.len() != nu + pu + 1
            || v.len() != nv + pv + 1
        {
            return Err(Error::invalid_argument("invalid B-surface dimensions"));
        }
        if u.iter().any(|k| !k.is_finite()) || v.iter().any(|k| !k.is_finite()) {
            return Err(Error::invalid_argument("surface knots must be finite"));
        }
        let mut rows = Vec::with_capacity(nu);
        for i in 0..nu {
            let mut row = Vec::with_capacity(nv);
            for j in 0..nv {
                let idx = i * nv + j;
                let p = surface.points()[idx];
                let w = surface.weights().map_or(1.0, |weights| weights[idx]);
                if !p.is_finite()
                    || !w.is_finite()
                    || w <= 0.0
                    || !(p.x * w).is_finite()
                    || !(p.y * w).is_finite()
                    || !(p.z * w).is_finite()
                {
                    return Err(Error::invalid_argument(
                        "surface control points and weights must be finite and positive",
                    ));
                }
                row.push(CPoint {
                    x: p.x * w,
                    y: p.y * w,
                    z: p.z * w,
                    w,
                });
            }
            rows.push(row);
        }
        let bs = BSplineSurface::from_data(
            nu - 1,
            nv - 1,
            pu,
            pv,
            u.len() - 1,
            v.len() - 1,
            &rows,
            u,
            v,
        )
        .map_err(map_error)?;
        let stored = bs
            .gw_nurb()
            .ok_or_else(|| Error::internal("constructed B-spline surface lacks NURBS data"))?;
        if !knots_match(u, stored.u_knots()) || !knots_match(v, stored.v_knots()) {
            return Err(Error::unsupported(
                "solidcore normalized input surface knots; representation was not preserved",
            ));
        }
        Ok(self.wrap_surface(ScSurface::from(bs)))
    }

    /// Finite analytic line segment over `[0, 1]`.
    pub fn line_segment(&self, start: Point3, end: Point3) -> Result<Curve> {
        let a = sc_point(start)?;
        let b = sc_point(end)?;
        positive(start.distance(end), self.precision())?;
        let line = api_curves::create_line_segment(a, b).map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(line)))
    }

    /// Rational circle in the XY plane centered at `center`.
    pub fn circle(&self, center: Point3, radius: f64) -> Result<Curve> {
        let c = sc_point(center)?;
        positive(radius, self.precision())?;
        let circle = api_curves::create_circle(c, radius).map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(circle)))
    }

    /// Analytic ellipse in the XY plane centered at `center`.
    pub fn ellipse(&self, center: Point3, radii: [f64; 2]) -> Result<Curve> {
        let c = sc_point(center)?;
        for r in radii {
            positive(r, self.precision())?;
        }
        let k = self.inner.kernel.borrow();
        let ellipse = api_curves::create_ellipse(&k, c, radii[0], radii[1]).map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(ellipse)))
    }

    /// Rational circular arc in the XY plane; `angles` are in radians.
    pub fn arc(&self, center: Point3, radius: f64, angles: [f64; 2]) -> Result<Curve> {
        let c = sc_point(center)?;
        positive(radius, self.precision())?;
        let a0 = require_finite(angles[0], "arc start angle")?.to_degrees();
        let a1 = require_finite(angles[1], "arc end angle")?.to_degrees();
        let arc = api_curves::create_arc(c, radius, a0, a1).map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(arc)))
    }

    /// Helix along +Z, fitted using session linear precision.
    pub fn helix(
        &self,
        origin: Point3,
        height: f64,
        radii: [f64; 2],
        turns: f64,
        right_handed: bool,
    ) -> Result<Curve> {
        let o = sc_point(origin)?;
        for x in [height, radii[0], radii[1]] {
            positive(x, self.precision())?;
        }
        let num_turns = require_finite(turns, "helix turns")?;
        let curve = api_curves::create_helix(&HelixParams {
            origin: o,
            height,
            radius_start: radii[0],
            radius_end: radii[1],
            num_turns,
            right_handed,
            tolerance: self.precision().linear(),
        })
        .map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(curve)))
    }

    /// Cubic interpolation through at least four points.
    pub fn interpolate_curve(
        &self,
        input: &[Point3],
        parameterization: Parameterization,
    ) -> Result<Curve> {
        at_least(input.len(), 4)?;
        let pts = sc_points(input)?;
        let curve = api_curve_ops::create_curve_interp_points(&pts, parameterization.to_sc())
            .map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(curve)))
    }

    /// Curve approximation using `solidcore`'s fitting policy.
    pub fn approximate_curve(&self, input: &[Point3]) -> Result<Curve> {
        at_least(input.len(), 4)?;
        let pts = sc_points(input)?;
        let curve = api_curves::create_curve_approx_points(&pts).map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(curve)))
    }

    /// Cubic B-spline curve from at least four control points.
    pub fn control_curve(&self, input: &[Point3]) -> Result<Curve> {
        at_least(input.len(), 4)?;
        let pts = sc_points(input)?;
        let curve = api_curves::create_curve(&pts).map_err(map_error)?;
        Ok(self.wrap_curve(ScCurve::from(curve)))
    }

    /// Four sides of an XY rectangle centered at `center`.
    pub fn rectangle(&self, center: Point3, dimensions: [f64; 2]) -> Result<Vec<Curve>> {
        let c = sc_point(center)?;
        for x in dimensions {
            positive(x, self.precision())?;
        }
        let lines =
            api_curves::create_rectangle(c, dimensions[0], dimensions[1]).map_err(map_error)?;
        Ok(lines
            .into_iter()
            .map(|l| self.wrap_curve(ScCurve::from(l)))
            .collect())
    }

    /// Sides of a regular XY polygon centered at `center`.
    pub fn regular_polygon(&self, center: Point3, sides: usize, radius: f64) -> Result<Vec<Curve>> {
        let c = sc_point(center)?;
        at_least(sides, 3)?;
        if sides > i32::MAX as usize {
            return Err(Error::new(
                ErrorCode::CapacityExceeded,
                "polygon side count exceeds indexing limit",
            ));
        }
        positive(radius, self.precision())?;
        let lines =
            api_curves::create_regular_polygon(c, sides as u64, radius).map_err(map_error)?;
        Ok(lines
            .into_iter()
            .map(|l| self.wrap_curve(ScCurve::from(l)))
            .collect())
    }

    /// Bilinear patch through four corners ordered `(umin,vmin)`, `(umax,vmin)`, `(umin,vmax)`, `(umax,vmax)`.
    pub fn corner_surface(&self, corners: [Point3; 4]) -> Result<Surface> {
        let c0 = sc_point(corners[0])?;
        let c1 = sc_point(corners[1])?;
        let c2 = sc_point(corners[2])?;
        let c3 = sc_point(corners[3])?;
        let s = api_surface_ops::create_surface_from_corner_points(c0, c1, c2, c3)
            .map_err(map_error)?;
        Ok(self.wrap_surface(ScSurface::from(s)))
    }

    /// Bicubic surface from a control grid of at least four rows and columns.
    pub fn control_surface(
        &self,
        input: &[Point3],
        rows: usize,
        columns: usize,
    ) -> Result<Surface> {
        grid(input, rows, columns)?;
        let pts = sc_points(input)?;
        let s = api_surface_ops::create_surface(&pts, rows as u64, columns as u64)
            .map_err(map_error)?;
        Ok(self.wrap_surface(ScSurface::from(s)))
    }

    /// Surface fitted to an ordered `rows × columns` point grid.
    pub fn interpolate_surface(
        &self,
        input: &[Point3],
        rows: usize,
        columns: usize,
    ) -> Result<Surface> {
        grid(input, rows, columns)?;
        let pts = sc_points(input)?;
        let s =
            api_surface_ops::create_surface_from_ordered_points(&pts, rows as u64, columns as u64)
                .map_err(map_error)?;
        Ok(self.wrap_surface(ScSurface::from(s)))
    }

    /// Fits a surface to an unordered point set using a tight oriented bounding frame
    /// (`PointSet3d::tight_pbox`) and least-squares B-spline approximation.
    pub fn approximate_surface(&self, input: &[Point3]) -> Result<Surface> {
        at_least(input.len(), 4)?;
        let pts = sc_points(input)?;
        let mut pset = PointSet3d::new(pts.clone(), self.precision().linear(), None);
        let pbox = pset.tight_pbox(None).map_err(map_error)?;
        if !pbox.assert_defined() {
            return Err(Error::degenerate(
                "point cloud does not span a valid surface domain",
            ));
        }
        let b0 = pbox
            .basis_vector(0)
            .ok_or_else(|| Error::internal("missing pseudo-box u basis"))?;
        let b1 = pbox
            .basis_vector(1)
            .ok_or_else(|| Error::internal("missing pseudo-box v basis"))?;
        let b2 = pbox
            .basis_vector(2)
            .ok_or_else(|| Error::internal("missing pseudo-box w basis"))?;
        let i0 = pbox
            .interval(0)
            .ok_or_else(|| Error::internal("missing pseudo-box u interval"))?;
        let i1 = pbox
            .interval(1)
            .ok_or_else(|| Error::internal("missing pseudo-box v interval"))?;
        if i0.length() <= self.precision().degenerate_length()
            || i1.length() <= self.precision().degenerate_length()
        {
            return Err(Error::degenerate(
                "unordered point set is collinear or degenerate",
            ));
        }
        // Fit a linear height plane z(u, v) = c0 + c1 * u + c2 * v in the oriented box frame,
        // then build a bicubic grid over [0, 1]^2 through `lst_sq_approx`.
        let mut ata = [[0.0_f64; 3]; 3];
        let mut atb = [0.0_f64; 3];
        for &p in &pts {
            let local = pbox.invert_absolute(p).map_err(map_error)?;
            let u = (local.x - i0.min()) / i0.length();
            let v = (local.y - i1.min()) / i1.length();
            let row = [1.0, u, v];
            for r in 0..3 {
                atb[r] += row[r] * local.z;
                for c in 0..3 {
                    ata[r][c] += row[r] * row[c];
                }
            }
        }
        let coeffs = solve_symmetric_3x3(ata, atb).unwrap_or([0.0, 0.0, 0.0]);
        let sample_grid: Vec<Vec<ScVec3>> = (0..=4)
            .map(|iu| {
                let u = iu as f64 / 4.0;
                let x = i0.min() + u * i0.length();
                (0..=4)
                    .map(|iv| {
                        let v = iv as f64 / 4.0;
                        let y = i1.min() + v * i1.length();
                        let z = coeffs[0] + coeffs[1] * u + coeffs[2] * v;
                        b0 * x + b1 * y + b2 * z
                    })
                    .collect()
            })
            .collect();
        let fitted =
            lst_sq_approx(&sample_grid, 3, 3, 3, 3, ParamMethod::Uniform).map_err(map_error)?;
        let bs = BSplineSurface::from_nurb(fitted, false);
        Ok(self.wrap_surface(ScSurface::from(bs)))
    }

    /// Extrudes a B-spline-backed profile curve by `direction`.
    pub fn extrude(&self, profile: &Curve, direction: Vec3) -> Result<Surface> {
        self.same(&profile.inner)?;
        let dir = sc_vec(direction)?;
        positive(direction.length(), self.precision())?;
        let mut bs = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &profile.curve)?
        };
        let s = api_surfaces::create_extrude(&mut bs, dir).map_err(map_error)?;
        Ok(self.wrap_surface(ScSurface::from(s)))
    }

    /// Sweeps a B-spline-backed profile along a B-spline-backed path using `solidcore` framing.
    pub fn sweep(&self, profile: &Curve, path: &Curve) -> Result<Surface> {
        self.same(&profile.inner)?;
        self.same(&path.inner)?;
        let (prof, pth) = {
            let k = self.inner.kernel.borrow();
            (
                bspline_of_curve(&k, &profile.curve)?,
                bspline_of_curve(&k, &path.curve)?,
            )
        };
        let s = BSplineSurface::create_sweep_surface(
            &prof,
            &pth,
            None,
            false,
            self.precision().linear(),
        )
        .map_err(map_error)?;
        Ok(self.wrap_surface(ScSurface::from(s)))
    }

    /// Ruled surface between two B-spline-backed curves.
    pub fn ruled(&self, first: &Curve, second: &Curve) -> Result<Surface> {
        self.same(&first.inner)?;
        self.same(&second.inner)?;
        let (a, b) = {
            let k = self.inner.kernel.borrow();
            (
                bspline_of_curve(&k, &first.curve)?,
                bspline_of_curve(&k, &second.curve)?,
            )
        };
        let s = api_surfaces::create_ruled_surface(&a, &b).map_err(map_error)?;
        Ok(self.wrap_surface(ScSurface::from(s)))
    }

    /// Revolves a B-spline-backed profile about `(origin, axis)` by `angle` radians.
    pub fn revolve(
        &self,
        profile: &Curve,
        origin: Point3,
        axis: Vec3,
        angle: f64,
    ) -> Result<Surface> {
        self.same(&profile.inner)?;
        let o = sc_point(origin)?;
        let normalized = axis
            .normalized(self.precision().angular())
            .ok_or_else(|| Error::degenerate("revolution axis is zero"))?;
        let ax = sc_vec(normalized)?;
        let angle_deg = require_finite(angle, "revolution angle")?.to_degrees();
        let prof = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &profile.curve)?
        };
        let mut k = self.inner.kernel.borrow_mut();
        let s = api_surfaces::create_surface_revolution(&mut k, Some(&prof), o, ax, angle_deg)
            .map_err(map_error)?;
        Ok(self.wrap_surface(s))
    }

    /// Surface skin (loft) through at least two B-spline-backed section curves.
    pub fn loft(&self, sections: &[Curve]) -> Result<Surface> {
        at_least(sections.len(), 2)?;
        let k = self.inner.kernel.borrow();
        let mut profiles = Vec::with_capacity(sections.len());
        for c in sections {
            self.same(&c.inner)?;
            profiles.push(bspline_of_curve(&k, &c.curve)?);
        }
        let s = api_surfaces::create_skin(&k, &profiles).map_err(map_error)?;
        Ok(self.wrap_surface(ScSurface::from(s)))
    }

    /// Makes copies of `curves` with compatible degrees and knot vectors.
    pub fn compatible_curves(&self, curves: &[Curve]) -> Result<Vec<Curve>> {
        at_least(curves.len(), 2)?;
        let mut copies = Vec::with_capacity(curves.len());
        {
            let k = self.inner.kernel.borrow();
            for c in curves {
                self.same(&c.inner)?;
                copies.push(bspline_of_curve(&k, &c.curve)?);
            }
        }
        api_curve_ops::make_curves_compatible(&mut copies).map_err(map_error)?;
        Ok(copies
            .into_iter()
            .map(|c| self.wrap_curve(ScCurve::from(c)))
            .collect())
    }

    /// Orders copies of `curves` using `solidcore` loop ordering.
    pub fn order_curves(&self, curves: &[Curve]) -> Result<Vec<Curve>> {
        at_least(curves.len(), 1)?;
        let mut raw = Vec::with_capacity(curves.len());
        for c in curves {
            self.same(&c.inner)?;
            raw.push((*c.curve).clone());
        }
        let k = self.inner.kernel.borrow();
        let indices = api_curve_ops::order_curves(&k, &raw).map_err(map_error)?;
        indices
            .into_iter()
            .map(|idx| {
                raw.get(idx)
                    .cloned()
                    .map(|c| self.wrap_curve(c))
                    .ok_or_else(|| Error::internal("invalid ordered curve index"))
            })
            .collect()
    }
}

fn solve_symmetric_3x3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for pivot in 0..3 {
        let mut best = pivot;
        for row in (pivot + 1)..3 {
            if a[row][pivot].abs() > a[best][pivot].abs() {
                best = row;
            }
        }
        a.swap(pivot, best);
        b.swap(pivot, best);
        let diag = a[pivot][pivot];
        if !diag.is_finite() || diag.abs() <= f64::EPSILON {
            return None;
        }
        for row in (pivot + 1)..3 {
            let factor = a[row][pivot] / diag;
            b[row] -= factor * b[pivot];
            let pivot_row = a[pivot];
            for (entry, &pv) in a[row][pivot..3].iter_mut().zip(&pivot_row[pivot..3]) {
                *entry -= factor * pv;
            }
        }
    }
    let mut x = [0.0; 3];
    for row in (0..3).rev() {
        let mut sum = b[row];
        for col in (row + 1)..3 {
            sum -= a[row][col] * x[col];
        }
        let val = sum / a[row][row];
        if !val.is_finite() {
            return None;
        }
        x[row] = val;
    }
    Some(x)
}

impl Curve {
    fn backend(&self) -> Backend {
        Backend {
            inner: self.inner.clone(),
        }
    }

    /// Natural parameter interval `[lo, hi]`.
    pub fn domain(&self) -> Result<[f64; 2]> {
        let ivl = self.curve.natural_interval();
        let (lo, hi) = (ivl.min(), ivl.max());
        if !lo.is_finite() || !hi.is_finite() {
            return Err(Error::internal("curve has non-finite natural interval"));
        }
        Ok([lo, hi])
    }

    /// Position, first derivative and second derivative at `t`.
    pub fn evaluate(&self, t: f64) -> Result<CurveEvaluation> {
        let t = require_finite(t, "curve parameter")?;
        let ivl = self.curve.natural_interval();
        if !ivl.contains_value(t, 0.0) {
            return Err(Error::invalid_argument(
                "curve parameter lies outside natural interval",
            ));
        }
        let k = self.inner.kernel.borrow();
        let derivs = self
            .curve
            .evaluate(&k, t, 2, false, false)
            .map_err(map_error)?;
        if derivs.len() < 3 {
            return Err(Error::internal("insufficient derivative output"));
        }
        let normalize_2d = |mut v: ScVec3| {
            if v.z.total_cmp(&NOZ).is_eq() {
                v.z = 0.0;
            }
            v
        };
        Ok(CurveEvaluation {
            point: gs_point(normalize_2d(derivs[0]))?,
            d1: gs_vec(normalize_2d(derivs[1]))?,
            d2: gs_vec(normalize_2d(derivs[2]))?,
        })
    }

    /// Exports NURBS curve data into an owned [`BCurve`].
    pub fn nurbs(&self) -> Result<BCurve> {
        let bs = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &self.curve)?
        };
        let nurb = bs
            .gw_nurb()
            .ok_or_else(|| Error::internal("B-spline curve lacks NURBS storage"))?;
        let knots = nurb.knots().to_vec();
        let cpts = nurb.cpts();
        let mut points = Vec::with_capacity(cpts.len());
        let mut weights = Vec::with_capacity(cpts.len());
        for &cp in cpts {
            let (p, w) = decode_cpoint(cp)?;
            points.push(p);
            weights.push(w);
        }
        BCurve::new(nurb.degree(), knots, points, Some(weights))
    }

    /// Numerical closest-point parameter solutions on the curve.
    pub fn project_point(&self, point: Point3) -> Result<Vec<Solution<1>>> {
        let pt = sc_point(point)?;
        let k = self.inner.kernel.borrow();
        let ivl = self.curve.natural_interval();
        let opts = solidcore::curve::curve::p1::GlobalPointSolveOptions::<UnnamedObject>::default();
        let arr = self
            .curve
            .global_point_solve(
                &k,
                &ivl,
                SolverOperation::Minimize,
                pt,
                self.inner.precision.linear(),
                SolutionRequested::All,
                &opts,
            )
            .map_err(map_error)?;
        extract_solutions(&arr)
    }

    /// Numerical curve-curve intersections, preserving overlap parameter intervals.
    pub fn intersect_curve(&self, other: &Curve) -> Result<Vec<Solution<2>>> {
        self.backend().same(&other.inner)?;
        if Rc::ptr_eq(&self.curve, &other.curve) {
            let [lo, hi] = self.domain()?;
            return Ok(vec![Solution {
                native_kind: SolutionType::RangeOfValues as usize,
                start: [lo, lo],
                end: [hi, hi],
            }]);
        }
        let k = self.inner.kernel.borrow();
        let opts = solidcore::curve::curve::p8::CurveIntersectOptions::<UnnamedObject>::default();
        let arr = self
            .curve
            .global_curve_intersect(
                &k,
                &self.curve.natural_interval(),
                &other.curve,
                &other.curve.natural_interval(),
                self.inner.precision.linear(),
                &opts,
            )
            .map_err(map_error)?;
        extract_solutions(&arr)
    }

    /// Numerical curve-surface intersections; solution parameters are `[t, u, v]`.
    pub fn intersect_surface(&self, other: &Surface) -> Result<Vec<Solution<3>>> {
        self.backend().same(&other.inner)?;
        let k = self.inner.kernel.borrow();
        let domain = other
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let arr = other
            .surface
            .global_curve_intersect(
                &k,
                [
                    CurveSurfaceOperand::Curve(&self.curve),
                    CurveSurfaceOperand::Surface(&other.surface),
                ],
                &domain,
                &self.curve,
                &self.curve.natural_interval(),
                self.inner.precision.linear(),
            )
            .map_err(map_error)?;
        extract_solutions(&arr)
    }

    /// Offsets a planar curve in its detected plane.
    pub fn offset(&self, distance: f64) -> Result<Curve> {
        let distance = require_finite(distance, "offset distance")?;
        let k = self.inner.kernel.borrow();
        let bs = bspline_of_curve(&k, &self.curve)?;
        let planar = bs
            .is_planar(self.inner.precision.linear())
            .map_err(map_error)?;
        if !planar.is_planar {
            return Err(Error::unsupported("curve is not planar within precision"));
        }
        let (out, _) = bs
            .create_simple_offset(
                &k,
                self.inner.precision.linear(),
                planar.normal,
                distance,
                false,
            )
            .map_err(map_error)?;
        Ok(self.backend().wrap_curve(ScCurve::from(out)))
    }

    /// Offsets a curve in the explicit plane with `normal`.
    pub fn offset_in_plane(&self, distance: f64, normal: Vec3) -> Result<Curve> {
        let distance = require_finite(distance, "offset distance")?;
        let n = normal
            .normalized(self.inner.precision.angular())
            .ok_or_else(|| Error::degenerate("offset plane normal is zero"))?;
        let sc_n = sc_vec(n)?;
        let k = self.inner.kernel.borrow();
        let bs = bspline_of_curve(&k, &self.curve)?;
        let (out, _) = bs
            .create_simple_offset(&k, self.inner.precision.linear(), sc_n, distance, false)
            .map_err(map_error)?;
        Ok(self.backend().wrap_curve(ScCurve::from(out)))
    }

    /// Lifts a UV curve (`x = u`, `y = v`) through `surface` into 3D world space.
    pub fn lift(&self, surface: &Surface) -> Result<Vec<Curve>> {
        self.backend().same(&surface.inner)?;
        let k = self.inner.kernel.borrow();
        let bs = bspline_of_curve(&k, &self.curve)?;
        let uv_curve = ScCurve::from(bs);
        let domain = surface
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let lifted = surface
            .surface
            .lift_curve(
                &k,
                &domain,
                &uv_curve,
                &uv_curve.natural_interval(),
                self.inner.precision.linear(),
                false,
            )
            .map_err(map_error)?;
        Ok(vec![self.backend().wrap_curve(lifted.curve)])
    }

    /// Drops a curve onto `surface`, returning paired UV and spatial pieces.
    pub fn drop_to(&self, surface: &Surface) -> Result<ProjectedCurves> {
        self.drop_or_project(surface, false)
    }

    /// Projects and trims a curve onto `surface`, returning paired UV and spatial pieces.
    pub fn project_to(&self, surface: &Surface) -> Result<ProjectedCurves> {
        self.drop_or_project(surface, true)
    }

    fn drop_or_project(&self, surface: &Surface, trim_to_domain: bool) -> Result<ProjectedCurves> {
        self.backend().same(&surface.inner)?;
        let k = self.inner.kernel.borrow();
        let bs = bspline_of_curve(&k, &self.curve)?;
        let input = ScCurve::from(bs);
        let domain = surface
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let tol = self.inner.precision.linear();
        let uv_curves = if trim_to_domain {
            surface
                .surface
                .drop_and_trim_curve::<UnnamedObject>(
                    &k,
                    &domain,
                    &input,
                    &input.natural_interval(),
                    tol,
                    true,
                    false,
                    None,
                )
                .map_err(map_error)?
                .uv_curves
        } else {
            surface
                .surface
                .drop_curve::<UnnamedObject>(
                    &k,
                    &domain,
                    &input,
                    &input.natural_interval(),
                    tol,
                    true,
                    None,
                    false,
                )
                .map_err(map_error)?
                .uv_curves
        };
        let mut uv = Vec::with_capacity(uv_curves.len());
        let mut spatial = Vec::with_capacity(uv_curves.len());
        for c in uv_curves {
            let lifted = surface
                .surface
                .lift_curve(&k, &domain, &c, &c.natural_interval(), tol, false)
                .map_err(map_error)?;
            uv.push(self.backend().wrap_curve(c));
            spatial.push(self.backend().wrap_curve(lifted.curve));
        }
        Ok(ProjectedCurves { uv, spatial })
    }

    /// Removes redundant interior knots from a copy within session linear precision.
    pub fn simplify(&self) -> Result<Curve> {
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &self.curve)?
        };
        copy.remove_extra_knots(self.inner.precision.linear())
            .map_err(map_error)?;
        Ok(self.backend().wrap_curve(ScCurve::from(copy)))
    }

    /// Decomposes a curve into line segments where certified within linear precision,
    /// retaining curved spans as B-spline copies.
    pub fn lines_and_arcs(&self) -> Result<Vec<Curve>> {
        let copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &self.curve)?
        };
        let seg_line = copy
            .is_segment_line(&copy.natural_interval(), 16, self.inner.precision.linear())
            .map_err(map_error)?;
        if seg_line.is_line {
            let ends = copy.ends().map_err(map_error)?;
            if ends.start.distance(ends.end) > self.inner.precision.degenerate_length() {
                let line = BSplineCurve::create_line_segment(3, ends.start, ends.end)
                    .map_err(map_error)?;
                return Ok(vec![self.backend().wrap_curve(ScCurve::from(line))]);
            }
        }
        Ok(vec![self.backend().wrap_curve(ScCurve::from(copy))])
    }

    /// Upstream `SmContinuityType` discriminant at parameter `t`.
    pub fn continuity(&self, t: f64) -> Result<usize> {
        let t = require_finite(t, "curve parameter")?;
        let k = self.inner.kernel.borrow();
        let ct = api_curve_ops::evaluate_continuity(&k, &self.curve, t).map_err(map_error)?;
        Ok(ct as usize)
    }

    /// Trims a NURBS copy to `interval`.
    pub fn trim(&self, interval: [f64; 2]) -> Result<Curve> {
        validate_interval(interval, self.domain()?)?;
        let k = self.inner.kernel.borrow();
        let mut copy = bspline_of_curve(&k, &self.curve)?;
        let mut ext = Extent1d::from_raw(interval[0], interval[1]);
        copy.trim(&k, &mut ext, TrimMode::Silent)
            .map_err(map_error)?;
        Ok(self.backend().wrap_curve(ScCurve::from(copy)))
    }

    /// Reverses a NURBS copy while preserving its natural parameter interval.
    pub fn reversed(&self) -> Result<Curve> {
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &self.curve)?
        };
        let ivl = copy.natural_interval();
        copy.reverse_parameterization(&ivl).map_err(map_error)?;
        Ok(self.backend().wrap_curve(ScCurve::from(copy)))
    }

    /// Elevates the degree of a NURBS copy to `degree`.
    pub fn elevate_degree(&self, degree: usize) -> Result<Curve> {
        let original = self.nurbs()?.degree();
        if degree < original || degree > i16::MAX as usize {
            return Err(Error::invalid_argument("invalid elevated degree"));
        }
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &self.curve)?
        };
        copy.degree_elevate(degree).map_err(map_error)?;
        Ok(self.backend().wrap_curve(ScCurve::from(copy)))
    }

    /// Inserts interior `knots` into a NURBS copy.
    pub fn insert_knots(&self, knots: &[f64]) -> Result<Curve> {
        let spline = self.nurbs()?;
        let sorted = validated_insertions(spline.knots(), spline.degree(), knots)?;
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_curve(&k, &self.curve)?
        };
        copy.insert_knots(&sorted).map_err(map_error)?;
        Ok(self.backend().wrap_curve(ScCurve::from(copy)))
    }
}

impl Surface {
    fn backend(&self) -> Backend {
        Backend {
            inner: self.inner.clone(),
        }
    }

    /// Natural parameter domain as `[u_min, u_max, v_min, v_max]`.
    pub fn domain(&self) -> Result<[f64; 4]> {
        let k = self.inner.kernel.borrow();
        let dom = self
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let out = [dom.u_min(), dom.u_max(), dom.v_min(), dom.v_max()];
        if out.iter().any(|x| !x.is_finite()) {
            return Err(Error::internal("surface has non-finite natural domain"));
        }
        Ok(out)
    }

    /// Position and all first/second partial derivatives at `uv`.
    pub fn evaluate(&self, uv: [f64; 2]) -> Result<SurfaceEvaluation> {
        let u = require_finite(uv[0], "surface u parameter")?;
        let v = require_finite(uv[1], "surface v parameter")?;
        let k = self.inner.kernel.borrow();
        let dom = self
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let pt = ScPoint2::new(u, v);
        if !dom.contains_point(pt, 0.0) {
            return Err(Error::invalid_argument(
                "surface parameters lie outside natural domain",
            ));
        }
        let d = self
            .surface
            .evaluate_2nd_derivatives(&k, pt, false, false)
            .map_err(map_error)?;
        // `evaluate_2nd_derivatives` returns `[point, du, dv, duv, duu, dvv]`.
        Ok(SurfaceEvaluation {
            point: gs_point(d[0])?,
            du: gs_vec(d[1])?,
            dv: gs_vec(d[2])?,
            duu: gs_vec(d[4])?,
            duv: gs_vec(d[3])?,
            dvv: gs_vec(d[5])?,
        })
    }

    /// Exports NURBS surface data in Gensolid row-major (`u` then `v`) order.
    pub fn nurbs(&self) -> Result<BSurface> {
        let bs = {
            let k = self.inner.kernel.borrow();
            bspline_of_surface(&k, &self.surface)?
        };
        let nurb = bs
            .gw_nurb()
            .ok_or_else(|| Error::internal("B-spline surface lacks NURBS storage"))?;
        let nu = nurb.n() + 1;
        let nv = nurb.m() + 1;
        let u_knots = nurb.u_knots().to_vec();
        let v_knots = nurb.v_knots().to_vec();
        let mut points = Vec::with_capacity(nu * nv);
        let mut weights = Vec::with_capacity(nu * nv);
        for row in nurb.rows() {
            for &cp in row {
                let (p, w) = decode_cpoint(cp)?;
                points.push(p);
                weights.push(w);
            }
        }
        BSurface::new(
            nurb.p(),
            nurb.q(),
            u_knots,
            v_knots,
            nu,
            nv,
            points,
            Some(weights),
        )
    }

    /// Approximates a procedural surface as a NURBS surface within a bounded subdivision budget.
    pub fn approximate_nurbs(&self, subdivision_levels: usize) -> Result<(BSurface, f64)> {
        if !(1..=8).contains(&subdivision_levels) {
            return Err(Error::invalid_argument(
                "surface approximation budget must be 1..=8",
            ));
        }
        let approx = {
            let k = self.inner.kernel.borrow();
            self.surface
                .approximate_surface(
                    &k,
                    solidcore::geom::tol::approx_tol3d(&k, self.inner.precision.linear()),
                    subdivision_levels,
                    true,
                )
                .map_err(map_error)?
        };
        let approx_surface = approx
            .surface
            .ok_or_else(|| Error::unsupported("surface approximation exceeded deviation budget"))?;
        let wrapped = self.backend().wrap_surface(approx_surface);
        Ok((wrapped.nurbs()?, approx.achieved_tol))
    }

    /// Numerical closest-point UV solutions on the surface.
    pub fn project_point(&self, point: Point3) -> Result<Vec<Solution<2>>> {
        let pt = sc_point(point)?;
        let k = self.inner.kernel.borrow();
        let domain = self
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let arr = self
            .surface
            .global_point_solve::<UnnamedObject>(
                &k,
                None,
                &domain,
                SolverOperation::Minimize,
                pt,
                self.inner.precision.linear(),
                None,
                SolutionRequested::All,
            )
            .map_err(map_error)?;
        extract_solutions(&arr)
    }

    /// Numerical surface-surface intersection branches with paired UV curves.
    pub fn intersect(&self, other: &Surface) -> Result<Vec<IntersectionBranch>> {
        self.backend().same(&other.inner)?;
        if Rc::ptr_eq(&self.surface, &other.surface) {
            return Err(Error::unsupported(
                "identical surface coincidence is not a curve branch",
            ));
        }
        let k = self.inner.kernel.borrow();
        let dom1 = self
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let dom2 = other
            .surface
            .natural_uv_domain_with_nurb(&k)
            .map_err(map_error)?;
        let hit = self
            .surface
            .global_surface_intersect(
                &k,
                &dom1,
                &other.surface,
                &dom2,
                [true, true],
                Some(solidcore::geom::tol::approx_tol3d(
                    &k,
                    self.inner.precision.linear(),
                )),
                None,
            )
            .map_err(map_error)?;
        let count = hit.curves_3d.len();
        if hit.surface1_uv_curves.len() != count
            || hit.surface2_uv_curves.len() != count
            || hit.curve_types.len() != count
            || hit.deviations.len() != count
        {
            return Err(Error::internal("mismatched surface intersection arrays"));
        }
        let mut branches = Vec::with_capacity(count);
        for ((((xyz, uv1), uv2), kind), dev) in hit
            .curves_3d
            .into_iter()
            .zip(hit.surface1_uv_curves)
            .zip(hit.surface2_uv_curves)
            .zip(hit.curve_types)
            .zip(hit.deviations)
        {
            let first_uv =
                uv1.ok_or_else(|| Error::internal("missing first UV intersection curve"))?;
            let second_uv =
                uv2.ok_or_else(|| Error::internal("missing second UV intersection curve"))?;
            if !dev.is_finite() {
                return Err(Error::internal("non-finite intersection deviation"));
            }
            branches.push(IntersectionBranch {
                curve: self.backend().wrap_curve(xyz),
                first_uv: self.backend().wrap_curve(first_uv),
                second_uv: self.backend().wrap_curve(second_uv),
                native_kind: kind as usize,
                reported_deviation: dev,
            });
        }
        Ok(branches)
    }

    /// Normal offset surface by signed `distance`.
    pub fn offset(&self, distance: f64) -> Result<Surface> {
        let distance = require_finite(distance, "offset distance")?;
        let k = self.inner.kernel.borrow();
        let out = self
            .surface
            .create_offset_surface(&k, distance, self.inner.precision.linear())
            .map_err(map_error)?
            .ok_or_else(|| Error::degenerate("surface offset collapsed"))?;
        Ok(self.backend().wrap_surface(out))
    }

    /// Isoparametric curve with `direction` held constant at `parameter`.
    pub fn iso_curve(&self, direction: Direction, parameter: f64) -> Result<Curve> {
        let domain = self.domain()?;
        let offset = 2 * direction.index();
        if !parameter.is_finite() || parameter < domain[offset] || parameter > domain[offset + 1] {
            return Err(Error::invalid_argument(
                "iso parameter lies outside the domain",
            ));
        }
        let k = self.inner.kernel.borrow();
        let iso = self
            .surface
            .create_iso_parametric_curve(
                &k,
                direction.to_sc(),
                parameter,
                self.inner.precision.linear(),
                None,
                false,
            )
            .map_err(map_error)?;
        Ok(self.backend().wrap_curve(iso.curve))
    }

    /// Elevates the degree of a NURBS copy along `direction` to `degree`.
    pub fn elevate_degree(&self, direction: Direction, degree: usize) -> Result<Surface> {
        let (u, v) = self.nurbs()?.degrees();
        let original = match direction {
            Direction::U => u,
            Direction::V => v,
        };
        if degree < original || degree > i16::MAX as usize {
            return Err(Error::invalid_argument("invalid elevated degree"));
        }
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_surface(&k, &self.surface)?
        };
        copy.degree_elevate(direction.to_sc(), degree)
            .map_err(map_error)?;
        Ok(self.backend().wrap_surface(ScSurface::from(copy)))
    }

    /// Inserts interior `knots` along `direction` into a NURBS copy.
    pub fn insert_knots(&self, direction: Direction, knots: &[f64]) -> Result<Surface> {
        let spline = self.nurbs()?;
        let (u, v) = spline.knots();
        let (pu, pv) = spline.degrees();
        let (original, degree) = match direction {
            Direction::U => (u, pu),
            Direction::V => (v, pv),
        };
        let sorted = validated_insertions(original, degree, knots)?;
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_surface(&k, &self.surface)?
        };
        copy.insert_knots(direction.to_sc(), &sorted)
            .map_err(map_error)?;
        Ok(self.backend().wrap_surface(ScSurface::from(copy)))
    }

    /// Splits a NURBS surface copy at an interior `parameter` along `direction`.
    pub fn split(&self, direction: Direction, parameter: f64) -> Result<(Surface, Surface)> {
        let domain = self.domain()?;
        let offset = 2 * direction.index();
        if !parameter.is_finite() || parameter <= domain[offset] || parameter >= domain[offset + 1]
        {
            return Err(Error::invalid_argument("split must be inside the domain"));
        }
        let copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_surface(&k, &self.surface)?
        };
        let (left, right) = copy
            .split_at(parameter, direction.to_sc())
            .map_err(map_error)?;
        Ok((
            self.backend().wrap_surface(ScSurface::from(left)),
            self.backend().wrap_surface(ScSurface::from(right)),
        ))
    }

    /// Rectangular parameter-domain trim `[u_min, u_max, v_min, v_max]` of a NURBS copy.
    pub fn trim(&self, domain: [f64; 4]) -> Result<Surface> {
        let [u0, u1, v0, v1] = self.domain()?;
        validate_interval([domain[0], domain[1]], [u0, u1])?;
        validate_interval([domain[2], domain[3]], [v0, v1])?;
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_surface(&k, &self.surface)?
        };
        let mut ext = Extent2d::from_bounds(domain[0], domain[2], domain[1], domain[3]);
        copy.trim_with_domain(&mut ext).map_err(map_error)?;
        Ok(self.backend().wrap_surface(ScSurface::from(copy)))
    }

    /// Reverses one parameter direction of a NURBS copy.
    pub fn reversed(&self, direction: Direction) -> Result<Surface> {
        let mut copy = {
            let k = self.inner.kernel.borrow();
            bspline_of_surface(&k, &self.surface)?
        };
        copy.reverse(direction.to_sc()).map_err(map_error)?;
        Ok(self.backend().wrap_surface(ScSurface::from(copy)))
    }
}

impl GeometryBackend for Backend {
    fn precision(&self) -> Precision {
        self.precision()
    }

    fn construct_surface(&self, request: &SurfaceConstruction) -> Result<BSurface> {
        match request {
            SurfaceConstruction::Loft(sections) => {
                let curves = sections
                    .iter()
                    .map(|c| self.import_curve(c))
                    .collect::<Result<Vec<_>>>()?;
                let compatible = self.compatible_curves(&curves)?;
                self.loft(&compatible)?.nurbs()
            }
            SurfaceConstruction::Extrude(profile, direction) => self
                .extrude(&self.import_curve(profile)?, *direction)?
                .nurbs(),
            SurfaceConstruction::Sweep(profile, path) => self
                .sweep(&self.import_curve(profile)?, &self.import_curve(path)?)?
                .nurbs(),
            SurfaceConstruction::TranslationSweep { profile, path } => {
                crate::translation_sweep::sweep_translation(profile, path, &self.precision())
            }
            SurfaceConstruction::Revolve {
                profile,
                origin,
                axis,
                angle,
            } => self
                .revolve(&self.import_curve(profile)?, *origin, *axis, *angle)?
                .nurbs(),
        }
    }

    fn intersect_surfaces(
        &self,
        first: &BSurface,
        second: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>> {
        self.import_surface(first)?
            .intersect(&self.import_surface(second)?)?
            .into_iter()
            .map(|b| {
                Ok(SurfaceIntersectionCandidate {
                    curve: b.curve.nurbs()?,
                    first_uv: b.first_uv.nurbs()?,
                    second_uv: b.second_uv.nurbs()?,
                    reported_deviation: b.reported_deviation,
                })
            })
            .collect()
    }

    fn offset_surface(&self, surface: &BSurface, distance: f64) -> Result<BSurface> {
        if let Ok(cylindrical) =
            crate::cylindrical_offset::offset_cylindrical(surface, distance, &self.precision())
        {
            return Ok(cylindrical);
        }
        let offset = self.import_surface(surface)?.offset(distance)?;
        match offset.nurbs() {
            Ok(surface) => Ok(surface),
            Err(error) if error.code() == ErrorCode::WrongEntityClass => {
                Ok(offset.approximate_nurbs(4)?.0)
            }
            Err(error) => Err(error),
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use gensolid_math::Frame;
    use proptest::prelude::*;

    #[test]
    fn incompatible_loft_and_curved_sweep_succeed_without_ffi() {
        let p = Precision::DEFAULT;
        let backend = Backend::new(p).unwrap();
        let c1 = BCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            None,
        )
        .unwrap();
        let c2 = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(1.0, 2.0, 0.5),
                Point3::new(2.0, 2.0, 0.0),
            ],
            None,
        )
        .unwrap();
        let loft = backend
            .construct_surface(&SurfaceConstruction::Loft(vec![c1.clone(), c2]))
            .unwrap();
        assert!(loft.eval(0.5, 0.5).point.is_finite());

        let path = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 1.0),
                Point3::new(0.0, 2.0, 2.0),
            ],
            None,
        )
        .unwrap();
        let swept = backend
            .construct_surface(&SurfaceConstruction::Sweep(c1, path))
            .unwrap();
        assert!(swept.eval(0.5, 0.5).point.is_finite());
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        #[test]
        fn solidcore_surface_intersection_candidates_lie_on_both_carriers(
            width in 0.5f64..10.0,
            depth in 0.5f64..10.0,
            bend in -0.5f64..0.5,
            cut_ratio in 0.2f64..0.8,
            angle in -2.5f64..2.5,
            offset in -10.0f64..10.0,
        ) {
            let p = Precision::DEFAULT;
            let frame = Frame::from_axis(
                Point3::new(offset, -offset * 0.5, offset),
                Vec3::new(angle.cos(), angle.sin(), 1.0),
                &p,
            )?;
            let backend = Backend::new(p)?;
            let profile = BCurve::new(
                2,
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                vec![
                    frame.to_world(Vec3::new(0.0, 0.0, 0.0)),
                    frame.to_world(Vec3::new(width * 0.5, 0.0, bend)),
                    frame.to_world(Vec3::new(width, 0.0, 0.0)),
                ],
                None,
            )?;
            let extruded = backend.construct_surface(&SurfaceConstruction::Extrude(
                profile,
                frame.y() * depth,
            ))?;
            let y_cut = depth * cut_ratio;
            let plane = BSurface::new(
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                2,
                2,
                vec![
                    frame.to_world(Vec3::new(-1.0, y_cut, -2.0)),
                    frame.to_world(Vec3::new(-1.0, y_cut, 2.0)),
                    frame.to_world(Vec3::new(width + 1.0, y_cut, -2.0)),
                    frame.to_world(Vec3::new(width + 1.0, y_cut, 2.0)),
                ],
                None,
            )?;
            let candidates = backend.intersect_surfaces(&extruded, &plane)?;
            prop_assert!(!candidates.is_empty());
            let tol = p.linear() * 128.0 * (1.0 + width + depth);
            for branch in candidates {
                let r = branch.curve.range();
                for frac in [0.0, 0.25, 0.5, 0.75, 1.0] {
                    let t = r.lo + (r.hi - r.lo) * frac;
                    let xyz = branch.curve.eval(t).point;
                    let uv1 = branch.first_uv.eval(t).point;
                    let uv2 = branch.second_uv.eval(t).point;
                    prop_assert!(xyz.distance(extruded.eval(uv1.x, uv1.y).point) <= tol);
                    prop_assert!(xyz.distance(plane.eval(uv2.x, uv2.y).point) <= tol);
                }
            }
        }
    }
}
