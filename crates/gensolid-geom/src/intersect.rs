//! Analytic curve-surface and surface-surface intersections.
//!
//! Closed-form pairs return every branch. Pairs without a complete classification return
//! [`ErrorCode::Unimplemented`] instead of a partial curve. Sampling an intersection for a
//! Boolean is not done here.

use std::f64::consts::TAU;

use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_math::{Frame, Point3, Sign, Vec3};

use crate::curve::{Circle, Curve, Ellipse, Line};
use crate::surface::{Cone, Cylinder, Plane, Sphere, Surface, Torus};

mod certified;
mod cone_plane;
mod sphere_torus;
mod spline_plane;
mod surface_plane;
use certified::{Scalar, Vector};

/// Intersection of two curves, a curve and a surface, or two surfaces.
#[derive(Clone, Debug, PartialEq)]
pub enum IntersectionSet {
    /// No real intersection.
    Empty,
    /// Isolated points.
    Points(Vec<Point3>),
    /// Intersection curves. Lines are infinite; circles and ellipses are full turns.
    Curves(Vec<Curve>),
    /// The two surfaces are the same set of points.
    Overlap,
}

impl IntersectionSet {
    /// Number of branches. An overlap is one branch; the empty set is zero.
    pub fn branch_count(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Points(points) => points.len(),
            Self::Curves(curves) => curves.len(),
            Self::Overlap => 1,
        }
    }
}

/// Branch agreement and sampled point error against independently supplied answers.
///
/// This score establishes coverage only for its input cases. In particular, the analytic
/// catalogue does not establish completeness for arbitrary NURBS or surface pairs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IntersectionScore {
    /// Pairs that were attempted.
    pub pairs: usize,
    /// Pairs whose branch count matched the known answer and whose points lie on both
    /// surfaces.
    pub matched: usize,
    /// Pairs explicitly rejected as unimplemented; these remain in `pairs`.
    pub unimplemented: usize,
    /// Pairs that failed for another reason; these remain in `pairs`.
    pub failed: usize,
    /// Largest measured distance from an isolated point or curve sample to its inputs.
    pub max_point_error: f64,
}

impl IntersectionScore {
    /// Fraction of pairs that matched the known branch count and sampled accuracy.
    /// An empty workload provides no evidence and scores zero.
    pub fn completeness(self) -> f64 {
        if self.pairs == 0 {
            0.0
        } else {
            self.matched as f64 / self.pairs as f64
        }
    }
}

/// A reference case whose branch count comes from an independent known answer.
#[derive(Clone, Debug)]
pub enum IntersectionCase {
    /// Curve-surface reference case.
    CurveSurface {
        /// Input curve.
        curve: Curve,
        /// Input surface.
        surface: Surface,
        /// Expected number of isolated points, curves, or overlap components.
        expected_branches: usize,
    },
    /// Surface-surface reference case.
    Surfaces {
        /// First input surface.
        a: Surface,
        /// Second input surface.
        b: Surface,
        /// Expected number of isolated points, curves, or overlap components.
        expected_branches: usize,
    },
}

/// Score reference cases without discarding failed or unimplemented input pairs.
///
/// Successful curves are sampled for residuals against both original inputs. Branch
/// counts must come from an independent reference; this routine cannot discover a
/// missing branch on an arbitrary pair with no reference answer.
pub fn intersection_score(cases: &[IntersectionCase], precision: &Precision) -> IntersectionScore {
    let mut score = empty_score();
    for case in cases {
        score.pairs += 1;
        let pair = match case {
            IntersectionCase::CurveSurface {
                curve,
                surface,
                expected_branches,
            } => score_curve_surface(curve, surface, *expected_branches, precision),
            IntersectionCase::Surfaces {
                a,
                b,
                expected_branches,
            } => score_surfaces(a, b, *expected_branches, precision),
        };
        let measurement = pair.and_then(|pair| {
            let error = pair_error(&pair, precision)?;
            Ok((pair, error))
        });
        match measurement {
            Ok((pair, error)) => {
                score.max_point_error = score.max_point_error.max(error);
                if pair.found == pair.expected_branches && error <= precision.linear() {
                    score.matched += 1;
                }
            }
            Err(error) if error.code() == ErrorCode::Unimplemented => score.unimplemented += 1,
            Err(_) => score.failed += 1,
        }
    }
    score
}

fn empty_score() -> IntersectionScore {
    IntersectionScore {
        pairs: 0,
        matched: 0,
        unimplemented: 0,
        failed: 0,
        max_point_error: 0.0,
    }
}

/// Intersect a curve with a surface.
pub fn intersect_curve_surface(
    curve: &Curve,
    surface: &Surface,
    precision: &Precision,
) -> Result<IntersectionSet> {
    match (curve, surface) {
        (Curve::Line(line), Surface::Plane(plane)) => line_plane(line, plane, precision),
        (Curve::Line(line), Surface::Sphere(sphere)) => line_sphere(line, sphere, precision),
        (Curve::Line(line), Surface::Cylinder(cylinder)) => {
            line_cylinder(line, cylinder, precision)
        }
        (Curve::Line(line), Surface::Cone(cone)) => line_cone(line, cone, precision),
        (Curve::Circle(circle), Surface::Plane(plane)) => circle_plane(circle, plane, precision),
        (Curve::Circle(circle), Surface::Sphere(sphere)) => {
            circle_sphere(circle, sphere, precision)
        }
        (Curve::Circle(circle), Surface::Cylinder(cylinder)) => {
            circle_cylinder(circle, cylinder, precision)
        }
        (Curve::Ellipse(ellipse), Surface::Plane(plane)) => {
            ellipse_plane(ellipse, plane, precision)
        }
        (Curve::Ellipse(ellipse), Surface::Sphere(sphere)) => {
            ellipse_sphere(ellipse, sphere, precision)
        }
        (Curve::BCurve(curve), Surface::Plane(plane)) => {
            spline_plane::intersect(curve, plane, precision)
        }
        (Curve::BCurve(_), _) | (_, Surface::BSurface(_)) => Err(Error::unimplemented(
            "NURBS intersection is not implemented",
        )),
        (_, Surface::Swept(_) | Surface::Spun(_)) => Err(Error::unimplemented(
            "intersection with a procedural surface is not implemented",
        )),
        _ => Err(Error::unimplemented(
            "this curve-surface pair has no complete intersection",
        )),
    }
}

/// Intersect two surfaces.
pub fn intersect_surfaces(
    a: &Surface,
    b: &Surface,
    precision: &Precision,
) -> Result<IntersectionSet> {
    match (a, b) {
        (Surface::Plane(a), Surface::Plane(b)) => plane_plane(a, b, precision),
        (Surface::Plane(plane), Surface::Sphere(sphere))
        | (Surface::Sphere(sphere), Surface::Plane(plane)) => {
            plane_sphere(plane, sphere, precision)
        }
        (Surface::Plane(plane), Surface::Cylinder(cylinder))
        | (Surface::Cylinder(cylinder), Surface::Plane(plane)) => {
            plane_cylinder(plane, cylinder, precision)
        }
        (Surface::Plane(plane), Surface::Torus(torus))
        | (Surface::Torus(torus), Surface::Plane(plane)) => plane_torus(plane, torus, precision),
        (Surface::Plane(plane), Surface::Cone(cone))
        | (Surface::Cone(cone), Surface::Plane(plane)) => plane_cone(plane, cone, precision),
        (Surface::Sphere(a), Surface::Sphere(b)) => sphere_sphere(a, b, precision),
        (Surface::Cylinder(a), Surface::Cylinder(b)) => cylinder_cylinder(a, b, precision),
        (Surface::Sphere(sphere), Surface::Cylinder(cylinder))
        | (Surface::Cylinder(cylinder), Surface::Sphere(sphere)) => {
            sphere_cylinder(sphere, cylinder, precision)
        }
        (Surface::Sphere(sphere), Surface::Torus(torus))
        | (Surface::Torus(torus), Surface::Sphere(sphere)) => {
            sphere_torus(sphere, torus, precision)
        }
        (Surface::BSurface(patch), Surface::Plane(plane))
        | (Surface::Plane(plane), Surface::BSurface(patch)) => {
            surface_plane::intersect(patch, plane, precision)
        }
        (Surface::BSurface(_), _) | (_, Surface::BSurface(_)) => Err(Error::unimplemented(
            "NURBS intersection is not implemented",
        )),
        (Surface::Swept(_) | Surface::Spun(_), _) | (_, Surface::Swept(_) | Surface::Spun(_)) => {
            Err(Error::unimplemented(
                "intersection with a procedural surface is not implemented",
            ))
        }
        _ => Err(Error::new(
            ErrorCode::Unimplemented,
            "this surface pair has no complete intersection",
        )),
    }
}

/// Score known analytic pairs placed in `frames`.
///
/// Each frame is a rigid motion of the same relative geometry, so the expected branches
/// stay known. `max_point_error` is measured with [`crate::GeomStore::distance_to_surface`].
pub fn analytic_intersection_score(
    frames: &[Frame],
    precision: &Precision,
) -> Result<IntersectionScore> {
    let mut score = empty_score();
    for frame in frames {
        for pair in pairs_in_frame(frame, precision)? {
            let error = pair_error(&pair, precision)?;
            score.max_point_error = score.max_point_error.max(error);
            score.pairs += 1;
            let accurate = error <= precision.linear();
            if pair.found == pair.expected_branches && accurate {
                score.matched += 1;
            }
        }
    }
    Ok(score)
}

struct ScoredPair {
    found: usize,
    expected_branches: usize,
    samples: Vec<Point3>,
    surfaces: Vec<Surface>,
    curves: Vec<Curve>,
}

#[allow(clippy::vec_init_then_push)]
fn pairs_in_frame(frame: &Frame, precision: &Precision) -> Result<Vec<ScoredPair>> {
    let origin = frame.origin();
    let (x, z) = (frame.x(), frame.z());
    let sphere = Surface::Sphere(Sphere::new(*frame, 2.0, precision)?);
    let plane_z = Surface::Plane(Plane::new(*frame));
    let cylinder = Surface::Cylinder(Cylinder::new(*frame, 1.0, precision)?);
    let torus = Surface::Torus(Torus::new(*frame, 3.0, 1.0, precision)?);
    let other = Frame::new(origin + z * 1.0, z, x, precision)?;
    let sphere_b = Surface::Sphere(Sphere::new(other, 1.0, precision)?);
    let far = Frame::new(origin + z * 5.0, z, x, precision)?;
    let sphere_far = Surface::Sphere(Sphere::new(far, 1.0, precision)?);
    let plane_x = Surface::Plane(Plane::new(Frame::new(origin, x, z, precision)?));
    let line = Curve::Line(Line::new(origin, x, precision)?);
    let circle = Curve::Circle(Circle::new(*frame, 2.0, precision)?);

    let mut pairs = Vec::new();
    pairs.push(score_surfaces(&plane_z, &sphere, 1, precision)?);
    pairs.push(score_surfaces(
        &Surface::Sphere(Sphere::new(*frame, 1.0, precision)?),
        &sphere_b,
        1,
        precision,
    )?);
    pairs.push(score_surfaces(&plane_z, &cylinder, 1, precision)?);
    pairs.push(score_surfaces(&plane_x, &cylinder, 2, precision)?);
    pairs.push(score_surfaces(&plane_z, &plane_x, 1, precision)?);
    let tilt = Surface::Plane(Plane::new(Frame::new(origin, x + z, frame.y(), precision)?));
    let cone = Surface::Cone(Cone::new(*frame, 1.0, 0.5, precision)?);
    pairs.push(score_surfaces(&tilt, &cylinder, 1, precision)?);
    pairs.push(score_surfaces(&plane_z, &cone, 1, precision)?);
    pairs.push(score_surfaces(&plane_z, &torus, 2, precision)?);
    pairs.push(score_surfaces(
        &Surface::Sphere(Sphere::new(*frame, 1.0, precision)?),
        &sphere_far,
        0,
        precision,
    )?);
    pairs.push(score_curve_surface(&line, &sphere, 2, precision)?);
    pairs.push(score_curve_surface(&line, &plane_z, 1, precision)?);
    pairs.push(score_curve_surface(&circle, &plane_x, 2, precision)?);
    Ok(pairs)
}

fn score_surfaces(
    a: &Surface,
    b: &Surface,
    expected_branches: usize,
    precision: &Precision,
) -> Result<ScoredPair> {
    let hit = intersect_surfaces(a, b, precision)?;
    let witnesses = if matches!(hit, IntersectionSet::Overlap) {
        let mut witnesses = surface_witnesses(a, precision)?;
        witnesses.extend(surface_witnesses(b, precision)?);
        witnesses
    } else {
        samples(&hit, precision)?
    };
    Ok(ScoredPair {
        found: hit.branch_count(),
        expected_branches,
        samples: witnesses,
        surfaces: vec![a.clone(), b.clone()],
        curves: Vec::new(),
    })
}

fn score_curve_surface(
    curve: &Curve,
    surface: &Surface,
    expected_branches: usize,
    precision: &Precision,
) -> Result<ScoredPair> {
    let hit = intersect_curve_surface(curve, surface, precision)?;
    Ok(ScoredPair {
        found: hit.branch_count(),
        expected_branches,
        samples: curve_hit_samples(&hit, curve, precision)?,
        surfaces: vec![surface.clone()],
        curves: vec![curve.clone()],
    })
}

fn curve_hit_samples(
    hit: &IntersectionSet,
    curve: &Curve,
    precision: &Precision,
) -> Result<Vec<Point3>> {
    if matches!(hit, IntersectionSet::Overlap) {
        samples(&IntersectionSet::Curves(vec![curve.clone()]), precision)
    } else {
        samples(hit, precision)
    }
}

/// Overlap has no explicit output geometry, so its accuracy needs finite witnesses
/// from the original surfaces. Both inputs are sampled by the scorer. Procedural
/// profiles cannot be recovered from a standalone Surface value without its store.
fn surface_witnesses(surface: &Surface, precision: &Precision) -> Result<Vec<Point3>> {
    if matches!(surface, Surface::Swept(_) | Surface::Spun(_)) {
        return Err(Error::unimplemented(
            "overlap surface witnesses require the procedural profile store",
        ));
    }
    let parameters = |range: crate::ParamRange| -> Result<Vec<f64>> {
        if range.is_bounded() {
            let width = range.hi - range.lo;
            if !(width.is_finite() && width > 0.0) {
                return Err(Error::unimplemented(
                    "overlap witness domain cannot be evaluated",
                ));
            }
            Ok((0..=4)
                .map(|i| range.lo + width * (i as f64 / 4.0))
                .collect())
        } else {
            let extent = precision.size() * 0.25;
            Ok(vec![-extent, 0.0, extent])
        }
    };
    let mut store = crate::GeomStore::new();
    let id = store.add_surface(surface.clone()).map_err(|error| {
        Error::unimplemented(format!(
            "overlap surface witnesses cannot be stored: {error}"
        ))
    })?;
    let (ur, vr) = store.surface_ranges(id)?;
    let us = parameters(ur)?;
    let vs = parameters(vr)?;
    let mut witnesses = Vec::new();
    for u in us {
        for &v in &vs {
            let point = store
                .eval_surface(id, u, v)
                .map_err(|error| {
                    Error::unimplemented(format!(
                        "overlap surface witness cannot be evaluated: {error}"
                    ))
                })?
                .point;
            witnesses.push(point);
        }
    }
    if witnesses.is_empty() {
        return Err(Error::unimplemented(
            "overlap surface has no finite witnesses",
        ));
    }
    Ok(witnesses)
}

fn samples(hit: &IntersectionSet, precision: &Precision) -> Result<Vec<Point3>> {
    let mut points = Vec::new();
    match hit {
        IntersectionSet::Empty => {}
        IntersectionSet::Overlap => {
            return Err(Error::unimplemented(
                "overlap sampling requires the original geometry",
            ));
        }
        IntersectionSet::Points(found) => points.extend(found.iter().copied()),
        IntersectionSet::Curves(curves) => {
            for curve in curves {
                let ts: Vec<f64> = match curve {
                    Curve::Line(_) => vec![-1.0, 0.0, 1.0],
                    Curve::Circle(_) | Curve::Ellipse(_) => {
                        (0..8).map(|i| TAU * i as f64 / 8.0).collect()
                    }
                    Curve::BCurve(curve) => {
                        let range = curve.range();
                        let mut parameters = Vec::new();
                        for pair in curve.knots().windows(2) {
                            let lo = pair[0].max(range.lo);
                            let hi = pair[1].min(range.hi);
                            if lo < hi {
                                parameters.push(lo);
                                parameters.push(lo + (hi - lo) * 0.5);
                                parameters.push(if hi < range.hi { hi.next_down() } else { hi });
                            }
                        }
                        parameters
                    }
                };
                for t in ts {
                    let point = curve.eval(t).point;
                    if !point.is_finite() {
                        return Err(Error::internal("intersection sample is not finite"));
                    }
                    let _ = precision;
                    points.push(point);
                }
            }
        }
    }
    Ok(points)
}

fn pair_error(pair: &ScoredPair, precision: &Precision) -> Result<f64> {
    let mut store = crate::GeomStore::new();
    let mut ids = Vec::new();
    for surface in &pair.surfaces {
        ids.push(store.add_surface(surface.clone())?);
    }
    let curve_ids = pair
        .curves
        .iter()
        .map(|curve| store.add_curve(curve.clone()))
        .collect::<Result<Vec<_>>>()?;
    let mut max_error: f64 = 0.0;
    for point in &pair.samples {
        for id in &ids {
            let distance = store.distance_to_surface(*id, *point, precision)?;
            max_error = max_error.max(distance);
        }
        for id in &curve_ids {
            let distance = store.project_to_curve(*id, *point, precision)?.distance;
            max_error = max_error.max(distance);
        }
    }
    Ok(max_error)
}

fn parallel(a: Vec3, b: Vec3, precision: &Precision) -> Result<bool> {
    Vector::new(a)
        .cross(Vector::new(b))
        .zero(precision.angular() * a.length() * b.length())
}

fn line_from(origin: Point3, direction: Vec3, precision: &Precision) -> Result<Curve> {
    Ok(Curve::Line(Line::new(origin, unit(direction)?, precision)?))
}

fn unit(direction: Vec3) -> Result<Vec3> {
    direction
        .normalized(0.0)
        .ok_or_else(|| Error::unsupported("intersection direction cannot be normalized"))
}

fn squared_length_tolerance(length: f64, precision: &Precision) -> f64 {
    precision.linear() * (2.0 * length.abs() + precision.linear())
}

fn circle_in(
    center: Point3,
    normal: Vec3,
    reference: Vec3,
    radius: f64,
    precision: &Precision,
) -> Result<Curve> {
    if radius <= precision.linear() {
        return Err(Error::degenerate("circle radius collapsed"));
    }
    let frame = Frame::new(center, normal, reference, precision)?;
    Ok(Curve::Circle(Circle::new(frame, radius, precision)?))
}

fn solve_quadratic(
    a: Scalar,
    b: Scalar,
    c: Scalar,
    precision: &Precision,
) -> Result<Option<Vec<f64>>> {
    if a.sign(precision.angular())? == Sign::Zero {
        if b.sign(precision.linear())? == Sign::Zero {
            return Ok(
                if c.sign(squared_length_tolerance(c.value.abs().sqrt(), precision))? == Sign::Zero
                {
                    None
                } else {
                    Some(Vec::new())
                },
            );
        }
        return Ok(Some(vec![-c.value / b.value]));
    }
    let disc = b.sqr() - Scalar::new(4.0) * a * c;
    let tolerance =
        squared_length_tolerance(b.value.abs() + (a.value * c.value).abs().sqrt(), precision);
    match disc.sign(tolerance)? {
        Sign::Negative => Ok(Some(Vec::new())),
        Sign::Zero => Ok(Some(vec![-b.value / (2.0 * a.value)])),
        Sign::Positive => {
            let root = disc.value.sqrt();
            let q = -0.5 * (b.value + root.copysign(b.value));
            Ok(Some(vec![q / a.value, c.value / q]))
        }
    }
}

fn line_plane(line: &Line, plane: &Plane, precision: &Precision) -> Result<IntersectionSet> {
    let normal = plane.frame().z();
    let denom = Vector::new(line.direction()).dot(Vector::new(normal));
    let signed = Vector::between(line.origin(), plane.frame().origin()).dot(Vector::new(normal));
    if denom.sign(precision.angular())? == Sign::Zero {
        return Ok(if signed.sign(precision.linear())? == Sign::Zero {
            IntersectionSet::Curves(vec![Curve::Line(*line)])
        } else {
            IntersectionSet::Empty
        });
    }
    let t = -signed.value / denom.value;
    Ok(IntersectionSet::Points(vec![line.eval(t).point]))
}

fn line_sphere(line: &Line, sphere: &Sphere, precision: &Precision) -> Result<IntersectionSet> {
    let offset = Vector::between(line.origin(), sphere.frame().origin());
    let direction = Vector::new(line.direction());
    circular_line_section(line, direction, offset, sphere.radius(), precision)
}

fn circular_line_section(
    line: &Line,
    direction: Vector,
    offset: Vector,
    radius: f64,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let a = direction.norm2();
    if a.sign(precision.angular() * precision.angular())? == Sign::Zero {
        let gap = offset.norm2() - Scalar::new(radius).sqr();
        return Ok(
            if gap.sign(squared_length_tolerance(radius, precision))? == Sign::Zero {
                IntersectionSet::Curves(vec![Curve::Line(*line)])
            } else {
                IntersectionSet::Empty
            },
        );
    }
    let center_t = -offset.dot(direction) / a;
    let nearest = offset - direction * -center_t;
    let radial2 = Scalar::new(radius).sqr() - nearest.norm2();
    match radial2.sign(squared_length_tolerance(radius, precision))? {
        Sign::Negative => Ok(IntersectionSet::Empty),
        Sign::Zero => Ok(points_from_line(line, &[center_t.value])),
        Sign::Positive => {
            let delta = (radial2.value / a.value).sqrt();
            Ok(points_from_line(
                line,
                &[center_t.value - delta, center_t.value + delta],
            ))
        }
    }
}

fn line_cylinder(
    line: &Line,
    cylinder: &Cylinder,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let axis = Vector::new(cylinder.frame().z());
    let offset = Vector::between(line.origin(), cylinder.frame().origin());
    circular_line_section(
        line,
        Vector::new(line.direction()).cross(axis),
        offset.cross(axis),
        cylinder.radius(),
        precision,
    )
}

fn line_cone(line: &Line, cone: &Cone, precision: &Precision) -> Result<IntersectionSet> {
    let axis = Vector::new(cone.frame().z());
    let offset = Vector::between(line.origin(), cone.frame().origin());
    let direction = Vector::new(line.direction());
    let (radius, angle) = cone.shape();
    let slope = Scalar::new(angle.tan());
    let h = Scalar::new(radius) + offset.dot(axis) * slope;
    let dh = direction.dot(axis) * slope;
    let radial = offset.cross(axis);
    let radial_direction = direction.cross(axis);
    let a = radial_direction.norm2() - dh.sqr();
    let b = Scalar::new(2.0) * (radial.dot(radial_direction) - h * dh);
    let c = radial.norm2() - h.sqr();
    match solve_quadratic(a, b, c, precision)? {
        Some(roots) => Ok(points_from_line(line, &roots)),
        None => Ok(IntersectionSet::Curves(vec![Curve::Line(*line)])),
    }
}

fn points_from_line(line: &Line, roots: &[f64]) -> IntersectionSet {
    if roots.is_empty() {
        IntersectionSet::Empty
    } else {
        IntersectionSet::Points(roots.iter().map(|t| line.eval(*t).point).collect())
    }
}

fn trig_roots(a: Scalar, b: Scalar, c: Scalar, precision: &Precision) -> Result<Option<Vec<f64>>> {
    let scale2 = a.sqr() + b.sqr();
    if a.sign(precision.linear())? == Sign::Zero && b.sign(precision.linear())? == Sign::Zero {
        return Ok(if c.sign(precision.linear())? == Sign::Zero {
            None
        } else {
            Some(Vec::new())
        });
    }
    let sign = (scale2 - c.sqr()).sign(squared_length_tolerance(scale2.value.sqrt(), precision))?;
    if sign == Sign::Negative {
        return Ok(Some(Vec::new()));
    }
    let ratio = c.value / scale2.value.sqrt();
    let clamped = ratio.clamp(-1.0, 1.0);
    let alpha = b.value.atan2(a.value);
    let delta = clamped.acos();
    Ok(Some(if sign == Sign::Zero {
        vec![(alpha + delta).rem_euclid(TAU)]
    } else {
        vec![
            (alpha + delta).rem_euclid(TAU),
            (alpha - delta).rem_euclid(TAU),
        ]
    }))
}

fn circle_plane(circle: &Circle, plane: &Plane, precision: &Precision) -> Result<IntersectionSet> {
    let center = circle.frame().origin();
    let normal = plane.frame().z();
    let normal = Vector::new(normal);
    let a = Scalar::new(circle.radius()) * Vector::new(circle.frame().x()).dot(normal);
    let b = Scalar::new(circle.radius()) * Vector::new(circle.frame().y()).dot(normal);
    let c = Vector::between(plane.frame().origin(), center).dot(normal);
    match trig_roots(a, b, c, precision)? {
        None => Ok(IntersectionSet::Curves(vec![Curve::Circle(*circle)])),
        Some(roots) if roots.is_empty() => Ok(IntersectionSet::Empty),
        Some(roots) => Ok(IntersectionSet::Points(
            roots.into_iter().map(|t| circle.eval(t).point).collect(),
        )),
    }
}

fn ellipse_plane(
    ellipse: &Ellipse,
    plane: &Plane,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let (major, minor) = ellipse.radii();
    let center = ellipse.frame().origin();
    let normal = plane.frame().z();
    let normal = Vector::new(normal);
    let a = Scalar::new(major) * Vector::new(ellipse.frame().x()).dot(normal);
    let b = Scalar::new(minor) * Vector::new(ellipse.frame().y()).dot(normal);
    let c = Vector::between(plane.frame().origin(), center).dot(normal);
    match trig_roots(a, b, c, precision)? {
        None => Ok(IntersectionSet::Curves(vec![Curve::Ellipse(*ellipse)])),
        Some(roots) if roots.is_empty() => Ok(IntersectionSet::Empty),
        Some(roots) => Ok(IntersectionSet::Points(
            roots.into_iter().map(|t| ellipse.eval(t).point).collect(),
        )),
    }
}

fn circle_sphere(
    circle: &Circle,
    sphere: &Sphere,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let plane = Plane::new(*circle.frame());
    let section = plane_sphere(&plane, sphere, precision)?;
    match section {
        IntersectionSet::Empty => Ok(IntersectionSet::Empty),
        IntersectionSet::Points(points) => {
            let mut on_circle = Vec::new();
            for point in points {
                let gap = Vector::between(point, circle.frame().origin()).norm2()
                    - Scalar::new(circle.radius()).sqr();
                if gap.sign(squared_length_tolerance(circle.radius(), precision))? == Sign::Zero {
                    on_circle.push(point);
                }
            }
            Ok(if on_circle.is_empty() {
                IntersectionSet::Empty
            } else {
                IntersectionSet::Points(on_circle)
            })
        }
        IntersectionSet::Curves(curves) => {
            let Some(Curve::Circle(other)) = curves.first() else {
                return Err(Error::internal("plane-sphere did not return a circle"));
            };
            circle_circle(circle, other, precision)
        }
        IntersectionSet::Overlap => Ok(IntersectionSet::Curves(vec![Curve::Circle(*circle)])),
    }
}

/// Every circle perpendicular to the cylinder axis reduces to a coplanar
/// circle-circle intersection, including tangencies and complete overlap.
fn circle_cylinder(
    circle: &Circle,
    cylinder: &Cylinder,
    precision: &Precision,
) -> Result<IntersectionSet> {
    if !parallel(circle.frame().z(), cylinder.frame().z(), precision)? {
        return Err(Error::unimplemented(
            "circle-cylinder intersection requires a circle perpendicular to the cylinder axis",
        ));
    }
    let section = plane_cylinder(&Plane::new(*circle.frame()), cylinder, precision)?;
    let IntersectionSet::Curves(curves) = section else {
        return Err(Error::internal(
            "perpendicular cylinder section is not a circle",
        ));
    };
    let Some(Curve::Circle(other)) = curves.first() else {
        return Err(Error::internal(
            "perpendicular cylinder section is not a circle",
        ));
    };
    circle_circle(circle, other, precision)
}

/// A sphere whose center lies on the ellipse axis cuts a concentric circle
/// from its plane. Comparing its radius to both ellipse axes gives all zero,
/// two, or four isolated points, or the whole ellipse when it is a circle.
fn ellipse_sphere(
    ellipse: &Ellipse,
    sphere: &Sphere,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let normal = Vector::new(ellipse.frame().z());
    let delta = Vector::between(sphere.frame().origin(), ellipse.frame().origin());
    if !delta.cross(normal).zero(precision.linear())? {
        return Err(Error::unimplemented(
            "ellipse-sphere intersection requires a sphere center on the ellipse axis",
        ));
    }
    let (major, minor) = ellipse.radii();
    let radius2 = Scalar::new(sphere.radius()).sqr() - delta.dot(normal).sqr() / normal.norm2();
    let major2 = Scalar::new(major).sqr();
    let minor2 = Scalar::new(minor).sqr();
    let low = radius2 - minor2;
    let high = major2 - radius2;
    let tolerance = squared_length_tolerance(sphere.radius().max(major), precision);
    let low_sign = low.sign(tolerance)?;
    let high_sign = high.sign(tolerance)?;
    if low_sign == Sign::Negative || high_sign == Sign::Negative {
        return Ok(IntersectionSet::Empty);
    }
    if (major2 - minor2).sign(squared_length_tolerance(major, precision))? == Sign::Zero {
        return Ok(IntersectionSet::Curves(vec![Curve::Ellipse(*ellipse)]));
    }
    let (x, y) = (ellipse.frame().x(), ellipse.frame().y());
    let center = ellipse.frame().origin();
    if low_sign == Sign::Zero {
        return Ok(IntersectionSet::Points(vec![
            center + y * minor,
            center - y * minor,
        ]));
    }
    if high_sign == Sign::Zero {
        return Ok(IntersectionSet::Points(vec![
            center + x * major,
            center - x * major,
        ]));
    }
    let difference = major2 - minor2;
    let x_offset = x * (major * (low.value / difference.value).sqrt());
    let y_offset = y * (minor * (high.value / difference.value).sqrt());
    Ok(IntersectionSet::Points(vec![
        center + x_offset + y_offset,
        center + x_offset - y_offset,
        center - x_offset + y_offset,
        center - x_offset - y_offset,
    ]))
}

fn circle_circle(a: &Circle, b: &Circle, precision: &Precision) -> Result<IntersectionSet> {
    match round_round(
        a.frame().origin(),
        a.radius(),
        b.frame().origin(),
        b.radius(),
        precision,
    )? {
        RoundSection::Empty => Ok(IntersectionSet::Empty),
        RoundSection::Overlap => Ok(IntersectionSet::Curves(vec![Curve::Circle(*a)])),
        RoundSection::Point(center) => Ok(IntersectionSet::Points(vec![center])),
        RoundSection::Circle {
            center,
            axis,
            radius,
        } => {
            let side = unit(a.frame().z().cross(axis))?;
            Ok(IntersectionSet::Points(vec![
                center + side * radius,
                center - side * radius,
            ]))
        }
    }
}

enum RoundSection {
    Empty,
    Overlap,
    Point(Point3),
    Circle {
        center: Point3,
        axis: Vec3,
        radius: f64,
    },
}

/// Common radical section of two circles or spheres. The factored height avoids
/// subtracting nearly equal squared radii at internal and external tangency.
fn round_round(
    a: Point3,
    r1: f64,
    b: Point3,
    r2: f64,
    precision: &Precision,
) -> Result<RoundSection> {
    let delta = Vector::between(b, a);
    let radii = Scalar::new(r1) - Scalar::new(r2);
    if delta.zero(precision.linear())? {
        return Ok(if radii.sign(precision.linear())? == Sign::Zero {
            RoundSection::Overlap
        } else {
            RoundSection::Empty
        });
    }
    let d2 = delta.norm2();
    let sum = Scalar::new(r1) + Scalar::new(r2);
    let outside = sum.sqr() - d2;
    let inside = d2 - radii.sqr();
    let outer_sign = outside.sign(squared_length_tolerance(r1 + r2, precision))?;
    let inner_sign = inside.sign(squared_length_tolerance(r1.max(r2), precision))?;
    if outer_sign == Sign::Negative || inner_sign == Sign::Negative {
        return Ok(RoundSection::Empty);
    }
    let distance = d2.sqrt();
    let along = (d2 + radii * sum) / (Scalar::new(2.0) * distance);
    let axis = unit(delta.value())?;
    let center = a + axis * along.value;
    if outer_sign == Sign::Zero || inner_sign == Sign::Zero {
        return Ok(RoundSection::Point(center));
    }
    let height2 = outside * inside / (Scalar::new(4.0) * d2);
    Ok(RoundSection::Circle {
        center,
        axis,
        radius: height2.value.sqrt(),
    })
}

fn plane_plane(a: &Plane, b: &Plane, precision: &Precision) -> Result<IntersectionSet> {
    let n1 = a.frame().z();
    let n2 = b.frame().z();
    let direction = n1.cross(n2);
    if parallel(n1, n2, precision)? {
        let signed = Vector::between(b.frame().origin(), a.frame().origin()).dot(Vector::new(n1));
        return Ok(if signed.sign(precision.linear())? == Sign::Zero {
            IntersectionSet::Overlap
        } else {
            IntersectionSet::Empty
        });
    }
    let h = (b.frame().origin() - a.frame().origin()).dot(n2);
    let scale = direction.length_squared();
    let point = a.frame().origin() + direction.cross(n1) * (h / scale);
    Ok(IntersectionSet::Curves(vec![line_from(
        point, direction, precision,
    )?]))
}

fn plane_sphere(plane: &Plane, sphere: &Sphere, precision: &Precision) -> Result<IntersectionSet> {
    let normal = plane.frame().z();
    let signed =
        Vector::between(sphere.frame().origin(), plane.frame().origin()).dot(Vector::new(normal));
    let normal2 = Vector::new(normal).norm2();
    let radius = sphere.radius();
    let radial2 = Scalar::new(radius).sqr() - signed.sqr() / normal2;
    let sign = radial2.sign(squared_length_tolerance(radius, precision))?;
    if sign == Sign::Negative {
        return Ok(IntersectionSet::Empty);
    }
    let center = sphere.frame().origin() - normal * (signed.value / normal2.value);
    if sign == Sign::Zero {
        return Ok(IntersectionSet::Points(vec![center]));
    }
    Ok(IntersectionSet::Curves(vec![circle_in(
        center,
        normal,
        plane.frame().x(),
        radial2.value.sqrt(),
        precision,
    )?]))
}

fn plane_cylinder(
    plane: &Plane,
    cylinder: &Cylinder,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let axis = cylinder.frame().z();
    let normal = plane.frame().z();
    let radius = cylinder.radius();
    let denom = Vector::new(axis).dot(Vector::new(normal));
    if parallel(axis, normal, precision)? {
        let t = (plane.frame().origin() - cylinder.frame().origin()).dot(normal) / denom.value;
        let center = cylinder.frame().origin() + axis * t;
        return Ok(IntersectionSet::Curves(vec![circle_in(
            center,
            normal,
            plane.frame().x(),
            radius,
            precision,
        )?]));
    }
    if denom.sign(precision.angular())? == Sign::Zero {
        let signed = Vector::between(cylinder.frame().origin(), plane.frame().origin())
            .dot(Vector::new(normal));
        let normal2 = Vector::new(normal).norm2();
        let half2 = Scalar::new(radius).sqr() - signed.sqr() / normal2;
        let sign = half2.sign(squared_length_tolerance(radius, precision))?;
        if sign == Sign::Negative {
            return Ok(IntersectionSet::Empty);
        }
        let lateral = axis.cross(normal);
        let on_plane = cylinder.frame().origin() - normal * (signed.value / normal2.value);
        if sign == Sign::Zero {
            return Ok(IntersectionSet::Curves(vec![line_from(
                on_plane, axis, precision,
            )?]));
        }
        let side = unit(lateral)?;
        let half = half2.value.sqrt();
        return Ok(IntersectionSet::Curves(vec![
            line_from(on_plane + side * half, axis, precision)?,
            line_from(on_plane - side * half, axis, precision)?,
        ]));
    }
    let t = (plane.frame().origin() - cylinder.frame().origin()).dot(normal) / denom.value;
    let center = cylinder.frame().origin() + axis * t;
    let minor_dir = axis.cross(normal);
    let major_dir = normal.cross(minor_dir);
    let major = radius / denom.value.abs();
    if major > precision.size() {
        return Err(Error::unsupported(
            "plane-cylinder ellipse exceeds the session size box",
        ));
    }
    let frame = Frame::new(center, normal, unit(major_dir)?, precision)?;
    let ellipse = Ellipse::new(frame, major, radius, precision)?;
    Ok(IntersectionSet::Curves(vec![Curve::Ellipse(ellipse)]))
}

fn plane_torus(plane: &Plane, torus: &Torus, precision: &Precision) -> Result<IntersectionSet> {
    let axis = torus.frame().z();
    let normal = plane.frame().z();
    if !parallel(axis, normal, precision)? {
        let direction = Vector::new(axis).dot(Vector::new(normal));
        let signed = Vector::between(torus.frame().origin(), plane.frame().origin())
            .dot(Vector::new(normal));
        if direction.sign(precision.angular())? == Sign::Zero
            && signed.sign(precision.linear())? == Sign::Zero
        {
            // A meridional plane through the axis contains the two opposite
            // generating tube circles. Ring tori keep those circles disjoint.
            let (major, minor) = torus.radii();
            let radial = unit(normal.cross(axis))?;
            let center = torus.frame().origin();
            return Ok(IntersectionSet::Curves(vec![
                circle_in(center + radial * major, normal, axis, minor, precision)?,
                circle_in(center - radial * major, normal, axis, minor, precision)?,
            ]));
        }
        return Err(Error::unimplemented(
            "plane-torus intersection requires a plane perpendicular to the axis or a meridional plane through it",
        ));
    }
    let (major, minor) = torus.radii();
    let height = Vector::between(plane.frame().origin(), torus.frame().origin())
        .dot(Vector::new(normal))
        / Vector::new(axis).dot(Vector::new(normal));
    let chord2 = Scalar::new(minor).sqr() - height.sqr();
    let sign = chord2.sign(squared_length_tolerance(minor, precision))?;
    if sign == Sign::Negative {
        return Ok(IntersectionSet::Empty);
    }
    let center = torus.frame().origin() + axis * height.value;
    if sign == Sign::Zero {
        return Ok(IntersectionSet::Curves(vec![circle_in(
            center,
            axis,
            torus.frame().x(),
            major,
            precision,
        )?]));
    }
    let chord = chord2.value.sqrt();
    Ok(IntersectionSet::Curves(vec![
        circle_in(center, axis, torus.frame().x(), major + chord, precision)?,
        circle_in(center, axis, torus.frame().x(), major - chord, precision)?,
    ]))
}

fn plane_cone(plane: &Plane, cone: &Cone, precision: &Precision) -> Result<IntersectionSet> {
    let axis = cone.frame().z();
    let normal = plane.frame().z();
    if !parallel(axis, normal, precision)? {
        return cone_plane::intersect(plane, cone, precision);
    }
    let (radius, semi_angle) = cone.shape();
    let height = Vector::between(plane.frame().origin(), cone.frame().origin())
        .dot(Vector::new(normal))
        / Vector::new(axis).dot(Vector::new(normal));
    let slice = Scalar::new(radius) + height * Scalar::new(semi_angle.tan());
    let center = cone.frame().origin() + axis * height.value;
    if slice.sign(precision.linear())? == Sign::Zero {
        return Ok(IntersectionSet::Points(vec![center]));
    }
    Ok(IntersectionSet::Curves(vec![circle_in(
        center,
        axis,
        cone.frame().x(),
        slice.value.abs(),
        precision,
    )?]))
}

fn sphere_sphere(a: &Sphere, b: &Sphere, precision: &Precision) -> Result<IntersectionSet> {
    match round_round(
        a.frame().origin(),
        a.radius(),
        b.frame().origin(),
        b.radius(),
        precision,
    )? {
        RoundSection::Empty => Ok(IntersectionSet::Empty),
        RoundSection::Overlap => Ok(IntersectionSet::Overlap),
        RoundSection::Point(center) => Ok(IntersectionSet::Points(vec![center])),
        RoundSection::Circle {
            center,
            axis,
            radius,
        } => {
            let frame = Frame::from_axis(center, axis, precision)?;
            Ok(IntersectionSet::Curves(vec![circle_in(
                center,
                axis,
                frame.x(),
                radius,
                precision,
            )?]))
        }
    }
}

/// Parallel cylinders reduce to a complete circle-circle section; skew axes need
/// intersection curves beyond the line/circle/ellipse representation.
fn cylinder_cylinder(a: &Cylinder, b: &Cylinder, precision: &Precision) -> Result<IntersectionSet> {
    let axis = a.frame().z();
    if !parallel(axis, b.frame().z(), precision)? {
        return Err(Error::unimplemented(
            "cylinder-cylinder intersection requires parallel axes",
        ));
    }
    let delta = Vector::between(b.frame().origin(), a.frame().origin());
    let axial = delta.dot(Vector::new(axis)) / Vector::new(axis).norm2();
    let center_b = b.frame().origin() - axis * axial.value;
    let circle_a = Circle::new(*a.frame(), a.radius(), precision)?;
    let circle_b = Circle::new(
        Frame::new(center_b, axis, a.frame().x(), precision)?,
        b.radius(),
        precision,
    )?;
    match circle_circle(&circle_a, &circle_b, precision)? {
        IntersectionSet::Empty => Ok(IntersectionSet::Empty),
        IntersectionSet::Points(points) => Ok(IntersectionSet::Curves(
            points
                .into_iter()
                .map(|point| line_from(point, axis, precision))
                .collect::<Result<Vec<_>>>()?,
        )),
        IntersectionSet::Curves(_) | IntersectionSet::Overlap => Ok(IntersectionSet::Overlap),
    }
}

/// A sphere centered on a cylinder axis has zero, one, or two latitude circles.
fn sphere_cylinder(
    sphere: &Sphere,
    cylinder: &Cylinder,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let axis = cylinder.frame().z();
    let delta = Vector::between(sphere.frame().origin(), cylinder.frame().origin());
    if !delta.cross(Vector::new(axis)).zero(precision.linear())? {
        return Err(Error::unimplemented(
            "sphere-cylinder intersection requires a sphere center on the cylinder axis",
        ));
    }
    let height2 = Scalar::new(sphere.radius()).sqr() - Scalar::new(cylinder.radius()).sqr();
    let sign = height2.sign(squared_length_tolerance(sphere.radius(), precision))?;
    match sign {
        Sign::Negative => Ok(IntersectionSet::Empty),
        Sign::Zero => Ok(IntersectionSet::Curves(vec![circle_in(
            sphere.frame().origin(),
            axis,
            cylinder.frame().x(),
            cylinder.radius(),
            precision,
        )?])),
        Sign::Positive => {
            let offset = axis * height2.value.sqrt();
            Ok(IntersectionSet::Curves(vec![
                circle_in(
                    sphere.frame().origin() - offset,
                    axis,
                    cylinder.frame().x(),
                    cylinder.radius(),
                    precision,
                )?,
                circle_in(
                    sphere.frame().origin() + offset,
                    axis,
                    cylinder.frame().x(),
                    cylinder.radius(),
                    precision,
                )?,
            ]))
        }
    }
}

/// Rotationally symmetric pairs reduce to two meridian circles. A ring torus's
/// generating circle lies strictly away from the axis, so each meridian point
/// generates a complete latitude circle and no negative-radius duplicate.
fn sphere_torus(sphere: &Sphere, torus: &Torus, precision: &Precision) -> Result<IntersectionSet> {
    let axis = Vector::new(torus.frame().z());
    let delta = Vector::between(sphere.frame().origin(), torus.frame().origin());
    if !delta.cross(axis).zero(precision.linear())? {
        return sphere_torus::off_axis(sphere, torus, precision);
    }
    let height = delta.dot(axis) / axis.norm2();
    let (major, minor) = torus.radii();
    let sphere_section = Circle::new(
        Frame::new(
            Point3::new(0.0, 0.0, height.value),
            Vec3::Y,
            Vec3::X,
            precision,
        )?,
        sphere.radius(),
        precision,
    )?;
    let tube_section = Circle::new(
        Frame::new(Point3::new(major, 0.0, 0.0), Vec3::Y, Vec3::X, precision)?,
        minor,
        precision,
    )?;
    match circle_circle(&sphere_section, &tube_section, precision)? {
        IntersectionSet::Empty => Ok(IntersectionSet::Empty),
        IntersectionSet::Points(points) => {
            let curves = points
                .into_iter()
                .map(|point| {
                    circle_in(
                        torus.frame().origin() + torus.frame().z() * point.z,
                        torus.frame().z(),
                        torus.frame().x(),
                        point.x,
                        precision,
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(IntersectionSet::Curves(curves))
        }
        IntersectionSet::Curves(_) | IntersectionSet::Overlap => Err(Error::internal(
            "a ring torus meridian cannot coincide with an axis-centered circle",
        )),
    }
}

#[cfg(test)]
mod tests {
    use gensolid_math::Frame;

    use super::*;

    fn frames() -> Vec<Frame> {
        let precision = Precision::DEFAULT;
        let specs = [
            (Point3::ORIGIN, Vec3::Z),
            (Point3::new(4.0, -2.0, 0.5), Vec3::new(0.2, 0.3, 0.9)),
            (Point3::new(-3.0, 1.0, 2.0), Vec3::new(-0.4, 0.8, 0.2)),
            (Point3::new(0.0, 5.0, -1.0), Vec3::new(1.0, -0.2, 0.3)),
        ];
        specs
            .into_iter()
            .map(|(origin, axis)| Frame::from_axis(origin, axis, &precision).unwrap())
            .collect()
    }

    #[test]
    #[allow(clippy::float_cmp)] // Completeness is a ratio of equal integer counts.
    fn analytic_pairs_are_complete_and_accurate() {
        let precision = Precision::DEFAULT;
        let score = analytic_intersection_score(&frames(), &precision).unwrap();
        assert_eq!(score.matched, score.pairs, "{score:?}");
        assert_eq!(score.completeness(), 1.0);
        assert!(score.max_point_error <= precision.linear(), "{score:?}");
        assert!(
            score.pairs >= 32,
            "expected several frames times the pair catalogue"
        );
    }

    #[test]
    fn refuses_incomplete_pairs() {
        let precision = Precision::DEFAULT;
        let plane = Surface::Plane(Plane::new(Frame::WORLD));
        let torus = Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0, &precision).unwrap());
        let tilted = Surface::Plane(Plane::new(
            Frame::new(
                Point3::ORIGIN,
                Vec3::new(0.0, 1.0, 1.0),
                Vec3::X,
                &precision,
            )
            .unwrap(),
        ));
        assert_eq!(
            intersect_surfaces(&tilted, &torus, &precision)
                .unwrap_err()
                .code(),
            ErrorCode::Unimplemented
        );
        assert!(intersect_surfaces(&plane, &torus, &precision).is_ok());
    }
}

#[cfg(test)]
mod robustness_tests;
