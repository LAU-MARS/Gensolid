//! Owned, native Rust surface construction for the geometry backend boundary.
//!
//! These are original implementations of tensor-product NURBS constructions,
//! rather than translations of a foreign kernel. They preserve rational geometry
//! and expose unsupported cases explicitly. No C++ session or pointer is needed.

use std::f64::consts::{FRAC_PI_2, TAU};

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Interval, Point3, Sign, Vec3};

use crate::backend::{GeometryBackend, SurfaceConstruction, SurfaceIntersectionCandidate};
use crate::{BCurve, BSurface};

/// Pure Rust provider of exact NURBS constructions with explicit capability limits.
#[derive(Clone, Copy, Debug)]
pub struct RustGeometryBackend {
    precision: Precision,
}

impl RustGeometryBackend {
    /// Creates an owned Rust backend; there is no foreign runtime to initialize.
    pub const fn new(precision: Precision) -> Self {
        Self { precision }
    }

    /// Precision used by this provider.
    pub const fn precision(&self) -> Precision {
        self.precision
    }
}

impl GeometryBackend for RustGeometryBackend {
    fn precision(&self) -> Precision {
        self.precision
    }

    fn construct_surface(&self, request: &SurfaceConstruction) -> Result<BSurface> {
        match request {
            SurfaceConstruction::Extrude(profile, direction) => {
                extrude_nurbs(profile, *direction, &self.precision)
            }
            SurfaceConstruction::Loft(sections) => loft_compatible(sections, &self.precision),
            SurfaceConstruction::Sweep(profile, path) => {
                sweep_linear(profile, path, &self.precision)
            }
            SurfaceConstruction::Revolve {
                profile,
                origin,
                axis,
                angle,
            } => revolve_nurbs(profile, *origin, *axis, *angle, &self.precision),
        }
    }

    fn intersect_surfaces(
        &self,
        first: &BSurface,
        second: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>> {
        crate::rust_surface_intersection::intersect_extrusion_plane(first, second, &self.precision)
    }

    fn offset_surface(&self, surface: &BSurface, distance: f64) -> Result<BSurface> {
        offset_planar(surface, distance, &self.precision)
    }
}

fn check_displacement(direction: Vec3, p: &Precision) -> Result<()> {
    if !direction.is_finite() {
        return Err(Error::invalid_argument("displacement must be finite"));
    }
    let length = direction.length();
    if !length.is_finite() || length <= p.degenerate_length() {
        return Err(Error::degenerate("displacement is degenerate at precision"));
    }
    Ok(())
}

/// Exactly translates the profile: `S(u,v) = C(u) + v * direction`, `v in [0,1]`.
/// The original U knot vector and rational weights are preserved.
pub fn extrude_nurbs(profile: &BCurve, direction: Vec3, p: &Precision) -> Result<BSurface> {
    check_displacement(direction, p)?;
    let points = profile
        .points()
        .iter()
        .flat_map(|&point| [point, point + direction])
        .collect();
    let weights = profile
        .weights()
        .map(|w| w.iter().flat_map(|&weight| [weight, weight]).collect());
    BSurface::new(
        profile.degree(),
        1,
        profile.knots().to_vec(),
        vec![0.0, 0.0, 1.0, 1.0],
        profile.points().len(),
        2,
        points,
        weights,
    )
}

/// Interpolates compatible sections at uniformly spaced parameters in `[0,1]`.
/// Degree in V is `min(3, section_count - 1)`. U degree, knots and control-point
/// correspondence must agree; this does not silently resample incompatible curves.
pub fn loft_compatible(sections: &[BCurve], p: &Precision) -> Result<BSurface> {
    if !(2..=64).contains(&sections.len()) {
        return Err(Error::invalid_argument("loft requires 2..64 sections"));
    }
    let parameters = (0..sections.len())
        .map(|i| i as f64 / (sections.len() - 1) as f64)
        .collect::<Vec<_>>();
    crate::loft::loft_nurbs_sections(sections, &parameters, 3.min(sections.len() - 1), p)
}

/// Creates a degree-one loft between compatible sections, with homogeneous
/// interpolation between them. Equal corresponding weights give a straight ruling.
pub fn ruled_compatible(first: &BCurve, second: &BCurve, p: &Precision) -> Result<BSurface> {
    crate::loft::loft_nurbs_sections(&[first.clone(), second.clone()], &[0.0, 1.0], 1, p)
}

/// Sweeps a profile by the exact displacement along a straight degree-one path.
/// `S(u,v) = C(u) + P(v) - P(v_min)` retains the path's parameterization, including
/// its positive rational weights. Curved paths require a framing policy and are
/// explicitly unimplemented here.
pub fn sweep_linear(profile: &BCurve, path: &BCurve, p: &Precision) -> Result<BSurface> {
    if path.degree() != 1 || path.points().len() != 2 {
        return Err(Error::unimplemented(
            "Rust sweep supports degree-one paths with two control points",
        ));
    }
    let direction = path.points()[1] - path.points()[0];
    check_displacement(direction, p)?;
    let pw = normalized_weights(profile.weights(), profile.points().len());
    let qw = normalized_weights(path.weights(), path.points().len());
    let points = profile
        .points()
        .iter()
        .flat_map(|&point| [point, point + direction])
        .collect();
    let weights = pw
        .iter()
        .flat_map(|w| qw.iter().map(move |q| w * q))
        .collect();
    BSurface::new(
        profile.degree(),
        path.degree(),
        profile.knots().to_vec(),
        path.knots().to_vec(),
        profile.points().len(),
        path.points().len(),
        points,
        Some(weights),
    )
}

fn normalized_weights(weights: Option<&[f64]>, count: usize) -> Vec<f64> {
    weights.map_or_else(
        || vec![1.0; count],
        |weights| {
            let scale = weights.iter().copied().fold(0.0, f64::max);
            weights.iter().map(|w| w / scale).collect()
        },
    )
}

/// Revolves a positive-weight NURBS profile about an arbitrary oriented axis.
/// V runs from 0 to 1; its rational quadratic parameter is not a linear angle.
/// Every arc span is at most a quarter turn, so its weights remain positive.
/// Signed rotations up to one full turn are supported. The U parameterization
/// and exact circular geometry are retained without sampling or tessellation.
pub fn revolve_nurbs(
    profile: &BCurve,
    origin: Point3,
    axis: Vec3,
    angle: f64,
    p: &Precision,
) -> Result<BSurface> {
    if !origin.is_finite() || !axis.is_finite() || !angle.is_finite() {
        return Err(Error::invalid_argument("revolution inputs must be finite"));
    }
    if angle.abs() <= p.angular() {
        return Err(Error::degenerate("revolution angle is degenerate"));
    }
    if angle.abs() > TAU {
        return Err(Error::unimplemented(
            "Rust revolution supports at most one full turn",
        ));
    }
    let axis = axis
        .normalized(p.degenerate_length())
        .ok_or_else(|| Error::degenerate("revolution axis is degenerate"))?;
    let spans = (angle.abs() / FRAC_PI_2).ceil() as usize;
    let step = angle / spans as f64;
    let middle_weight = (step / 2.0).cos();
    let count = 2 * spans + 1;
    let mut knots = vec![0.0; 3];
    for span in 1..spans {
        knots.extend([span as f64 / spans as f64; 2]);
    }
    knots.extend([1.0; 3]);
    let profile_weights = normalized_weights(profile.weights(), profile.points().len());
    let mut points = Vec::with_capacity(profile.points().len() * count);
    let mut weights = Vec::with_capacity(points.capacity());
    for (&point, weight) in profile.points().iter().zip(profile_weights) {
        let relative = point - origin;
        let center = origin + axis * axis.dot(relative);
        let radial = point - center;
        let tangent = axis.cross(radial);
        for j in 0..count {
            let theta = step * j as f64 / 2.0;
            let (sin, cos) = theta.sin_cos();
            let arc_weight = if j % 2 == 0 { 1.0 } else { middle_weight };
            points.push(center + (radial * cos + tangent * sin) / arc_weight);
            weights.push(weight * arc_weight);
        }
    }
    BSurface::new(
        profile.degree(),
        2,
        profile.knots().to_vec(),
        knots,
        profile.points().len(),
        count,
        points,
        Some(weights),
    )
}

fn cross_dot_interval(
    u_start: Point3,
    u_end: Point3,
    v_start: Point3,
    v_end: Point3,
    normal: Vec3,
) -> Interval {
    let difference = |a: Point3, b: Point3| {
        [
            Interval::point(b.x) - Interval::point(a.x),
            Interval::point(b.y) - Interval::point(a.y),
            Interval::point(b.z) - Interval::point(a.z),
        ]
    };
    let u = difference(u_start, u_end);
    let v = difference(v_start, v_end);
    (u[1] * v[2] - u[2] * v[1]) * Interval::point(normal.x)
        + (u[2] * v[0] - u[0] * v[2]) * Interval::point(normal.y)
        + (u[0] * v[1] - u[1] * v[0]) * Interval::point(normal.z)
}

/// Translates a non-folded bilinear planar patch along its oriented normal.
/// Only polynomial patches (or a common rational weight) are accepted. The four
/// corner cross products enclose all bilinear normal directions: their dot
/// products with the reference normal must have certified positive signs.
/// Coordinate round-off is bounded using session precision; varying normals must
/// remain within the solver budget after the requested offset. Curved, folded,
/// higher-degree and nonuniform-weight surfaces are explicitly unimplemented.
pub fn offset_planar(surface: &BSurface, distance: f64, p: &Precision) -> Result<BSurface> {
    if !distance.is_finite() {
        return Err(Error::invalid_argument("offset distance must be finite"));
    }
    if surface.degrees() != (1, 1) || surface.counts() != (2, 2) {
        return Err(Error::unimplemented(
            "Rust planar offset requires a bilinear control net",
        ));
    }
    if let Some(weights) = surface.weights()
        && weights.iter().any(|w| !w.total_cmp(&weights[0]).is_eq())
    {
        return Err(Error::unimplemented(
            "Rust planar offset requires a common rational weight",
        ));
    }
    let points = surface.points();
    let a = points[2] - points[0];
    let b = points[1] - points[0];
    let cross = a.cross(b);
    let normal = cross
        .normalized(p.angular() * a.length() * b.length())
        .ok_or_else(|| Error::degenerate("planar offset has a singular corner"))?;
    if (points[3] - points[0]).dot(normal).abs() > p.solver() {
        return Err(Error::unimplemented("Rust offset requires a planar patch"));
    }
    for u in 0..2 {
        for v in 0..2 {
            let du = points[2 + v] - points[v];
            let dv = points[2 * u + 1] - points[2 * u];
            let corner = du.cross(dv);
            // Enclose subtraction and cross-product rounding from the original
            // control coordinates, not only a dot product of rounded vectors.
            if cross_dot_interval(
                points[v],
                points[2 + v],
                points[2 * u],
                points[2 * u + 1],
                normal,
            )
            .sign()
                != Some(Sign::Positive)
            {
                return Err(Error::degenerate(
                    "planar offset orientation is folded or unresolved",
                ));
            }
            let unit = corner
                .normalized(p.angular() * du.length() * dv.length())
                .ok_or_else(|| Error::degenerate("planar offset has a singular corner"))?;
            if (unit - normal).length() * distance.abs() > p.solver() {
                return Err(Error::unimplemented(
                    "patch normal varies beyond the planar offset error budget",
                ));
            }
        }
    }
    let (u_knots, v_knots) = surface.knots();
    BSurface::new(
        1,
        1,
        u_knots.to_vec(),
        v_knots.to_vec(),
        2,
        2,
        points.iter().map(|&q| q + normal * distance).collect(),
        surface.weights().map(<[f64]>::to_vec),
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use gensolid_core::ErrorCode;
    use gensolid_math::Frame;
    use proptest::prelude::*;

    fn curve(frame: Frame, size: f64, weight: f64) -> BCurve {
        BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            [
                Vec3::new(size, 0.0, 0.0),
                Vec3::new(size, size, size / 3.0),
                Vec3::new(0.0, size, size),
            ]
            .map(|v| frame.to_world(v))
            .to_vec(),
            Some(vec![1.0, weight, 1.0]),
        )
        .unwrap()
    }

    fn line(a: Point3, b: Point3) -> BCurve {
        BCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![a, b], None).unwrap()
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(48))]

        #[test]
        fn rational_extrusion_and_sweep_are_exact_in_random_frames(
            size in 0.1f64..20.0,
            angle in -3.0f64..3.0,
            translation in -20.0f64..20.0,
            weight in 0.2f64..5.0,
            u in 0.0f64..1.0,
            v in 0.0f64..1.0,
        ) {
            let p = Precision::DEFAULT;
            let frame = Frame::from_axis(Point3::new(translation, -translation, translation),
                Vec3::new(angle.cos(), angle.sin(), 1.0), &p)?;
            let profile = curve(frame, size, weight);
            let direction = frame.z() * size;
            let extrusion = extrude_nurbs(&profile, direction, &p)?;
            let actual = extrusion.eval(u, v);
            let expected = profile.eval(u);
            prop_assert!(actual.point.distance(expected.point + direction * v) <= p.linear());
            prop_assert!((actual.du - expected.d1).length() <= p.linear());
            prop_assert!((actual.dv - direction).length() <= p.linear());
            let path = BCurve::new(1, vec![2.0, 2.0, 5.0, 5.0],
                vec![frame.origin(), frame.origin() + direction], Some(vec![weight, 1.0]))?;
            let swept = sweep_linear(&profile, &path, &p)?;
            let parameter = 2.0 + 3.0 * v;
            prop_assert!(swept.eval(u, parameter).point.distance(
                expected.point + (path.eval(parameter).point - path.eval(2.0).point)) <= p.linear());
        }

        #[test]
        fn revolution_preserves_axis_distance_height_and_endpoints(
            size in 0.1f64..20.0,
            tilt in -3.0f64..3.0,
            translation in -20.0f64..20.0,
            weight in 0.2f64..5.0,
            angle in -6.0f64..6.0,
            u in 0.0f64..1.0,
            v in 0.0f64..1.0,
        ) {
            let p = Precision::DEFAULT;
            prop_assume!(angle.abs() > p.angular());
            let frame = Frame::from_axis(Point3::new(translation, -translation, translation),
                Vec3::new(tilt.cos(), tilt.sin(), 1.0), &p)?;
            let profile = curve(frame, size, weight);
            let surface = revolve_nurbs(&profile, frame.origin(), frame.z() * 7.0, angle, &p)?;
            let original = frame.to_local(profile.eval(u).point);
            let actual = frame.to_local(surface.eval(u, v).point);
            prop_assert!((original.z - actual.z).abs() <= p.linear());
            prop_assert!((original.x.hypot(original.y) - actual.x.hypot(actual.y)).abs() <= p.linear());
            prop_assert!(surface.eval(u, 0.0).point.distance(profile.eval(u).point) <= p.linear());
            let (sin, cos) = angle.sin_cos();
            let end = frame.to_world(Vec3::new(original.x*cos - original.y*sin,
                original.x*sin + original.y*cos, original.z));
            prop_assert!(surface.eval(u, 1.0).point.distance(end) <= p.linear());
        }

        #[test]
        fn planar_offset_preserves_signed_distance_under_rotation(
            width in 0.1f64..20.0,
            height in 0.1f64..20.0,
            tilt in -3.0f64..3.0,
            translation in -20.0f64..20.0,
            distance in -5.0f64..5.0,
            u in 0.0f64..1.0,
            v in 0.0f64..1.0,
        ) {
            let p = Precision::DEFAULT;
            let frame = Frame::from_axis(Point3::new(translation, -translation, translation),
                Vec3::new(tilt.cos(), tilt.sin(), 1.0), &p)?;
            let profile = line(frame.origin(), frame.origin() + frame.x()*width);
            let plane = extrude_nurbs(&profile, frame.y()*height, &p)?;
            let offset = offset_planar(&plane, distance, &p)?;
            let actual = offset.eval(u,v).point - plane.eval(u,v).point;
            prop_assert!((actual - frame.z()*distance).length() <= p.linear());
            prop_assert_eq!(offset.knots(), plane.knots());
        }

        #[test]
        fn precision_boundary_extrusions_and_revolutions_fail(
            factor in 0.01f64..0.99,
            shift in -20.0f64..20.0,
        ) {
            let p = Precision::DEFAULT;
            let profile = line(Point3::new(shift, 0.0, 0.0), Point3::new(shift+1.0, 0.0, 0.0));
            prop_assert_eq!(extrude_nurbs(&profile, Vec3::Z*p.linear()*factor, &p)
                .err().map(|e|e.code()), Some(ErrorCode::Degenerate));
            prop_assert_eq!(revolve_nurbs(&profile, Point3::ORIGIN, Vec3::Z, p.angular()*factor, &p)
                .err().map(|e|e.code()), Some(ErrorCode::Degenerate));
            prop_assert_eq!(revolve_nurbs(&profile, Point3::ORIGIN, Vec3::Z*p.linear()*factor, 1.0, &p)
                .err().map(|e|e.code()), Some(ErrorCode::Degenerate));
        }
    }

    #[test]
    fn compatible_loft_interpolates_all_sections_and_rejects_mismatches() -> Result<()> {
        let p = Precision::DEFAULT;
        let sections = (0..4)
            .map(|i| {
                line(
                    Point3::new(0.0, f64::from(i), 0.0),
                    Point3::new(2.0, f64::from(i), 1.0),
                )
            })
            .collect::<Vec<_>>();
        let surface = loft_compatible(&sections, &p)?;
        for (i, section) in sections.iter().enumerate() {
            for u in [0.0, 0.2, 0.7, 1.0] {
                assert!(
                    surface
                        .eval(u, i as f64 / 3.0)
                        .point
                        .distance(section.eval(u).point)
                        <= p.linear()
                );
            }
        }
        let incompatible = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::ORIGIN,
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ],
            None,
        )?;
        assert!(loft_compatible(&[sections[0].clone(), incompatible.clone()], &p).is_err());
        assert!(sweep_linear(&sections[0], &incompatible, &p).is_err());
        Ok(())
    }

    #[test]
    fn full_turn_closes_and_invalid_inputs_never_produce_a_surface() -> Result<()> {
        let p = Precision::DEFAULT;
        let profile = line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 3.0));
        for angle in [TAU, -TAU] {
            let surface = revolve_nurbs(&profile, Point3::ORIGIN, Vec3::Z, angle, &p)?;
            for u in [0.0, 0.5, 1.0] {
                assert!(
                    surface
                        .eval(u, 0.0)
                        .point
                        .distance(surface.eval(u, 1.0).point)
                        <= p.linear()
                );
            }
        }
        assert!(revolve_nurbs(&profile, Point3::ORIGIN, Vec3::ZERO, 1.0, &p).is_err());
        assert!(revolve_nurbs(&profile, Point3::ORIGIN, Vec3::Z, f64::NAN, &p).is_err());
        assert!(revolve_nurbs(&profile, Point3::ORIGIN, Vec3::Z, TAU * 2.0, &p).is_err());
        assert!(extrude_nurbs(&profile, Vec3::new(f64::INFINITY, 0.0, 0.0), &p).is_err());
        Ok(())
    }

    #[test]
    fn folded_and_curved_offsets_are_explicitly_rejected() -> Result<()> {
        let p = Precision::DEFAULT;
        let patch = |last: Point3| {
            BSurface::new(
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                2,
                2,
                vec![
                    Point3::ORIGIN,
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    last,
                ],
                None,
            )
        };
        assert_eq!(
            offset_planar(&patch(Point3::new(-1.0, -1.0, 0.0))?, 1.0, &p)
                .err()
                .map(|e| e.code()),
            Some(ErrorCode::Degenerate)
        );
        assert_eq!(
            offset_planar(&patch(Point3::new(1.0, 1.0, 1.0))?, 1.0, &p)
                .err()
                .map(|e| e.code()),
            Some(ErrorCode::Unimplemented)
        );
        assert!(offset_planar(&patch(Point3::new(1.0, 1.0, 0.0))?, f64::INFINITY, &p).is_err());
        Ok(())
    }

    #[test]
    fn unresolved_nearly_parallel_corner_is_not_assumed_regular() -> Result<()> {
        let base = Precision::DEFAULT;
        let p = Precision::new(base.linear(), base.angular() * base.angular(), base.size())?;
        // The rounded cross products are positive, but full outward-rounded
        // arithmetic cannot certify their signs. Dotting an already rounded
        // cross product with an Interval would incorrectly accept this patch.
        let patch = BSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            2,
            2,
            vec![
                Point3::ORIGIN,
                Point3::new(1.0, 1.0_f64.next_up(), 0.0),
                Point3::new(0.5, 0.5, 0.0),
                Point3::new(1.5, 1.5_f64.next_up(), 0.0),
            ],
            None,
        )?;
        assert_eq!(
            offset_planar(&patch, 1.0, &p).err().map(|e| e.code()),
            Some(ErrorCode::Degenerate)
        );
        Ok(())
    }
}
