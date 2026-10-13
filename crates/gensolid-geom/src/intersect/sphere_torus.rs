//! Complete off-axis sphere/ring-torus cases without tracing a quartic curve.
//!
//! A ring torus is the radius-r tube about its major circle of radius R. For a
//! center at cylindrical coordinates (rho, h), the extreme distances to the
//! major circle are d_minus = hypot(rho - R, h) and
//! d_plus = hypot(rho + R, h). The torus surface distances range continuously
//! from |d_minus - r| to d_plus + r. This proves disjointness and isolated
//! extremal contacts independently of any curve sampling.
//!
//! One additional complete family consists of meridian circles: h = 0 and
//! s^2 = rho^2 - R^2 + r^2. Substituting the torus parameterization into the
//! sphere equation factors its residual into
//! 2 (R + r cos(v)) (R - rho cos(u)). Since R > r, the first factor is strictly
//! positive. Thus rho = R gives one full circle, and rho > R gives exactly two.

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Sign, Vec3};

use super::certified::{Scalar, Vector};
use super::{IntersectionSet, circle_in, unit};
use crate::surface::{Sphere, Torus};

/// Called only after the rotationally symmetric implementation excludes the axis.
pub(super) fn off_axis(
    sphere: &Sphere,
    torus: &Torus,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let axis = Vector::new(torus.frame().z());
    let delta = Vector::between(sphere.frame().origin(), torus.frame().origin());
    let axis_squared = axis.norm2();
    let projection = delta.dot(axis);
    let height = projection / axis_squared.sqrt();
    let radial_delta = delta - axis * (projection / axis_squared);
    let rho = radial_delta.norm2().sqrt();
    let (major_value, minor_value) = torus.radii();
    let (major, minor, radius) = (
        Scalar::new(major_value),
        Scalar::new(minor_value),
        Scalar::new(sphere.radius()),
    );
    let near_generator = ((rho - major).sqr() + height.sqr()).sqrt();
    let far_generator = ((rho + major).sqr() + height.sqr()).sqrt();
    // Scalar::sign first uses a deterministic interval sign. Only an interval
    // containing zero uses its explicit tolerance fallback. Consequently every
    // Positive below certifies strict separation, including sub-tolerance gaps.
    let contact_tolerance = precision.linear() * 0.5;
    let outside = (near_generator - minor - radius).sign(contact_tolerance)?;
    let inside = (minor - near_generator - radius).sign(contact_tolerance)?;
    let enclosing = (radius - far_generator - minor).sign(contact_tolerance)?;
    if outside == Sign::Positive || inside == Sign::Positive || enclosing == Sign::Positive {
        return Ok(IntersectionSet::Empty);
    }

    let radial_direction = unit(radial_delta.value())?;
    let axis_direction = unit(torus.frame().z())?;
    if let Some(circles) = meridian_circles(
        torus,
        radial_direction,
        axis_direction,
        rho,
        height,
        radius,
        precision,
    )? {
        return Ok(circles);
    }
    if enclosing == Sign::Zero {
        let generator = torus.frame().origin() - radial_direction * major_value;
        let point =
            generator + (generator - sphere.frame().origin()) * (minor_value / far_generator.value);
        return Ok(IntersectionSet::Points(vec![point]));
    }
    if outside == Sign::Zero || inside == Sign::Zero {
        if (near_generator - Scalar::new(precision.linear())).sign(contact_tolerance)?
            != Sign::Positive
        {
            return Err(Error::unimplemented(
                "unresolved sphere-torus contact at a major-circle center",
            ));
        }
        let generator = torus.frame().origin() + radial_direction * major_value;
        let point = generator
            + (sphere.frame().origin() - generator) * (minor_value / near_generator.value);
        return Ok(IntersectionSet::Points(vec![point]));
    }
    Err(Error::unimplemented(
        "off-axis sphere-torus intersection requires a general quartic curve representation",
    ))
}

#[allow(clippy::too_many_arguments)]
fn meridian_circles(
    torus: &Torus,
    radial: Vec3,
    axis: Vec3,
    rho: Scalar,
    height: Scalar,
    sphere_radius: Scalar,
    precision: &Precision,
) -> Result<Option<IntersectionSet>> {
    let tolerance = precision.linear() * 0.5;
    if height.sign(tolerance)? != Sign::Zero {
        return Ok(None);
    }
    let (major_value, minor_value) = torus.radii();
    let major = Scalar::new(major_value);
    let excess = rho - major;
    let position = excess.sign(tolerance)?;
    if position == Sign::Negative {
        return Ok(None);
    }
    // Factoring rho^2 - R^2 preserves small positive excess near the one-circle
    // limit. The two half-tolerance budgets bound the axial and radius residuals.
    let expected_radius = (excess * (rho + major) + Scalar::new(minor_value).sqr()).sqrt();
    if (sphere_radius - expected_radius).sign(tolerance)? != Sign::Zero {
        return Ok(None);
    }
    let meridian = |direction: Vec3| {
        circle_in(
            torus.frame().origin() + direction * major_value,
            axis.cross(direction),
            direction,
            minor_value,
            precision,
        )
    };
    if position == Sign::Zero {
        return Ok(Some(IntersectionSet::Curves(vec![meridian(radial)?])));
    }
    let cosine = (major / rho).value;
    let sine = ((excess * (rho + major)).sqrt() / rho).value;
    let side = axis.cross(radial);
    Ok(Some(IntersectionSet::Curves(vec![
        meridian(radial * cosine + side * sine)?,
        meridian(radial * cosine - side * sine)?,
    ])))
}

#[cfg(test)]
mod tests {
    use std::f64::consts::TAU;

    use gensolid_math::{Frame, Point3};
    use proptest::prelude::*;

    use super::*;
    use crate::{GeomStore, Surface, intersect_surfaces};

    fn pair(frame: Frame, center: Vec3, radius: f64, major: f64, minor: f64) -> (Surface, Surface) {
        let precision = Precision::DEFAULT;
        let sphere_frame =
            Frame::new(frame.to_world(center), frame.z(), frame.x(), &precision).unwrap();
        (
            Surface::Sphere(Sphere::new(sphere_frame, radius, &precision).unwrap()),
            Surface::Torus(Torus::new(frame, major, minor, &precision).unwrap()),
        )
    }

    fn check_residuals(hit: &IntersectionSet, sphere: &Surface, torus: &Surface) {
        let precision = Precision::DEFAULT;
        let mut store = GeomStore::new();
        let ids = [
            store.add_surface(sphere.clone()).unwrap(),
            store.add_surface(torus.clone()).unwrap(),
        ];
        let points = match hit {
            IntersectionSet::Points(points) => points.clone(),
            IntersectionSet::Curves(curves) => curves
                .iter()
                .flat_map(|curve| (0..32).map(move |i| curve.eval(TAU * i as f64 / 32.0).point))
                .collect(),
            _ => vec![],
        };
        for point in points {
            for id in ids {
                let residual = store.distance_to_surface(id, point, &precision).unwrap();
                assert!(
                    residual <= precision.linear(),
                    "residual={residual}, hit={hit:?}"
                );
            }
        }
    }

    fn frames() -> impl Strategy<Value = Frame> {
        (
            prop::array::uniform3(-10.0f64..10.0),
            prop::array::uniform3(-1.0f64..1.0),
        )
            .prop_filter_map("finite nonzero axis", |(origin, axis)| {
                Frame::from_axis(
                    Point3::new(origin[0], origin[1], origin[2]),
                    Vec3::new(axis[0], axis[1], axis[2]),
                    &Precision::DEFAULT,
                )
                .ok()
            })
    }

    /// Exact signed-permutation bases keep the special family's axial equation
    /// exactly zero. Arbitrary floating rotations do not retain the provenance
    /// of a constructed center and may correctly become a general quartic pair.
    fn meridian_frames() -> impl Strategy<Value = Frame> {
        (
            prop::array::uniform3(-10.0f64..10.0),
            0usize..3,
            any::<bool>(),
            any::<bool>(),
        )
            .prop_map(|(origin, axis_index, flip_axis, flip_reference)| {
                let directions = [Vec3::X, Vec3::Y, Vec3::Z];
                let axis = directions[axis_index] * if flip_axis { -1.0 } else { 1.0 };
                let reference = directions[(axis_index + 1) % directions.len()]
                    * if flip_reference { -1.0 } else { 1.0 };
                Frame::new(
                    Point3::new(origin[0], origin[1], origin[2]),
                    axis,
                    reference,
                    &Precision::DEFAULT,
                )
                .unwrap()
            })
    }

    #[test]
    fn strict_distance_bounds_cover_hole_tube_exterior_and_enclosing_spheres() {
        let precision = Precision::DEFAULT;
        for (center, radius) in [
            (Vec3::new(0.5, 0.0, 0.0), 0.5),
            (Vec3::new(6.0, 0.0, 4.0), 0.5),
            (Vec3::new(3.0, 0.0, 0.0), 0.5),
            (Vec3::new(0.5, 0.0, 0.0), 6.0),
        ] {
            let (sphere, torus) = pair(Frame::WORLD, center, radius, 3.0, 1.0);
            assert_eq!(
                intersect_surfaces(&sphere, &torus, &precision).unwrap(),
                IntersectionSet::Empty
            );
            assert_eq!(
                intersect_surfaces(&torus, &sphere, &precision).unwrap(),
                IntersectionSet::Empty
            );
        }
    }

    #[test]
    fn nearest_and_farthest_contacts_are_complete_isolated_points() {
        let precision = Precision::DEFAULT;
        for (center, radius) in [
            (Vec3::new(1.0, 0.0, 0.0), 1.0),
            (Vec3::new(6.0, 0.0, 0.0), 2.0),
            (Vec3::new(3.5, 0.0, 0.0), 0.5),
            (Vec3::new(1.0, 0.0, 0.0), 5.0),
        ] {
            let (sphere, torus) = pair(Frame::WORLD, center, radius, 3.0, 1.0);
            let hit = intersect_surfaces(&sphere, &torus, &precision).unwrap();
            assert!(matches!(&hit, IntersectionSet::Points(points) if points.len() == 1));
            check_residuals(&hit, &sphere, &torus);
        }
    }

    #[test]
    fn quartic_branches_are_not_replaced_with_a_partial_circle() {
        let precision = Precision::DEFAULT;
        let (sphere, torus) = pair(Frame::WORLD, Vec3::new(1.0, 0.0, 0.0), 2.0, 3.0, 1.0);
        assert_eq!(
            intersect_surfaces(&sphere, &torus, &precision)
                .unwrap_err()
                .code(),
            gensolid_core::ErrorCode::Unimplemented
        );
        let (sphere, torus) = pair(
            Frame::WORLD,
            Vec3::new(3.0, 0.0, precision.linear() * 0.25),
            1.0,
            3.0,
            1.0,
        );
        assert_eq!(
            intersect_surfaces(&sphere, &torus, &precision)
                .unwrap_err()
                .code(),
            gensolid_core::ErrorCode::Unimplemented
        );
    }

    #[test]
    fn meridian_circle_family_keeps_its_one_and_two_branch_limits() {
        let precision = Precision::DEFAULT;
        for (rho, branches) in [(3.0, 1), (3.0 + precision.linear() * 4.0, 2), (4.0, 2)] {
            let radius = ((rho - 3.0) * (rho + 3.0) + 1.0).sqrt();
            let (sphere, torus) = pair(Frame::WORLD, Vec3::new(rho, 0.0, 0.0), radius, 3.0, 1.0);
            let hit = intersect_surfaces(&sphere, &torus, &precision).unwrap();
            assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == branches));
            check_residuals(&hit, &sphere, &torus);
        }
    }

    #[test]
    fn near_extremal_boundaries_keep_certified_gaps() {
        let precision = Precision::DEFAULT;
        for (radius, empty) in [
            (1.0 - precision.linear() * 4.0, true),
            (1.0 + precision.linear() * 4.0, false),
            (5.0 - precision.linear() * 4.0, false),
            (5.0 + precision.linear() * 4.0, true),
        ] {
            let (sphere, torus) = pair(Frame::WORLD, Vec3::new(1.0, 0.0, 0.0), radius, 3.0, 1.0);
            let result = intersect_surfaces(&sphere, &torus, &precision);
            if empty {
                assert_eq!(result.unwrap(), IntersectionSet::Empty);
            } else {
                assert_eq!(
                    result.unwrap_err().code(),
                    gensolid_core::ErrorCode::Unimplemented
                );
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]
        #[test]
        fn random_frames_preserve_hole_separation(frame in frames(), major in 1.0f64..20.0, minor_fraction in 0.1f64..0.4, axial in -2.0f64..2.0) {
            let minor = major * minor_fraction;
            let hole = major - minor;
            let (sphere, torus) = pair(frame, Vec3::new(hole * 0.25, 0.0, axial), hole * 0.25, major, minor);
            let hit = intersect_surfaces(&sphere, &torus, &Precision::DEFAULT).unwrap();
            prop_assert_eq!(hit, IntersectionSet::Empty);
        }

        #[test]
        fn meridian_family_returns_both_complete_circles(frame in meridian_frames(), major in 1.0f64..10.0, minor_fraction in 0.1f64..0.4, radial_factor in 1.1f64..3.0) {
            let minor = major * minor_fraction;
            let rho = major * radial_factor;
            let stored_rho = frame.to_local(frame.to_world(Vec3::new(rho, 0.0, 0.0))).x;
            let radius = ((stored_rho - major) * (stored_rho + major) + minor * minor).sqrt();
            let (sphere, torus) = pair(frame, Vec3::new(rho, 0.0, 0.0), radius, major, minor);
            let hit = intersect_surfaces(&sphere, &torus, &Precision::DEFAULT).unwrap();
            prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == 2));
            check_residuals(&hit, &sphere, &torus);
        }

        #[test]
        fn near_hole_contact_in_random_frames_does_not_become_false_empty(frame in frames(), boundary_factor in 2.0f64..32.0) {
            let precision = Precision::DEFAULT;
            let gap = precision.linear() * boundary_factor;
            for (radius, empty) in [(1.0 - gap, true), (1.0 + gap, false)] {
                let (sphere, torus) = pair(frame, Vec3::new(1.0, 0.0, 0.0), radius, 3.0, 1.0);
                let result = intersect_surfaces(&sphere, &torus, &precision);
                if empty { prop_assert_eq!(result.unwrap(), IntersectionSet::Empty); }
                else { prop_assert_eq!(result.unwrap_err().code(), gensolid_core::ErrorCode::Unimplemented); }
            }
        }
    }
}
