//! Complete representable sections of a full double cone by an oblique plane.
//!
//! In an orthonormal plane basis with x along the projected cone axis, the
//! cone equation is A x^2 + y^2 + B x + C = 0, where
//! A = q^2 - k^2 t^2, q is the axis/normal cosine, t its sine, and k the cone
//! slope. A > 0 is a complete ellipse on one nappe. Through the apex the same
//! equation gives one point, one generator, or two generators according to A.
//! Nondegenerate parabolas and hyperbolas need an unbounded conic representation.

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Frame, Sign};

use super::certified::{Scalar, Vector};
use super::{IntersectionSet, circle_in, line_from, unit};
use crate::{Cone, Curve, Ellipse, Plane};

pub(super) fn intersect(
    plane: &Plane,
    cone: &Cone,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let normal = Vector::new(plane.frame().z());
    let axis = Vector::new(cone.frame().z());
    let normal_norm = normal.norm2().sqrt();
    let axis_norm = axis.norm2().sqrt();
    let q = axis.dot(normal) / (axis_norm * normal_norm);
    let t = axis.cross(normal).norm2().sqrt() / (axis_norm * normal_norm);
    let (radius, semi_angle) = cone.shape();
    let slope = Scalar::new(semi_angle.tan());
    let apex_offset = Scalar::new(radius) / slope;
    // Keep the apex construction in the arithmetic expression: rebuilding a
    // Scalar from a rounded apex coordinate would erase cancellation bounds.
    let to_plane =
        Vector::between(plane.frame().origin(), cone.frame().origin()) - axis * (-apex_offset);
    let height = to_plane.dot(normal) / normal_norm;
    let apex = cone.frame().origin() - cone.frame().z() * apex_offset.value;
    let coefficient = q.sqr() - slope.sqr() * t.sqr();
    let angle_tolerance = precision.angular() * (1.0 + slope.value.abs()).powi(2);
    let conic = coefficient.sign(angle_tolerance)?;
    let through_apex = height.sign(precision.linear())? == Sign::Zero;
    let normal_direction = unit(plane.frame().z())?;
    let axis_direction = unit(cone.frame().z())?;
    if t.sign(precision.angular())? == Sign::Zero {
        if through_apex {
            return Ok(IntersectionSet::Points(vec![apex]));
        }
        let center = apex + normal_direction * height.value;
        return Ok(IntersectionSet::Curves(vec![circle_in(
            center,
            normal_direction,
            cone.frame().x(),
            (slope * height).value.abs(),
            precision,
        )?]));
    }
    let major_direction = unit(
        normal_direction
            .cross(axis_direction)
            .cross(normal_direction),
    )?;
    let minor_direction = unit(normal_direction.cross(major_direction))?;
    if through_apex {
        return match conic {
            Sign::Positive => Ok(IntersectionSet::Points(vec![apex])),
            Sign::Zero => Ok(IntersectionSet::Curves(vec![line_from(
                apex,
                major_direction,
                precision,
            )?])),
            Sign::Negative => {
                let scale = (-coefficient).sqrt().value;
                Ok(IntersectionSet::Curves(vec![
                    line_from(apex, major_direction + minor_direction * scale, precision)?,
                    line_from(apex, major_direction - minor_direction * scale, precision)?,
                ]))
            }
        };
    }
    if conic != Sign::Positive {
        return Err(Error::unimplemented(
            "nondegenerate oblique cone section requires an unbounded parabola or hyperbola carrier",
        ));
    }
    let displacement = (Scalar::new(1.0) + slope.sqr()) * q * height * t / coefficient;
    let center = apex + normal_direction * height.value + major_direction * displacement.value;
    let numerator = (slope * height).value.abs();
    let major = numerator / coefficient.value;
    let minor = numerator / coefficient.sqrt().value;
    // Near a perpendicular section, independent rounding can reverse equal
    // semi-axis labels. Relabeling both axes retains the same complete conic.
    let (major, minor, reference) =
        if (Scalar::new(major) - Scalar::new(minor)).sign(precision.linear())? == Sign::Negative {
            (minor, major, minor_direction)
        } else {
            (major, minor, major_direction)
        };
    let frame = Frame::new(center, normal_direction, reference, precision)?;
    Ok(IntersectionSet::Curves(vec![Curve::Ellipse(Ellipse::new(
        frame, major, minor, precision,
    )?)]))
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_4, TAU};

    use gensolid_math::{Point3, Vec3};
    use proptest::prelude::*;

    use super::*;
    use crate::{GeomStore, Surface, intersect_surfaces};

    fn plane(frame: Frame, height: f64, tilt: f64) -> Plane {
        let precision = Precision::DEFAULT;
        Plane::new(
            Frame::new(
                frame.origin() + frame.z() * height,
                frame.z() + frame.x() * tilt,
                frame.y(),
                &precision,
            )
            .unwrap(),
        )
    }

    fn check(hit: &IntersectionSet, cone: &Cone, plane: &Plane) {
        let precision = Precision::DEFAULT;
        let mut store = GeomStore::new();
        let ids = [
            store.add_surface(Surface::Cone(*cone)).unwrap(),
            store.add_surface(Surface::Plane(*plane)).unwrap(),
        ];
        let points = match hit {
            IntersectionSet::Points(points) => points.clone(),
            IntersectionSet::Curves(curves) => curves
                .iter()
                .flat_map(|curve| {
                    let bounded = curve.range().is_bounded();
                    (0..32).map(move |i| {
                        curve
                            .eval(if bounded {
                                TAU * i as f64 / 32.0
                            } else {
                                i as f64 - 16.0
                            })
                            .point
                    })
                })
                .collect(),
            _ => vec![],
        };
        for point in points {
            for id in ids {
                let residual = store.distance_to_surface(id, point, &precision).unwrap();
                assert!(
                    residual <= precision.linear(),
                    "residual={residual}, point={point:?}, hit={hit:?}"
                );
            }
        }
    }

    #[test]
    fn both_nappes_have_complete_closed_ellipses() {
        let precision = Precision::DEFAULT;
        let cone = Cone::new(Frame::WORLD, 0.0, FRAC_PI_4, &precision).unwrap();
        for height in [-2.0, 2.0] {
            let plane = plane(Frame::WORLD, height, 0.1);
            let hit = intersect(&plane, &cone, &precision).unwrap();
            assert!(
                matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Ellipse(_)]))
            );
            check(&hit, &cone, &plane);
        }
    }

    #[test]
    fn apex_sections_have_one_point_one_generator_or_two_generators() {
        let precision = Precision::DEFAULT;
        let cone = Cone::new(Frame::WORLD, 2.0, FRAC_PI_4, &precision).unwrap();
        let apex_height = -2.0 / FRAC_PI_4.tan();
        for (tilt, branches, point) in [(0.1, 1, true), (1.0, 1, false), (2.0, 2, false)] {
            let plane = plane(Frame::WORLD, apex_height, tilt);
            let hit = intersect(&plane, &cone, &precision).unwrap();
            assert_eq!(hit.branch_count(), branches);
            assert_eq!(matches!(&hit, IntersectionSet::Points(_)), point);
            check(&hit, &cone, &plane);
        }
    }

    #[test]
    fn missing_unbounded_conics_never_return_partial_branches() {
        let precision = Precision::DEFAULT;
        let cone = Cone::new(Frame::WORLD, 0.0, FRAC_PI_4, &precision).unwrap();
        for tilt in [1.0, 2.0] {
            let plane = plane(Frame::WORLD, 2.0, tilt);
            assert_eq!(
                intersect(&plane, &cone, &precision).unwrap_err().code(),
                gensolid_core::ErrorCode::Unimplemented
            );
        }
    }

    #[test]
    fn nonzero_apex_offsets_retain_their_ellipse_center() {
        let precision = Precision::DEFAULT;
        let cone = Cone::new(Frame::WORLD, 1.0, FRAC_PI_4, &precision).unwrap();
        let plane = plane(Frame::WORLD, 0.0, 0.25);
        let hit = intersect(&plane, &cone, &precision).unwrap();
        assert!(
            matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Ellipse(_)]))
        );
        check(&hit, &cone, &plane);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]
        #[test]
        fn random_frames_scales_and_nappes_have_complete_elliptic_sections(
            origin in prop::array::uniform3(-25.0f64..25.0),
            axis in prop::array::uniform3(-1.0f64..1.0),
            scale in 0.1f64..10.0,
            tilt in 0.05f64..0.7,
            negative_nappe in any::<bool>(),
        ) {
            let precision = Precision::DEFAULT;
            let Some(frame) = Frame::from_axis(Point3::new(origin[0], origin[1], origin[2]), Vec3::new(axis[0], axis[1], axis[2]), &precision).ok() else { return Ok(()); };
            let cone = Cone::new(frame, 0.0, FRAC_PI_4, &precision).unwrap();
            let plane = plane(frame, if negative_nappe { -2.0 * scale } else { 2.0 * scale }, tilt);
            let hit = intersect_surfaces(&Surface::Cone(cone), &Surface::Plane(plane), &precision).unwrap();
            prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Ellipse(_)])));
            check(&hit, &cone, &plane);
        }

        #[test]
        fn near_apex_ellipses_do_not_disappear(scale in 16.0f64..1000.0, negative_nappe in any::<bool>()) {
            let precision = Precision::DEFAULT;
            let cone = Cone::new(Frame::WORLD, 0.0, FRAC_PI_4, &precision).unwrap();
            let height = precision.linear() * scale * if negative_nappe { -1.0 } else { 1.0 };
            let plane = plane(Frame::WORLD, height, 0.1);
            let hit = intersect(&plane, &cone, &precision).unwrap();
            prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Ellipse(_)])));
            check(&hit, &cone, &plane);
        }

        #[test]
        fn near_generator_angles_keep_ellipse_and_hyperbola_classes_distinct(factor in 2.0f64..32.0) {
            let precision = Precision::DEFAULT;
            let cone = Cone::new(Frame::WORLD, 0.0, FRAC_PI_4, &precision).unwrap();
            let critical_tilt = 1.0 / FRAC_PI_4.tan();
            for (tilt, ellipse) in [(critical_tilt - factor * precision.angular(), true), (critical_tilt + factor * precision.angular(), false)] {
                let plane = plane(Frame::WORLD, precision.angular() * 64.0, tilt);
                let result = intersect(&plane, &cone, &precision);
                if ellipse {
                    let hit = result.unwrap();
                    prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Ellipse(_)])));
                    check(&hit, &cone, &plane);
                } else { prop_assert_eq!(result.unwrap_err().code(), gensolid_core::ErrorCode::Unimplemented); }
            }
        }
    }
}
