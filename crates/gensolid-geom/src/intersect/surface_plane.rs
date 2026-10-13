//! Complete plane sections of a single positive-weight rational bilinear patch.
//!
//! The homogeneous plane residual is bilinear. Between its boundary roots, a
//! section is the graph of a linear-fractional parameter function. Substituting
//! that function into the patch produces a rational quadratic BCurve. Its three
//! positive Bernstein weights follow from products of nonnegative linear domain
//! basis functions. Domain classification and branch endpoints are algebraic;
//! geometric samples are used only by tests to measure residuals.

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Point3, Sign, Vec3};

use super::IntersectionSet;
use super::certified::{Scalar, Vector};
use crate::{BCurve, BSurface, Curve, Plane};

/// The bounded patch is never replaced by an infinite intersection line.
pub(super) fn intersect(
    surface: &BSurface,
    plane: &Plane,
    precision: &Precision,
) -> Result<IntersectionSet> {
    let normal = Vector::new(plane.frame().z());
    let heights = surface
        .points()
        .iter()
        .map(|&point| Vector::between(point, plane.frame().origin()).dot(normal))
        .collect::<Vec<_>>();
    let signs = heights
        .iter()
        .map(|height| height.sign(precision.linear()))
        .collect::<Result<Vec<_>>>()?;
    // Positive rational basis functions certify the same half-space for the
    // entire patch, including higher degrees unsupported by the section solver.
    if signs.iter().all(|sign| *sign == Sign::Positive)
        || signs.iter().all(|sign| *sign == Sign::Negative)
    {
        return Ok(IntersectionSet::Empty);
    }
    if surface.degrees() != (1, 1) || surface.counts() != (2, 2) {
        return Err(Error::unimplemented(
            "plane section requires a single degree-(1,1) rational patch",
        ));
    }
    if signs.iter().all(|sign| *sign == Sign::Zero) {
        return Err(Error::unimplemented(
            "bounded coplanar surface-patch overlap requires an area representation",
        ));
    }
    let weight_scale = surface
        .weights()
        .map_or(1.0, |weights| weights.iter().copied().fold(0.0, f64::max));
    let weights: [Scalar; 4] = std::array::from_fn(|i| {
        Scalar::new(surface.weights().map_or(1.0, |weights| weights[i])) / Scalar::new(weight_scale)
    });
    for weight in weights {
        if weight.strict_sign().ok() != Some(Sign::Positive) {
            return Err(Error::unimplemented(
                "bilinear homogeneous weights are not representable at this numeric scale",
            ));
        }
    }
    let residuals: [Scalar; 4] = std::array::from_fn(|i| heights[i] * weights[i]);
    let signs: [Sign; 4] = std::array::from_fn(|i| signs[i]);
    if let Some(hit) = boundary_factors(surface, residuals, signs, precision)? {
        return Ok(hit);
    }

    // Prefer the graph direction whose denominator is uniformly nonzero. This
    // also retains nearly vertical sections of floating-frame planar patches.
    let u_slopes = [residuals[2] - residuals[0], residuals[3] - residuals[1]];
    let v_slopes = [residuals[1] - residuals[0], residuals[3] - residuals[2]];
    let stable = |slopes: [Scalar; 2]| -> Result<bool> {
        let a = slopes[0].sign(precision.linear())?;
        let b = slopes[1].sign(precision.linear())?;
        Ok(a == b && a != Sign::Zero)
    };
    let transpose = !stable(u_slopes)? && stable(v_slopes)?;
    let order = if transpose {
        [0, 2, 1, 3]
    } else {
        [0, 1, 2, 3]
    };
    let coefficients: [Scalar; 4] = std::array::from_fn(|i| residuals[order[i]]);
    let weights: [Scalar; 4] = std::array::from_fn(|i| weights[order[i]]);
    let points: [Point3; 4] = std::array::from_fn(|i| surface.points()[order[i]]);
    let mut boundaries = vec![Scalar::new(0.0), Scalar::new(1.0)];
    for ends in [
        [coefficients[0], coefficients[1]],
        [coefficients[2], coefficients[3]],
    ] {
        if let Some(root) = linear_root(ends, precision)? {
            match root.location {
                RootLocation::Interior => insert_boundary(&mut boundaries, root.parameter)?,
                RootLocation::Start | RootLocation::End => {}
            }
        }
    }
    let mut curves = Vec::new();
    let mut covered_end = [false; 2];
    for (index, interval) in boundaries.windows(2).enumerate() {
        let middle = (interval[0] + interval[1]) * Scalar::new(0.5);
        let a = linear_value([coefficients[0], coefficients[1]], middle);
        let b = linear_value([coefficients[2], coefficients[3]], middle);
        let (sa, sb) = (a.sign(precision.linear())?, b.sign(precision.linear())?);
        if matches!(
            (sa, sb),
            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
        ) {
            curves.push(section_curve(
                points,
                weights,
                coefficients,
                [interval[0], interval[1]],
                precision,
            )?);
            if index == 0 {
                covered_end[0] = true;
            }
            if index + 2 == boundaries.len() {
                covered_end[1] = true;
            }
        } else if sa == Sign::Zero || sb == Sign::Zero {
            return Err(Error::unimplemented(
                "bilinear plane-section interval or common zero is unresolved",
            ));
        }
    }
    let mut isolated = Vec::new();
    for i in 0..4 {
        if signs[order[i]] == Sign::Zero && !covered_end[i % 2] {
            let mut duplicate = false;
            for &existing in &isolated {
                let delta = Vector::between(points[i], existing);
                if delta.zero(precision.linear())? {
                    duplicate = true;
                    break;
                }
                if (delta.norm2() - Scalar::new(precision.linear()).sqr())
                    .sign(precision.linear() * precision.linear())?
                    != Sign::Positive
                {
                    return Err(Error::unimplemented(
                        "bilinear isolated contacts are below session resolution",
                    ));
                }
            }
            if !duplicate {
                isolated.push(points[i]);
            }
        }
    }
    if !curves.is_empty() && !isolated.is_empty() {
        return Err(Error::unimplemented(
            "mixed point and curve sections require a mixed intersection-set representation",
        ));
    }
    if !curves.is_empty() {
        Ok(IntersectionSet::Curves(curves))
    } else if !isolated.is_empty() {
        Ok(IntersectionSet::Points(isolated))
    } else {
        Ok(IntersectionSet::Empty)
    }
}

#[derive(Clone, Copy)]
enum RootLocation {
    Start,
    End,
    Interior,
}
struct Root {
    parameter: Scalar,
    location: RootLocation,
}

fn linear_value(ends: [Scalar; 2], parameter: Scalar) -> Scalar {
    ends[0] * (Scalar::new(1.0) - parameter) + ends[1] * parameter
}

fn linear_root(ends: [Scalar; 2], precision: &Precision) -> Result<Option<Root>> {
    let a = ends[0].sign(precision.linear())?;
    let b = ends[1].sign(precision.linear())?;
    match (a, b) {
        (Sign::Zero, Sign::Zero) => Err(Error::unimplemented(
            "identically zero boundary needs a factored bilinear section",
        )),
        (Sign::Zero, _) => Ok(Some(Root {
            parameter: Scalar::new(0.0),
            location: RootLocation::Start,
        })),
        (_, Sign::Zero) => Ok(Some(Root {
            parameter: Scalar::new(1.0),
            location: RootLocation::End,
        })),
        (Sign::Negative, Sign::Positive) | (Sign::Positive, Sign::Negative) => {
            let root = ends[0] / (ends[0] - ends[1]);
            if root.strict_sign()? != Sign::Positive
                || (Scalar::new(1.0) - root).strict_sign()? != Sign::Positive
            {
                return Err(Error::unimplemented(
                    "bilinear boundary root domain is unresolved",
                ));
            }
            Ok(Some(Root {
                parameter: root,
                location: RootLocation::Interior,
            }))
        }
        _ => Ok(None),
    }
}

fn insert_boundary(boundaries: &mut Vec<Scalar>, root: Scalar) -> Result<()> {
    for i in 1..boundaries.len() {
        let sign = (root - boundaries[i]).strict_sign().map_err(|_| {
            Error::unimplemented("bilinear boundary roots coincide or their order is unresolved")
        })?;
        if sign == Sign::Negative {
            boundaries.insert(i, root);
            return Ok(());
        }
        if sign == Sign::Zero {
            return Err(Error::unimplemented(
                "bilinear section has a common boundary zero",
            ));
        }
    }
    Err(Error::internal(
        "bilinear root escaped its certified parameter domain",
    ))
}

fn segment(points: [Point3; 2], weights: [f64; 2]) -> Result<Curve> {
    Ok(Curve::BCurve(BCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        points.to_vec(),
        Some(weights.to_vec()),
    )?))
}

/// A zero boundary factors F into one parameter basis and a linear function of
/// the other parameter. Its full edge and any transverse branch are retained.
fn boundary_factors(
    surface: &BSurface,
    residuals: [Scalar; 4],
    signs: [Sign; 4],
    precision: &Precision,
) -> Result<Option<IntersectionSet>> {
    let weights = surface
        .weights()
        .map_or([1.0; 4], |values| std::array::from_fn(|i| values[i]));
    for (edge, opposite, fixed_u) in [
        ([0, 1], [2, 3], true),
        ([2, 3], [0, 1], true),
        ([0, 2], [1, 3], false),
        ([1, 3], [0, 2], false),
    ] {
        if signs[edge[0]] != Sign::Zero || signs[edge[1]] != Sign::Zero {
            continue;
        }
        let mut curves = vec![segment(
            [surface.points()[edge[0]], surface.points()[edge[1]]],
            [weights[edge[0]], weights[edge[1]]],
        )?];
        if let Some(root) =
            linear_root([residuals[opposite[0]], residuals[opposite[1]]], precision)?
        {
            let t = root.parameter.value;
            let (ur, vr) = surface.ranges();
            let u = ur.lo + (ur.hi - ur.lo) * t;
            let v = vr.lo + (vr.hi - vr.lo) * t;
            let points = if fixed_u {
                [surface.eval(ur.lo, v).point, surface.eval(ur.hi, v).point]
            } else {
                [surface.eval(u, vr.lo).point, surface.eval(u, vr.hi).point]
            };
            let ids = if fixed_u {
                [[0, 1], [2, 3]]
            } else {
                [[0, 2], [1, 3]]
            };
            let w = ids.map(|pair| weights[pair[0]] * (1.0 - t) + weights[pair[1]] * t);
            curves.push(segment(points, w)?);
        }
        return Ok(Some(IntersectionSet::Curves(curves)));
    }
    Ok(None)
}

fn section_curve(
    points: [Point3; 4],
    weights: [Scalar; 4],
    residuals: [Scalar; 4],
    interval: [Scalar; 2],
    precision: &Precision,
) -> Result<Curve> {
    let first: [Scalar; 2] = interval.map(|t| linear_value([residuals[0], residuals[1]], t));
    let second: [Scalar; 2] = interval.map(|t| linear_value([residuals[2], residuals[3]], t));
    let middle = (interval[0] + interval[1]) * Scalar::new(0.5);
    let denominator = linear_value(
        [residuals[2] - residuals[0], residuals[3] - residuals[1]],
        middle,
    );
    let flip = match denominator.sign(precision.linear())? {
        Sign::Positive => Scalar::new(1.0),
        Sign::Negative => Scalar::new(-1.0),
        Sign::Zero => {
            return Err(Error::unimplemented(
                "bilinear section denominator is unresolved",
            ));
        }
    };
    let clean = |value: Scalar| -> Result<Scalar> {
        match value.sign(precision.linear())? {
            Sign::Positive => Ok(value),
            // Explicit residual fallback is limited to an undecidable zero sign;
            // a certified negative domain basis is never snapped into the patch.
            Sign::Zero => Ok(Scalar::new(0.0)),
            Sign::Negative => Err(Error::unimplemented(
                "bilinear endpoint leaves its parameter domain",
            )),
        }
    };
    let u = [
        [clean(second[0] * flip)?, clean(second[1] * flip)?],
        [clean(-first[0] * flip)?, clean(-first[1] * flip)?],
    ];
    let v = [
        [
            Scalar::new(1.0) - interval[0],
            Scalar::new(1.0) - interval[1],
        ],
        interval,
    ];
    let mut homogeneous_weights = [Scalar::new(0.0); 3];
    let mut offsets = [Vec3::ZERO; 3];
    for i in 0..2 {
        for j in 0..2 {
            let coefficients = [
                u[i][0] * v[j][0],
                (u[i][0] * v[j][1] + u[i][1] * v[j][0]) * Scalar::new(0.5),
                u[i][1] * v[j][1],
            ];
            for k in 0..3 {
                let weighted = coefficients[k] * weights[i * 2 + j];
                homogeneous_weights[k] = homogeneous_weights[k] + weighted;
                offsets[k] += (points[i * 2 + j] - points[0]) * weighted.value;
            }
        }
    }
    let mut poles = Vec::with_capacity(3);
    let mut output_weights = Vec::with_capacity(3);
    for k in 0..3 {
        if homogeneous_weights[k].strict_sign()? != Sign::Positive {
            return Err(Error::unimplemented(
                "bilinear section has an unrepresentable non-positive homogeneous weight",
            ));
        }
        poles.push(points[0] + offsets[k] / homogeneous_weights[k].value);
        output_weights.push(homogeneous_weights[k].value);
    }
    Ok(Curve::BCurve(BCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        poles,
        Some(output_weights),
    )?))
}

#[cfg(test)]
mod tests {
    use gensolid_math::{Frame, Point3};
    use proptest::prelude::*;

    use super::*;
    use crate::{GeomStore, Surface, intersect_surfaces};

    fn patch(points: Vec<Point3>, weights: Option<Vec<f64>>) -> BSurface {
        BSurface::new(
            1,
            1,
            vec![-2.0, -2.0, 5.0, 5.0],
            vec![1.0, 1.0, 4.0, 4.0],
            2,
            2,
            points,
            weights,
        )
        .unwrap()
    }

    fn rectangle(frame: Frame, radius: f64, weights: Option<Vec<f64>>) -> BSurface {
        patch(
            [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)]
                .into_iter()
                .map(|(x, y)| frame.to_world(Vec3::new(x * radius, y * radius, 0.0)))
                .collect(),
            weights,
        )
    }

    fn points(hit: &IntersectionSet) -> Vec<Point3> {
        match hit {
            IntersectionSet::Points(points) => points.clone(),
            IntersectionSet::Curves(curves) => curves
                .iter()
                .flat_map(|curve| (0..=32).map(move |i| curve.eval(i as f64 / 32.0).point))
                .collect(),
            _ => vec![],
        }
    }

    fn check(hit: &IntersectionSet, surface: &BSurface, plane: &Plane) {
        let precision = Precision::DEFAULT;
        let mut store = GeomStore::new();
        let ids = [
            store
                .add_surface(Surface::BSurface(surface.clone()))
                .unwrap(),
            store.add_surface(Surface::Plane(*plane)).unwrap(),
        ];
        for point in points(hit) {
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
    fn clipped_patch_section_retains_both_domain_endpoints() {
        let precision = Precision::DEFAULT;
        let surface = rectangle(Frame::WORLD, 1.0, Some(vec![1.0, 2.0, 3.0, 0.5]));
        let plane = Plane::new(Frame::new(Point3::ORIGIN, Vec3::X, Vec3::Z, &precision).unwrap());
        let hit = intersect(&surface, &plane, &precision).unwrap();
        assert_eq!(hit.branch_count(), 1);
        let IntersectionSet::Curves(curves) = &hit else {
            panic!("expected bounded section")
        };
        for (parameter, y) in [(0.0, -1.0), (1.0, 1.0)] {
            assert!(
                curves[0]
                    .eval(parameter)
                    .point
                    .distance(Point3::new(0.0, y, 0.0))
                    <= precision.linear()
            );
        }
        check(&hit, &surface, &plane);
    }

    #[test]
    fn isolated_corner_contacts_are_complete() {
        let precision = Precision::DEFAULT;
        let surface = rectangle(Frame::WORLD, 1.0, None);
        let plane = Plane::new(
            Frame::new(
                Point3::new(-1.0, -1.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::Z,
                &precision,
            )
            .unwrap(),
        );
        let hit = intersect(&surface, &plane, &precision).unwrap();
        assert!(matches!(&hit, IntersectionSet::Points(points) if points.len() == 1));
        check(&hit, &surface, &plane);
        let surface = patch(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 1.0),
                Point3::new(1.0, 0.0, 1.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            None,
        );
        let plane = Plane::new(Frame::WORLD);
        let hit = intersect(&surface, &plane, &precision).unwrap();
        assert!(matches!(&hit, IntersectionSet::Points(points) if points.len() == 2));
        check(&hit, &surface, &plane);
    }

    #[test]
    fn a_zero_boundary_and_its_transverse_branch_are_both_retained() {
        let precision = Precision::DEFAULT;
        let surface = patch(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, -1.0),
                Point3::new(1.0, 1.0, 1.0),
            ],
            Some(vec![1.0, 2.0, 3.0, 4.0]),
        );
        let plane = Plane::new(Frame::WORLD);
        let hit = intersect(&surface, &plane, &precision).unwrap();
        assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == 2));
        check(&hit, &surface, &plane);
    }

    #[test]
    fn nonplanar_saddle_section_keeps_both_rational_branches() {
        let precision = Precision::DEFAULT;
        let surface = patch(
            vec![
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(0.0, 1.0, -1.0),
                Point3::new(1.0, 0.0, -1.0),
                Point3::new(1.0, 1.0, 2.0),
            ],
            None,
        );
        let plane = Plane::new(Frame::WORLD);
        let hit = intersect(&surface, &plane, &precision).unwrap();
        assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == 2));
        check(&hit, &surface, &plane);
    }

    #[test]
    fn overlap_and_unresolved_common_zero_are_explicit() {
        let precision = Precision::DEFAULT;
        let plane = Plane::new(Frame::WORLD);
        let coplanar = rectangle(Frame::WORLD, 1.0, None);
        assert_eq!(
            intersect(&coplanar, &plane, &precision).unwrap_err().code(),
            gensolid_core::ErrorCode::Unimplemented
        );
        let saddle = patch(
            vec![
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(0.0, 1.0, -1.0),
                Point3::new(1.0, 0.0, -1.0),
                Point3::new(1.0, 1.0, 1.0),
            ],
            None,
        );
        assert_eq!(
            intersect(&saddle, &plane, &precision).unwrap_err().code(),
            gensolid_core::ErrorCode::Unimplemented
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]
        #[test]
        fn random_positive_weight_patches_have_complete_bounded_center_sections(
            origin in prop::array::uniform3(-25.0f64..25.0),
            axis in prop::array::uniform3(-1.0f64..1.0),
            radius in 0.1f64..10.0,
            weights in prop::array::uniform4(0.5f64..2.0),
        ) {
            let precision = Precision::DEFAULT;
            let Some(frame) = Frame::from_axis(Point3::new(origin[0], origin[1], origin[2]), Vec3::new(axis[0], axis[1], axis[2]), &precision).ok() else { return Ok(()); };
            let surface = rectangle(frame, radius, Some(weights.to_vec()));
            let plane = Plane::new(Frame::new(frame.origin(), frame.x(), frame.z(), &precision).unwrap());
            let hit = intersect_surfaces(&Surface::BSurface(surface.clone()), &Surface::Plane(plane), &precision).unwrap();
            prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == 1));
            check(&hit, &surface, &plane);
        }

        #[test]
        fn near_boundary_cuts_respect_the_finite_patch_domain(factor in 2.0f64..100.0) {
            let precision = Precision::DEFAULT;
            let surface = rectangle(Frame::WORLD, 1.0, None);
            for (x, empty) in [(1.0 - factor * precision.linear(), false), (1.0 + factor * precision.linear(), true)] {
                let plane = Plane::new(Frame::new(Point3::new(x, 0.0, 0.0), Vec3::X, Vec3::Z, &precision).unwrap());
                let hit = intersect(&surface, &plane, &precision).unwrap();
                if empty { prop_assert_eq!(hit, IntersectionSet::Empty); }
                else { prop_assert_eq!(hit.branch_count(), 1); check(&hit, &surface, &plane); }
            }
        }
    }
}
