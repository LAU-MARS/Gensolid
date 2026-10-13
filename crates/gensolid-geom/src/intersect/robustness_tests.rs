use proptest::prelude::*;

use super::*;

fn check_on_surfaces(hit: &IntersectionSet, surfaces: &[Surface], precision: &Precision) {
    let mut store = crate::GeomStore::new();
    let ids: Vec<_> = surfaces
        .iter()
        .map(|surface| store.add_surface(surface.clone()).unwrap())
        .collect();
    for point in samples(hit, precision).unwrap() {
        assert!(point.is_finite());
        for &id in &ids {
            let error = store.distance_to_surface(id, point, precision).unwrap();
            assert!(
                error <= precision.linear(),
                "point={point:?}, error={error}, hit={hit:?}"
            );
        }
    }
}

#[test]
fn overlapping_surfaces_have_finite_witnesses_from_both_parameterizations() {
    let p = Precision::DEFAULT;
    let pairs = [
        (
            Surface::Plane(Plane::new(Frame::WORLD)),
            Surface::Plane(Plane::new(
                Frame::new(Point3::new(2.0, 3.0, 0.0), Vec3::Z, Vec3::Y, &p).unwrap(),
            )),
        ),
        (
            Surface::Sphere(Sphere::new(Frame::WORLD, 2.0, &p).unwrap()),
            Surface::Sphere(
                Sphere::new(
                    Frame::from_axis(Point3::ORIGIN, Vec3::X, &p).unwrap(),
                    2.0,
                    &p,
                )
                .unwrap(),
            ),
        ),
        (
            Surface::Cylinder(Cylinder::new(Frame::WORLD, 1.0, &p).unwrap()),
            Surface::Cylinder(
                Cylinder::new(
                    Frame::new(Point3::new(0.0, 0.0, 3.0), -Vec3::Z, Vec3::Y, &p).unwrap(),
                    1.0,
                    &p,
                )
                .unwrap(),
            ),
        ),
    ];
    let mut cases = Vec::new();
    for (a, b) in pairs {
        assert!(matches!(
            intersect_surfaces(&a, &b, &p).unwrap(),
            IntersectionSet::Overlap
        ));
        let measured = score_surfaces(&a, &b, 1, &p).unwrap();
        assert!(
            measured.samples.len() >= 18,
            "both inputs need finite witnesses"
        );
        assert!(measured.samples.iter().all(|point| point.is_finite()));
        assert!(pair_error(&measured, &p).unwrap() <= p.linear());
        cases.push(IntersectionCase::Surfaces {
            a,
            b,
            expected_branches: 1,
        });
    }
    let score = intersection_score(&cases, &p);
    assert_eq!(score.matched, 3, "{score:?}");
    assert_eq!(score.failed, 0);
    assert_eq!(score.unimplemented, 0);
}

#[test]
fn unresolved_overlap_offset_still_records_its_nonzero_residual() {
    let p = Precision::DEFAULT;
    let a = Surface::Plane(Plane::new(Frame::WORLD));
    let b = Surface::Plane(Plane::new(
        Frame::new(
            Point3::new(0.0, 0.0, f64::from_bits(1)),
            Vec3::Z,
            Vec3::X,
            &p,
        )
        .unwrap(),
    ));
    // The subnormal separation's interval contains zero, so the existing explicit
    // linear fallback classifies overlap. Its measured residual must remain visible.
    assert!(matches!(
        intersect_surfaces(&a, &b, &p).unwrap(),
        IntersectionSet::Overlap
    ));
    let score = intersection_score(
        &[IntersectionCase::Surfaces {
            a,
            b,
            expected_branches: 1,
        }],
        &p,
    );
    assert_eq!(score.matched, 1, "{score:?}");
    assert!(score.max_point_error > 0.0, "{score:?}");
    assert!(score.max_point_error <= p.linear());
}

#[test]
fn overlap_requires_original_geometry_and_empty_requires_no_witnesses() {
    let p = Precision::DEFAULT;
    assert_eq!(
        samples(&IntersectionSet::Overlap, &p).unwrap_err().code(),
        ErrorCode::Unimplemented
    );
    assert!(samples(&IntersectionSet::Empty, &p).unwrap().is_empty());
    let curve = Curve::Line(Line::new(Point3::ORIGIN, Vec3::X, &p).unwrap());
    let witnesses = curve_hit_samples(&IntersectionSet::Overlap, &curve, &p).unwrap();
    assert_eq!(witnesses.len(), 3);
    let mut store = crate::GeomStore::new();
    let plane = store
        .add_surface(Surface::Plane(Plane::new(Frame::WORLD)))
        .unwrap();
    for point in witnesses {
        assert!(store.distance_to_surface(plane, point, &p).unwrap() <= p.linear());
    }
    let profile = store.add_curve(curve).unwrap();
    let swept = Surface::Swept(crate::Swept::new(profile, Vec3::Y, &p).unwrap());
    assert_eq!(
        surface_witnesses(&swept, &p).unwrap_err().code(),
        ErrorCode::Unimplemented
    );
    let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, 1.0, &p).unwrap());
    let far_plane = Surface::Plane(Plane::new(
        Frame::new(Point3::new(0.0, 0.0, 2.0), Vec3::Z, Vec3::X, &p).unwrap(),
    ));
    let empty = score_surfaces(&sphere, &far_plane, 0, &p).unwrap();
    assert!(empty.samples.is_empty());
    assert_eq!(empty.found, 0);
}

#[test]
fn small_center_cuts_keep_both_roots_and_equatorial_circle() {
    let precision = Precision::DEFAULT;
    let radius = 1000.0 * precision.linear();
    let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, radius, &precision).unwrap());
    let plane = Surface::Plane(Plane::new(Frame::WORLD));
    let line = Curve::Line(Line::new(Point3::ORIGIN, Vec3::X, &precision).unwrap());
    let hit = intersect_curve_surface(&line, &sphere, &precision).unwrap();
    assert!(matches!(&hit, IntersectionSet::Points(points) if points.len() == 2));
    check_on_surfaces(&hit, std::slice::from_ref(&sphere), &precision);
    let hit = intersect_surfaces(&plane, &sphere, &precision).unwrap();
    assert!(
        matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Circle(_)]))
    );
    check_on_surfaces(&hit, &[plane, sphere], &precision);
}

#[test]
fn tangent_disjoint_and_overlapping_pairs_have_distinct_results() {
    let p = Precision::DEFAULT;
    let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, 1.0, &p).unwrap());
    for (height, branches) in [
        (1.0 - 4.0 * p.linear(), 2),
        (1.0, 1),
        (1.0 + 4.0 * p.linear(), 0),
    ] {
        let line = Curve::Line(Line::new(Point3::new(0.0, height, 0.0), Vec3::X, &p).unwrap());
        let hit = intersect_curve_surface(&line, &sphere, &p).unwrap();
        assert_eq!(hit.branch_count(), branches);
        check_on_surfaces(&hit, std::slice::from_ref(&sphere), &p);
    }
    assert_eq!(
        intersect_surfaces(&sphere, &sphere, &p).unwrap(),
        IntersectionSet::Overlap
    );
    let smaller = Surface::Sphere(Sphere::new(Frame::WORLD, 0.5, &p).unwrap());
    assert_eq!(
        intersect_surfaces(&sphere, &smaller, &p).unwrap(),
        IntersectionSet::Empty
    );
    let cylinder = Surface::Cylinder(Cylinder::new(Frame::WORLD, 1.0, &p).unwrap());
    assert_eq!(
        intersect_surfaces(&cylinder, &cylinder, &p).unwrap(),
        IntersectionSet::Overlap
    );
    let hit = intersect_surfaces(&sphere, &cylinder, &p).unwrap();
    assert_eq!(hit.branch_count(), 1);
    check_on_surfaces(&hit, &[sphere, cylinder], &p);
}

#[test]
fn double_cone_includes_negative_nappe_and_generator_overlap() {
    let p = Precision::DEFAULT;
    let angle = std::f64::consts::FRAC_PI_4;
    let cone = Surface::Cone(Cone::new(Frame::WORLD, 0.0, angle, &p).unwrap());
    let line =
        Curve::Line(Line::new(Point3::ORIGIN, Vec3::new(angle.tan(), 0.0, 1.0), &p).unwrap());
    let hit = intersect_curve_surface(&line, &cone, &p).unwrap();
    assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == 1));
    check_on_surfaces(&hit, std::slice::from_ref(&cone), &p);
    let frame = Frame::new(Point3::new(0.0, 0.0, -2.0), Vec3::Z, Vec3::X, &p).unwrap();
    let plane = Surface::Plane(Plane::new(frame));
    let hit = intersect_surfaces(&plane, &cone, &p).unwrap();
    assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == 1));
    check_on_surfaces(&hit, &[plane, cone], &p);
}

#[test]
fn non_special_cylinder_pairs_are_explicitly_unimplemented() {
    let p = Precision::DEFAULT;
    let cylinder = Surface::Cylinder(Cylinder::new(Frame::WORLD, 1.0, &p).unwrap());
    let tilted = Surface::Cylinder(
        Cylinder::new(
            Frame::from_axis(Point3::ORIGIN, Vec3::X, &p).unwrap(),
            1.0,
            &p,
        )
        .unwrap(),
    );
    let off_axis = Surface::Sphere(
        Sphere::new(
            Frame::from_axis(Point3::new(1.0, 0.0, 0.0), Vec3::Z, &p).unwrap(),
            2.0,
            &p,
        )
        .unwrap(),
    );
    for other in [tilted, off_axis] {
        assert_eq!(
            intersect_surfaces(&cylinder, &other, &p)
                .unwrap_err()
                .code(),
            ErrorCode::Unimplemented
        );
    }
}

#[test]
#[allow(clippy::float_cmp)] // Empty workloads and integer count ratios are exact.
fn scoring_keeps_unimplemented_cases_and_empty_workloads_visible() {
    let p = Precision::DEFAULT;
    let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, 2.0, &p).unwrap());
    let plane = Surface::Plane(Plane::new(Frame::WORLD));
    let spline = Curve::BCurve(
        crate::BCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![Point3::ORIGIN, Point3::new(1.0, 0.0, 0.0)],
            None,
        )
        .unwrap(),
    );
    let cases = vec![
        IntersectionCase::Surfaces {
            a: plane,
            b: sphere.clone(),
            expected_branches: 1,
        },
        IntersectionCase::CurveSurface {
            curve: spline,
            surface: sphere,
            expected_branches: 0,
        },
    ];
    let score = intersection_score(&cases, &p);
    assert_eq!(score.pairs, 2);
    assert_eq!(score.matched, 1);
    assert_eq!(score.unimplemented, 1);
    assert_eq!(score.failed, 0);
    assert_eq!(score.completeness(), 0.5);
    assert_eq!(intersection_score(&[], &p).completeness(), 0.0);
    assert_eq!(
        analytic_intersection_score(&[], &p).unwrap().completeness(),
        0.0
    );
}

#[test]
fn circle_cylinder_and_ellipse_sphere_classify_degenerate_branches() {
    let p = Precision::DEFAULT;
    let cylinder = Surface::Cylinder(Cylinder::new(Frame::WORLD, 1.0, &p).unwrap());
    let circle = Curve::Circle(Circle::new(Frame::WORLD, 1.0, &p).unwrap());
    assert!(matches!(
        intersect_curve_surface(&circle, &cylinder, &p).unwrap(),
        IntersectionSet::Curves(_)
    ));
    for (distance, branches) in [
        (2.0 - 4.0 * p.linear(), 2),
        (2.0, 1),
        (2.0 + 4.0 * p.linear(), 0),
    ] {
        let circle = Curve::Circle(
            Circle::new(
                Frame::new(Point3::new(distance, 0.0, 3.0), -Vec3::Z, Vec3::X, &p).unwrap(),
                1.0,
                &p,
            )
            .unwrap(),
        );
        let hit = intersect_curve_surface(&circle, &cylinder, &p).unwrap();
        assert_eq!(hit.branch_count(), branches);
        check_on_surfaces(&hit, std::slice::from_ref(&cylinder), &p);
    }
    let ellipse = Curve::Ellipse(Ellipse::new(Frame::WORLD, 2.0, 1.0, &p).unwrap());
    for (radius, branches) in [(0.5, 0), (1.0, 2), (1.5, 4), (2.0, 2), (3.0, 0)] {
        let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, radius, &p).unwrap());
        let hit = intersect_curve_surface(&ellipse, &sphere, &p).unwrap();
        assert_eq!(hit.branch_count(), branches);
        check_on_surfaces(&hit, std::slice::from_ref(&sphere), &p);
    }
    let ellipse_circle = Curve::Ellipse(Ellipse::new(Frame::WORLD, 1.0, 1.0, &p).unwrap());
    let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, 1.0, &p).unwrap());
    assert!(matches!(
        intersect_curve_surface(&ellipse_circle, &sphere, &p).unwrap(),
        IntersectionSet::Curves(_)
    ));
    let tilted_circle = Curve::Circle(
        Circle::new(
            Frame::from_axis(Point3::ORIGIN, Vec3::X, &p).unwrap(),
            1.0,
            &p,
        )
        .unwrap(),
    );
    assert_eq!(
        intersect_curve_surface(&tilted_circle, &cylinder, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let off_axis_sphere = Surface::Sphere(
        Sphere::new(
            Frame::from_axis(Point3::new(1.0, 0.0, 0.0), Vec3::Z, &p).unwrap(),
            1.0,
            &p,
        )
        .unwrap(),
    );
    assert_eq!(
        intersect_curve_surface(&ellipse, &off_axis_sphere, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
}

#[test]
fn linear_rational_spline_plane_handles_vertices_overlap_and_discontinuity() {
    let p = Precision::DEFAULT;
    let plane = Surface::Plane(Plane::new(Frame::WORLD));
    let make = |knots, points| Curve::BCurve(crate::BCurve::new(1, knots, points, None).unwrap());
    let vertex = make(
        vec![0.0, 0.0, 0.5, 1.0, 1.0],
        vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 1.0),
        ],
    );
    let hit = intersect_curve_surface(&vertex, &plane, &p).unwrap();
    assert!(matches!(hit, IntersectionSet::Points(points) if points.len() == 1));
    let coplanar = make(
        vec![0.0, 0.0, 0.5, 1.0, 1.0],
        vec![
            Point3::ORIGIN,
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
    );
    assert!(matches!(
        intersect_curve_surface(&coplanar, &plane, &p).unwrap(),
        IntersectionSet::Curves(_)
    ));
    let overlap_score = intersection_score(
        &[IntersectionCase::CurveSurface {
            curve: coplanar,
            surface: plane.clone(),
            expected_branches: 1,
        }],
        &p,
    );
    assert_eq!(overlap_score.matched, 1, "{overlap_score:?}");
    let partial = make(
        vec![0.0, 0.0, 0.5, 1.0, 1.0],
        vec![
            Point3::ORIGIN,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 1.0),
        ],
    );
    assert_eq!(
        intersect_curve_surface(&partial, &plane, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let discontinuous = make(
        vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0],
        vec![
            Point3::new(0.0, 0.0, -1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, -1.0),
            Point3::new(3.0, 0.0, 1.0),
        ],
    );
    let hit = intersect_curve_surface(&discontinuous, &plane, &p).unwrap();
    assert!(matches!(hit, IntersectionSet::Points(points) if points.len() == 2));
}

fn frame_strategy() -> impl Strategy<Value = Frame> {
    (
        -100.0f64..100.0,
        -100.0f64..100.0,
        -100.0f64..100.0,
        -1.0f64..1.0,
        -1.0f64..1.0,
        0.1f64..1.0,
    )
        .prop_map(|(x, y, z, ax, ay, az)| {
            Frame::from_axis(
                Point3::new(x, y, z),
                Vec3::new(ax, ay, az),
                &Precision::DEFAULT,
            )
            .unwrap()
        })
}

#[test]
fn seeded_frames_preserve_orthogonal_cylinder_section_classification() {
    let p = Precision::DEFAULT;
    let mut seed = 20261013_u64;
    let mut random = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1_u64 << 53) as f64
    };
    for _ in 0..64 {
        let origin = Point3::new(
            (2.0 * random() - 1.0) * 25.0,
            (2.0 * random() - 1.0) * 25.0,
            (2.0 * random() - 1.0) * 25.0,
        );
        let axis = Vec3::new(
            2.0 * random() - 1.0,
            2.0 * random() - 1.0,
            2.0 * random() - 1.0,
        );
        let frame = Frame::from_axis(origin, axis, &p).unwrap();
        let radius = 0.1 + 10.0 * random();
        let cylinder = Surface::Cylinder(Cylinder::new(frame, radius * 0.5, &p).unwrap());
        let side = Surface::Plane(Plane::new(
            Frame::new(
                origin + frame.x() * (radius * 0.25),
                frame.x(),
                frame.z(),
                &p,
            )
            .unwrap(),
        ));
        let hit = intersect_surfaces(&cylinder, &side, &p).unwrap();
        assert!(
            matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Line(_), Curve::Line(_)])),
            "{hit:?}"
        );
        check_on_surfaces(&hit, &[cylinder, side], &p);
    }
}

#[test]
fn coaxial_sphere_torus_distinguishes_both_tangencies_and_disjoint_sides() {
    let p = Precision::DEFAULT;
    let torus = Surface::Torus(Torus::new(Frame::WORLD, 3.0, 1.0, &p).unwrap());
    for (radius, branches) in [
        (1.0, 0),
        (2.0 - 4.0 * p.linear(), 0),
        (2.0, 1),
        (2.0 + 4.0 * p.linear(), 2),
        (3.0, 2),
        (4.0 - 4.0 * p.linear(), 2),
        (4.0, 1),
        (4.0 + 4.0 * p.linear(), 0),
        (5.0, 0),
    ] {
        let sphere = Surface::Sphere(Sphere::new(Frame::WORLD, radius, &p).unwrap());
        let hit = intersect_surfaces(&sphere, &torus, &p).unwrap();
        assert_eq!(hit.branch_count(), branches);
        check_on_surfaces(&hit, &[sphere, torus.clone()], &p);
    }
    let off_axis = Surface::Sphere(
        Sphere::new(
            Frame::from_axis(Point3::new(1.0, 0.0, 0.0), Vec3::Z, &p).unwrap(),
            2.0,
            &p,
        )
        .unwrap(),
    );
    assert_eq!(
        intersect_surfaces(&off_axis, &torus, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
}

proptest! {
    #[test]
    fn coaxial_sphere_torus_has_both_complete_latitude_circles(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        major_ratio in 1.5f64..4.0,
    ) {
        let p = Precision::DEFAULT;
        let minor = 10.0f64.powf(log_radius);
        let major = minor * major_ratio;
        let torus = Surface::Torus(Torus::new(frame, major, minor, &p).unwrap());
        let sphere = Surface::Sphere(Sphere::new(frame, major, &p).unwrap());
        let hit = intersect_surfaces(&sphere, &torus, &p).unwrap();
        prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Circle(_), Curve::Circle(_)])));
        check_on_surfaces(&hit, &[sphere, torus], &p);
    }

    #[test]
    fn meridional_torus_section_has_both_complete_tube_circles(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        major_ratio in 1.5f64..4.0,
    ) {
        let p = Precision::DEFAULT;
        let minor = 10.0f64.powf(log_radius);
        let torus = Surface::Torus(Torus::new(frame, minor * major_ratio, minor, &p).unwrap());
        let plane = Surface::Plane(Plane::new(Frame::new(frame.origin(), frame.x(), frame.z(), &p).unwrap()));
        let hit = intersect_surfaces(&plane, &torus, &p).unwrap();
        prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if matches!(curves.as_slice(), [Curve::Circle(_), Curve::Circle(_)])));
        check_on_surfaces(&hit, &[plane, torus], &p);
    }

    #[test]
    fn linear_rational_spline_plane_finds_every_span_root(
        frame in frame_strategy(),
        log_scale in -5.0f64..2.0,
        weights in prop::array::uniform5(0.2f64..5.0),
        gap_factor in 4.0f64..64.0,
    ) {
        let p = Precision::DEFAULT;
        let size = 10.0f64.powf(log_scale);
        for height in [size, gap_factor * p.linear()] {
            let points = (0..5).map(|i| frame.origin() + frame.x() * (size * i as f64) + frame.z() * (if i % 2 == 0 { -height } else { height })).collect();
            let curve = Curve::BCurve(crate::BCurve::new(1, vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0], points, Some(weights.to_vec())).unwrap());
            let plane = Surface::Plane(Plane::new(frame));
            let score = intersection_score(&[IntersectionCase::CurveSurface { curve, surface: plane, expected_branches: 4 }], &p);
            prop_assert_eq!(score.matched, 1, "{:?}", score);
            prop_assert!(score.max_point_error <= p.linear(), "{score:?}");
        }
    }

    #[test]
    fn circle_cylinder_and_ellipse_sphere_find_all_random_frame_branches(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        separation in 0.25f64..1.75,
        radius_ratio in 1.1f64..1.9,
    ) {
        let p = Precision::DEFAULT;
        let radius = 10.0f64.powf(log_radius);
        let cylinder = Surface::Cylinder(Cylinder::new(frame, radius, &p).unwrap());
        let circle = Curve::Circle(Circle::new(Frame::new(frame.origin() + frame.x() * (radius * separation), frame.z(), frame.x(), &p).unwrap(), radius, &p).unwrap());
        let ellipse = Curve::Ellipse(Ellipse::new(frame, radius * 2.0, radius, &p).unwrap());
        let sphere = Surface::Sphere(Sphere::new(frame, radius * radius_ratio, &p).unwrap());
        let cases = [
            IntersectionCase::CurveSurface { curve: circle, surface: cylinder, expected_branches: 2 },
            IntersectionCase::CurveSurface { curve: ellipse, surface: sphere, expected_branches: 4 },
        ];
        let score = intersection_score(&cases, &p);
        prop_assert_eq!(score.matched, 2, "{:?}", score);
        prop_assert!(score.max_point_error <= p.linear(), "{score:?}");
    }

    #[test]
    fn circle_sphere_and_ellipse_plane_keep_both_points(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        separation in 0.25f64..1.75,
        fraction in -0.9f64..0.9,
    ) {
        let p = Precision::DEFAULT;
        let radius = 10.0f64.powf(log_radius);
        let circle = Curve::Circle(Circle::new(frame, radius, &p).unwrap());
        let sphere = Surface::Sphere(Sphere::new(Frame::new(frame.origin() + frame.x() * (radius * separation), frame.z(), frame.x(), &p).unwrap(), radius, &p).unwrap());
        let ellipse = Curve::Ellipse(Ellipse::new(frame, radius * 2.0, radius, &p).unwrap());
        let plane = Surface::Plane(Plane::new(Frame::new(frame.origin() + frame.x() * (radius * fraction), frame.x(), frame.z(), &p).unwrap()));
        let cases = [
            IntersectionCase::CurveSurface { curve: circle, surface: sphere, expected_branches: 2 },
            IntersectionCase::CurveSurface { curve: ellipse, surface: plane, expected_branches: 2 },
        ];
        let score = intersection_score(&cases, &p);
        prop_assert_eq!(score.pairs, 2);
        prop_assert_eq!(score.matched, 2, "{:?}", score);
        prop_assert!(score.max_point_error <= p.linear(), "{score:?}");
    }

    #[test]
    fn random_scale_and_frame_sections_are_complete(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        fraction in -0.9f64..0.9,
    ) {
        let p = Precision::DEFAULT;
        let radius = 10.0f64.powf(log_radius);
        let sphere = Surface::Sphere(Sphere::new(frame, radius, &p).unwrap());
        let cylinder = Surface::Cylinder(Cylinder::new(frame, radius, &p).unwrap());
        let line = Curve::Line(Line::new(frame.origin() + frame.y() * (radius * fraction), frame.x(), &p).unwrap());
        for surface in [&sphere, &cylinder] {
            let hit = intersect_curve_surface(&line, surface, &p).unwrap();
            prop_assert!(matches!(&hit, IntersectionSet::Points(points) if points.len() == 2));
            check_on_surfaces(&hit, std::slice::from_ref(surface), &p);
            for point in samples(&hit, &p).unwrap() {
                let Curve::Line(line) = &line else { unreachable!() };
                prop_assert!(line.eval(line.closest_param(point)).point.distance(point) <= p.linear());
            }
        }
        let plane = Surface::Plane(Plane::new(Frame::new(frame.origin() + frame.z() * (radius * fraction), frame.z(), frame.x(), &p).unwrap()));
        let hit = intersect_surfaces(&plane, &sphere, &p).unwrap();
        prop_assert!(matches!(&hit, IntersectionSet::Curves(curves) if curves.len() == 1));
        check_on_surfaces(&hit, &[plane, sphere], &p);
    }

    #[test]
    fn near_tangency_resolves_both_sides_at_random_scales(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        gap_factor in 4.0f64..64.0,
    ) {
        let p = Precision::DEFAULT;
        let radius = 10.0f64.powf(log_radius);
        let sphere = Surface::Sphere(Sphere::new(frame, radius, &p).unwrap());
        for side in [-1.0, 1.0] {
            let offset = radius + side * gap_factor * p.linear();
            let line = Curve::Line(Line::new(frame.origin() + frame.y() * offset, frame.x(), &p).unwrap());
            let hit = intersect_curve_surface(&line, &sphere, &p).unwrap();
            prop_assert_eq!(hit.branch_count(), if side < 0.0 { 2 } else { 0 });
            check_on_surfaces(&hit, std::slice::from_ref(&sphere), &p);
            let plane = Surface::Plane(Plane::new(Frame::new(frame.origin() + frame.z() * offset, frame.z(), frame.x(), &p).unwrap()));
            let hit = intersect_surfaces(&plane, &sphere, &p).unwrap();
            prop_assert_eq!(hit.branch_count(), if side < 0.0 { 1 } else { 0 });
            check_on_surfaces(&hit, &[plane, sphere.clone()], &p);
        }
    }

    #[test]
    fn round_surface_pairs_keep_complete_branches(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        separation in 0.25f64..1.75,
    ) {
        let p = Precision::DEFAULT;
        let radius = 10.0f64.powf(log_radius);
        let other = Frame::new(frame.origin() + frame.x() * (radius * separation), frame.z(), frame.x(), &p).unwrap();
        let sphere_a = Surface::Sphere(Sphere::new(frame, radius, &p).unwrap());
        let sphere_b = Surface::Sphere(Sphere::new(other, radius, &p).unwrap());
        let hit = intersect_surfaces(&sphere_a, &sphere_b, &p).unwrap();
        prop_assert_eq!(hit.branch_count(), 1);
        check_on_surfaces(&hit, &[sphere_a, sphere_b], &p);
        let cylinder_a = Surface::Cylinder(Cylinder::new(frame, radius, &p).unwrap());
        let cylinder_b = Surface::Cylinder(Cylinder::new(other, radius, &p).unwrap());
        let hit = intersect_surfaces(&cylinder_a, &cylinder_b, &p).unwrap();
        prop_assert_eq!(hit.branch_count(), 2);
        check_on_surfaces(&hit, &[cylinder_a.clone(), cylinder_b], &p);
        let large_sphere = Surface::Sphere(Sphere::new(frame, radius * 2.0, &p).unwrap());
        let hit = intersect_surfaces(&large_sphere, &cylinder_a, &p).unwrap();
        prop_assert_eq!(hit.branch_count(), 2);
        check_on_surfaces(&hit, &[large_sphere, cylinder_a], &p);
    }

    #[test]
    fn double_cone_line_cuts_work_in_both_nappes(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        angle in 0.2f64..1.1,
    ) {
        let p = Precision::DEFAULT;
        let radius = 10.0f64.powf(log_radius);
        let cone = Surface::Cone(Cone::new(frame, radius, angle, &p).unwrap());
        for height in [0.0, -2.0 * radius / angle.tan()] {
            let center = frame.origin() + frame.z() * height;
            let line = Curve::Line(Line::new(center, frame.x(), &p).unwrap());
            let hit = intersect_curve_surface(&line, &cone, &p).unwrap();
            prop_assert_eq!(hit.branch_count(), 2);
            check_on_surfaces(&hit, std::slice::from_ref(&cone), &p);
        }
    }
}
