//! Exchange regressions against independent OCCT exports and random geometry.

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_step::{StepStatus, read_solids, validity_report, write_block_step};
use proptest::prelude::*;

#[test]
fn occt_exports_import_with_no_skipped_checks() {
    for (name, text, expected_faces) in [
        ("block", include_str!("fixtures/block.step"), 6),
        ("cylinder", include_str!("fixtures/cylinder.step"), 3),
        ("cone", include_str!("fixtures/cone.step"), 3),
        ("torus", include_str!("fixtures/torus.step"), 1),
        (
            "block_with_hole",
            include_str!("fixtures/block_with_hole.step"),
            10,
        ),
    ] {
        let model =
            read_solids(text, Precision::DEFAULT).unwrap_or_else(|error| panic!("{name}: {error}"));
        let bodies: Vec<_> = model.bodies().collect();
        assert_eq!(bodies.len(), 1, "{name}");
        assert_eq!(
            model.body_faces(bodies[0]).unwrap().len(),
            expected_faces,
            "{name}"
        );
        let checked = check_body(&model, bodies[0], &CheckOptions::default()).unwrap();
        assert!(checked.is_valid(), "{name}: {checked:?}");
        let report = validity_report(text, Precision::DEFAULT);
        assert_eq!(report.status, StepStatus::Valid, "{name}: {report:?}");
        assert_eq!(report.checked_bodies, 1);
    }
}

#[test]
fn collapsed_spherical_vertex_loop_is_explicitly_unsupported() {
    let report = validity_report(include_str!("fixtures/sphere.step"), Precision::DEFAULT);
    assert_eq!(report.status, StepStatus::Unsupported, "{report:?}");
    assert_eq!(report.checked_bodies, 0);
}

#[test]
fn malformed_and_invalid_are_distinct() {
    let precision = Precision::DEFAULT;
    let source = write_block_step(2.0, 3.0, 4.0, &Frame::WORLD, precision).unwrap();
    let malformed = source.replace("#1=CARTESIAN_POINT", "#2=CARTESIAN_POINT");
    assert_eq!(
        validity_report(&malformed, precision).status,
        StepStatus::Parse
    );
    let marker = "#1=CARTESIAN_POINT('',(";
    let start = source.find(marker).unwrap() + marker.len();
    let end = start + source[start..].find(")").unwrap();
    let invalid = format!("{}0.5,0.,0.{}", &source[..start], &source[end..]);
    assert_eq!(
        validity_report(&invalid, precision).status,
        StepStatus::Invalid
    );
    assert_eq!(
        read_solids(&invalid, precision).err().unwrap().code(),
        ErrorCode::InvalidArgument
    );
}

#[test]
fn quoted_section_names_and_comments_are_not_delimiters() {
    let precision = Precision::DEFAULT;
    let source = write_block_step(2.0, 3.0, 4.0, &Frame::WORLD, precision).unwrap();
    let source = source
        .replace("DATA;", "/* DATA; ENDSEC; */ DATA ;")
        .replace("VERTEX_POINT('',", "VERTEX_POINT('ENDSEC; 数据',");
    assert!(read_solids(&source, precision).is_ok());
}

proptest! {
    #[test]
    fn random_scaled_and_rotated_step_round_trip(
        exponent in -4i32..=2,
        factors in (1.0f64..3.0, 1.0f64..3.0, 1.0f64..3.0),
        translation in (-5.0f64..5.0, -5.0f64..5.0, -5.0f64..5.0),
        axis in (-1.0f64..1.0, -1.0f64..1.0, -1.0f64..1.0),
    ) {
        let precision = Precision::DEFAULT;
        let scale = 10.0f64.powi(exponent);
        let axis = Vec3::new(axis.0, axis.1, axis.2);
        prop_assume!(axis.length() > precision.linear());
        let frame = Frame::from_axis(Point3::new(translation.0, translation.1, translation.2), axis, &precision).unwrap();
        let sizes = (scale * factors.0, scale * factors.1, scale * factors.2);
        let text = write_block_step(sizes.0, sizes.1, sizes.2, &frame, precision).unwrap();
        let model = read_solids(&text, precision).unwrap();
        let body = model.bodies().next().unwrap();
        prop_assert!(check_body(&model, body, &CheckOptions::default()).unwrap().is_valid());
        let volume = gensolid_topo::body_volume(&model, body).unwrap();
        let area = 2.0 * (sizes.0 * sizes.1 + sizes.1 * sizes.2 + sizes.2 * sizes.0);
        prop_assert!((volume - sizes.0 * sizes.1 * sizes.2).abs() <= precision.linear() * area);
    }

    #[test]
    fn near_degenerate_blocks_survive_file_round_trip(multiplier in 2.0f64..20.0, theta in -1.0f64..1.0) {
        let precision = Precision::DEFAULT;
        let size = precision.linear() * multiplier;
        let frame = Frame::from_axis(Point3::ORIGIN, Vec3::new(theta, 1.0, 1.0), &precision).unwrap();
        let text = write_block_step(size, size, size, &frame, precision).unwrap();
        let model = read_solids(&text, precision).unwrap();
        let body = model.bodies().next().unwrap();
        prop_assert!(check_body(&model, body, &CheckOptions::default()).unwrap().is_valid());
    }

    #[test]
    fn arbitrary_unicode_returns_a_report(text in ".{0,512}") {
        let report = validity_report(&text, Precision::DEFAULT);
        prop_assert_eq!(report.status, StepStatus::Parse);
    }
}
