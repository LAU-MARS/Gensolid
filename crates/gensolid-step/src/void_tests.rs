//! Closed shell, orientation and cavity-containment regressions.

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

use crate::{StepStatus, read_solids, validity_report, write_block_step};

#[path = "../examples/support/void_export_assertions.rs"]
mod export_contract;

fn records(source: &str) -> &str {
    let start = source.find("DATA;\n").unwrap() + "DATA;\n".len();
    let end = source.rfind("ENDSEC;").unwrap();
    &source[start..end]
}

fn shifted_ids(source: &str, shift: u32) -> String {
    let mut result = String::new();
    let mut remaining = source;
    while let Some(mark) = remaining.find('#') {
        result.push_str(&remaining[..mark + 1]);
        remaining = &remaining[mark + 1..];
        let end = remaining.bytes().take_while(u8::is_ascii_digit).count();
        let id = remaining[..end].parse::<u32>().unwrap();
        result.push_str(&(id + shift).to_string());
        remaining = &remaining[end..];
    }
    result.push_str(remaining);
    result
}

fn hollow_step(sizes: Vec3, frame: &Frame, cavity_origins: &[Point3], complex: bool) -> String {
    let precision = Precision::DEFAULT;
    let outer = write_block_step(sizes.x, sizes.y, sizes.z, frame, precision).unwrap();
    let data = crate::parse::parse_document(&outer).unwrap();
    let body = *data
        .iter()
        .find(|(_, entity)| entity.type_name == "MANIFOLD_SOLID_BREP")
        .unwrap()
        .0;
    let outer_shell =
        crate::parse::ref_of(&crate::parse::name_skipped(&data[&body].args).unwrap()[0]).unwrap();
    let mut all = records(&outer)
        .lines()
        .filter(|line| !line.contains("=MANIFOLD_SOLID_BREP("))
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let mut shift = *data.keys().max().unwrap();
    let mut voids = Vec::new();
    for origin in cavity_origins {
        let local = Frame::new(*origin, frame.z(), frame.x(), &precision).unwrap();
        let inner = write_block_step(
            sizes.x * 0.25,
            sizes.y * 0.25,
            sizes.z * 0.25,
            &local,
            precision,
        )
        .unwrap();
        let data = crate::parse::parse_document(&inner).unwrap();
        let shell = *data
            .iter()
            .find(|(_, entity)| entity.type_name == "CLOSED_SHELL")
            .unwrap()
            .0
            + shift;
        let text = records(&inner)
            .lines()
            .filter(|line| !line.contains("=MANIFOLD_SOLID_BREP("))
            .collect::<Vec<_>>()
            .join("\n");
        all.push(shifted_ids(&text, shift));
        shift += data.keys().max().unwrap();
        shift += 1;
        all.push(format!(
            "#{shift}=ORIENTED_CLOSED_SHELL('',*,#{shell},.F.);"
        ));
        voids.push(format!("#{shift}"));
    }
    let voids = voids.join(",");
    all.push(if complex {
        format!("#{body}=(BREP_WITH_VOIDS(({voids})) GEOMETRIC_REPRESENTATION_ITEM() MANIFOLD_SOLID_BREP(#{outer_shell}) REPRESENTATION_ITEM('') SOLID_MODEL());")
    } else { format!("#{body}=BREP_WITH_VOIDS('',#{outer_shell},({voids}));") });
    let header = &outer[..outer.find("DATA;").unwrap()];
    format!(
        "{header}DATA;\n{}\nENDSEC;\nEND-ISO-10303-21;",
        all.join("\n")
    )
}

#[test]
fn independent_occt_void_export_preserves_regions_and_volume() {
    let text = include_str!("../tests/fixtures/hollow_block.step");
    let precision = Precision::DEFAULT;
    let model = read_solids(text, precision).unwrap();
    let bodies: Vec<_> = model.bodies().collect();
    assert_eq!(bodies.len(), 1);
    let body = bodies[0];
    let regions = &model.body(body).unwrap().regions;
    assert_eq!(regions.len(), 3);
    assert!(!model.region(regions[0]).unwrap().solid);
    assert!(model.region(regions[1]).unwrap().solid);
    assert!(!model.region(regions[2]).unwrap().solid);
    assert_eq!(model.region(regions[1]).unwrap().shells.len(), 2);
    assert_eq!(model.body_faces(body).unwrap().len(), 12);
    assert!(
        check_body(&model, body, &CheckOptions::default())
            .unwrap()
            .is_valid()
    );
    assert!(
        (gensolid_topo::body_volume(&model, body).unwrap() - 168.0).abs()
            <= precision.linear() * 264.0
    );
    assert_eq!(validity_report(text, precision).status, StepStatus::Valid);
    // A maintainer-owned capability transition: lack of support stays explicit;
    // successful export must preserve the independently specified finite cavity.
    // The frozen executable probe separately requires success to award capability.
    match crate::write_planar_solid(&model, body) {
        Err(error) => assert_eq!(error.code(), ErrorCode::Unimplemented),
        Ok(text) => {
            let report = validity_report(&text, precision);
            assert_eq!(report.status, StepStatus::Valid);
            assert_eq!(report.decoded_bodies, 1);
            assert_eq!(report.checked_bodies, 1);
            assert_eq!(report.faults, 0);
            assert_eq!(report.unchecked, 0);
            let output = read_solids(&text, precision).unwrap();
            let body = output.bodies().next().unwrap();
            export_contract::certify(&output, body, &export_contract::occt_boxes()).unwrap();
        }
    }
}

#[test]
fn uncertified_cavity_relationships_never_become_valid() {
    let precision = Precision::DEFAULT;
    let sizes = Vec3::new(4.0, 6.0, 8.0);
    let outside = hollow_step(sizes, &Frame::WORLD, &[Point3::new(4.0, 0.0, 0.0)], false);
    assert_eq!(
        validity_report(&outside, precision).status,
        StepStatus::Invalid
    );
    let touching = hollow_step(sizes, &Frame::WORLD, &[Point3::new(1.5, 0.0, 0.0)], false);
    assert_eq!(
        read_solids(&touching, precision).err().unwrap().code(),
        ErrorCode::Unimplemented
    );
    let overlap = hollow_step(
        sizes,
        &Frame::WORLD,
        &[Point3::ORIGIN, Point3::ORIGIN],
        false,
    );
    assert_eq!(
        read_solids(&overlap, precision).err().unwrap().code(),
        ErrorCode::Unimplemented
    );
    let reversed = hollow_step(sizes, &Frame::WORLD, &[Point3::ORIGIN], false)
        .lines()
        .map(|line| {
            if line.contains("=ORIENTED_CLOSED_SHELL(") {
                line.replace(".F.", ".T.")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_ne!(
        validity_report(&reversed, precision).status,
        StepStatus::Valid
    );
}

#[test]
fn nonconvex_and_curved_void_exports_remain_unimplemented() {
    for text in [
        include_str!("../tests/fixtures/hollow_concave_outer.step"),
        include_str!("../tests/fixtures/hollow_cylinder.step"),
    ] {
        assert_eq!(
            read_solids(text, Precision::DEFAULT).err().unwrap().code(),
            ErrorCode::Unimplemented
        );
        assert_eq!(
            validity_report(text, Precision::DEFAULT).status,
            StepStatus::Unsupported
        );
        let report = validity_report(text, Precision::DEFAULT);
        assert_eq!(report.decoded_bodies, 1);
        assert_eq!(report.checked_bodies, 1);
        assert_eq!(report.faults, 0);
    }
}

#[test]
fn malformed_void_shells_are_not_silently_imported() {
    let precision = Precision::DEFAULT;
    let original = include_str!("../tests/fixtures/hollow_block.step");
    let open = original.replace("#346 = CLOSED_SHELL('',(#347,", "#346 = CLOSED_SHELL('',(");
    assert_eq!(
        validity_report(&open, precision).status,
        StepStatus::Invalid
    );
    let report = validity_report(&open, precision);
    assert_eq!(report.decoded_bodies, 1);
    assert_eq!(report.checked_bodies, 1);
    assert!(report.faults > 0);
    let duplicate = original.replace(
        "BREP_WITH_VOIDS('',#16,(#345))",
        "BREP_WITH_VOIDS('',#16,(#345,#345))",
    );
    assert_eq!(
        validity_report(&duplicate, precision).status,
        StepStatus::Parse
    );
    let empty = original.replace(
        "BREP_WITH_VOIDS('',#16,(#345))",
        "BREP_WITH_VOIDS('',#16,())",
    );
    assert_eq!(validity_report(&empty, precision).status, StepStatus::Parse);
    let nested = original.replace(
        "ORIENTED_CLOSED_SHELL('',*,#346,.F.)",
        "ORIENTED_CLOSED_SHELL('',*,#345,.F.)",
    );
    assert_eq!(
        validity_report(&nested, precision).status,
        StepStatus::Parse
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn cavities_close_to_the_containment_margin_are_checked_or_refused(
        multiplier in 0.25f64..4.0,
        theta in -1.0f64..1.0,
    ) {
        let precision = Precision::DEFAULT;
        let frame = Frame::from_axis(Point3::new(1.0, 2.0, 3.0), Vec3::new(theta, 1.0, 1.0), &precision).unwrap();
        let sizes = Vec3::new(4.0, 6.0, 8.0);
        let gap = precision.linear() * multiplier;
        let origin = frame.to_world(Vec3::new(sizes.x * 0.375 - gap, 0.0, 0.0));
        let text = hollow_step(sizes, &frame, &[origin], false);
        let report = validity_report(&text, precision);
        prop_assert!(matches!(report.status, StepStatus::Valid | StepStatus::Unsupported));
        if report.status == StepStatus::Valid {
            let model = read_solids(&text, precision).unwrap();
            let body = model.bodies().next().unwrap();
            prop_assert!(check_body(&model, body, &CheckOptions::default()).unwrap().is_valid());
        }
    }

    #[test]
    fn random_hollow_blocks_are_checked_and_subtract_cavity_volume(
        factors in (1.0f64..3.0, 1.0f64..3.0, 1.0f64..3.0),
        exponent in -3i32..=2,
        translation in (-10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0),
        axis in (-1.0f64..1.0, -1.0f64..1.0, -1.0f64..1.0),
        near_degenerate in any::<bool>(),
        complex in any::<bool>(),
        multiple in any::<bool>(),
    ) {
        let precision = Precision::DEFAULT;
        let axis = Vec3::new(axis.0, axis.1, axis.2);
        prop_assume!(axis.length() > precision.linear());
        let frame = Frame::from_axis(Point3::new(translation.0, translation.1, translation.2), axis, &precision).unwrap();
        let scale = if near_degenerate { precision.linear() * 20.0 } else { 10.0f64.powi(exponent) };
        let sizes = Vec3::new(scale * factors.0, scale * factors.1, scale * factors.2);
        let origins = if multiple {
            vec![frame.to_world(Vec3::new(-sizes.x * 0.25, 0.0, 0.0)), frame.to_world(Vec3::new(sizes.x * 0.25, 0.0, 0.0))]
        } else { vec![frame.origin()] };
        let text = hollow_step(sizes, &frame, &origins, complex);
        let model = read_solids(&text, precision).unwrap();
        let body = model.bodies().next().unwrap();
        prop_assert!(check_body(&model, body, &CheckOptions::default()).unwrap().is_valid());
        prop_assert_eq!(model.body(body).unwrap().regions.len(), 2 + origins.len());
        let expected = sizes.x * sizes.y * sizes.z * (1.0 - origins.len() as f64 / 64.0);
        let area = 2.0 * (sizes.x * sizes.y + sizes.x * sizes.z + sizes.y * sizes.z) * (1.0 + origins.len() as f64 / 16.0);
        prop_assert!((gensolid_topo::body_volume(&model, body).unwrap() - expected).abs() <= precision.linear() * area);
    }
}
