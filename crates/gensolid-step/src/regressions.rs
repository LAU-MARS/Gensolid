//! Import regressions for physical-file structure, topology and rational surfaces.

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

use crate::{StepStatus, read_solids, validity_report, write_block_step};

fn block() -> String {
    write_block_step(2.0, 3.0, 4.0, &Frame::WORLD, Precision::DEFAULT).unwrap()
}

fn document(data: &str) -> String {
    let source = block();
    let start = source.find("DATA;").unwrap();
    format!(
        "{}DATA;\n{data}\nENDSEC;\nEND-ISO-10303-21;",
        &source[..start]
    )
}

#[test]
fn header_and_all_references_are_validated() {
    let source = block();
    for malformed in [
        source.replace("HEADER;\n", ""),
        source.replace(
            "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));",
            "FILE_SCHEMA(());",
        ),
        source.replace("FILE_DESCRIPTION", "FILE_NAME"),
        source.replace("END-ISO-10303-21;", "END-ISO-10303-21; EXTRA;"),
        source.replace("DATA;", "DATA; #9999=STYLED_ITEM('',(),#9998);"),
        source.replace("ORIENTED_EDGE('',*,*,", "ORIENTED_EDGE('',#1,*,"),
    ] {
        assert_eq!(
            validity_report(&malformed, Precision::DEFAULT).status,
            StepStatus::Parse
        );
    }
    let decorated = format!("/* physical file */\n{source}\n/* ENDSEC; */");
    assert_eq!(
        validity_report(&decorated, Precision::DEFAULT).status,
        StepStatus::Valid
    );
}

#[test]
fn named_data_sections_share_one_entity_namespace() {
    let source = block();
    let start = source.find("DATA;\n").unwrap();
    let end = source.rfind("ENDSEC;").unwrap();
    let records: Vec<_> = source[start + "DATA;\n".len()..end].lines().collect();
    let midpoint = records.len() / 2;
    let split = format!(
        "{}DATA('first',('CONFIG_CONTROL_DESIGN'));\n{}\nENDSEC;\nDATA('second',('CONFIG_CONTROL_DESIGN'));\n{}\nENDSEC;\nEND-ISO-10303-21;",
        &source[..start],
        records[..midpoint].join("\n"),
        records[midpoint..].join("\n")
    );
    let report = validity_report(&split, Precision::DEFAULT);
    assert_eq!(report.status, StepStatus::Valid, "{report:?}");
    assert_eq!(report.checked_bodies, 1);
    for malformed in [
        split.replace("'second'", "'first'"),
        split.replace("DATA('second',('CONFIG_CONTROL_DESIGN'))", "DATA"),
        split.replace(
            "DATA('second',('CONFIG_CONTROL_DESIGN'))",
            "DATA('second',('MISSING_SCHEMA'))",
        ),
        split.replace(
            "DATA('second',('CONFIG_CONTROL_DESIGN'));",
            "DATA('second',('CONFIG_CONTROL_DESIGN')); #1=CARTESIAN_POINT('',(0.,0.,0.));",
        ),
    ] {
        assert_eq!(
            validity_report(&malformed, Precision::DEFAULT).status,
            StepStatus::Parse
        );
    }
}

#[test]
fn external_reference_sections_are_explicitly_unimplemented() {
    let source = block().replace("DATA;", "REFERENCE; ENDSEC; DATA;");
    assert_eq!(
        read_solids(&source, Precision::DEFAULT)
            .err()
            .unwrap()
            .code(),
        ErrorCode::Unimplemented
    );
    assert_eq!(
        validity_report(&source, Precision::DEFAULT).status,
        StepStatus::Unsupported
    );
}

#[test]
fn empty_closed_shell_is_invalid_geometry() {
    let source = block();
    let empty = source
        .lines()
        .map(|line| {
            if line.contains("=CLOSED_SHELL(") {
                format!("{}=CLOSED_SHELL('',());", line.split('=').next().unwrap())
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        validity_report(&empty, Precision::DEFAULT).status,
        StepStatus::Invalid
    );
}

#[test]
fn reversed_edge_carriers_preserve_loop_orientation() {
    let source = block();
    let reversed = source
        .lines()
        .map(|line| {
            if let Some((prefix, rest)) = line.split_once("=EDGE_CURVE('',") {
                let args: Vec<_> = rest.trim_end_matches(");").split(',').collect();
                format!(
                    "{prefix}=EDGE_CURVE('',{},{},{},.F.);",
                    args[1], args[0], args[2]
                )
            } else if line.contains("=ORIENTED_EDGE(") {
                if line.contains(".T.") {
                    line.replace(".T.", ".F.")
                } else {
                    line.replace(".F.", ".T.")
                }
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let model = read_solids(&reversed, Precision::DEFAULT).unwrap();
    let body = model.bodies().next().unwrap();
    assert!(
        check_body(&model, body, &CheckOptions::default())
            .unwrap()
            .is_valid()
    );
    assert!(
        (gensolid_topo::body_volume(&model, body).unwrap() - 24.0).abs()
            <= Precision::DEFAULT.linear()
    );
}

#[test]
fn rational_curve_and_surface_preserve_nonuniform_weights() {
    let weight = 0.5f64.sqrt();
    let text = document(&format!(
        "#1=CARTESIAN_POINT('',(1.,0.,0.));
         #2=CARTESIAN_POINT('',(1.,1.,0.));
         #3=CARTESIAN_POINT('',(0.,1.,0.));
         #4=(BOUNDED_CURVE() B_SPLINE_CURVE(2,(#1,#2,#3),.CIRCULAR_ARC.,.F.,.F.)
             B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.) CURVE()
             GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,{weight},1.)) REPRESENTATION_ITEM(''));
         #5=CARTESIAN_POINT('',(0.,0.,0.));
         #6=(BOUNDED_SURFACE() B_SPLINE_SURFACE(1,1,((#5,#3),(#1,#2)),.PLANE_SURF.,.F.,.F.,.U.)
             B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)
             GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(((1.,2.),(1.,2.))) REPRESENTATION_ITEM('') SURFACE());"
    ));
    let data = crate::parse::parse_document(&text).unwrap();
    let precision = Precision::DEFAULT;
    let gensolid_geom::Curve::BCurve(curve) =
        crate::geometry::read_curve(&data, 4, &precision).unwrap()
    else {
        panic!("expected rational B-spline curve");
    };
    let point = curve.eval(0.5).point;
    assert!(point.distance(Point3::new(weight, weight, 0.0)) <= precision.linear());
    let gensolid_geom::Surface::BSurface(surface) =
        crate::geometry::read_surface(&data, 6, &precision).unwrap()
    else {
        panic!("expected rational B-spline surface");
    };
    let evaluated = surface.eval(0.5, 0.5);
    assert!(evaluated.point.distance(Point3::new(0.5, 2.0 / 3.0, 0.0)) <= precision.linear());
    assert!((evaluated.du - Vec3::X).length() <= precision.linear());
    assert!((evaluated.dv - Vec3::Y * (8.0 / 9.0)).length() <= precision.linear());
    for malformed in [
        text.replace(".CIRCULAR_ARC.", ".INVALID_FORM."),
        text.replace(".PLANE_SURF.", ".INVALID_FORM."),
        text.replace(".F.,.F.,.U.", ".F.,.F.,$"),
        text.replace("((1.,2.),(1.,2.))", "((1.,2.),(1.))"),
    ] {
        let data = crate::parse::parse_document(&malformed).unwrap();
        assert!(
            crate::geometry::read_curve(&data, 4, &precision).is_err()
                || crate::geometry::read_surface(&data, 6, &precision).is_err()
        );
    }
}

#[test]
fn cyclic_surface_curve_carriers_are_parse_errors() {
    let text = document(
        "#1=SURFACE_CURVE('',#1,(#2),.CURVE_3D.); #2=PLANE('',#3); #3=AXIS2_PLACEMENT_3D('',#4,$,$); #4=CARTESIAN_POINT('',(0.,0.,0.));",
    );
    let data = crate::parse::parse_document(&text).unwrap();
    assert_eq!(
        crate::geometry::read_curve(&data, 1, &Precision::DEFAULT)
            .err()
            .unwrap()
            .code(),
        ErrorCode::Parse
    );
}

fn rational_planar_block(source: &str, extent: f64) -> String {
    let precision = Precision::DEFAULT;
    let data = crate::parse::parse_document(source).unwrap();
    let mut next = data.keys().copied().max().unwrap() + 1;
    let mut replacements = std::collections::HashMap::new();
    let mut appended = Vec::new();
    for (id, entity) in &data {
        if entity.type_name != "PLANE" {
            continue;
        }
        let args = crate::parse::name_skipped(&entity.args).unwrap();
        let frame =
            crate::placement(&data, crate::parse::ref_of(&args[0]).unwrap(), &precision).unwrap();
        let mut point_ids = Vec::new();
        for u in [-extent, extent] {
            for v in [-extent, extent] {
                let point = frame.to_world(Vec3::new(u, v, 0.0));
                appended.push(format!(
                    "#{next}=CARTESIAN_POINT('',({},{},{}));",
                    crate::format_number(point.x),
                    crate::format_number(point.y),
                    crate::format_number(point.z)
                ));
                point_ids.push(next);
                next += 1;
            }
        }
        replacements.insert(format!("#{id}="), format!("#{id}=(BOUNDED_SURFACE() B_SPLINE_SURFACE(1,1,((#{},#{}),(#{},#{})),.PLANE_SURF.,.F.,.F.,.U.) B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(((1.,1.),(1.,1.))) REPRESENTATION_ITEM('') SURFACE());", point_ids[0], point_ids[1], point_ids[2], point_ids[3]));
    }
    source
        .lines()
        .map(|line| {
            replacements
                .iter()
                .find_map(|(prefix, replacement)| {
                    line.starts_with(prefix).then(|| replacement.clone())
                })
                .unwrap_or_else(|| line.to_owned())
        })
        .collect::<Vec<_>>()
        .join("\n")
        .replace(
            "ENDSEC;\nEND-ISO-10303-21;",
            &format!("{}\nENDSEC;\nEND-ISO-10303-21;", appended.join("\n")),
        )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn rational_bilinear_faces_import_in_random_frames(
        factors in (1.0f64..3.0, 1.0f64..3.0, 1.0f64..3.0),
        exponent in -3i32..=2,
        translation in (-5.0f64..5.0, -5.0f64..5.0, -5.0f64..5.0),
        axis in (-1.0f64..1.0, -1.0f64..1.0, -1.0f64..1.0),
        near_degenerate in any::<bool>(),
    ) {
        let precision = Precision::DEFAULT;
        let axis = Vec3::new(axis.0, axis.1, axis.2);
        prop_assume!(axis.length() > precision.linear());
        let frame = Frame::from_axis(Point3::new(translation.0, translation.1, translation.2), axis, &precision).unwrap();
        let scale = if near_degenerate { precision.linear() * 2.0 } else { 10.0f64.powi(exponent) };
        let sizes = (scale * factors.0, scale * factors.1, scale * factors.2);
        let source = write_block_step(sizes.0, sizes.1, sizes.2, &frame, precision).unwrap();
        let text = rational_planar_block(&source, sizes.0.max(sizes.1).max(sizes.2) * 2.0);
        let model = read_solids(&text, precision).unwrap();
        let body = model.bodies().next().unwrap();
        prop_assert!(check_body(&model, body, &CheckOptions::default()).unwrap().is_valid());
        prop_assert_eq!(model.body_faces(body).unwrap().len(), 6);
    }

    #[test]
    fn arbitrary_framed_unicode_never_panics(payload in ".{0,1024}") {
        let text = document(&payload);
        let report = validity_report(&text, Precision::DEFAULT);
        prop_assert!(!report.messages.is_empty());
    }
}
