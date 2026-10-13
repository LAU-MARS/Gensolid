//! Schema framing, adapter fidelity, and kernel validity regression tests.
mod support;

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::body_volume;
use gensolid_xt::{
    TextHeader, XtStatus, read_bytes, validity_report, validity_report_bytes, write_text_header,
};
use proptest::prelude::*;

#[test]
fn text_and_neutral_binary_blocks_pass_kernel_checks() {
    for binary in [false, true] {
        let bytes = support::block([1.0, 2.0, 3.0], Frame::WORLD, binary);
        let report = validity_report_bytes(&bytes);
        assert_eq!(report.status, XtStatus::Checked, "{report:?}");
        assert_eq!(report.decoded_bodies, 1);
        let imported = read_bytes(&bytes, Precision::DEFAULT).unwrap();
        let body = imported.bodies[0];
        assert_eq!(imported.model.body_faces(body).unwrap().len(), 6);
        assert_eq!(imported.model.body_edges(body).unwrap().len(), 12);
        assert_eq!(imported.model.body_vertices(body).unwrap().len(), 8);
        assert!(
            (body_volume(&imported.model, body).unwrap() - 6.0).abs() < Precision::DEFAULT.linear()
        );
    }
}

#[test]
fn truncated_and_trailing_streams_never_claim_validity() {
    for binary in [false, true] {
        let bytes = support::block([1.0, 2.0, 3.0], Frame::WORLD, binary);
        for offset in (0..bytes.len()).step_by(31) {
            assert_ne!(
                validity_report_bytes(&bytes[..offset]).status,
                XtStatus::Checked
            );
        }
        let mut trailing = bytes;
        trailing.push(b'X');
        assert_eq!(validity_report_bytes(&trailing).status, XtStatus::Parse);
    }
}

#[test]
fn source_and_kernel_invalidity_have_invalid_reports() {
    for binary in [false, true] {
        for bytes in [
            support::block_with_invalid_geometry(binary),
            support::block_with_invalid_topology(binary),
        ] {
            let report = validity_report_bytes(&bytes);
            assert_eq!(report.status, XtStatus::Invalid, "{report:?}");
            assert!(read_bytes(&bytes, Precision::DEFAULT).is_err());
        }
    }
}

#[test]
fn unchecked_acorn_geometry_is_explicitly_unsupported() {
    for binary in [false, true] {
        let bytes = support::acorn(binary);
        let report = validity_report_bytes(&bytes);
        assert_eq!(report.status, XtStatus::UnsupportedGeometry, "{report:?}");
        assert_eq!(report.decoded_bodies, 0);
        assert!(
            matches!(read_bytes(&bytes, Precision::DEFAULT), Err(error) if error.code() == ErrorCode::Unimplemented)
        );
    }
}

#[test]
fn common_header_charset_damage_cannot_be_checked() {
    for binary in [false, true] {
        let header = write_text_header(&TextHeader::default()).unwrap();
        let mut bytes = header.as_bytes().to_vec();
        bytes.extend_from_slice(&support::block([1.0, 2.0, 3.0], Frame::WORLD, binary));
        assert_eq!(validity_report_bytes(&bytes).status, XtStatus::Checked);
        bytes[2] = b'Z';
        assert_eq!(validity_report_bytes(&bytes).status, XtStatus::Parse);
    }
}

#[test]
fn exact_schema_selection_does_not_guess_adjacent_versions() {
    let bytes = support::block([1.0, 2.0, 3.0], Frame::WORLD, false);
    let text = String::from_utf8(bytes)
        .unwrap()
        .replace("SCH_3000000_30000", "SCH_3000001_30000");
    let report = validity_report_bytes(text.as_bytes());
    assert_eq!(report.status, XtStatus::UnsupportedSchema);
    assert!(
        matches!(read_bytes(text.as_bytes(), Precision::DEFAULT), Err(error) if error.code() == ErrorCode::Unimplemented)
    );
}

#[test]
fn shared_geometry_senses_preserve_boundary_direction_and_material_side() {
    for binary in [false, true] {
        let bytes = support::block_with_negative_carrier_senses(binary);
        let imported = read_bytes(&bytes, Precision::DEFAULT).unwrap();
        let body = imported.bodies[0];
        assert!(
            check_body(&imported.model, body, &CheckOptions::default())
                .unwrap()
                .is_valid()
        );
        assert!(
            (body_volume(&imported.model, body).unwrap() - 6.0).abs() < Precision::DEFAULT.linear()
        );
    }
}

/// The public ABC smoke corpus is fetched separately, with origin and hashes in
/// its manifest. The regular synthetic suite does not require network access.
#[test]
fn downloaded_abc_parts_preserve_all_boundary_topology() {
    let corpus = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../corpus/abc-smoke");
    if !corpus.is_dir() {
        return;
    }
    for part in 1..=3 {
        let bytes = std::fs::read(corpus.join(format!("Part Studio 1 - Part {part}.x_t"))).unwrap();
        let text = std::str::from_utf8(&bytes).unwrap();
        let header = gensolid_xt::read_text(text).unwrap();
        assert_eq!(header.schema.as_deref(), Some("SCH_3001168_30100_13006"));
        assert_eq!(header.modeller_version.as_deref(), Some("3001168"));
        let report = validity_report(text);
        assert_eq!(report.status, XtStatus::Checked, "part {part}: {report:?}");
        assert_eq!(report.schema.as_deref(), Some("SCH_3001168_30100_13006"));
        let imported = read_bytes(&bytes, Precision::DEFAULT).unwrap();
        let body = imported.bodies[0];
        assert!(!imported.model.body_faces(body).unwrap().is_empty());
        assert!(!imported.model.body_edges(body).unwrap().is_empty());
        // Ring edges have no endpoint vertices; ABC Part 3 exercises that case.
        if part == 3 {
            assert!(imported.model.body_vertices(body).unwrap().is_empty());
        }
        assert!(
            check_body(&imported.model, body, &CheckOptions::default())
                .unwrap()
                .is_valid()
        );
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(40))]
    #[test]
    fn arbitrary_bytes_have_a_report(input in prop::collection::vec(any::<u8>(), 0..2048)) {
        let report = validity_report_bytes(&input);
        prop_assert_ne!(report.status, XtStatus::Checked);
    }

    #[test]
    fn imported_blocks_preserve_random_sizes_and_frames(
        dimensions in prop::array::uniform3(0.05f64..20.0),
        origin in prop::array::uniform3(-20.0f64..20.0),
        axis in prop::array::uniform3(-1.0f64..1.0),
        binary in any::<bool>(),
    ) {
        let precision = Precision::DEFAULT;
        let Ok(frame) = Frame::from_axis(Point3::new(origin[0], origin[1], origin[2]), Vec3::new(axis[0], axis[1], axis[2]), &precision) else { return Ok(()); };
        let bytes = support::block(dimensions, frame, binary);
        let imported = read_bytes(&bytes, precision).unwrap();
        let body = imported.bodies[0];
        prop_assert!(check_body(&imported.model, body, &CheckOptions::default()).unwrap().is_valid());
        let expected: f64 = dimensions.iter().product();
        prop_assert!((body_volume(&imported.model, body).unwrap() - expected).abs() <= precision.linear() * expected.max(1.0));
    }

    #[test]
    fn blocks_close_to_degenerate_resolution_have_explicit_outcomes(
        factors in prop::array::uniform3(2.0f64..1000.0),
        binary in any::<bool>(),
    ) {
        let precision = Precision::DEFAULT;
        let dimensions = factors.map(|f| f * precision.linear());
        let bytes = support::block(dimensions, Frame::WORLD, binary);
        let imported = read_bytes(&bytes, precision).unwrap();
        prop_assert!(check_body(&imported.model, imported.bodies[0], &CheckOptions::default()).unwrap().is_valid());
    }
}
