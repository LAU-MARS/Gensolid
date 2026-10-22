//! Frozen, deterministic capability probe. Unsupported is observable, never a success.
//! The RSI evaluator pins this file, its assertion module and the independent fixture.

#[path = "support/void_export_assertions.rs"]
mod contract;

use contract::{BoxSpec, certify, occt_boxes, require};
use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_geom::{Curve, Line, Plane, Surface};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_step::{
    StepStatus, read_solids, validity_report, write_block_step, write_planar_solid,
};
use gensolid_topo::{BodyId, Model, ShellId};
use std::collections::HashSet;

const IDS: [&str; 16] = [
    "occt_hollow_block",
    "one_void_world_simple",
    "one_void_transformed_complex",
    "two_void_world_simple",
    "two_void_transformed_complex",
    "one_void_scaled_small",
    "one_void_scaled_large",
    "one_void_near_precision",
    "single_block_control",
    "reject_vertex_outside",
    "reject_fin_owner",
    "reject_empty_shell",
    "reject_duplicate_shell",
    "reject_region_backpointer",
    "reject_closed_cavity_outside",
    "reject_closed_cavities_overlap",
];

fn first_body(model: &Model) -> Result<BodyId> {
    model
        .bodies()
        .next()
        .ok_or_else(|| Error::invalid_argument("missing body"))
}

fn references(text: &str) -> Result<Vec<u32>> {
    let mut result = Vec::new();
    let mut tail = text;
    while let Some(mark) = tail.find('#') {
        tail = &tail[mark + 1..];
        let end = tail.bytes().take_while(u8::is_ascii_digit).count();
        result.push(
            tail[..end]
                .parse()
                .map_err(|_| Error::parse("invalid generated entity identifier"))?,
        );
        tail = &tail[end..];
    }
    Ok(result)
}

fn shifted_ids(source: &str, shift: u32) -> Result<String> {
    let mut result = String::new();
    let mut tail = source;
    while let Some(mark) = tail.find('#') {
        result.push_str(&tail[..mark + 1]);
        tail = &tail[mark + 1..];
        let end = tail.bytes().take_while(u8::is_ascii_digit).count();
        let id: u32 = tail[..end]
            .parse()
            .map_err(|_| Error::parse("invalid generated identifier"))?;
        result.push_str(&(id + shift).to_string());
        tail = &tail[end..];
    }
    result.push_str(tail);
    Ok(result)
}

fn records(text: &str) -> Result<&str> {
    let start = text
        .find("DATA;\n")
        .ok_or_else(|| Error::parse("generated DATA missing"))?
        + "DATA;\n".len();
    let end = text
        .rfind("ENDSEC;")
        .ok_or_else(|| Error::parse("generated ENDSEC missing"))?;
    Ok(&text[start..end])
}

fn source(boxes: &[BoxSpec], complex: bool) -> Result<String> {
    let mut lines = Vec::new();
    let mut shift = 0;
    let mut outer_shell = 0;
    let mut voids = Vec::new();
    let mut header = String::new();
    for (index, spec) in boxes.iter().enumerate() {
        let text = write_block_step(
            spec.sizes.x,
            spec.sizes.y,
            spec.sizes.z,
            &spec.frame,
            Precision::DEFAULT,
        )?;
        if index == 0 {
            let at = text
                .find("DATA;")
                .ok_or_else(|| Error::parse("generated header missing"))?;
            header.push_str(&text[..at]);
        }
        let data = records(&text)?;
        let body_line = data
            .lines()
            .find(|line| line.contains("=MANIFOLD_SOLID_BREP("))
            .ok_or_else(|| Error::parse("generated body missing"))?;
        let refs = references(body_line)?;
        require(refs.len() == 2, "generated body reference layout changed")?;
        let shell = refs[1] + shift;
        let max_id = references(data)?
            .into_iter()
            .max()
            .ok_or_else(|| Error::parse("generated identifiers missing"))?;
        for line in data
            .lines()
            .filter(|line| !line.contains("=MANIFOLD_SOLID_BREP("))
        {
            lines.push(shifted_ids(line, shift)?);
        }
        shift += max_id;
        if index == 0 {
            outer_shell = shell;
        } else {
            shift += 1;
            lines.push(format!(
                "#{shift}=ORIENTED_CLOSED_SHELL('',*,#{shell},.F.);"
            ));
            voids.push(format!("#{shift}"));
        }
    }
    shift += 1;
    let voids = voids.join(",");
    lines.push(if boxes.len() == 1 {
        format!("#{shift}=MANIFOLD_SOLID_BREP('',#{outer_shell});")
    } else if complex {
        format!("#{shift}=(BREP_WITH_VOIDS(({voids})) GEOMETRIC_REPRESENTATION_ITEM() MANIFOLD_SOLID_BREP(#{outer_shell}) REPRESENTATION_ITEM('') SOLID_MODEL());")
    } else {
        format!("#{shift}=BREP_WITH_VOIDS('',#{outer_shell},({voids}));")
    });
    Ok(format!(
        "{header}DATA;\n{}\nENDSEC;\nEND-ISO-10303-21;",
        lines.join("\n")
    ))
}

fn boxes(scale: f64, transformed: bool, multiple: bool) -> Result<Vec<BoxSpec>> {
    let p = Precision::DEFAULT;
    let frame = if transformed {
        Frame::new(
            Point3::new(3.0, -2.0, 5.0),
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(2.0, -1.0, 0.0),
            &p,
        )?
    } else {
        Frame::WORLD
    };
    let sizes = Vec3::new(4.0 * scale, 6.0 * scale, 8.0 * scale);
    let mut result = vec![BoxSpec { frame, sizes }];
    for offset in if multiple {
        vec![-0.25, 0.25]
    } else {
        vec![0.0]
    } {
        let frame = Frame::new(
            frame.to_world(Vec3::new(sizes.x * offset, 0.0, 0.0)),
            frame.z(),
            frame.x(),
            &p,
        )?;
        result.push(BoxSpec {
            frame,
            sizes: sizes * 0.25,
        });
    }
    Ok(result)
}

fn verify_text(text: &str, boxes: &[BoxSpec]) -> Result<()> {
    let p = Precision::DEFAULT;
    let report = validity_report(text, p);
    require(
        report.status == StepStatus::Valid
            && report.decoded_bodies == 1
            && report.checked_bodies == 1
            && report.faults == 0
            && report.unchecked == 0,
        "export failed complete STEP validity contract",
    )?;
    let model = read_solids(text, p)?;
    certify(&model, first_body(&model)?, boxes)
}

fn positive(index: usize) -> Result<bool> {
    let (text, expected) = if index == 0 {
        (
            include_str!("../tests/fixtures/hollow_block.step").to_owned(),
            occt_boxes(),
        )
    } else {
        let scale = match index {
            5 => 0.001,
            6 => 100.0,
            7 => Precision::DEFAULT.linear() * 20.0,
            _ => 1.0,
        };
        let mut expected = boxes(scale, matches!(index, 2 | 4), matches!(index, 3 | 4))?;
        if index == 8 {
            expected.truncate(1);
        }
        (source(&expected, matches!(index, 2 | 4))?, expected)
    };
    verify_text(&text, &expected)?;
    let model = read_solids(&text, Precision::DEFAULT)?;
    match write_planar_solid(&model, first_body(&model)?) {
        Ok(exported) => {
            verify_text(&exported, &expected)?;
            Ok(true)
        }
        Err(error) if error.code() == ErrorCode::Unimplemented && index != 8 => Ok(false),
        Err(error) => Err(error),
    }
}

/// Move the entire closed boundary, including its exact line/plane carriers.
/// This deliberately preserves all checker invariants while violating relationships
/// between shells, which the immutable body checker does not certify.
fn translate_shell(model: &mut Model, shell: ShellId, delta: Vec3) -> Result<()> {
    let precision = *model.precision();
    let faces = model.shell(shell)?.faces.clone();
    let mut vertices = HashSet::new();
    let mut edges = HashSet::new();
    for face in &faces {
        for lp in &model.face(*face)?.loops {
            for fin in &model.loop_(*lp)?.fins {
                let edge_id = model.fin(*fin)?.edge;
                edges.insert(edge_id);
                let edge = model.edge(edge_id)?;
                for vertex in [edge.start, edge.end] {
                    vertices.insert(vertex.ok_or_else(|| {
                        Error::invalid_argument("translation needs bounded edges")
                    })?);
                }
            }
        }
    }
    for vertex in vertices {
        let point = model.vertex(vertex)?.point;
        model.raw().vertex_mut(vertex)?.point = point + delta;
    }
    for edge_id in edges {
        let edge = model.edge(edge_id)?;
        let start = edge
            .start
            .ok_or_else(|| Error::invalid_argument("translation start missing"))?;
        let end = edge
            .end
            .ok_or_else(|| Error::invalid_argument("translation end missing"))?;
        let line = Line::through(
            model.vertex(start)?.point,
            model.vertex(end)?.point,
            &precision,
        )?;
        let curve = model.raw().geom_mut().add_curve(Curve::Line(line))?;
        model.raw().edge_mut(edge_id)?.curve = Some(curve);
    }
    for face_id in faces {
        let surface = model
            .face(face_id)?
            .surface
            .ok_or_else(|| Error::invalid_argument("translation surface missing"))?;
        let Surface::Plane(plane) = model.geom().surface(surface)? else {
            return Err(Error::invalid_argument("translation needs planes"));
        };
        let frame = Frame::new(
            plane.frame().origin() + delta,
            plane.frame().z(),
            plane.frame().x(),
            &precision,
        )?;
        let surface = model
            .raw()
            .geom_mut()
            .add_surface(Surface::Plane(Plane::new(frame)))?;
        model.raw().face_mut(face_id)?.surface = Some(surface);
    }
    Ok(())
}

fn negative(index: usize) -> Result<bool> {
    let mut expected = if index == 15 {
        boxes(1.0, false, true)?
    } else {
        occt_boxes()
    };
    let text = if index == 15 {
        source(&expected, false)?
    } else {
        include_str!("../tests/fixtures/hollow_block.step").to_owned()
    };
    let mut model = read_solids(&text, Precision::DEFAULT)?;
    let body = first_body(&model)?;
    let regions = model.body(body)?.regions.clone();
    let shells = model.region(regions[1])?.shells.clone();
    match index {
        9 => {
            let face = model.shell(shells[1])?.faces[0];
            let lp = model.face(face)?.loops[0];
            let fin = model.loop_(lp)?.fins[0];
            let vertex = model
                .fin_vertices(fin)?
                .0
                .ok_or_else(|| Error::invalid_argument("negative fixture missing vertex"))?;
            let point = model.vertex(vertex)?.point;
            model.raw().vertex_mut(vertex)?.point = point + Vec3::new(100.0, 0.0, 0.0);
        }
        10 => {
            let faces = model.shell(shells[0])?.faces.clone();
            let lp = model.face(faces[0])?.loops[0];
            let other = model.face(faces[1])?.loops[0];
            let fin = model.loop_(lp)?.fins[0];
            model.raw().fin_mut(fin)?.owner = other;
        }
        11 => model.raw().shell_mut(shells[1])?.faces.clear(),
        12 => model.raw().region_mut(regions[1])?.shells.push(shells[1]),
        13 => {
            let shell = model.region(regions[2])?.shells[0];
            model.raw().shell_mut(shell)?.region = regions[1];
        }
        14 | 15 => {
            // Outside: cavity x range [9,11], while the outer box ends at x=2.
            // Overlap: two unit-width cavities centered at -1 and -0.5, each
            // still strictly inside the outer box but overlapping by width 0.5.
            let (which, delta) = if index == 14 {
                (1, Vec3::new(10.0, 0.0, 0.0))
            } else {
                (2, Vec3::new(-1.5, 0.0, 0.0))
            };
            translate_shell(&mut model, shells[which], delta)?;
            let spec = &mut expected[which];
            spec.frame = Frame::new(
                spec.frame.origin() + delta,
                spec.frame.z(),
                spec.frame.x(),
                model.precision(),
            )?;
            // These assertions establish the intended analytic translated shells;
            // they do not certify containment or nonintersection between them.
            certify(&model, body, &expected)?;
            let report = check_body(&model, body, &CheckOptions::default())?;
            require(
                report.faults.is_empty() && report.unchecked == 0,
                "relationship negative must remain individually checker-valid",
            )?;
        }
        _ => return Err(Error::invalid_argument("unknown negative contract case")),
    }
    match write_planar_solid(&model, body) {
        Err(error)
            if matches!(
                error.code(),
                ErrorCode::Unimplemented
                    | ErrorCode::Unsupported
                    | ErrorCode::InvalidArgument
                    | ErrorCode::Degenerate
                    | ErrorCode::NotFound
            ) =>
        {
            Ok(true)
        }
        Err(error) => Err(error),
        Ok(_) => Err(Error::invalid_argument("writer accepted a damaged model")),
    }
}

fn quoted(text: &str) -> String {
    let mut result = String::from("\"");
    for ch in text.chars() {
        match ch {
            '"' => result.push_str("\\\""),
            '\\' => result.push_str("\\\\"),
            '\n' => result.push_str("\\n"),
            '\r' => result.push_str("\\r"),
            '\t' => result.push_str("\\t"),
            ch if ch.is_control() => result.push_str(&format!("\\u{:04x}", ch as u32)),
            ch => result.push(ch),
        }
    }
    result.push('"');
    result
}

fn main() {
    let (mut passed, mut unsupported, mut failed) = (0, 0, 0);
    let mut rows = Vec::new();
    for (index, id) in IDS.iter().enumerate() {
        let result = if index < 9 {
            positive(index)
        } else {
            negative(index)
        };
        let (status, detail) = match result {
            Ok(true) => {
                passed += 1;
                ("passed", "all fixed assertions passed".to_owned())
            }
            Ok(false) => {
                unsupported += 1;
                (
                    "unsupported",
                    "writer returned exactly ErrorCode::Unimplemented".to_owned(),
                )
            }
            Err(error) => {
                failed += 1;
                ("failed", error.to_string())
            }
        };
        rows.push(format!(
            "{{\"id\":{},\"status\":{},\"detail\":{}}}",
            quoted(id),
            quoted(status),
            quoted(&detail)
        ));
    }
    println!(
        "{{\"schema_version\":1,\"contract\":\"step_planar_void_export_v1\",\"cases\":{},\"passed\":{passed},\"unsupported\":{unsupported},\"failed\":{failed},\"case_results\":[{}]}}",
        IDS.len(),
        rows.join(",")
    );
    if failed != 0 {
        std::process::exit(1);
    }
}
