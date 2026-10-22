//! Writer regressions with independent box coordinates and material volumes.

#![cfg(test)]

use std::collections::HashSet;

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision};
use gensolid_geom::{Curve, Line, Plane, Surface};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_step::{read_solids, write_block_step, write_planar_solid};
use gensolid_topo::{BodyId, Model, ShellId, body_volume};
use proptest::prelude::*;

#[derive(Clone, Copy)]
struct BoxSpec {
    frame: Frame,
    sizes: Vec3,
}

impl BoxSpec {
    fn volume(self) -> f64 {
        self.sizes.x * self.sizes.y * self.sizes.z
    }

    fn area(self) -> f64 {
        2.0 * (self.sizes.x * self.sizes.y
            + self.sizes.x * self.sizes.z
            + self.sizes.y * self.sizes.z)
    }

    fn corners(self) -> Vec<Point3> {
        let mut points = Vec::new();
        for x in [-0.5, 0.5] {
            for y in [-0.5, 0.5] {
                for z in [-0.5, 0.5] {
                    points.push(self.frame.to_world(Vec3::new(
                        x * self.sizes.x,
                        y * self.sizes.y,
                        z * self.sizes.z,
                    )));
                }
            }
        }
        points
    }
}

fn body(model: &Model) -> BodyId {
    assert_eq!(model.bodies().count(), 1);
    model.bodies().next().unwrap()
}

fn references(text: &str) -> Vec<u32> {
    text.split('#')
        .skip(1)
        .map(|tail| {
            let end = tail.bytes().take_while(u8::is_ascii_digit).count();
            tail[..end].parse().unwrap()
        })
        .collect()
}

fn shift_references(text: &str, offset: u32) -> String {
    let mut parts = text.split('#');
    let mut result = parts.next().unwrap().to_owned();
    for tail in parts {
        let end = tail.bytes().take_while(u8::is_ascii_digit).count();
        let id: u32 = tail[..end].parse().unwrap();
        result.push_str(&format!("#{}{}", id + offset, &tail[end..]));
    }
    result
}

// Compose independently placed solid boxes into a conventional STEP void input.
// Every cavity shell initially has outward faces and a false oriented-shell use.
fn input(boxes: &[BoxSpec]) -> Model {
    let mut lines = Vec::new();
    let mut offset = 0;
    let mut outer = 0;
    let mut cavities = Vec::new();
    let mut header = String::new();
    for (index, spec) in boxes.iter().enumerate() {
        let step = write_block_step(
            spec.sizes.x,
            spec.sizes.y,
            spec.sizes.z,
            &spec.frame,
            Precision::DEFAULT,
        )
        .unwrap();
        if index == 0 {
            header = step.split("DATA;").next().unwrap().to_owned();
        }
        let data = step.split("DATA;\n").nth(1).unwrap();
        let data = data.split("ENDSEC;").next().unwrap();
        let shell = data
            .lines()
            .find(|line| line.contains("=CLOSED_SHELL("))
            .map(|line| references(line)[0] + offset)
            .unwrap();
        for line in data
            .lines()
            .filter(|line| !line.contains("=MANIFOLD_SOLID_BREP("))
        {
            lines.push(shift_references(line, offset));
        }
        offset += references(data).into_iter().max().unwrap();
        if index == 0 {
            outer = shell;
        } else {
            offset += 1;
            lines.push(format!(
                "#{offset}=ORIENTED_CLOSED_SHELL('',*,#{shell},.F.);"
            ));
            cavities.push(format!("#{offset}"));
        }
    }
    offset += 1;
    lines.push(format!(
        "#{offset}=BREP_WITH_VOIDS('',#{outer},({}));",
        cavities.join(",")
    ));
    read_solids(
        &format!(
            "{header}DATA;\n{}\nENDSEC;\nEND-ISO-10303-21;\n",
            lines.join("\n")
        ),
        Precision::DEFAULT,
    )
    .unwrap()
}

fn checked(model: &Model) {
    let report = check_body(model, body(model), &CheckOptions::default()).unwrap();
    assert!(report.is_valid(), "{report:?}");
}

fn assert_geometry(model: &Model, boxes: &[BoxSpec]) {
    checked(model);
    let id = body(model);
    let regions = &model.body(id).unwrap().regions;
    assert_eq!(regions.len(), boxes.len() + 1);
    assert!(!model.region(regions[0]).unwrap().solid);
    let material = model.region(regions[1]).unwrap();
    assert!(material.solid);
    assert_eq!(material.shells.len(), boxes.len());
    assert_eq!(model.body_faces(id).unwrap().len(), boxes.len() * 6);
    assert_eq!(model.body_edges(id).unwrap().len(), boxes.len() * 12);
    assert_eq!(model.body_vertices(id).unwrap().len(), boxes.len() * 8);
    let tolerance = model.precision().linear() * 4.0;
    for (index, (shell, spec)) in material.shells.iter().zip(boxes).enumerate() {
        let mut vertices = HashSet::new();
        let shell_value = model.shell(*shell).unwrap();
        assert_eq!(shell_value.faces.len(), 6);
        let empty = model
            .region(regions[if index == 0 { 0 } else { index + 1 }])
            .unwrap();
        assert!(!empty.solid);
        assert_eq!(empty.shells.len(), 1);
        for face in &shell_value.faces {
            let face = model.face(*face).unwrap();
            assert_eq!(face.back, Some(*shell));
            assert_eq!(face.front, Some(empty.shells[0]));
            for lp in &face.loops {
                for fin in &model.loop_(*lp).unwrap().fins {
                    vertices.insert(model.fin_vertices(*fin).unwrap().0.unwrap());
                }
            }
        }
        let mut expected = spec.corners();
        assert_eq!(vertices.len(), expected.len());
        for vertex in vertices {
            let point = model.vertex(vertex).unwrap().point;
            let match_index = expected
                .iter()
                .position(|corner| (point - *corner).length() <= tolerance)
                .unwrap();
            expected.swap_remove(match_index);
        }
        assert!(expected.is_empty());
    }
    let expected = boxes[0].volume() - boxes[1..].iter().map(|b| b.volume()).sum::<f64>();
    let budget = model.precision().linear() * boxes.iter().map(|b| b.area()).sum::<f64>();
    assert!((body_volume(model, id).unwrap() - expected).abs() <= budget);
}

fn round_trip(boxes: &[BoxSpec]) {
    let model = input(boxes);
    assert_geometry(&model, boxes);
    let text = write_planar_solid(&model, body(&model)).unwrap();
    assert_eq!(text.matches("=BREP_WITH_VOIDS(").count(), 1);
    assert_eq!(text.matches("=CLOSED_SHELL(").count(), boxes.len());
    let output = read_solids(&text, *model.precision()).unwrap();
    assert_geometry(&output, boxes);
    // A second export/import must preserve all shell and analytic invariants too.
    let text = write_planar_solid(&output, body(&output)).unwrap();
    assert_geometry(&read_solids(&text, *model.precision()).unwrap(), boxes);
}

fn placed(frame: Frame, sizes: Vec3, center: Vec3) -> BoxSpec {
    BoxSpec {
        frame: Frame::new(
            frame.to_world(center),
            frame.z(),
            frame.x(),
            &Precision::DEFAULT,
        )
        .unwrap(),
        sizes,
    }
}

// Move all carriers along with the closed boundary. These inputs still pass the
// immutable checker, so rejection must come from shell-relationship certification.
fn translate_shell(model: &mut Model, shell: ShellId, delta: Vec3) {
    let precision = *model.precision();
    let faces = model.shell(shell).unwrap().faces.clone();
    let mut vertices = HashSet::new();
    let mut edges = HashSet::new();
    for face in &faces {
        for lp in &model.face(*face).unwrap().loops {
            for fin in &model.loop_(*lp).unwrap().fins {
                let edge = model.fin(*fin).unwrap().edge;
                edges.insert(edge);
                let edge = model.edge(edge).unwrap();
                vertices.extend([edge.start.unwrap(), edge.end.unwrap()]);
            }
        }
    }
    for vertex in vertices {
        let point = model.vertex(vertex).unwrap().point;
        model.raw().vertex_mut(vertex).unwrap().point = point + delta;
    }
    for edge in edges {
        let value = model.edge(edge).unwrap();
        let line = Line::through(
            model.vertex(value.start.unwrap()).unwrap().point,
            model.vertex(value.end.unwrap()).unwrap().point,
            &precision,
        )
        .unwrap();
        let curve = model.raw().geom_mut().add_curve(Curve::Line(line)).unwrap();
        model.raw().edge_mut(edge).unwrap().curve = Some(curve);
    }
    for face in faces {
        let surface = model.face(face).unwrap().surface.unwrap();
        let Surface::Plane(plane) = model.geom().surface(surface).unwrap() else {
            panic!("test input must be planar");
        };
        let frame = Frame::new(
            plane.frame().origin() + delta,
            plane.frame().z(),
            plane.frame().x(),
            &precision,
        )
        .unwrap();
        let surface = model
            .raw()
            .geom_mut()
            .add_surface(Surface::Plane(Plane::new(frame)))
            .unwrap();
        model.raw().face_mut(face).unwrap().surface = Some(surface);
    }
}

fn frame_at(translation: (f64, f64, f64), axis: (f64, f64), roll: f64) -> Frame {
    let p = Precision::DEFAULT;
    let frame = Frame::from_axis(
        Point3::new(translation.0, translation.1, translation.2),
        Vec3::new(axis.0, axis.1, 1.0),
        &p,
    )
    .unwrap();
    Frame::new(
        frame.origin(),
        frame.z(),
        frame.x() * roll.cos() + frame.y() * roll.sin(),
        &p,
    )
    .unwrap()
}

#[test]
fn independent_occt_fixture_exports_all_boundaries() {
    let model = read_solids(
        include_str!("fixtures/hollow_block.step"),
        Precision::DEFAULT,
    )
    .unwrap();
    let boxes = [
        BoxSpec {
            frame: Frame::WORLD,
            sizes: Vec3::new(4.0, 6.0, 8.0),
        },
        BoxSpec {
            frame: Frame::WORLD,
            sizes: Vec3::new(2.0, 3.0, 4.0),
        },
    ];
    let text = write_planar_solid(&model, body(&model)).unwrap();
    assert_geometry(&read_solids(&text, Precision::DEFAULT).unwrap(), &boxes);
}

#[test]
fn inconsistent_region_and_shell_memberships_are_rejected() {
    let original = read_solids(
        include_str!("fixtures/hollow_block.step"),
        Precision::DEFAULT,
    )
    .unwrap();
    let id = body(&original);
    let regions = original.body(id).unwrap().regions.clone();
    let boundary = original.region(regions[1]).unwrap().shells[1];
    for corruption in 0..6 {
        let mut model = original.clone();
        let face = model.shell(boundary).unwrap().faces[0];
        match corruption {
            0 => model.raw().body_mut(id).unwrap().regions.push(regions[2]),
            1 => model
                .raw()
                .region_mut(regions[1])
                .unwrap()
                .shells
                .push(boundary),
            2 => model.raw().shell_mut(boundary).unwrap().faces.push(face),
            3 => model.raw().face_mut(face).unwrap().front = None,
            4 => {
                let front = model.face(face).unwrap().front.unwrap();
                model.raw().shell_mut(front).unwrap().faces.pop();
            }
            _ => {
                let vertex = model.body_vertices(id).unwrap()[0];
                model
                    .raw()
                    .shell_mut(boundary)
                    .unwrap()
                    .acorn_vertices
                    .push(vertex);
            }
        }
        assert!(
            write_planar_solid(&model, id).is_err(),
            "corruption {corruption}"
        );
    }
}

#[test]
fn touching_overlapping_and_nested_closed_cavities_are_rejected() {
    let frame = Frame::WORLD;
    let original = input(&[
        BoxSpec {
            frame,
            sizes: Vec3::new(12.0, 8.0, 6.0),
        },
        placed(frame, Vec3::new(3.0, 3.0, 3.0), Vec3::new(-2.0, 0.0, 0.0)),
        placed(frame, Vec3::new(2.0, 2.0, 2.0), Vec3::new(2.0, 0.0, 0.0)),
    ]);
    let id = body(&original);
    let shell = original
        .region(original.body(id).unwrap().regions[1])
        .unwrap()
        .shells[2];
    for (target, expected) in [
        (5.0, ErrorCode::Unimplemented),   // cavity touches the outer face
        (9.0, ErrorCode::InvalidArgument), // cavity is wholly outside
        (0.5, ErrorCode::Unimplemented),   // cavity faces touch each other
        (0.0, ErrorCode::Unimplemented),   // cavity interiors overlap
        (-2.0, ErrorCode::Unimplemented),  // smaller cavity lies inside the other
    ] {
        let mut model = original.clone();
        translate_shell(&mut model, shell, Vec3::X * (target - 2.0));
        checked(&model);
        assert_eq!(write_planar_solid(&model, id).unwrap_err().code(), expected);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn separated_cavities_survive_scales_frames_and_repeated_export(
        factors in (1.0f64..3.0, 1.0f64..3.0, 1.0f64..3.0),
        exponent in -4i32..=2,
        translation in (-10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0),
        axis in (-2.0f64..2.0, -2.0f64..2.0),
        roll in -3.0f64..3.0,
        count in 1usize..=5,
        near_precision in any::<bool>(),
    ) {
        let scale = if near_precision { Precision::DEFAULT.linear() * 128.0 } else { 10.0f64.powi(exponent) };
        let frame = frame_at(translation, axis, roll);
        let sizes = Vec3::new(scale * factors.0 * 2.0, scale * factors.1 * 2.0, scale * factors.2 * 2.0);
        let mut boxes = vec![BoxSpec { frame, sizes }];
        let slot = sizes.x / count as f64;
        for index in 0..count {
            let center = Vec3::new(-sizes.x * 0.5 + slot * (index as f64 + 0.5), 0.0, 0.0);
            boxes.push(placed(frame, Vec3::new(slot * 0.25, sizes.y * 0.25, sizes.z * 0.25), center));
        }
        round_trip(&boxes);
    }

    #[test]
    fn strictly_contained_cavities_near_session_margin_export(
        gap_factor in 4.0f64..32.0,
        axis in (-1.0f64..1.0, -1.0f64..1.0),
        roll in -3.0f64..3.0,
        scale in 0.25f64..4.0,
    ) {
        let frame = frame_at((1.0, 2.0, 3.0), axis, roll);
        let gap = Precision::DEFAULT.linear() * gap_factor;
        let boxes = [
            BoxSpec { frame, sizes: Vec3::new(8.0 * scale, 6.0 * scale, 4.0 * scale) },
            placed(frame, Vec3::new(2.0 * scale, 2.0 * scale, 2.0 * scale), Vec3::new(3.0 * scale - gap, 0.0, 0.0)),
        ];
        round_trip(&boxes);
    }

    #[test]
    fn containment_and_separation_boundaries_are_refused_after_valid_carrier_moves(
        axis in (-1.0f64..1.0, -1.0f64..1.0),
        roll in -3.0f64..3.0,
        gap_factor in -2.0f64..0.75,
        outside in any::<bool>(),
    ) {
        let frame = frame_at((1.0, -2.0, 3.0), axis, roll);
        let boxes = [
            BoxSpec { frame, sizes: Vec3::new(12.0, 8.0, 6.0) },
            placed(frame, Vec3::new(2.0, 2.0, 2.0), Vec3::new(-2.0, 0.0, 0.0)),
            placed(frame, Vec3::new(2.0, 2.0, 2.0), Vec3::new(2.0, 0.0, 0.0)),
        ];
        let original = input(&boxes);
        let id = body(&original);
        let shells = original.region(original.body(id).unwrap().regions[1]).unwrap().shells.clone();
        let gap = Precision::DEFAULT.linear() * gap_factor;
        let mut containment = original.clone();
        let target = if outside { 9.0 } else { 5.0 - gap };
        translate_shell(&mut containment, shells[2], frame.x() * (target - 2.0));
        checked(&containment);
        let code = write_planar_solid(&containment, id).unwrap_err().code();
        prop_assert!(matches!(code, ErrorCode::InvalidArgument | ErrorCode::Unimplemented));
        let mut separation = original;
        translate_shell(&mut separation, shells[2], frame.x() * (-2.0 + gap));
        checked(&separation);
        prop_assert_eq!(write_planar_solid(&separation, id).unwrap_err().code(), ErrorCode::Unimplemented);
    }
}
