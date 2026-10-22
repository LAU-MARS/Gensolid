//! Raw STEP conformance checks, independent of shell interpretation by the reader.
//!
//! References: ISO 10303-42 closed_shell IP10 and brep_with_voids 6.4.3:
//! https://www.steptools.com/stds/smrl/data/resource_docs/geometric_and_topological_representation/sys/5_schema.htm
//! https://www.steptools.com/stds/smrl/data/resource_docs/geometric_and_topological_representation/sys/6_schema.htm

#![cfg(test)]

use std::collections::{HashMap, HashSet};

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::Precision;
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{Model, body_volume};
use proptest::prelude::*;

use crate::parse::{Instance, Value, bool_of, name_skipped, nums, parse_document, ref_of, refs};
use crate::{Writer, read_solids, write_block_step, write_planar_solid};

#[derive(Clone, Copy)]
struct AnalyticBox {
    frame: Frame,
    sizes: Vec3,
}

impl AnalyticBox {
    fn volume(self) -> f64 {
        self.sizes.x * self.sizes.y * self.sizes.z
    }

    fn area(self) -> f64 {
        2.0 * (self.sizes.x * self.sizes.y
            + self.sizes.y * self.sizes.z
            + self.sizes.z * self.sizes.x)
    }

    fn corners(self) -> Vec<Point3> {
        let mut result = Vec::new();
        for x in [-0.5, 0.5] {
            for y in [-0.5, 0.5] {
                for z in [-0.5, 0.5] {
                    result.push(self.frame.to_world(Vec3::new(
                        x * self.sizes.x,
                        y * self.sizes.y,
                        z * self.sizes.z,
                    )));
                }
            }
        }
        result
    }
}

type Entities = HashMap<u32, Instance>;

fn args<'a>(data: &'a Entities, id: u32, kind: &str) -> &'a [Value] {
    assert_eq!(data[&id].type_name, kind);
    name_skipped(&data[&id].args).unwrap()
}

fn vector(data: &Entities, id: u32, kind: &str) -> Vec3 {
    let v = nums(&args(data, id, kind)[0]).unwrap();
    assert_eq!(v.len(), 3);
    Vec3::new(v[0], v[1], v[2])
}

fn point(data: &Entities, vertex: u32) -> Point3 {
    let id = ref_of(&args(data, vertex, "VERTEX_POINT")[0]).unwrap();
    let v = vector(data, id, "CARTESIAN_POINT");
    Point3::new(v.x, v.y, v.z)
}

struct ShellMeasure {
    volume: f64,
    points: Vec<Point3>,
    normal_projections: Vec<f64>,
}

// Inspect ordinary CLOSED_SHELL geometry without resolving any oriented-shell use
// and without asking the STEP reader to interpret normals or loops.
fn outward_base_shell(data: &Entities, shell: u32, spec: AnalyticBox) -> ShellMeasure {
    let p = Precision::DEFAULT;
    let faces = refs(&args(data, shell, "CLOSED_SHELL")[0]).unwrap();
    assert_eq!(faces.len(), 6);
    assert_eq!(faces.iter().collect::<HashSet<_>>().len(), 6);
    let mut vertices = HashSet::new();
    let mut edge_uses: HashMap<u32, Vec<bool>> = HashMap::new();
    let mut volume = 0.0;
    let mut normal_projections = Vec::new();
    let mut expected_normals = vec![
        spec.frame.x(),
        -spec.frame.x(),
        spec.frame.y(),
        -spec.frame.y(),
        spec.frame.z(),
        -spec.frame.z(),
    ];
    for face in faces {
        let face = args(data, face, "ADVANCED_FACE");
        let bounds = refs(&face[0]).unwrap();
        assert_eq!(bounds.len(), 1);
        let bound = args(data, bounds[0], "FACE_OUTER_BOUND");
        assert!(bool_of(&bound[1]).unwrap());
        let coedges = refs(&args(data, ref_of(&bound[0]).unwrap(), "EDGE_LOOP")[0]).unwrap();
        assert_eq!(coedges.len(), 4);
        let mut ends = Vec::new();
        for coedge in coedges {
            let coedge = args(data, coedge, "ORIENTED_EDGE");
            let edge_id = ref_of(&coedge[2]).unwrap();
            let edge = args(data, edge_id, "EDGE_CURVE");
            let forward = bool_of(&coedge[3]).unwrap();
            edge_uses.entry(edge_id).or_default().push(forward);
            let a = ref_of(&edge[0]).unwrap();
            let b = ref_of(&edge[1]).unwrap();
            ends.push(if forward { (a, b) } else { (b, a) });
            vertices.extend([a, b]);
        }
        for i in 0..ends.len() {
            assert_eq!(
                ends[i].1,
                ends[(i + 1) % ends.len()].0,
                "broken co-edge cycle"
            );
        }
        let points: Vec<_> = ends.iter().map(|(a, _)| point(data, *a)).collect();
        let plane = args(data, ref_of(&face[1]).unwrap(), "PLANE");
        let placement = args(data, ref_of(&plane[0]).unwrap(), "AXIS2_PLACEMENT_3D");
        let axis = vector(data, ref_of(&placement[1]).unwrap(), "DIRECTION");
        let normal = axis
            * if bool_of(&face[2]).unwrap() {
                1.0
            } else {
                -1.0
            };
        let projection = (points[0] - spec.frame.origin()).dot(normal);
        assert!(
            projection > p.linear(),
            "raw CLOSED_SHELL face points into its finite region"
        );
        normal_projections.push(projection);
        let matching = expected_normals
            .iter()
            .position(|n| (*n - normal).length() <= p.angular() * 8.0)
            .unwrap();
        expected_normals.swap_remove(matching);

        // Exact planar polygon boundary integrals, with the box center as origin
        // to avoid cancellation under translation. No mesh or reader volume is used.
        let mut area_vector = Vec3::ZERO;
        for i in 0..points.len() {
            area_vector += (points[i] - spec.frame.origin())
                .cross(points[(i + 1) % points.len()] - spec.frame.origin())
                * 0.5;
        }
        assert!(
            area_vector.dot(normal) > p.linear() * p.linear(),
            "loop winding disagrees with face normal"
        );
        assert!(
            (area_vector - normal * area_vector.length()).length()
                <= p.linear() * spec.sizes.length() * 8.0
        );
        volume += (points[0] - spec.frame.origin()).dot(area_vector) / 3.0;
    }
    assert!(expected_normals.is_empty());
    assert_eq!(edge_uses.len(), 12);
    for uses in edge_uses.values() {
        assert_eq!(uses.len(), 2);
        assert_ne!(uses[0], uses[1], "shared edge uses must oppose each other");
    }
    assert_eq!(vertices.len(), 8);
    let points: Vec<_> = vertices.into_iter().map(|v| point(data, v)).collect();
    let mut corners = spec.corners();
    for actual in &points {
        let matching = corners
            .iter()
            .position(|c| (*actual - *c).length() <= p.linear() * 8.0)
            .unwrap();
        corners.swap_remove(matching);
    }
    assert!(corners.is_empty());
    assert!(
        (volume - spec.volume()).abs() <= p.linear() * spec.area() * 8.0,
        "raw shell signed volume: {volume}"
    );
    ShellMeasure {
        volume,
        points,
        normal_projections,
    }
}

fn checked_volume(model: &Model, expected: f64, boxes: &[AnalyticBox]) {
    let bodies: Vec<_> = model.bodies().collect();
    assert_eq!(bodies.len(), 1);
    let id = bodies[0];
    let report = check_body(model, id, &CheckOptions::default()).unwrap();
    assert!(report.is_valid(), "{report:?}");
    assert_eq!(report.unchecked, 0);
    assert_eq!(model.body(id).unwrap().regions.len(), boxes.len() + 1);
    assert_eq!(model.body_faces(id).unwrap().len(), boxes.len() * 6);
    assert_eq!(model.body_edges(id).unwrap().len(), boxes.len() * 12);
    assert_eq!(model.body_vertices(id).unwrap().len(), boxes.len() * 8);
    let budget = model.precision().linear() * boxes.iter().map(|b| b.area()).sum::<f64>() * 8.0;
    assert!((body_volume(model, id).unwrap() - expected).abs() <= budget);
}

fn verify_export(model: &Model, boxes: &[AnalyticBox]) {
    let p = Precision::DEFAULT;
    let expected = boxes[0].volume() - boxes[1..].iter().map(|b| b.volume()).sum::<f64>();
    checked_volume(model, expected, boxes);
    let text = write_planar_solid(model, model.bodies().next().unwrap()).unwrap();
    let data = parse_document(&text).unwrap();
    let count = |kind| data.values().filter(|e| e.type_name == kind).count();
    assert_eq!(count("CLOSED_SHELL"), boxes.len());
    assert_eq!(count("ORIENTED_CLOSED_SHELL"), boxes.len() - 1);
    assert_eq!(count("ADVANCED_FACE"), boxes.len() * 6);
    assert_eq!(count("EDGE_LOOP"), boxes.len() * 6);
    assert_eq!(count("ORIENTED_EDGE"), boxes.len() * 24);
    assert_eq!(count("EDGE_CURVE"), boxes.len() * 12);
    assert_eq!(count("VERTEX_POINT"), boxes.len() * 8);
    let (outer, wrappers) = if boxes.len() == 1 {
        assert_eq!(count("BREP_WITH_VOIDS"), 0);
        assert_eq!(count("MANIFOLD_SOLID_BREP"), 1);
        let entity = data
            .values()
            .find(|e| e.type_name == "MANIFOLD_SOLID_BREP")
            .unwrap();
        (
            ref_of(&name_skipped(&entity.args).unwrap()[0]).unwrap(),
            Vec::new(),
        )
    } else {
        assert_eq!(count("BREP_WITH_VOIDS"), 1);
        assert_eq!(count("MANIFOLD_SOLID_BREP"), 0);
        let entity = data
            .values()
            .find(|e| e.type_name == "BREP_WITH_VOIDS")
            .unwrap();
        let values = name_skipped(&entity.args).unwrap();
        (ref_of(&values[0]).unwrap(), refs(&values[1]).unwrap())
    };
    assert_eq!(wrappers.len(), boxes.len() - 1);
    assert_eq!(
        wrappers.iter().collect::<HashSet<_>>().len(),
        wrappers.len()
    );
    let mut base_shells = HashSet::from([outer]);
    let mut effective_volume = outward_base_shell(&data, outer, boxes[0]).volume;
    let mut cavity_intervals: Vec<(f64, f64)> = Vec::new();
    for (wrapper, spec) in wrappers.iter().zip(&boxes[1..]) {
        let wrapper = args(&data, *wrapper, "ORIENTED_CLOSED_SHELL");
        let shell = ref_of(&wrapper[1]).unwrap();
        assert!(base_shells.insert(shell));
        // Check the base FIRST: a true wrapper around an inward base preserves
        // geometry but violates the CLOSED_SHELL requirement independently.
        let raw = outward_base_shell(&data, shell, *spec);
        let orientation = bool_of(&wrapper[2]).unwrap();
        assert!(!orientation, "BREP_WITH_VOIDS requires false void uses");
        let sign = if orientation { 1.0 } else { -1.0 };
        for projection in raw.normal_projections {
            assert!(
                sign * projection < -p.linear(),
                "effective void normal must point into cavity"
            );
        }
        effective_volume += sign * raw.volume;
        let mut interval = (f64::INFINITY, f64::NEG_INFINITY);
        for point in raw.points {
            let local = boxes[0].frame.to_local(point);
            let half = boxes[0].sizes * 0.5;
            assert!(half.x - local.x.abs() > p.linear());
            assert!(half.y - local.y.abs() > p.linear());
            assert!(half.z - local.z.abs() > p.linear());
            interval.0 = interval.0.min(local.x);
            interval.1 = interval.1.max(local.x);
        }
        for prior in &cavity_intervals {
            assert!(
                interval.0 - prior.1 > p.linear() || prior.0 - interval.1 > p.linear(),
                "analytic cavity projections must be disjoint"
            );
        }
        cavity_intervals.push(interval);
    }
    assert_eq!(base_shells.len(), count("CLOSED_SHELL"));
    let budget = p.linear() * boxes.iter().map(|b| b.area()).sum::<f64>() * 8.0;
    assert!((effective_volume - expected).abs() <= budget);
    // Secondary regression only: raw conformance above does not depend on this reader.
    checked_volume(&read_solids(&text, p).unwrap(), expected, boxes);
}

// Assemble conventional outward box shells with false void uses as input only.
// Output conformance is measured from coordinates and analytic boxes above.
fn box_input(boxes: &[AnalyticBox]) -> Model {
    let mut combined = Writer::default();
    let mut shells = Vec::new();
    for spec in boxes {
        let step = write_block_step(
            spec.sizes.x,
            spec.sizes.y,
            spec.sizes.z,
            &spec.frame,
            Precision::DEFAULT,
        )
        .unwrap();
        let offset = combined.next;
        let records = step
            .split("DATA;\n")
            .nth(1)
            .unwrap()
            .split("ENDSEC;")
            .next()
            .unwrap();
        for record in records
            .lines()
            .filter(|s| !s.contains("=MANIFOLD_SOLID_BREP("))
        {
            let mut shifted = String::new();
            let mut parts = record.split('#');
            shifted.push_str(parts.next().unwrap());
            for part in parts {
                let end = part.bytes().take_while(u8::is_ascii_digit).count();
                let id = part[..end].parse::<u32>().unwrap() + offset;
                combined.next = combined.next.max(id);
                shifted.push_str(&format!("#{id}{}", &part[end..]));
            }
            if record.contains("=CLOSED_SHELL(") {
                shells.push(combined.next);
            }
            combined.lines.push(shifted);
        }
    }
    if shells.len() == 1 {
        combined.entity(format!("MANIFOLD_SOLID_BREP('',#{})", shells[0]));
    } else {
        let mut uses = Vec::new();
        for shell in &shells[1..] {
            let id = combined.entity(format!("ORIENTED_CLOSED_SHELL('',*,#{shell},.F.)"));
            uses.push(format!("#{id}"));
        }
        combined.entity(format!(
            "BREP_WITH_VOIDS('',#{},({}))",
            shells[0],
            uses.join(",")
        ));
    }
    read_solids(&combined.finish(), Precision::DEFAULT).unwrap()
}

#[test]
fn single_shell_keeps_outward_faces_without_void_wrappers() {
    let spec = AnalyticBox {
        frame: Frame::new(
            Point3::new(7.0, -4.0, 2.0),
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::X,
            &Precision::DEFAULT,
        )
        .unwrap(),
        sizes: Vec3::new(2.0, 3.0, 4.0),
    };
    verify_export(&box_input(&[spec]), &[spec]);
}

#[test]
fn independent_fixture_has_outward_base_and_false_void_use() {
    let model = read_solids(
        include_str!("../tests/fixtures/hollow_block.step"),
        Precision::DEFAULT,
    )
    .unwrap();
    verify_export(
        &model,
        &[
            AnalyticBox {
                frame: Frame::WORLD,
                sizes: Vec3::new(4.0, 6.0, 8.0),
            },
            AnalyticBox {
                frame: Frame::WORLD,
                sizes: Vec3::new(2.0, 3.0, 4.0),
            },
        ],
    );
}

#[test]
fn close_separated_cavities_and_thin_outer_wall_keep_orientation() {
    let p = Precision::DEFAULT;
    let gap = p.linear() * 128.0;
    let boxes = [
        AnalyticBox {
            frame: Frame::WORLD,
            sizes: Vec3::new(8.0, 6.0, 4.0),
        },
        AnalyticBox {
            frame: Frame::from_axis(Point3::new(-1.0 - gap * 0.5, 0.0, 0.0), Vec3::Z, &p).unwrap(),
            sizes: Vec3::new(2.0, 2.0, 2.0),
        },
        AnalyticBox {
            frame: Frame::from_axis(Point3::new(1.0 + gap * 0.5, 0.0, 0.0), Vec3::Z, &p).unwrap(),
            sizes: Vec3::new(2.0, 2.0, 2.0),
        },
    ];
    verify_export(&box_input(&boxes), &boxes);
    let thin = [
        boxes[0],
        AnalyticBox {
            frame: Frame::from_axis(Point3::new(3.0 - gap, 0.0, 0.0), Vec3::Z, &p).unwrap(),
            sizes: Vec3::new(2.0, 2.0, 2.0),
        },
    ];
    verify_export(&box_input(&thin), &thin);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn raw_shells_conform_under_random_sizes_rigid_frames_and_cavity_counts(
        factors in (1.0f64..2.0, 1.0f64..2.0, 1.0f64..2.0),
        translation in (-20.0f64..20.0, -20.0f64..20.0, -20.0f64..20.0),
        axes in (-2.0f64..2.0, -2.0f64..2.0, -2.0f64..2.0, -2.0f64..2.0),
        roll in -3.0f64..3.0,
        exponent in -3i32..=1,
        count in 0usize..=4,
        tiny in any::<bool>(),
    ) {
        let p = Precision::DEFAULT;
        let scale = if tiny { p.linear() * 256.0 } else { 10.0f64.powi(exponent) };
        let axis_frame = Frame::from_axis(Point3::new(translation.0, translation.1, translation.2), Vec3::new(axes.0, axes.1, 1.0), &p).unwrap();
        let frame = Frame::new(axis_frame.origin(), axis_frame.z(), axis_frame.x() * roll.cos() + axis_frame.y() * roll.sin(), &p).unwrap();
        let outer = AnalyticBox { frame, sizes: Vec3::new(8.0 * count.max(1) as f64 * scale * factors.0, 10.0 * scale * factors.1, 12.0 * scale * factors.2) };
        let mut boxes = vec![outer];
        for i in 0..count {
            let slot = outer.sizes.x / count as f64;
            let origin = frame.to_world(Vec3::new(-outer.sizes.x * 0.5 + slot * (i as f64 + 0.5), scale * roll.sin(), scale * roll.cos()));
            let rotated = Frame::from_axis(origin, frame.vector_to_world(Vec3::new(axes.2, axes.3, 1.0)), &p).unwrap();
            let cavity = Frame::new(origin, rotated.z(), rotated.x() * roll.cos() + rotated.y() * roll.sin(), &p).unwrap();
            boxes.push(AnalyticBox { frame: cavity, sizes: Vec3::new(scale * factors.2, scale * factors.0, scale * factors.1) });
        }
        verify_export(&box_input(&boxes), &boxes);
    }
}
