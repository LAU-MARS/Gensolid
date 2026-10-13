//! Original test-only XT encoder for planar blocks. This is a framing/adapter fixture,
//! not a Parasolid-produced file and not an independent CAD oracle.
#![allow(clippy::unwrap_used)] // Failed fixture setup must fail the test immediately.
use std::collections::BTreeMap;

use gensolid_core::Precision;
use gensolid_geom::{Curve, Surface};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{EntityRef, Model, Sense, make_block};
use parasolid_core::{BuiltinProfileRegistry, FieldType, SchemaKey};

#[derive(Clone, Debug)]
enum Value {
    Int(i32),
    Real(f64),
    Char(u8),
    Vector([f64; 3]),
}

struct Node {
    kind: u16,
    index: i32,
    values: BTreeMap<&'static str, Value>,
}

impl Node {
    fn new(kind: u16, index: i32) -> Self {
        Self {
            kind,
            index,
            values: BTreeMap::new(),
        }
    }
    fn int(&mut self, key: &'static str, value: i32) {
        self.values.insert(key, Value::Int(value));
    }
    fn real(&mut self, key: &'static str, value: f64) {
        self.values.insert(key, Value::Real(value));
    }
    fn char(&mut self, key: &'static str, value: u8) {
        self.values.insert(key, Value::Char(value));
    }
    fn vector(&mut self, key: &'static str, value: Vec3) {
        self.values
            .insert(key, Value::Vector([value.x, value.y, value.z]));
    }
    fn point(&mut self, key: &'static str, value: Point3) {
        self.vector(key, value.to_vec());
    }
}

/// Serialize a kernel block to a deliberately small, schema-driven test stream.
pub fn block(dimensions: [f64; 3], frame: Frame, binary: bool) -> Vec<u8> {
    encode(&block_nodes(dimensions, frame), binary)
}

fn block_nodes(dimensions: [f64; 3], frame: Frame) -> Vec<Node> {
    let precision = Precision::DEFAULT;
    let mut model = Model::new(precision);
    let body = make_block(
        &mut model,
        dimensions[0],
        dimensions[1],
        dimensions[2],
        &frame,
    )
    .unwrap();
    let regions = model.body(body).unwrap().regions.clone();
    let shells: Vec<_> = regions
        .iter()
        .flat_map(|&id| model.region(id).unwrap().shells.iter().copied())
        .collect();
    let faces = model.body_faces(body).unwrap();
    let loops: Vec<_> = faces
        .iter()
        .flat_map(|&id| model.face(id).unwrap().loops.iter().copied())
        .collect();
    let fins: Vec<_> = loops
        .iter()
        .flat_map(|&id| model.loop_(id).unwrap().fins.iter().copied())
        .collect();
    let edges = model.body_edges(body).unwrap();
    let vertices = model.body_vertices(body).unwrap();
    let entities = std::iter::once(EntityRef::Body(body))
        .chain(regions.iter().copied().map(EntityRef::Region))
        .chain(shells.iter().copied().map(EntityRef::Shell))
        .chain(faces.iter().copied().map(EntityRef::Face))
        .chain(loops.iter().copied().map(EntityRef::Loop))
        .chain(fins.iter().copied().map(EntityRef::Fin))
        .chain(edges.iter().copied().map(EntityRef::Edge))
        .chain(vertices.iter().copied().map(EntityRef::Vertex))
        .collect::<Vec<_>>();
    let index =
        |e| i32::try_from(entities.iter().position(|&item| item == e).unwrap() + 1).unwrap();
    let point_start = i32::try_from(entities.len() + 1).unwrap();
    let curve_start = point_start + i32::try_from(vertices.len()).unwrap();
    let surface_start = curve_start + i32::try_from(edges.len()).unwrap();
    let neighbour = |values: &[EntityRef], i: usize, delta: isize| {
        i.checked_add_signed(delta)
            .and_then(|j| values.get(j))
            .map_or(0, |&id| index(id))
    };
    let vertex_refs: Vec<_> = vertices.iter().copied().map(EntityRef::Vertex).collect();
    let edge_refs: Vec<_> = edges.iter().copied().map(EntityRef::Edge).collect();
    let face_refs: Vec<_> = faces.iter().copied().map(EntityRef::Face).collect();
    let region_refs: Vec<_> = regions.iter().copied().map(EntityRef::Region).collect();
    let mut nodes = Vec::new();
    let mut b = Node::new(12, 1);
    b.int(
        "max_local_id",
        surface_start + i32::try_from(faces.len()).unwrap(),
    );
    b.real("size_precision", precision.size());
    b.real("linear_precision", precision.linear());
    b.int("body_kind", 1);
    b.int("geometry_state", 1);
    b.int("region_head", index(region_refs[0]));
    b.int("edge_head", index(edge_refs[0]));
    b.int("vertex_head", index(vertex_refs[0]));
    b.int("boundary_surfaces", surface_start);
    b.int("boundary_curves", curve_start);
    b.int("boundary_points", point_start);
    nodes.push(b);
    for (i, &id) in regions.iter().enumerate() {
        let r = model.region(id).unwrap();
        let mut n = Node::new(19, index(EntityRef::Region(id)));
        n.int("body_ref", 1);
        n.int("owner_ref", 1);
        n.int("shell_head", index(EntityRef::Shell(r.shells[0])));
        n.int("next_region", neighbour(&region_refs, i, 1));
        n.int("previous_region", neighbour(&region_refs, i, -1));
        n.char("region_kind", if r.solid { b'S' } else { b'V' });
        nodes.push(n);
    }
    for &id in &shells {
        let s = model.shell(id).unwrap();
        let mut n = Node::new(13, index(EntityRef::Shell(id)));
        n.int("region_ref", index(EntityRef::Region(s.region)));
        n.int(
            if model.region(s.region).unwrap().solid {
                "back_faces"
            } else {
                "front_faces"
            },
            index(face_refs[0]),
        );
        nodes.push(n);
    }
    for (i, &id) in faces.iter().enumerate() {
        let f = model.face(id).unwrap();
        let mut n = Node::new(14, index(EntityRef::Face(id)));
        n.int("loop_head", index(EntityRef::Loop(f.loops[0])));
        n.int("back_shell", index(EntityRef::Shell(f.back.unwrap())));
        n.int("front_shell", index(EntityRef::Shell(f.front.unwrap())));
        n.int("surface_ref", surface_start + i32::try_from(i).unwrap());
        n.char("orientation", b'+');
        for key in ["next_back", "next_front"] {
            n.int(key, neighbour(&face_refs, i, 1));
        }
        for key in ["previous_back", "previous_front"] {
            n.int(key, neighbour(&face_refs, i, -1));
        }
        nodes.push(n);
    }
    for &id in &loops {
        let lp = model.loop_(id).unwrap();
        let mut n = Node::new(15, index(EntityRef::Loop(id)));
        n.int("fin_ref", index(EntityRef::Fin(lp.fins[0])));
        n.int("face_ref", index(EntityRef::Face(lp.face)));
        nodes.push(n);
    }
    for &id in &fins {
        let f = model.fin(id).unwrap();
        let lp = model.loop_(f.owner).unwrap();
        let edge = model.edge(f.edge).unwrap();
        let fi = lp.fins.iter().position(|&item| item == id).unwrap();
        let end = model.fin_vertices(id).unwrap().1.unwrap();
        let vertex_fins: Vec<_> = fins
            .iter()
            .copied()
            .filter(|&id| model.fin_vertices(id).unwrap().1 == Some(end))
            .collect();
        let vi = vertex_fins.iter().position(|&item| item == id).unwrap();
        let mut n = Node::new(17, index(EntityRef::Fin(id)));
        n.int("loop_ref", index(EntityRef::Loop(f.owner)));
        n.int(
            "forward_fin",
            index(EntityRef::Fin(lp.fins[(fi + 1) % lp.fins.len()])),
        );
        n.int(
            "backward_fin",
            index(EntityRef::Fin(
                lp.fins[(fi + lp.fins.len() - 1) % lp.fins.len()],
            )),
        );
        n.int("end_vertex", index(EntityRef::Vertex(end)));
        n.int(
            "radial_fin",
            index(EntityRef::Fin(
                *edge.fins.iter().find(|&&other| other != id).unwrap(),
            )),
        );
        n.int("edge_ref", index(EntityRef::Edge(f.edge)));
        n.int(
            "next_vertex_fin",
            index(EntityRef::Fin(vertex_fins[(vi + 1) % vertex_fins.len()])),
        );
        n.char(
            "orientation",
            if f.sense == Sense::Positive {
                b'+'
            } else {
                b'-'
            },
        );
        nodes.push(n);
    }
    for (i, &id) in edges.iter().enumerate() {
        let e = model.edge(id).unwrap();
        let mut n = Node::new(16, index(EntityRef::Edge(id)));
        n.int("fin_ref", index(EntityRef::Fin(e.fins[0])));
        n.int("owner_ref", 1);
        n.int("curve_ref", curve_start + i32::try_from(i).unwrap());
        n.int("next_edge", neighbour(&edge_refs, i, 1));
        n.int("previous_edge", neighbour(&edge_refs, i, -1));
        nodes.push(n);
    }
    for (i, &id) in vertices.iter().enumerate() {
        let mut n = Node::new(18, index(EntityRef::Vertex(id)));
        n.int("owner_ref", 1);
        n.int("point_ref", point_start + i32::try_from(i).unwrap());
        let fin = fins
            .iter()
            .copied()
            .find(|&f| model.fin_vertices(f).unwrap().1 == Some(id))
            .unwrap();
        n.int("fin_head", index(EntityRef::Fin(fin)));
        n.int("next_vertex", neighbour(&vertex_refs, i, 1));
        n.int("previous_vertex", neighbour(&vertex_refs, i, -1));
        nodes.push(n);
        let mut p = Node::new(29, point_start + i32::try_from(i).unwrap());
        p.int("owner_ref", index(EntityRef::Vertex(id)));
        p.point("position", model.vertex(id).unwrap().point);
        p.int(
            "next_point",
            if i + 1 == vertices.len() {
                0
            } else {
                p.index + 1
            },
        );
        p.int("previous_point", if i == 0 { 0 } else { p.index - 1 });
        nodes.push(p);
    }
    for (i, &id) in edges.iter().enumerate() {
        let Curve::Line(line) = model
            .geom()
            .curve(model.edge(id).unwrap().curve.unwrap())
            .unwrap()
        else {
            unreachable!()
        };
        let mut n = Node::new(30, curve_start + i32::try_from(i).unwrap());
        n.int("owner_ref", index(EntityRef::Edge(id)));
        n.char("orientation", b'+');
        n.point("origin", line.origin());
        n.vector("tangent", line.direction());
        n.int(
            "next_curve",
            if i + 1 == edges.len() { 0 } else { n.index + 1 },
        );
        n.int("previous_curve", if i == 0 { 0 } else { n.index - 1 });
        nodes.push(n);
    }
    for (i, &id) in faces.iter().enumerate() {
        let Surface::Plane(plane) = model
            .geom()
            .surface(model.face(id).unwrap().surface.unwrap())
            .unwrap()
        else {
            unreachable!()
        };
        let mut n = Node::new(50, surface_start + i32::try_from(i).unwrap());
        n.int("owner_ref", index(EntityRef::Face(id)));
        n.char("orientation", b'+');
        n.point("origin", plane.frame().origin());
        n.vector("normal", plane.frame().z());
        n.vector("x_direction", plane.frame().x());
        n.int(
            "next_surface",
            if i + 1 == faces.len() { 0 } else { n.index + 1 },
        );
        n.int("previous_surface", if i == 0 { 0 } else { n.index - 1 });
        nodes.push(n);
    }
    nodes
}

/// A well-framed block whose first vertex no longer lies on its edge carriers.
pub fn block_with_invalid_geometry(binary: bool) -> Vec<u8> {
    let mut nodes = block_nodes([1.0, 2.0, 3.0], Frame::WORLD);
    let point = nodes.iter_mut().find(|node| node.kind == 29).unwrap();
    point.point("position", Point3::new(100.0, 100.0, 100.0));
    encode(&nodes, binary)
}

/// A well-framed block with a broken source radial fin ring.
pub fn block_with_invalid_topology(binary: bool) -> Vec<u8> {
    let mut nodes = block_nodes([1.0, 2.0, 3.0], Frame::WORLD);
    let fin = nodes.iter_mut().find(|node| node.kind == 17).unwrap();
    fin.int("radial_fin", 0);
    encode(&nodes, binary)
}

/// The same block with reversed shared-carrier senses and compensating use senses.
pub fn block_with_negative_carrier_senses(binary: bool) -> Vec<u8> {
    let mut nodes = block_nodes([1.0, 2.0, 3.0], Frame::WORLD);
    for node in &mut nodes {
        if node.kind == 17 {
            let Value::Char(value) = node.values["orientation"] else {
                unreachable!()
            };
            node.char("orientation", if value == b'+' { b'-' } else { b'+' });
        } else if matches!(node.kind, 14 | 30 | 50) {
            node.char("orientation", b'-');
        }
    }
    encode(&nodes, binary)
}

/// A consistent XT acorn whose geometry the immutable body checker does not visit.
pub fn acorn(binary: bool) -> Vec<u8> {
    let precision = Precision::DEFAULT;
    let mut body = Node::new(12, 1);
    body.real("size_precision", precision.size());
    body.real("linear_precision", precision.linear());
    body.int("body_kind", 2);
    body.int("region_head", 2);
    let mut region = Node::new(19, 2);
    region.int("body_ref", 1);
    region.int("owner_ref", 1);
    region.int("shell_head", 3);
    region.char("region_kind", b'V');
    let mut shell = Node::new(13, 3);
    shell.int("region_ref", 2);
    shell.int("isolated_vertex", 4);
    let mut vertex = Node::new(18, 4);
    vertex.int("point_ref", 5);
    vertex.int("owner_ref", 3);
    let mut point = Node::new(29, 5);
    point.int("owner_ref", 4);
    point.point("position", Point3::new(1.0, 2.0, 3.0));
    encode(&[body, region, shell, vertex, point], binary)
}

fn encode(nodes: &[Node], binary: bool) -> Vec<u8> {
    let key = SchemaKey::parse("SCH_3000000_30000").unwrap();
    let registry = BuiltinProfileRegistry::compiled().unwrap();
    let profile = registry.profile_for_key(&key).unwrap();
    let description = ": TRANSMIT FILE created by modeller version 3000000";
    let mut out = Vec::new();
    if binary {
        out.extend_from_slice(b"PS\0\0");
        out.extend_from_slice(&u16::try_from(description.len()).unwrap().to_be_bytes());
        out.extend_from_slice(description.as_bytes());
        out.extend_from_slice(&i32::try_from(key.raw().len()).unwrap().to_be_bytes());
        out.extend_from_slice(key.raw().as_bytes());
        out.extend_from_slice(&0_i32.to_be_bytes());
    } else {
        out.extend_from_slice(
            format!(
                "T{} {description}{} {}0 ",
                description.len(),
                key.raw().len(),
                key.raw()
            )
            .as_bytes(),
        );
    }
    for node in nodes {
        if binary {
            out.extend_from_slice(&node.kind.to_be_bytes());
            out.extend_from_slice(&i16::try_from(node.index + 1).unwrap().to_be_bytes());
        } else {
            out.extend_from_slice(format!("{} {} ", node.kind, node.index).as_bytes());
        }
        for field in &profile.definition(node.kind).unwrap().fields {
            assert!(field.transmitted && field.element_count == 0);
            let value =
                node.values
                    .get(field.name.as_str())
                    .cloned()
                    .unwrap_or(match field.field_type {
                        FieldType::Double => Value::Real(0.0),
                        _ => Value::Int(0),
                    });
            if binary {
                match (field.field_type, value) {
                    (FieldType::PointerIndex, Value::Int(v)) => {
                        out.extend_from_slice(&i16::try_from(v + 1).unwrap().to_be_bytes())
                    }
                    (FieldType::Integer, Value::Int(v)) => out.extend_from_slice(&v.to_be_bytes()),
                    (FieldType::UnsignedByte | FieldType::Logical, Value::Int(v)) => {
                        out.push(u8::try_from(v).unwrap())
                    }
                    (FieldType::Character, Value::Char(v)) => out.push(v),
                    (FieldType::Double, Value::Real(v)) => out.extend_from_slice(&v.to_be_bytes()),
                    (FieldType::Vector, Value::Vector(v)) => {
                        for component in v {
                            out.extend_from_slice(&component.to_be_bytes());
                        }
                    }
                    _ => unreachable!(),
                }
            } else {
                let text = match value {
                    Value::Int(v) => format!("{v} "),
                    Value::Real(v) => format!("{v} "),
                    Value::Char(v) => format!("{}", char::from(v)),
                    Value::Vector(v) => format!("{} {} {} ", v[0], v[1], v[2]),
                };
                out.extend_from_slice(text.as_bytes());
            }
        }
    }
    if binary {
        out.extend_from_slice(&[0, 1, 0, 1]);
    } else {
        out.extend_from_slice(b"1 0 ");
    }
    out
}
