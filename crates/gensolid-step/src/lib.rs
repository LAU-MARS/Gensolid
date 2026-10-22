//! STEP exchange for manifold solids with analytic or B-spline geometry.
//!
//! The reader preserves numeric file coordinates as model units and checks every imported
//! body. Planar convex solids with finite voids also require certified containment
//! and separation. Assemblies and unsupported classes are refused explicitly. The
//! AP203 writer emits planar solids with straight edges and no face holes, including
//! finite voids whose convex containment and separation can be certified.

mod geometry;
mod parse;
mod report;
mod voids;

pub use report::{StepStatus, StepValidityReport, validity_report};

use std::collections::{HashMap, HashSet};

use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{Curve, Line, Surface};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{
    Body, BodyId, BodyKind, Edge, Face, Fin, Loop, Model, Region, Sense, Shell, Vertex, make_block,
};
use parse::{Instance, Value, bool_of, name_skipped, nums, parse_document, ref_of, refs};

/// Writes `body` as an AP203 file. Faces must be planes and edges must be lines, with one
/// loop and no holes. Finite voids require connected convex planar shells with certified
/// strict containment and separation; curved or unproven relationships are refused.
pub fn write_planar_solid(model: &Model, body: BodyId) -> Result<String> {
    let boundaries = if model.body(body)?.regions.len() == 2 {
        None
    } else {
        Some(planar_void_export_shells(model, body)?)
    };
    let precision = *model.precision();
    let mut writer = Writer::default();
    let mut points = HashMap::new();
    let mut vertices = HashMap::new();
    for vertex_id in model.body_vertices(body)? {
        let point = model.vertex(vertex_id)?.point;
        let point_id = writer.cartesian(point);
        let vertex = writer.entity(format!("VERTEX_POINT('',#{point_id})"));
        points.insert(vertex_id, point_id);
        vertices.insert(vertex_id, vertex);
    }
    let mut edges = HashMap::new();
    for edge_id in model.body_edges(body)? {
        let edge = model.edge(edge_id)?;
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Err(Error::unsupported(
                "STEP export requires every edge to have two vertices",
            ));
        };
        let curve = edge
            .curve
            .ok_or_else(|| Error::unsupported("STEP export requires every edge to have a curve"))?;
        if !matches!(model.geom().curve(curve)?, Curve::Line(_)) {
            return Err(Error::unsupported(
                "STEP export currently writes line edges only",
            ));
        }
        let start_point = model.vertex(start)?.point;
        let end_point = model.vertex(end)?.point;
        let line = Line::through(start_point, end_point, &precision)?;
        let origin = *points
            .get(&start)
            .ok_or_else(|| Error::internal("missing vertex point"))?;
        let direction = writer.direction(line.direction());
        let vector = writer.entity(format!(
            "VECTOR('',#{direction},{length})",
            length = format_number((end_point - start_point).length())
        ));
        let curve = writer.entity(format!("LINE('',#{origin},#{vector})"));
        let edge_step = writer.entity(format!(
            "EDGE_CURVE('',#{start_v},#{end_v},#{curve},.T.)",
            start_v = vertices
                .get(&start)
                .ok_or_else(|| Error::internal("missing start vertex"))?,
            end_v = vertices
                .get(&end)
                .ok_or_else(|| Error::internal("missing end vertex"))?
        ));
        edges.insert(edge_id, edge_step);
    }
    let mut face_ids = Vec::new();
    let mut face_steps = HashMap::new();
    for face_id in model.body_faces(body)? {
        let face = model.face(face_id)?;
        if face.loops.len() != 1 {
            return Err(Error::unsupported(
                "STEP export currently writes faces with exactly one loop",
            ));
        }
        let surface = face.surface.ok_or_else(|| {
            Error::unsupported("STEP export requires every face to have a surface")
        })?;
        let Surface::Plane(plane) = model.geom().surface(surface)? else {
            return Err(Error::unsupported(
                "STEP export currently writes planar faces only",
            ));
        };
        let loop_id = face.loops[0];
        let mut oriented = Vec::new();
        for fin_id in &model.loop_(loop_id)?.fins {
            let fin = model.fin(*fin_id)?;
            let edge_step = edges
                .get(&fin.edge)
                .ok_or_else(|| Error::internal("missing edge"))?;
            let flag = if fin.sense == Sense::Positive {
                ".T."
            } else {
                ".F."
            };
            oriented.push(writer.entity(format!("ORIENTED_EDGE('',*,*,#{edge_step},{flag})")));
        }
        let list = oriented
            .iter()
            .map(|id| format!("#{id}"))
            .collect::<Vec<_>>()
            .join(",");
        let edge_loop = writer.entity(format!("EDGE_LOOP('',({list}))"));
        let bound = writer.entity(format!("FACE_OUTER_BOUND('',#{edge_loop},.T.)"));
        let placement = writer.placement(plane.frame());
        let plane_id = writer.entity(format!("PLANE('',#{placement})"));
        let sense = if face.sense == Sense::Positive {
            ".T."
        } else {
            ".F."
        };
        let step = writer.entity(format!("ADVANCED_FACE('',(#{bound}),#{plane_id},{sense})"));
        face_ids.push(step);
        face_steps.insert(face_id, step);
    }
    if let Some(boundaries) = boundaries {
        let mut shells = Vec::new();
        for boundary in boundaries {
            let faces = model
                .shell(boundary)?
                .faces
                .iter()
                .map(|face| {
                    face_steps
                        .get(face)
                        .map(|id| format!("#{id}"))
                        .ok_or_else(|| Error::internal("missing STEP boundary face"))
                })
                .collect::<Result<Vec<_>>>()?
                .join(",");
            shells.push(writer.entity(format!("CLOSED_SHELL('',({faces}))")));
        }
        let outer = shells
            .first()
            .ok_or_else(|| Error::internal("missing certified outer shell"))?;
        let mut cavities = Vec::new();
        for shell in shells.iter().skip(1) {
            // The stored face normals already point out of material and into the
            // cavity. Preserve both face senses and loop traversal with a true use.
            let oriented = writer.entity(format!("ORIENTED_CLOSED_SHELL('',*,#{shell},.T.)"));
            cavities.push(format!("#{oriented}"));
        }
        writer.entity(format!(
            "BREP_WITH_VOIDS('',#{outer},({}))",
            cavities.join(",")
        ));
        return Ok(writer.finish());
    }
    let faces = face_ids
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let shell = writer.entity(format!("CLOSED_SHELL('',({faces}))"));
    writer.entity(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    Ok(writer.finish())
}

/// Establish the region/shell pairing required by the existing void certifier.
/// The checker validates carriers and backpointers, but does not require distinct
/// boundary components or account for every side of every material face.
fn planar_void_export_shells(model: &Model, body: BodyId) -> Result<Vec<gensolid_topo::ShellId>> {
    let value = model.body(body)?;
    if value.kind != BodyKind::Solid || value.regions.len() < 3 {
        return Err(Error::unimplemented(
            "STEP void export needs one solid with finite empty regions",
        ));
    }
    let report = gensolid_check::check_body(model, body, &gensolid_check::CheckOptions::default())?;
    if !report.faults.is_empty() {
        return Err(Error::invalid_argument(format!(
            "STEP export body failed check_body: {:?}",
            report.faults
        )));
    }
    if report.unchecked != 0 {
        return Err(Error::unimplemented(
            "STEP void export requires complete geometric checks",
        ));
    }
    let mut material = None;
    let mut regions = HashSet::new();
    for id in &value.regions {
        if !regions.insert(*id) {
            return Err(Error::invalid_argument("STEP body repeats a region"));
        }
        if model.region(*id)?.solid && material.replace(*id).is_some() {
            return Err(Error::unimplemented(
                "STEP void export needs one connected material region",
            ));
        }
    }
    let material =
        material.ok_or_else(|| Error::invalid_argument("STEP void body has no material region"))?;
    let boundaries = &model.region(material)?.shells;
    if boundaries.len() != value.regions.len() - 1 {
        return Err(Error::invalid_argument("STEP void boundary count mismatch"));
    }
    let mut seen_shells = HashSet::new();
    let mut seen_faces = HashSet::new();
    let mut empty_regions = HashSet::new();
    for (index, boundary) in boundaries.iter().enumerate() {
        let shell = model.shell(*boundary)?;
        if !seen_shells.insert(*boundary)
            || shell.faces.is_empty()
            || !shell.wire_edges.is_empty()
            || !shell.acorn_vertices.is_empty()
        {
            return Err(Error::invalid_argument("invalid STEP material boundary"));
        }
        let mut front_shell = None;
        for face in &shell.faces {
            let face = model.face(*face)?;
            if face.back != Some(*boundary) {
                return Err(Error::unimplemented(
                    "STEP void export needs face normals pointing out of material",
                ));
            }
            let front = face
                .front
                .ok_or_else(|| Error::invalid_argument("STEP void face has no empty side"))?;
            if front_shell
                .replace(front)
                .is_some_and(|previous| previous != front)
            {
                return Err(Error::invalid_argument(
                    "STEP boundary has mixed empty sides",
                ));
            }
        }
        let front = front_shell.ok_or_else(|| Error::internal("empty STEP boundary"))?;
        let empty_shell = model.shell(front)?;
        let empty = model.region(empty_shell.region)?;
        if !seen_shells.insert(front)
            || !regions.contains(&empty_shell.region)
            || empty.solid
            || !empty_regions.insert(empty_shell.region)
            || (empty_shell.region == value.regions[0]) != (index == 0)
            || empty.shells.as_slice() != [front]
            || !empty_shell.wire_edges.is_empty()
            || !empty_shell.acorn_vertices.is_empty()
        {
            return Err(Error::invalid_argument(
                "invalid STEP empty boundary pairing",
            ));
        }
        let faces: HashSet<_> = shell.faces.iter().copied().collect();
        let front_faces: HashSet<_> = empty_shell.faces.iter().copied().collect();
        if faces.len() != shell.faces.len()
            || front_faces.len() != empty_shell.faces.len()
            || faces != front_faces
            || faces.iter().any(|face| !seen_faces.insert(*face))
        {
            return Err(Error::invalid_argument(
                "invalid STEP boundary face membership",
            ));
        }
    }
    voids::validate(model, body)?;
    Ok(boundaries.clone())
}

#[derive(Default)]
struct Writer {
    next: u32,
    lines: Vec<String>,
}

impl Writer {
    fn id(&mut self) -> u32 {
        self.next += 1;
        self.next
    }

    fn entity(&mut self, body: String) -> u32 {
        let id = self.id();
        self.lines.push(format!("#{id}={body};"));
        id
    }

    fn cartesian(&mut self, point: Point3) -> u32 {
        self.entity(format!(
            "CARTESIAN_POINT('',({x},{y},{z}))",
            x = format_number(point.x),
            y = format_number(point.y),
            z = format_number(point.z)
        ))
    }

    fn placement(&mut self, frame: &Frame) -> u32 {
        let origin = self.cartesian(frame.origin());
        let axis = self.direction(frame.z());
        let reference = self.direction(frame.x());
        self.entity(format!(
            "AXIS2_PLACEMENT_3D('',#{origin},#{axis},#{reference})"
        ))
    }

    fn direction(&mut self, vector: Vec3) -> u32 {
        self.entity(format!(
            "DIRECTION('',({x},{y},{z}))",
            x = format_number(vector.x),
            y = format_number(vector.y),
            z = format_number(vector.z)
        ))
    }

    fn finish(self) -> String {
        let mut text = String::from(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('Gensolid planar solid'),'2;1');\n",
        );
        text.push_str(
            "FILE_NAME('body.step','',('Gensolid'),('Gensolid'),'Gensolid','Gensolid','');\n",
        );
        text.push_str("FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\nENDSEC;\nDATA;\n");
        for line in self.lines {
            text.push_str(&line);
            text.push('\n');
        }
        text.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        text
    }
}

fn format_number(value: f64) -> String {
    // Rust's shortest round-trip representation retains small features and large offsets.
    let text = value.to_string();
    if text.contains('.') || text.contains(['e', 'E']) {
        text
    } else {
        format!("{text}.")
    }
}

/// Compatibility entry point for [`read_solids`], now also accepting analytic/NURBS faces.
pub fn read_planar_solids(text: &str, precision: Precision) -> Result<Model> {
    read_solids(text, precision)
}

/// Imports standalone manifold solids, preserving the file's numeric length unit.
///
/// `BREP_WITH_VOIDS` is accepted for connected convex planar shells with line edges
/// when interval halfspaces certify cavity containment and mutual separation.
/// Curved, nonconvex, touching or otherwise unproven cavity relationships are refused.
///
/// Each output body must pass `check_body`. Unsupported assembly transforms, body
/// classes, geometry, and incomplete checks return an error rather than a partial model.
pub fn read_solids(text: &str, precision: Precision) -> Result<Model> {
    let model = read_unchecked(text, precision)?;
    for body in model.bodies() {
        let report =
            gensolid_check::check_body(&model, body, &gensolid_check::CheckOptions::default())?;
        if report.unchecked != 0 {
            return Err(Error::unimplemented(format!(
                "STEP body has {} unchecked geometric queries",
                report.unchecked
            )));
        }
        if !report.faults.is_empty() {
            return Err(Error::invalid_argument(format!(
                "STEP body failed check_body: {:?}",
                report.faults
            )));
        }
        voids::validate(&model, body)?;
    }
    Ok(model)
}

fn read_unchecked(text: &str, precision: Precision) -> Result<Model> {
    let instances = parse_document(text)?;
    let mut ids: Vec<_> = instances.keys().copied().collect();
    ids.sort_unstable();
    let unsupported = [
        "MAPPED_ITEM",
        "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
        "ASSEMBLY_COMPONENT_USAGE",
        "SPECIFIED_HIGHER_USAGE_OCCURRENCE",
        "FACETED_BREP",
        "FACETED_BREP_AND_BREP_WITH_VOIDS",
        "SHELL_BASED_SURFACE_MODEL",
        "FACE_BASED_SURFACE_MODEL",
        "EDGE_BASED_WIREFRAME_MODEL",
        "CSG_SOLID",
        "BOOLEAN_RESULT",
        "BOOLEAN_CLIPPING_RESULT",
        "SOLID_REPLICA",
        "TESSELLATED_SOLID",
        "EXTRUDED_AREA_SOLID",
        "REVOLVED_AREA_SOLID",
        "SWEPT_DISK_SOLID",
    ];
    let mut solids = Vec::new();
    for id in ids {
        let entity = instance(&instances, id)?;
        for name in unsupported {
            if entity.type_name == name || entity.components.contains_key(name) {
                return Err(Error::unimplemented(format!(
                    "STEP #{id} {name} needs unsupported body or assembly semantics"
                )));
            }
        }
        if matches!(
            entity.type_name.as_str(),
            "MANIFOLD_SOLID_BREP" | "BREP_WITH_VOIDS"
        ) || entity.components.contains_key("MANIFOLD_SOLID_BREP")
        {
            solids.push(id);
        }
    }
    if solids.is_empty() {
        return Err(Error::unsupported("STEP file has no MANIFOLD_SOLID_BREP"));
    }
    let mut model = Model::new(precision);
    for solid in solids {
        read_solid(&mut model, &instances, solid)?;
    }
    Ok(model)
}

fn instance(instances: &HashMap<u32, Instance>, id: u32) -> Result<&Instance> {
    instances
        .get(&id)
        .ok_or_else(|| Error::parse(format!("STEP reference #{id} is missing")))
}

#[allow(clippy::needless_lifetimes)]
fn expect_type<'a>(
    instances: &'a HashMap<u32, Instance>,
    id: u32,
    type_name: &str,
) -> Result<&'a [Value]> {
    let found = instance(instances, id)?;
    if found.type_name != type_name {
        return Err(Error::unsupported(format!(
            "STEP #{} is {}, expected {type_name}",
            id, found.type_name
        )));
    }
    let args = name_skipped(&found.args)?;
    let fields = match type_name {
        "CARTESIAN_POINT"
        | "DIRECTION"
        | "MANIFOLD_SOLID_BREP"
        | "CLOSED_SHELL"
        | "EDGE_LOOP"
        | "VERTEX_POINT" => Some(1),
        "VECTOR" | "BREP_WITH_VOIDS" => Some(2),
        "AXIS2_PLACEMENT_3D" | "ADVANCED_FACE" | "ORIENTED_CLOSED_SHELL" => Some(3),
        "ORIENTED_EDGE" | "EDGE_CURVE" => Some(4),
        _ => None,
    };
    if let Some(fields) = fields
        && args.len() != fields
    {
        return Err(Error::parse(format!(
            "STEP #{id} {type_name} has {} fields, expected {fields}",
            args.len()
        )));
    }
    Ok(args)
}

fn cartesian(instances: &HashMap<u32, Instance>, id: u32) -> Result<Point3> {
    let args = expect_type(instances, id, "CARTESIAN_POINT")?;
    let coords = nums(
        args.first()
            .ok_or_else(|| Error::parse("CARTESIAN_POINT has no coordinates"))?,
    )?;
    if coords.len() != 3 {
        return Err(Error::unsupported(
            "only 3D CARTESIAN_POINT entities are supported",
        ));
    }
    let point = Point3::new(coords[0], coords[1], coords[2]);
    if !point.is_finite() {
        return Err(Error::parse("STEP point is not finite"));
    }
    Ok(point)
}

fn direction(instances: &HashMap<u32, Instance>, id: u32) -> Result<Vec3> {
    let args = expect_type(instances, id, "DIRECTION")?;
    let coords = nums(
        args.first()
            .ok_or_else(|| Error::parse("DIRECTION has no components"))?,
    )?;
    if coords.len() != 3 {
        return Err(Error::unsupported(
            "only 3D DIRECTION entities are supported",
        ));
    }
    Ok(Vec3::new(coords[0], coords[1], coords[2]))
}

fn placement(instances: &HashMap<u32, Instance>, id: u32, precision: &Precision) -> Result<Frame> {
    let args = expect_type(instances, id, "AXIS2_PLACEMENT_3D")?;
    if args.len() < 3 {
        return Err(Error::parse("AXIS2_PLACEMENT_3D is missing its axes"));
    }
    let origin = cartesian(instances, ref_of(&args[0])?)?;
    let axis = if matches!(args[1], Value::Unset) {
        Vec3::Z
    } else {
        direction(instances, ref_of(&args[1])?)?
    };
    if matches!(args[2], Value::Unset) {
        Frame::from_axis(origin, axis, precision)
    } else {
        Frame::new(
            origin,
            axis,
            direction(instances, ref_of(&args[2])?)?,
            precision,
        )
    }
}

fn read_solid(model: &mut Model, instances: &HashMap<u32, Instance>, id: u32) -> Result<BodyId> {
    let (shell_id, cavities) = solid_shells(instances, id)?;
    if cavities.is_empty() {
        // Standalone manifold imports retain their original closed-shell scope.
        // Oriented wrappers require the relationship certification used for voids.
        expect_type(instances, shell_id, "CLOSED_SHELL")?;
    }
    let (shell_id, reverse_outer) = closed_shell(instances, shell_id)?;
    let mut cavity_shells = Vec::new();
    let mut seen_shells = std::collections::HashSet::from([shell_id]);
    for cavity in cavities {
        if instance(instances, cavity)?.type_name != "ORIENTED_CLOSED_SHELL" {
            return Err(Error::parse(format!(
                "STEP #{id} void #{cavity} is not an ORIENTED_CLOSED_SHELL"
            )));
        }
        let (shell, reversed) = closed_shell(instances, cavity)?;
        if !seen_shells.insert(shell) {
            return Err(Error::parse(format!("STEP #{id} repeats shell #{shell}")));
        }
        cavity_shells.push((shell, reversed));
    }
    let precision = *model.precision();
    let (body, infinite, solid) = {
        let mut raw = model.raw();
        let body = raw.insert_body(Body {
            kind: BodyKind::Solid,
            regions: Vec::new(),
        })?;
        let infinite = raw.insert_region(Region {
            body,
            solid: false,
            shells: Vec::new(),
        })?;
        let solid = raw.insert_region(Region {
            body,
            solid: true,
            shells: Vec::new(),
        })?;
        raw.body_mut(body)?.regions = vec![infinite, solid];
        (body, infinite, solid)
    };
    let mut vertex_ids = HashMap::new();
    let mut edge_ids = HashMap::new();
    for (shell, reversed, cavity) in std::iter::once((shell_id, reverse_outer, false)).chain(
        cavity_shells
            .into_iter()
            .map(|(shell, reversed)| (shell, reversed, true)),
    ) {
        let empty = if cavity {
            let mut raw = model.raw();
            let region = raw.insert_region(Region {
                body,
                solid: false,
                shells: Vec::new(),
            })?;
            raw.body_mut(body)?.regions.push(region);
            region
        } else {
            infinite
        };
        let (front, back) = insert_shell_pair(model, empty, solid)?;
        for face in shell_faces(instances, shell)? {
            read_face(
                model,
                instances,
                face,
                front,
                back,
                &precision,
                &mut vertex_ids,
                &mut edge_ids,
                reversed,
            )?;
        }
    }
    Ok(body)
}

fn solid_shells(instances: &HashMap<u32, Instance>, id: u32) -> Result<(u32, Vec<u32>)> {
    let found = instance(instances, id)?;
    let (outer, cavities) = if found.components.is_empty() {
        let args = expect_type(instances, id, &found.type_name)?;
        let outer = args
            .first()
            .ok_or_else(|| Error::parse("STEP solid has no outer shell"))?;
        (
            outer,
            if found.type_name == "BREP_WITH_VOIDS" {
                args.get(1)
            } else {
                None
            },
        )
    } else {
        for component in found.components.keys() {
            if !matches!(
                component.as_str(),
                "MANIFOLD_SOLID_BREP"
                    | "BREP_WITH_VOIDS"
                    | "REPRESENTATION_ITEM"
                    | "GEOMETRIC_REPRESENTATION_ITEM"
                    | "SOLID_MODEL"
            ) {
                return Err(Error::unimplemented(format!(
                    "STEP solid #{id} complex component {component} is not implemented"
                )));
            }
            if matches!(
                component.as_str(),
                "GEOMETRIC_REPRESENTATION_ITEM" | "SOLID_MODEL"
            ) && found
                .components
                .get(component)
                .is_some_and(|args| !args.is_empty())
            {
                return Err(Error::parse(format!(
                    "complex STEP solid {component} must have no local fields"
                )));
            }
        }
        let name = found
            .components
            .get("REPRESENTATION_ITEM")
            .ok_or_else(|| Error::parse("complex STEP solid is missing its representation name"))?;
        if !matches!(name.as_slice(), [Value::Str(_)]) {
            return Err(Error::parse("invalid complex STEP solid name"));
        }
        let args = found
            .components
            .get("MANIFOLD_SOLID_BREP")
            .ok_or_else(|| Error::parse("complex STEP solid has no outer shell"))?;
        let [outer] = args.as_slice() else {
            return Err(Error::parse("invalid complex MANIFOLD_SOLID_BREP fields"));
        };
        let cavities = found
            .components
            .get("BREP_WITH_VOIDS")
            .map(|args| match args.as_slice() {
                [value] => Ok(value),
                _ => Err(Error::parse("invalid complex BREP_WITH_VOIDS fields")),
            })
            .transpose()?;
        (outer, cavities)
    };
    let cavities = if let Some(Value::List(values)) = cavities {
        if values.is_empty() {
            return Err(Error::parse("BREP_WITH_VOIDS requires at least one cavity"));
        }
        values.iter().map(ref_of).collect::<Result<Vec<_>>>()?
    } else if cavities.is_some() {
        return Err(Error::parse(
            "BREP_WITH_VOIDS voids must be a reference list",
        ));
    } else {
        Vec::new()
    };
    Ok((ref_of(outer)?, cavities))
}

fn closed_shell(instances: &HashMap<u32, Instance>, id: u32) -> Result<(u32, bool)> {
    let found = instance(instances, id)?;
    if found.type_name == "CLOSED_SHELL" {
        return Ok((id, false));
    }
    let args = expect_type(instances, id, "ORIENTED_CLOSED_SHELL")?;
    if !matches!(args[0], Value::Star) {
        return Err(Error::parse(
            "ORIENTED_CLOSED_SHELL faces must be derived '*'",
        ));
    }
    let shell = ref_of(&args[1])?;
    if instance(instances, shell)?.type_name == "ORIENTED_CLOSED_SHELL" {
        return Err(Error::parse(
            "ORIENTED_CLOSED_SHELL must not reference another oriented shell",
        ));
    }
    expect_type(instances, shell, "CLOSED_SHELL")?;
    Ok((shell, !bool_of(&args[2])?))
}

fn shell_faces(instances: &HashMap<u32, Instance>, shell_id: u32) -> Result<Vec<u32>> {
    let shell_args = expect_type(instances, shell_id, "CLOSED_SHELL")?;
    let face_refs = refs(
        shell_args
            .first()
            .ok_or_else(|| Error::parse("shell has no faces"))?,
    )?;
    if face_refs.is_empty() {
        return Err(Error::invalid_argument(format!(
            "STEP #{shell_id} CLOSED_SHELL has no faces"
        )));
    }
    let mut seen_faces = std::collections::HashSet::new();
    if face_refs.iter().any(|id| !seen_faces.insert(*id)) {
        return Err(Error::parse(format!(
            "STEP #{shell_id} CLOSED_SHELL repeats a face"
        )));
    }
    Ok(face_refs)
}

fn insert_shell_pair(
    model: &mut Model,
    empty: gensolid_topo::RegionId,
    solid: gensolid_topo::RegionId,
) -> Result<(gensolid_topo::ShellId, gensolid_topo::ShellId)> {
    let mut raw = model.raw();
    let outer = raw.insert_shell(Shell {
        region: empty,
        faces: Vec::new(),
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    })?;
    let inner = raw.insert_shell(Shell {
        region: solid,
        faces: Vec::new(),
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    })?;
    raw.region_mut(empty)?.shells.push(outer);
    raw.region_mut(solid)?.shells.push(inner);
    Ok((outer, inner))
}

#[allow(clippy::too_many_arguments)]
fn read_face(
    model: &mut Model,
    instances: &HashMap<u32, Instance>,
    id: u32,
    outer: gensolid_topo::ShellId,
    inner: gensolid_topo::ShellId,
    precision: &Precision,
    vertex_ids: &mut HashMap<u32, gensolid_topo::VertexId>,
    edge_ids: &mut HashMap<u32, gensolid_topo::EdgeId>,
    reversed: bool,
) -> Result<()> {
    let args = expect_type(instances, id, "ADVANCED_FACE")?;
    if args.len() < 3 {
        return Err(Error::parse("ADVANCED_FACE is missing its surface"));
    }
    let bounds = refs(&args[0])?;
    if bounds.is_empty() {
        return Err(Error::unimplemented(
            "unbounded STEP faces are not implemented",
        ));
    }
    let same_sense = bool_of(&args[2])? ^ reversed;
    let surface_id = ref_of(&args[1])?;
    let surface = geometry::read_surface(instances, surface_id, precision)?;
    let surface = model.raw().geom_mut().add_surface(surface)?;
    let face = model.raw().insert_face(Face {
        front: Some(outer),
        back: Some(inner),
        loops: Vec::new(),
        surface: Some(surface),
        sense: if same_sense {
            Sense::Positive
        } else {
            Sense::Negative
        },
    })?;
    for bound in bounds {
        read_bound(
            model, instances, bound, face, precision, vertex_ids, edge_ids, reversed,
        )?;
    }
    model.raw().shell_mut(outer)?.faces.push(face);
    model.raw().shell_mut(inner)?.faces.push(face);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn read_bound(
    model: &mut Model,
    instances: &HashMap<u32, Instance>,
    id: u32,
    face: gensolid_topo::FaceId,
    precision: &Precision,
    vertex_ids: &mut HashMap<u32, gensolid_topo::VertexId>,
    edge_ids: &mut HashMap<u32, gensolid_topo::EdgeId>,
    reversed: bool,
) -> Result<()> {
    let found = instance(instances, id)?;
    if found.type_name != "FACE_OUTER_BOUND" && found.type_name != "FACE_BOUND" {
        return Err(Error::unsupported(format!(
            "STEP bound {} is not implemented",
            found.type_name
        )));
    }
    let args = name_skipped(&found.args)?;
    if args.len() != 2 {
        return Err(Error::parse("face bound is missing its loop"));
    }
    let orientation = bool_of(&args[1])? ^ reversed;
    let loop_args = expect_type(instances, ref_of(&args[0])?, "EDGE_LOOP")?;
    let mut oriented = refs(
        loop_args
            .first()
            .ok_or_else(|| Error::parse("EDGE_LOOP is empty"))?,
    )?;
    if !orientation {
        oriented.reverse();
    }
    let loop_id = model.raw().insert_loop(Loop {
        face,
        fins: Vec::new(),
    })?;
    for oriented_id in oriented {
        read_oriented_edge(
            model,
            instances,
            oriented_id,
            loop_id,
            precision,
            vertex_ids,
            edge_ids,
            !orientation,
        )?;
    }
    model.raw().face_mut(face)?.loops.push(loop_id);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn read_oriented_edge(
    model: &mut Model,
    instances: &HashMap<u32, Instance>,
    id: u32,
    loop_id: gensolid_topo::LoopId,
    precision: &Precision,
    vertex_ids: &mut HashMap<u32, gensolid_topo::VertexId>,
    edge_ids: &mut HashMap<u32, gensolid_topo::EdgeId>,
    flip: bool,
) -> Result<()> {
    let args = expect_type(instances, id, "ORIENTED_EDGE")?;
    if args.len() < 4 {
        return Err(Error::parse("ORIENTED_EDGE is missing its edge"));
    }
    if !matches!((&args[0], &args[1]), (Value::Star, Value::Star)) {
        return Err(Error::parse(format!(
            "STEP #{id} ORIENTED_EDGE endpoints must be derived '*' values"
        )));
    }
    let mut sense = if bool_of(&args[3])? {
        Sense::Positive
    } else {
        Sense::Negative
    };
    if flip {
        sense = sense.reversed();
    }
    let edge_args = expect_type(instances, ref_of(&args[2])?, "EDGE_CURVE")?;
    if !bool_of(
        edge_args
            .get(3)
            .ok_or_else(|| Error::parse("EDGE_CURVE has no same_sense"))?,
    )? {
        sense = sense.reversed();
    }
    let edge = read_edge_curve(
        model,
        instances,
        ref_of(&args[2])?,
        precision,
        vertex_ids,
        edge_ids,
    )?;
    let fin = model.raw().insert_fin(Fin {
        owner: loop_id,
        edge,
        sense,
    })?;
    model.raw().edge_mut(edge)?.fins.push(fin);
    model.raw().loop_mut(loop_id)?.fins.push(fin);
    Ok(())
}

fn read_edge_curve(
    model: &mut Model,
    instances: &HashMap<u32, Instance>,
    id: u32,
    precision: &Precision,
    vertex_ids: &mut HashMap<u32, gensolid_topo::VertexId>,
    edge_ids: &mut HashMap<u32, gensolid_topo::EdgeId>,
) -> Result<gensolid_topo::EdgeId> {
    if let Some(existing) = edge_ids.get(&id) {
        return Ok(*existing);
    }
    let args = expect_type(instances, id, "EDGE_CURVE")?;
    if args.len() < 4 {
        return Err(Error::parse("EDGE_CURVE is missing its curve"));
    }
    let mut start = read_vertex(model, instances, ref_of(&args[0])?, vertex_ids)?;
    let mut end = read_vertex(model, instances, ref_of(&args[1])?, vertex_ids)?;
    if !bool_of(&args[3])? {
        std::mem::swap(&mut start, &mut end);
    }
    let curve_id = ref_of(&args[2])?;
    let curve = geometry::read_curve(instances, curve_id, precision)?;
    let curve = model.raw().geom_mut().add_curve(curve)?;
    let edge = model.raw().insert_edge(Edge {
        fins: Vec::new(),
        start: Some(start),
        end: Some(end),
        curve: Some(curve),
        tolerance: None,
    })?;
    edge_ids.insert(id, edge);
    Ok(edge)
}

fn read_vertex(
    model: &mut Model,
    instances: &HashMap<u32, Instance>,
    id: u32,
    vertex_ids: &mut HashMap<u32, gensolid_topo::VertexId>,
) -> Result<gensolid_topo::VertexId> {
    if let Some(existing) = vertex_ids.get(&id) {
        return Ok(*existing);
    }
    let args = expect_type(instances, id, "VERTEX_POINT")?;
    let point = cartesian(
        instances,
        ref_of(
            args.first()
                .ok_or_else(|| Error::parse("VERTEX_POINT has no point"))?,
        )?,
    )?;
    let vertex = model.raw().insert_vertex(Vertex {
        point,
        tolerance: None,
    })?;
    vertex_ids.insert(id, vertex);
    Ok(vertex)
}

/// Builds a block and writes it. Used by the verifier to hand the same solid to OCCT.
pub fn write_block_step(
    x: f64,
    y: f64,
    z: f64,
    frame: &Frame,
    precision: Precision,
) -> Result<String> {
    let mut model = Model::new(precision);
    let body = make_block(&mut model, x, y, z, frame)?;
    write_planar_solid(&model, body)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::float_cmp)]
    fn parses_a_point_entity() {
        let text = "ISO-10303-21;DATA;\n#1=CARTESIAN_POINT('',(1.5,2.,-3.));\nENDSEC;\n";
        let data = parse::parse_data(text).unwrap();
        let point = cartesian(&data, 1).unwrap();
        assert_eq!(point, Point3::new(1.5, 2.0, -3.0));
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn block_round_trips_and_passes_check() {
        use gensolid_check::{CheckOptions, check_body};

        let precision = Precision::DEFAULT;
        let step = write_block_step(2.0, 3.0, 4.0, &Frame::WORLD, precision).unwrap();
        let model = read_planar_solids(&step, precision).unwrap();
        let body = model.bodies().next().unwrap();
        let report = check_body(&model, body, &CheckOptions::default()).unwrap();
        assert!(report.is_valid(), "{report:?}");
        assert_eq!(model.body_faces(body).unwrap().len(), 6);
        assert_eq!(model.body_edges(body).unwrap().len(), 12);
        assert_eq!(model.body_vertices(body).unwrap().len(), 8);
        let volume = gensolid_topo::body_volume(&model, body).unwrap();
        assert_eq!(volume, 24.0);
    }
}

#[cfg(test)]
mod regressions;
#[cfg(test)]
mod void_tests;
