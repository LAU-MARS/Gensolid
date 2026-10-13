//! Schema-aware XT import into the kernel model.
//!
//! Format authority: Siemens, Parasolid XT Format Reference (April 2008),
//! <https://www.plm.automation.siemens.com/en_us/Images/XT_Format_April_2008_tcm1023-62642.pdf>.
//! Public official alternative: Siemens JT v10 Rev-B, Appendix E,
//! <https://www.plm.automation.siemens.com/en_us/Images/JT-v10-file-format-reference-rev-B_tcm1023-233786.pdf>.
//! Raw decoding uses the separately distributed `parasolid-core` 0.3.7 dependency
//! (MIT AND Apache-2.0, <https://github.com/monozukuri-ai/parasolid-kit>).
//! This adapter is original Gensolid code; no OCCT code or manual prose is included.

use std::collections::{BTreeMap, HashSet};

use gensolid_check::{CheckOptions, Report, check_body};
use gensolid_core::{Error, ErrorCode, Precision, Result, Tolerance};
use gensolid_geom::{
    BCurve, BSurface, Circle, Cone, Curve, Cylinder, Ellipse, Line, Plane, Sphere, Surface, Torus,
};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{
    Body, BodyId, BodyKind, Edge, Face, Fin, Loop, Model, Region, Sense, Shell, Vertex,
};
use parasolid_core::brep::{self as source, BrepId, BrepModel};
use parasolid_core::{
    BuiltinProfileRegistry, DocumentLimits, ErrorKind, InspectionLimits, ParseError, SchemaKey,
    SchemaProvider,
};

/// A successfully imported model. Every body has passed `check_body` without skips.
pub struct XtImport {
    /// Kernel topology and geometry, in the units of the source file.
    pub model: Model,
    /// Imported bodies, in source order.
    pub bodies: Vec<BodyId>,
    /// Exact internal schema identifier used to decode the file.
    pub schema: String,
}

pub(crate) struct Decoded {
    pub model: Model,
    pub bodies: Vec<BodyId>,
    pub reports: Vec<Report>,
}

pub(crate) struct Header {
    pub schema: String,
    pub modeller: String,
    pub binary: bool,
    pub data_start: usize,
}

/// Inspect either text or neutral-binary XT without guessing the schema.
pub(crate) fn inspect(input: &[u8]) -> Result<Header> {
    // The binary reader recognizes both bare PS streams and the common text header.
    match parasolid_core::inspect_xb(input, InspectionLimits::default()) {
        Ok(header) => {
            validate_common_header(input, header.text_header_range)?;
            Ok(Header {
                schema: header.schema_key,
                modeller: header.modeller_version,
                binary: true,
                data_start: header.header_range.end,
            })
        }
        Err(binary_error) => match parasolid_core::inspect_xt(input, InspectionLimits::default()) {
            Ok(header) => {
                validate_common_header(input, header.common_header_range)?;
                Ok(Header {
                    schema: header.schema_key,
                    modeller: header.modeller_version,
                    binary: false,
                    data_start: header.header_range.end,
                })
            }
            Err(text_error) => {
                if binary_error.kind() == ErrorKind::InvalidSignature {
                    Err(parser_error(text_error))
                } else {
                    Err(parser_error(binary_error))
                }
            }
        },
    }
}

fn validate_common_header(input: &[u8], range: Option<std::ops::Range<usize>>) -> Result<()> {
    if let Some(range) = range {
        let text = std::str::from_utf8(&input[range])
            .map_err(|_| Error::parse("XT common header is not UTF-8/ASCII text"))?;
        crate::header::read_text(text)?;
    }
    Ok(())
}

pub(crate) fn parser_error(error: ParseError) -> Error {
    let code = match error.kind() {
        ErrorKind::MissingBaseSchema
        | ErrorKind::MissingSchemaType
        | ErrorKind::UnsupportedBaseSchemaType
        | ErrorKind::UnknownBaseSchemaType
        | ErrorKind::UnsupportedBuiltinSchemaKey
        | ErrorKind::BuiltinProfileUncoveredType
        | ErrorKind::UnsupportedSchemaFieldType
        | ErrorKind::UnsupportedBinaryFormat
        | ErrorKind::UnsupportedUserFields => ErrorCode::Unimplemented,
        ErrorKind::LimitExceeded => ErrorCode::CapacityExceeded,
        ErrorKind::MissingBrepBody
        | ErrorKind::InvalidBrepField
        | ErrorKind::InvalidBrepReference
        | ErrorKind::InvalidBrepTopology
        | ErrorKind::InvalidGeometryParameter => ErrorCode::InvalidArgument,
        _ => ErrorCode::Parse,
    };
    Error::new(
        code,
        format!(
            "{} at byte {}: {error}",
            error.kind().code(),
            error.offset()
        ),
    )
}

pub(crate) fn source_with_provider<P: SchemaProvider>(
    input: &[u8],
    header: &Header,
    provider: &P,
) -> Result<BrepModel> {
    if header.binary {
        let document = parasolid_core::parse_xb(input, provider, DocumentLimits::default())
            .map_err(parser_error)?;
        source::map_xb_brep(&document).map_err(parser_error)
    } else {
        let document = parasolid_core::parse_xt(input, provider, DocumentLimits::default())
            .map_err(parser_error)?;
        source::map_xt_brep(&document).map_err(parser_error)
    }
}

pub(crate) fn source_builtin(input: &[u8], header: &Header) -> Result<BrepModel> {
    if header.schema == crate::schema::ABC_KEY {
        let provider = crate::schema::AbcSchema::new()?;
        if header.binary {
            let mut document =
                parasolid_core::parse_xb(input, &provider, DocumentLimits::default())
                    .map_err(parser_error)?;
            provider.audit(&document.schemas)?;
            provider.audit_nodes(&document.nodes)?;
            for node in &mut document.nodes {
                if node.definition.source == parasolid_core::SchemaSource::EmbeddedFull {
                    // A newly declared metadata node cannot acquire a base B-rep role
                    // merely by using a name such as BODY or FACE in the input.
                    node.definition.name = format!("XT_EXTENSION_{}", node.node_type);
                }
            }
            return source::map_xb_brep(&document).map_err(parser_error);
        }
        let mut document = parasolid_core::parse_xt(input, &provider, DocumentLimits::default())
            .map_err(parser_error)?;
        provider.audit(&document.schemas)?;
        provider.audit_nodes(&document.nodes)?;
        for node in &mut document.nodes {
            if node.definition.source == parasolid_core::SchemaSource::EmbeddedFull {
                node.definition.name = format!("XT_EXTENSION_{}", node.node_type);
            }
        }
        return source::map_xt_brep(&document).map_err(parser_error);
    }
    let key = SchemaKey::parse(&header.schema).map_err(parser_error)?;
    let registry = BuiltinProfileRegistry::compiled().map_err(parser_error)?;
    let provider = registry.provider_for_key(&key).ok_or_else(|| Error::unimplemented(format!(
        "schema {} has no registered decoder; no body was built or checked; supply an exact schema provider", header.schema
    )))?;
    source_with_provider(input, header, &provider)
}

/// Imports text (`.x_t`) or neutral binary (`.x_b`) with a reviewed exact built-in schema.
///
/// Unsupported schemas and unrepresentable geometry return `Unimplemented`.
/// Invalid topology or a failed/unfinished kernel check never returns a model.
/// No units are guessed or converted; `precision` is expressed in source units.
pub fn read_bytes(input: &[u8], precision: Precision) -> Result<XtImport> {
    let header = inspect(input)?;
    let source = source_builtin(input, &header)?;
    finish(source, precision)
}

/// Imports using caller-supplied exact schema definitions instead of the built-in subset.
///
/// The provider must supply the schema named by the internal stream, including its base
/// schema for embedded definitions. There is no fallback to a nearby version.
pub fn read_bytes_with_schema<P: SchemaProvider>(
    input: &[u8],
    precision: Precision,
    provider: &P,
) -> Result<XtImport> {
    let header = inspect(input)?;
    let source = source_with_provider(input, &header, provider)?;
    finish(source, precision)
}

fn finish(source: BrepModel, precision: Precision) -> Result<XtImport> {
    let decoded = map_model(&source, precision)?;
    for (index, report) in decoded.reports.iter().enumerate() {
        if !report.is_valid() {
            let message = format!(
                "body {index}: {} faults, {} unchecked geometric checks: {:?}",
                report.faults.len(),
                report.unchecked,
                report.faults
            );
            return Err(if report.unchecked != 0 && report.faults.is_empty() {
                Error::unimplemented(message)
            } else {
                Error::parse(message)
            });
        }
    }
    Ok(XtImport {
        model: decoded.model,
        bodies: decoded.bodies,
        schema: source.schema_key,
    })
}

fn lookup<T: Copy>(map: &BTreeMap<BrepId, T>, id: BrepId, what: &str) -> Result<T> {
    map.get(&id)
        .copied()
        .ok_or_else(|| Error::parse(format!("unmapped {what} source id {id}")))
}

fn ids<T: Copy>(map: &BTreeMap<BrepId, T>, values: &[BrepId], what: &str) -> Result<Vec<T>> {
    values.iter().map(|&id| lookup(map, id, what)).collect()
}

fn optional<T: Copy>(
    map: &BTreeMap<BrepId, T>,
    id: Option<BrepId>,
    what: &str,
) -> Result<Option<T>> {
    id.map(|id| lookup(map, id, what)).transpose()
}

fn sense(value: source::Sense) -> Result<Sense> {
    match value {
        source::Sense::Positive => Ok(Sense::Positive),
        source::Sense::Negative => Ok(Sense::Negative),
        source::Sense::Unknown => Err(Error::unimplemented(
            "unknown XT orientation cannot be represented",
        )),
    }
}

fn combine(a: Sense, b: Sense) -> Sense {
    if a == b {
        Sense::Positive
    } else {
        Sense::Negative
    }
}

fn point(v: source::Vector3) -> Point3 {
    Point3::new(v.x, v.y, v.z)
}
fn vector(v: source::Vector3) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}
fn frame(
    origin: source::Vector3,
    axis: source::Vector3,
    reference: source::Vector3,
    precision: &Precision,
) -> Result<Frame> {
    Frame::new(point(origin), vector(axis), vector(reference), precision)
}

fn tolerance(value: Option<f64>, precision: &Precision) -> Result<Option<Tolerance>> {
    // Smaller declared radii remain exact at this session's resolution; this is a
    // representational floor, never a geometric branch-selection tolerance.
    value
        .map(|v| {
            if !v.is_finite() || v.is_sign_negative() {
                return Err(Error::parse("invalid XT entity tolerance"));
            }
            Tolerance::new(v.max(precision.default_local()), precision)
        })
        .transpose()
}

fn expand_knots(
    knots: &[f64],
    multiplicities: &[u16],
    count: usize,
    degree: usize,
) -> Result<Vec<f64>> {
    if knots.len() != multiplicities.len() {
        return Err(Error::parse("XT knot multiplicities do not match knots"));
    }
    let expected = count
        .checked_add(degree)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| {
            Error::new(
                ErrorCode::CapacityExceeded,
                "XT expanded knot count overflow",
            )
        })?;
    let actual = multiplicities
        .iter()
        .try_fold(0usize, |sum, &value| sum.checked_add(usize::from(value)))
        .ok_or_else(|| {
            Error::new(
                ErrorCode::CapacityExceeded,
                "XT knot multiplicity sum overflow",
            )
        })?;
    if actual != expected {
        return Err(Error::parse(
            "XT expanded knot count does not match degree and control count",
        ));
    }
    let mut out = Vec::new();
    out.try_reserve_exact(expected)
        .map_err(|_| Error::new(ErrorCode::CapacityExceeded, "XT knot allocation failed"))?;
    for (&knot, &count) in knots.iter().zip(multiplicities) {
        out.extend(std::iter::repeat_n(knot, usize::from(count)));
    }
    Ok(out)
}

fn poles(vertices: &[Vec<f64>], rational: bool) -> Result<(Vec<Point3>, Option<Vec<f64>>)> {
    let mut points = Vec::with_capacity(vertices.len());
    let mut weights = rational.then(Vec::new);
    for vertex in vertices {
        match (rational, vertex.as_slice()) {
            (false, &[x, y, z]) => points.push(Point3::new(x, y, z)),
            (true, &[x, y, z, w]) if w.is_finite() && w > 0.0 => {
                points.push(Point3::new(x / w, y / w, z / w));
                if let Some(weights) = &mut weights {
                    weights.push(w);
                }
            }
            _ => {
                return Err(Error::unimplemented(
                    "XT spline requires finite 3D poles with positive weights; parameter-space splines are not kernel curves",
                ));
            }
        }
    }
    Ok((points, weights))
}

fn curve(source: &source::CurveGeometry, precision: &Precision) -> Result<Curve> {
    use source::CurveKind as K;
    match &source.kind {
        K::Line {
            point: p,
            direction,
        } => Ok(Curve::Line(Line::new(
            point(*p),
            vector(*direction),
            precision,
        )?)),
        K::Circle {
            center,
            normal,
            x_axis,
            radius,
        } => Ok(Curve::Circle(Circle::new(
            frame(*center, *normal, *x_axis, precision)?,
            *radius,
            precision,
        )?)),
        K::Ellipse {
            center,
            normal,
            x_axis,
            major_radius,
            minor_radius,
        } => Ok(Curve::Ellipse(Ellipse::new(
            frame(*center, *normal, *x_axis, precision)?,
            *major_radius,
            *minor_radius,
            precision,
        )?)),
        K::Nurbs(n) => {
            if n.periodic {
                return Err(Error::unimplemented(
                    "periodic XT spline parameter wrapping is not represented by BCurve",
                ));
            }
            let (points, weights) = poles(&n.control_vertices, n.rational)?;
            let degree = usize::from(n.degree);
            Ok(Curve::BCurve(BCurve::new(
                degree,
                expand_knots(&n.knots, &n.knot_multiplicities, points.len(), degree)?,
                points,
                weights,
            )?))
        }
        _ => Err(Error::unimplemented(format!(
            "XT curve node {} ({}) is not representable by the kernel",
            source.source.node_index,
            source.kind.as_str()
        ))),
    }
}

fn surface(source: &source::SurfaceGeometry, precision: &Precision) -> Result<Surface> {
    use source::SurfaceKind as K;
    match &source.kind {
        K::Plane {
            point,
            normal,
            x_axis,
        } => Ok(Surface::Plane(Plane::new(frame(
            *point, *normal, *x_axis, precision,
        )?))),
        K::Cylinder {
            point,
            axis,
            radius,
            x_axis,
        } => Ok(Surface::Cylinder(Cylinder::new(
            frame(*point, *axis, *x_axis, precision)?,
            *radius,
            precision,
        )?)),
        K::Cone {
            point,
            axis,
            radius,
            sin_half_angle,
            cos_half_angle,
            x_axis,
        } => {
            // A negative cone slope is represented by reversing the axial direction.
            // Both parameterizations retain the outward natural normal.
            let axis = if sin_half_angle.is_sign_negative() {
                source::Vector3::from_array([-axis.x, -axis.y, -axis.z])
            } else {
                *axis
            };
            Ok(Surface::Cone(Cone::new(
                frame(*point, axis, *x_axis, precision)?,
                *radius,
                sin_half_angle.abs().atan2(*cos_half_angle),
                precision,
            )?))
        }
        K::Sphere {
            center,
            radius,
            axis,
            x_axis,
        } => Ok(Surface::Sphere(Sphere::new(
            frame(*center, *axis, *x_axis, precision)?,
            *radius,
            precision,
        )?)),
        K::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
            x_axis,
        } => Ok(Surface::Torus(Torus::new(
            frame(*center, *axis, *x_axis, precision)?,
            *major_radius,
            *minor_radius,
            precision,
        )?)),
        K::Nurbs(n) => {
            if n.u_periodic || n.v_periodic {
                return Err(Error::unimplemented(
                    "periodic XT surface parameter wrapping is not represented by BSurface",
                ));
            }
            let (points, weights) = poles(&n.control_vertices, n.rational)?;
            let (u_degree, v_degree) = (usize::from(n.u_degree), usize::from(n.v_degree));
            let u_count = usize::try_from(n.u_control_vertex_count).map_err(|_| {
                Error::new(ErrorCode::CapacityExceeded, "XT U control count overflow")
            })?;
            let v_count = usize::try_from(n.v_control_vertex_count).map_err(|_| {
                Error::new(ErrorCode::CapacityExceeded, "XT V control count overflow")
            })?;
            Ok(Surface::BSurface(BSurface::new(
                u_degree,
                v_degree,
                expand_knots(&n.u_knots, &n.u_knot_multiplicities, u_count, u_degree)?,
                expand_knots(&n.v_knots, &n.v_knot_multiplicities, v_count, v_degree)?,
                u_count,
                v_count,
                points,
                weights,
            )?))
        }
        _ => Err(Error::unimplemented(format!(
            "XT surface node {} ({}) is not representable by the kernel",
            source.source.node_index,
            source.kind.as_str()
        ))),
    }
}

pub(crate) fn map_model(source: &BrepModel, precision: Precision) -> Result<Decoded> {
    if !source.complete {
        return Err(Error::unimplemented(format!(
            "incomplete XT source mapping: {:?}",
            source.diagnostics
        )));
    }
    if source.bodies.is_empty() {
        return Err(Error::parse("XT contains no bodies"));
    }
    if source.bodies.iter().any(|body| {
        !(body.size_resolution.is_finite()
            && body.size_resolution > 0.0
            && body.linear_resolution.is_finite()
            && body.linear_resolution > 0.0)
    }) {
        return Err(Error::parse(
            "XT body declares non-positive or non-finite modeller resolution",
        ));
    }
    // The immutable checker visits geometry through face loops. Wire/acorn
    // geometry is not currently visited, so an empty check is not evidence of
    // validity for these parts, even if their source topology is consistent.
    if source
        .bodies
        .iter()
        .any(|body| body.kind == source::BodyKind::Wire)
        || source
            .shells
            .iter()
            .any(|shell| !shell.wire_edges.is_empty() || shell.isolated_vertex.is_some())
    {
        return Err(Error::unimplemented(
            "XT wire/acorn geometry is not covered by the kernel body checker",
        ));
    }
    // Gensolid has no per-fin carrier/pcurve. Accepting these would silently erase
    // tolerant topology even if the coarse carrier happened to pass sample checks.
    if source.half_edges.iter().any(|fin| fin.curve.is_some()) {
        return Err(Error::unimplemented(
            "XT per-fin curves require a tolerant fin geometry representation",
        ));
    }
    let mut model = Model::new(precision);
    let mut raw = model.raw();
    let mut curves = BTreeMap::new();
    let mut curve_senses = BTreeMap::new();
    for geometry in &source.curves {
        let id = raw.geom_mut().add_curve(curve(geometry, &precision)?)?;
        curves.insert(geometry.id, id);
        curve_senses.insert(geometry.id, sense(geometry.sense)?);
    }
    let mut surfaces = BTreeMap::new();
    let mut surface_senses = BTreeMap::new();
    for geometry in &source.surfaces {
        let id = raw.geom_mut().add_surface(surface(geometry, &precision)?)?;
        surfaces.insert(geometry.id, id);
        surface_senses.insert(geometry.id, sense(geometry.sense)?);
    }
    let points: BTreeMap<_, _> = source
        .points
        .iter()
        .map(|p| (p.id, point(p.position)))
        .collect();
    let mut vertices = BTreeMap::new();
    for vertex in &source.vertices {
        let id = raw.insert_vertex(Vertex {
            point: lookup(&points, vertex.point, "point")?,
            tolerance: tolerance(vertex.tolerance, &precision)?,
        })?;
        vertices.insert(vertex.id, id);
    }
    let mut bodies = BTreeMap::new();
    for body in &source.bodies {
        let kind = match body.kind {
            source::BodyKind::Solid => BodyKind::Solid,
            source::BodyKind::Wire => BodyKind::Wire,
            source::BodyKind::Sheet => BodyKind::Sheet,
            source::BodyKind::General => BodyKind::General,
        };
        bodies.insert(
            body.id,
            raw.insert_body(Body {
                kind,
                regions: Vec::new(),
            })?,
        );
    }
    let mut regions = BTreeMap::new();
    for region in &source.regions {
        regions.insert(
            region.id,
            raw.insert_region(Region {
                body: lookup(&bodies, region.body, "body")?,
                solid: region.kind == source::RegionKind::Solid,
                shells: Vec::new(),
            })?,
        );
    }
    let mut shells = BTreeMap::new();
    for shell in &source.shells {
        shells.insert(
            shell.id,
            raw.insert_shell(Shell {
                region: lookup(&regions, shell.region, "region")?,
                faces: Vec::new(),
                wire_edges: Vec::new(),
                acorn_vertices: optional(&vertices, shell.isolated_vertex, "vertex")?
                    .into_iter()
                    .collect(),
            })?,
        );
    }
    let mut faces = BTreeMap::new();
    for face in &source.faces {
        let surface = face
            .surface
            .ok_or_else(|| Error::parse(format!("face {} has no surface", face.id)))?;
        faces.insert(
            face.id,
            raw.insert_face(Face {
                // These are document-local BRep IDs: zero is a valid ID, not an XT null pointer.
                front: Some(lookup(&shells, face.front_shell, "front shell")?),
                back: Some(lookup(&shells, face.back_shell, "back shell")?),
                loops: Vec::new(),
                surface: Some(lookup(&surfaces, surface, "surface")?),
                sense: combine(
                    sense(face.sense)?,
                    lookup(&surface_senses, surface, "surface sense")?,
                ),
            })?,
        );
    }
    let mut loops = BTreeMap::new();
    for lp in &source.loops {
        loops.insert(
            lp.id,
            raw.insert_loop(Loop {
                face: lookup(&faces, lp.face, "face")?,
                fins: Vec::new(),
            })?,
        );
    }
    let mut edges = BTreeMap::new();
    let mut edge_senses = BTreeMap::new();
    let source_fins: BTreeMap<_, _> = source.half_edges.iter().map(|fin| (fin.id, fin)).collect();
    for edge in &source.edges {
        let curve = edge
            .curve
            .ok_or_else(|| Error::unimplemented(format!("edge {} has no shared curve", edge.id)))?;
        let orientation = lookup(&curve_senses, curve, "curve sense")?;
        // The source mapper orders endpoints by the first fin in the radial ring.
        // Normalize that ordering to the carrier's natural parameter direction.
        let primary_id = edge
            .half_edges
            .first()
            .ok_or_else(|| Error::parse("XT edge has no fin ring"))?;
        let primary = lookup(&source_fins, *primary_id, "primary fin")?;
        let endpoint_orientation = combine(orientation, sense(primary.sense)?);
        let (start, end) = if endpoint_orientation == Sense::Positive {
            (edge.start_vertex, edge.end_vertex)
        } else {
            (edge.end_vertex, edge.start_vertex)
        };
        edges.insert(
            edge.id,
            raw.insert_edge(Edge {
                fins: Vec::new(),
                start: optional(&vertices, start, "vertex")?,
                end: optional(&vertices, end, "vertex")?,
                curve: Some(lookup(&curves, curve, "curve")?),
                tolerance: tolerance(edge.tolerance, &precision)?,
            })?,
        );
        edge_senses.insert(edge.id, orientation);
    }
    let mut fins = BTreeMap::new();
    for fin in &source.half_edges {
        if fin.dummy {
            continue;
        }
        let edge = fin.edge.ok_or_else(|| {
            Error::unimplemented("isolated-loop XT fins have no kernel representation")
        })?;
        let lp = fin
            .loop_id
            .ok_or_else(|| Error::parse("non-dummy XT fin lacks a loop"))?;
        fins.insert(
            fin.id,
            raw.insert_fin(Fin {
                owner: lookup(&loops, lp, "loop")?,
                edge: lookup(&edges, edge, "edge")?,
                sense: combine(sense(fin.sense)?, lookup(&edge_senses, edge, "edge sense")?),
            })?,
        );
    }
    for body in &source.bodies {
        raw.body_mut(lookup(&bodies, body.id, "body")?)?.regions =
            ids(&regions, &body.regions, "region")?;
    }
    for region in &source.regions {
        raw.region_mut(lookup(&regions, region.id, "region")?)?
            .shells = ids(&shells, &region.shells, "shell")?;
    }
    for shell in &source.shells {
        let item = raw.shell_mut(lookup(&shells, shell.id, "shell")?)?;
        item.faces = ids(&faces, &shell.front_faces, "face")?;
        for id in ids(&faces, &shell.back_faces, "face")? {
            if !item.faces.contains(&id) {
                item.faces.push(id);
            }
        }
        item.wire_edges = ids(&edges, &shell.wire_edges, "edge")?;
    }
    for face in &source.faces {
        raw.face_mut(lookup(&faces, face.id, "face")?)?.loops = ids(&loops, &face.loops, "loop")?;
    }
    for lp in &source.loops {
        raw.loop_mut(lookup(&loops, lp.id, "loop")?)?.fins = ids(&fins, &lp.half_edges, "fin")?;
    }
    for edge in &source.edges {
        raw.edge_mut(lookup(&edges, edge.id, "edge")?)?.fins = edge
            .half_edges
            .iter()
            .filter_map(|id| fins.get(id).copied())
            .collect();
    }
    // Check the source's complete owner chains against the traversable target.
    // Inserting an entity into an arena is insufficient if it is then omitted
    // from the body graph used by the checker and public inquiry functions.
    for body in &source.bodies {
        let id = lookup(&bodies, body.id, "body")?;
        let actual_edges: HashSet<_> = model.body_edges(id)?.into_iter().collect();
        let expected_edges: HashSet<_> = ids(&edges, &body.edges, "owned edge")?
            .into_iter()
            .collect();
        let actual_vertices: HashSet<_> = model.body_vertices(id)?.into_iter().collect();
        let expected_vertices: HashSet<_> = ids(&vertices, &body.vertices, "owned vertex")?
            .into_iter()
            .collect();
        if actual_edges != expected_edges || actual_vertices != expected_vertices {
            return Err(Error::unimplemented(format!(
                "XT body {} has source entities outside the mapped face boundary graph",
                body.source.node_index
            )));
        }
    }
    let bodies = source
        .bodies
        .iter()
        .map(|body| lookup(&bodies, body.id, "body"))
        .collect::<Result<Vec<_>>>()?;
    let reports = bodies
        .iter()
        .map(|&body| check_body(&model, body, &CheckOptions::default()))
        .collect::<Result<Vec<_>>>()?;
    Ok(Decoded {
        model,
        bodies,
        reports,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn homogeneous_poles_retain_rational_weights_and_small_common_scales() {
        // A positive subnormal common scale is valid homogeneous data. Its
        // magnitude is not a geometric degeneracy or tolerance decision.
        let scale = f64::MIN_POSITIVE / 2.0;
        let vertices = vec![
            vec![scale, scale * 2.0, 0.0, scale],
            vec![scale * 6.0, scale * 8.0, 0.0, scale * 2.0],
        ];
        let (points, weights) = poles(&vertices, true).unwrap();
        assert!(points[0].distance(Point3::new(1.0, 2.0, 0.0)) <= Precision::DEFAULT.linear());
        assert!(points[1].distance(Point3::new(3.0, 4.0, 0.0)) <= Precision::DEFAULT.linear());
        assert_eq!(weights.unwrap(), vec![scale, scale * 2.0]);
    }

    #[test]
    fn oversized_knot_multiplicities_fail_before_allocation() {
        assert_eq!(
            expand_knots(&[0.0, 1.0], &[u16::MAX, u16::MAX], 2, 1)
                .unwrap_err()
                .code(),
            ErrorCode::Parse
        );
        assert_eq!(
            expand_knots(&[0.0], &[1], usize::MAX, 1)
                .unwrap_err()
                .code(),
            ErrorCode::CapacityExceeded
        );
    }
}
