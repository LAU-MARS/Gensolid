//! PK-style application interface.
//!
//! Applications address entities by integer *tags*, as in Parasolid's PK interface, and
//! every function is named `<class>_<verb>_<object>` (`body_create_solid_block`,
//! `body_ask_faces`). The C ABI (`gensolid-ffi`), Python and WASM bindings are thin
//! wrappers over this crate, so all of them share one set of semantics.

use std::collections::HashMap;

pub use gensolid_check::FaultKind;
use gensolid_check::{CheckOptions, check_body};
pub use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_math::{Frame, Point2, Point3, Vec2, Vec3};
pub use gensolid_modeling::MassProperties;
use gensolid_topo::{BodyId, EdgeId, EntityRef, Model, body_bounding_box, body_volume};

/// Integer handle to an entity. Tags are never reused within a session.
pub type Tag = i32;

/// The null tag.
pub const NULL_TAG: Tag = 0;

/// Entity class of a tag.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum Class {
    /// Body.
    Body = 1,
    /// Region.
    Region = 2,
    /// Shell.
    Shell = 3,
    /// Face.
    Face = 4,
    /// Loop.
    Loop = 5,
    /// Fin.
    Fin = 6,
    /// Edge.
    Edge = 7,
    /// Vertex.
    Vertex = 8,
}

/// Session start options.
#[derive(Clone, Copy, Debug, Default)]
pub struct SessionOptions {
    /// Precision for the session's model.
    pub precision: Precision,
}

/// Coordinate system argument (PK's `AXIS2_sf`): location, principal axis, reference
/// direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Axis2 {
    /// Origin.
    pub location: [f64; 3],
    /// Principal (z) axis.
    pub axis: [f64; 3],
    /// Reference (x) direction.
    pub ref_direction: [f64; 3],
}

/// One fault reported by [`Session::body_check`].
#[derive(Clone, Debug, PartialEq)]
pub struct CheckFault {
    /// Fault category.
    pub state: FaultKind,
    /// Entity at fault.
    pub entity: Tag,
    /// Explanation.
    pub message: String,
}

/// Stable session tags returned by a checked line-edge subdivision.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EdgeSplit {
    /// Original edge, retained as the first segment.
    pub original_edge: Tag,
    /// New second segment.
    pub new_edge: Tag,
    /// New vertex shared by the two segments.
    pub vertex: Tag,
}

/// A corresponding convex polygon at a local Z height for ruled solid lofting.
#[derive(Clone, Debug, PartialEq)]
pub struct LoftSection {
    /// Height in the operation's coordinate frame; strictly increasing across sections.
    pub height: f64,
    /// Local XY vertices; sections must have the same vertex count and correspondence.
    pub profile: Vec<[f64; 2]>,
}

/// A modelling session: one model plus its tag table.
pub struct Session {
    model: Model,
    entities: Vec<EntityRef>,
    tags: HashMap<EntityRef, Tag>,
}

fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}

fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

fn profile_points(profile: &[[f64; 2]]) -> Vec<Point2> {
    profile
        .iter()
        .map(|point| Point2::new(point[0], point[1]))
        .collect()
}

impl Session {
    /// Starts a session.
    pub fn start(options: SessionOptions) -> Self {
        Self {
            model: Model::new(options.precision),
            entities: Vec::new(),
            tags: HashMap::new(),
        }
    }

    /// The underlying model, for read-only inspection.
    pub const fn model(&self) -> &Model {
        &self.model
    }

    fn ensure_tag_capacity(&self, additional: usize) -> Result<()> {
        let count = self
            .entities
            .len()
            .checked_add(additional)
            .ok_or_else(|| Error::new(ErrorCode::CapacityExceeded, "tag space exhausted"))?;
        Tag::try_from(count)
            .map_err(|_| Error::new(ErrorCode::CapacityExceeded, "tag space exhausted"))?;
        Ok(())
    }

    fn basis_frame(&self, basis: Option<&Axis2>) -> Result<Frame> {
        match basis {
            None => Ok(Frame::WORLD),
            Some(basis) => Frame::new(
                point(basis.location),
                vector(basis.axis),
                vector(basis.ref_direction),
                self.model.precision(),
            ),
        }
    }

    fn tag(&mut self, entity: EntityRef) -> Result<Tag> {
        if let Some(&tag) = self.tags.get(&entity) {
            return Ok(tag);
        }
        let tag = Tag::try_from(self.entities.len() + 1)
            .map_err(|_| Error::new(ErrorCode::CapacityExceeded, "tag space exhausted"))?;
        self.entities.push(entity);
        self.tags.insert(entity, tag);
        Ok(tag)
    }

    fn tags(&mut self, entities: impl IntoIterator<Item = EntityRef>) -> Result<Vec<Tag>> {
        entities.into_iter().map(|e| self.tag(e)).collect()
    }

    fn entity(&self, tag: Tag) -> Result<EntityRef> {
        let entity = usize::try_from(tag)
            .ok()
            .and_then(|t| t.checked_sub(1))
            .and_then(|i| self.entities.get(i))
            .copied()
            .ok_or_else(|| Error::not_found(format!("tag {tag} was never issued")))?;
        if !self.model.contains(entity) {
            return Err(Error::not_found(format!(
                "tag {tag} refers to a deleted entity"
            )));
        }
        Ok(entity)
    }

    fn body_id(&self, tag: Tag) -> Result<BodyId> {
        match self.entity(tag)? {
            EntityRef::Body(id) => Ok(id),
            other => Err(Error::new(
                ErrorCode::WrongEntityClass,
                format!("tag {tag} is {other:?}, not a body"),
            )),
        }
    }

    fn edge_id(&self, tag: Tag) -> Result<EdgeId> {
        match self.entity(tag)? {
            EntityRef::Edge(id) => Ok(id),
            other => Err(Error::new(
                ErrorCode::WrongEntityClass,
                format!("tag {tag} is {other:?}, not an edge"),
            )),
        }
    }

    /// Class of the entity a tag refers to.
    pub fn entity_ask_class(&self, tag: Tag) -> Result<Class> {
        Ok(match self.entity(tag)? {
            EntityRef::Body(_) => Class::Body,
            EntityRef::Region(_) => Class::Region,
            EntityRef::Shell(_) => Class::Shell,
            EntityRef::Face(_) => Class::Face,
            EntityRef::Loop(_) => Class::Loop,
            EntityRef::Fin(_) => Class::Fin,
            EntityRef::Edge(_) => Class::Edge,
            EntityRef::Vertex(_) => Class::Vertex,
        })
    }

    /// Creates a solid block of size `x × y × z` centred on `basis.location` (the world
    /// origin when `None`), with edges along the basis axes.
    pub fn body_create_solid_block(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        basis: Option<&Axis2>,
    ) -> Result<Tag> {
        self.ensure_tag_capacity(1)?;
        let frame = self.basis_frame(basis)?;
        let body = gensolid_modeling::body_create_solid_block(&mut self.model, x, y, z, &frame)?;
        self.tag(EntityRef::Body(body))
    }

    /// Create a checked solid prism from a convex local XY polygon and a positive
    /// height. The polygon may have either winding and must not repeat its first point.
    pub fn body_create_solid_prism(
        &mut self,
        profile: &[[f64; 2]],
        height: f64,
        basis: Option<&Axis2>,
    ) -> Result<Tag> {
        self.ensure_tag_capacity(1)?;
        let frame = self.basis_frame(basis)?;
        let body = gensolid_modeling::body_create_solid_prism(
            &mut self.model,
            &profile_points(profile),
            height,
            &frame,
        )?;
        self.tag(EntityRef::Body(body))
    }

    /// Create a checked linear sweep of a convex local XY polygon by a finite
    /// world-coordinate displacement. Either normal direction is supported.
    pub fn body_create_solid_linear_sweep(
        &mut self,
        profile: &[[f64; 2]],
        displacement: [f64; 3],
        basis: Option<&Axis2>,
    ) -> Result<Tag> {
        self.ensure_tag_capacity(1)?;
        let frame = self.basis_frame(basis)?;
        let body = gensolid_modeling::body_create_solid_linear_sweep(
            &mut self.model,
            &profile_points(profile),
            vector(displacement),
            &frame,
        )?;
        self.tag(EntityRef::Body(body))
    }

    /// Create a checked planar homothetic loft to a positively scaled copy of a
    /// convex local XY polygon. The top translation uses the basis XY axes.
    pub fn body_create_solid_homothetic_loft(
        &mut self,
        profile: &[[f64; 2]],
        height: f64,
        top_scale: f64,
        top_translation: [f64; 2],
        basis: Option<&Axis2>,
    ) -> Result<Tag> {
        self.ensure_tag_capacity(1)?;
        let frame = self.basis_frame(basis)?;
        let body = gensolid_modeling::body_create_solid_homothetic_loft(
            &mut self.model,
            &profile_points(profile),
            height,
            top_scale,
            Vec2::new(top_translation[0], top_translation[1]),
            &frame,
        )?;
        self.tag(EntityRef::Body(body))
    }

    /// Revolve a convex (radius, height) polygon about local Z by `angle` radians.
    /// Supports partial and full turns; touching/crossing the axis is unsupported.
    pub fn body_create_solid_revolution(
        &mut self,
        profile: &[[f64; 2]],
        angle: f64,
        basis: Option<&Axis2>,
    ) -> Result<Tag> {
        self.ensure_tag_capacity(1)?;
        let frame = self.basis_frame(basis)?;
        let result = gensolid_modeling::body_create_solid_revolution(
            &mut self.model,
            &profile_points(profile),
            angle,
            &frame,
        )?;
        self.tag(EntityRef::Body(result.body))
    }

    /// Create a checked ruled solid through corresponding convex polygon sections.
    /// Intermediate sections must remain strictly convex; invalid correspondence is rejected.
    pub fn body_create_solid_loft_sections(
        &mut self,
        sections: &[LoftSection],
        basis: Option<&Axis2>,
    ) -> Result<Tag> {
        self.ensure_tag_capacity(1)?;
        let frame = self.basis_frame(basis)?;
        let sections = sections
            .iter()
            .map(|s| gensolid_topo::LoftSection {
                height: s.height,
                profile: profile_points(&s.profile),
            })
            .collect::<Vec<_>>();
        let result =
            gensolid_modeling::body_create_solid_loft_sections(&mut self.model, &sections, &frame)?;
        self.tag(EntityRef::Body(result.body))
    }

    /// Create a checked convex polyhedron from world-coordinate vertices and
    /// outward counter-clockwise faces using zero-based vertex indices.
    pub fn body_create_solid_convex_polyhedron(
        &mut self,
        vertices: &[[f64; 3]],
        faces: &[Vec<usize>],
    ) -> Result<Tag> {
        self.ensure_tag_capacity(1)?;
        let points: Vec<_> = vertices.iter().copied().map(point).collect();
        let body = gensolid_modeling::body_create_solid_convex_polyhedron(
            &mut self.model,
            &points,
            faces,
        )?;
        self.tag(EntityRef::Body(body))
    }

    /// Split a bounded line edge, checking every body sharing it before and after
    /// subdivision. Existing tags stay valid; the new edge and vertex get fresh tags.
    pub fn edge_split_line(&mut self, edge: Tag, fraction: f64) -> Result<EdgeSplit> {
        let id = self.edge_id(edge)?;
        self.ensure_tag_capacity(2)?;
        let split = gensolid_modeling::edge_split_line(&mut self.model, id, fraction)?;
        Ok(EdgeSplit {
            original_edge: edge,
            new_edge: self.tag(EntityRef::Edge(split.new_edge))?,
            vertex: self.tag(EntityRef::Vertex(split.vertex))?,
        })
    }

    /// Faces of a body.
    pub fn body_ask_faces(&mut self, body: Tag) -> Result<Vec<Tag>> {
        let faces = self.model.body_faces(self.body_id(body)?)?;
        self.tags(faces.into_iter().map(EntityRef::Face))
    }

    /// Edges of a body.
    pub fn body_ask_edges(&mut self, body: Tag) -> Result<Vec<Tag>> {
        let edges = self.model.body_edges(self.body_id(body)?)?;
        self.tags(edges.into_iter().map(EntityRef::Edge))
    }

    /// Vertices of a body.
    pub fn body_ask_vertices(&mut self, body: Tag) -> Result<Vec<Tag>> {
        let vertices = self.model.body_vertices(self.body_id(body)?)?;
        self.tags(vertices.into_iter().map(EntityRef::Vertex))
    }

    /// Checks a body; an empty result means the body is valid. Checks that could not be
    /// run are reported as an [`ErrorCode::Unimplemented`] error rather than silently
    /// passing.
    pub fn body_check(&mut self, body: Tag) -> Result<Vec<CheckFault>> {
        let report = check_body(&self.model, self.body_id(body)?, &CheckOptions::default())?;
        if report.faults.is_empty() && report.unchecked > 0 {
            return Err(Error::unimplemented(format!(
                "{} geometric checks are not available for this body's geometry",
                report.unchecked
            )));
        }
        report
            .faults
            .into_iter()
            .map(|f| {
                Ok(CheckFault {
                    state: f.kind,
                    entity: self.tag(f.entity)?,
                    message: f.message,
                })
            })
            .collect()
    }

    /// Axis-aligned bounding box of a body, as `(min_xyz, max_xyz)`.
    pub fn body_ask_bounding_box(&self, body: Tag) -> Result<([f64; 3], [f64; 3])> {
        let box_ = body_bounding_box(&self.model, self.body_id(body)?)?;
        let (min, max) = (box_.min(), box_.max());
        Ok(([min.x, min.y, min.z], [max.x, max.y, max.z]))
    }

    /// Volume of the solid regions of a body.
    pub fn body_ask_volume(&self, body: Tag) -> Result<f64> {
        body_volume(&self.model, self.body_id(body)?)
    }

    /// Checked unit-density volume, centroid, boundary area and centroidal inertia
    /// for planar straight-edged material boundaries, including holes and cavities.
    pub fn body_ask_mass_properties(&self, body: Tag) -> Result<MassProperties> {
        gensolid_modeling::body_ask_mass_properties(&self.model, self.body_id(body)?)
    }

    /// Position of a vertex.
    pub fn vertex_ask_point(&self, vertex: Tag) -> Result<[f64; 3]> {
        match self.entity(vertex)? {
            EntityRef::Vertex(id) => {
                let p = self.model.vertex(id)?.point;
                Ok([p.x, p.y, p.z])
            }
            other => Err(Error::new(
                ErrorCode::WrongEntityClass,
                format!("tag {vertex} is {other:?}, not a vertex"),
            )),
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::float_cmp)] // Half-sizes of the block are exactly representable.
    fn block_round_trip_through_tags() {
        let mut s = Session::start(SessionOptions::default());
        let body = s.body_create_solid_block(1.0, 2.0, 3.0, None).unwrap();
        assert_eq!(s.entity_ask_class(body).unwrap(), Class::Body);
        let faces = s.body_ask_faces(body).unwrap();
        assert_eq!(faces.len(), 6);
        assert_eq!(s.body_ask_faces(body).unwrap(), faces, "tags are stable");
        assert_eq!(s.body_ask_edges(body).unwrap().len(), 12);
        let vertices = s.body_ask_vertices(body).unwrap();
        assert_eq!(vertices.len(), 8);
        assert_eq!(s.vertex_ask_point(vertices[0]).unwrap(), [-0.5, -1.0, -1.5]);
        assert!(s.body_check(body).unwrap().is_empty());
        assert_eq!(s.body_ask_volume(body).unwrap(), 6.0);
        assert_eq!(
            s.body_ask_bounding_box(body).unwrap(),
            ([-0.5, -1.0, -1.5], [0.5, 1.0, 1.5])
        );
    }

    #[test]
    fn rejects_bad_tags() {
        let mut s = Session::start(SessionOptions::default());
        assert_eq!(
            s.entity_ask_class(NULL_TAG).unwrap_err().code(),
            ErrorCode::NotFound
        );
        assert_eq!(
            s.entity_ask_class(-3).unwrap_err().code(),
            ErrorCode::NotFound
        );
        let body = s.body_create_solid_block(1.0, 1.0, 1.0, None).unwrap();
        let face = s.body_ask_faces(body).unwrap()[0];
        assert_eq!(
            s.body_ask_faces(face).unwrap_err().code(),
            ErrorCode::WrongEntityClass
        );
    }

    #[test]
    fn honours_basis_set() {
        let mut s = Session::start(SessionOptions::default());
        let basis = Axis2 {
            location: [10.0, 0.0, 0.0],
            axis: [0.0, 0.0, 1.0],
            ref_direction: [0.0, 1.0, 0.0],
        };
        let body = s
            .body_create_solid_block(2.0, 2.0, 2.0, Some(&basis))
            .unwrap();
        let tol = s.model().precision().linear();
        for v in s.body_ask_vertices(body).unwrap() {
            let [x, y, z] = s.vertex_ask_point(v).unwrap();
            assert!(
                ((x - 10.0).abs() - 1.0).abs() <= tol
                    && (y.abs() - 1.0).abs() <= tol
                    && (z.abs() - 1.0).abs() <= tol
            );
        }
    }

    #[test]
    fn checked_polyhedral_operations_and_mass_are_available_through_tags() {
        let mut session = Session::start(SessionOptions::default());
        let p = *session.model().precision();
        let profile = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        let basis = Axis2 {
            location: [3.0, 2.0, 1.0],
            axis: [0.0, 1.0, 1.0],
            ref_direction: [1.0, 0.0, 0.0],
        };
        let frame = session.basis_frame(Some(&basis)).unwrap();
        let prism = session
            .body_create_solid_prism(&profile, 3.0, Some(&basis))
            .unwrap();
        let sweep = session
            .body_create_solid_linear_sweep(&profile, [1.0, 0.0, -3.0], None)
            .unwrap();
        let loft = session
            .body_create_solid_homothetic_loft(&profile, 3.0, 2.0, [0.5, -0.25], None)
            .unwrap();
        let tetrahedron = session
            .body_create_solid_convex_polyhedron(
                &[
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 2.0, 0.0],
                    [0.0, 0.0, 3.0],
                ],
                &[vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 3]],
            )
            .unwrap();
        for (body, volume) in [
            (prism, 12.0),
            (sweep, 12.0),
            (loft, 28.0),
            (tetrahedron, 1.0),
        ] {
            assert_eq!(session.entity_ask_class(body).unwrap(), Class::Body);
            assert!(session.body_check(body).unwrap().is_empty());
            let mass = session.body_ask_mass_properties(body).unwrap();
            assert!((mass.volume - volume).abs() <= p.linear() * mass.surface_area);
        }
        assert!(
            session
                .body_ask_mass_properties(prism)
                .unwrap()
                .centroid
                .distance(frame.origin() + frame.z() * 1.5)
                <= p.linear()
        );
        assert!(
            session
                .body_ask_mass_properties(sweep)
                .unwrap()
                .centroid
                .distance(Point3::new(0.5, 0.0, -1.5))
                <= p.linear()
        );
        assert!(
            session
                .body_ask_mass_properties(tetrahedron)
                .unwrap()
                .centroid
                .distance(Point3::new(0.25, 0.5, 0.75))
                <= p.linear()
        );
    }

    #[test]
    fn checked_euler_split_preserves_old_tags_and_publishes_new_entities() {
        let mut session = Session::start(SessionOptions::default());
        let body = session
            .body_create_solid_block(1.0, 2.0, 3.0, None)
            .unwrap();
        let old_edges = session.body_ask_edges(body).unwrap();
        let old_vertices = session.body_ask_vertices(body).unwrap();
        let before = session.body_ask_mass_properties(body).unwrap();
        let edge = old_edges[0];
        let id = session.edge_id(edge).unwrap();
        let start = session
            .model
            .vertex(session.model.edge(id).unwrap().start.unwrap())
            .unwrap()
            .point;
        let end = session
            .model
            .vertex(session.model.edge(id).unwrap().end.unwrap())
            .unwrap()
            .point;
        let split = session.edge_split_line(edge, 0.4).unwrap();
        assert_eq!(split.original_edge, edge);
        assert_eq!(
            session.entity_ask_class(split.new_edge).unwrap(),
            Class::Edge
        );
        assert_eq!(
            session.entity_ask_class(split.vertex).unwrap(),
            Class::Vertex
        );
        let p = session.model.precision().linear();
        assert!(
            point(session.vertex_ask_point(split.vertex).unwrap())
                .distance(start + (end - start) * 0.4)
                <= p
        );
        let edges = session.body_ask_edges(body).unwrap();
        let vertices = session.body_ask_vertices(body).unwrap();
        assert_eq!(edges.len(), 13);
        assert_eq!(vertices.len(), 9);
        assert!(old_edges.iter().all(|edge| edges.contains(edge)));
        assert!(old_vertices.iter().all(|vertex| vertices.contains(vertex)));
        assert!(edges.contains(&split.new_edge));
        assert!(vertices.contains(&split.vertex));
        assert!(session.body_check(body).unwrap().is_empty());
        let after = session.body_ask_mass_properties(body).unwrap();
        assert!((after.volume - before.volume).abs() <= p * after.surface_area);
        assert_eq!(
            session.edge_split_line(body, 0.5).unwrap_err().code(),
            ErrorCode::WrongEntityClass
        );
        assert_eq!(
            session.body_ask_mass_properties(edge).unwrap_err().code(),
            ErrorCode::WrongEntityClass
        );
    }

    #[test]
    fn invalid_checked_constructors_do_not_publish_bodies_or_tags() {
        let mut session = Session::start(SessionOptions::default());
        let concave = [[0.0, 0.0], [2.0, 0.0], [1.0, 0.25], [2.0, 2.0], [0.0, 2.0]];
        assert_eq!(
            session
                .body_create_solid_prism(&concave, 1.0, None)
                .unwrap_err()
                .code(),
            ErrorCode::InvalidArgument
        );
        let profile = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        assert_eq!(
            session
                .body_create_solid_linear_sweep(&profile, [1.0, 0.0, 0.0], None)
                .unwrap_err()
                .code(),
            ErrorCode::Degenerate
        );
        assert_eq!(
            session
                .body_create_solid_homothetic_loft(&profile, 1.0, 0.0, [0.0, 0.0], None)
                .unwrap_err()
                .code(),
            ErrorCode::InvalidArgument
        );
        assert_eq!(session.model.bodies().count(), 0);
        assert!(session.entities.is_empty());
    }
}
