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
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{BodyId, EntityRef, Model, make_block};

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
        let precision = *self.model.precision();
        let frame = match basis {
            None => Frame::WORLD,
            Some(b) => Frame::new(
                point(b.location),
                vector(b.axis),
                vector(b.ref_direction),
                &precision,
            )?,
        };
        let body = make_block(&mut self.model, x, y, z, &frame)?;
        self.tag(EntityRef::Body(body))
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
}
