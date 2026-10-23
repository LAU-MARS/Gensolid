//! Transactional session entry points for owned surfaces and checked topology edits.
use super::*;
use gensolid_geom::BSurface;
use gensolid_topo::FaceId;

/// Tags created by a checked bounded-sheet cut.
#[derive(Clone, Debug)]
pub struct FaceSplit {
    /// Newly created body; the source is retained.
    pub body: Tag,
    /// Replacement faces.
    pub faces: [Tag; 2],
    /// Shared intersection edge.
    pub intersection: Tag,
}

impl Session {
    fn staged<T>(&mut self, edit: impl FnOnce(&mut Self) -> Result<T>) -> Result<T> {
        let mut staged = self.clone();
        let result = edit(&mut staged)?;
        *self = staged;
        Ok(result)
    }
    fn face_id(&self, tag: Tag) -> Result<FaceId> {
        match self.entity(tag)? {
            EntityRef::Face(id) => Ok(id),
            _ => Err(Error::new(
                ErrorCode::WrongEntityClass,
                "expected a face tag",
            )),
        }
    }
    fn body_tags(&mut self, bodies: Vec<BodyId>) -> Result<Vec<Tag>> {
        self.ensure_tag_capacity(bodies.len())?;
        bodies
            .into_iter()
            .map(|b| self.tag(EntityRef::Body(b)))
            .collect()
    }
    /// Create a checked regular graph sheet from owned NURBS data.
    pub fn body_create_sheet_surface(&mut self, surface: &BSurface) -> Result<Tag> {
        self.staged(|s| {
            let b = gensolid_modeling::body_create_sheet_surface(&mut s.model, surface)?;
            s.tag(EntityRef::Body(b))
        })
    }
    /// Construct a surface using a borrowed backend and publish a checked sheet.
    pub fn body_create_sheet_backend(
        &mut self,
        backend: &dyn GeometryBackend,
        request: &SurfaceConstruction,
    ) -> Result<Tag> {
        self.staged(|s| {
            let b = gensolid_modeling::body_create_sheet_backend(&mut s.model, backend, request)?;
            s.tag(EntityRef::Body(b))
        })
    }
    /// Make a directional solid prism from a regular graph patch, retaining curved faces.
    pub fn body_create_solid_surface_prism(
        &mut self,
        surface: &BSurface,
        height: f64,
    ) -> Result<Tag> {
        self.staged(|s| {
            let b =
                gensolid_modeling::body_create_solid_surface_prism(&mut s.model, surface, height)?;
            s.tag(EntityRef::Body(b))
        })
    }
    /// Construct a regular surface via a borrowed provider, then make a checked prism.
    pub fn body_create_solid_backend_prism(
        &mut self,
        backend: &dyn GeometryBackend,
        request: &SurfaceConstruction,
        height: f64,
    ) -> Result<Tag> {
        self.staged(|s| {
            if backend.precision() != *s.model.precision() {
                return Err(Error::invalid_argument(
                    "backend and session precisions differ",
                ));
            }
            let surface = backend.construct_surface(request)?;
            let b =
                gensolid_modeling::body_create_solid_surface_prism(&mut s.model, &surface, height)?;
            s.tag(EntityRef::Body(b))
        })
    }
    /// Copy all faces into separately checked sheets.
    pub fn body_extract_sheets(&mut self, body: Tag) -> Result<Vec<Tag>> {
        self.staged(|s| {
            let id = s.body_id(body)?;
            let out = gensolid_modeling::body_extract_sheets(&mut s.model, id)?;
            s.body_tags(out)
        })
    }
    /// Sew oriented planar sheets into a convex solid; source sheets remain valid.
    pub fn body_sew_sheets(&mut self, sheets: &[Tag]) -> Result<Tag> {
        self.staged(|s| {
            let ids = sheets
                .iter()
                .map(|&b| s.body_id(b))
                .collect::<Result<Vec<_>>>()?;
            let b = gensolid_modeling::body_sew_sheets(&mut s.model, &ids)?;
            s.tag(EntityRef::Body(b))
        })
    }
    /// Create a subdivided copy of a body, propagating edge splits to neighboring faces.
    pub fn body_split_face_isoparametric(
        &mut self,
        body: Tag,
        face: Tag,
        constant_u: bool,
        fraction: f64,
    ) -> Result<Tag> {
        self.staged(|s| {
            let b = s.body_id(body)?;
            let f = s.face_id(face)?;
            let out = gensolid_modeling::body_split_face_isoparametric(
                &mut s.model,
                b,
                f,
                constant_u,
                fraction,
            )?;
            s.tag(EntityRef::Body(out))
        })
    }
    /// Cut a ruled NURBS face with a bounded planar face. Uses a borrowed backend
    /// for candidate generation when present, with independent kernel validation.
    pub fn body_split_face_by_sheet(
        &mut self,
        body: Tag,
        face: Tag,
        tool: Tag,
        backend: Option<&dyn GeometryBackend>,
    ) -> Result<FaceSplit> {
        self.staged(|s| {
            let b = s.body_id(body)?;
            let f = s.face_id(face)?;
            let tool = s.face_id(tool)?;
            let out =
                gensolid_modeling::body_split_face_by_sheet(&mut s.model, b, f, tool, backend)?;
            s.ensure_tag_capacity(4)?;
            Ok(FaceSplit {
                body: s.tag(EntityRef::Body(out.body))?,
                faces: [
                    s.tag(EntityRef::Face(out.faces[0]))?,
                    s.tag(EntityRef::Face(out.faces[1]))?,
                ],
                intersection: s.tag(EntityRef::Edge(out.intersection))?,
            })
        })
    }
    /// Convex planar Boolean with internal interfaces removed, one body per component.
    pub fn body_boolean(
        &mut self,
        a: Tag,
        b: Tag,
        operation: BooleanOperation,
    ) -> Result<Vec<Tag>> {
        self.staged(|s| {
            let a = s.body_id(a)?;
            let b = s.body_id(b)?;
            let out = gensolid_modeling::body_boolean(&mut s.model, a, b, operation)?;
            s.body_tags(out)
        })
    }
    /// Convex planar regularized Boolean, returned as bodies with disjoint interiors.
    /// Adjacent cells retain internal interfaces; this is not a fused boundary result.
    pub fn body_boolean_cells(
        &mut self,
        a: Tag,
        b: Tag,
        operation: BooleanOperation,
    ) -> Result<Vec<Tag>> {
        self.staged(|s| {
            let a = s.body_id(a)?;
            let b = s.body_id(b)?;
            let out = gensolid_modeling::body_boolean_cells(&mut s.model, a, b, operation)?;
            s.body_tags(out)
        })
    }
    /// Offset a planar or supported circular cylindrical NURBS sheet through a provider.
    /// The complete output is independently certified before checked topology is published.
    pub fn body_offset_sheet_backend(
        &mut self,
        body: Tag,
        backend: &dyn GeometryBackend,
        distance: f64,
    ) -> Result<Tag> {
        self.staged(|s| {
            let b = s.body_id(body)?;
            let out =
                gensolid_modeling::body_offset_sheet_backend(&mut s.model, b, backend, distance)?;
            s.tag(EntityRef::Body(out))
        })
    }
    /// Offset a simple convex planar solid with unchanged face incidence.
    pub fn body_offset_planar(&mut self, body: Tag, distance: f64) -> Result<Tag> {
        self.staged(|s| {
            let b = s.body_id(body)?;
            let out = gensolid_modeling::body_offset_planar(&mut s.model, b, distance)?;
            s.tag(EntityRef::Body(out))
        })
    }
    /// Remove one planar face and offset the remaining faces inward with a connected rim.
    pub fn body_shell_planar(&mut self, body: Tag, remove: Tag, thickness: f64) -> Result<Tag> {
        self.staged(|s| {
            let b = s.body_id(body)?;
            let f = s.face_id(remove)?;
            let out = gensolid_modeling::body_shell_planar(&mut s.model, b, f, thickness)?;
            s.tag(EntityRef::Body(out))
        })
    }
    /// Create a closed hollow planar solid with a finite empty cavity region.
    pub fn body_hollow_planar(&mut self, body: Tag, thickness: f64) -> Result<Tag> {
        self.staged(|s| {
            let b = s.body_id(body)?;
            let out = gensolid_modeling::body_hollow_planar(&mut s.model, b, thickness)?;
            s.tag(EntityRef::Body(out))
        })
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn failed_staged_edit_restores_geometry_topology_and_tags() {
        fn require_send<T: Send>() {}
        require_send::<Session>();
        let mut session = Session::start(SessionOptions::default());
        let original = session.body_create_solid_block(2., 3., 4., None).unwrap();
        let original_faces = session.body_ask_faces(original).unwrap();
        let mut control = session.clone();
        let result: Result<()> = session.staged(|s| {
            let body = s.body_create_solid_block(3., 4., 5., None)?;
            s.body_ask_faces(body)?;
            Err(Error::unimplemented(
                "injected failure after topology and tag allocation",
            ))
        });
        assert_eq!(result.unwrap_err().code(), ErrorCode::Unimplemented);
        assert_eq!(session.entities, control.entities);
        assert_eq!(session.tags, control.tags);
        assert_eq!(
            session.model.bodies().count(),
            control.model.bodies().count()
        );
        assert_eq!(session.body_ask_faces(original).unwrap(), original_faces);
        assert!(session.body_check(original).unwrap().is_empty());

        // A successful follow-up receives exactly the IDs of a session that never
        // attempted the failed operation, including geometry and arena allocation.
        let next = session.body_create_solid_block(1., 2., 3., None).unwrap();
        let expected = control.body_create_solid_block(1., 2., 3., None).unwrap();
        assert_eq!(next, expected);
        assert_eq!(
            session.body_ask_faces(next).unwrap(),
            control.body_ask_faces(expected).unwrap()
        );
        assert_eq!(session.entities, control.entities);
        for &entity in &session.entities {
            let EntityRef::Face(id) = entity else {
                continue;
            };
            assert_eq!(
                session.model.face(id).unwrap().surface,
                control.model.face(id).unwrap().surface
            );
        }
    }
}
