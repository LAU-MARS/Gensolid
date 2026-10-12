//! The model (partition): all topology and geometry of a session.

use std::collections::HashSet;
use std::hash::Hash;

use gensolid_core::{Arena, Error, Id, Precision, Result};
use gensolid_geom::GeomStore;

use crate::entity::*;

/// Container for topology and geometry sharing one precision.
pub struct Model {
    precision: Precision,
    geom: GeomStore,
    bodies: Arena<Body>,
    regions: Arena<Region>,
    shells: Arena<Shell>,
    faces: Arena<Face>,
    loops: Arena<Loop>,
    fins: Arena<Fin>,
    edges: Arena<Edge>,
    vertices: Arena<Vertex>,
}

fn lookup<'a, T>(arena: &'a Arena<T>, id: Id<T>, what: &str) -> Result<&'a T> {
    arena
        .get(id)
        .ok_or_else(|| Error::not_found(format!("{what} {id:?}")))
}

fn lookup_mut<'a, T>(arena: &'a mut Arena<T>, id: Id<T>, what: &str) -> Result<&'a mut T> {
    arena
        .get_mut(id)
        .ok_or_else(|| Error::not_found(format!("{what} {id:?}")))
}

fn push_unique<T: Copy + Eq + Hash>(seen: &mut HashSet<T>, out: &mut Vec<T>, item: T) {
    if seen.insert(item) {
        out.push(item);
    }
}

impl Model {
    /// Creates an empty model.
    pub fn new(precision: Precision) -> Self {
        Self {
            precision,
            geom: GeomStore::new(),
            bodies: Arena::new(),
            regions: Arena::new(),
            shells: Arena::new(),
            faces: Arena::new(),
            loops: Arena::new(),
            fins: Arena::new(),
            edges: Arena::new(),
            vertices: Arena::new(),
        }
    }

    /// Session precision.
    pub const fn precision(&self) -> &Precision {
        &self.precision
    }

    /// Geometry.
    pub const fn geom(&self) -> &GeomStore {
        &self.geom
    }

    /// Low-level editor that bypasses all validity rules.
    ///
    /// Used by modelling operations and file readers; anything built with it must pass
    /// `gensolid-check` before it is handed to callers.
    pub fn raw(&mut self) -> RawEditor<'_> {
        RawEditor { model: self }
    }

    /// Looks up a body.
    pub fn body(&self, id: BodyId) -> Result<&Body> {
        lookup(&self.bodies, id, "body")
    }

    /// Looks up a region.
    pub fn region(&self, id: RegionId) -> Result<&Region> {
        lookup(&self.regions, id, "region")
    }

    /// Looks up a shell.
    pub fn shell(&self, id: ShellId) -> Result<&Shell> {
        lookup(&self.shells, id, "shell")
    }

    /// Looks up a face.
    pub fn face(&self, id: FaceId) -> Result<&Face> {
        lookup(&self.faces, id, "face")
    }

    /// Looks up a loop.
    pub fn loop_(&self, id: LoopId) -> Result<&Loop> {
        lookup(&self.loops, id, "loop")
    }

    /// Looks up a fin.
    pub fn fin(&self, id: FinId) -> Result<&Fin> {
        lookup(&self.fins, id, "fin")
    }

    /// Looks up an edge.
    pub fn edge(&self, id: EdgeId) -> Result<&Edge> {
        lookup(&self.edges, id, "edge")
    }

    /// Looks up a vertex.
    pub fn vertex(&self, id: VertexId) -> Result<&Vertex> {
        lookup(&self.vertices, id, "vertex")
    }

    /// Whether `entity` exists.
    pub fn contains(&self, entity: EntityRef) -> bool {
        match entity {
            EntityRef::Body(id) => self.bodies.contains(id),
            EntityRef::Region(id) => self.regions.contains(id),
            EntityRef::Shell(id) => self.shells.contains(id),
            EntityRef::Face(id) => self.faces.contains(id),
            EntityRef::Loop(id) => self.loops.contains(id),
            EntityRef::Fin(id) => self.fins.contains(id),
            EntityRef::Edge(id) => self.edges.contains(id),
            EntityRef::Vertex(id) => self.vertices.contains(id),
        }
    }

    /// All bodies.
    pub fn bodies(&self) -> impl Iterator<Item = BodyId> + '_ {
        self.bodies.iter().map(|(id, _)| id)
    }

    /// Faces of a body, each once, in shell order.
    pub fn body_faces(&self, body: BodyId) -> Result<Vec<FaceId>> {
        let (mut seen, mut out) = (HashSet::new(), Vec::new());
        for &region in &self.body(body)?.regions {
            for &shell in &self.region(region)?.shells {
                for &face in &self.shell(shell)?.faces {
                    push_unique(&mut seen, &mut out, face);
                }
            }
        }
        Ok(out)
    }

    /// Edges of a body (face edges and wire edges), each once.
    pub fn body_edges(&self, body: BodyId) -> Result<Vec<EdgeId>> {
        let (mut seen, mut out) = (HashSet::new(), Vec::new());
        for face in self.body_faces(body)? {
            for &lp in &self.face(face)?.loops {
                for &fin in &self.loop_(lp)?.fins {
                    push_unique(&mut seen, &mut out, self.fin(fin)?.edge);
                }
            }
        }
        for &region in &self.body(body)?.regions {
            for &shell in &self.region(region)?.shells {
                for &edge in &self.shell(shell)?.wire_edges {
                    push_unique(&mut seen, &mut out, edge);
                }
            }
        }
        Ok(out)
    }

    /// Vertices of a body (edge vertices and acorn vertices), each once.
    pub fn body_vertices(&self, body: BodyId) -> Result<Vec<VertexId>> {
        let (mut seen, mut out) = (HashSet::new(), Vec::new());
        for edge in self.body_edges(body)? {
            let e = self.edge(edge)?;
            for vertex in [e.start, e.end].into_iter().flatten() {
                push_unique(&mut seen, &mut out, vertex);
            }
        }
        for &region in &self.body(body)?.regions {
            for &shell in &self.region(region)?.shells {
                for &vertex in &self.shell(shell)?.acorn_vertices {
                    push_unique(&mut seen, &mut out, vertex);
                }
            }
        }
        Ok(out)
    }

    /// Start and end vertices of a fin in the fin's direction.
    pub fn fin_vertices(&self, fin: FinId) -> Result<(Option<VertexId>, Option<VertexId>)> {
        let f = self.fin(fin)?;
        let e = self.edge(f.edge)?;
        Ok(match f.sense {
            Sense::Positive => (e.start, e.end),
            Sense::Negative => (e.end, e.start),
        })
    }
}

/// Unchecked mutation of a [`Model`]; see [`Model::raw`].
pub struct RawEditor<'a> {
    model: &'a mut Model,
}

macro_rules! raw_accessors {
    ($($insert:ident, $get_mut:ident, $field:ident, $ty:ty, $what:literal;)*) => {
        $(
            #[doc = concat!("Inserts a ", $what, " without checking its references.")]
            pub fn $insert(&mut self, value: $ty) -> Result<Id<$ty>> {
                self.model.$field.insert(value)
            }

            #[doc = concat!("Mutable access to a ", $what, ".")]
            pub fn $get_mut(&mut self, id: Id<$ty>) -> Result<&mut $ty> {
                lookup_mut(&mut self.model.$field, id, $what)
            }
        )*
    };
}

impl RawEditor<'_> {
    raw_accessors! {
        insert_body, body_mut, bodies, Body, "body";
        insert_region, region_mut, regions, Region, "region";
        insert_shell, shell_mut, shells, Shell, "shell";
        insert_face, face_mut, faces, Face, "face";
        insert_loop, loop_mut, loops, Loop, "loop";
        insert_fin, fin_mut, fins, Fin, "fin";
        insert_edge, edge_mut, edges, Edge, "edge";
        insert_vertex, vertex_mut, vertices, Vertex, "vertex";
    }

    /// Mutable geometry.
    pub fn geom_mut(&mut self) -> &mut GeomStore {
        &mut self.model.geom
    }

    /// Session precision.
    pub fn precision(&self) -> Precision {
        self.model.precision
    }
}
