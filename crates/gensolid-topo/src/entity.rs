//! Topological entities.

use gensolid_core::{Id, Tolerance};
use gensolid_geom::{CurveId, SurfaceId};
use gensolid_math::Point3;

/// Handle to a [`Body`].
pub type BodyId = Id<Body>;
/// Handle to a [`Region`].
pub type RegionId = Id<Region>;
/// Handle to a [`Shell`].
pub type ShellId = Id<Shell>;
/// Handle to a [`Face`].
pub type FaceId = Id<Face>;
/// Handle to a [`Loop`].
pub type LoopId = Id<Loop>;
/// Handle to a [`Fin`].
pub type FinId = Id<Fin>;
/// Handle to an [`Edge`].
pub type EdgeId = Id<Edge>;
/// Handle to a [`Vertex`].
pub type VertexId = Id<Vertex>;

/// Orientation of an entity relative to the geometry or entity it uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Sense {
    /// Same direction.
    Positive,
    /// Opposite direction.
    Negative,
}

impl Sense {
    /// The opposite sense.
    pub const fn reversed(self) -> Self {
        match self {
            Self::Positive => Self::Negative,
            Self::Negative => Self::Positive,
        }
    }
}

/// What a body is allowed to contain. Operations keep the body within its kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BodyKind {
    /// A single isolated vertex.
    Acorn,
    /// Edges only.
    Wire,
    /// Faces without enclosed volume.
    Sheet,
    /// Manifold solid: one or more solid regions bounded by closed shells.
    Solid,
    /// Anything, including mixed-dimensional and non-manifold cellular models.
    General,
}

/// Top-level topological entity.
#[derive(Clone, Debug, PartialEq)]
pub struct Body {
    /// Kind constraint.
    pub kind: BodyKind,
    /// Regions; `regions[0]` is the infinite void region.
    pub regions: Vec<RegionId>,
}

/// Connected volume of space, either material (`solid`) or void.
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    /// Owning body.
    pub body: BodyId,
    /// Whether the region is material.
    pub solid: bool,
    /// Boundary components.
    pub shells: Vec<ShellId>,
}

/// Connected component of the boundary of a region.
#[derive(Clone, Debug, PartialEq)]
pub struct Shell {
    /// Region this shell bounds.
    pub region: RegionId,
    /// Faces; each face records on which of its sides this shell lies.
    pub faces: Vec<FaceId>,
    /// Edges with no faces (wireframe).
    pub wire_edges: Vec<EdgeId>,
    /// Isolated vertices.
    pub acorn_vertices: Vec<VertexId>,
}

/// Bounded portion of a surface.
///
/// The face normal is the surface normal flipped when `sense` is negative. The *front*
/// shell lies on the side the face normal points to, the *back* shell on the other.
#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    /// Shell on the side the face normal points into.
    pub front: Option<ShellId>,
    /// Shell on the opposite side.
    pub back: Option<ShellId>,
    /// Boundary loops.
    pub loops: Vec<LoopId>,
    /// Carrier surface.
    pub surface: Option<SurfaceId>,
    /// Face normal relative to the surface normal.
    pub sense: Sense,
}

/// Closed boundary component of a face. Fins run anticlockwise about the face normal
/// around the material of the face.
#[derive(Clone, Debug, PartialEq)]
pub struct Loop {
    /// Owning face.
    pub face: FaceId,
    /// Fins in order.
    pub fins: Vec<FinId>,
}

/// Use of an edge by a loop (a co-edge / half-edge).
#[derive(Clone, Debug, PartialEq)]
pub struct Fin {
    /// Owning loop.
    pub owner: LoopId,
    /// Edge used.
    pub edge: EdgeId,
    /// Direction of the fin along the loop relative to the edge direction.
    pub sense: Sense,
}

/// Bounded portion of a curve. Runs from `start` to `end` along the curve direction.
/// Ring edges (closed, vertex-less) have neither.
#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
    /// Fins using this edge, in radial order around the edge.
    pub fins: Vec<FinId>,
    /// Start vertex.
    pub start: Option<VertexId>,
    /// End vertex.
    pub end: Option<VertexId>,
    /// Carrier curve.
    pub curve: Option<CurveId>,
    /// Local tolerance; `None` means exact at session precision.
    pub tolerance: Option<Tolerance>,
}

/// Point of the topology.
#[derive(Clone, Debug, PartialEq)]
pub struct Vertex {
    /// Position.
    pub point: Point3,
    /// Local tolerance; `None` means exact at session precision.
    pub tolerance: Option<Tolerance>,
}

/// Class-tagged reference to any topological entity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EntityRef {
    /// A body.
    Body(BodyId),
    /// A region.
    Region(RegionId),
    /// A shell.
    Shell(ShellId),
    /// A face.
    Face(FaceId),
    /// A loop.
    Loop(LoopId),
    /// A fin.
    Fin(FinId),
    /// An edge.
    Edge(EdgeId),
    /// A vertex.
    Vertex(VertexId),
}
