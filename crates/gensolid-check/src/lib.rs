//! Body validity checking.
//!
//! The checker is the kernel's ground truth: every modelling operation and every imported
//! file must produce bodies that pass it, and the verifier scores generated code by it.
//! It reports faults instead of stopping at the first one, so a single run gives a full
//! picture of what is wrong.

use gensolid_core::{ErrorCode, Result, effective_tolerance};
use gensolid_math::Point3;
use gensolid_topo::*;

/// Category of a validity fault. Numeric values are part of the C API.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(i32)]
#[non_exhaustive]
pub enum FaultKind {
    /// A referenced entity or geometry does not exist.
    MissingEntity = 1,
    /// A child does not point back to its parent, or a parent does not list its child.
    BackPointer = 2,
    /// The body's region list is malformed (empty, or region 0 is not void).
    RegionLayout = 3,
    /// A loop has no fins.
    EmptyLoop = 4,
    /// Consecutive fins of a loop do not share a vertex.
    OpenLoop = 5,
    /// An edge of a solid body does not have exactly two fins.
    NonManifoldEdge = 6,
    /// The two fins of a manifold edge run in the same direction.
    FinSenses = 7,
    /// A face has no surface or an edge has no curve.
    MissingGeometry = 8,
    /// A vertex is further from an edge curve than the tolerances allow.
    VertexOffCurve = 9,
    /// A vertex is further from a face surface than its tolerance allows.
    VertexOffSurface = 10,
    /// An edge is further from a face surface than its tolerance allows.
    EdgeOffSurface = 11,
}

/// One problem found by the checker.
#[derive(Clone, Debug, PartialEq)]
pub struct Fault {
    /// Category.
    pub kind: FaultKind,
    /// Entity at fault.
    pub entity: EntityRef,
    /// Explanation.
    pub message: String,
}

/// Which checks to run.
#[derive(Clone, Copy, Debug)]
pub struct CheckOptions {
    /// Whether to run geometric consistency checks in addition to topological ones.
    pub geometry: bool,
    /// Interior samples per edge for the edge-on-surface check.
    pub edge_samples: usize,
}

impl Default for CheckOptions {
    fn default() -> Self {
        Self {
            geometry: true,
            edge_samples: 3,
        }
    }
}

/// Result of checking a body.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Report {
    /// Faults found; empty for a valid body.
    pub faults: Vec<Fault>,
    /// Geometric checks skipped because an algorithm is not available for that geometry.
    /// A body with skipped checks is not proven valid.
    pub unchecked: usize,
}

impl Report {
    /// No faults and nothing skipped.
    pub fn is_valid(&self) -> bool {
        self.faults.is_empty() && self.unchecked == 0
    }
}

struct Checker<'a> {
    model: &'a Model,
    options: CheckOptions,
    report: Report,
}

impl Checker<'_> {
    fn fault(&mut self, kind: FaultKind, entity: EntityRef, message: impl Into<String>) {
        self.report.faults.push(Fault {
            kind,
            entity,
            message: message.into(),
        });
    }

    /// Records a geometric query failure: missing geometry is a fault, a missing algorithm
    /// is an unchecked result, anything else is a fault too.
    fn geometry_error(&mut self, entity: EntityRef, error: &gensolid_core::Error) {
        match error.code() {
            ErrorCode::Unsupported | ErrorCode::Unimplemented => self.report.unchecked += 1,
            ErrorCode::NotFound => self.fault(FaultKind::MissingEntity, entity, error.message()),
            _ => self.fault(FaultKind::MissingGeometry, entity, error.to_string()),
        }
    }

    fn body(&mut self, id: BodyId) {
        let Ok(body) = self.model.body(id) else {
            return;
        };
        let regions = body.regions.clone();
        match regions.first().map(|&r| self.model.region(r)) {
            None => self.fault(
                FaultKind::RegionLayout,
                EntityRef::Body(id),
                "body has no regions",
            ),
            Some(Ok(r)) if r.solid => self.fault(
                FaultKind::RegionLayout,
                EntityRef::Body(id),
                "region 0 must be the infinite void",
            ),
            _ => {}
        }
        for region in regions {
            self.region(id, region);
        }
        if body.kind == BodyKind::Solid
            && let Ok(edges) = self.model.body_edges(id)
        {
            for edge in edges {
                self.manifold_edge(edge);
            }
        }
    }

    fn region(&mut self, body: BodyId, id: RegionId) {
        let Ok(region) = self.model.region(id) else {
            return self.fault(
                FaultKind::MissingEntity,
                EntityRef::Body(body),
                format!("region {id:?}"),
            );
        };
        if region.body != body {
            self.fault(
                FaultKind::BackPointer,
                EntityRef::Region(id),
                "region does not point to its body",
            );
        }
        for &shell in &region.shells {
            self.shell(id, shell);
        }
    }

    fn shell(&mut self, region: RegionId, id: ShellId) {
        let Ok(shell) = self.model.shell(id) else {
            return self.fault(
                FaultKind::MissingEntity,
                EntityRef::Region(region),
                format!("shell {id:?}"),
            );
        };
        if shell.region != region {
            self.fault(
                FaultKind::BackPointer,
                EntityRef::Shell(id),
                "shell does not point to its region",
            );
        }
        for &face in &shell.faces {
            self.face(id, face);
        }
    }

    fn face(&mut self, shell: ShellId, id: FaceId) {
        let Ok(face) = self.model.face(id) else {
            return self.fault(
                FaultKind::MissingEntity,
                EntityRef::Shell(shell),
                format!("face {id:?}"),
            );
        };
        if face.front != Some(shell) && face.back != Some(shell) {
            self.fault(
                FaultKind::BackPointer,
                EntityRef::Face(id),
                "face does not reference the shell listing it",
            );
        }
        // Loops are shared by both shells of a face; check them only from the front side
        // (or the back side for faces with no front shell).
        if face.front.is_some() && face.front != Some(shell) {
            return;
        }
        if face.surface.is_none() {
            self.fault(
                FaultKind::MissingGeometry,
                EntityRef::Face(id),
                "face has no surface",
            );
        }
        for &lp in &face.loops {
            self.loop_(id, lp);
        }
    }

    fn loop_(&mut self, face: FaceId, id: LoopId) {
        let Ok(lp) = self.model.loop_(id) else {
            return self.fault(
                FaultKind::MissingEntity,
                EntityRef::Face(face),
                format!("loop {id:?}"),
            );
        };
        if lp.face != face {
            self.fault(
                FaultKind::BackPointer,
                EntityRef::Loop(id),
                "loop does not point to its face",
            );
        }
        if lp.fins.is_empty() {
            return self.fault(
                FaultKind::EmptyLoop,
                EntityRef::Loop(id),
                "loop has no fins",
            );
        }
        let mut ends = Vec::with_capacity(lp.fins.len());
        for &fin in &lp.fins {
            self.fin(id, fin);
            ends.push(self.model.fin_vertices(fin).ok());
        }
        for (k, (&fin, current)) in lp.fins.iter().zip(&ends).enumerate() {
            let next = &ends[(k + 1) % ends.len()];
            if let (Some((_, end)), Some((start, _))) = (current, next)
                && end != start
            {
                self.fault(
                    FaultKind::OpenLoop,
                    EntityRef::Fin(fin),
                    "fin does not end where the next fin starts",
                );
            }
        }
        if self.options.geometry {
            for &fin in &lp.fins {
                self.fin_geometry(face, fin);
            }
        }
    }

    fn fin(&mut self, lp: LoopId, id: FinId) {
        let Ok(fin) = self.model.fin(id) else {
            return self.fault(
                FaultKind::MissingEntity,
                EntityRef::Loop(lp),
                format!("fin {id:?}"),
            );
        };
        if fin.owner != lp {
            self.fault(
                FaultKind::BackPointer,
                EntityRef::Fin(id),
                "fin does not point to its loop",
            );
        }
        match self.model.edge(fin.edge) {
            Err(_) => self.fault(
                FaultKind::MissingEntity,
                EntityRef::Fin(id),
                format!("edge {:?}", fin.edge),
            ),
            Ok(edge) => {
                if !edge.fins.contains(&id) {
                    self.fault(
                        FaultKind::BackPointer,
                        EntityRef::Edge(fin.edge),
                        "edge does not list its fin",
                    );
                }
                for v in [edge.start, edge.end].into_iter().flatten() {
                    if self.model.vertex(v).is_err() {
                        self.fault(
                            FaultKind::MissingEntity,
                            EntityRef::Edge(fin.edge),
                            format!("vertex {v:?}"),
                        );
                    }
                }
            }
        }
    }

    fn manifold_edge(&mut self, id: EdgeId) {
        let Ok(edge) = self.model.edge(id) else {
            return;
        };
        if edge.fins.len() != 2 {
            return self.fault(
                FaultKind::NonManifoldEdge,
                EntityRef::Edge(id),
                format!("solid edge has {} fins, expected 2", edge.fins.len()),
            );
        }
        let senses: Vec<_> = edge
            .fins
            .iter()
            .filter_map(|&f| self.model.fin(f).ok())
            .map(|f| f.sense)
            .collect();
        if senses.len() == 2 && senses[0] == senses[1] {
            self.fault(
                FaultKind::FinSenses,
                EntityRef::Edge(id),
                "both fins run in the same direction",
            );
        }
    }

    /// Vertices of the fin's edge lie on the edge curve and the face surface; edge
    /// interior samples lie on the face surface.
    fn fin_geometry(&mut self, face: FaceId, fin: FinId) {
        let model = self.model;
        let precision = model.precision();
        let geom = model.geom();
        let (Ok(f), Ok(fc)) = (model.fin(fin), model.face(face)) else {
            return;
        };
        let Ok(edge) = model.edge(f.edge) else { return };
        let edge_ref = EntityRef::Edge(f.edge);
        let edge_tol = effective_tolerance(edge.tolerance, precision);
        let Some(curve) = edge.curve else {
            return self.fault(FaultKind::MissingGeometry, edge_ref, "edge has no curve");
        };
        let Some(surface) = fc.surface else { return };

        let mut params = Vec::new();
        for v in [edge.start, edge.end].into_iter().flatten() {
            let Ok(vertex) = model.vertex(v) else {
                continue;
            };
            let vertex_tol = effective_tolerance(vertex.tolerance, precision);
            match geom.project_to_curve(curve, vertex.point, precision) {
                Ok(proj) => {
                    params.push(proj.param);
                    if proj.distance > vertex_tol + edge_tol {
                        self.fault(
                            FaultKind::VertexOffCurve,
                            EntityRef::Vertex(v),
                            format!("vertex is {:e} from its edge curve", proj.distance),
                        );
                    }
                }
                Err(e) => self.geometry_error(edge_ref, &e),
            }
            self.point_on_surface(
                surface,
                vertex.point,
                vertex_tol,
                FaultKind::VertexOffSurface,
                EntityRef::Vertex(v),
            );
        }

        let Ok(range) = geom.curve_range(curve) else {
            return;
        };
        let (t0, mut t1) = match params.as_slice() {
            [a, b] => (*a, *b),
            [] if range.is_bounded() => (range.lo, range.hi),
            _ => return,
        };
        if let Some(period) = range.period()
            && t1 <= t0
        {
            t1 += period;
        }
        let n = self.options.edge_samples;
        for k in 1..=n {
            let t = t0 + (t1 - t0) * k as f64 / (n + 1) as f64;
            match geom.eval_curve(curve, t) {
                Ok(d) => self.point_on_surface(
                    surface,
                    d.point,
                    edge_tol,
                    FaultKind::EdgeOffSurface,
                    edge_ref,
                ),
                Err(e) => self.geometry_error(edge_ref, &e),
            }
        }
    }

    fn point_on_surface(
        &mut self,
        surface: gensolid_geom::SurfaceId,
        p: Point3,
        tol: f64,
        kind: FaultKind,
        entity: EntityRef,
    ) {
        match self
            .model
            .geom()
            .distance_to_surface(surface, p, self.model.precision())
        {
            Ok(d) if d > tol => self.fault(
                kind,
                entity,
                format!("{d:e} from face surface, tolerance {tol:e}"),
            ),
            Ok(_) => {}
            Err(e) => self.geometry_error(entity, &e),
        }
    }
}

/// Checks a body. Fails only if the body itself does not exist; every other problem is
/// reported as a [`Fault`].
pub fn check_body(model: &Model, body: BodyId, options: &CheckOptions) -> Result<Report> {
    model.body(body)?;
    let mut checker = Checker {
        model,
        options: *options,
        report: Report::default(),
    };
    checker.body(body);
    Ok(checker.report)
}

#[cfg(test)]
mod tests {
    use gensolid_core::{Precision, Tolerance};
    use gensolid_math::{Frame, Point3, Vec3};
    use proptest::prelude::*;

    use super::*;

    fn block() -> (Model, BodyId) {
        let mut model = Model::new(Precision::DEFAULT);
        let body = make_block(&mut model, 2.0, 3.0, 4.0, &Frame::WORLD).unwrap();
        (model, body)
    }

    fn kinds(report: &Report) -> Vec<FaultKind> {
        report.faults.iter().map(|f| f.kind).collect()
    }

    #[test]
    fn block_is_valid() {
        let (model, body) = block();
        let report = check_body(&model, body, &CheckOptions::default()).unwrap();
        assert!(report.is_valid(), "{report:?}");
    }

    #[test]
    fn detects_moved_vertex_and_respects_its_tolerance() {
        let (mut model, body) = block();
        let v = model.body_vertices(body).unwrap()[0];
        let precision = *model.precision();
        let offset = precision.linear() * 100.0;
        model.raw().vertex_mut(v).unwrap().point =
            model.vertex(v).unwrap().point + Vec3::new(offset, offset, offset);
        let report = check_body(&model, body, &CheckOptions::default()).unwrap();
        assert!(
            kinds(&report).contains(&FaultKind::VertexOffCurve),
            "{report:?}"
        );
        assert!(
            kinds(&report).contains(&FaultKind::VertexOffSurface),
            "{report:?}"
        );

        // Tolerant modelling: a vertex tolerance covering the displacement makes it valid.
        let tolerance = Tolerance::new(offset * 2.0, &precision).unwrap();
        model.raw().vertex_mut(v).unwrap().tolerance = Some(tolerance);
        let report = check_body(&model, body, &CheckOptions::default()).unwrap();
        assert!(report.is_valid(), "{report:?}");
    }

    #[test]
    fn detects_broken_loop_and_back_pointer() {
        let (mut model, body) = block();
        let face = model.body_faces(body).unwrap()[0];
        let lp = model.face(face).unwrap().loops[0];
        let fin = model.loop_(lp).unwrap().fins[0];
        let mut raw = model.raw();
        let f = raw.fin_mut(fin).unwrap();
        f.sense = f.sense.reversed();
        let report = check_body(
            &model,
            body,
            &CheckOptions {
                geometry: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(kinds(&report).contains(&FaultKind::OpenLoop), "{report:?}");
        assert!(kinds(&report).contains(&FaultKind::FinSenses), "{report:?}");

        let other_loop = model
            .face(model.body_faces(body).unwrap()[1])
            .unwrap()
            .loops[0];
        model.raw().fin_mut(fin).unwrap().owner = other_loop;
        let report = check_body(
            &model,
            body,
            &CheckOptions {
                geometry: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            kinds(&report).contains(&FaultKind::BackPointer),
            "{report:?}"
        );
    }

    #[test]
    fn detects_non_manifold_edge_in_solid() {
        let (mut model, body) = block();
        let edge = model.body_edges(body).unwrap()[0];
        let extra = model.edge(edge).unwrap().fins[0];
        model.raw().edge_mut(edge).unwrap().fins.push(extra);
        let report = check_body(&model, body, &CheckOptions::default()).unwrap();
        assert!(
            kinds(&report).contains(&FaultKind::NonManifoldEdge),
            "{report:?}"
        );
    }

    proptest! {
        #[test]
        fn arbitrary_blocks_are_valid(
            x in 0.01f64..100.0, y in 0.01f64..100.0, z in 0.01f64..100.0,
            ox in -100.0f64..100.0, oy in -100.0f64..100.0, oz in -100.0f64..100.0,
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in 0.1f64..1.0,
        ) {
            let precision = Precision::DEFAULT;
            let frame = Frame::from_axis(Point3::new(ox, oy, oz), Vec3::new(ax, ay, az), &precision).unwrap();
            let mut model = Model::new(precision);
            let body = make_block(&mut model, x, y, z, &frame).unwrap();
            let report = check_body(&model, body, &CheckOptions::default()).unwrap();
            prop_assert!(report.is_valid(), "{report:?}");
        }
    }
}
