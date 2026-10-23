//! Rational oblique cuts must change actual checked topology, not merely return candidates.
#![allow(clippy::unwrap_used)]
use gensolid_api::{Session, SessionOptions, SheetCutOptions, SurfaceConstruction, Tag};
use gensolid_core::{ErrorCode, Precision, Result};
use gensolid_geom::{
    BCurve, BSurface,
    backend::{GeometryBackend, SurfaceIntersectionCandidate},
    rust_surface::RustGeometryBackend,
};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

const OPTIONS: SheetCutOptions = SheetCutOptions {
    allow_rational_oblique: true,
};
fn profile(frame: Frame, width: f64, weight: f64) -> BCurve {
    BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        [(0., 0.), (0.5, 0.1), (1., 0.)]
            .into_iter()
            .map(|(x, z)| frame.to_world(Vec3::new(x * width, 0., z)))
            .collect(),
        Some(vec![1., weight, 1.]),
    )
    .unwrap()
}
fn tool(frame: Frame, width: f64, depth: f64, slope: f64, partial: bool) -> BSurface {
    let xs = if partial {
        [width * 0.25, width * 0.75]
    } else {
        [-5., width + 5.]
    };
    BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        xs.into_iter()
            .flat_map(|x| {
                [-5., 5.]
                    .map(move |z| frame.to_world(Vec3::new(x, depth * 0.35 + slope * (x + z), z)))
            })
            .collect(),
        None,
    )
    .unwrap()
}
#[derive(Debug, PartialEq)]
struct Snapshot {
    faces: Vec<Tag>,
    edges: Vec<Tag>,
    vertices: Vec<(Tag, [u64; 3])>,
}
fn snapshot(s: &mut Session, body: Tag) -> Snapshot {
    Snapshot {
        faces: s.body_ask_faces(body).unwrap(),
        edges: s.body_ask_edges(body).unwrap(),
        vertices: s
            .body_ask_vertices(body)
            .unwrap()
            .into_iter()
            .map(|v| (v, s.vertex_ask_point(v).unwrap().map(f64::to_bits)))
            .collect(),
    }
}
fn solid(frame: Frame, width: f64, depth: f64, weight: f64) -> (Session, Tag) {
    let mut s = Session::start(SessionOptions::default());
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
    let request = SurfaceConstruction::Extrude(profile(frame, width, weight), frame.y() * depth);
    let body = s
        .body_create_solid_backend_prism(&backend, &request, 0.4)
        .unwrap();
    (s, body)
}
fn verify_cut(
    frame: Frame,
    width: f64,
    depth: f64,
    weight: f64,
    slope: f64,
    backend: Option<&dyn GeometryBackend>,
) {
    let (mut s, body) = solid(frame, width, depth, weight);
    let before = snapshot(&mut s, body);
    let target = s.body_ask_faces(body).unwrap()[1];
    let cutter = s
        .body_create_sheet_surface(&tool(frame, width, depth, slope, false))
        .unwrap();
    let tool_face = s.body_ask_faces(cutter).unwrap()[0];
    let previous: Vec<_> = s.model().bodies().collect();
    let result = s
        .body_split_face_by_sheet_with_options(body, target, tool_face, backend, OPTIONS)
        .unwrap();
    assert!(s.body_check(result.body).unwrap().is_empty());
    assert!(s.body_check(body).unwrap().is_empty());
    assert_eq!(s.body_ask_faces(result.body).unwrap().len(), 7);
    assert_eq!(s.body_ask_edges(result.body).unwrap().len(), 15);
    assert_eq!(s.body_ask_vertices(result.body).unwrap().len(), 10);
    assert_eq!(snapshot(&mut s, body), before);
    let model = s.model();
    let created = model.bodies().find(|b| !previous.contains(b)).unwrap();
    for edge in model.body_edges(created).unwrap() {
        assert_eq!(model.edge(edge).unwrap().fins.len(), 2);
    }
    // The new edge and every propagated adjacent subdivision retain paired UV.
    for face in model.body_faces(created).unwrap() {
        let surface = model.face(face).unwrap().surface.unwrap();
        for &lp in &model.face(face).unwrap().loops {
            for &fin in &model.loop_(lp).unwrap().fins {
                let geometry = model.fin_geometry(fin).unwrap().unwrap();
                let edge = model.edge(model.fin(fin).unwrap().edge).unwrap();
                for i in 0..=12 {
                    let fraction = i as f64 / 12.;
                    let range = geometry.pcurve.range();
                    let uv = geometry
                        .pcurve
                        .eval(range.lo + (range.hi - range.lo) * fraction)
                        .point;
                    let t = geometry.edge_range.lo
                        + (geometry.edge_range.hi - geometry.edge_range.lo) * fraction;
                    let xyz = model
                        .geom()
                        .eval_curve(edge.curve.unwrap(), t)
                        .unwrap()
                        .point;
                    assert!(
                        xyz.distance(
                            model
                                .geom()
                                .eval_surface(surface, uv.x, uv.y)
                                .unwrap()
                                .point
                        ) <= Precision::DEFAULT.linear()
                    );
                }
            }
        }
    }
}

#[test]
fn explicit_rational_oblique_cut_preserves_solid_adjacency_with_and_without_provider() {
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
    verify_cut(Frame::WORLD, 3., 2., 0.8, 0.12, Some(&backend));
    verify_cut(Frame::WORLD, 3., 2., 0.8, 0.12, None);
}

struct FalseProvider {
    empty: bool,
}
impl GeometryBackend for FalseProvider {
    fn precision(&self) -> Precision {
        Precision::DEFAULT
    }
    fn construct_surface(&self, request: &SurfaceConstruction) -> Result<BSurface> {
        RustGeometryBackend::new(self.precision()).construct_surface(request)
    }
    fn offset_surface(&self, s: &BSurface, d: f64) -> Result<BSurface> {
        RustGeometryBackend::new(self.precision()).offset_surface(s, d)
    }
    fn intersect_surfaces(
        &self,
        a: &BSurface,
        b: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>> {
        if self.empty {
            return Ok(vec![]);
        }
        let mut results = RustGeometryBackend::new(self.precision()).intersect_surfaces(a, b)?;
        for result in &mut results {
            result.curve = BCurve::new(
                1,
                vec![0., 0., 1., 1.],
                vec![Point3::new(100., 100., 100.), Point3::new(101., 100., 100.)],
                None,
            )?;
            result.first_uv = BCurve::new(
                1,
                vec![0., 0., 1., 1.],
                vec![Point3::new(0., 0., 0.), Point3::new(1., 0., 0.)],
                None,
            )?;
            result.reported_deviation = 0.;
        }
        Ok(results)
    }
}
#[test]
fn false_provider_geometry_and_zero_deviation_cannot_certify_the_cut() {
    verify_cut(
        Frame::WORLD,
        3.,
        2.,
        0.8,
        0.12,
        Some(&FalseProvider { empty: false }),
    );
}

#[test]
fn legacy_empty_and_partial_operations_keep_existing_state() {
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
    let (mut s, body) = solid(Frame::WORLD, 3., 2., 0.8);
    let before = snapshot(&mut s, body);
    let target = before.faces[1];
    for partial in [false, true] {
        let cutter = s
            .body_create_sheet_surface(&tool(Frame::WORLD, 3., 2., 0.12, partial))
            .unwrap();
        let face = s.body_ask_faces(cutter).unwrap()[0];
        let count = s.model().bodies().count();
        assert_eq!(
            s.body_split_face_by_sheet(body, target, face, Some(&backend))
                .unwrap_err()
                .code(),
            ErrorCode::Unimplemented
        );
        assert_eq!(
            s.body_split_face_by_sheet_with_options(
                body,
                target,
                face,
                Some(&FalseProvider { empty: true }),
                OPTIONS
            )
            .unwrap_err()
            .code(),
            ErrorCode::Unimplemented
        );
        if partial {
            assert!(
                s.body_split_face_by_sheet_with_options(body, target, face, None, OPTIONS)
                    .is_err()
            );
        }
        assert_eq!(s.model().bodies().count(), count);
        assert_eq!(snapshot(&mut s, body), before);
        assert!(s.body_check(body).unwrap().is_empty());
    }
}

#[test]
fn nonmonotone_uv_control_ordinates_do_not_reject_a_monotone_rational_identity() {
    let p = Precision::DEFAULT;
    let curve = BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(1.5, 0., 0.1),
            Point3::new(3., 0., 0.),
        ],
        Some(vec![1., 0.49, 1.]),
    )
    .unwrap();
    let backend = RustGeometryBackend::new(p);
    let request = SurfaceConstruction::Extrude(curve, Vec3::Y * 2.);
    let source = backend.construct_surface(&request).unwrap();
    let candidate = backend
        .intersect_surfaces(&source, &tool(Frame::WORLD, 3., 2., 0.12, false))
        .unwrap()
        .remove(0);
    assert!(
        candidate
            .first_uv
            .points()
            .windows(2)
            .any(|pair| pair[1].x < pair[0].x)
    );
    let mut s = Session::start(SessionOptions::default());
    let sheet = s.body_create_sheet_backend(&backend, &request).unwrap();
    let cutter = s
        .body_create_sheet_surface(&tool(Frame::WORLD, 3., 2., 0.12, false))
        .unwrap();
    let face = s.body_ask_faces(sheet).unwrap()[0];
    let tool_face = s.body_ask_faces(cutter).unwrap()[0];
    let result = s
        .body_split_face_by_sheet_with_options(sheet, face, tool_face, Some(&backend), OPTIONS)
        .unwrap();
    assert!(s.body_check(result.body).unwrap().is_empty());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn rational_oblique_solid_cuts_pass_in_random_frames(
        width in 2.0f64..4.,depth in 1.5f64..3.,weight in 0.6f64..1.4,
        slope in 0.03f64..0.12,ox in -5.0f64..5.,ax in -0.3f64..0.3,ay in -0.3f64..0.3,
    ) {
        let frame=Frame::from_axis(Point3::new(ox,1.,2.),Vec3::new(ax,ay,1.),&Precision::DEFAULT).unwrap();
        verify_cut(frame,width,depth,weight,slope,Some(&RustGeometryBackend::new(Precision::DEFAULT)));
    }

    #[test]
    fn near_boundary_crossings_are_refused_without_publishing_topology(factor in 0.01f64..0.9) {
        let p=Precision::DEFAULT;
        let (mut s,body)=solid(Frame::WORLD,3.,2.,0.8);
        let before=snapshot(&mut s,body);
        let cutter=BSurface::new(1,1,vec![0.,0.,1.,1.],vec![0.,0.,1.,1.],2,2,
            [(-5.,-5.),(-5.,5.),(8.,-5.),(8.,5.)].into_iter()
                .map(|(x,z)|Point3::new(x,p.linear()*factor+0.1*x,z)).collect(),None).unwrap();
        let cutter=s.body_create_sheet_surface(&cutter).unwrap();
        let tool_face=s.body_ask_faces(cutter).unwrap()[0];
        let count=s.model().bodies().count();
        let error=s.body_split_face_by_sheet_with_options(body,before.faces[1],tool_face,None,OPTIONS).unwrap_err();
        prop_assert_eq!(error.code(),ErrorCode::Unimplemented);
        prop_assert_eq!(s.model().bodies().count(),count);
        prop_assert_eq!(snapshot(&mut s,body),before);
        prop_assert!(s.body_check(body).unwrap().is_empty());
    }
}
