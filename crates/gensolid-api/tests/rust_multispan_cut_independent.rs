//! Independent end-to-end oracles for complete multi-span rational cuts.
//! Samples check published geometry; production acceptance must use domain bounds.
#![allow(clippy::unwrap_used)]

use gensolid_api::{Session, SessionOptions, SheetCutOptions, SurfaceConstruction, Tag};
use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_geom::{
    BCurve, BSurface, Curve, Surface,
    backend::{GeometryBackend, SurfaceIntersectionCandidate},
    nurbs_curve_ops::arc,
    rust_surface::{RustGeometryBackend, extrude_nurbs},
};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{BodyId, Model};
use proptest::prelude::*;

fn cylinder(frame: Frame, radius: f64, refined: bool, reverse: bool) -> BSurface {
    let p = Precision::DEFAULT;
    let arc = arc(Point3::ORIGIN, radius, [-1.2, 1.2], &p).unwrap();
    let mut profile = BCurve::new(
        2,
        arc.knots()
            .iter()
            .map(|t| -2. + 7. * ((t + 1.2) / 2.4))
            .collect(),
        arc.points().to_vec(),
        arc.weights().map(<[f64]>::to_vec),
    )
    .unwrap();
    if refined {
        // Binary-inexact simple knots coexist with the original double seam.
        profile = profile.insert_knots(&[-1.3, -0.2, 3.25]).unwrap();
    }
    if reverse {
        profile = profile.reversed().unwrap();
    }
    let profile = BCurve::new(
        profile.degree(),
        profile.knots().to_vec(),
        profile
            .points()
            .iter()
            .map(|q| frame.to_world(q.to_vec()))
            .collect(),
        profile.weights().map(<[f64]>::to_vec),
    )
    .unwrap();
    extrude_nurbs(&profile, frame.z() * 2., &p).unwrap()
}

fn cutter(frame: Frame, slope: f64, clipped: bool) -> BSurface {
    let ys = if clipped { [-0.25, 0.25] } else { [-8., 8.] };
    BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [-8., 8.]
            .into_iter()
            .flat_map(|x| ys.map(move |y| frame.to_world(Vec3::new(x, y, 0.9 + slope * y))))
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
fn snapshot(session: &mut Session, body: Tag) -> Snapshot {
    Snapshot {
        faces: session.body_ask_faces(body).unwrap(),
        edges: session.body_ask_edges(body).unwrap(),
        vertices: session
            .body_ask_vertices(body)
            .unwrap()
            .into_iter()
            .map(|tag| {
                (
                    tag,
                    session.vertex_ask_point(tag).unwrap().map(f64::to_bits),
                )
            })
            .collect(),
    }
}

fn shared_curve(model: &Model, body: BodyId) -> gensolid_topo::EdgeId {
    let faces = model.body_faces(body).unwrap();
    let face_edges = |face| {
        model
            .face(face)
            .unwrap()
            .loops
            .iter()
            .flat_map(|&lp| model.loop_(lp).unwrap().fins.iter())
            .map(|&fin| model.fin(fin).unwrap().edge)
            .collect::<Vec<_>>()
    };
    let first = face_edges(faces[0]);
    let second = face_edges(faces[1]);
    let shared: Vec<_> = first.into_iter().filter(|e| second.contains(e)).collect();
    assert_eq!(
        shared.len(),
        1,
        "internal spline seams must not create extra topology"
    );
    shared[0]
}

#[allow(clippy::panic)] // A wrong published carrier is an unconditional test failure.
fn inspect_section(session: &Session, body: BodyId, frame: Frame, radius: f64, slope: f64) {
    let p = Precision::DEFAULT;
    let model = session.model();
    let edge = model.edge(shared_curve(model, body)).unwrap();
    let curve = model.geom().curve(edge.curve.unwrap()).unwrap();
    let Curve::BCurve(curve) = curve else {
        panic!("curved section must retain its rational carrier")
    };
    let range = curve.range();
    let span_count = curve.knots().windows(2).filter(|k| k[0] < k[1]).count();
    assert!(span_count > 1);
    assert_eq!(
        [range.lo, range.hi].map(f64::to_bits),
        [0., span_count as f64].map(f64::to_bits)
    );
    let mut parameters: Vec<_> = (0..=31)
        .map(|i| range.lo + (range.hi - range.lo) * (i as f64 / 31.))
        .collect();
    for &seam in curve
        .knots()
        .iter()
        .filter(|&&u| u > range.lo && u < range.hi)
    {
        parameters.extend([seam.next_down(), seam, seam.next_up()]);
        assert!(
            curve
                .eval(seam.next_down())
                .point
                .distance(curve.eval(seam.next_up()).point)
                <= p.linear()
        );
    }
    for t in parameters {
        let point = curve.eval(t).point;
        let q = frame.to_local(point);
        assert!(
            (q.x.hypot(q.y) - radius).abs() <= p.linear(),
            "wrong cylinder at {t}"
        );
        assert!(
            (q.z - 0.9 - slope * q.y).abs() <= p.linear(),
            "wrong plane at {t}"
        );
        for &fin in &edge.fins {
            let fg = model.fin_geometry(fin).unwrap().unwrap();
            let fraction = (t - fg.edge_range.lo) / (fg.edge_range.hi - fg.edge_range.lo);
            let ur = fg.pcurve.range();
            let uv = fg.pcurve.eval(ur.lo + (ur.hi - ur.lo) * fraction).point;
            let a = fg.pcurve.points()[0].x;
            let b = fg.pcurve.points().last().unwrap().x;
            let face = model.loop_(model.fin(fin).unwrap().owner).unwrap().face;
            let surface = model.face(face).unwrap().surface.unwrap();
            let Surface::BSurface(carrier) = model.geom().surface(surface).unwrap() else {
                panic!("cut must retain the original NURBS surface")
            };
            let mut source_breaks = carrier.knots().0.to_vec();
            source_breaks.dedup();
            assert_eq!(source_breaks.len(), span_count + 1);
            if a > b {
                source_breaks.reverse();
            }
            assert_eq!(a.to_bits(), source_breaks[0].to_bits());
            assert_eq!(b.to_bits(), source_breaks[span_count].to_bits());
            // The edge uses integer spans. Its local parameter maps back to the
            // corresponding original, nonuniform carrier interval via the UV.
            let span = (t.floor() as usize).min(span_count - 1);
            let local = t - span as f64;
            let expected_profile =
                source_breaks[span] + (source_breaks[span + 1] - source_breaks[span]) * local;
            assert!(
                (uv.x - expected_profile).abs() <= p.solver(),
                "profile parameter lost at seam {t}"
            );
            assert!(
                model
                    .geom()
                    .eval_surface(surface, uv.x, uv.y)
                    .unwrap()
                    .point
                    .distance(point)
                    <= p.linear(),
                "UV lost at {t}"
            );
        }
    }
}

fn offset_cut(
    frame: Frame,
    radius: f64,
    distance: f64,
    refined: bool,
    reverse: bool,
    backend: Option<&dyn GeometryBackend>,
) {
    let p = Precision::DEFAULT;
    let native = RustGeometryBackend::new(p);
    let mut session = Session::start(SessionOptions::default());
    let original = session
        .body_create_sheet_surface(&cylinder(frame, radius, refined, reverse))
        .unwrap();
    let original_state = snapshot(&mut session, original);
    let offset = session
        .body_offset_sheet_backend(original, &native, distance)
        .unwrap();
    let before = snapshot(&mut session, offset);
    let tool = session
        .body_create_sheet_surface(&cutter(frame, 0.12, false))
        .unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let old = session.body_split_face_by_sheet_with_options(
        offset,
        before.faces[0],
        tool_face,
        backend,
        SheetCutOptions {
            allow_rational_oblique: true,
        },
    );
    assert_eq!(old.unwrap_err().code(), ErrorCode::Unimplemented);
    let previous: Vec<_> = session.model().bodies().collect();
    let result = session
        .body_split_face_by_sheet_multispan(offset, before.faces[0], tool_face, backend)
        .unwrap();
    assert_eq!(session.body_ask_faces(result.body).unwrap().len(), 2);
    assert_eq!(session.body_ask_edges(result.body).unwrap().len(), 7);
    assert_eq!(session.body_ask_vertices(result.body).unwrap().len(), 6);
    assert!(session.body_check(result.body).unwrap().is_empty());
    assert_eq!(snapshot(&mut session, original), original_state);
    assert_eq!(snapshot(&mut session, offset), before);
    assert!(session.body_check(original).unwrap().is_empty());
    assert!(session.body_check(offset).unwrap().is_empty());
    let created = session
        .model()
        .bodies()
        .find(|b| !previous.contains(b))
        .unwrap();
    inspect_section(
        &session,
        created,
        frame,
        radius + if reverse { -distance } else { distance },
        0.12,
    );
}

#[test]
fn real_offset_then_multispan_oblique_cut_retains_one_continuous_shared_edge() {
    let native = RustGeometryBackend::new(Precision::DEFAULT);
    for refined in [false, true] {
        for reverse in [false, true] {
            offset_cut(Frame::WORLD, 2., 0.25, refined, reverse, Some(&native));
            offset_cut(Frame::WORLD, 2., 0.25, refined, reverse, None);
        }
    }
}

#[test]
fn negative_sense_solid_face_cut_propagates_both_endpoint_splits_to_adjacent_faces() {
    let p = Precision::DEFAULT;
    let mut session = Session::start(SessionOptions::default());
    // Use the already supported two-span prism constructor. Binary-inexact
    // refined prism sewing has its own explicit, transactional refusal below.
    let source = cylinder(Frame::WORLD, 2., false, false);
    let solid = session
        .body_create_solid_surface_prism(&source, 0.3)
        .unwrap();
    let before = snapshot(&mut session, solid);
    // The bottom face has negative sense. Its crossing edge is traversed in
    // the opposite direction from the same carrier on a standalone sheet.
    let face = before.faces[0];
    let tool = session
        .body_create_sheet_surface(&cutter(Frame::WORLD, 0.12, false))
        .unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let previous: Vec<_> = session.model().bodies().collect();
    let result = session
        .body_split_face_by_sheet_multispan(
            solid,
            face,
            tool_face,
            Some(&RustGeometryBackend::new(p)),
        )
        .unwrap();
    assert_eq!(session.body_ask_faces(result.body).unwrap().len(), 7);
    assert_eq!(session.body_ask_edges(result.body).unwrap().len(), 15);
    assert_eq!(session.body_ask_vertices(result.body).unwrap().len(), 10);
    assert!(session.body_check(result.body).unwrap().is_empty());
    assert_eq!(snapshot(&mut session, solid), before);
    let model = session.model();
    let created = model.bodies().find(|b| !previous.contains(b)).unwrap();
    for id in model.body_edges(created).unwrap() {
        let edge = model.edge(id).unwrap();
        assert_eq!(edge.fins.len(), 2);
        for &fin in &edge.fins {
            let fg = model.fin_geometry(fin).unwrap().unwrap();
            let face = model.loop_(model.fin(fin).unwrap().owner).unwrap().face;
            let surface = model.face(face).unwrap().surface.unwrap();
            for i in 0..=11 {
                let fraction = i as f64 / 11.;
                let ur = fg.pcurve.range();
                let uv = fg.pcurve.eval(ur.lo + (ur.hi - ur.lo) * fraction).point;
                let t = fg.edge_range.lo + (fg.edge_range.hi - fg.edge_range.lo) * fraction;
                let point = model
                    .geom()
                    .eval_curve(edge.curve.unwrap(), t)
                    .unwrap()
                    .point;
                assert!(
                    model
                        .geom()
                        .eval_surface(surface, uv.x, uv.y)
                        .unwrap()
                        .point
                        .distance(point)
                        <= p.linear()
                );
            }
        }
    }
    inspect_section(&session, created, Frame::WORLD, 2., 0.12);
}

#[test]
fn refined_prism_sewing_refuses_incompatible_boundary_bases_transactionally() {
    let source = cylinder(Frame::WORLD, 2., true, false);
    let mut session = Session::start(SessionOptions::default());
    let mut control = session.clone();
    let result = session.body_create_solid_surface_prism(&source, 0.3);
    assert_eq!(result.unwrap_err().code(), ErrorCode::Unimplemented);
    assert_eq!(session.model().bodies().count(), 0);
    assert_eq!(
        session.body_create_solid_block(1., 2., 3., None).unwrap(),
        control.body_create_solid_block(1., 2., 3., None).unwrap()
    );
}

#[derive(Clone, Copy)]
enum Attack {
    FalseGeometry,
    OneSpan,
    ReverseOrder,
    Empty,
}
struct Adversary(Attack);
impl GeometryBackend for Adversary {
    fn precision(&self) -> Precision {
        Precision::DEFAULT
    }
    fn construct_surface(&self, _: &SurfaceConstruction) -> Result<BSurface> {
        Err(Error::unimplemented("cut-only provider"))
    }
    fn offset_surface(&self, _: &BSurface, _: f64) -> Result<BSurface> {
        Err(Error::unimplemented("cut-only provider"))
    }
    fn intersect_surfaces(
        &self,
        a: &BSurface,
        b: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>> {
        if matches!(self.0, Attack::Empty) {
            return Ok(Vec::new());
        }
        let mut candidates = RustGeometryBackend::new(self.precision()).intersect_surfaces(a, b)?;
        match self.0 {
            Attack::OneSpan => candidates.truncate(1),
            Attack::ReverseOrder => candidates.reverse(),
            Attack::FalseGeometry => {
                for candidate in &mut candidates {
                    let wrong = BCurve::new(
                        1,
                        vec![-2., -2., 5., 5.],
                        vec![Point3::new(40., 50., 60.), Point3::new(45., 55., 65.)],
                        None,
                    )?;
                    candidate.curve = wrong.clone();
                    candidate.first_uv = wrong.clone();
                    candidate.second_uv = wrong;
                    candidate.reported_deviation = 0.;
                }
            }
            Attack::Empty => {}
        }
        Ok(candidates)
    }
}

#[test]
fn incomplete_reordered_and_false_provider_proposals_cannot_replace_complete_geometry() {
    for attack in [Attack::OneSpan, Attack::ReverseOrder, Attack::FalseGeometry] {
        offset_cut(
            Frame::WORLD,
            2.,
            0.25,
            true,
            false,
            Some(&Adversary(attack)),
        );
    }
}

fn rejected_cut(source: &BSurface, tool: &BSurface, backend: Option<&dyn GeometryBackend>) {
    let mut session = Session::start(SessionOptions::default());
    let body = session.body_create_sheet_surface(source).unwrap();
    let before = snapshot(&mut session, body);
    let cutter = session.body_create_sheet_surface(tool).unwrap();
    let face = session.body_ask_faces(cutter).unwrap()[0];
    let count = session.model().bodies().count();
    let mut control = session.clone();
    let result = session.body_split_face_by_sheet_multispan(body, before.faces[0], face, backend);
    assert_eq!(result.unwrap_err().code(), ErrorCode::Unimplemented);
    assert_eq!(session.model().bodies().count(), count);
    assert_eq!(snapshot(&mut session, body), before);
    assert!(session.body_check(body).unwrap().is_empty());
    let expected = control.body_create_solid_block(1., 2., 3., None).unwrap();
    let next = session.body_create_solid_block(1., 2., 3., None).unwrap();
    assert_eq!(next, expected);
    assert_eq!(
        snapshot(&mut session, next),
        snapshot(&mut control, expected)
    );
}

#[test]
fn missing_provider_result_partial_branch_and_clipped_tool_roll_back_every_tag() {
    let source = cylinder(Frame::WORLD, 2., true, false);
    rejected_cut(
        &source,
        &cutter(Frame::WORLD, 0.12, false),
        Some(&Adversary(Attack::Empty)),
    );
    rejected_cut(&source, &cutter(Frame::WORLD, 0.8, false), None);
    rejected_cut(&source, &cutter(Frame::WORLD, 0.12, true), None);
}

#[test]
fn tangent_and_near_tangent_cases_remain_explicit_in_the_new_entry() {
    let p = Precision::DEFAULT;
    let source = cylinder(Frame::WORLD, 2., true, false);
    let tangent = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        vec![
            Point3::new(2., -8., -8.),
            Point3::new(2., -8., 8.),
            Point3::new(2., 8., -8.),
            Point3::new(2., 8., 8.),
        ],
        None,
    )
    .unwrap();
    rejected_cut(&source, &tangent, None);
    let epsilon = p.linear() * 10_000.;
    let profile = BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![
            Point3::ORIGIN,
            Point3::new(epsilon * 0.1, 0., 1.5),
            Point3::new(epsilon * 0.4, 0., 3.),
        ],
        Some(vec![1., 0.7, 1.]),
    )
    .unwrap()
    .insert_knots(&[0.21, 0.65, 0.65])
    .unwrap();
    let source = extrude_nurbs(&profile, Vec3::Y * 2., &p).unwrap();
    let plane = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-5., -5.), (-5., 5.), (5., -5.), (5., 5.)]
            .into_iter()
            .map(|(y, z)| Point3::new(epsilon * (1. - y), y, z))
            .collect(),
        None,
    )
    .unwrap();
    rejected_cut(&source, &plane, None);
}

#[test]
fn discontinuous_surface_is_rejected_before_cutting_without_allocating_tags() {
    let p = Precision::DEFAULT;
    let curve = BCurve::new(
        2,
        vec![0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1.],
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(0.5, 0., 0.1),
            Point3::new(1., 0., 0.),
            Point3::new(2., 0., 0.),
            Point3::new(2.5, 0., 0.1),
            Point3::new(3., 0., 0.),
        ],
        Some(vec![1., 0.8, 1., 1., 0.8, 1.]),
    )
    .unwrap();
    let source = extrude_nurbs(&curve, Vec3::Y * 2., &p).unwrap();
    let mut session = Session::start(SessionOptions::default());
    let mut control = session.clone();
    assert_eq!(
        session
            .body_create_sheet_surface(&source)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    assert_eq!(session.model().bodies().count(), 0);
    assert_eq!(
        session.body_create_solid_block(1., 2., 3., None).unwrap(),
        control.body_create_solid_block(1., 2., 3., None).unwrap()
    );
}

#[test]
fn common_homogeneous_weight_scales_do_not_change_multispan_sections() {
    let p = Precision::DEFAULT;
    let source = cylinder(Frame::WORLD, 2., true, false);
    for factor in [2.0_f64.powi(-400), 1., 2.0_f64.powi(400)] {
        let scaled = BSurface::new(
            2,
            1,
            source.knots().0.to_vec(),
            source.knots().1.to_vec(),
            source.counts().0,
            2,
            source.points().to_vec(),
            Some(
                source
                    .weights()
                    .unwrap()
                    .iter()
                    .map(|w| w * factor)
                    .collect(),
            ),
        )
        .unwrap();
        let mut session = Session::start(SessionOptions::default());
        let body = session.body_create_sheet_surface(&scaled).unwrap();
        let face = session.body_ask_faces(body).unwrap()[0];
        let tool = session
            .body_create_sheet_surface(&cutter(Frame::WORLD, 0.12, false))
            .unwrap();
        let tool_face = session.body_ask_faces(tool).unwrap()[0];
        let previous: Vec<_> = session.model().bodies().collect();
        let result = session
            .body_split_face_by_sheet_multispan(
                body,
                face,
                tool_face,
                Some(&RustGeometryBackend::new(p)),
            )
            .unwrap();
        assert!(session.body_check(result.body).unwrap().is_empty());
        let created = session
            .model()
            .bodies()
            .find(|b| !previous.contains(b))
            .unwrap();
        inspect_section(&session, created, Frame::WORLD, 2., 0.12);
    }
}

#[test]
fn unresolved_weights_or_coarse_parameters_cannot_publish_a_multispan_cut() {
    let source = cylinder(Frame::WORLD, 2., false, false);
    let mut weights = source.weights().unwrap().to_vec();
    weights[0] = f64::from_bits(1);
    weights[1] = f64::from_bits(1);
    let bad_weights = BSurface::new(
        2,
        1,
        source.knots().0.to_vec(),
        source.knots().1.to_vec(),
        source.counts().0,
        2,
        source.points().to_vec(),
        Some(weights),
    )
    .unwrap();
    let large = 1_000_000_000_000_000.;
    let bad_domain = BSurface::new(
        2,
        1,
        source
            .knots()
            .0
            .iter()
            .map(|k| large + (k + 2.) / 7.)
            .collect(),
        source.knots().1.to_vec(),
        source.counts().0,
        2,
        source.points().to_vec(),
        source.weights().map(<[f64]>::to_vec),
    )
    .unwrap();
    let mut session = Session::start(SessionOptions::default());
    let mut control = session.clone();
    // The underflowing weights fail at creation, before any cut can run.
    assert!(session.body_create_sheet_surface(&bad_weights).is_err());
    assert_eq!(session.model().bodies().count(), 0);
    assert_eq!(
        session.body_create_solid_block(1., 2., 3., None).unwrap(),
        control.body_create_solid_block(1., 2., 3., None).unwrap()
    );
    // The large parameter origin does enter the model, but its coarse parameter
    // resolution cannot represent this oblique branch to the spatial tolerance.
    // Integer edge parameters must not hide the original surface's quantization.
    rejected_cut(&bad_domain, &cutter(Frame::WORLD, 0.12, false), None);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn offset_and_multispan_cut_compose_in_random_frames(
        radius in 1.0f64..3., distance in -0.2f64..0.4,
        tx in -3.0f64..3., tilt in -0.4f64..0.4, reverse in any::<bool>(),
    ) {
        let p=Precision::DEFAULT;
        let frame=Frame::from_axis(Point3::new(tx,1.,-2.),Vec3::new(tilt,0.3,1.),&p).unwrap();
        offset_cut(frame,radius,distance,true,reverse,Some(&RustGeometryBackend::new(p)));
    }
}
