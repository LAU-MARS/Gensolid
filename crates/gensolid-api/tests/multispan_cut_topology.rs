//! Certified multi-span cuts retain original carriers and real shared topology.
#![allow(clippy::unwrap_used)]

use gensolid_api::{Session, SessionOptions, SheetCutOptions, Tag};
use gensolid_core::Precision;
use gensolid_geom::{
    BCurve, BSurface, Surface,
    backend::GeometryBackend,
    rust_surface::{RustGeometryBackend, extrude_nurbs},
};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::Sense;
use proptest::prelude::*;

#[derive(Clone, Copy)]
enum KnotLayout {
    Dyadic,
    Nonbinary,
    Translated,
}

fn source(frame: Frame, width: f64, depth: f64, linear_u: bool, layout: KnotLayout) -> BSurface {
    let (domain, knots) = match layout {
        KnotLayout::Dyadic => ([0.0, 1.0], vec![0.25, 0.5, 0.5, 0.75]),
        KnotLayout::Nonbinary => ([0.0, 1.0], vec![0.3, 0.6, 0.6, 0.8]),
        KnotLayout::Translated => ([-2.0, 5.0], vec![-1.3, -0.2, -0.2, 3.25]),
    };
    let profile = BCurve::new(
        2,
        vec![
            domain[0], domain[0], domain[0], domain[1], domain[1], domain[1],
        ],
        [(0.0, 0.0), (0.5, 0.12), (1.0, 0.0)]
            .into_iter()
            .map(|(x, z)| frame.to_world(Vec3::new(x * width, 0.0, z)))
            .collect(),
        Some(vec![1.0, 0.8, 1.0]),
    )
    .unwrap()
    .insert_knots(&knots)
    .unwrap();
    let s = extrude_nurbs(&profile, frame.y() * depth, &Precision::DEFAULT).unwrap();
    if !linear_u {
        return s;
    }
    let points = (0..2)
        .flat_map(|j| (0..s.counts().0).map(move |i| i * 2 + j))
        .map(|i| s.points()[i])
        .collect();
    let weights = (0..2)
        .flat_map(|j| (0..s.counts().0).map(move |i| i * 2 + j))
        .map(|i| s.weights().unwrap()[i])
        .collect();
    BSurface::new(
        1,
        2,
        s.knots().1.to_vec(),
        s.knots().0.to_vec(),
        2,
        s.counts().0,
        points,
        Some(weights),
    )
    .unwrap()
}

fn tool(frame: Frame, width: f64, depth: f64, slope: f64) -> BSurface {
    BSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        2,
        2,
        [-5.0, width + 5.0]
            .into_iter()
            .flat_map(|x| {
                [-5.0, 5.0]
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

// The two explicit panics below report a test failure if a required NURBS
// carrier is replaced with another geometry type; this is not library code.
#[allow(clippy::too_many_arguments, clippy::panic)]
fn checked_cut(
    frame: Frame,
    width: f64,
    depth: f64,
    slope: f64,
    linear_u: bool,
    solid: bool,
    target_index: usize,
    layout: KnotLayout,
    provider: bool,
) {
    let p = Precision::DEFAULT;
    let surface = source(frame, width, depth, linear_u, layout);
    let mut session = Session::start(SessionOptions::default());
    let body = if solid {
        session
            .body_create_solid_surface_prism(&surface, 0.4)
            .unwrap()
    } else {
        session.body_create_sheet_surface(&surface).unwrap()
    };
    let original_id = session.model().bodies().next().unwrap();
    let original_target = session.model().body_faces(original_id).unwrap()[target_index];
    let original_carrier = session
        .model()
        .geom()
        .surface(
            session
                .model()
                .face(original_target)
                .unwrap()
                .surface
                .unwrap(),
        )
        .unwrap()
        .clone();
    let original_sense = session.model().face(original_target).unwrap().sense;
    let before = snapshot(&mut session, body);
    let target = before.faces[target_index];
    let cutter = session
        .body_create_sheet_surface(&tool(frame, width, depth, slope))
        .unwrap();
    let tool_face = session.body_ask_faces(cutter).unwrap()[0];
    let old_bodies = session.model().bodies().collect::<Vec<_>>();
    let backend = RustGeometryBackend::new(p);
    let provider = if provider {
        Some(&backend as &dyn GeometryBackend)
    } else {
        None
    };
    // The old public literal and single-span capability keep their contract.
    assert!(
        session
            .body_split_face_by_sheet_with_options(
                body,
                target,
                tool_face,
                provider,
                SheetCutOptions {
                    allow_rational_oblique: true
                }
            )
            .is_err()
    );
    assert_eq!(snapshot(&mut session, body), before);
    assert_eq!(session.model().bodies().count(), old_bodies.len());
    let cut = session
        .body_split_face_by_sheet_multispan(body, target, tool_face, provider)
        .unwrap();
    assert_eq!(session.model().bodies().count(), old_bodies.len() + 1);
    assert!(session.body_check(body).unwrap().is_empty());
    assert!(session.body_check(cut.body).unwrap().is_empty());
    assert_eq!(
        session.body_ask_faces(cut.body).unwrap().len(),
        if solid { 7 } else { 2 }
    );
    assert_eq!(
        session.body_ask_edges(cut.body).unwrap().len(),
        if solid { 15 } else { 7 }
    );
    assert_eq!(
        session.body_ask_vertices(cut.body).unwrap().len(),
        if solid { 10 } else { 6 }
    );
    assert!(
        session
            .body_ask_edges(cut.body)
            .unwrap()
            .contains(&cut.intersection)
    );
    assert_eq!(snapshot(&mut session, body), before);
    let model = session.model();
    let result = model.bodies().find(|b| !old_bodies.contains(b)).unwrap();
    let faces = model.body_faces(result).unwrap();
    let replacements = [faces[target_index], faces[target_index + 1]];
    let edges = |face| {
        model
            .face(face)
            .unwrap()
            .loops
            .iter()
            .flat_map(|&lp| model.loop_(lp).unwrap().fins.iter())
            .map(|&f| model.fin(f).unwrap().edge)
            .collect::<Vec<_>>()
    };
    let first = edges(replacements[0]);
    let second = edges(replacements[1]);
    let shared = first
        .into_iter()
        .filter(|e| second.contains(e))
        .collect::<Vec<_>>();
    assert_eq!(shared.len(), 1);
    let edge = model.edge(shared[0]).unwrap();
    assert_eq!(edge.fins.len(), 2);
    assert_ne!(
        model.fin(edge.fins[0]).unwrap().sense,
        model.fin(edge.fins[1]).unwrap().sense
    );
    let gensolid_geom::Curve::BCurve(curve) = model.geom().curve(edge.curve.unwrap()).unwrap()
    else {
        panic!("multi-span crossing lost its NURBS carrier")
    };
    let spans = curve
        .knots()
        .windows(2)
        .filter(|w| w[0] < w[1])
        .collect::<Vec<_>>();
    assert!(spans.len() > 1);
    for face in replacements {
        let face = model.face(face).unwrap();
        assert_eq!(face.sense, original_sense);
        assert_eq!(
            model.geom().surface(face.surface.unwrap()).unwrap(),
            &original_carrier
        );
        for &lp in &face.loops {
            for &fin in &model.loop_(lp).unwrap().fins {
                if model.fin(fin).unwrap().edge != shared[0] {
                    continue;
                }
                let geometry = model.fin_geometry(fin).unwrap().unwrap();
                let uv_range = geometry.pcurve.range();
                for span in &spans {
                    for i in 0..=8 {
                        let t = span[0] + (span[1] - span[0]) * (i as f64 / 8.0);
                        let fraction = (t - geometry.edge_range.lo)
                            / (geometry.edge_range.hi - geometry.edge_range.lo);
                        let uv = geometry
                            .pcurve
                            .eval(uv_range.lo + (uv_range.hi - uv_range.lo) * fraction)
                            .point;
                        let xyz = curve.eval(t).point;
                        assert!(
                            xyz.distance(
                                model
                                    .geom()
                                    .eval_surface(face.surface.unwrap(), uv.x, uv.y)
                                    .unwrap()
                                    .point
                            ) <= p.linear()
                        );
                    }
                }
            }
        }
    }
    if solid {
        assert!(
            model.body_edges(result).unwrap().into_iter().all(|e| model
                .edge(e)
                .unwrap()
                .fins
                .len()
                == 2)
        );
        assert_eq!(
            original_sense,
            if target_index == 0 {
                Sense::Negative
            } else {
                Sense::Positive
            }
        );
    }
    let Surface::BSurface(stored) = original_carrier else {
        panic!("source is not NURBS")
    };
    assert_eq!(stored.knots(), surface.knots());
}

#[test]
fn multispan_sheet_cuts_preserve_carriers_and_both_parameter_directions() {
    for linear_u in [false, true] {
        for provider in [false, true] {
            checked_cut(
                Frame::WORLD,
                3.0,
                2.0,
                0.12,
                linear_u,
                false,
                0,
                KnotLayout::Dyadic,
                provider,
            );
        }
    }
}

#[test]
fn multispan_solid_cuts_preserve_reversed_faces_and_neighbor_adjacency() {
    for linear_u in [false, true] {
        for target in [0, 1] {
            checked_cut(
                Frame::WORLD,
                3.0,
                2.0,
                0.12,
                linear_u,
                true,
                target,
                KnotLayout::Dyadic,
                true,
            );
        }
    }
}

#[test]
fn nonbinary_seams_and_translated_domains_retain_both_fin_parameter_curves() {
    for linear_u in [false, true] {
        for layout in [KnotLayout::Nonbinary, KnotLayout::Translated] {
            checked_cut(
                Frame::WORLD,
                3.0,
                2.0,
                0.12,
                linear_u,
                false,
                0,
                layout,
                true,
            );
        }
    }
}

#[test]
fn the_ruled_graph_fast_path_rejects_a_fold_across_profile_spans() {
    let p = Precision::DEFAULT;
    let profile = BCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![
            Point3::ORIGIN,
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        Some(vec![1.0, 0.8, 1.0]),
    )
    .unwrap()
    .insert_knots(&[0.3, 0.7, 0.7])
    .unwrap();
    let source = extrude_nurbs(&profile, Vec3::Y * 2.0, &p).unwrap();
    let mut session = Session::start(SessionOptions::default());
    let original = session
        .body_create_solid_block(1.0, 1.0, 1.0, None)
        .unwrap();
    let before = snapshot(&mut session, original);
    let count = session.model().bodies().count();
    assert!(session.body_create_sheet_surface(&source).is_err());
    assert!(
        session
            .body_create_solid_surface_prism(&source, 0.4)
            .is_err()
    );
    assert_eq!(session.model().bodies().count(), count);
    assert_eq!(snapshot(&mut session, original), before);
    assert!(session.body_check(original).unwrap().is_empty());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn multispan_topology_survives_random_frames_sizes_and_face_orientation(
        width in 1.0f64..4.0,depth in 1.5f64..3.0,slope in 0.04f64..0.15,
        tilt in -2.0f64..2.0,translation in -3.0f64..3.0,
        linear_u in any::<bool>(),negative in any::<bool>(),
    ) {
        let frame=Frame::from_axis(Point3::new(translation,-translation,translation),Vec3::new(tilt,1.0,1.0),&Precision::DEFAULT)?;
        checked_cut(frame,width,depth,slope,linear_u,true,usize::from(!negative),KnotLayout::Dyadic,true);
    }
}
