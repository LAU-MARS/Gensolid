//! End-to-end modelling regressions and independent volume identities.
#![allow(clippy::unwrap_used)]
use gensolid_core::Precision;
use gensolid_geom::BSurface;
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_modeling::*;
use gensolid_topo::Model;
use proptest::prelude::*;

fn graph(frame: Frame, x: f64, y: f64, bend: f64) -> BSurface {
    let points = (0..3)
        .flat_map(|i| {
            (0..2).map(move |j| {
                frame.to_world(Vec3::new(
                    i as f64 * x / 2.0,
                    j as f64 * y,
                    if i == 1 { bend } else { 0.0 },
                ))
            })
        })
        .collect();
    BSurface::new(
        2,
        1,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        3,
        2,
        points,
        None,
    )
    .unwrap()
}
#[test]
fn graph_sheet_prism_split_and_rollback() {
    let mut m = Model::new(Precision::DEFAULT);
    let s = graph(Frame::WORLD, 3.0, 2.0, 0.4);
    let sheet = body_create_sheet_surface(&mut m, &s).unwrap();
    body_check_complete(&m, sheet).unwrap();
    let body = body_create_solid_surface_prism(&mut m, &s, 1.0).unwrap();
    let faces = m.body_faces(body).unwrap();
    let split = body_split_face_isoparametric(&mut m, body, faces[0], true, 0.4).unwrap();
    assert_eq!(m.body_faces(split).unwrap().len(), 7);
    body_check_complete(&m, split).unwrap();
    body_check_complete(&m, body).unwrap();
    let count = m.bodies().count();
    assert!(body_split_face_isoparametric(&mut m, body, faces[0], true, 0.0).is_err());
    assert_eq!(m.bodies().count(), count);
}
#[test]
fn sew_extracted_block() {
    let mut m = Model::new(Precision::DEFAULT);
    let b = body_create_solid_block(&mut m, 2.0, 3.0, 4.0, &Frame::WORLD).unwrap();
    let sheets = body_extract_sheets(&mut m, b).unwrap();
    let s = body_sew_sheets(&mut m, &sheets).unwrap();
    assert!(
        (body_ask_mass_properties(&m, s).unwrap().volume - 24.0).abs() < m.precision().linear()
    );
    let count = m.bodies().count();
    assert!(body_sew_sheets(&mut m, &sheets[..5]).is_err());
    assert_eq!(m.bodies().count(), count);
}
proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn random_graph_prisms(x in 0.5f64..5.0,y in 0.5f64..5.0,h in 0.1f64..3.0,bend in -0.1f64..0.1,ox in -10.0f64..10.0,ay in -0.5f64..0.5,t in 0.2f64..0.8) {
        let p=Precision::DEFAULT;
        let f=Frame::from_axis(Point3::new(ox,0.5,1.0),Vec3::new(0.2,ay,1.0),&p).unwrap();
        let s=graph(f,x,y,bend);let mut m=Model::new(p);
        let b=body_create_solid_surface_prism(&mut m,&s,h).unwrap();
        let face=m.body_faces(b).unwrap()[1];
        let result=body_split_face_isoparametric(&mut m,b,face,true,t).unwrap();
        body_check_complete(&m,result).unwrap();
    }
}

#[test]
fn planar_boolean_offset_hollow_volume() {
    let p = Precision::DEFAULT;
    let mut m = Model::new(p);
    let a = body_create_solid_block(&mut m, 4.0, 4.0, 4.0, &Frame::WORLD).unwrap();
    let f = Frame::from_axis(Point3::new(1.0, 0.0, 0.0), Vec3::Z, &p).unwrap();
    let b = body_create_solid_block(&mut m, 4.0, 4.0, 4.0, &f).unwrap();
    for (operation, expected) in [
        (BooleanOperation::Intersection, 48.0),
        (BooleanOperation::Difference, 16.0),
        (BooleanOperation::Union, 80.0),
    ] {
        let result = body_boolean_cells(&mut m, a, b, operation).unwrap();
        let mut volume = 0.0;
        for body in result {
            body_check_complete(&m, body).unwrap();
            volume += body_ask_mass_properties(&m, body).unwrap().volume;
        }
        assert!((volume - expected).abs() < p.linear() * expected);
    }
    let offset = body_offset_planar(&mut m, a, 0.5).unwrap();
    assert!(
        (body_ask_mass_properties(&m, offset).unwrap().volume - 125.0).abs() < p.linear() * 125.0
    );
    let hollow = body_hollow_planar(&mut m, a, 0.5).unwrap();
    assert_eq!(m.body(hollow).unwrap().regions.len(), 3);
    assert!(
        (body_ask_mass_properties(&m, hollow).unwrap().volume - 37.0).abs() < p.linear() * 64.0
    );
    let count = m.bodies().count();
    assert!(body_hollow_planar(&mut m, a, 2.0).is_err());
    assert_eq!(count, m.bodies().count());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn rotated_planar_csg(x in 1.0f64..4.0,y in 1.0f64..4.0,z in 1.0f64..4.0,ox in -5.0f64..5.0,ax in -0.5f64..0.5,ay in -0.5f64..0.5,shift in 0.1f64..0.8) {
        let p=Precision::DEFAULT;let mut m=Model::new(p);
        let f=Frame::from_axis(Point3::new(ox,1.0,2.0),Vec3::new(ax,ay,1.0),&p).unwrap();
        let g=Frame::new(f.origin()+f.x()*(x*shift),f.z(),f.x(),&p).unwrap();
        let a=body_create_solid_block(&mut m,x,y,z,&f).unwrap();let b=body_create_solid_block(&mut m,x,y,z,&g).unwrap();
        for (op,expected) in [(BooleanOperation::Intersection,x*y*z*(1.0-shift)),(BooleanOperation::Difference,x*y*z*shift),(BooleanOperation::Union,x*y*z*(1.0+shift))] {
            let cells=body_boolean_cells(&mut m,a,b,op).unwrap();
            let volume: f64=cells.iter().map(|&b|body_ask_mass_properties(&m,b).unwrap().volume).sum();
            prop_assert!((volume-expected).abs()<p.linear()*100.0);
        }
        let h=x.min(y).min(z)*0.1;
        let hollow=body_hollow_planar(&mut m,a,h).unwrap();let volume=body_ask_mass_properties(&m,hollow).unwrap().volume;
        prop_assert!((volume-(x*y*z-(x-2.0*h)*(y-2.0*h)*(z-2.0*h))).abs()<p.linear()*100.0);
    }
}

fn cutting_surface(y: f64, xmax: f64) -> BSurface {
    BSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        2,
        2,
        vec![
            Point3::new(-1.0, y, -2.0),
            Point3::new(-1.0, y, 3.0),
            Point3::new(xmax, y, -2.0),
            Point3::new(xmax, y, 3.0),
        ],
        None,
    )
    .unwrap()
}
#[test]
fn bounded_intersection_drives_shared_topology() {
    let mut m = Model::new(Precision::DEFAULT);
    let s = graph(Frame::WORLD, 3.0, 2.0, 0.4);
    let body = body_create_solid_surface_prism(&mut m, &s, 1.0).unwrap();
    let target = m.body_faces(body).unwrap()[1];
    let tool = body_create_sheet_surface(&mut m, &cutting_surface(0.7, 4.0)).unwrap();
    let tool_face = m.body_faces(tool).unwrap()[0];
    let result = body_split_face_by_sheet(&mut m, body, target, tool_face, None).unwrap();
    body_check_complete(&m, result.body).unwrap();
    assert_eq!(m.edge(result.intersection).unwrap().fins.len(), 2);
    for &fin in &m.edge(result.intersection).unwrap().fins {
        assert!(m.fin_geometry(fin).unwrap().is_some());
    }
    let too_small = body_create_sheet_surface(&mut m, &cutting_surface(0.7, 1.0)).unwrap();
    let face = m.body_faces(too_small).unwrap()[0];
    let count = m.bodies().count();
    assert!(body_split_face_by_sheet(&mut m, body, target, face, None).is_err());
    assert_eq!(count, m.bodies().count());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]
    #[test]
    fn planar_pipeline_in_random_frames(
        x in 1.0f64..4.0, y in 1.0f64..4.0, z in 1.0f64..4.0,
        ox in -5.0f64..5.0, ax in -0.5f64..0.5, ay in -0.5f64..0.5,
        opening in 0usize..6, boundary in 0.0f64..1.0,
    ) {
        let p = Precision::DEFAULT;
        let mut m = Model::new(p);
        let frame = Frame::from_axis(Point3::new(ox, 1., 2.), Vec3::new(ax, ay, 1.), &p).unwrap();
        let body = body_create_solid_block(&mut m, x, y, z, &frame).unwrap();
        let sheets = body_extract_sheets(&mut m, body).unwrap();
        let sewn = body_sew_sheets(&mut m, &sheets).unwrap();
        let h = x.min(y).min(z) * 0.1;
        let offset = body_offset_planar(&mut m, body, h).unwrap();
        let faces = m.body_faces(body).unwrap();
        let shell = body_shell_planar(&mut m, body, faces[opening], h).unwrap();
        let surface = m.geom().eval_surface(m.face(faces[opening]).unwrap().surface.unwrap(), 0., 0.).unwrap();
        let normal = surface.du.cross(surface.dv);
        let axes = [frame.x(), frame.y(), frame.z()];
        let dimensions = [x, y, z];
        let mut inner = [x - 2. * h, y - 2. * h, z - 2. * h];
        let axis = (0..3).max_by(|&a, &b| normal.dot(axes[a]).abs().total_cmp(&normal.dot(axes[b]).abs())).unwrap();
        inner[axis] = dimensions[axis] - h;
        let mass = |id| body_ask_mass_properties(&m, id).unwrap().volume;
        prop_assert!((mass(sewn) - x*y*z).abs() < p.linear()*100.);
        prop_assert!((mass(offset) - (x+2.*h)*(y+2.*h)*(z+2.*h)).abs() < p.linear()*100.);
        prop_assert!((mass(shell) - (x*y*z-inner.iter().product::<f64>())).abs() < p.linear()*100.);
        for result in [sewn, offset, shell] { body_check_complete(&m, result).unwrap(); }
        let count = m.bodies().count();
        prop_assert!(body_shell_planar(&mut m, body, faces[opening], p.linear()*boundary).is_err());
        prop_assert!(body_offset_planar(&mut m, body, -x.min(y).min(z)*0.5).is_err());
        prop_assert_eq!(m.bodies().count(), count);
    }

    #[test]
    fn bounded_curved_cuts_in_random_frames(
        x in 1.0f64..4.0, y in 1.0f64..4.0, fraction in 0.2f64..0.8,
        ox in -5.0f64..5.0, ax in -0.5f64..0.5, ay in -0.5f64..0.5,
        bend in -0.1f64..0.1,
    ) {
        let p = Precision::DEFAULT;
        let mut m = Model::new(p);
        let frame = Frame::from_axis(Point3::new(ox, 1., 2.), Vec3::new(ax, ay, 1.), &p).unwrap();
        let surface = graph(frame, x, y, bend);
        let body = body_create_solid_surface_prism(&mut m, &surface, 1.).unwrap();
        let cutter = BSurface::new(1,1,vec![0.,0.,1.,1.],vec![0.,0.,1.,1.],2,2,
            [Vec3::new(-1.,y*fraction,-2.),Vec3::new(-1.,y*fraction,3.),Vec3::new(x+1.,y*fraction,-2.),Vec3::new(x+1.,y*fraction,3.)]
                .into_iter().map(|q|frame.to_world(q)).collect(),None).unwrap();
        let tool = body_create_sheet_surface(&mut m, &cutter).unwrap();
        let target = m.body_faces(body).unwrap()[1];
        let tool_face = m.body_faces(tool).unwrap()[0];
        let split = body_split_face_by_sheet(&mut m, body, target, tool_face, None).unwrap();
        body_check_complete(&m, split.body).unwrap();
        prop_assert_eq!(m.body_faces(split.body).unwrap().len(), 7);
        let edge = m.edge(split.intersection).unwrap();
        prop_assert_eq!(edge.fins.len(), 2);
        for id in [edge.start.unwrap(),edge.end.unwrap()] {
            let point = frame.to_local(m.vertex(id).unwrap().point);
            prop_assert!((point.y - y*fraction).abs() < p.default_local());
        }
    }
}

#[test]
fn shell_opening_retains_constant_wall_thickness() {
    let mut m = Model::new(Precision::DEFAULT);
    let b = body_create_solid_block(&mut m, 4.0, 4.0, 4.0, &Frame::WORLD).unwrap();
    for face in m.body_faces(b).unwrap() {
        let shell = body_shell_planar(&mut m, b, face, 0.5).unwrap();
        let volume = body_ask_mass_properties(&m, shell).unwrap().volume;
        assert!((volume - (64.0 - 3.0 * 3.0 * 3.5)).abs() < m.precision().linear() * 64.0);
        assert_eq!(m.body(shell).unwrap().regions.len(), 2);
    }
}

#[test]
fn fused_boolean_boundary_regions_and_area() {
    let p = Precision::DEFAULT;
    let mut m = Model::new(p);
    let a = body_create_solid_block(&mut m, 4.0, 4.0, 4.0, &Frame::WORLD).unwrap();
    let f = Frame::from_axis(Point3::new(1.0, 0.0, 0.0), Vec3::Z, &p).unwrap();
    let b = body_create_solid_block(&mut m, 4.0, 4.0, 4.0, &f).unwrap();
    let union = body_boolean(&mut m, a, b, BooleanOperation::Union).unwrap();
    assert_eq!(union.len(), 1);
    let mass = body_ask_mass_properties(&m, union[0]).unwrap();
    assert!((mass.volume - 80.0).abs() < p.linear() * 100.0);
    assert!((mass.surface_area - 112.0).abs() < p.linear() * 100.0);
    let inner = body_create_solid_block(&mut m, 2.0, 2.0, 2.0, &Frame::WORLD).unwrap();
    let cavity = body_boolean(&mut m, a, inner, BooleanOperation::Difference).unwrap();
    assert_eq!(cavity.len(), 1);
    assert_eq!(m.body(cavity[0]).unwrap().regions.len(), 3);
    assert!(
        (body_ask_mass_properties(&m, cavity[0]).unwrap().volume - 56.0).abs() < p.linear() * 100.0
    );
    let far = Frame::from_axis(Point3::new(8.0, 0.0, 0.0), Vec3::Z, &p).unwrap();
    let b = body_create_solid_block(&mut m, 4.0, 4.0, 4.0, &far).unwrap();
    assert_eq!(
        body_boolean(&mut m, a, b, BooleanOperation::Union)
            .unwrap()
            .len(),
        2
    );
    assert!(
        body_boolean(&mut m, a, b, BooleanOperation::Intersection)
            .unwrap()
            .is_empty()
    );
    assert!(
        body_boolean(&mut m, a, a, BooleanOperation::Difference)
            .unwrap()
            .is_empty()
    );
}
proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]
    #[test]
    fn independently_rotated_boolean_conservation(ax in -0.3f64..0.3,ay in -0.3f64..0.3,bx in -0.3f64..0.3,by in -0.3f64..0.3,dx in 0.1f64..0.5) {
        let p=Precision::DEFAULT;let mut m=Model::new(p);
        let f=Frame::from_axis(Point3::new(0.0,0.0,0.0),Vec3::new(ax,ay,1.0),&p).unwrap();let g=Frame::from_axis(Point3::new(dx,0.1,0.2),Vec3::new(bx,by,1.0),&p).unwrap();
        let a=body_create_solid_block(&mut m,2.0,2.0,2.0,&f).unwrap();let b=body_create_solid_block(&mut m,2.0,2.0,2.0,&g).unwrap();
        let union=body_boolean(&mut m,a,b,BooleanOperation::Union).unwrap();let intersection=body_boolean(&mut m,a,b,BooleanOperation::Intersection).unwrap();let difference=body_boolean(&mut m,a,b,BooleanOperation::Difference).unwrap();
        let volume=|ids:&[gensolid_topo::BodyId]|ids.iter().map(|&id|body_ask_mass_properties(&m,id).unwrap().volume).sum::<f64>();
        let vi=volume(&intersection);let vu=volume(&union);let vd=volume(&difference);
        prop_assert!((vi+vu-16.0).abs()<p.linear()*100.0);prop_assert!((vd+vi-8.0).abs()<p.linear()*100.0);
        for id in union.into_iter().chain(intersection).chain(difference) {body_check_complete(&m,id).unwrap();}
    }
}

#[test]
fn existing_revolution_and_loft_can_be_subdivided() {
    use gensolid_math::Point2;
    let mut m = Model::new(Precision::DEFAULT);
    let profile = vec![
        Point2::new(1.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(2.0, 3.0),
        Point2::new(1.0, 3.0),
    ];
    let revolution =
        body_create_solid_revolution(&mut m, &profile, std::f64::consts::PI, &Frame::WORLD)
            .unwrap();
    for &face in &revolution.patches {
        let split =
            body_split_face_isoparametric(&mut m, revolution.body, face, false, 0.4).unwrap();
        body_check_complete(&m, split).unwrap();
    }
    let polygon = vec![
        Point2::new(-1.0, -1.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 1.0),
    ];
    let sections = vec![
        gensolid_topo::LoftSection {
            height: 0.0,
            profile: polygon.clone(),
        },
        gensolid_topo::LoftSection {
            height: 3.0,
            profile: polygon
                .iter()
                .map(|p| Point2::new(p.x * 1.2 + 0.2, p.y * 0.8))
                .collect(),
        },
    ];
    let loft = body_create_solid_loft_sections(&mut m, &sections, &Frame::WORLD).unwrap();
    for &face in &loft.patches {
        let split = body_split_face_isoparametric(&mut m, loft.body, face, false, 0.6).unwrap();
        body_check_complete(&m, split).unwrap();
    }
}

#[test]
fn rejects_offset_inversion_discontinuous_sheets_and_near_boundary_cuts() {
    let p = Precision::DEFAULT;
    let mut m = Model::new(p);
    let body = body_create_solid_block(&mut m, 2.0, 2.0, 10.0, &Frame::WORLD).unwrap();
    let count = m.bodies().count();
    assert!(body_offset_planar(&mut m, body, -2.0).is_err());
    assert_eq!(count, m.bodies().count());
    let points = (0..4)
        .flat_map(|i| (0..2).map(move |j| Point3::new(i as f64, j as f64, 0.0)))
        .collect();
    let discontinuous = BSurface::new(
        1,
        1,
        vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        4,
        2,
        points,
        None,
    )
    .unwrap();
    assert!(body_create_sheet_surface(&mut m, &discontinuous).is_err());
    assert_eq!(count, m.bodies().count());
    let s = graph(Frame::WORLD, 3.0, 2.0, 0.4);
    let body = body_create_solid_surface_prism(&mut m, &s, 1.0).unwrap();
    let face = m.body_faces(body).unwrap()[1];
    let count = m.bodies().count();
    assert!(body_split_face_isoparametric(&mut m, body, face, true, p.linear() * 0.01).is_err());
    assert_eq!(count, m.bodies().count());
}

fn assert_pcurves(model: &Model, body: gensolid_topo::BodyId) {
    for face in model.body_faces(body).unwrap() {
        let f = model.face(face).unwrap();
        let surface = f.surface.unwrap();
        for &lp in &f.loops {
            for &fin in &model.loop_(lp).unwrap().fins {
                let g = model.fin_geometry(fin).unwrap().unwrap();
                let edge = model.edge(model.fin(fin).unwrap().edge).unwrap();
                for fraction in [0.0, 0.2, 0.5, 0.8, 1.0] {
                    let r = g.pcurve.range();
                    let uv = g.pcurve.eval(r.lo + (r.hi - r.lo) * fraction).point;
                    let a = model
                        .geom()
                        .eval_surface(surface, uv.x, uv.y)
                        .unwrap()
                        .point;
                    let t = g.edge_range.lo + (g.edge_range.hi - g.edge_range.lo) * fraction;
                    let b = model
                        .geom()
                        .eval_curve(edge.curve.unwrap(), t)
                        .unwrap()
                        .point;
                    assert!(
                        (a - b).length() <= model.precision().linear(),
                        "UV carrier mismatch at {fin:?}"
                    );
                }
            }
        }
    }
}
#[test]
fn parameter_curves_survive_repeated_edits_and_shared_edges() {
    let mut m = Model::new(Precision::DEFAULT);
    let s = graph(Frame::WORLD, 3.0, 2.0, 0.4);
    let original = body_create_solid_surface_prism(&mut m, &s, 1.0).unwrap();
    assert_pcurves(&m, original);
    let f = m.body_faces(original).unwrap()[1];
    let a = body_split_face_isoparametric(&mut m, original, f, true, 0.4).unwrap();
    assert_pcurves(&m, a);
    let f = m.body_faces(a).unwrap()[1];
    let b = body_split_face_isoparametric(&mut m, a, f, false, 0.6).unwrap();
    assert_pcurves(&m, b);
    let block = body_create_solid_block(&mut m, 2.0, 3.0, 4.0, &Frame::WORLD).unwrap();
    for sheet in body_extract_sheets(&mut m, block).unwrap() {
        assert_pcurves(&m, sheet);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]
    #[test]
    fn non_isoparametric_sections_keep_trim_loops_and_shared_edges(
        x in 1.0f64..4.0, y in 2.0f64..4.0, ox in -5.0f64..5.0,
        ax in -0.4f64..0.4, ay in -0.4f64..0.4, bend in 0.1f64..0.4,
        tilt in 0.1f64..0.2, transpose in any::<bool>(),
    ) {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(Point3::new(ox,1.,2.),Vec3::new(ax,ay,1.),&p).unwrap();
        let mut surface = graph(frame,x,y,bend);
        if transpose {
            surface = BSurface::new(1,2,vec![0.,0.,1.,1.],vec![0.,0.,0.,1.,1.,1.],2,3,
                [0,2,4,1,3,5].into_iter().map(|i|surface.points()[i]).collect(),None).unwrap();
        }
        let mut m = Model::new(p);
        let body = body_create_solid_surface_prism(&mut m,&surface,1.).unwrap();
        let cutter = BSurface::new(1,1,vec![0.,0.,1.,1.],vec![0.,0.,1.,1.],2,2,
            [(-5.,-5.),(-5.,5.),(8.,-5.),(8.,5.)].into_iter()
                .map(|(a,b)|frame.to_world(Vec3::new(a,0.7+tilt*a+0.2*b,b))).collect(),None).unwrap();
        let tool = body_create_sheet_surface(&mut m,&cutter).unwrap();
        let tool_face = m.body_faces(tool).unwrap()[0];
        for target in m.body_faces(body).unwrap()[..2].to_vec() {
            let original_surface = m.geom().surface(m.face(target).unwrap().surface.unwrap()).unwrap().clone();
            let split = body_split_face_by_sheet(&mut m,body,target,tool_face,None).unwrap();
            body_check_complete(&m,split.body).unwrap();
            assert_pcurves(&m,split.body);
            for face in split.faces {
                let carrier = m.geom().surface(m.face(face).unwrap().surface.unwrap()).unwrap();
                // The curved trim lives on the original carrier; no rectangle
                // surface subdivision can represent this boundary.
                prop_assert_eq!(carrier,&original_surface);
            }
            let edge = m.edge(split.intersection).unwrap();
            prop_assert_eq!(edge.fins.len(),2);
            let curve = edge.curve.unwrap();
            let range = m.geom().curve_range(curve).unwrap();
            for j in 0..=16 {
                let t = range.lo+(range.hi-range.lo)*j as f64/16.;
                let point = frame.to_local(m.geom().eval_curve(curve,t).unwrap().point);
                prop_assert!((point.y-(0.7+tilt*point.x+0.2*point.z)).abs()<p.default_local());
            }
        }
        body_check_complete(&m,body).unwrap();
    }
}

#[test]
fn oblique_cubic_section_preserves_non_unit_parameter_domains() {
    let p = Precision::DEFAULT;
    let mut m = Model::new(p);
    let surface = BSurface::new(
        3,
        1,
        vec![2., 2., 2., 2., 5., 5., 5., 5.],
        vec![-3., -3., 2., 2.],
        4,
        2,
        (0..4)
            .flat_map(|i| {
                (0..2).map(move |j| Point3::new(i as f64, j as f64 * 2., [0., 0.3, 0.2, 0.][i]))
            })
            .collect(),
        None,
    )
    .unwrap();
    let body = body_create_solid_surface_prism(&mut m, &surface, 1.).unwrap();
    let cutter = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-5., -5.), (-5., 5.), (8., -5.), (8., 5.)]
            .into_iter()
            .map(|(x, z)| Point3::new(x, 0.7 + 0.2 * x + 0.2 * z, z))
            .collect(),
        None,
    )
    .unwrap();
    let tool = body_create_sheet_surface(&mut m, &cutter).unwrap();
    let target = m.body_faces(body).unwrap()[1];
    let tool_face = m.body_faces(tool).unwrap()[0];
    let result = body_split_face_by_sheet(&mut m, body, target, tool_face, None).unwrap();
    body_check_complete(&m, result.body).unwrap();
    assert_pcurves(&m, result.body);
    let count = m.bodies().count();
    // A nonrectangular trimmed face must not silently lose its trim in a later
    // operation that only supports complete parameter rectangles.
    assert!(
        body_split_face_isoparametric(&mut m, result.body, result.faces[0], false, 0.5).is_err()
    );
    assert_eq!(m.bodies().count(), count);
}

#[test]
fn backend_empty_intersection_and_incorrect_offset_are_not_certificates() {
    use gensolid_core::{Error, ErrorCode, Result};
    use gensolid_geom::backend::{
        GeometryBackend, SurfaceConstruction, SurfaceIntersectionCandidate,
    };
    struct IncorrectBackend;
    impl GeometryBackend for IncorrectBackend {
        fn precision(&self) -> Precision {
            Precision::DEFAULT
        }
        fn construct_surface(&self, _: &SurfaceConstruction) -> Result<BSurface> {
            Err(Error::unimplemented("unused construction"))
        }
        fn intersect_surfaces(
            &self,
            _: &BSurface,
            _: &BSurface,
        ) -> Result<Vec<SurfaceIntersectionCandidate>> {
            Ok(Vec::new())
        }
        fn offset_surface(&self, surface: &BSurface, _: f64) -> Result<BSurface> {
            // Deliberately return the unmodified source for a nonzero offset.
            Ok(surface.clone())
        }
    }
    let mut m = Model::new(Precision::DEFAULT);
    let source = graph(Frame::WORLD, 3., 2., 0.4);
    let body = body_create_solid_surface_prism(&mut m, &source, 1.).unwrap();
    let tool = body_create_sheet_surface(&mut m, &cutting_surface(0.7, 4.)).unwrap();
    let target = m.body_faces(body).unwrap()[1];
    let tool_face = m.body_faces(tool).unwrap()[0];
    let count = m.bodies().count();
    let error = body_split_face_by_sheet(&mut m, body, target, tool_face, Some(&IncorrectBackend))
        .unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
    assert_eq!(m.bodies().count(), count);
    let planar = body_create_sheet_surface(&mut m, &graph(Frame::WORLD, 3., 2., 0.)).unwrap();
    let count = m.bodies().count();
    let error = body_offset_sheet_backend(&mut m, planar, &IncorrectBackend, 0.5).unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
    assert_eq!(m.bodies().count(), count);
    body_check_complete(&m, body).unwrap();
    body_check_complete(&m, planar).unwrap();
}

#[test]
fn near_tangent_plane_residual_does_not_certify_spatial_section_error() {
    let p = Precision::DEFAULT;
    let mut m = Model::new(p);
    let height = p.linear() * 4.;
    let surface = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(0., 1., height),
            Point3::new(1., 0., 0.),
            Point3::new(1., 1., height + p.solver() * 2.),
        ],
        None,
    )
    .unwrap();
    let body = body_create_sheet_surface(&mut m, &surface).unwrap();
    let target = m.body_faces(body).unwrap()[0];
    let plane = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        vec![
            Point3::new(-1., -1., height * 0.5),
            Point3::new(-1., 2., height * 0.5),
            Point3::new(2., -1., height * 0.5),
            Point3::new(2., 2., height * 0.5),
        ],
        None,
    )
    .unwrap();
    let tool = body_create_sheet_surface(&mut m, &plane).unwrap();
    let tool_face = m.body_faces(tool).unwrap()[0];
    let count = m.bodies().count();
    assert!(body_split_face_by_sheet(&mut m, body, target, tool_face, None).is_err());
    assert_eq!(m.bodies().count(), count);
}
