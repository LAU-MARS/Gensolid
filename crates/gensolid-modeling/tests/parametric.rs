//! Checked curved modelling regression tests.
use gensolid_core::Precision;
use gensolid_math::{Frame, Point2, Point3, Vec3};
use gensolid_modeling::{
    body_check_complete, body_create_solid_loft_sections, body_create_solid_revolution,
};
use gensolid_topo::{LoftSection, Model};
use proptest::prelude::*;
fn profile(r: f64, h: f64) -> Vec<Point2> {
    vec![
        Point2::new(r, -h),
        Point2::new(2.0 * r, -h),
        Point2::new(2.0 * r, h),
        Point2::new(r, h),
    ]
}
#[test]
fn full_and_partial_revolution_keep_rational_surfaces_and_closed_incidence() {
    let p = Precision::DEFAULT;
    for angle in [std::f64::consts::TAU, std::f64::consts::PI, 0.7] {
        let mut m = Model::new(p);
        let result =
            body_create_solid_revolution(&mut m, &profile(1.0, 1.0), angle, &Frame::WORLD).unwrap();
        body_check_complete(&m, result.body).unwrap();
        for e in m.body_edges(result.body).unwrap() {
            assert_eq!(m.edge(e).unwrap().fins.len(), 2);
        }
        assert!(!result.patches.is_empty());
    }
}
proptest! {
 #![proptest_config(ProptestConfig::with_cases(12))]
 #[test]
 fn parametric_solids_pass_complete_checks_in_random_frames(r in 0.2f64..3.0,h in 0.2f64..3.0,angle in 0.3f64..6.0,ox in -5.0f64..5.0,oy in -5.0f64..5.0,oz in -5.0f64..5.0,ax in -1.0f64..1.0,ay in -1.0f64..1.0){
  let p=Precision::DEFAULT;let frame=Frame::from_axis(Point3::new(ox,oy,oz),Vec3::new(ax,ay,1.0),&p).unwrap();
  let mut m=Model::new(p);let a=body_create_solid_revolution(&mut m,&profile(r,h),angle,&frame).unwrap();body_check_complete(&m,a.body).unwrap();
  let sections=(0..3).map(|j|{let a=j as f64*0.12;LoftSection{height:j as f64*h,profile:(0..5).map(|i|{let t=a+std::f64::consts::TAU*i as f64/5.0;Point2::new(r*t.cos()*(1.0+j as f64*0.2),r*t.sin()+j as f64*0.1)}).collect()}}).collect::<Vec<_>>();
  let a=body_create_solid_loft_sections(&mut m,&sections,&frame).unwrap();body_check_complete(&m,a.body).unwrap();
  prop_assert_eq!(m.body_faces(a.body).unwrap().len(),12);
 }
}
#[test]
fn rejected_axis_and_foldover_do_not_publish_bodies() {
    let p = Precision::DEFAULT;
    let mut m = Model::new(p);
    assert!(
        body_create_solid_revolution(&mut m, &profile(p.linear() * 0.5, 1.0), 1.0, &Frame::WORLD)
            .is_err()
    );
    let a = vec![
        Point2::new(-1.0, -1.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 1.0),
    ];
    let mut b = a.clone();
    b.rotate_left(2);
    assert!(
        body_create_solid_loft_sections(
            &mut m,
            &[
                LoftSection {
                    height: 0.0,
                    profile: a
                },
                LoftSection {
                    height: 1.0,
                    profile: b
                }
            ],
            &Frame::WORLD
        )
        .is_err()
    );
    assert_eq!(m.bodies().count(), 0);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]
    #[test]
    fn near_degenerate_profiles_are_rejected_without_publishing(factor in 0.1f64..0.99, ox in -5.0f64..5.0, ax in -1.0f64..1.0) {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(Point3::new(ox,0.0,0.0),Vec3::new(ax,1.0,1.0),&p).unwrap();
        let mut m = Model::new(p);
        prop_assert!(body_create_solid_revolution(&mut m,&profile(p.linear()*factor,1.0),1.0,&frame).is_err());
        let polygon = profile(1.0,1.0);
        let sections = [LoftSection {height:0.0,profile:polygon.clone()},LoftSection{height:p.linear()*factor,profile:polygon}];
        prop_assert!(body_create_solid_loft_sections(&mut m,&sections,&frame).is_err());
        prop_assert_eq!(m.bodies().count(),0);
    }
}
