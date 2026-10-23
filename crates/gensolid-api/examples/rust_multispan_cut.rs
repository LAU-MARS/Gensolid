//! A continuous multi-span cylinder, normal offset, and certified oblique cut.
use gensolid_api::{
    Error, Precision, Result, RustGeometryBackend, Session, SessionOptions, SurfaceConstruction,
    Tag,
};
use gensolid_geom::{BSurface, nurbs_curve_ops::arc};
use gensolid_math::{Point3, Vec3};

fn checked(session: &mut Session, body: Tag, label: &str) -> Result<()> {
    let faults = session.body_check(body)?;
    if !faults.is_empty() {
        return Err(Error::internal(format!("{label}: {} faults", faults.len())));
    }
    println!(
        "{label}: {} faces, {} edges, {} vertices, 0 faults",
        session.body_ask_faces(body)?.len(),
        session.body_ask_edges(body)?.len(),
        session.body_ask_vertices(body)?.len()
    );
    Ok(())
}

fn first_face(session: &mut Session, body: Tag) -> Result<Tag> {
    session
        .body_ask_faces(body)?
        .first()
        .copied()
        .ok_or_else(|| Error::internal("example body has no face"))
}

fn main() -> Result<()> {
    let p = Precision::DEFAULT;
    let backend = RustGeometryBackend::new(p);
    let mut session = Session::start(SessionOptions::default());
    // The original double knot and three simple knots yield five unequal spans.
    let profile = arc(Point3::ORIGIN, 2., [-1.2, 1.2], &p)?.insert_knots(&[-0.8, -0.2, 0.6])?;
    let spans = profile.knots().windows(2).filter(|k| k[0] < k[1]).count();
    println!("Profile: {spans} rational spans");
    let source = session.body_create_sheet_backend(
        &backend,
        &SurfaceConstruction::Extrude(profile, Vec3::Z * 1.5),
    )?;
    let offset = session.body_offset_sheet_backend(source, &backend, 0.25)?;
    let plane = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-5., -5.), (-5., 5.), (5., -5.), (5., 5.)]
            .into_iter()
            .map(|(x, y)| Point3::new(x, y, 0.65 + 0.1 * y))
            .collect(),
        None,
    )?;
    let tool = session.body_create_sheet_surface(&plane)?;
    let target_face = first_face(&mut session, offset)?;
    let tool_face = first_face(&mut session, tool)?;
    let cut = session.body_split_face_by_sheet_multispan(
        offset,
        target_face,
        tool_face,
        Some(&backend),
    )?;
    checked(&mut session, cut.body, "Multi-span offset / oblique cut")?;
    checked(&mut session, source, "Original cylinder retained")?;
    checked(&mut session, offset, "Original offset retained")?;
    Ok(())
}
