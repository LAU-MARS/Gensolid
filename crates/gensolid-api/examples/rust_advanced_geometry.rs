//! A complete Rust-only modelling round: spatial sweep, cylinder offset and cut.
use gensolid_api::{
    Error, Precision, Result, RustGeometryBackend, Session, SessionOptions, SurfaceConstruction,
    Tag,
};
use gensolid_geom::{BCurve, BSurface, nurbs_curve_ops::arc};
use gensolid_math::{Point3, Vec3};

fn checked(session: &mut Session, body: Tag, label: &str) -> Result<()> {
    let faults = session.body_check(body)?;
    if !faults.is_empty() {
        return Err(Error::internal(format!(
            "{label}: {} body faults",
            faults.len()
        )));
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

    let profile = BCurve::new(
        1,
        vec![0., 0., 1., 1.],
        vec![Point3::ORIGIN, Point3::new(3., 0., 0.)],
        None,
    )?;
    let path = BCurve::new(
        3,
        vec![0., 0., 0., 0., 1., 1., 1., 1.],
        vec![
            Point3::ORIGIN,
            Point3::new(0.1, 0.7, 0.1),
            Point3::new(-0.1, 1.3, 0.2),
            Point3::new(0., 2., 0.),
        ],
        None,
    )?;
    let sweep = SurfaceConstruction::TranslationSweep { profile, path };
    let sheet = session.body_create_sheet_backend(&backend, &sweep)?;
    let solid = session.body_create_solid_backend_prism(&backend, &sweep, 0.5)?;
    checked(&mut session, sheet, "Fixed-orientation spatial sweep sheet")?;
    checked(&mut session, solid, "Spatial sweep prism")?;

    let circular = arc(Point3::ORIGIN, 2., [-0.5, 0.5], &p)?;
    let cylinder = session.body_create_sheet_backend(
        &backend,
        &SurfaceConstruction::Extrude(circular, Vec3::Z * 1.5),
    )?;
    let offset = session.body_offset_sheet_backend(cylinder, &backend, 0.2)?;
    checked(&mut session, offset, "Normal offset of cylindrical sheet")?;

    // Cut the offset result itself, exercising the two new capabilities together.
    let plane = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        vec![
            Point3::new(-5., -5., 0.7),
            Point3::new(-5., 5., 0.7),
            Point3::new(5., -5., 0.7),
            Point3::new(5., 5., 0.7),
        ],
        None,
    )?;
    let tool = session.body_create_sheet_surface(&plane)?;
    let face = first_face(&mut session, offset)?;
    let tool_face = first_face(&mut session, tool)?;
    let split = session.body_split_face_by_sheet(offset, face, tool_face, Some(&backend))?;
    checked(&mut session, split.body, "Rational cylinder / plane cut")?;
    checked(
        &mut session,
        cylinder,
        "Original cylindrical sheet retained",
    )?;
    checked(&mut session, offset, "Original offset sheet retained")?;
    Ok(())
}
