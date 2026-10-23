//! Certified rational oblique cuts and multi-span cylindrical offsets in Rust.
use gensolid_api::{
    Error, Precision, Result, RustGeometryBackend, Session, SessionOptions, SheetCutOptions,
    SurfaceConstruction, Tag,
};
use gensolid_geom::{BSurface, nurbs_curve_ops::arc};
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

fn cutting_sheet(session: &mut Session, slope: f64) -> Result<Tag> {
    let plane = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-5., -5.), (-5., 5.), (5., -5.), (5., 5.)]
            .into_iter()
            .map(|(x, y)| Point3::new(x, y, 0.65 + slope * x))
            .collect(),
        None,
    )?;
    session.body_create_sheet_surface(&plane)
}

fn main() -> Result<()> {
    let p = Precision::DEFAULT;
    let backend = RustGeometryBackend::new(p);
    let mut session = Session::start(SessionOptions::default());
    let source = session.body_create_sheet_backend(
        &backend,
        &SurfaceConstruction::Extrude(arc(Point3::ORIGIN, 2., [-0.5, 0.5], &p)?, Vec3::Z * 1.5),
    )?;
    let offset = session.body_offset_sheet_backend(source, &backend, 0.15)?;
    let tool = cutting_sheet(&mut session, 0.1)?;
    let face = first_face(&mut session, offset)?;
    let tool_face = first_face(&mut session, tool)?;
    let cut = session.body_split_face_by_sheet_with_options(
        offset,
        face,
        tool_face,
        Some(&backend),
        SheetCutOptions {
            allow_rational_oblique: true,
        },
    )?;
    checked(
        &mut session,
        cut.body,
        "Offset cylinder / oblique plane cut",
    )?;
    checked(
        &mut session,
        source,
        "Original single-span cylinder retained",
    )?;
    checked(&mut session, offset, "Original single-span offset retained")?;

    // This arc crosses an internal NURBS knot and uses two rational spans.
    let multi = session.body_create_sheet_backend(
        &backend,
        &SurfaceConstruction::Extrude(arc(Point3::ORIGIN, 2., [-0.8, 0.8], &p)?, Vec3::Z * 1.5),
    )?;
    let multi_offset = session.body_offset_sheet_backend(multi, &backend, -0.25)?;
    checked(
        &mut session,
        multi_offset,
        "Two-span cylindrical normal offset",
    )?;
    let tool = cutting_sheet(&mut session, 0.)?;
    let face = first_face(&mut session, multi_offset)?;
    let tool_face = first_face(&mut session, tool)?;
    let cut = session.body_split_face_by_sheet(multi_offset, face, tool_face, Some(&backend))?;
    checked(
        &mut session,
        cut.body,
        "Multi-span offset / isoparametric plane cut",
    )?;
    checked(&mut session, multi, "Original multi-span cylinder retained")?;
    checked(
        &mut session,
        multi_offset,
        "Original multi-span offset retained",
    )?;
    Ok(())
}
