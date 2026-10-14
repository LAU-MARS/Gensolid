//! Native geometry candidates enter checked Gensolid modelling operations.
use gensolid::smlib::Backend;
use gensolid_api::{Axis2, BooleanOperation, Session, SessionOptions, SurfaceConstruction, Tag};
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{BCurve, BSurface};
use gensolid_math::{Point3, Vec3};

fn face(session: &mut Session, body: Tag, index: usize) -> Result<Tag> {
    session
        .body_ask_faces(body)?
        .get(index)
        .copied()
        .ok_or_else(|| Error::internal("example body lacks the expected face"))
}

fn report(session: &mut Session, name: &str, body: Tag) -> Result<()> {
    let faults = session.body_check(body)?;
    if !faults.is_empty() {
        return Err(Error::internal(format!("{name}: {faults:?}")));
    }
    println!(
        "{name}: {} faces, check passed",
        session.body_ask_faces(body)?.len()
    );
    Ok(())
}

fn rectangle(points: Vec<Point3>) -> Result<BSurface> {
    BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        points,
        None,
    )
}

fn main() -> Result<()> {
    let backend = Backend::new(Precision::DEFAULT)?;
    let mut s = Session::start(SessionOptions::default());
    let profile = BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(1.5, 0., 0.4),
            Point3::new(3., 0., 0.),
        ],
        None,
    )?;
    let request = SurfaceConstruction::Extrude(profile, Vec3::new(0., 2., 0.));
    let sheet = s.body_create_sheet_backend(&backend, &request)?;
    report(&mut s, "SMLib curved sheet", sheet)?;
    let prism = s.body_create_solid_backend_prism(&backend, &request, 1.)?;
    report(&mut s, "Directional curved prism", prism)?;
    let target = face(&mut s, prism, 1)?;
    let cutter = rectangle(vec![
        Point3::new(-1., 0.7, -2.),
        Point3::new(-1., 0.7, 3.),
        Point3::new(4., 0.7, -2.),
        Point3::new(4., 0.7, 3.),
    ])?;
    let tool = s.body_create_sheet_surface(&cutter)?;
    let tool_face = face(&mut s, tool, 0)?;
    let split = s.body_split_face_by_sheet(prism, target, tool_face, Some(&backend))?;
    report(&mut s, "Bounded intersection face split", split.body)?;
    println!(
        "Replacement face tags {:?}; shared edge tag {}",
        split.faces, split.intersection
    );

    let oblique = rectangle(
        [(-5., -5.), (-5., 5.), (8., -5.), (8., 5.)]
            .into_iter()
            .map(|(x, z)| Point3::new(x, 0.7 + 0.2 * x + 0.2 * z, z))
            .collect(),
    )?;
    let tool = s.body_create_sheet_surface(&oblique)?;
    let tool_face = face(&mut s, tool, 0)?;
    let trimmed = s.body_split_face_by_sheet(prism, target, tool_face, Some(&backend))?;
    report(
        &mut s,
        "Oblique section with non-isoparametric UV loops",
        trimmed.body,
    )?;

    let plane = rectangle(vec![
        Point3::new(0., 0., 0.),
        Point3::new(0., 2., 0.),
        Point3::new(3., 0., 0.),
        Point3::new(3., 2., 0.),
    ])?;
    let plane_sheet = s.body_create_sheet_surface(&plane)?;
    let offset_sheet = s.body_offset_sheet_backend(plane_sheet, &backend, 0.5)?;
    drop(backend);
    report(&mut s, "Certified planar sheet offset", offset_sheet)?;
    report(&mut s, "Split body after backend release", split.body)?;

    // These operations currently accept convex planar solids; the curved prism
    // above is deliberately not presented as a supported Boolean/shell input.
    let a = s.body_create_solid_block(4., 4., 4., None)?;
    let basis = Axis2 {
        location: [1., 0., 0.],
        axis: [0., 0., 1.],
        ref_direction: [1., 0., 0.],
    };
    let b = s.body_create_solid_block(4., 4., 4., Some(&basis))?;
    let sheets = s.body_extract_sheets(a)?;
    let sewn = s.body_sew_sheets(&sheets)?;
    report(&mut s, "Sewn planar solid", sewn)?;
    for result in s.body_boolean(a, b, BooleanOperation::Union)? {
        report(&mut s, "Fused planar union", result)?;
        println!(
            "Union volume: {} (expected 80)",
            s.body_ask_mass_properties(result)?.volume
        );
    }
    let offset = s.body_offset_planar(a, 0.5)?;
    report(&mut s, "Planar solid offset", offset)?;
    let hollow = s.body_hollow_planar(a, 0.5)?;
    report(&mut s, "Closed hollow with finite cavity", hollow)?;
    let opening = face(&mut s, a, 0)?;
    let shell = s.body_shell_planar(a, opening, 0.5)?;
    report(&mut s, "Open planar shell", shell)?;
    Ok(())
}
