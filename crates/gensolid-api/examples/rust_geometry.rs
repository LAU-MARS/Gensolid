//! Real checked modelling through the built-in Rust provider, with no SMLib feature.
use gensolid_api::{
    Precision, Result, RustGeometryBackend, Session, SessionOptions, SurfaceConstruction,
};
use gensolid_geom::BCurve;
use gensolid_math::{Point3, Vec3};

fn main() -> Result<()> {
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
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
    let mut session = Session::start(SessionOptions::default());
    let request = SurfaceConstruction::Extrude(profile, Vec3::new(0., 2., 0.));
    let sheet = session.body_create_sheet_backend(&backend, &request)?;
    let solid = session.body_create_solid_backend_prism(&backend, &request, 1.)?;
    println!(
        "Rust backend: sheet faults={}, solid faults={}, solid faces={}",
        session.body_check(sheet)?.len(),
        session.body_check(solid)?.len(),
        session.body_ask_faces(solid)?.len()
    );
    Ok(())
}
