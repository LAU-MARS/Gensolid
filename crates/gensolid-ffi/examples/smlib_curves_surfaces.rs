//! Run a native surface intersection without constructing a B-rep body.
use gensolid::smlib::Backend;
use gensolid_core::{Precision, Result};
use gensolid_math::Point3;

fn main() -> Result<()> {
    let backend = Backend::new(Precision::DEFAULT)?;
    let first = backend.corner_surface([
        Point3::new(0., 0., 0.),
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
        Point3::new(1., 1., 0.),
    ])?;
    let second = backend.corner_surface([
        Point3::new(0., 0.5, -1.),
        Point3::new(1., 0.5, -1.),
        Point3::new(0., 0.5, 1.),
        Point3::new(1., 0.5, 1.),
    ])?;
    let branches = first.intersect(&second)?;
    println!(
        "{} native numerical intersection branch(es)",
        branches.len()
    );
    for branch in branches {
        let [lo, hi] = branch.curve.domain()?;
        let parameter = (lo + hi) / 2.;
        let point = branch.curve.evaluate(parameter)?.point;
        let uv = branch.first_uv.evaluate(parameter)?.point;
        let surface_point = first.evaluate([uv.x, uv.y])?.point;
        println!(
            "midpoint {point:?}; first-surface residual {}",
            point.distance(surface_point)
        );
    }
    Ok(())
}
