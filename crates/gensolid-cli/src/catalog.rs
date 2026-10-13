//! Deterministic development cases with closed-form reference branch counts.
//! Random scales and common rigid frames preserve the independently known answers.

use gensolid_core::{Precision, Result};
use gensolid_geom::{
    BCurve, BSurface, Circle, Cone, Curve, Cylinder, Ellipse, IntersectionCase, IntersectionScore,
    Line, Plane, Sphere, Surface, Torus, intersection_score,
};
use gensolid_math::{Frame, Point3, Vec3};

struct Generator(u64);

impl Generator {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn signed(&mut self) -> f64 {
        2.0 * self.unit() - 1.0
    }
}

pub(crate) fn score(seed: u64, frames: usize, precision: &Precision) -> Result<IntersectionScore> {
    let mut random = Generator(seed);
    let mut cases = Vec::new();
    for _ in 0..frames {
        let origin = Point3::new(
            random.signed() * 25.0,
            random.signed() * 25.0,
            random.signed() * 25.0,
        );
        let axis = Vec3::new(random.signed(), random.signed(), random.signed());
        let frame = Frame::from_axis(origin, axis, precision)?;
        let radius = 0.1 + 10.0 * random.unit();
        let sphere = Surface::Sphere(Sphere::new(frame, radius, precision)?);
        let cylinder = Surface::Cylinder(Cylinder::new(frame, radius * 0.5, precision)?);
        let equator = Surface::Plane(Plane::new(frame));
        cases.push(IntersectionCase::Surfaces {
            a: sphere.clone(),
            b: equator.clone(),
            expected_branches: 1,
        });
        cases.push(IntersectionCase::Surfaces {
            a: sphere.clone(),
            b: cylinder.clone(),
            expected_branches: 2,
        });
        cases.push(IntersectionCase::Surfaces {
            a: cylinder.clone(),
            b: equator,
            expected_branches: 1,
        });
        let side = Frame::new(
            origin + frame.x() * (radius * 0.25),
            frame.x(),
            frame.z(),
            precision,
        )?;
        cases.push(IntersectionCase::Surfaces {
            a: cylinder.clone(),
            b: Surface::Plane(Plane::new(side)),
            expected_branches: 2,
        });
        let offset = Frame::new(origin + frame.z() * radius, frame.z(), frame.x(), precision)?;
        let other_sphere = Surface::Sphere(Sphere::new(offset, radius, precision)?);
        cases.push(IntersectionCase::Surfaces {
            a: sphere.clone(),
            b: other_sphere,
            expected_branches: 1,
        });
        cases.push(IntersectionCase::CurveSurface {
            curve: Curve::Line(Line::new(origin, frame.x(), precision)?),
            surface: sphere.clone(),
            expected_branches: 2,
        });
        // The equatorial circle is contained in the sphere: one overlap component.
        cases.push(IntersectionCase::CurveSurface {
            curve: Curve::Circle(Circle::new(frame, radius, precision)?),
            surface: sphere.clone(),
            expected_branches: 1,
        });
        let circle_frame = Frame::new(
            origin + frame.x() * (radius * 0.75),
            frame.z(),
            frame.x(),
            precision,
        )?;
        cases.push(IntersectionCase::CurveSurface {
            curve: Curve::Circle(Circle::new(circle_frame, radius, precision)?),
            surface: cylinder,
            expected_branches: 2,
        });
        cases.push(IntersectionCase::CurveSurface {
            curve: Curve::Ellipse(Ellipse::new(frame, radius * 1.5, radius * 0.5, precision)?),
            surface: sphere.clone(),
            expected_branches: 4,
        });
        let segment = BCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![origin - frame.z() * radius, origin + frame.z() * radius],
            Some(vec![1.0, 2.0]),
        )?;
        cases.push(IntersectionCase::CurveSurface {
            curve: Curve::BCurve(segment),
            surface: Surface::Plane(Plane::new(frame)),
            expected_branches: 1,
        });
        let torus = Surface::Torus(Torus::new(frame, radius * 3.0, radius, precision)?);
        // A meridional plane cuts a ring torus into two circles.
        let meridian = Frame::new(origin, frame.x(), frame.z(), precision)?;
        cases.push(IntersectionCase::Surfaces {
            a: torus.clone(),
            b: Surface::Plane(Plane::new(meridian)),
            expected_branches: 2,
        });
        // The sphere lies inside the torus hole because r < R - tube_radius.
        cases.push(IntersectionCase::Surfaces {
            a: torus.clone(),
            b: sphere,
            expected_branches: 0,
        });
        let off_axis = Frame::new(
            origin + frame.x() * (radius * 0.25),
            frame.z(),
            frame.x(),
            precision,
        )?;
        // Even after this offset, sphere extent 1.25r stays inside torus hole 2r.
        cases.push(IntersectionCase::Surfaces {
            a: torus,
            b: Surface::Sphere(Sphere::new(off_axis, radius, precision)?),
            expected_branches: 0,
        });
        // For equal weights the plane residual is r(1 - 8t + 8t²): two roots
        // strictly inside (0,1). Higher-degree NURBS intersection remains visible.
        let quadratic = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                origin - frame.x() * radius + frame.z() * radius,
                origin - frame.z() * (radius * 3.0),
                origin + frame.x() * radius + frame.z() * radius,
            ],
            Some(vec![1.0, 1.0, 1.0]),
        )?;
        cases.push(IntersectionCase::CurveSurface {
            curve: Curve::BCurve(quadratic),
            surface: Surface::Plane(Plane::new(frame)),
            expected_branches: 2,
        });
    }
    Ok(intersection_score(&cases, precision))
}

/// A separate fixed workload for classes beyond the original fourteen configurations.
/// Keep this score separate so implementing a baseline class improves a comparable run.
pub(crate) fn extended_score(
    seed: u64,
    frames: usize,
    precision: &Precision,
) -> Result<IntersectionScore> {
    let mut random = Generator(seed);
    let mut cases = Vec::new();
    for _ in 0..frames {
        let origin = Point3::new(
            random.signed() * 25.0,
            random.signed() * 25.0,
            random.signed() * 25.0,
        );
        let frame = Frame::from_axis(
            origin,
            Vec3::new(random.signed(), random.signed(), random.signed()),
            precision,
        )?;
        let r = 0.1 + 10.0 * random.unit();
        // Degree elevation preserves the quadratic residual r(1 - 8t + 8t²).
        let cubic = BCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![
                origin - frame.x() * r + frame.z() * r,
                origin - frame.x() * (r / 3.0) - frame.z() * (5.0 * r / 3.0),
                origin + frame.x() * (r / 3.0) - frame.z() * (5.0 * r / 3.0),
                origin + frame.x() * r + frame.z() * r,
            ],
            None,
        )?;
        cases.push(IntersectionCase::CurveSurface {
            curve: Curve::BCurve(cubic),
            surface: Surface::Plane(Plane::new(frame)),
            expected_branches: 2,
        });
        let patch = BSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            2,
            2,
            vec![
                origin - frame.x() * r - frame.y() * r,
                origin - frame.x() * r + frame.y() * r,
                origin + frame.x() * r - frame.y() * r,
                origin + frame.x() * r + frame.y() * r,
            ],
            None,
        )?;
        let vertical = Frame::new(origin, frame.x(), frame.z(), precision)?;
        cases.push(IntersectionCase::Surfaces {
            a: Surface::BSurface(patch),
            b: Surface::Plane(Plane::new(vertical)),
            expected_branches: 1,
        });
        let cone = Surface::Cone(Cone::new(
            frame,
            0.0,
            std::f64::consts::FRAC_PI_4,
            precision,
        )?);
        // Small tilt is below the cone generator angle, so the elevated plane cuts
        // one closed ellipse on the positive nappe.
        let oblique = Frame::new(
            origin + frame.z() * (2.0 * r),
            frame.z() + frame.x() * 0.1,
            frame.y(),
            precision,
        )?;
        cases.push(IntersectionCase::Surfaces {
            a: cone,
            b: Surface::Plane(Plane::new(oblique)),
            expected_branches: 1,
        });
    }
    Ok(intersection_score(&cases, precision))
}
