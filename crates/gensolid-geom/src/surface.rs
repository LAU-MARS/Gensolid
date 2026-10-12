//! Surfaces.
//!
//! Parameterisations put the natural outward normal along `du × dv`: away from the axis
//! for cylinders, cones, spheres and tori, along `z` for planes.

use std::f64::consts::{FRAC_PI_2, TAU};

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Frame, Point3, Vec3};

use crate::nurbs::BSurface;
use crate::param::ParamRange;
use crate::store::CurveId;

/// Point and first partial derivatives of a surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDerivs {
    /// Position.
    pub point: Point3,
    /// Partial derivative along `u`.
    pub du: Vec3,
    /// Partial derivative along `v`.
    pub dv: Vec3,
}

impl SurfaceDerivs {
    /// Unit normal `du × dv`, or `None` at a singular point.
    pub fn normal(&self, precision: &Precision) -> Option<Vec3> {
        let n = self.du.cross(self.dv);
        let scale = self.du.length() * self.dv.length();
        n.normalized(precision.angular() * scale)
    }
}

/// Plane through the frame origin spanned by `x` (u) and `y` (v); normal `z`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    frame: Frame,
}

impl Plane {
    /// Creates a plane.
    pub const fn new(frame: Frame) -> Self {
        Self { frame }
    }

    /// Frame.
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// Evaluates at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> SurfaceDerivs {
        let f = &self.frame;
        SurfaceDerivs {
            point: f.origin() + f.x() * u + f.y() * v,
            du: f.x(),
            dv: f.y(),
        }
    }

    fn distance(&self, p: Point3) -> f64 {
        self.frame.to_local(p).z.abs()
    }
}

fn radial(frame: &Frame, u: f64) -> (Vec3, Vec3) {
    let (s, c) = u.sin_cos();
    (frame.x() * c + frame.y() * s, frame.y() * c - frame.x() * s)
}

fn check_radius(radius: f64, precision: &Precision, what: &str) -> Result<()> {
    if radius.is_finite() && radius > precision.linear() {
        Ok(())
    } else {
        Err(Error::degenerate(format!(
            "{what} must exceed linear precision"
        )))
    }
}

/// Cylinder of given radius about the frame `z` axis; `u` angle, `v` height.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cylinder {
    frame: Frame,
    radius: f64,
}

impl Cylinder {
    /// Creates a cylinder.
    pub fn new(frame: Frame, radius: f64, precision: &Precision) -> Result<Self> {
        check_radius(radius, precision, "cylinder radius")?;
        Ok(Self { frame, radius })
    }

    /// Frame.
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// Radius.
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// Evaluates at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> SurfaceDerivs {
        let (r, t) = radial(&self.frame, u);
        SurfaceDerivs {
            point: self.frame.origin() + r * self.radius + self.frame.z() * v,
            du: t * self.radius,
            dv: self.frame.z(),
        }
    }

    fn distance(&self, p: Point3) -> f64 {
        let l = self.frame.to_local(p);
        (l.x.hypot(l.y) - self.radius).abs()
    }
}

/// Cone about the frame `z` axis with `radius` at `v = 0` and half-angle `semi_angle`;
/// the radius grows with `v` as `radius + v * tan(semi_angle)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cone {
    frame: Frame,
    radius: f64,
    semi_angle: f64,
}

impl Cone {
    /// Creates a cone with `radius >= 0` and `semi_angle` strictly inside `(0, π/2)`.
    pub fn new(frame: Frame, radius: f64, semi_angle: f64, precision: &Precision) -> Result<Self> {
        if !(radius.is_finite() && radius >= 0.0) {
            return Err(Error::invalid_argument(
                "cone radius must be finite and non-negative",
            ));
        }
        let a = precision.angular();
        if !(semi_angle > a && semi_angle < FRAC_PI_2 - a) {
            return Err(Error::degenerate(
                "cone semi-angle must lie strictly inside (0, π/2)",
            ));
        }
        Ok(Self {
            frame,
            radius,
            semi_angle,
        })
    }

    /// Frame.
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// Radius at `v = 0` and half-angle.
    pub const fn shape(&self) -> (f64, f64) {
        (self.radius, self.semi_angle)
    }

    /// Evaluates at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> SurfaceDerivs {
        let (r, t) = radial(&self.frame, u);
        let tan = self.semi_angle.tan();
        let rv = self.radius + v * tan;
        SurfaceDerivs {
            point: self.frame.origin() + r * rv + self.frame.z() * v,
            du: t * rv,
            dv: r * tan + self.frame.z(),
        }
    }

    /// Distance to the full double cone: in the half-plane through the axis, the surface is
    /// the generator line and its mirror image in the axis.
    fn distance(&self, p: Point3) -> f64 {
        let l = self.frame.to_local(p);
        let rho = l.x.hypot(l.y);
        let (s, c) = self.semi_angle.sin_cos();
        let to_line = |rho: f64| ((rho - self.radius) * c - l.z * s).abs();
        to_line(rho).min(to_line(-rho))
    }
}

/// Sphere centred on the frame origin; `u` longitude, `v` latitude in `[-π/2, π/2]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere {
    frame: Frame,
    radius: f64,
}

impl Sphere {
    /// Creates a sphere.
    pub fn new(frame: Frame, radius: f64, precision: &Precision) -> Result<Self> {
        check_radius(radius, precision, "sphere radius")?;
        Ok(Self { frame, radius })
    }

    /// Frame.
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// Radius.
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// Evaluates at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> SurfaceDerivs {
        let (r, t) = radial(&self.frame, u);
        let (sv, cv) = v.sin_cos();
        let z = self.frame.z();
        SurfaceDerivs {
            point: self.frame.origin() + (r * cv + z * sv) * self.radius,
            du: t * (self.radius * cv),
            dv: (z * cv - r * sv) * self.radius,
        }
    }

    fn distance(&self, p: Point3) -> f64 {
        (p.distance(self.frame.origin()) - self.radius).abs()
    }
}

/// Torus about the frame `z` axis; `u` around the axis, `v` around the tube.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Torus {
    frame: Frame,
    major: f64,
    minor: f64,
}

impl Torus {
    /// Creates a ring torus (`major > minor`). Apple and lemon tori are not supported yet.
    pub fn new(frame: Frame, major: f64, minor: f64, precision: &Precision) -> Result<Self> {
        check_radius(minor, precision, "torus minor radius")?;
        if !(major.is_finite() && major > minor + precision.linear()) {
            return Err(Error::unsupported(
                "only ring tori (major > minor) are supported",
            ));
        }
        Ok(Self {
            frame,
            major,
            minor,
        })
    }

    /// Frame.
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// Major and minor radii.
    pub const fn radii(&self) -> (f64, f64) {
        (self.major, self.minor)
    }

    /// Evaluates at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> SurfaceDerivs {
        let (r, t) = radial(&self.frame, u);
        let (sv, cv) = v.sin_cos();
        let z = self.frame.z();
        let ring = self.major + self.minor * cv;
        SurfaceDerivs {
            point: self.frame.origin() + r * ring + z * (self.minor * sv),
            du: t * ring,
            dv: (z * cv - r * sv) * self.minor,
        }
    }

    fn distance(&self, p: Point3) -> f64 {
        let l = self.frame.to_local(p);
        ((l.x.hypot(l.y) - self.major).hypot(l.z) - self.minor).abs()
    }
}

/// Profile curve translated along `direction`: `S(u, v) = C(u) + v * direction`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Swept {
    /// Profile curve.
    pub profile: CurveId,
    /// Unit sweep direction.
    pub direction: Vec3,
}

impl Swept {
    /// Creates a swept surface; `direction` is normalised.
    pub fn new(profile: CurveId, direction: Vec3, precision: &Precision) -> Result<Self> {
        let direction = direction
            .normalized(precision.degenerate_length())
            .ok_or_else(|| Error::degenerate("sweep direction has zero length"))?;
        Ok(Self { profile, direction })
    }
}

/// Profile curve rotated about an axis: `S(u, v) = Rot(axis, v) C(u)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spun {
    /// Profile curve.
    pub profile: CurveId,
    /// Point on the axis.
    pub axis_origin: Point3,
    /// Unit axis direction.
    pub axis_direction: Vec3,
}

impl Spun {
    /// Creates a spun surface; the axis direction is normalised.
    pub fn new(
        profile: CurveId,
        axis_origin: Point3,
        axis_direction: Vec3,
        precision: &Precision,
    ) -> Result<Self> {
        if !axis_origin.is_finite() {
            return Err(Error::invalid_argument("spin axis origin is not finite"));
        }
        let axis_direction = axis_direction
            .normalized(precision.degenerate_length())
            .ok_or_else(|| Error::degenerate("spin axis has zero length"))?;
        Ok(Self {
            profile,
            axis_origin,
            axis_direction,
        })
    }

    /// Rotates `v` by `angle` about the unit axis (Rodrigues' formula).
    pub(crate) fn rotate(&self, v: Vec3, angle: f64) -> Vec3 {
        let k = self.axis_direction;
        let (s, c) = angle.sin_cos();
        v * c + k.cross(v) * s + k * (k.dot(v) * (1.0 - c))
    }
}

/// Every surface class of the kernel.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Surface {
    /// Plane.
    Plane(Plane),
    /// Circular cylinder.
    Cylinder(Cylinder),
    /// Circular cone.
    Cone(Cone),
    /// Sphere.
    Sphere(Sphere),
    /// Ring torus.
    Torus(Torus),
    /// NURBS surface.
    BSurface(BSurface),
    /// Linear extrusion of a curve.
    Swept(Swept),
    /// Revolution of a curve.
    Spun(Spun),
}

impl Surface {
    /// Distance from `p` to an analytic surface, or `None` for classes that need the store.
    pub(crate) fn analytic_distance(&self, p: Point3) -> Option<f64> {
        match self {
            Self::Plane(s) => Some(s.distance(p)),
            Self::Cylinder(s) => Some(s.distance(p)),
            Self::Cone(s) => Some(s.distance(p)),
            Self::Sphere(s) => Some(s.distance(p)),
            Self::Torus(s) => Some(s.distance(p)),
            Self::BSurface(_) | Self::Swept(_) | Self::Spun(_) => None,
        }
    }

    /// Ranges of analytic and NURBS surfaces; `None` for procedural ones.
    pub(crate) fn intrinsic_ranges(&self) -> Option<(ParamRange, ParamRange)> {
        let latitude = ParamRange::bounded(-FRAC_PI_2, FRAC_PI_2);
        match self {
            Self::Plane(_) => Some((ParamRange::UNBOUNDED, ParamRange::UNBOUNDED)),
            Self::Cylinder(_) | Self::Cone(_) => {
                Some((ParamRange::FULL_TURN, ParamRange::UNBOUNDED))
            }
            Self::Sphere(_) => Some((ParamRange::FULL_TURN, latitude)),
            Self::Torus(_) => Some((
                ParamRange::FULL_TURN,
                ParamRange {
                    lo: 0.0,
                    hi: TAU,
                    periodic: true,
                },
            )),
            Self::BSurface(s) => Some(s.ranges()),
            Self::Swept(_) | Self::Spun(_) => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn frame() -> Frame {
        Frame::new(
            Point3::new(0.5, -1.0, 2.0),
            Vec3::new(0.0, 1.0, 1.0),
            Vec3::X,
            &Precision::DEFAULT,
        )
        .unwrap()
    }

    fn analytic() -> Vec<Surface> {
        let p = Precision::DEFAULT;
        vec![
            Surface::Plane(Plane::new(frame())),
            Surface::Cylinder(Cylinder::new(frame(), 2.0, &p).unwrap()),
            Surface::Cone(Cone::new(frame(), 1.5, 0.4, &p).unwrap()),
            Surface::Sphere(Sphere::new(frame(), 3.0, &p).unwrap()),
            Surface::Torus(Torus::new(frame(), 4.0, 1.0, &p).unwrap()),
        ]
    }

    fn eval(s: &Surface, u: f64, v: f64) -> SurfaceDerivs {
        match s {
            Surface::Plane(s) => s.eval(u, v),
            Surface::Cylinder(s) => s.eval(u, v),
            Surface::Cone(s) => s.eval(u, v),
            Surface::Sphere(s) => s.eval(u, v),
            Surface::Torus(s) => s.eval(u, v),
            _ => unreachable!(),
        }
    }

    proptest! {
        #[test]
        fn evaluated_points_have_zero_distance(u in 0.0f64..TAU, v in -1.2f64..1.2) {
            let tol = Precision::DEFAULT.linear();
            for s in analytic() {
                let d = eval(&s, u, v);
                prop_assert!(s.analytic_distance(d.point).unwrap() <= tol, "{s:?}");
            }
        }

        #[test]
        fn partials_match_finite_differences(u in 0.0f64..TAU, v in -1.2f64..1.2) {
            let h = Precision::DEFAULT.linear().sqrt();
            for s in analytic() {
                let d = eval(&s, u, v);
                let fu = (eval(&s, u + h, v).point - eval(&s, u - h, v).point) / (2.0 * h);
                let fv = (eval(&s, u, v + h).point - eval(&s, u, v - h).point) / (2.0 * h);
                prop_assert!((fu - d.du).length() <= h * (1.0 + d.du.length()), "{s:?}");
                prop_assert!((fv - d.dv).length() <= h * (1.0 + d.dv.length()), "{s:?}");
            }
        }

        #[test]
        fn offset_points_have_offset_distance(u in 0.0f64..TAU, v in -1.2f64..1.2, off in -0.5f64..0.5) {
            let p = Precision::DEFAULT;
            for s in analytic() {
                let d = eval(&s, u, v);
                let Some(n) = d.normal(&p) else { continue };
                let q = d.point + n * off;
                prop_assert!((s.analytic_distance(q).unwrap() - off.abs()).abs() <= p.linear(), "{s:?}");
            }
        }
    }

    #[test]
    fn normals_point_outward() {
        let p = Precision::DEFAULT;
        for s in analytic().into_iter().skip(1) {
            let d = eval(&s, 0.3, 0.2);
            let n = d.normal(&p).unwrap();
            let outside = d.point + n * 0.01;
            let inside = d.point - n * 0.01;
            let axis_dist = |q: Point3| {
                let l = frame().to_local(q);
                l.x.hypot(l.y)
            };
            match s {
                Surface::Sphere(_) => {
                    assert!(outside.distance(frame().origin()) > inside.distance(frame().origin()));
                }
                _ => assert!(axis_dist(outside) > axis_dist(inside), "{s:?}"),
            }
        }
    }
}
