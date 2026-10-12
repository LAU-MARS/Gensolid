//! Geometry storage and evaluation.

use gensolid_core::{Arena, Error, Id, Precision, Result};
use gensolid_math::{Point3, Vec3};

use crate::curve::{Curve, CurveDerivs};
use crate::param::ParamRange;
use crate::project;
use crate::surface::{Surface, SurfaceDerivs};

/// Handle to a curve in a [`GeomStore`].
pub type CurveId = Id<Curve>;
/// Handle to a surface in a [`GeomStore`].
pub type SurfaceId = Id<Surface>;

/// Closest point on a curve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurveProjection {
    /// Curve parameter of the closest point.
    pub param: f64,
    /// The closest point.
    pub point: Point3,
    /// Distance from the query point.
    pub distance: f64,
}

/// Owns all curves and surfaces of a model.
///
/// Procedural surfaces may only reference geometry that already exists, so the reference
/// graph is acyclic by construction and evaluation always terminates.
#[derive(Default)]
pub struct GeomStore {
    curves: Arena<Curve>,
    surfaces: Arena<Surface>,
}

impl GeomStore {
    /// Creates an empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a curve.
    pub fn add_curve(&mut self, curve: Curve) -> Result<CurveId> {
        self.curves.insert(curve)
    }

    /// Adds a surface after checking that every referenced curve exists.
    pub fn add_surface(&mut self, surface: Surface) -> Result<SurfaceId> {
        match &surface {
            Surface::Swept(s) => {
                self.curve(s.profile)?;
            }
            Surface::Spun(s) => {
                self.curve(s.profile)?;
            }
            _ => {}
        }
        self.surfaces.insert(surface)
    }

    /// Looks up a curve.
    pub fn curve(&self, id: CurveId) -> Result<&Curve> {
        self.curves
            .get(id)
            .ok_or_else(|| Error::not_found(format!("curve {id:?}")))
    }

    /// Looks up a surface.
    pub fn surface(&self, id: SurfaceId) -> Result<&Surface> {
        self.surfaces
            .get(id)
            .ok_or_else(|| Error::not_found(format!("surface {id:?}")))
    }

    /// Number of curves.
    pub fn curve_count(&self) -> usize {
        self.curves.len()
    }

    /// Number of surfaces.
    pub fn surface_count(&self) -> usize {
        self.surfaces.len()
    }

    /// Parameter range of a curve.
    pub fn curve_range(&self, id: CurveId) -> Result<ParamRange> {
        Ok(self.curve(id)?.range())
    }

    /// Parameter ranges `(u, v)` of a surface.
    pub fn surface_ranges(&self, id: SurfaceId) -> Result<(ParamRange, ParamRange)> {
        let surface = self.surface(id)?;
        if let Some(ranges) = surface.intrinsic_ranges() {
            return Ok(ranges);
        }
        match surface {
            Surface::Swept(s) => Ok((self.curve_range(s.profile)?, ParamRange::UNBOUNDED)),
            Surface::Spun(s) => Ok((self.curve_range(s.profile)?, ParamRange::FULL_TURN)),
            _ => Err(Error::internal("surface without ranges")),
        }
    }

    /// Point and derivative of a curve at `t`.
    pub fn eval_curve(&self, id: CurveId, t: f64) -> Result<CurveDerivs> {
        Ok(self.curve(id)?.eval(t))
    }

    /// Point and partial derivatives of a surface at `(u, v)`.
    pub fn eval_surface(&self, id: SurfaceId, u: f64, v: f64) -> Result<SurfaceDerivs> {
        Ok(match self.surface(id)? {
            Surface::Plane(s) => s.eval(u, v),
            Surface::Cylinder(s) => s.eval(u, v),
            Surface::Cone(s) => s.eval(u, v),
            Surface::Sphere(s) => s.eval(u, v),
            Surface::Torus(s) => s.eval(u, v),
            Surface::BSurface(s) => s.eval(u, v),
            Surface::Swept(s) => {
                let c = self.eval_curve(s.profile, u)?;
                SurfaceDerivs {
                    point: c.point + s.direction * v,
                    du: c.d1,
                    dv: s.direction,
                }
            }
            Surface::Spun(s) => {
                let c = self.eval_curve(s.profile, u)?;
                let rel = s.rotate(c.point - s.axis_origin, v);
                SurfaceDerivs {
                    point: s.axis_origin + rel,
                    du: s.rotate(c.d1, v),
                    dv: s.axis_direction.cross(rel),
                }
            }
        })
    }

    /// Closest point on a curve to `p`.
    ///
    /// Exact for lines and circles. Other curves use sampling plus Gauss-Newton refinement,
    /// which finds the global minimum when the sampling resolves the curve's features.
    pub fn project_to_curve(
        &self,
        id: CurveId,
        p: Point3,
        precision: &Precision,
    ) -> Result<CurveProjection> {
        let curve = self.curve(id)?;
        let exact = |param: f64| {
            let point = curve.eval(param).point;
            CurveProjection {
                param,
                point,
                distance: point.distance(p),
            }
        };
        match curve {
            Curve::Line(c) => Ok(exact(c.closest_param(p))),
            Curve::Circle(c) => Ok(exact(c.closest_param(p))),
            _ => project::closest_on_curve(|t| curve.eval(t), curve.range(), p, precision),
        }
    }

    /// Distance from `p` to a surface (ignoring any face boundary).
    ///
    /// Closed form for analytic surfaces; numerical for NURBS and procedural surfaces,
    /// with the same caveat as [`GeomStore::project_to_curve`].
    pub fn distance_to_surface(
        &self,
        id: SurfaceId,
        p: Point3,
        precision: &Precision,
    ) -> Result<f64> {
        let surface = self.surface(id)?;
        if let Some(d) = surface.analytic_distance(p) {
            return Ok(d);
        }
        match surface {
            Surface::BSurface(s) => {
                let (ur, vr) = s.ranges();
                project::distance_to_patch(|u, v| s.eval(u, v), ur, vr, p, precision)
            }
            Surface::Swept(s) => {
                // For fixed u the optimal v is exact, so only the profile projected onto
                // the plane perpendicular to the sweep remains to be searched.
                let d = s.direction;
                let flatten = |v: Vec3| v - d * v.dot(d);
                let target = Point3::from_vec(flatten(p.to_vec()));
                let range = self.curve_range(s.profile)?;
                let curve = self.curve(s.profile)?;
                let flat = |t: f64| {
                    let c = curve.eval(t);
                    CurveDerivs {
                        point: Point3::from_vec(flatten(c.point.to_vec())),
                        d1: flatten(c.d1),
                    }
                };
                Ok(project::closest_on_curve(flat, range, target, precision)?.distance)
            }
            Surface::Spun(s) => {
                // Rotation sweeps each profile point over a full circle, so the distance
                // is the planar distance in (radius, height) coordinates about the axis.
                let (o, k) = (s.axis_origin, s.axis_direction);
                let half_plane = |q: Point3| {
                    let rel = q - o;
                    let h = rel.dot(k);
                    (rel - k * h, h)
                };
                let (radial, height) = half_plane(p);
                let target = Point3::new(radial.length(), height, 0.0);
                let range = self.curve_range(s.profile)?;
                let curve = self.curve(s.profile)?;
                let planar = |t: f64| {
                    let c = curve.eval(t);
                    let (r, h) = half_plane(c.point);
                    let rho = r.length();
                    let d1_axial = c.d1.dot(k);
                    let d1_radial = c.d1 - k * d1_axial;
                    let d_rho = if rho > precision.degenerate_length() {
                        r.dot(d1_radial) / rho
                    } else {
                        d1_radial.length()
                    };
                    CurveDerivs {
                        point: Point3::new(rho, h, 0.0),
                        d1: Vec3::new(d_rho, d1_axial, 0.0),
                    }
                };
                Ok(project::closest_on_curve(planar, range, target, precision)?.distance)
            }
            _ => Err(Error::internal("surface without distance evaluation")),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use gensolid_math::Frame;

    use super::*;
    use crate::curve::{Circle, Line};
    use crate::nurbs::BCurve;
    use crate::surface::{Cylinder, Sphere, Spun, Swept};

    #[test]
    fn swept_circle_matches_cylinder() {
        let p = Precision::DEFAULT;
        let mut store = GeomStore::new();
        let circle = Circle::new(Frame::WORLD, 2.0, &p).unwrap();
        let profile = store.add_curve(Curve::Circle(circle)).unwrap();
        let swept = store
            .add_surface(Surface::Swept(Swept::new(profile, Vec3::Z, &p).unwrap()))
            .unwrap();
        let cyl = store
            .add_surface(Surface::Cylinder(
                Cylinder::new(Frame::WORLD, 2.0, &p).unwrap(),
            ))
            .unwrap();
        for (u, v) in [(0.0, 0.0), (1.0, 2.5), (4.0, -3.0)] {
            let a = store.eval_surface(swept, u, v).unwrap();
            let b = store.eval_surface(cyl, u, v).unwrap();
            assert!(a.point.distance(b.point) <= p.linear());
            assert!((a.du - b.du).length() <= p.linear());
        }
        let q = Point3::new(3.0, 1.0, 7.0);
        let expected = (10.0f64.sqrt() - 2.0).abs();
        assert!((store.distance_to_surface(swept, q, &p).unwrap() - expected).abs() <= p.linear());
    }

    #[test]
    fn spun_semicircle_matches_sphere() {
        let p = Precision::DEFAULT;
        let mut store = GeomStore::new();
        // Meridian of a unit sphere: a half circle in the x-z plane from the south pole up.
        let frame = Frame::new(Point3::ORIGIN, -Vec3::Y, Vec3::X, &p).unwrap();
        let circle = store
            .add_curve(Curve::Circle(Circle::new(frame, 1.0, &p).unwrap()))
            .unwrap();
        let spun = store
            .add_surface(Surface::Spun(
                Spun::new(circle, Point3::ORIGIN, Vec3::Z, &p).unwrap(),
            ))
            .unwrap();
        let sphere = store
            .add_surface(Surface::Sphere(Sphere::new(Frame::WORLD, 1.0, &p).unwrap()))
            .unwrap();
        let d = store.eval_surface(spun, FRAC_PI_2, 1.0).unwrap();
        assert!(store.distance_to_surface(sphere, d.point, &p).unwrap() <= p.linear());
        for q in [
            Point3::new(2.0, 0.5, 0.3),
            Point3::new(0.1, -0.2, 0.4),
            Point3::new(-1.0, 1.0, -1.0),
        ] {
            let a = store.distance_to_surface(spun, q, &p).unwrap();
            let b = store.distance_to_surface(sphere, q, &p).unwrap();
            assert!((a - b).abs() <= p.linear(), "{q:?}: {a} vs {b}");
        }
    }

    #[test]
    fn projects_onto_bcurve() {
        let p = Precision::DEFAULT;
        let mut store = GeomStore::new();
        let curve = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ],
            None,
        )
        .unwrap();
        let id = store.add_curve(Curve::BCurve(curve)).unwrap();
        let on = store.eval_curve(id, 0.3).unwrap().point;
        let proj = store.project_to_curve(id, on, &p).unwrap();
        assert!(proj.distance <= p.linear());
        assert!((proj.param - 0.3).abs() <= p.linear().sqrt());
    }

    #[test]
    fn rejects_dangling_profile() {
        let p = Precision::DEFAULT;
        let mut empty = GeomStore::new();
        let mut other = GeomStore::new();
        let line = Line::new(Point3::ORIGIN, Vec3::X, &p).unwrap();
        let foreign = other.add_curve(Curve::Line(line)).unwrap();
        let swept = Surface::Swept(Swept::new(foreign, Vec3::Z, &p).unwrap());
        assert!(empty.add_surface(swept).is_err());
    }
}
