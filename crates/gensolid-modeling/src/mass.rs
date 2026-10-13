//! Original analytic integration of planar B-rep loops at uniform unit density.
//! Signed simplex terms are boundary integration identities, not a mesh operation.
//! Formula context: Mirtich, JGT 1(2), 1996,
//! https://people.eecs.berkeley.edu/~jfc/mirtich/massProps.html . No reference code is used.

use gensolid_core::{Error, Result};
use gensolid_geom::{Curve, Surface};
use gensolid_math::{Interval, Point3, Sign, Vec3};
use gensolid_topo::{BodyId, BodyKind, Model, Sense};

/// Uniform-unit-density properties in source model units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties {
    /// Material volume (and mass at unit density).
    pub volume: f64,
    /// Area of the material boundary, including cavity boundaries.
    pub surface_area: f64,
    /// Centre of mass in the model coordinate frame.
    pub centroid: Point3,
    /// Symmetric inertia tensor about the centroid, expressed in world axes.
    pub inertia: [[f64; 3]; 3],
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let next = self.value + x;
        // This branch chooses a compensated arithmetic expression, not topology.
        self.correction += if self.value.abs() >= x.abs() {
            (self.value - next) + x
        } else {
            (x - next) + self.value
        };
        self.value = next;
    }
    fn total(&self) -> f64 {
        self.value + self.correction
    }
}

fn vector_interval(a: Vec3) -> [Interval; 3] {
    [
        Interval::point(a.x),
        Interval::point(a.y),
        Interval::point(a.z),
    ]
}
fn difference(a: Point3, b: Point3) -> [Interval; 3] {
    [
        Interval::point(a.x) - Interval::point(b.x),
        Interval::point(a.y) - Interval::point(b.y),
        Interval::point(a.z) - Interval::point(b.z),
    ]
}
fn determinant(a: [Interval; 3], b: [Interval; 3], c: [Interval; 3]) -> Interval {
    a[0] * (b[1] * c[2] - b[2] * c[1])
        + a[1] * (b[2] * c[0] - b[0] * c[2])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// Integrate checked planar, straight-edged B-rep boundaries, including oriented holes
/// and finite cavities. Curved carriers are explicitly unimplemented.
pub fn body_ask_mass_properties(model: &Model, body: BodyId) -> Result<MassProperties> {
    crate::body_check_complete(model, body)?;
    if model.body(body)?.kind != BodyKind::Solid {
        return Err(Error::unsupported("mass properties require a solid body"));
    }
    let vertices = model.body_vertices(body)?;
    let first = vertices
        .first()
        .ok_or_else(|| Error::degenerate("body has no vertices"))?;
    let reference = model.vertex(*first)?.point;
    let mut volume = Sum::default();
    let mut area = Sum::default();
    let mut moment: [Sum; 3] = std::array::from_fn(|_| Sum::default());
    let mut second: [[Sum; 3]; 3] =
        std::array::from_fn(|_| std::array::from_fn(|_| Sum::default()));
    let mut volume_interval = Interval::point(0.0);
    let mut face_count = 0;
    for face_id in model.body_faces(body)? {
        let face = model.face(face_id)?;
        let material = |shell: Option<gensolid_topo::ShellId>| -> Result<bool> {
            match shell {
                Some(shell) => Ok(model.region(model.shell(shell)?.region)?.solid),
                None => Ok(false),
            }
        };
        let out = material(face.back)?;
        if out == material(face.front)? {
            continue;
        }
        face_count += 1;
        let surface = face
            .surface
            .ok_or_else(|| Error::invalid_argument("face has no surface"))?;
        let Surface::Plane(plane) = model.geom().surface(surface)? else {
            return Err(Error::unimplemented(
                "curved B-rep mass integration is not implemented",
            ));
        };
        let sign = if (face.sense == Sense::Positive) == out {
            1.0
        } else {
            -1.0
        };
        let normal = plane.frame().z() * sign;
        let normal = normal
            .normalized(0.0)
            .ok_or_else(|| Error::degenerate("face normal is singular"))?;
        let mut face_area = Sum::default();
        let mut face_area_interval = Interval::point(0.0);
        for loop_id in &face.loops {
            let mut points = Vec::new();
            for fin in &model.loop_(*loop_id)?.fins {
                let fin_value = model.fin(*fin)?;
                let edge = model.edge(fin_value.edge)?;
                let curve = edge
                    .curve
                    .ok_or_else(|| Error::invalid_argument("boundary edge has no curve"))?;
                if !matches!(model.geom().curve(curve)?, Curve::Line(_)) {
                    return Err(Error::unimplemented(
                        "curved-edge mass integration is not implemented",
                    ));
                }
                let (start, _) = model.fin_vertices(*fin)?;
                let vertex = start.ok_or_else(|| {
                    Error::unimplemented("mass integration requires bounded edges")
                })?;
                points.push(model.vertex(vertex)?.point);
            }
            if points.len() < 3 {
                return Err(Error::degenerate("mass integration needs polygon loops"));
            }
            if !out {
                points.reverse();
            }
            let anchor = points[0];
            for pair in points[1..].windows(2) {
                let a = anchor - reference;
                let b = pair[0] - reference;
                let c = pair[1] - reference;
                let det = a.dot(b.cross(c));
                let v = det / 6.0;
                volume.add(v);
                let sixth: f64 = 1.0 / 6.0;
                let sixth = Interval::new(sixth.next_down(), sixth.next_up())
                    .ok_or_else(|| Error::internal("invalid exact integration coefficient"))?;
                volume_interval = volume_interval
                    + determinant(
                        difference(anchor, reference),
                        difference(pair[0], reference),
                        difference(pair[1], reference),
                    ) * sixth;
                let vectors = [[a.x, a.y, a.z], [b.x, b.y, b.z], [c.x, c.y, c.z]];
                let sums = [a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z];
                for i in 0..3 {
                    moment[i].add(v * sums[i] / 4.0);
                    for j in 0..3 {
                        let diagonal = vectors
                            .iter()
                            .map(|vector| vector[i] * vector[j])
                            .sum::<f64>();
                        second[i][j].add(v * (sums[i] * sums[j] + diagonal) / 20.0);
                    }
                }
                face_area.add((pair[0] - anchor).cross(pair[1] - anchor).dot(normal) / 2.0);
                face_area_interval = face_area_interval
                    + determinant(
                        difference(pair[0], anchor),
                        difference(pair[1], anchor),
                        vector_interval(normal),
                    ) * Interval::point(0.5);
            }
        }
        if !face_area.total().is_finite() || face_area_interval.sign() != Some(Sign::Positive) {
            return Err(Error::degenerate(
                "positive oriented boundary area cannot be certified",
            ));
        }
        area.add(face_area.total());
    }
    if face_count == 0 || volume_interval.sign() != Some(Sign::Positive) {
        return Err(Error::degenerate(
            "positive material volume cannot be certified",
        ));
    }
    let volume = volume.total();
    let offset = [
        moment[0].total() / volume,
        moment[1].total() / volume,
        moment[2].total() / volume,
    ];
    let centroid = reference + Vec3::new(offset[0], offset[1], offset[2]);
    let centred: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| second[i][j].total() - volume * offset[i] * offset[j])
    });
    let inertia = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i == j {
                centred[(i + 1) % 3][(i + 1) % 3] + centred[(i + 2) % 3][(i + 2) % 3]
            } else {
                -centred[i][j]
            }
        })
    });
    let result = MassProperties {
        volume,
        surface_area: area.total(),
        centroid,
        inertia,
    };
    if !result.volume.is_finite()
        || !result.surface_area.is_finite()
        || !result.centroid.is_finite()
        || result.inertia.iter().flatten().any(|x| !x.is_finite())
    {
        return Err(Error::degenerate("mass properties are not finite"));
    }
    Ok(result)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use gensolid_core::Precision;
    use gensolid_math::Frame;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn blocks_retain_closed_form_mass_after_random_rigid_motion(
            x in 0.01f64..20.0,y in 0.01f64..20.0,z in 0.01f64..20.0,
            ox in -30.0f64..30.0,oy in -30.0f64..30.0,oz in -30.0f64..30.0,
            ax in -1.0f64..1.0,ay in -1.0f64..1.0,az in 0.1f64..1.0,
        ) {
            let p=Precision::DEFAULT;
            let f=Frame::from_axis(Point3::new(ox,oy,oz),Vec3::new(ax,ay,az),&p).unwrap();
            let mut model=Model::new(p);
            let body=crate::body_create_solid_block(&mut model,x,y,z,&f).unwrap();
            let mass=body_ask_mass_properties(&model,body).unwrap();
            let expected=x*y*z;
            let area=2.0*(x*y+x*z+y*z);
            prop_assert!((mass.volume-expected).abs()<=p.linear()*area);
            prop_assert!((mass.surface_area-area).abs()<=p.linear()*(x+y+z)*model.body_edges(body).unwrap().len() as f64);
            prop_assert!((mass.centroid-f.origin()).length()<=p.linear());
            let axes=[f.x(),f.y(),f.z()];
            let principal=[expected*(y*y+z*z)/12.0,expected*(x*x+z*z)/12.0,expected*(x*x+y*y)/12.0];
            for i in 0..3 {for j in 0..3 {
                let component=|v:Vec3,index:usize|match index {0=>v.x,1=>v.y,_=>v.z};
                let expected=axes.iter().zip(principal).map(|(axis,moment)|component(*axis,i)*component(*axis,j)*moment).sum::<f64>();
                prop_assert!((mass.inertia[i][j]-expected).abs()<=p.linear()*mass.volume*(x+y+z));
            }}
        }
    }

    #[test]
    fn thin_rod_keeps_small_axial_inertia() {
        let p = Precision::DEFAULT;
        let mut model = Model::new(p);
        let width = p.linear() * 10.0;
        let body =
            crate::body_create_solid_block(&mut model, 100.0, width, width, &Frame::WORLD).unwrap();
        let mass = body_ask_mass_properties(&model, body).unwrap();
        let expected = mass.volume * width * width / 6.0;
        assert!((mass.inertia[0][0] - expected).abs() <= expected * p.angular());
    }

    #[test]
    fn inverted_face_sense_is_not_positive_material_area() {
        let p = Precision::DEFAULT;
        let mut model = Model::new(p);
        let body =
            crate::body_create_solid_block(&mut model, 2.0, 3.0, 4.0, &Frame::WORLD).unwrap();
        let face = model.body_faces(body).unwrap()[0];
        let reversed = model.face(face).unwrap().sense.reversed();
        model.raw().face_mut(face).unwrap().sense = reversed;
        assert!(body_ask_mass_properties(&model, body).is_err());
    }

    #[test]
    fn cavities_and_face_holes_have_material_mass_without_meshes() {
        let p = Precision::DEFAULT;
        let cavity = gensolid_step::read_solids(
            include_str!("../../gensolid-step/tests/fixtures/hollow_block.step"),
            p,
        )
        .unwrap();
        let body = cavity.bodies().next().unwrap();
        let mass = body_ask_mass_properties(&cavity, body).unwrap();
        assert!((mass.volume - 168.0).abs() <= p.linear() * mass.surface_area);
        let holes = gensolid_step::read_solids(
            include_str!("../../gensolid-step/tests/fixtures/block_with_hole.step"),
            p,
        )
        .unwrap();
        let body = holes.bodies().next().unwrap();
        let mass = body_ask_mass_properties(&holes, body).unwrap();
        assert!((mass.volume - 20.0).abs() <= p.linear() * mass.surface_area);
    }
}
