//! Planar B-rep output tessellation.

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{Error, Precision, Result, Tolerance, effective_tolerance};
use gensolid_geom::{Curve, Line, Plane, Surface};
use gensolid_math::{Interval, Point3, Sign, Vec3};
use gensolid_topo::{BodyId, BodyKind, Model, Sense};

use crate::Tessellation;

#[derive(Clone, Copy)]
struct Projected {
    x: Interval,
    y: Interval,
}

fn coordinates(point: Point3) -> [Interval; 3] {
    [
        Interval::point(point.x),
        Interval::point(point.y),
        Interval::point(point.z),
    ]
}
fn components(vector: Vec3) -> [Interval; 3] {
    [
        Interval::point(vector.x),
        Interval::point(vector.y),
        Interval::point(vector.z),
    ]
}
fn dot(a: [Interval; 3], b: [Interval; 3]) -> Interval {
    (0..3).fold(Interval::point(0.0), |sum, i| sum + a[i] * b[i])
}
fn reciprocal_positive(value: Interval) -> Result<Interval> {
    if value.sign() != Some(Sign::Positive) {
        return Err(Error::degenerate(
            "projection norm has no certified positive bound",
        ));
    }
    Interval::new((1.0 / value.hi()).next_down(), (1.0 / value.lo()).next_up())
        .ok_or_else(|| Error::degenerate("projection reciprocal is unbounded"))
}
fn projected(point: Point3, plane: &Plane) -> Result<Projected> {
    let (point, origin, normal) = (
        coordinates(point),
        coordinates(plane.frame().origin()),
        components(plane.frame().z()),
    );
    let delta = std::array::from_fn(|i| point[i] - origin[i]);
    let parameter = dot(delta, normal) * reciprocal_positive(dot(normal, normal))?;
    let shadow = std::array::from_fn(|i| delta[i] - normal[i] * parameter);
    Ok(Projected {
        x: dot(shadow, components(plane.frame().x())),
        y: dot(shadow, components(plane.frame().y())),
    })
}
fn orientation(a: Projected, b: Projected, c: Projected) -> Result<Sign> {
    let determinant = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    // The predicates concern the exact plane projections, not their rounded UV
    // approximations. An undecidable orientation cannot certify an ear or a
    // removable collinear vertex, so it is rejected instead of being snapped.
    determinant
        .sign()
        .ok_or_else(|| Error::unimplemented("polygon projection orientation cannot be certified"))
}
fn plane_basis(plane: &Plane) -> Result<()> {
    let (a, b, c) = (
        components(plane.frame().x()),
        components(plane.frame().y()),
        components(plane.frame().z()),
    );
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    if dot(cross, c).sign() != Some(Sign::Positive) {
        return Err(Error::unimplemented(
            "plane projection handedness cannot be certified",
        ));
    }
    Ok(())
}
fn line_endpoint_bound(line: &Line, point: Point3) -> Result<f64> {
    let (point, origin, direction) = (
        coordinates(point),
        coordinates(line.origin()),
        components(line.direction()),
    );
    let delta = std::array::from_fn(|i| point[i] - origin[i]);
    let parameter = dot(delta, direction);
    // Enclose both the kernel's unit-direction closest_param expression and the
    // exact projection for the stored direction's slight norm roundoff.
    let exact = parameter * reciprocal_positive(dot(direction, direction))?;
    let parameter = parameter.hull(exact);
    let residual: [Interval; 3] = std::array::from_fn(|i| delta[i] - direction[i] * parameter);
    let squared = residual
        .into_iter()
        .fold(Interval::point(0.0), |sum, value| sum + value.sqr());
    let bound = squared.hi().sqrt().next_up();
    if !bound.is_finite() {
        return Err(Error::degenerate(
            "trim carrier endpoint bound is not finite",
        ));
    }
    Ok(bound)
}
fn opposite(a: Sign, b: Sign) -> bool {
    matches!(
        (a, b),
        (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
    )
}
fn contains_range(a: Interval, b: Interval, x: Interval, p: &Precision) -> Result<bool> {
    let product = (x - a) * (x - b);
    match product.sign() {
        Some(Sign::Negative | Sign::Zero) => Ok(true),
        Some(Sign::Positive) => Ok(false),
        None => {
            // Only an undecidable interval uses the explicit precision margin.
            let width = b - a;
            let margin = p.solver() * (width.lo().abs().max(width.hi().abs()) + p.solver());
            if product.lo() >= -margin && product.hi() <= margin {
                Ok(true)
            } else {
                Err(Error::unsupported("segment bound is unresolved"))
            }
        }
    }
}
fn simple(points: &[Projected], p: &Precision) -> Result<()> {
    let n = points.len();
    for i in 0..n {
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            let (a, b, c, d) = (
                points[i],
                points[(i + 1) % n],
                points[j],
                points[(j + 1) % n],
            );
            let signs = [
                orientation(a, b, c)?,
                orientation(a, b, d)?,
                orientation(c, d, a)?,
                orientation(c, d, b)?,
            ];
            if opposite(signs[0], signs[1]) && opposite(signs[2], signs[3]) {
                return Err(Error::invalid_argument("self-intersecting polygon loop"));
            }
            for (s, x, y, z) in [
                (signs[0], a, b, c),
                (signs[1], a, b, d),
                (signs[2], c, d, a),
                (signs[3], c, d, b),
            ] {
                if s == Sign::Zero
                    && contains_range(x.x, y.x, z.x, p)?
                    && contains_range(x.y, y.y, z.y, p)?
                {
                    return Err(Error::invalid_argument(
                        "polygon has overlapping nonadjacent edges",
                    ));
                }
            }
        }
    }
    Ok(())
}
fn triangles(points: &[Projected], expected: Sign, p: &Precision) -> Result<Vec<[usize; 3]>> {
    simple(points, p)?;
    let mut active: Vec<_> = (0..points.len()).collect();
    let mut result = Vec::new();
    while active.len() > 3 {
        let mut found = false;
        for i in 0..active.len() {
            let prev = active[(i + active.len() - 1) % active.len()];
            let here = active[i];
            let next = active[(i + 1) % active.len()];
            let turn = orientation(points[prev], points[here], points[next])?;
            if turn == Sign::Zero {
                active.remove(i);
                found = true;
                break;
            }
            if turn != expected {
                continue;
            }
            let mut blocked = false;
            for &other in &active {
                if [prev, here, next].contains(&other) {
                    continue;
                }
                let signs = [
                    orientation(points[prev], points[here], points[other])?,
                    orientation(points[here], points[next], points[other])?,
                    orientation(points[next], points[prev], points[other])?,
                ];
                if signs
                    .iter()
                    .all(|sign| *sign == expected || *sign == Sign::Zero)
                {
                    blocked = true;
                    break;
                }
            }
            if !blocked {
                result.push([prev, here, next]);
                active.remove(i);
                found = true;
                break;
            }
        }
        if !found {
            return Err(Error::unsupported(
                "polygon has no certified ear at its required orientation",
            ));
        }
    }
    if active.len() != 3
        || orientation(points[active[0]], points[active[1]], points[active[2]])? != expected
    {
        return Err(Error::degenerate("polygon residual triangle is degenerate"));
    }
    result.push([active[0], active[1], active[2]]);
    Ok(result)
}

/// Tessellate simple planar straight-edged material faces. Finite cavity boundaries
/// are supported; face holes and curved trims remain explicit unimplemented cases.
/// The bound includes the plane shadow displacement, each bounded line carrier's
/// endpoint projection error plus its edge tube, and every vertex tolerance ball.
/// The exact projected polygon must have certified topology; ambiguous projections,
/// including unresolved collinear ears, are rejected without altering the B-rep.
pub fn body_tessellate_planar(
    model: &Model,
    body: BodyId,
    chord: Tolerance,
) -> Result<Tessellation> {
    let chord = Tolerance::new(chord.value(), model.precision())?.value();
    if model.body(body)?.kind != BodyKind::Solid {
        return Err(Error::unsupported("body tessellation requires a solid"));
    }
    if !check_body(model, body, &CheckOptions::default())?.is_valid() {
        return Err(Error::invalid_argument(
            "body does not pass complete validity checks",
        ));
    }
    let mut mesh = Tessellation {
        vertices: Vec::new(),
        triangles: Vec::new(),
        interpolation_bound: 0.0,
        vertex_roundoff_bound: 0.0,
    };
    for id in model.body_faces(body)? {
        let face = model.face(id)?;
        let solid = |shell: Option<gensolid_topo::ShellId>| -> Result<bool> {
            match shell {
                Some(shell) => Ok(model.region(model.shell(shell)?.region)?.solid),
                None => Ok(false),
            }
        };
        let (front, back) = (solid(face.front)?, solid(face.back)?);
        if front == back {
            continue;
        }
        if face.loops.len() != 1 {
            return Err(Error::unimplemented(
                "planar mesh faces with holes are not implemented",
            ));
        }
        let surface = face
            .surface
            .ok_or_else(|| Error::invalid_argument("face has no carrier"))?;
        let Surface::Plane(plane) = model.geom().surface(surface)? else {
            return Err(Error::unimplemented(
                "trimmed curved-body tessellation is not implemented",
            ));
        };
        plane_basis(plane)?;
        let mut points = Vec::new();
        let mut local_tolerance = model.precision().default_local();
        for fin in &model.loop_(face.loops[0])?.fins {
            let edge = model.edge(model.fin(*fin)?.edge)?;
            let curve = edge
                .curve
                .ok_or_else(|| Error::invalid_argument("edge has no carrier"))?;
            let Curve::Line(line) = model.geom().curve(curve)? else {
                return Err(Error::unimplemented(
                    "curved planar face trims are not implemented",
                ));
            };
            let edge_tolerance = effective_tolerance(edge.tolerance, model.precision());
            let (start, end) = model.fin_vertices(*fin)?;
            for endpoint in [start, end] {
                let endpoint = model.vertex(endpoint.ok_or_else(|| {
                    Error::unimplemented("mesh trim carrier needs two endpoints")
                })?)?;
                let vertex_tolerance = effective_tolerance(endpoint.tolerance, model.precision());
                let carrier =
                    (line_endpoint_bound(line, endpoint.point)? + edge_tolerance).next_up();
                local_tolerance = local_tolerance.max(carrier).max(vertex_tolerance);
            }
            let vertex =
                model
                    .vertex(start.ok_or_else(|| {
                        Error::unimplemented("mesh polygon needs bounded edges")
                    })?)?;
            if let Some(tol) = vertex.tolerance {
                local_tolerance = local_tolerance.max(tol.value());
            }
            points.push(vertex.point);
        }
        if points.len() < 3 {
            return Err(Error::degenerate("mesh face has fewer than three vertices"));
        }
        if front {
            points.reverse();
        }
        let expected = if (face.sense == Sense::Positive) == back {
            Sign::Positive
        } else {
            Sign::Negative
        };
        let projected = points
            .iter()
            .map(|&point| projected(point, plane))
            .collect::<Result<Vec<_>>>()?;
        let tris = triangles(&projected, expected, model.precision())?;
        let offset = mesh.vertices.len();
        for point in &points {
            let origin = plane.frame().origin();
            let n = plane.frame().z();
            let components = [
                (point.x, origin.x, n.x),
                (point.y, origin.y, n.y),
                (point.z, origin.z, n.z),
            ];
            let signed = components
                .into_iter()
                .fold(Interval::point(0.0), |acc, (x, o, n)| {
                    acc + (Interval::point(x) - Interval::point(o)) * Interval::point(n)
                });
            let norm_squared = Interval::point(n.x).sqr()
                + Interval::point(n.y).sqr()
                + Interval::point(n.z).sqr();
            let norm_lower = norm_squared.lo().sqrt().next_down();
            if Interval::point(norm_lower).sign() != Some(Sign::Positive) {
                return Err(Error::degenerate("plane normal has no positive norm bound"));
            }
            let displacement = (signed.lo().abs().max(signed.hi().abs()) / norm_lower).next_up();
            mesh.interpolation_bound = mesh
                .interpolation_bound
                .max((displacement + local_tolerance).next_up());
        }
        mesh.vertices.extend(points);
        mesh.triangles.extend(
            tris.into_iter()
                .map(|triangle| triangle.map(|i| i + offset)),
        );
    }
    if mesh.triangles.is_empty() {
        return Err(Error::degenerate("body has no tessellated boundary"));
    }
    if mesh.deviation_bound() > chord {
        return Err(Error::unsupported(
            "local B-rep tolerance exceeds mesh chord budget",
        ));
    }
    Ok(mesh)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use gensolid_core::Precision;
    use gensolid_math::{Frame, Point3, Vec3};
    use proptest::prelude::*;

    #[test]
    fn displaced_shared_edge_tubes_are_not_bounded_by_the_largest_single_tolerance() {
        let p = Precision::DEFAULT;
        let mut model = Model::new(p);
        let body = gensolid_topo::make_block(&mut model, 2.0, 2.0, 2.0, &Frame::WORLD).unwrap();
        let tolerance = Tolerance::new(p.linear() * 100.0, &p).unwrap();
        let edge_id = model
            .body_edges(body)
            .unwrap()
            .into_iter()
            .find(|&id| {
                let edge = model.edge(id).unwrap();
                let a = model.vertex(edge.start.unwrap()).unwrap().point;
                let b = model.vertex(edge.end.unwrap()).unwrap().point;
                a.x.is_sign_negative()
                    && a.y.is_sign_negative()
                    && b.x.is_sign_negative()
                    && b.y.is_sign_negative()
            })
            .unwrap();
        let edge = model.edge(edge_id).unwrap().clone();
        let Curve::Line(line) = model.geom().curve(edge.curve.unwrap()).unwrap() else {
            panic!("block edge must be a line")
        };
        let displacement = Vec3::new(-0.9 * tolerance.value(), -0.9 * tolerance.value(), 0.0);
        let displaced = Line::new(line.origin() + displacement, line.direction(), &p).unwrap();
        let midpoint = model.vertex(edge.start.unwrap()).unwrap().point
            + (model.vertex(edge.end.unwrap()).unwrap().point
                - model.vertex(edge.start.unwrap()).unwrap().point)
                * 0.5
            + displacement;
        let mut raw = model.raw();
        let curve = raw.geom_mut().add_curve(Curve::Line(displaced)).unwrap();
        let mutable = raw.edge_mut(edge_id).unwrap();
        mutable.curve = Some(curve);
        mutable.tolerance = Some(tolerance);
        for vertex in [edge.start.unwrap(), edge.end.unwrap()] {
            raw.vertex_mut(vertex).unwrap().tolerance = Some(tolerance);
        }
        assert!(
            check_body(&model, body, &CheckOptions::default())
                .unwrap()
                .is_valid()
        );
        let narrow = Tolerance::new(tolerance.value() * 1.1, &p).unwrap();
        assert!(
            matches!(body_tessellate_planar(&model,body,narrow),Err(error) if error.code()==gensolid_core::ErrorCode::Unsupported)
        );
        let wide = Tolerance::new(tolerance.value() * 4.0, &p).unwrap();
        let mesh = body_tessellate_planar(&model, body, wide).unwrap();
        // All output triangles remain on the original cube. This carrier point
        // is outside both adjacent supporting planes, so its exact nearest mesh
        // point is (-1,-1,0), farther than the old max(edge_tol,vertex_tol).
        let measured = midpoint.distance(Point3::new(-1.0, -1.0, 0.0));
        assert!(measured > narrow.value());
        assert!(measured + tolerance.value() <= mesh.deviation_bound());
        assert!(mesh.deviation_bound() <= wide.value());
    }

    #[test]
    fn collinear_projected_ears_need_an_actual_orientation_certificate() {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(
            Point3::new(100.0, -75.0, 25.0),
            Vec3::new(1.0, 2.0, 3.0),
            &p,
        )
        .unwrap();
        let plane = Plane::new(frame);
        let projected = [
            Point3::new(-1.0, 0.0, 0.0),
            Point3::ORIGIN,
            Point3::new(1.0, 0.0, 0.0),
        ]
        .map(|point| projected(point, &plane).unwrap());
        assert!(
            matches!(orientation(projected[0],projected[1],projected[2]),Err(error) if error.code()==gensolid_core::ErrorCode::Unimplemented)
        );
    }
    proptest! {
        #[test]
        fn random_block_mesh_is_closed_oriented_and_within_chord_budget(
            x in 0.01f64..20.0,y in 0.01f64..20.0,z in 0.01f64..20.0,
            ox in -30.0f64..30.0,oy in -30.0f64..30.0,oz in -30.0f64..30.0,
            ax in -1.0f64..1.0,ay in -1.0f64..1.0,az in 0.1f64..1.0,
        ) {
            let p=Precision::DEFAULT;let f=Frame::from_axis(Point3::new(ox,oy,oz),Vec3::new(ax,ay,az),&p).unwrap();
            let mut model=Model::new(p);let body=gensolid_topo::make_block(&mut model,x,y,z,&f).unwrap();
            let chord=Tolerance::new(p.linear()*64.0,&p).unwrap();let mesh=body_tessellate_planar(&model,body,chord).unwrap();
            prop_assert_eq!(mesh.triangles.len(),12);
            prop_assert!(mesh.deviation_bound()<=chord.value());
            for triangle in mesh.triangles {let a=mesh.vertices[triangle[0]];let b=mesh.vertices[triangle[1]];let c=mesh.vertices[triangle[2]];
                prop_assert!((b-a).cross(c-a).dot(a-f.origin())>0.0);
            }
        }
    }
}
