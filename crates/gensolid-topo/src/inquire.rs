//! Inquiries over an existing body: bounding box and solid volume.
//!
//! Volume uses the divergence theorem on the face loops. Each loop is fanned into
//! triangles; a closed boundary oriented out of the solid contributes
//! `(1/6) Σ p0 · ((pi − p0) × (p{i+1} − p0))`. Faces with holes, and edges without
//! vertices, are refused rather than approximated.

use gensolid_core::{Error, Result};
use gensolid_math::{Aabb, Point3};

use crate::entity::{BodyId, BodyKind, FaceId, ShellId};
use crate::model::Model;

/// Exact vertex-extrema box for planar bodies with straight edges.
/// Curved bodies are rejected until curve and trimmed-face extrema are implemented.
pub fn body_bounding_box(model: &Model, body: BodyId) -> Result<Aabb> {
    require_polyhedral(model, body)?;
    let points = model
        .body_vertices(body)?
        .into_iter()
        .map(|id| model.vertex(id).map(|vertex| vertex.point))
        .collect::<Result<Vec<_>>>()?;
    Aabb::from_points(points).ok_or_else(|| Error::degenerate("body has no vertices"))
}

/// Volume of the solid regions of `body`.
///
/// The boundary of each solid region is walked once. A face whose back side is that
/// region keeps its loop order (the face normal points out of the solid); a face whose
/// front side is that region is reversed.
pub fn body_volume(model: &Model, body: BodyId) -> Result<f64> {
    require_polyhedral(model, body)?;
    let kind = model.body(body)?.kind;
    if kind != BodyKind::Solid {
        return Err(Error::unsupported(
            "volume is implemented for solid bodies only",
        ));
    }
    let regions = model.body(body)?.regions.clone();
    let mut volume = 0.0;
    let mut faces = 0usize;
    for region_id in regions {
        let region = model.region(region_id)?;
        if !region.solid {
            continue;
        }
        let shells = region.shells.clone();
        for shell_id in shells {
            let face_ids = model.shell(shell_id)?.faces.clone();
            for face_id in face_ids {
                volume += oriented_face_volume(model, face_id, shell_id)?;
                faces += 1;
            }
        }
    }
    if faces == 0 {
        return Err(Error::degenerate("solid body has no faces"));
    }
    Ok(volume)
}

fn require_polyhedral(model: &Model, body: BodyId) -> Result<()> {
    for face in model.body_faces(body)? {
        let surface = model
            .face(face)?
            .surface
            .ok_or_else(|| Error::invalid_argument("face has no surface"))?;
        if !matches!(
            model.geom().surface(surface)?,
            gensolid_geom::Surface::Plane(_)
        ) {
            return Err(Error::unimplemented(
                "mass properties and bounding boxes of curved faces are not implemented",
            ));
        }
    }
    for edge in model.body_edges(body)? {
        let curve = model
            .edge(edge)?
            .curve
            .ok_or_else(|| Error::invalid_argument("edge has no curve"))?;
        if !matches!(model.geom().curve(curve)?, gensolid_geom::Curve::Line(_)) {
            return Err(Error::unimplemented(
                "mass properties and bounding boxes of curved edges are not implemented",
            ));
        }
    }
    Ok(())
}

fn oriented_face_volume(model: &Model, face_id: FaceId, shell_id: ShellId) -> Result<f64> {
    let face = model.face(face_id)?;
    let outward = face.back == Some(shell_id);
    let inward = face.front == Some(shell_id);
    if outward == inward {
        return Err(Error::new(
            gensolid_core::ErrorCode::InvalidArgument,
            "face does not bound this shell on exactly one side",
        ));
    }
    if face.loops.len() != 1 {
        return Err(Error::unimplemented(
            "volume of a face with holes, or with no loop, is not implemented",
        ));
    }
    let loop_id = face.loops[0];
    let fins = model.loop_(loop_id)?.fins.clone();
    if fins.len() < 3 {
        return Err(Error::degenerate("face loop has fewer than 3 fins"));
    }
    let mut points = Vec::with_capacity(fins.len());
    for fin_id in fins {
        let (start, _) = model.fin_vertices(fin_id)?;
        let Some(vertex_id) = start else {
            return Err(Error::unimplemented(
                "volume of an edge without vertices is not implemented",
            ));
        };
        points.push(model.vertex(vertex_id)?.point);
    }
    if inward {
        points.reverse();
    }
    Ok(polygon_volume(&points))
}

/// Signed volume of the cone from the origin to the polygon.
fn polygon_volume(points: &[Point3]) -> f64 {
    let Some(&origin) = points.first() else {
        return 0.0;
    };
    let mut acc = 0.0;
    for window in points[1..].windows(2) {
        let a = window[0] - origin;
        let b = window[1] - origin;
        acc += origin.to_vec().dot(a.cross(b));
    }
    acc / 6.0
}

#[cfg(test)]
mod tests {
    use gensolid_core::Precision;
    use gensolid_math::{Frame, Point3, Vec3};
    use proptest::prelude::*;

    use super::*;
    use crate::make_block;
    use crate::model::Model;

    fn volume_tolerance(x: f64, y: f64, z: f64, precision: &Precision) -> f64 {
        let surface_area = 2.0 * (x * y + y * z + z * x);
        precision.linear() * surface_area
    }

    #[test]
    fn axis_aligned_block_has_exact_box_and_volume() {
        let precision = Precision::DEFAULT;
        let mut model = Model::new(precision);
        let body = make_block(&mut model, 2.0, 4.0, 8.0, &Frame::WORLD).unwrap();
        let box_ = body_bounding_box(&model, body).unwrap();
        let tol = precision.linear();
        assert!((box_.min().x + 1.0).abs() <= tol);
        assert!((box_.min().y + 2.0).abs() <= tol);
        assert!((box_.min().z + 4.0).abs() <= tol);
        assert!((box_.max().x - 1.0).abs() <= tol);
        assert!((box_.max().y - 2.0).abs() <= tol);
        assert!((box_.max().z - 4.0).abs() <= tol);
        let volume = body_volume(&model, body).unwrap();
        assert!((volume - 64.0).abs() <= volume_tolerance(2.0, 4.0, 8.0, &precision));
    }

    proptest! {
        #[test]
        fn rotated_block_volume_is_the_size_product(
            x in 0.01f64..100.0,
            y in 0.01f64..100.0,
            z in 0.01f64..100.0,
            ox in -50.0f64..50.0,
            oy in -50.0f64..50.0,
            oz in -50.0f64..50.0,
            ax in -1.0f64..1.0,
            ay in -1.0f64..1.0,
            az in 0.1f64..1.0,
        ) {
            let precision = Precision::DEFAULT;
            let frame = Frame::from_axis(Point3::new(ox, oy, oz), Vec3::new(ax, ay, az), &precision).unwrap();
            let mut model = Model::new(precision);
            let body = make_block(&mut model, x, y, z, &frame).unwrap();
            let volume = body_volume(&model, body).unwrap();
            let expected = x * y * z;
            prop_assert!((volume - expected).abs() <= volume_tolerance(x, y, z, &precision));
            let box_ = body_bounding_box(&model, body).unwrap();
            prop_assert!(box_.diagonal().x >= 0.0 && box_.diagonal().y >= 0.0 && box_.diagonal().z >= 0.0);
        }
    }
}
