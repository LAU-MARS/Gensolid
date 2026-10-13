//! Polygon extrusion, linear sweep and homothetic planar loft.

use gensolid_core::{Error, Result};
use gensolid_math::{Frame, Interval, Point2, Sign, Vec2, Vec3, predicates::orient2d};

use crate::polyhedron::dot;
use crate::{BodyId, Model, make_convex_polyhedron};

fn counterclockwise(profile: &[Point2]) -> Result<Vec<Point2>> {
    if profile.len() < 3 || profile.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
        return Err(Error::invalid_argument(
            "polygon profile needs at least three finite points",
        ));
    }
    let mut orientation = None;
    for point in &profile[2..] {
        match orient2d(profile[0], profile[1], *point) {
            Some(Sign::Positive | Sign::Negative) => {
                orientation = orient2d(profile[0], profile[1], *point);
                break;
            }
            Some(Sign::Zero) => (),
            None => return Err(Error::degenerate("polygon orientation is not finite")),
        }
    }
    let orientation =
        orientation.ok_or_else(|| Error::degenerate("polygon profile is collinear"))?;
    for (index, start) in profile.iter().enumerate() {
        let end = profile[(index + 1) % profile.len()];
        for point in profile {
            match orient2d(*start, end, *point) {
                Some(Sign::Zero) => (),
                Some(sign) if sign == orientation => (),
                _ => {
                    return Err(Error::invalid_argument(
                        "polygon profile must be simple and convex",
                    ));
                }
            }
        }
    }
    let mut profile = profile.to_vec();
    if orientation == Sign::Negative {
        profile.reverse();
    }
    // Remove redundant, exactly collinear local corners before rounding them into
    // a world frame. Otherwise a rotated straight profile side can acquire a tiny
    // artificial reflex corner from independent floating-point transformations.
    loop {
        let redundant = (0..profile.len()).find(|index| {
            let previous = profile[(*index + profile.len() - 1) % profile.len()];
            let next = profile[(*index + 1) % profile.len()];
            orient2d(previous, profile[*index], next) == Some(Sign::Zero)
        });
        let Some(index) = redundant else {
            break;
        };
        profile.remove(index);
        if profile.len() < 3 {
            return Err(Error::degenerate(
                "polygon profile collapses after removing collinear corners",
            ));
        }
    }
    Ok(profile)
}

fn prism_faces(count: usize, reversed: bool) -> Result<Vec<Vec<usize>>> {
    let top = |index: usize| {
        index
            .checked_add(count)
            .ok_or_else(|| Error::invalid_argument("polygon vertex count overflows"))
    };
    let mut faces = vec![
        (0..count).rev().collect(),
        (0..count).map(top).collect::<Result<Vec<_>>>()?,
    ];
    for i in 0..count {
        let next = (i + 1) % count;
        faces.push(vec![i, next, top(next)?, top(i)?]);
    }
    if reversed {
        for face in &mut faces {
            face.reverse();
        }
    }
    Ok(faces)
}

/// Extrudes a convex polygon from local `z = 0` to `z = height` in `frame`.
/// The profile may be clockwise or counter-clockwise; it must not repeat its first
/// point at the end. Exactly collinear intermediate corners are removed. A higher
/// layer must validate the returned body with `check_body`.
pub fn make_prism(
    model: &mut Model,
    profile: &[Point2],
    height: f64,
    frame: &Frame,
) -> Result<BodyId> {
    if !(height.is_finite()
        && height > model.precision().linear()
        && height < model.precision().size())
    {
        return Err(Error::invalid_argument(
            "prism height must exceed linear precision and lie within the size box",
        ));
    }
    make_linear_sweep(model, profile, frame.z() * height, frame)
}

/// Sweeps a convex polygon by a finite world-coordinate displacement.
/// Its normal component must exceed linear precision. Both positive and negative
/// sweeps are supported, with outward cap and side-loop orientations preserved.
/// This is an unchecked topology-layer constructor; callers must run `check_body`.
pub fn make_linear_sweep(
    model: &mut Model,
    profile: &[Point2],
    displacement: Vec3,
    frame: &Frame,
) -> Result<BodyId> {
    if !displacement.is_finite() {
        return Err(Error::invalid_argument("sweep displacement must be finite"));
    }
    let profile = counterclockwise(profile)?;
    let normal = frame.z();
    let height = dot(
        [
            Interval::point(displacement.x),
            Interval::point(displacement.y),
            Interval::point(displacement.z),
        ],
        [
            Interval::point(normal.x),
            Interval::point(normal.y),
            Interval::point(normal.z),
        ],
    );
    let tolerance = Interval::point(model.precision().linear());
    let reversed = if (height - tolerance).sign() == Some(Sign::Positive) {
        false
    } else if (height + tolerance).sign() == Some(Sign::Negative) {
        true
    } else {
        return Err(Error::degenerate(
            "sweep has no certified normal displacement above linear precision",
        ));
    };
    let mut points = profile
        .iter()
        .map(|p| frame.to_world(Vec3::new(p.x, p.y, 0.0)))
        .collect::<Vec<_>>();
    let top = points.iter().map(|p| *p + displacement).collect::<Vec<_>>();
    points.extend(top);
    let faces = prism_faces(profile.len(), reversed)?;
    make_convex_polyhedron(model, &points, &faces)
}

/// Lofts a convex polygon to a positive uniformly scaled copy in a parallel plane.
///
/// `top_translation` is measured in the local frame's XY axes and `height` along Z.
/// Side faces remain planar because corresponding profile edges are parallel. This
/// is a polyhedral frustum, not a general loft between unrelated profiles. Callers
/// must run `check_body` on the returned topology-layer body.
pub fn make_homothetic_loft(
    model: &mut Model,
    profile: &[Point2],
    height: f64,
    top_scale: f64,
    top_translation: Vec2,
    frame: &Frame,
) -> Result<BodyId> {
    if !(height.is_finite()
        && height > model.precision().linear()
        && height < model.precision().size())
        || !(top_scale.is_finite() && top_scale > 0.0)
        || !top_translation.x.is_finite()
        || !top_translation.y.is_finite()
    {
        return Err(Error::invalid_argument(
            "loft needs a positive height and scale and finite top translation",
        ));
    }
    let profile = counterclockwise(profile)?;
    let mut points = profile
        .iter()
        .map(|p| frame.to_world(Vec3::new(p.x, p.y, 0.0)))
        .collect::<Vec<_>>();
    points.extend(profile.iter().map(|p| {
        frame.to_world(Vec3::new(
            p.x * top_scale + top_translation.x,
            p.y * top_scale + top_translation.y,
            height,
        ))
    }));
    let faces = prism_faces(profile.len(), false)?;
    make_convex_polyhedron(model, &points, &faces)
}
