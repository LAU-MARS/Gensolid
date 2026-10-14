//! Checked surface edits and topology reconstruction.
use crate::{
    body_check_complete,
    boundary::{self, BoundaryEdge, BoundaryFace},
    patch,
};
use gensolid_core::{Error, Result};
use gensolid_geom::{BSurface, Surface};
use gensolid_math::{Interval, Point3, Sign, Vec3};
use gensolid_topo::{BodyId, BodyKind, EdgeId, FaceId, Model, Sense};

fn translate(s: &BSurface, d: Vec3) -> Result<BSurface> {
    let (pu, pv) = s.degrees();
    let (nu, nv) = s.counts();
    let (uk, vk) = s.knots();
    BSurface::new(
        pu,
        pv,
        uk.to_vec(),
        vk.to_vec(),
        nu,
        nv,
        s.points().iter().map(|&p| p + d).collect(),
        s.weights().map(|w| w.to_vec()),
    )
}

/// Solid between a regular graph patch and its translated copy along the graph
/// projection normal. This is a directional prism, not a normal surface offset.
pub fn body_create_solid_surface_prism(
    model: &mut Model,
    surface: &BSurface,
    height: f64,
) -> Result<BodyId> {
    let p = *model.precision();
    if !height.is_finite()
        || (Interval::point(height) - Interval::point(p.linear())).sign() != Some(Sign::Positive)
    {
        return Err(Error::degenerate("prism height must exceed precision"));
    }
    patch::certify_graph(surface, &p)?;
    let (u, v) = surface.ranges();
    let mid = surface.eval((u.lo + u.hi) * 0.5, (v.lo + v.hi) * 0.5);
    let direction = mid
        .du
        .cross(mid.dv)
        .normalized(p.angular())
        .ok_or_else(|| Error::degenerate("singular prism normal"))?
        * height;
    let top = translate(surface, direction)?;
    patch::certify_graph(&top, &p)?;
    let bottom = patch::patch_face(surface)?;
    let upper = patch::patch_face(&top)?;
    let mut faces = vec![boundary::reverse(&bottom)?, upper];
    for edge in &bottom.edges {
        let c = &edge.curve;
        let mut points = Vec::new();
        let mut weights = Vec::new();
        for (i, &point) in c.points().iter().enumerate() {
            points.extend([point, point + direction]);
            let w = c.weights().map_or(1.0, |w| w[i]);
            weights.extend([w, w]);
        }
        let side = BSurface::new(
            c.degree(),
            1,
            c.knots().to_vec(),
            vec![0.0, 0.0, 1.0, 1.0],
            c.points().len(),
            2,
            points,
            Some(weights),
        )?;
        faces.push(patch::patch_face(&side)?);
    }
    model.transaction(|m| {
        let body = boundary::assemble(m, &faces, true)?;
        body_check_complete(m, body)?;
        Ok(body)
    })
}

/// Copy a body's faces to separately checked sheets; source topology is retained.
pub fn body_extract_sheets(model: &mut Model, body: BodyId) -> Result<Vec<BodyId>> {
    body_check_complete(model, body)?;
    let data = model
        .body_faces(body)?
        .into_iter()
        .map(|f| boundary::extract(model, f))
        .collect::<Result<Vec<_>>>()?;
    model.transaction(|m| {
        data.into_iter()
            .map(|f| {
                let b = boundary::assemble(m, &[f], false)?;
                body_check_complete(m, b)?;
                Ok(b)
            })
            .collect()
    })
}

/// Split a single-span NURBS face along an isoparameter, propagating both boundary
/// subdivisions into every adjacent face. Returns a new body; the source is unchanged.
/// Topological relations, not endpoint-only proximity, select the adjacent uses.
pub fn body_split_face_isoparametric(
    model: &mut Model,
    body: BodyId,
    face: FaceId,
    constant_u: bool,
    fraction: f64,
) -> Result<BodyId> {
    let (original, canonical, correspondence) = rectangle_boundary(model, body, face)?;
    let Surface::BSurface(surface) = &original.surface else {
        return Err(Error::internal("rectangle surface missing"));
    };
    let (left, right) = patch::split_surface(surface, constant_u, fraction)?;
    let replacements = [patch::patch_face(&left)?, patch::patch_face(&right)?];
    let lp = model.face(face)?.loops[0];
    let source_fins = &model.loop_(lp)?.fins;
    let mut cuts = Vec::new();
    // Work in the actual loop order and infer each crossing from UV isoparameters.
    for (i, &k) in correspondence.iter().enumerate() {
        let uv = canonical.edges[k]
            .uv
            .as_ref()
            .ok_or_else(|| Error::internal("missing canonical UV"))?;
        let a = uv.eval(uv.range().lo).point;
        let b = uv.eval(uv.range().hi).point;
        let (x, y) = if constant_u { (a.x, b.x) } else { (a.y, b.y) };
        // Canonical UV endpoints reuse the very same parameter bound.
        if x.to_bits() == y.to_bits() {
            continue;
        }
        let domain = if constant_u {
            surface.ranges().0
        } else {
            surface.ranges().1
        };
        let value = domain.lo + (domain.hi - domain.lo) * fraction;
        let f = (value - x) / (y - x);
        let fin = model.fin(source_fins[i])?;
        cuts.push((
            fin.edge,
            if fin.sense == Sense::Positive {
                f
            } else {
                1.0 - f
            },
        ));
    }
    if cuts.len() != 2 {
        return Err(Error::internal(
            "rectangle cut needs two boundary intersections",
        ));
    }
    replace_face_boundary(model, body, face, original.sense, replacements, cuts)
}

fn rectangle_boundary(
    model: &Model,
    body: BodyId,
    face: FaceId,
) -> Result<(BoundaryFace, BoundaryFace, Vec<usize>)> {
    body_check_complete(model, body)?;
    let source_body = model.body(body)?;
    let expected_regions = match source_body.kind {
        BodyKind::Sheet => 1,
        BodyKind::Solid => 2,
        _ => {
            return Err(Error::unimplemented(
                "face subdivision supports sheets and manifold solids",
            ));
        }
    };
    if source_body.regions.len() != expected_regions
        || source_body
            .regions
            .iter()
            .any(|&r| model.region(r).is_ok_and(|r| r.shells.len() != 1))
    {
        return Err(Error::unimplemented(
            "face subdivision of multiple shells or finite cavities",
        ));
    }
    if !model.body_faces(body)?.contains(&face) {
        return Err(Error::invalid_argument("face does not belong to body"));
    }
    let original = boundary::extract(model, face)?;
    let Surface::BSurface(surface) = &original.surface else {
        return Err(Error::unimplemented("isoparametric split needs NURBS face"));
    };
    let canonical = patch::patch_face(surface)?;
    let canonical = if original.sense == Sense::Negative {
        boundary::reverse(&canonical)?
    } else {
        canonical
    };
    if original.edges.len() != 4 {
        return Err(Error::unimplemented(
            "already trimmed face needs general trim intersection",
        ));
    }
    // Ensure the source is the complete rectangle, allowing a cyclic loop start.
    let p = *model.precision();
    let mut correspondence = Vec::new();
    for edge in &original.edges {
        let a = edge.curve.eval(edge.curve.range().lo).point;
        let b = edge.curve.eval(edge.curve.range().hi).point;
        let mut index = None;
        for (k, e) in canonical.edges.iter().enumerate() {
            if boundary::near(a, e.curve.eval(e.curve.range().lo).point, &p)?
                && boundary::near(b, e.curve.eval(e.curve.range().hi).point, &p)?
            {
                index = Some(k);
            }
        }
        let index = index
            .ok_or_else(|| Error::unimplemented("face boundary is not the full NURBS rectangle"))?;
        if !boundary::same_curve(&edge.curve, &canonical.edges[index].curve, &p)? {
            return Err(Error::unimplemented(
                "face boundary carrier differs from the NURBS rectangle",
            ));
        }
        correspondence.push(index);
    }
    let unique = correspondence
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    if unique.len() != 4 {
        return Err(Error::invalid_argument(
            "rectangle loop repeats a boundary edge",
        ));
    }
    Ok((original, canonical, correspondence))
}

fn replace_face_boundary(
    model: &mut Model,
    body: BodyId,
    face: FaceId,
    sense: Sense,
    replacements: [BoundaryFace; 2],
    cuts: Vec<(EdgeId, f64)>,
) -> Result<BodyId> {
    let mut result = Vec::new();
    for id in model.body_faces(body)? {
        if id == face {
            for f in &replacements {
                result.push(if sense == Sense::Negative {
                    boundary::reverse(f)?
                } else {
                    f.clone()
                });
            }
            continue;
        }
        let mut data = boundary::extract(model, id)?;
        let lp = model.face(id)?.loops[0];
        let fins = &model.loop_(lp)?.fins;
        let mut edges = Vec::new();
        for (entry, &fin_id) in data.edges.iter().zip(fins) {
            let fin = model.fin(fin_id)?;
            if let Some((_, t)) = cuts.iter().find(|c| c.0 == fin.edge) {
                let t = if fin.sense == Sense::Positive {
                    *t
                } else {
                    1.0 - *t
                };
                let (a, b) = patch::split_curve(&entry.curve, t)?;
                let uv = entry
                    .uv
                    .as_ref()
                    .map(|u| patch::split_curve(u, t))
                    .transpose()?;
                edges.push(BoundaryEdge {
                    curve: a,
                    uv: uv.as_ref().map(|u| u.0.clone()),
                });
                edges.push(BoundaryEdge {
                    curve: b,
                    uv: uv.map(|u| u.1),
                });
            } else {
                edges.push(entry.clone());
            }
        }
        data.edges = edges;
        result.push(data);
    }
    let solid = model.body(body)?.kind == BodyKind::Solid;
    model.transaction(|m| {
        let b = boundary::assemble(m, &result, solid)?;
        body_check_complete(m, b)?;
        Ok(b)
    })
}

/// Replace a complete rectangle by two trimmed loops on its original carrier.
/// The caller certifies a simple crossing curve and its boundary fractions.
pub(crate) fn split_rectangle_crossing(
    model: &mut Model,
    body: BodyId,
    face: FaceId,
    linear_u: bool,
    crossing: BoundaryEdge,
    fractions: [f64; 2],
) -> Result<BodyId> {
    let (original, canonical, correspondence) = rectangle_boundary(model, body, face)?;
    let positive = if original.sense == Sense::Negative {
        boundary::reverse(&canonical)?
    } else {
        canonical.clone()
    };
    let edge_cuts = if linear_u {
        [(0, fractions[0]), (2, 1.0 - fractions[1])]
    } else {
        [(1, fractions[1]), (3, 1.0 - fractions[0])]
    };
    let split = |index: usize, fraction| -> Result<[BoundaryEdge; 2]> {
        let e = &positive.edges[index];
        let (a, b) = patch::split_curve(&e.curve, fraction)?;
        let uv =
            e.uv.as_ref()
                .ok_or_else(|| Error::internal("missing rectangle UV"))?;
        let (u, v) = patch::split_curve(uv, fraction)?;
        Ok([
            BoundaryEdge {
                curve: a,
                uv: Some(u),
            },
            BoundaryEdge {
                curve: b,
                uv: Some(v),
            },
        ])
    };
    let [a, b] = split(edge_cuts[0].0, edge_cuts[0].1)?;
    let [c, d] = split(edge_cuts[1].0, edge_cuts[1].1)?;
    let reverse = BoundaryEdge {
        curve: crossing.curve.reversed()?,
        uv: crossing.uv.as_ref().map(|c| c.reversed()).transpose()?,
    };
    let [first, second] = if linear_u {
        [
            vec![a, crossing, d, positive.edges[3].clone()],
            vec![b, positive.edges[1].clone(), c, reverse],
        ]
    } else {
        [
            vec![positive.edges[0].clone(), a, reverse, d],
            vec![crossing, b, positive.edges[2].clone(), c],
        ]
    };
    let fins = &model.loop_(model.face(face)?.loops[0])?.fins;
    let mut cuts = Vec::new();
    for (index, fraction) in edge_cuts {
        // Negative faces reverse both loop order and each boundary curve.
        let oriented_index = if original.sense == Sense::Negative {
            3 - index
        } else {
            index
        };
        let local = correspondence
            .iter()
            .position(|&k| k == oriented_index)
            .ok_or_else(|| Error::internal("rectangle boundary correspondence missing"))?;
        let fin = model.fin(fins[local])?;
        let traversal = if original.sense == Sense::Negative {
            1.0 - fraction
        } else {
            fraction
        };
        cuts.push((
            fin.edge,
            if fin.sense == Sense::Negative {
                1.0 - traversal
            } else {
                traversal
            },
        ));
    }
    let replacements = [first, second].map(|edges| BoundaryFace {
        surface: positive.surface.clone(),
        sense: Sense::Positive,
        edges,
    });
    replace_face_boundary(model, body, face, original.sense, replacements, cuts)
}

pub(crate) fn planar_polygon(model: &Model, face: FaceId) -> Result<Vec<Point3>> {
    let data = boundary::extract(model, face)?;
    if !matches!(data.surface, Surface::Plane(_)) {
        return Err(Error::unimplemented("operation requires planar carriers"));
    }
    if data.edges.len() < 3 {
        return Err(Error::degenerate("planar face needs at least three edges"));
    }
    data.edges
        .iter()
        .map(|e| {
            if e.curve.degree() != 1 || e.curve.points().len() != 2 {
                return Err(Error::unimplemented(
                    "operation requires straight boundary edges",
                ));
            }
            Ok(e.curve.eval(e.curve.range().lo).point)
        })
        .collect()
}

/// Sew independent oriented planar sheets into a convex solid. All faces must be
/// supplied; gaps, incompatible shared curves and inward faces are rejected.
/// Curved sewing currently uses constructive operations with known shell geometry.
pub fn body_sew_sheets(model: &mut Model, sheets: &[BodyId]) -> Result<BodyId> {
    if sheets.is_empty() {
        return Err(Error::invalid_argument("no sheets supplied"));
    }
    let mut points = Vec::new();
    let mut polygons = Vec::new();
    for &sheet in sheets {
        if model.body(sheet)?.kind != BodyKind::Sheet {
            return Err(Error::invalid_argument("sewing input must be sheets"));
        }
        body_check_complete(model, sheet)?;
        for face in model.body_faces(sheet)? {
            let mut ids = Vec::new();
            for point in planar_polygon(model, face)? {
                let mut found = None;
                for (i, &p) in points.iter().enumerate() {
                    if boundary::near(point, p, model.precision())? {
                        found = Some(i);
                        break;
                    }
                }
                ids.push(match found {
                    Some(i) => i,
                    None => {
                        points.push(point);
                        points.len() - 1
                    }
                });
            }
            polygons.push(ids);
        }
    }
    model.transaction(|m| {
        let b = gensolid_topo::make_convex_polyhedron(m, &points, &polygons)?;
        body_check_complete(m, b)?;
        Ok(b)
    })
}
