//! Analytic boundary assembly; polygons here are B-rep loops, never mesh facets.
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{BCurve, Curve, Line, Surface};
use gensolid_math::{Interval, Point3, Sign};
use gensolid_topo::*;

#[derive(Clone)]
pub(crate) struct BoundaryEdge {
    pub curve: BCurve,
    pub uv: Option<BCurve>,
}
#[derive(Clone)]
pub(crate) struct BoundaryFace {
    pub surface: Surface,
    pub sense: Sense,
    pub edges: Vec<BoundaryEdge>,
}

pub(crate) fn line(a: Point3, b: Point3) -> Result<BCurve> {
    BCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![a, b], None)
}
pub(crate) fn near(a: Point3, b: Point3, p: &Precision) -> Result<bool> {
    let d = [a.x, a.y, a.z]
        .into_iter()
        .zip([b.x, b.y, b.z])
        .map(|(a, b)| (Interval::point(a) - Interval::point(b)).sqr())
        .fold(Interval::point(0.0), |a, b| a + b);
    match (d - Interval::point(p.default_local()).sqr()).sign() {
        Some(Sign::Negative | Sign::Zero) => Ok(true),
        Some(Sign::Positive) => Ok(false),
        None => Err(Error::unimplemented(
            "vertex coincidence unresolved at local tolerance",
        )),
    }
}
pub(crate) fn same_curve(a: &BCurve, b: &BCurve, p: &Precision) -> Result<bool> {
    // A common positive rational basis gives a convex-hull error bound. Different
    // bases are refused rather than equating edges from their endpoints alone.
    if a.degree() != b.degree() || a.points().len() != b.points().len() {
        return Ok(false);
    }
    let ar = a.range();
    let br = b.range();
    let ak: Vec<_> = a
        .knots()
        .iter()
        .map(|t| (t - ar.lo) / (ar.hi - ar.lo))
        .collect();
    let bk: Vec<_> = b
        .knots()
        .iter()
        .map(|t| (t - br.lo) / (br.hi - br.lo))
        .collect();
    if ak != bk {
        return Ok(false);
    }
    let aw: Vec<_> = (0..a.points().len())
        .map(|i| a.weights().map_or(1.0, |w| w[i]))
        .collect();
    let bw: Vec<_> = (0..b.points().len())
        .map(|i| b.weights().map_or(1.0, |w| w[i]))
        .collect();
    // Representation equality, not a geometric sign decision.
    if aw.iter().map(|w| w / aw[0]).collect::<Vec<_>>()
        != bw.iter().map(|w| w / bw[0]).collect::<Vec<_>>()
    {
        return Ok(false);
    }
    for (&a, &b) in a.points().iter().zip(b.points()) {
        if !near(a, b, p)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn carrier(c: &BCurve, p: &Precision) -> Result<Curve> {
    // A degree-one, equal-weight span is an exact straight segment. Preserve its
    // analytic type so mass integration and subsequent planar operations can use it.
    if c.degree() == 1
        && c.points().len() == 2
        && c.weights()
            .is_none_or(|w| w.windows(2).all(|v| v[0].to_bits() == v[1].to_bits()))
    {
        Ok(Curve::Line(Line::through(c.points()[0], c.points()[1], p)?))
    } else {
        Ok(Curve::BCurve(c.clone()))
    }
}
fn carrier_range(c: &BCurve, p: &Precision) -> Result<gensolid_geom::ParamRange> {
    match carrier(c, p)? {
        Curve::Line(line) => Ok(gensolid_geom::ParamRange::bounded(
            0.0,
            line.closest_param(c.points()[1]),
        )),
        _ => Ok(c.range()),
    }
}

pub(crate) fn assemble(model: &mut Model, faces: &[BoundaryFace], solid: bool) -> Result<BodyId> {
    if faces.is_empty() {
        return Err(Error::invalid_argument("empty boundary"));
    }
    let p = *model.precision();
    let mut points = Vec::new();
    let mut corners = Vec::new();
    for face in faces {
        if face.edges.len() < 3 {
            return Err(Error::degenerate("face requires three boundary segments"));
        }
        let mut ids = Vec::new();
        for (i, e) in face.edges.iter().enumerate() {
            let a = e.curve.eval(e.curve.range().lo).point;
            let b = e.curve.eval(e.curve.range().hi).point;
            let next = &face.edges[(i + 1) % face.edges.len()].curve;
            if !near(b, next.eval(next.range().lo).point, &p)? {
                return Err(Error::invalid_argument("open geometric boundary"));
            }
            if near(a, b, &p)? {
                return Err(Error::degenerate("collapsed boundary edge"));
            }
            if !a.is_finite() || [a.x, a.y, a.z].iter().any(|x| x.abs() > p.size() * 0.5) {
                return Err(Error::invalid_argument("boundary outside session size"));
            }
            let mut found = None;
            for (j, &point) in points.iter().enumerate() {
                if near(a, point, &p)? {
                    if found.is_some() {
                        return Err(Error::unimplemented("ambiguous sewing vertex cluster"));
                    }
                    found = Some(j);
                }
            }
            let index = match found {
                Some(j) => j,
                None => {
                    points.push(a);
                    points.len() - 1
                }
            };
            if ids.contains(&index) {
                return Err(Error::unimplemented("periodic or repeated boundary vertex"));
            }
            ids.push(index);
        }
        corners.push(ids);
    }
    let mut raw = model.raw();
    let body = raw.insert_body(Body {
        kind: if solid {
            BodyKind::Solid
        } else {
            BodyKind::Sheet
        },
        regions: Vec::new(),
    })?;
    let empty = raw.insert_region(Region {
        body,
        solid: false,
        shells: Vec::new(),
    })?;
    let front = raw.insert_shell(Shell {
        region: empty,
        faces: Vec::new(),
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    })?;
    raw.region_mut(empty)?.shells.push(front);
    raw.body_mut(body)?.regions.push(empty);
    let back = if solid {
        let material = raw.insert_region(Region {
            body,
            solid: true,
            shells: Vec::new(),
        })?;
        let shell = raw.insert_shell(Shell {
            region: material,
            faces: Vec::new(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        })?;
        raw.region_mut(material)?.shells.push(shell);
        raw.body_mut(body)?.regions.push(material);
        Some(shell)
    } else {
        None
    };
    let vertices = points
        .iter()
        .map(|&point| {
            raw.insert_vertex(Vertex {
                point,
                tolerance: None,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut edges: Vec<(usize, usize, BCurve, EdgeId, usize)> = Vec::new();
    for (data, ids) in faces.iter().zip(&corners) {
        let surface = raw.geom_mut().add_surface(data.surface.clone())?;
        let face = raw.insert_face(Face {
            front: Some(front),
            back,
            loops: Vec::new(),
            surface: Some(surface),
            sense: data.sense,
        })?;
        let lp = raw.insert_loop(Loop {
            face,
            fins: Vec::new(),
        })?;
        for (k, entry) in data.edges.iter().enumerate() {
            let a = ids[k];
            let b = ids[(k + 1) % ids.len()];
            let mut found = None;
            for (index, (c, d, curve, _, _)) in edges.iter().enumerate() {
                if (a == *c && b == *d) || (a == *d && b == *c) {
                    let reverse = a == *d;
                    let candidate = if reverse {
                        entry.curve.reversed()?
                    } else {
                        entry.curve.clone()
                    };
                    if !same_curve(curve, &candidate, &p)? {
                        return Err(Error::unimplemented(
                            "sewing curves have incompatible representations",
                        ));
                    }
                    if found.is_some() {
                        return Err(Error::invalid_argument("ambiguous shared edge"));
                    }
                    found = Some((index, reverse));
                }
            }
            let (edge, sense, range) = if let Some((index, reverse)) = found {
                let (_, _, curve, edge, uses) = &mut edges[index];
                if *uses >= 2 || !reverse {
                    return Err(Error::invalid_argument(
                        "nonmanifold or inconsistently oriented sewing",
                    ));
                }
                *uses += 1;
                (*edge, Sense::Negative, carrier_range(curve, &p)?)
            } else {
                let curve = raw.geom_mut().add_curve(carrier(&entry.curve, &p)?)?;
                let edge = raw.insert_edge(Edge {
                    fins: Vec::new(),
                    start: Some(vertices[a]),
                    end: Some(vertices[b]),
                    curve: Some(curve),
                    tolerance: None,
                })?;
                edges.push((a, b, entry.curve.clone(), edge, 1));
                (edge, Sense::Positive, carrier_range(&entry.curve, &p)?)
            };
            let fin = raw.insert_fin(Fin {
                owner: lp,
                edge,
                sense,
            })?;
            raw.edge_mut(edge)?.fins.push(fin);
            raw.loop_mut(lp)?.fins.push(fin);
            if let Some(uv) = &entry.uv {
                let uv = if sense == Sense::Negative {
                    uv.reversed()?
                } else {
                    uv.clone()
                };
                raw.set_fin_geometry(
                    fin,
                    FinGeometry {
                        pcurve: uv,
                        edge_range: range,
                    },
                )?;
            }
        }
        raw.face_mut(face)?.loops.push(lp);
        raw.shell_mut(front)?.faces.push(face);
        if let Some(back) = back {
            raw.shell_mut(back)?.faces.push(face);
        }
    }
    if solid && edges.iter().any(|e| e.4 != 2) {
        return Err(Error::invalid_argument(
            "sewn shell has free boundary edges",
        ));
    }
    Ok(body)
}

pub(crate) fn extract(model: &Model, face: FaceId) -> Result<BoundaryFace> {
    let f = model.face(face)?;
    let [lp] = f.loops.as_slice() else {
        return Err(Error::unimplemented(
            "boundary extraction requires one loop",
        ));
    };
    let surface = model
        .geom()
        .surface(
            f.surface
                .ok_or_else(|| Error::invalid_argument("missing face surface"))?,
        )?
        .clone();
    let mut edges = Vec::new();
    for &fin in &model.loop_(*lp)?.fins {
        let value = model.fin(fin)?;
        let edge = model.edge(value.edge)?;
        let curve = model.geom().curve(
            edge.curve
                .ok_or_else(|| Error::invalid_argument("missing edge curve"))?,
        )?;
        let (Some(a), Some(b)) = (edge.start, edge.end) else {
            return Err(Error::unimplemented("ring boundary extraction"));
        };
        let mut curve = match curve {
            Curve::BCurve(c) => {
                if !near(
                    c.eval(c.range().lo).point,
                    model.vertex(a)?.point,
                    model.precision(),
                )? || !near(
                    c.eval(c.range().hi).point,
                    model.vertex(b)?.point,
                    model.precision(),
                )? {
                    return Err(Error::unimplemented(
                        "B-curve edge is not its full carrier domain",
                    ));
                }
                c.clone()
            }
            Curve::Line(_) => line(model.vertex(a)?.point, model.vertex(b)?.point)?,
            _ => return Err(Error::unimplemented("analytic boundary conversion")),
        };
        let mut uv = model.fin_geometry(fin)?.map(|g| g.pcurve.clone());
        if value.sense == Sense::Negative {
            curve = curve.reversed()?;
            uv = uv.map(|c| c.reversed()).transpose()?;
        }
        if uv.is_none() {
            match &surface {
                Surface::Plane(plane) => {
                    let points = curve
                        .points()
                        .iter()
                        .map(|&p| {
                            let q = plane.frame().to_local(p);
                            Point3::new(q.x, q.y, 0.0)
                        })
                        .collect();
                    uv = Some(BCurve::new(
                        curve.degree(),
                        curve.knots().to_vec(),
                        points,
                        curve.weights().map(|w| w.to_vec()),
                    )?);
                }
                Surface::BSurface(s) => {
                    let canonical = crate::patch::patch_face(s)?;
                    for edge in canonical.edges {
                        if same_curve(&curve, &edge.curve, model.precision())? {
                            uv = edge.uv;
                            break;
                        }
                        if same_curve(&curve, &edge.curve.reversed()?, model.precision())? {
                            uv = edge.uv.map(|c| c.reversed()).transpose()?;
                            break;
                        }
                    }
                }
                _ => (),
            }
        }
        edges.push(BoundaryEdge { curve, uv });
    }
    Ok(BoundaryFace {
        surface,
        sense: f.sense,
        edges,
    })
}

pub(crate) fn reverse(face: &BoundaryFace) -> Result<BoundaryFace> {
    Ok(BoundaryFace {
        surface: face.surface.clone(),
        sense: face.sense.reversed(),
        edges: face
            .edges
            .iter()
            .rev()
            .map(|e| {
                Ok(BoundaryEdge {
                    curve: e.curve.reversed()?,
                    uv: e.uv.as_ref().map(|c| c.reversed()).transpose()?,
                })
            })
            .collect::<Result<_>>()?,
    })
}
