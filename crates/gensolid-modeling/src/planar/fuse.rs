//! Remove internal interfaces from a finite arrangement of convex planar cells.
use super::*;
use gensolid_topo::{FaceId, Sense};
use std::collections::{BTreeMap, BTreeSet};

fn split_cells(mut cells: Vec<Poly>, planes: &[Plane], p: &Precision) -> Result<Vec<Poly>> {
    for &plane in planes {
        let reverse = opposite(plane, p)?;
        let mut next = Vec::new();
        for cell in cells {
            if let Some(a) = clip(&cell, plane, p)? {
                next.push(a);
            }
            if let Some(b) = clip(&cell, reverse, p)? {
                next.push(b);
            }
            if next.len() > 512 {
                return Err(Error::unimplemented(
                    "planar arrangement cell budget exceeded",
                ));
            }
        }
        cells = next;
    }
    Ok(cells)
}
fn same_cycle(a: &[usize], b: &[usize]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    (0..b.len()).any(|shift| {
        a.iter()
            .enumerate()
            .all(|(i, &v)| v == b[(shift + b.len() - i) % b.len()])
    })
}
struct Use {
    cell: usize,
    plane: Plane,
    vertices: Vec<usize>,
}
fn conform(points: &[Point3], polygon: &[usize], p: &Precision) -> Result<Vec<usize>> {
    let mut out = Vec::new();
    for (i, &start) in polygon.iter().enumerate() {
        let end = polygon[(i + 1) % polygon.len()];
        let a = points[start];
        let b = points[end];
        let d = b - a;
        let length = d.length();
        let component = |b, a| Interval::point(b) - Interval::point(a);
        let edge = [
            component(b.x, a.x),
            component(b.y, a.y),
            component(b.z, a.z),
        ];
        let l2 = edge[0].sqr() + edge[1].sqr() + edge[2].sqr();
        if (l2 - Interval::point(p.linear()).sqr()).sign() != Some(Sign::Positive) {
            return Err(Error::degenerate("arrangement edge collapsed"));
        }
        let mut inside = Vec::new();
        for (k, &q) in points.iter().enumerate() {
            if k == start || k == end {
                continue;
            }
            let delta = [
                Interval::point(q.x) - Interval::point(a.x),
                Interval::point(q.y) - Interval::point(a.y),
                Interval::point(q.z) - Interval::point(a.z),
            ];
            let product =
                |x: [Interval; 3], y: [Interval; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
            let projection = product(delta, edge);
            if (projection - Interval::point(p.linear() * length)).sign() != Some(Sign::Positive)
                || (l2 - projection - Interval::point(p.linear() * length)).sign()
                    != Some(Sign::Positive)
            {
                continue;
            }
            let cross = [
                delta[1] * edge[2] - delta[2] * edge[1],
                delta[2] * edge[0] - delta[0] * edge[2],
                delta[0] * edge[1] - delta[1] * edge[0],
            ];
            let separation = product(cross, cross) - Interval::point(p.default_local()).sqr() * l2;
            match separation.sign() {
                Some(Sign::Negative | Sign::Zero) => inside.push((projection, k)),
                Some(Sign::Positive) => (),
                None => {
                    return Err(Error::unimplemented(
                        "arrangement edge incidence unresolved",
                    ));
                }
            }
        }
        // Sorting proposes an order; disjoint interval projections must certify
        // every neighboring pair before that order is used in the boundary loop.
        inside.sort_by(|a, b| a.0.lo().total_cmp(&b.0.lo()));
        for pair in inside.windows(2) {
            if (pair[1].0 - pair[0].0).sign() != Some(Sign::Positive) {
                return Err(Error::unimplemented("arrangement vertex order unresolved"));
            }
        }
        out.push(start);
        out.extend(inside.into_iter().map(|(_, k)| k));
    }
    Ok(out)
}
fn components(adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut seen = BTreeSet::new();
    let mut groups = Vec::new();
    for root in 0..adjacency.len() {
        if seen.contains(&root) {
            continue;
        }
        let mut work = vec![root];
        let mut group = Vec::new();
        seen.insert(root);
        while let Some(i) = work.pop() {
            group.push(i);
            for &j in &adjacency[i] {
                if seen.insert(j) {
                    work.push(j);
                }
            }
        }
        groups.push(group);
    }
    groups
}
fn volume_sign(faces: &[&Use], points: &[Point3]) -> Result<Sign> {
    let origin = points[faces[0].vertices[0]];
    let vector = |p: Point3| {
        [
            Interval::point(p.x) - Interval::point(origin.x),
            Interval::point(p.y) - Interval::point(origin.y),
            Interval::point(p.z) - Interval::point(origin.z),
        ]
    };
    let mut volume = Interval::point(0.0);
    for f in faces {
        let a = vector(points[f.vertices[0]]);
        for i in 1..f.vertices.len() - 1 {
            let b = vector(points[f.vertices[i]]);
            let c = vector(points[f.vertices[i + 1]]);
            volume = volume
                + a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
        }
    }
    match volume.sign() {
        Some(Sign::Positive) => Ok(Sign::Positive),
        Some(Sign::Negative) => Ok(Sign::Negative),
        _ => Err(Error::unimplemented(
            "boundary shell orientation is not certified",
        )),
    }
}
fn make_face(f: &Use, points: &[Point3]) -> Result<boundary::BoundaryFace> {
    let frame = f.plane.frame();
    let mut edges = Vec::new();
    for (i, &a) in f.vertices.iter().enumerate() {
        let b = f.vertices[(i + 1) % f.vertices.len()];
        let u = frame.to_local(points[a]);
        let v = frame.to_local(points[b]);
        edges.push(boundary::BoundaryEdge {
            curve: boundary::line(points[a], points[b])?,
            uv: Some(boundary::line(
                Point3::new(u.x, u.y, 0.0),
                Point3::new(v.x, v.y, 0.0),
            )?),
        });
    }
    Ok(boundary::BoundaryFace {
        surface: gensolid_geom::Surface::Plane(f.plane),
        sense: Sense::Positive,
        edges,
    })
}
fn cavity_shells(
    m: &mut Model,
    body: BodyId,
    outer_faces: &[FaceId],
    cavities: &[Vec<FaceId>],
) -> Result<()> {
    let regions = m.body(body)?.regions.clone();
    let outer_front = m.region(regions[0])?.shells[0];
    let outer_back = m.region(regions[1])?.shells[0];
    let mut raw = m.raw();
    raw.shell_mut(outer_front)?.faces = outer_faces.to_vec();
    raw.shell_mut(outer_back)?.faces = outer_faces.to_vec();
    for faces in cavities {
        let region = raw.insert_region(Region {
            body,
            solid: false,
            shells: Vec::new(),
        })?;
        let front = raw.insert_shell(Shell {
            region,
            faces: faces.clone(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        })?;
        let back = raw.insert_shell(Shell {
            region: regions[1],
            faces: faces.clone(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        })?;
        raw.body_mut(body)?.regions.push(region);
        raw.region_mut(region)?.shells.push(front);
        raw.region_mut(regions[1])?.shells.push(back);
        for &face in faces {
            let f = raw.face_mut(face)?;
            f.front = Some(front);
            f.back = Some(back);
        }
    }
    Ok(())
}
fn fuse(m: &mut Model, cells: &[Poly]) -> Result<Vec<BodyId>> {
    if cells.is_empty() {
        return Ok(Vec::new());
    }
    let p = *m.precision();
    let mut points = Vec::new();
    let mut uses = Vec::new();
    for (cell, poly) in cells.iter().enumerate() {
        for face in &poly.faces {
            let mut vertices = Vec::new();
            for &q in &face.points {
                let mut found = None;
                for (i, &r) in points.iter().enumerate() {
                    if boundary::near(q, r, &p)? {
                        if found.is_some() {
                            return Err(Error::unimplemented("ambiguous arrangement vertex"));
                        }
                        found = Some(i);
                    }
                }
                vertices.push(match found {
                    Some(i) => i,
                    None => {
                        points.push(q);
                        points.len() - 1
                    }
                });
            }
            uses.push(Use {
                cell,
                plane: face.plane,
                vertices,
            });
        }
    }
    if points.len() > 4096 {
        return Err(Error::unimplemented("arrangement vertex budget exceeded"));
    }
    let mut facets: BTreeMap<Vec<usize>, Vec<usize>> = BTreeMap::new();
    for (i, f) in uses.iter_mut().enumerate() {
        f.vertices = conform(&points, &f.vertices, &p)?;
        let mut key = f.vertices.clone();
        key.sort_unstable();
        facets.entry(key).or_default().push(i);
    }
    let mut adjacency = vec![Vec::new(); cells.len()];
    let mut exterior = Vec::new();
    for indices in facets.values() {
        match indices.as_slice() {
            [a] => exterior.push(*a),
            [a, b] => {
                let a = &uses[*a];
                let b = &uses[*b];
                if a.cell == b.cell || !same_cycle(&a.vertices, &b.vertices) {
                    return Err(Error::unimplemented(
                        "arrangement interface orientation mismatch",
                    ));
                }
                adjacency[a.cell].push(b.cell);
                adjacency[b.cell].push(a.cell);
            }
            _ => return Err(Error::unimplemented("nonmanifold arrangement interface")),
        }
    }
    let mut results = Vec::new();
    for cells in components(&adjacency) {
        let faces: Vec<_> = exterior
            .iter()
            .map(|&i| &uses[i])
            .filter(|f| cells.contains(&f.cell))
            .collect();
        let mut edge_uses: BTreeMap<(usize, usize), Vec<(usize, bool)>> = BTreeMap::new();
        for (i, f) in faces.iter().enumerate() {
            for (k, &a) in f.vertices.iter().enumerate() {
                let b = f.vertices[(k + 1) % f.vertices.len()];
                edge_uses
                    .entry((a.min(b), a.max(b)))
                    .or_default()
                    .push((i, a < b));
            }
        }
        let mut neighbors = vec![Vec::new(); faces.len()];
        for uses in edge_uses.values() {
            let [(a, sa), (b, sb)] = uses.as_slice() else {
                return Err(Error::unimplemented(
                    "Boolean boundary is not a closed two-manifold",
                ));
            };
            if a == b || sa == sb {
                return Err(Error::unimplemented(
                    "Boolean boundary edge winding mismatch",
                ));
            }
            neighbors[*a].push(*b);
            neighbors[*b].push(*a);
        }
        let mut outer = None;
        let mut cavities = Vec::new();
        for group in components(&neighbors) {
            let shell = group.iter().map(|&i| faces[i]).collect::<Vec<_>>();
            if volume_sign(&shell, &points)? == Sign::Positive {
                if outer.is_some() {
                    return Err(Error::unimplemented(
                        "multiple outer shells in one material component",
                    ));
                }
                outer = Some(group);
            } else {
                cavities.push(group);
            }
        }
        let outer =
            outer.ok_or_else(|| Error::unimplemented("Boolean component has no outward shell"))?;
        let mut groups = vec![outer];
        groups.extend(cavities);
        let ordered = groups
            .iter()
            .flatten()
            .map(|&i| make_face(faces[i], &points))
            .collect::<Result<Vec<_>>>()?;
        let body = boundary::assemble(m, &ordered, true)?;
        let ids = m.body_faces(body)?;
        let mut offset = groups[0].len();
        let mut inner_ids = Vec::new();
        for group in &groups[1..] {
            inner_ids.push(ids[offset..offset + group.len()].to_vec());
            offset += group.len();
        }
        cavity_shells(m, body, &ids[..groups[0].len()], &inner_ids)?;
        body_check_complete(m, body)?;
        results.push(body);
    }
    Ok(results)
}

/// Convex planar input Boolean with internal cell interfaces removed. Returns one
/// body per material component, including finite cavity regions when necessary.
/// Coplanar exterior subdivisions are retained as faces. Edge/point-only nonmanifold
/// connections and unresolved arrangements are rejected; curved inputs are unsupported.
pub fn body_boolean(
    model: &mut Model,
    first: BodyId,
    second: BodyId,
    operation: BooleanOperation,
) -> Result<Vec<BodyId>> {
    let a = read(model, first)?;
    let b = read(model, second)?;
    let p = *model.precision();
    let planes: Vec<_> = a.faces.iter().chain(&b.faces).map(|f| f.plane).collect();
    let cells = match operation {
        BooleanOperation::Difference => difference(&a, &b, &p)?,
        BooleanOperation::Union => {
            let mut cells = difference(&a, &b, &p)?;
            cells.push(b);
            cells
        }
        BooleanOperation::Intersection => {
            let mut cell = Some(a);
            for face in b.faces {
                if let Some(c) = cell {
                    cell = clip(&c, face.plane, &p)?;
                } else {
                    break;
                }
            }
            cell.into_iter().collect()
        }
    };
    let arrangement = split_cells(cells, &planes, &p)?;
    model.transaction(|m| fuse(m, &arrangement))
}
