//! Direct NURBS parameter-domain intersection. Interval exclusion and Krawczyk
//! contraction classify regular isolated curve/surface roots without tessellation.
//! Unresolved domains remain in the report; they never become a complete result.
use crate::nurbs_bounds::{bound_curve, bound_surface, knot_spans};
use crate::{BCurve, BSurface, ParamRange};
use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Interval, Point3, Sign, Vec3};

/// Work limits; these are resource controls, not geometric tolerances.
#[derive(Clone, Copy, Debug)]
pub struct IntersectionOptions {
    /// Maximum interval cells examined.
    pub max_cells: usize,
    /// Maximum subdivision depth per knot cell.
    pub max_depth: usize,
}
impl Default for IntersectionOptions {
    fn default() -> Self {
        Self {
            max_cells: 20000,
            max_depth: 64,
        }
    }
}
/// A regular isolated root with an existence/uniqueness interval certificate.
#[derive(Clone, Debug)]
pub struct NurbsRoot {
    /// Curve parameter followed by surface U and V.
    pub parameters: [f64; 3],
    /// Enclosing parameter box where Krawczyk proves a unique root.
    pub enclosure: [ParamRange; 3],
    /// Evaluated curve point; the certified root lies within `point_error_bound`.
    pub point: Point3,
    /// Residual between original curve and original surface evaluations.
    pub residual: f64,
    /// Upper bound on distance from `point` to the exact isolated root.
    pub point_error_bound: f64,
}
/// Whole-domain curve/surface result, including unresolved singular/boundary cells.
#[derive(Clone, Debug)]
pub struct NurbsIntersectionReport {
    /// Certified isolated roots; duplicate certificates are merged only by uniqueness.
    pub roots: Vec<NurbsRoot>,
    /// Boxes not excluded or uniquely resolved, including any exhausted work queue.
    pub unresolved: Vec<[ParamRange; 3]>,
    /// Examined interval cells.
    pub visited_cells: usize,
    /// Interval-excluded cells.
    pub excluded_cells: usize,
}
impl NurbsIntersectionReport {
    /// True only when the entire input parameter domain is accounted for.
    pub fn complete(&self) -> bool {
        self.unresolved.is_empty()
    }
}
fn interval(r: ParamRange) -> Result<Interval> {
    Interval::new(r.lo, r.hi).ok_or_else(|| Error::invalid_argument("invalid parameter box"))
}
fn middle(r: ParamRange) -> f64 {
    r.lo + (r.hi - r.lo) * 0.5
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn magnitude(i: Interval) -> f64 {
    i.lo().abs().max(i.hi().abs())
}
pub(crate) fn inverse(a: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let mut aug = [[0.0; 6]; 3];
    for i in 0..3 {
        for j in 0..3 {
            aug[i][j] = a[i][j];
        }
        aug[i][i + 3] = 1.0;
    }
    for i in 0..3 {
        // Pivot selection only chooses a numerical preconditioner. All root decisions
        // below still require outward-rounded interval inclusion and contraction.
        let pivot = (i..3).max_by(|&a, &b| aug[a][i].abs().total_cmp(&aug[b][i].abs()))?;
        aug.swap(i, pivot);
        let divisor = aug[i][i];
        for value in &mut aug[i] {
            *value /= divisor;
        }
        let pivot_row = aug[i];
        for (row, values) in aug.iter_mut().enumerate() {
            if row == i {
                continue;
            }
            let factor = values[i];
            for (value, pivot) in values.iter_mut().zip(pivot_row) {
                *value -= factor * pivot;
            }
        }
    }
    let result: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| aug[i][j + 3]));
    result
        .iter()
        .flatten()
        .all(|v| v.is_finite())
        .then_some(result)
}
fn numeric(c: &BCurve, s: &BSurface, x: [f64; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let a = c.eval(x[0]);
    let b = s.eval(x[1], x[2]);
    let columns = [components(a.d1), components(-b.du), components(-b.dv)];
    (
        components(a.point - b.point),
        std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i])),
    )
}
fn functions(
    c: &BCurve,
    s: &BSurface,
    x: [ParamRange; 3],
) -> Result<([Interval; 3], [[Interval; 3]; 3])> {
    let a = bound_curve(c, x[0])?;
    let b = bound_surface(s, x[1], x[2])?;
    Ok((
        std::array::from_fn(|i| a.point[i] - b.point[i]),
        std::array::from_fn(|i| [a.du[i], -b.du[i], -b.dv[i]]),
    ))
}
struct Contract {
    image: [Interval; 3],
    unique: bool,
}
fn contract(c: &BCurve, s: &BSurface, x: [ParamRange; 3]) -> Result<Option<Contract>> {
    let center = x.map(middle);
    let (_, j) = numeric(c, s, center);
    let Some(pre) = inverse(j) else {
        return Ok(None);
    };
    let (_, bounds) = functions(c, s, x)?;
    let (f, _) = functions(c, s, center.map(|t| ParamRange::bounded(t, t)))?;
    let mut image = [Interval::point(0.0); 3];
    let mut contraction = true;
    for i in 0..3 {
        image[i] = Interval::point(center[i]);
        for k in 0..3 {
            image[i] = image[i] - Interval::point(pre[i][k]) * f[k];
        }
        let mut norm = Interval::point(0.0);
        for j in 0..3 {
            let mut a = Interval::point(if i == j { 1.0 } else { 0.0 });
            for k in 0..3 {
                a = a - Interval::point(pre[i][k]) * bounds[k][j];
            }
            norm = norm + Interval::point(magnitude(a));
            image[i] = image[i] + a * (interval(x[j])? - Interval::point(center[j]));
        }
        contraction &= (Interval::point(1.0) - norm).sign() == Some(Sign::Positive);
    }
    let inclusion = (0..3).all(|i| image[i].lo() > x[i].lo && image[i].hi() < x[i].hi);
    Ok(Some(Contract {
        image,
        unique: inclusion && contraction,
    }))
}
fn refined(
    c: &BCurve,
    s: &BSurface,
    mut x: [f64; 3],
    cell: [ParamRange; 3],
    p: &Precision,
) -> Option<([f64; 3], Point3, f64)> {
    for _ in 0..32 {
        let (f, j) = numeric(c, s, x);
        let a = c.eval(x[0]).point;
        let b = s.eval(x[1], x[2]).point;
        let residual = a.distance(b);
        if residual <= p.solver() {
            return Some((x, a, residual));
        }
        let pre = inverse(j)?;
        let next =
            std::array::from_fn(|i| x[i] - pre[i].iter().zip(f).map(|(a, b)| a * b).sum::<f64>());
        if !(0..3).all(|i| next[i].is_finite() && next[i] >= cell[i].lo && next[i] <= cell[i].hi) {
            return None;
        }
        x = next;
    }
    None
}
fn location_bound(
    c: &BCurve,
    s: &BSurface,
    mut cell: [ParamRange; 3],
    point: Point3,
) -> Result<f64> {
    for _ in 0..16 {
        let Some(k) = contract(c, s, cell)? else {
            break;
        };
        for (r, image) in cell.iter_mut().zip(k.image) {
            r.lo = r.lo.max(image.lo());
            r.hi = r.hi.min(image.hi());
        }
    }
    let b = bound_curve(c, cell[0])?;
    let q = [point.x, point.y, point.z];
    let squared = (0..3).fold(Interval::point(0.0), |sum, i| {
        sum + (b.point[i] - Interval::point(q[i])).sqr()
    });
    Ok(squared.hi().sqrt().next_up())
}
fn overlaps(a: [ParamRange; 3], b: [ParamRange; 3]) -> bool {
    (0..3).all(|i| a[i].lo <= b[i].hi && b[i].lo <= a[i].hi)
}
/// Intersect arbitrary-degree positive-weight NURBS curve and surface on all knot cells.
/// Complete reports certify exclusion or a unique root in every region. Tangencies,
/// overlap, boundary roots or exhausted budgets may leave explicit unresolved cells.
pub fn intersect_nurbs_curve_surface(
    c: &BCurve,
    s: &BSurface,
    p: &Precision,
    options: IntersectionOptions,
) -> Result<NurbsIntersectionReport> {
    if options.max_cells == 0 || options.max_depth == 0 {
        return Err(Error::invalid_argument(
            "intersection budgets must be nonzero",
        ));
    }
    let (a, b) = s.degrees();
    let (n, m) = s.counts();
    let (uk, vk) = s.knots();
    let mut report = NurbsIntersectionReport {
        roots: Vec::new(),
        unresolved: Vec::new(),
        visited_cells: 0,
        excluded_cells: 0,
    };
    let domains = [c.range(), s.ranges().0, s.ranges().1];
    let spans = [
        knot_spans(c.knots(), c.degree(), c.points().len()),
        knot_spans(uk, a, n),
        knot_spans(vk, b, m),
    ];
    let count = spans.iter().try_fold(1usize, |n, s| n.checked_mul(s.len()));
    if count.is_none_or(|n| n > options.max_cells) {
        report.unresolved.push(domains);
        return Ok(report);
    }
    let mut queue = Vec::new();
    for &t in &spans[0] {
        for &u in &spans[1] {
            for &v in &spans[2] {
                queue.push(([t, u, v], 0usize));
            }
        }
    }
    while let Some((mut cell, depth)) = queue.pop() {
        if report.visited_cells >= options.max_cells {
            report.unresolved.push(cell);
            report.unresolved.extend(queue.into_iter().map(|x| x.0));
            break;
        }
        report.visited_cells += 1;
        let (f, _) = functions(c, s, cell)?;
        if f.iter().any(|x| !x.contains(0.0)) {
            report.excluded_cells += 1;
            continue;
        }
        if let Some(k) = contract(c, s, cell)? {
            if (0..3).any(|i| k.image[i].hi() < cell[i].lo || k.image[i].lo() > cell[i].hi) {
                report.excluded_cells += 1;
                continue;
            }
            if k.unique
                && let Some((parameters, point, residual)) = refined(
                    c,
                    s,
                    k.image.map(|i| i.lo() + (i.hi() - i.lo()) * 0.5),
                    cell,
                    p,
                )
            {
                let point_error_bound = location_bound(c, s, cell, point)?;
                if point_error_bound > p.linear() {
                    report.unresolved.push(cell);
                    continue;
                }
                let root = NurbsRoot {
                    parameters,
                    enclosure: cell,
                    point,
                    residual,
                    point_error_bound,
                };
                let mut duplicate = false;
                let mut ambiguous = false;
                for old in &report.roots {
                    if overlaps(old.enclosure, cell) {
                        let hull = std::array::from_fn(|i| {
                            ParamRange::bounded(
                                old.enclosure[i].lo.min(cell[i].lo),
                                old.enclosure[i].hi.max(cell[i].hi),
                            )
                        });
                        if contract(c, s, hull)
                            .ok()
                            .flatten()
                            .is_some_and(|k| k.unique)
                        {
                            duplicate = true;
                        } else {
                            ambiguous = true;
                        }
                    }
                }
                if !ambiguous {
                    if !duplicate {
                        report.roots.push(root);
                    }
                    continue;
                }
            }
            // Root-preserving interval contraction; excluded portions contain no root.
            for (i, r) in cell.iter_mut().enumerate() {
                r.lo = r.lo.max(k.image[i].lo());
                r.hi = r.hi.min(k.image[i].hi());
            }
        }
        if depth >= options.max_depth {
            report.unresolved.push(cell);
            continue;
        }
        let axis = (0..3)
            .max_by(|&i, &j| {
                ((cell[i].hi - cell[i].lo) / (domains[i].hi - domains[i].lo))
                    .total_cmp(&((cell[j].hi - cell[j].lo) / (domains[j].hi - domains[j].lo)))
            })
            .ok_or_else(|| Error::internal("no split axis"))?;
        let r = cell[axis];
        let mid = r.lo + (r.hi - r.lo) * 0.47;
        if mid <= r.lo || mid >= r.hi {
            report.unresolved.push(cell);
            continue;
        }
        let mut right = cell;
        cell[axis].hi = mid;
        right[axis].lo = mid;
        queue.push((right, depth + 1));
        queue.push((cell, depth + 1));
    }
    report
        .roots
        .sort_by(|a, b| a.parameters[0].total_cmp(&b.parameters[0]));
    Ok(report)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    fn patch() -> BSurface {
        BSurface::new(
            3,
            3,
            [vec![0.0; 4], vec![1.0; 4]].concat(),
            [vec![0.0; 4], vec![1.0; 4]].concat(),
            4,
            4,
            (0..4)
                .flat_map(|i| (0..4).map(move |j| Point3::new(i as f64 / 3.0, j as f64 / 3.0, 0.0)))
                .collect(),
            None,
        )
        .unwrap()
    }
    #[test]
    fn arbitrary_degree_curve_surface_reports_three_certified_roots() {
        let z = [-0.08, 0.14, -0.14, 0.08];
        let c = BCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            (0..4)
                .map(|i| Point3::new(i as f64 / 3.0, 0.4, z[i]))
                .collect(),
            None,
        )
        .unwrap();
        let r = intersect_nurbs_curve_surface(
            &c,
            &patch(),
            &Precision::DEFAULT,
            IntersectionOptions::default(),
        )
        .unwrap();
        assert!(
            r.complete(),
            "visited {}, unresolved {}",
            r.visited_cells,
            r.unresolved.len()
        );
        assert_eq!(r.roots.len(), 3);
        for (root, t) in r.roots.iter().zip([0.2, 0.5, 0.8]) {
            assert!((root.parameters[0] - t).abs() <= Precision::DEFAULT.linear());
            assert!(root.point_error_bound <= Precision::DEFAULT.linear());
            assert!(root.point.distance(Point3::new(t, 0.4, 0.0)) <= root.point_error_bound);
        }
    }
    #[test]
    fn exhausted_and_overlap_queries_never_claim_completeness() {
        let c = BCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![Point3::new(0.2, 0.4, 0.0), Point3::new(0.8, 0.4, 0.0)],
            None,
        )
        .unwrap();
        let r = intersect_nurbs_curve_surface(
            &c,
            &patch(),
            &Precision::DEFAULT,
            IntersectionOptions {
                max_cells: 16,
                max_depth: 8,
            },
        )
        .unwrap();
        assert!(!r.complete());
        assert!(r.roots.is_empty());
    }

    proptest::proptest! {
        #![proptest_config(proptest::test_runner::Config::with_cases(12))]
        #[test]
        fn rational_multispan_roots_in_random_frames(weight in 0.5f64..2.0, ox in -2.0f64..2.0, ax in -1.0f64..1.0) {
            let p=Precision::DEFAULT;
            let frame=gensolid_math::Frame::from_axis(Point3::new(ox,0.0,0.0),Vec3::new(ax,1.0,1.0),&p).unwrap();
            let transform=|q:Point3|frame.origin() + frame.x()*q.x + frame.y()*q.y + frame.z()*q.z;
            // Two linear rational knot spans cross the plane, each at a known point.
            let c=BCurve::new(1,vec![0.0,0.0,0.5,1.0,1.0],vec![Point3::new(0.2,0.4,-1.0),Point3::new(0.5,0.4,1.0),Point3::new(0.8,0.4,-1.0)].into_iter().map(transform).collect(),Some(vec![1.0,weight,1.0])).unwrap();
            let base=patch();let(uk,vk)=base.knots();
            let s=BSurface::new(3,3,uk.to_vec(),vk.to_vec(),4,4,base.points().iter().copied().map(transform).collect(),None).unwrap();
            let result=intersect_nurbs_curve_surface(&c,&s,&p,Default::default()).unwrap();
            proptest::prop_assert!(result.complete());proptest::prop_assert_eq!(result.roots.len(),2);
            for (root,x) in result.roots.iter().zip([0.35,0.65]) {
                proptest::prop_assert!(root.point.distance(transform(Point3::new(x,0.4,0.0)))<=p.linear());
                proptest::prop_assert!(root.point_error_bound<=p.linear());
            }
            let budget=intersect_nurbs_curve_surface(&c,&s,&p,IntersectionOptions{max_cells:1,max_depth:1}).unwrap();
            proptest::prop_assert!(!budget.complete());proptest::prop_assert_eq!(budget.unresolved.len(),1);
        }
    }
}
