//! JSON scene bridge. Display tessellation is output-only and never used by modelling.
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::nurbs_intersect::{IntersectionOptions, intersect_nurbs_curve_surface};
use gensolid_geom::nurbs_surface_intersect::intersect_nurbs_surfaces;
use gensolid_geom::{BCurve, BSurface, Curve, Surface};
use gensolid_math::{Frame, Point2, Point3};
use gensolid_topo::{BodyId, LoftSection, Model};
use serde_json::{Value, json};
use std::io::{self, Read};
fn xyz(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn number(v: &Value, key: &str, default: f64, lo: f64, hi: f64) -> Result<f64> {
    let x = match v.get(key) {
        None => default,
        Some(value) => value
            .as_f64()
            .ok_or_else(|| Error::invalid_argument(format!("{key} must be numeric")))?,
    };
    if !x.is_finite() || x < lo || x > hi {
        return Err(Error::invalid_argument(format!(
            "{key} must lie in [{lo}, {hi}]"
        )));
    }
    Ok(x)
}
fn mesh(surface: &BSurface, name: &str, color: &str) -> Value {
    let (u, v) = surface.ranges();
    let n = 24;
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            vertices.push(xyz(surface
                .eval(
                    u.lo + (u.hi - u.lo) * i as f64 / n as f64,
                    v.lo + (v.hi - v.lo) * j as f64 / n as f64,
                )
                .point));
        }
    }
    for i in 0..n {
        for j in 0..n {
            let a = i * (n + 1) + j;
            triangles.extend([[a, a + n + 1, a + n + 2], [a, a + n + 2, a + 1]]);
        }
    }
    json!({"name":name,"color":color,"vertices":vertices,"triangles":triangles})
}
fn line(points: Vec<[f64; 3]>, color: &str) -> Value {
    json!({"points":points,"color":color})
}
fn body_scene(model: &Model, body: BodyId, name: &str) -> Result<Value> {
    let mut objects = Vec::new();
    let mut lines = Vec::new();
    for face in model.body_faces(body)? {
        let f = model.face(face)?;
        let surface = model.geom().surface(
            f.surface
                .ok_or_else(|| Error::internal("missing face surface"))?,
        )?;
        if let Surface::BSurface(s) = surface {
            objects.push(mesh(s, name, "#5bd5c5"));
        } else {
            let mut vertices = Vec::new();
            for fin in &model.loop_(f.loops[0])?.fins {
                let id = model
                    .fin_vertices(*fin)?
                    .0
                    .ok_or_else(|| Error::internal("missing vertex"))?;
                vertices.push(xyz(model.vertex(id)?.point));
            }
            let triangles = (1..vertices.len() - 1)
                .map(|i| [0, i, i + 1])
                .collect::<Vec<_>>();
            objects.push(
                json!({"name":"cap","color":"#91b9ee","vertices":vertices,"triangles":triangles}),
            );
        }
    }
    for edge in model.body_edges(body)? {
        let e = model.edge(edge)?;
        let c = model
            .geom()
            .curve(e.curve.ok_or_else(|| Error::internal("missing curve"))?)?;
        let points = match c {
            Curve::BCurve(c) => (0..=32)
                .map(|i| {
                    xyz(c
                        .eval(c.range().lo + (c.range().hi - c.range().lo) * i as f64 / 32.0)
                        .point)
                })
                .collect(),
            _ => [e.start, e.end]
                .into_iter()
                .flatten()
                .map(|v| model.vertex(v).map(|p| xyz(p.point)))
                .collect::<Result<Vec<_>>>()?,
        };
        lines.push(line(points, "#163d50"));
    }
    let checked = gensolid_check::check_body(model, body, &Default::default())?;
    Ok(
        json!({"objects":objects,"lines":lines,"points":[],"status":"checked","summary":"完整 B-rep 检查通过","metrics":{"faces":model.body_faces(body)?.len(),"edges":model.body_edges(body)?.len(),"vertices":model.body_vertices(body)?.len(),"checker_faults":checked.faults.len(),"unchecked":checked.unchecked},"notes":["实体保存精确 NURBS 几何；显示网格仅用于交互预览，没有弦高认证。"]}),
    )
}
fn plane_patch() -> Result<BSurface> {
    BSurface::new(
        3,
        3,
        [vec![0.0; 4], vec![1.0; 4]].concat(),
        [vec![0.0; 4], vec![1.0; 4]].concat(),
        4,
        4,
        (0..4)
            .flat_map(|i| {
                (0..4).map(move |j| {
                    Point3::new(
                        i as f64 / 3.0,
                        j as f64 / 3.0,
                        if i == 3 {
                            0.3 * (j as f64 / 3.0 - 0.4)
                        } else {
                            0.0
                        },
                    )
                })
            })
            .collect(),
        None,
    )
}
fn choose(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |a, i| a * (n - i) as f64 / (i + 1) as f64)
}
fn curve_case() -> Result<BCurve> {
    // z=10(t-.2)(t-.5)(t-.8)(1+.2t^2): genuinely degree five, three known roots.
    let c = [-0.8, 6.6, -15.16, 11.32, -3.0, 2.0];
    let points = (0..6)
        .map(|i| {
            let z = (0..=i).map(|k| c[k] * choose(i, k) / choose(5, k)).sum();
            Point3::new(i as f64 / 5.0, 0.4, z)
        })
        .collect();
    BCurve::new(5, [vec![0.0; 6], vec![1.0; 6]].concat(), points, None)
}
fn scene(input: &Value) -> Result<Value> {
    if !input.is_object() {
        return Err(Error::invalid_argument("scene request must be an object"));
    }
    let p = Precision::DEFAULT;
    let mode = input["scene"].as_str().unwrap_or("revolve");
    let mut value = match mode {
        "revolve" => {
            let r = number(input, "radius", 1.0, 0.1, 3.0)?;
            let h = number(input, "height", 2.0, 0.2, 5.0)?;
            let angle = number(input, "angle", 270.0, 10.0, 360.0)?.to_radians();
            let profile = vec![
                Point2::new(r, -h * 0.5),
                Point2::new(r + 0.7, -h * 0.5),
                Point2::new(r + 1.0, h * 0.2),
                Point2::new(r + 0.5, h * 0.5),
                Point2::new(r, h * 0.5),
            ];
            let mut model = Model::new(p);
            let result = gensolid_modeling::body_create_solid_revolution(
                &mut model,
                &profile,
                angle,
                &Frame::WORLD,
            )?;
            body_scene(&model, result.body, "revolution")?
        }
        "loft" => {
            let twist = number(input, "twist", 25.0, -65.0, 65.0)?.to_radians();
            let scale = number(input, "scale", 0.6, 0.3, 1.7)?;
            let count = number(input, "sections", 4.0, 2.0, 6.0)?;
            if count.fract() != 0.0 {
                return Err(Error::invalid_argument("section count must be an integer"));
            }
            let count = count as usize;
            let sections = (0..count)
                .map(|j| {
                    let t = j as f64 / (count - 1) as f64;
                    LoftSection {
                        height: t * 3.0 - 1.5,
                        profile: (0..6)
                            .map(|i| {
                                let a = std::f64::consts::TAU * i as f64 / 6.0 + t * twist;
                                let rx = 1.0 + (scale - 1.0) * t;
                                let ry = 1.0 + 0.25 * (t * std::f64::consts::PI).sin();
                                Point2::new(rx * a.cos() + t * 0.6, ry * a.sin())
                            })
                            .collect(),
                    }
                })
                .collect::<Vec<_>>();
            let mut model = Model::new(p);
            let result = gensolid_modeling::body_create_solid_loft_sections(
                &mut model,
                &sections,
                &Frame::WORLD,
            )?;
            body_scene(&model, result.body, "multi-section loft")?
        }
        "smooth_loft" => {
            let bend = number(input, "bend", 1.0, -2.0, 2.0)?;
            let sections = (0..4)
                .map(|j| {
                    BCurve::new(
                        3,
                        [vec![0.0; 4], vec![1.0; 4]].concat(),
                        (0..4)
                            .map(|i| {
                                Point3::new(
                                    i as f64 - 1.5,
                                    j as f64 - 1.5,
                                    bend * ((i + j) as f64).sin(),
                                )
                            })
                            .collect(),
                        Some(vec![1.0, 0.85, 0.85, 1.0]),
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            let parameters = [0.0, 0.2, 0.7, 1.0];
            let surface = gensolid_geom::loft::loft_nurbs_sections(&sections, &parameters, 3, &p)?;
            let mut residual: f64 = 0.0;
            for (c, v) in sections.iter().zip(parameters) {
                for i in 0..=32 {
                    let u = i as f64 / 32.0;
                    residual = residual.max(c.eval(u).point.distance(surface.eval(u, v).point));
                }
            }
            json!({"objects":[mesh(&surface,"smooth rational loft","#5bd5c5")],"lines":sections.iter().map(|c|line((0..=64).map(|i|xyz(c.eval(i as f64/64.0).point)).collect(),"#ffd48a")).collect::<Vec<_>>(),"points":[],"status":"geometry","summary":"三次有理 NURBS 插值放样曲面","metrics":{"sections":4,"u_degree":3,"v_degree":3,"section_residual":residual},"notes":["黄色线为原始截面；曲面通过全部截面。","这是曲面几何，尚未封闭为实体；不显示 body checker 通过标记。"]})
        }
        "curve_surface" => {
            let curve = curve_case()?;
            let surface = plane_patch()?;
            let report = intersect_nurbs_curve_surface(
                &curve,
                &surface,
                &p,
                IntersectionOptions::default(),
            )?;
            json!({"objects":[mesh(&surface,"bicubic patch","#91b9ee")],"lines":[line((0..=128).map(|i|xyz(curve.eval(i as f64/128.0).point)).collect(),"#ffd48a")],"points":report.roots.iter().map(|r|xyz(r.point)).collect::<Vec<_>>(),"status":if report.complete(){"certified"}else{"unresolved"},"summary":if report.complete(){"全参数域求交已认证"}else{"仍有未解决参数区域"},"metrics":{"roots":report.roots.len(),"expected_roots":3,"unresolved_cells":report.unresolved.len(),"visited_cells":report.visited_cells,"max_residual":report.roots.iter().map(|r|r.residual).fold(0.0,f64::max),"root_position_bound":report.roots.iter().map(|r|r.point_error_bound).fold(0.0,f64::max)},"notes":["五次曲线 × 双三次 NURBS 曲面；区间排除与 Krawczyk 唯一根认证。"],"root_parameters":report.roots.iter().map(|r|r.parameters).collect::<Vec<_>>()})
        }
        "surface_surface" => {
            let radius = number(input, "radius", 0.65, 0.2, 0.9)?;
            let coefficients = [1.0, -1.0, 1.0];
            let a = BSurface::new(
                2,
                2,
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                3,
                3,
                (0..3)
                    .flat_map(|i| {
                        (0..3).map(move |j| {
                            Point3::new(
                                i as f64 - 1.0,
                                j as f64 - 1.0,
                                0.7 * (coefficients[i] + coefficients[j] - radius * radius),
                            )
                        })
                    })
                    .collect(),
                None,
            )?;
            let b = BSurface::new(
                1,
                1,
                vec![0.0, 0.0, 1.0, 1.0],
                vec![0.0, 0.0, 1.0, 1.0],
                2,
                2,
                vec![
                    Point3::new(-1.0, -1.0, 0.0),
                    Point3::new(-1.0, 1.0, 0.0),
                    Point3::new(1.0, -1.0, 0.0),
                    Point3::new(1.0, 1.0, 0.0),
                ],
                None,
            )?;
            let report = intersect_nurbs_surfaces(
                &a,
                &b,
                &p,
                IntersectionOptions {
                    max_cells: 8000,
                    max_depth: 12,
                },
            )?;
            json!({"objects":[mesh(&a,"curved NURBS patch","#5bd5c5"),mesh(&b,"planar NURBS patch","#91b9ee")],"lines":report.curves.iter().map(|c|line(c.samples.iter().map(|s|xyz(s.point)).collect(),"#ffbd69")).collect::<Vec<_>>(),"points":[],"status":"candidate","summary":"数值交线 · 全局分支连接尚未认证","metrics":{"candidate_curves":report.curves.len(),"possible_cells":report.possible_cells.len(),"visited_cells":report.visited_cells,"budget_exhausted":report.budget_exhausted,"max_residual":report.curves.iter().flat_map(|c|c.samples.iter().map(|s|s.residual)).fold(0.0,f64::max)},"notes":["黄色交线直接由原始 NURBS 参数方程求解；未用三角网格求交。","区间可能区域覆盖未排除的交集；小环、相切和奇异点的连接完整性未认证。"]})
        }
        _ => return Err(Error::invalid_argument("unknown viewer scene")),
    };
    value["scene"] = json!(mode);
    value["precision"] = json!({"linear":p.linear(),"angular":p.angular()});
    value["request"] = input.clone();
    Ok(value)
}
fn main() {
    let mut input = String::new();
    let result = io::stdin()
        .take(65537)
        .read_to_string(&mut input)
        .map_err(|e| Error::parse(e.to_string()))
        .and_then(|_| {
            if input.len() > 65536 {
                return Err(Error::invalid_argument("scene request too large"));
            }
            let input: Value =
                serde_json::from_str(&input).map_err(|e| Error::parse(e.to_string()))?;
            scene(&input)
        });
    match result {
        Ok(value) => println!("{value}"),
        Err(error) => {
            println!(
                "{}",
                json!({"status":"error","error_code":format!("{:?}",error.code()),"message":error.message()})
            );
            std::process::exit(1);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scene_bridge_uses_real_checked_geometry() -> Result<()> {
        let r = scene(&json!({"scene":"revolve","angle":180.0}))?;
        assert_eq!(r["status"], "checked");
        assert_eq!(r["metrics"]["unchecked"], 0);
        let r = scene(&json!({"scene":"curve_surface"}))?;
        assert_eq!(r["status"], "certified");
        assert_eq!(r["metrics"]["roots"], 3);
        Ok(())
    }
}
