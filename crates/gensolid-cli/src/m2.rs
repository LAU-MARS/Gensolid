//! Frozen modelling inputs, analytic volume answers and independent OCCT integration.

use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

use gensolid_core::{Error, ErrorCode, Precision, Result, Tolerance};
use gensolid_geom::{Cylinder, ParamRange, Sphere, Surface, Torus};
use gensolid_math::{Frame, Point2, Point3, Vec2, Vec3};
use gensolid_modeling::{MassProperties, body_ask_mass_properties};
use gensolid_topo::{BodyId, Model};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

const SEED: u64 = 20261013;
const FRAMES: usize = 16;

struct Generator(u64);
impl Generator {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
    fn signed(&mut self) -> f64 {
        self.unit() * 2.0 - 1.0
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct Properties {
    volume: f64,
    surface_area: f64,
    centroid: [f64; 3],
    inertia: [[f64; 3]; 3],
}
impl From<MassProperties> for Properties {
    fn from(value: MassProperties) -> Self {
        Self {
            volume: value.volume,
            surface_area: value.surface_area,
            centroid: [value.centroid.x, value.centroid.y, value.centroid.z],
            inertia: value.inertia,
        }
    }
}

#[derive(Serialize)]
struct Boundary {
    id: String,
    points: Vec<[f64; 3]>,
    faces: Vec<Vec<usize>>,
}

struct Work {
    row: usize,
    boundary: Boundary,
    properties: Properties,
    edge_length: f64,
    diameter: f64,
}

#[derive(Deserialize)]
struct OracleCase {
    id: String,
    status: String,
    valid: bool,
    closed: bool,
    properties: Option<Properties>,
    message: Option<String>,
}
#[derive(Deserialize)]
struct OracleBatch {
    available: bool,
    status: String,
    backend: Option<String>,
    message: Option<String>,
    cases: Vec<OracleCase>,
}

fn boundary(model: &Model, body: BodyId, id: String) -> Result<(Boundary, f64, f64)> {
    let vertices = model.body_vertices(body)?;
    let positions = vertices
        .iter()
        .map(|id| model.vertex(*id).map(|vertex| vertex.point))
        .collect::<Result<Vec<_>>>()?;
    let indices: HashMap<_, _> = vertices
        .iter()
        .enumerate()
        .map(|(index, id)| (*id, index))
        .collect();
    let material_shells = model
        .body(body)?
        .regions
        .iter()
        .map(|id| model.region(*id))
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .filter(|region| region.solid)
        .flat_map(|region| region.shells.iter().copied())
        .collect::<Vec<_>>();
    let mut faces = Vec::new();
    for id in model.body_faces(body)? {
        let face = model.face(id)?;
        if face.loops.len() != 1 {
            return Err(Error::unimplemented(
                "M2 same-carrier OCCT reconstruction requires one polygon loop per face",
            ));
        }
        let out = face
            .back
            .is_some_and(|shell| material_shells.contains(&shell));
        let inward = face
            .front
            .is_some_and(|shell| material_shells.contains(&shell));
        if out == inward {
            return Err(Error::invalid_argument(
                "face does not bound one material shell",
            ));
        }
        let lp = face
            .loops
            .first()
            .ok_or_else(|| Error::internal("face loop disappeared"))?;
        let mut polygon = Vec::new();
        for fin in &model.loop_(*lp)?.fins {
            let (start, _) = model.fin_vertices(*fin)?;
            let start = start.ok_or_else(|| {
                Error::unimplemented("M2 reconstruction requires bounded straight edges")
            })?;
            polygon.push(
                *indices
                    .get(&start)
                    .ok_or_else(|| Error::internal("face vertex is absent from body vertices"))?,
            );
        }
        if inward {
            polygon.reverse();
        }
        faces.push(polygon);
    }
    let mut edge_length = 0.0;
    for id in model.body_edges(body)? {
        let edge = model.edge(id)?;
        let (Some(a), Some(b)) = (edge.start, edge.end) else {
            return Err(Error::unimplemented(
                "M2 reconstruction requires two edge vertices",
            ));
        };
        edge_length += model.vertex(a)?.point.distance(model.vertex(b)?.point);
    }
    let mut diameter: f64 = 0.0;
    for a in &positions {
        for b in &positions {
            diameter = diameter.max(a.distance(*b));
        }
    }
    Ok((
        Boundary {
            id,
            points: positions.iter().map(|p| [p.x, p.y, p.z]).collect(),
            faces,
        },
        edge_length,
        diameter,
    ))
}

fn area(profile: &[Point2]) -> f64 {
    profile
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let q = profile[(i + 1) % profile.len()];
            p.x * q.y - p.y * q.x
        })
        .sum::<f64>()
        .abs()
        * 0.5
}

fn invoke(work: &[Work], precision: &Precision) -> Result<OracleBatch> {
    let python = ["python", "python3"]
        .into_iter()
        .find(|name| {
            Command::new(name)
                .arg("--version")
                .output()
                .is_ok_and(|output| output.status.success())
        })
        .ok_or_else(|| Error::unsupported("Python executable not found"))?;
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occt-m2.py");
    let mut child = Command::new(python)
        .arg(script)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|error| Error::internal(format!("starting M2 OCCT helper: {error}")))?;
    let request = json!({"precision":{"linear":precision.linear()}, "cases":work.iter().map(|case| &case.boundary).collect::<Vec<_>>()});
    let bytes = serde_json::to_vec(&request).map_err(|error| Error::internal(error.to_string()))?;
    {
        let mut input = child
            .stdin
            .take()
            .ok_or_else(|| Error::internal("OCCT helper has no stdin pipe"))?;
        input
            .write_all(&bytes)
            .map_err(|error| Error::internal(format!("writing OCCT batch: {error}")))?;
    }
    let output = child
        .wait_with_output()
        .map_err(|error| Error::internal(format!("waiting for OCCT: {error}")))?;
    if !output.status.success() {
        return Err(Error::internal(format!(
            "M2 OCCT helper failed ({}): {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    serde_json::from_slice(&output.stdout)
        .map_err(|error| Error::parse(format!("M2 OCCT returned malformed JSON: {error}")))
}

fn comparison(
    ours: &Properties,
    theirs: &Properties,
    edge_length: f64,
    diameter: f64,
    precision: &Precision,
) -> Value {
    let linear = precision.linear();
    let area = ours.surface_area.abs().max(theirs.surface_area.abs());
    let volume = ours.volume.abs().max(theirs.volume.abs());
    // One linear boundary displacement supplies the volume budget. Area is bounded
    // by displaced polygon edges; centroid and inertia budgets propagate volume and
    // first-moment perturbations over the body's diameter at unit density.
    let volume_budget = linear * area;
    let area_budget = linear * edge_length * 2.0;
    let centroid_budget = linear * (1.0 + area * diameter / volume.max(linear.powi(3)));
    let inertia_budget =
        volume_budget * diameter * diameter + 2.0 * volume * diameter * centroid_budget;
    let volume_error = (ours.volume - theirs.volume).abs();
    let area_error = (ours.surface_area - theirs.surface_area).abs();
    let centroid_error = ours
        .centroid
        .iter()
        .zip(theirs.centroid)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt();
    let mut inertia_error: f64 = 0.0;
    let mut inertia_norm: f64 = 0.0;
    for (a, b) in ours
        .inertia
        .iter()
        .flatten()
        .zip(theirs.inertia.iter().flatten())
    {
        inertia_error = inertia_error.max((a - b).abs());
        inertia_norm = inertia_norm.max(b.abs());
    }
    json!({"absolute_errors":{"volume":volume_error,"surface_area":area_error,"centroid_distance":centroid_error,"inertia_max_component":inertia_error},
        "relative_errors":{"volume":volume_error/theirs.volume.abs().max(linear.powi(3)),"surface_area":area_error/theirs.surface_area.abs().max(linear.powi(2)),
            "centroid_over_body_diameter":centroid_error/diameter.max(linear),"inertia_max_norm":inertia_error/inertia_norm.max(linear.powi(5))},
        "precision_budgets":{"volume":volume_budget,"surface_area":area_budget,"centroid_distance":centroid_budget,"inertia_max_component":inertia_budget},
        "passed":volume_error<=volume_budget && area_error<=area_budget && centroid_error<=centroid_budget && inertia_error<=inertia_budget})
}

fn mesh_samples(precision: &Precision) -> Result<Value> {
    let chord = Tolerance::new(
        (precision.linear() * 65536.0)
            .min(precision.size() * 0.125)
            .max(precision.default_local()),
        precision,
    )?;
    let u = ParamRange::bounded(0.0, std::f64::consts::FRAC_PI_2);
    let v = ParamRange::bounded(-0.5, 0.5);
    let constructors = [
        (
            "sphere",
            Sphere::new(Frame::WORLD, 2.0, precision).map(Surface::Sphere),
            v,
        ),
        (
            "cylinder",
            Cylinder::new(Frame::WORLD, 2.0, precision).map(Surface::Cylinder),
            ParamRange::bounded(0.0, 2.0),
        ),
        (
            "torus",
            Torus::new(Frame::WORLD, 3.0, 1.0, precision).map(Surface::Torus),
            v,
        ),
    ];
    let mut rows = Vec::new();
    let mut passed = 0;
    for (name, surface, v) in constructors {
        match surface.and_then(|surface| {
            gensolid_mesh::surface_tessellate_patch(&surface, u, v, chord, precision)
        }) {
            Ok(mesh) => {
                let pass = mesh.deviation_bound() <= chord.value();
                passed += usize::from(pass);
                rows.push(json!({"carrier":name,"status":"generated","vertices":mesh.vertices.len(),"triangles":mesh.triangles.len(),
                    "u":[u.lo,u.hi],"v":[v.lo,v.hi],"chord_tolerance":chord.value(),"interpolation_bound":mesh.interpolation_bound,
                    "vertex_roundoff_bound":mesh.vertex_roundoff_bound,"guaranteed_deviation_bound":mesh.deviation_bound(),"passed":pass}));
            }
            Err(error) => rows.push(
                json!({"carrier":name,"status":"error","message":error.to_string(),"passed":false}),
            ),
        }
    }
    Ok(
        json!({"scope":"three analytic parameter patches; full trimmed curved bodies and unproven carriers are outside this catalog",
        "cases":3,"passed":passed,"rows":rows}),
    )
}

pub(crate) fn run(precision: &Precision) -> Result<Value> {
    let mut random = Generator(SEED);
    let mut rows = Vec::new();
    let mut work = Vec::new();
    let mut known_answers_passed = 0;
    let mut planar_mesh_passed = 0;
    let planar_chord = Tolerance::new(precision.linear() * 64.0, precision)?;
    for frame_index in 0..FRAMES {
        let scale = 10.0_f64.powi((frame_index % 5) as i32 - 3);
        let origin = Point3::new(
            random.signed() * 20.0,
            random.signed() * 20.0,
            random.signed() * 20.0,
        );
        let frame = Frame::from_axis(
            origin,
            Vec3::new(random.signed(), random.signed(), random.signed()),
            precision,
        )?;
        let count = 3 + (random.unit() * 6.0) as usize;
        let rx = scale * (0.8 + random.unit() * 1.2);
        let ry = scale * (0.8 + random.unit() * 1.2);
        let theta = random.signed() * std::f64::consts::PI;
        let mut profile = (0..count)
            .map(|i| {
                let angle = theta
                    + std::f64::consts::TAU * (i as f64 + random.signed() * 0.2) / count as f64;
                Point2::new(rx * angle.cos(), ry * angle.sin())
            })
            .collect::<Vec<_>>();
        if frame_index % 2 != 0 {
            profile.reverse();
        }
        let height = scale * (0.8 + random.unit() * 1.2);
        let loft_scale = 0.75 + random.unit();
        let translation = Vec2::new(scale * random.signed() * 0.5, scale * random.signed() * 0.5);
        let sign = if frame_index % 2 == 0 { 1.0 } else { -1.0 };
        let displacement =
            frame.vector_to_world(Vec3::new(translation.x, translation.y, height * sign));
        for operation in ["prism", "linear_sweep", "homothetic_loft"] {
            let id = format!("frame-{frame_index:02}-{operation}");
            let input = json!({"operation":operation,"frame":{"origin":[origin.x,origin.y,origin.z],"x":[frame.x().x,frame.x().y,frame.x().z],
                "y":[frame.y().x,frame.y().y,frame.y().z],"z":[frame.z().x,frame.z().y,frame.z().z]},
                "profile":profile.iter().map(|p|[p.x,p.y]).collect::<Vec<_>>(),"height":height,"loft_scale":loft_scale,
                "top_translation_local":[translation.x,translation.y],"sweep_displacement_world":[displacement.x,displacement.y,displacement.z]});
            let expected = area(&profile)
                * height
                * if operation == "homothetic_loft" {
                    (1.0 + loft_scale + loft_scale * loft_scale) / 3.0
                } else {
                    1.0
                };
            let mut model = Model::new(*precision);
            let body = match operation {
                "prism" => {
                    gensolid_modeling::body_create_solid_prism(&mut model, &profile, height, &frame)
                }
                "linear_sweep" => gensolid_modeling::body_create_solid_linear_sweep(
                    &mut model,
                    &profile,
                    displacement,
                    &frame,
                ),
                _ => gensolid_modeling::body_create_solid_homothetic_loft(
                    &mut model,
                    &profile,
                    height,
                    loft_scale,
                    translation,
                    &frame,
                ),
            };
            let result = body.and_then(|body| {
                // Meshing is a separate observation: its failure must never remove
                // an otherwise checked body from the mass comparison denominator.
                let mesh = match gensolid_mesh::body_tessellate_planar(&model, body, planar_chord) {
                    Ok(mesh) => {
                        let passed = mesh.deviation_bound() <= planar_chord.value();
                        planar_mesh_passed += usize::from(passed);
                        json!({"status":"generated","vertices":mesh.vertices.len(),"triangles":mesh.triangles.len(),
                            "chord_tolerance":planar_chord.value(),"interpolation_bound":mesh.interpolation_bound,
                            "vertex_roundoff_bound":mesh.vertex_roundoff_bound,"guaranteed_deviation_bound":mesh.deviation_bound(),"passed":passed})
                    }
                    Err(error) => json!({"status":"error","chord_tolerance":planar_chord.value(),"message":error.to_string(),"passed":false}),
                };
                let properties:Properties=body_ask_mass_properties(&model,body)?.into();
                let (boundary,edge_length,diameter)=boundary(&model,body,id.clone())?;
                Ok((properties,boundary,edge_length,diameter,mesh))
            });
            match result {
                Ok((properties,boundary,edge_length,diameter,mesh))=>{
                    let budget=precision.linear()*properties.surface_area;
                    let error=(properties.volume-expected).abs();
                    let passed=error<=budget;
                    known_answers_passed+=usize::from(passed);
                    rows.push(json!({"id":id,"input":input,"status":"checked","check_body":{"valid":true,"unchecked":0},"gensolid":properties,"planar_tessellation":mesh,
                        "independent_expected_volume":expected,"known_answer":{"absolute_error":error,"relative_error":error/expected.abs().max(precision.linear().powi(3)),
                            "precision_budget":budget,"passed":passed},"occt":{"status":"pending"}}));
                    work.push(Work { row:rows.len()-1,boundary,properties,edge_length,diameter });
                }
                Err(error)=>rows.push(json!({"id":id,"input":input,"status":"gensolid_error","message":error.to_string(),"occt":{"status":"not_generated"}})),
            }
        }
    }
    let total = FRAMES * 3;
    let mut compared = 0;
    let mut agreed = 0;
    let oracle = match invoke(&work, precision) {
        Ok(batch) if batch.available => {
            let mut observed = BTreeMap::new();
            for item in batch.cases {
                if observed.insert(item.id.clone(), item).is_some() {
                    return Err(Error::parse("M2 OCCT batch repeats a case identifier"));
                }
            }
            for case in &work {
                let row = rows
                    .get_mut(case.row)
                    .ok_or_else(|| Error::internal("M2 case row is missing"))?;
                let observation = if let Some(item) = observed.remove(&case.boundary.id) {
                    if item.valid && item.closed && item.status == "ok" {
                        if let Some(properties) = item.properties {
                            let measured = comparison(
                                &case.properties,
                                &properties,
                                case.edge_length,
                                case.diameter,
                                precision,
                            );
                            compared += 1;
                            agreed += usize::from(measured["passed"].as_bool() == Some(true));
                            json!({"status":"compared","valid":true,"closed":true,"properties":properties,"comparison":measured})
                        } else {
                            json!({"status":"malformed_observation","message":"valid OCCT result lacks properties"})
                        }
                    } else {
                        json!({"status":item.status,"valid":item.valid,"closed":item.closed,"message":item.message})
                    }
                } else {
                    json!({"status":"missing_observation"})
                };
                row["occt"] = observation;
            }
            if !observed.is_empty() {
                return Err(Error::parse("M2 OCCT returned unexpected case identifiers"));
            }
            json!({"available":true,"status":batch.status,"backend":batch.backend})
        }
        Ok(batch) => {
            for case in &work {
                if let Some(row) = rows.get_mut(case.row) {
                    row["occt"] = json!({"status":"oracle_unavailable","message":batch.message});
                }
            }
            json!({"available":false,"status":batch.status,"message":batch.message})
        }
        Err(error) => {
            let status = if error.code() == ErrorCode::Unsupported {
                "oracle_unavailable"
            } else {
                "oracle_error"
            };
            for case in &work {
                if let Some(row) = rows.get_mut(case.row) {
                    row["occt"] = json!({"status":status,"message":error.to_string()});
                }
            }
            json!({"available":false,"status":status,"message":error.to_string()})
        }
    };
    Ok(
        json!({"scope":"convex polygon prism, linear sweep and homothetic planar loft; not general M2 or L0 acceptance",
        "seed":SEED,"frames":FRAMES,"cases":total,"checked_bodies":work.len(),"gensolid_errors":total-work.len(),
        "known_answer_volume":{"source":"independent input-profile shoelace area times normal height; frustum factor (1+s+s^2)/3",
            "passed":known_answers_passed,"inputs":total,"pass_rate":known_answers_passed as f64/total as f64},
        "planar_tessellation":{"inputs":total,"passed":planar_mesh_passed,"pass_rate":planar_mesh_passed as f64/total as f64,
            "chord_tolerance":planar_chord.value()},
        "occt_mass":{"oracle":oracle,"scope":"OCCT reconstructs the same world vertices and outward planar B-rep loops; independent integration and native validity, not independent modelling intent",
            "use_triangulation":false,"inertia_reference":"centroid, world-coordinate axes, unit density","inputs":total,"compared":compared,"agreed":agreed,
            "comparison_coverage":compared as f64/total as f64,"agreement":if compared>0 {Some(agreed as f64/compared as f64)} else {None}},
        "rows":rows,"analytic_patch_tessellation":mesh_samples(precision)?}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_m2_inputs_are_checked_and_keep_the_full_denominator() -> Result<()> {
        let report = run(&Precision::DEFAULT)?;
        assert_eq!(report["cases"], 48);
        assert_eq!(report["checked_bodies"], 48);
        assert_eq!(report["known_answer_volume"]["passed"], 48);
        assert_eq!(report["rows"].as_array().map(Vec::len), Some(48));
        assert_eq!(report["analytic_patch_tessellation"]["passed"], 3);
        assert_eq!(report["planar_tessellation"]["inputs"], 48);
        assert_eq!(report["planar_tessellation"]["passed"], 48);
        if report["occt_mass"]["oracle"]["available"].as_bool() == Some(true) {
            assert_eq!(report["occt_mass"]["compared"], 48);
            assert_eq!(report["occt_mass"]["agreed"], 48);
        }
        Ok(())
    }
}
