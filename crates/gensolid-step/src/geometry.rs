//! Independent STEP geometry mapping from the published EXPRESS entity definitions.
//! References: https://www.steptools.com/docs/stp_aim/html/t_edge_curve.html and
//! https://www.steptools.com/docs/stp_aim/html/t_b_spline_surface_with_knots.html .

use std::collections::{HashMap, HashSet};

use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{
    BCurve, BSurface, Circle, Cone, Curve, Cylinder, Ellipse, Line, Plane, Sphere, Surface, Torus,
};

use crate::parse::{Instance, Value, name_skipped, nums, ref_of, refs};
use crate::{cartesian, direction, expect_type, instance, placement};

fn arg(args: &[Value], index: usize) -> Result<&Value> {
    args.get(index)
        .ok_or_else(|| Error::parse(format!("missing STEP geometry argument {index}")))
}

fn number(value: &Value) -> Result<f64> {
    match value {
        Value::Num(n) if n.is_finite() => Ok(*n),
        _ => Err(Error::parse("expected a finite STEP geometry number")),
    }
}

fn integer(value: &Value) -> Result<usize> {
    let n = number(value)?;
    if n < 0.0 || n.fract() != 0.0 || n >= usize::MAX as f64 {
        return Err(Error::parse("expected a non-negative STEP integer"));
    }
    Ok(n as usize)
}

fn list(value: &Value) -> Result<&[Value]> {
    match value {
        Value::List(items) => Ok(items),
        _ => Err(Error::parse("expected a STEP geometry list")),
    }
}

fn logical(value: &Value) -> Result<()> {
    match value {
        Value::Bool(_) => Ok(()),
        Value::Enum(name) if name == "U" => Ok(()),
        _ => Err(Error::parse("expected a STEP logical (.T., .F. or .U.)")),
    }
}

fn enumeration(value: &Value, allowed: &[&str]) -> Result<()> {
    match value {
        Value::Enum(name) if allowed.contains(&name.as_str()) => Ok(()),
        _ => Err(Error::parse("invalid STEP geometry enumeration")),
    }
}

fn fields(args: &[Value], count: usize, name: &str) -> Result<()> {
    if args.len() != count {
        return Err(Error::parse(format!(
            "STEP {name} has {} fields, expected {count}",
            args.len()
        )));
    }
    Ok(())
}

fn knots(multiplicities: &Value, values: &Value, degree: usize, count: usize) -> Result<Vec<f64>> {
    if degree >= count {
        return Err(Error::parse(
            "B-spline degree exceeds its control point count",
        ));
    }
    let multiplicities = list(multiplicities)?;
    let values = nums(values)?;
    if multiplicities.len() != values.len() || values.windows(2).any(|w| w[0] >= w[1]) {
        return Err(Error::parse("invalid STEP B-spline knot table"));
    }
    let expected = count
        .checked_add(degree)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| Error::parse("B-spline knot count overflows"))?;
    let mut result = Vec::with_capacity(expected);
    for (multiplicity, knot) in multiplicities.iter().zip(values) {
        let n = integer(multiplicity)?;
        if n == 0 || n > degree + 1 || n > expected - result.len() {
            return Err(Error::parse("invalid STEP knot multiplicity"));
        }
        result.extend(std::iter::repeat_n(knot, n));
    }
    if result.len() != expected {
        return Err(Error::parse("inconsistent STEP knot count"));
    }
    Ok(result)
}

pub(crate) fn read_curve(
    instances: &HashMap<u32, Instance>,
    id: u32,
    precision: &Precision,
) -> Result<Curve> {
    let mut current = id;
    let mut seen = HashSet::new();
    loop {
        if !seen.insert(current) {
            return Err(Error::parse(format!(
                "cyclic STEP curve reference #{current}"
            )));
        }
        let entity = instance(instances, current)?;
        if entity.components.contains_key("B_SPLINE_CURVE")
            || entity.type_name == "B_SPLINE_CURVE_WITH_KNOTS"
        {
            return read_bcurve(instances, entity).map(Curve::BCurve);
        }
        if !matches!(
            entity.type_name.as_str(),
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" | "LINE" | "CIRCLE" | "ELLIPSE"
        ) {
            return Err(Error::unimplemented(format!(
                "STEP curve #{current} {} is not implemented",
                entity.type_name
            )));
        }
        let args = name_skipped(&entity.args)?;
        match entity.type_name.as_str() {
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" => {
                fields(args, 3, &entity.type_name)?;
                let associated = refs(arg(args, 1)?)?;
                if !(1..=2).contains(&associated.len()) {
                    return Err(Error::parse(
                        "STEP SURFACE_CURVE needs one or two associated geometries",
                    ));
                }
                enumeration(arg(args, 2)?, &["CURVE_3D", "PCURVE_S1", "PCURVE_S2"])?;
                // The existing model stores 3D carrier curves. PCURVE trimming is not
                // substituted with a mesh or an invented spatial approximation.
                current = ref_of(arg(args, 0)?)?;
            }
            "LINE" => {
                fields(args, 2, "LINE")?;
                let origin = cartesian(instances, ref_of(arg(args, 0)?)?)?;
                let vector = expect_type(instances, ref_of(arg(args, 1)?)?, "VECTOR")?;
                let d = direction(instances, ref_of(arg(vector, 0)?)?)?;
                let magnitude = number(arg(vector, 1)?)?;
                if magnitude <= 0.0 {
                    return Err(Error::parse("STEP VECTOR magnitude must be positive"));
                }
                return Line::new(origin, d * magnitude, precision).map(Curve::Line);
            }
            "CIRCLE" => {
                fields(args, 2, "CIRCLE")?;
                return Circle::new(
                    placement(instances, ref_of(arg(args, 0)?)?, precision)?,
                    number(arg(args, 1)?)?,
                    precision,
                )
                .map(Curve::Circle);
            }
            "ELLIPSE" => {
                fields(args, 3, "ELLIPSE")?;
                return Ellipse::new(
                    placement(instances, ref_of(arg(args, 0)?)?, precision)?,
                    number(arg(args, 1)?)?,
                    number(arg(args, 2)?)?,
                    precision,
                )
                .map(Curve::Ellipse);
            }
            _ => {
                return Err(Error::unimplemented(format!(
                    "STEP curve #{current} {} is not implemented",
                    entity.type_name
                )));
            }
        }
    }
}

pub(crate) fn read_surface(
    instances: &HashMap<u32, Instance>,
    id: u32,
    precision: &Precision,
) -> Result<Surface> {
    let entity = instance(instances, id)?;
    if entity.components.contains_key("B_SPLINE_SURFACE")
        || entity.type_name == "B_SPLINE_SURFACE_WITH_KNOTS"
    {
        return read_bsurface(instances, entity).map(Surface::BSurface);
    }
    if !matches!(
        entity.type_name.as_str(),
        "PLANE"
            | "CYLINDRICAL_SURFACE"
            | "CONICAL_SURFACE"
            | "SPHERICAL_SURFACE"
            | "TOROIDAL_SURFACE"
    ) {
        return Err(Error::unimplemented(format!(
            "STEP surface #{id} {} is not implemented",
            entity.type_name
        )));
    }
    let args = name_skipped(&entity.args)?;
    fields(
        args,
        match entity.type_name.as_str() {
            "PLANE" => 1,
            "CONICAL_SURFACE" | "TOROIDAL_SURFACE" => 3,
            _ => 2,
        },
        &entity.type_name,
    )?;
    let frame = placement(instances, ref_of(arg(args, 0)?)?, precision)?;
    match entity.type_name.as_str() {
        "PLANE" => Ok(Surface::Plane(Plane::new(frame))),
        "CYLINDRICAL_SURFACE" => {
            Cylinder::new(frame, number(arg(args, 1)?)?, precision).map(Surface::Cylinder)
        }
        "CONICAL_SURFACE" => Cone::new(
            frame,
            number(arg(args, 1)?)?,
            number(arg(args, 2)?)?,
            precision,
        )
        .map(Surface::Cone),
        "SPHERICAL_SURFACE" => {
            Sphere::new(frame, number(arg(args, 1)?)?, precision).map(Surface::Sphere)
        }
        "TOROIDAL_SURFACE" => Torus::new(
            frame,
            number(arg(args, 1)?)?,
            number(arg(args, 2)?)?,
            precision,
        )
        .map(Surface::Torus),
        _ => Err(Error::internal("unhandled supported STEP surface")),
    }
}

fn component<'a>(entity: &'a Instance, name: &str) -> Result<&'a [Value]> {
    entity
        .components
        .get(name)
        .map(Vec::as_slice)
        .ok_or_else(|| Error::unimplemented(format!("STEP complex entity requires {name}")))
}

fn read_bcurve(instances: &HashMap<u32, Instance>, entity: &Instance) -> Result<BCurve> {
    let (base, knot_args) = if entity.components.is_empty() {
        let args = name_skipped(&entity.args)?;
        fields(args, 8, "B_SPLINE_CURVE_WITH_KNOTS")?;
        (&args[..5], &args[5..])
    } else {
        (
            component(entity, "B_SPLINE_CURVE")?,
            component(entity, "B_SPLINE_CURVE_WITH_KNOTS")?,
        )
    };
    fields(base, 5, "B_SPLINE_CURVE")?;
    fields(knot_args, 3, "B_SPLINE_CURVE_WITH_KNOTS component")?;
    enumeration(
        arg(base, 2)?,
        &[
            "POLYLINE_FORM",
            "CIRCULAR_ARC",
            "ELLIPTIC_ARC",
            "PARABOLIC_ARC",
            "HYPERBOLIC_ARC",
            "UNSPECIFIED",
        ],
    )?;
    logical(arg(base, 3)?)?;
    logical(arg(base, 4)?)?;
    enumeration(
        arg(knot_args, 2)?,
        &[
            "UNIFORM_KNOTS",
            "QUASI_UNIFORM_KNOTS",
            "PIECEWISE_BEZIER_KNOTS",
            "UNSPECIFIED",
        ],
    )?;
    let degree = integer(arg(base, 0)?)?;
    let points = refs(arg(base, 1)?)?
        .into_iter()
        .map(|id| cartesian(instances, id))
        .collect::<Result<Vec<_>>>()?;
    let knots = knots(arg(knot_args, 0)?, arg(knot_args, 1)?, degree, points.len())?;
    let weights = entity
        .components
        .get("RATIONAL_B_SPLINE_CURVE")
        .map(|args| {
            fields(args, 1, "RATIONAL_B_SPLINE_CURVE")?;
            nums(arg(args, 0)?)
        })
        .transpose()?;
    BCurve::new(degree, knots, points, weights)
}

fn read_bsurface(instances: &HashMap<u32, Instance>, entity: &Instance) -> Result<BSurface> {
    let (base, knot_args) = if entity.components.is_empty() {
        let args = name_skipped(&entity.args)?;
        fields(args, 12, "B_SPLINE_SURFACE_WITH_KNOTS")?;
        (&args[..7], &args[7..])
    } else {
        (
            component(entity, "B_SPLINE_SURFACE")?,
            component(entity, "B_SPLINE_SURFACE_WITH_KNOTS")?,
        )
    };
    fields(base, 7, "B_SPLINE_SURFACE")?;
    fields(knot_args, 5, "B_SPLINE_SURFACE_WITH_KNOTS component")?;
    enumeration(
        arg(base, 3)?,
        &[
            "PLANE_SURF",
            "CYLINDRICAL_SURF",
            "CONICAL_SURF",
            "SPHERICAL_SURF",
            "TOROIDAL_SURF",
            "SURF_OF_REVOLUTION",
            "RULED_SURF",
            "GENERALISED_CONE",
            "QUADRIC_SURF",
            "SURF_OF_LINEAR_EXTRUSION",
            "UNSPECIFIED",
        ],
    )?;
    for value in &base[4..] {
        logical(value)?;
    }
    enumeration(
        arg(knot_args, 4)?,
        &[
            "UNIFORM_KNOTS",
            "QUASI_UNIFORM_KNOTS",
            "PIECEWISE_BEZIER_KNOTS",
            "UNSPECIFIED",
        ],
    )?;
    let u_degree = integer(arg(base, 0)?)?;
    let v_degree = integer(arg(base, 1)?)?;
    let rows = list(arg(base, 2)?)?;
    let u_count = rows.len();
    let mut v_count = None;
    let mut points = Vec::new();
    for row in rows {
        let row = refs(row)?;
        if let Some(count) = v_count
            && count != row.len()
        {
            return Err(Error::parse("ragged B-spline control net"));
        }
        v_count = Some(row.len());
        for id in row {
            points.push(cartesian(instances, id)?);
        }
    }
    let v_count = v_count.ok_or_else(|| Error::parse("empty B-spline control net"))?;
    let u_knots = knots(arg(knot_args, 0)?, arg(knot_args, 2)?, u_degree, u_count)?;
    let v_knots = knots(arg(knot_args, 1)?, arg(knot_args, 3)?, v_degree, v_count)?;
    let weights = if let Some(args) = entity.components.get("RATIONAL_B_SPLINE_SURFACE") {
        fields(args, 1, "RATIONAL_B_SPLINE_SURFACE")?;
        let rows = list(arg(args, 0)?)?;
        if rows.len() != u_count {
            return Err(Error::parse(
                "B-spline weight row count differs from control net",
            ));
        }
        let mut values = Vec::new();
        for row in rows {
            let row = nums(row)?;
            if row.len() != v_count {
                return Err(Error::parse(
                    "B-spline weight column count differs from control net",
                ));
            }
            values.extend(row);
        }
        Some(values)
    } else {
        None
    };
    BSurface::new(
        u_degree, v_degree, u_knots, v_knots, u_count, v_count, points, weights,
    )
}
