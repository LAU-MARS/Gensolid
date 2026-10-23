//! Fixed-orientation translation sweeps along arbitrary positive-weight NURBS paths.
//!
//! This original tensor-product construction uses the identity
//! `S(u,v) = C(u) + P(v) - P(v_min)`. It preserves the two input parameterizations
//! and never fits, tessellates, or implicitly rotates the profile along the path.

use crate::{BCurve, BSurface};
use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_math::{Interval, Point3};

fn normalized_weights(curve: &BCurve) -> Result<Vec<f64>> {
    let Some(weights) = curve.weights() else {
        return Ok(vec![1.0; curve.points().len()]);
    };
    let scale = weights.iter().copied().fold(0.0_f64, f64::max);
    let normalized: Vec<_> = weights.iter().map(|w| w / scale).collect();
    if normalized.iter().any(|&w| !w.is_normal() || w <= 0.0) {
        return Err(Error::degenerate(
            "translation sweep weight range is not representable",
        ));
    }
    Ok(normalized)
}

fn coordinates(point: Point3) -> [f64; 3] {
    [point.x, point.y, point.z]
}

fn check_coordinate_resolution(points: &[Point3], p: &Precision) -> Result<()> {
    // Positive NURBS bases confine all evaluations to the control hull. Check
    // its largest coordinate spacing, not only rounding at the control points:
    // exactly representable controls can still have unresolved interior values.
    let mut spacing = [0.0_f64; 3];
    for point in points {
        for (k, coordinate) in coordinates(*point).into_iter().enumerate() {
            let scale = coordinate.abs();
            spacing[k] = spacing[k].max(scale.next_up() - scale);
        }
    }
    let resolution = spacing[0].hypot(spacing[1]).hypot(spacing[2]);
    if !resolution.is_finite() || resolution > p.solver() {
        return Err(Error::unimplemented(
            "translation sweep coordinates cannot resolve session precision across the control hull",
        ));
    }
    Ok(())
}

fn check_extent(curve: &BCurve, p: &Precision) -> Result<()> {
    let origin = curve.points()[0];
    let mut extent = 0.0_f64;
    for point in curve.points() {
        let length = point.distance(origin);
        if !length.is_finite() {
            return Err(Error::invalid_argument(
                "translation sweep displacement is not representable",
            ));
        }
        extent = extent.max(length);
    }
    // Positive-weight NURBS lie inside their control hull. A hull contained in
    // a session-resolution ball cannot provide a resolved profile/path extent.
    if extent <= p.degenerate_length() {
        return Err(Error::degenerate(
            "translation sweep profile or path is below session precision",
        ));
    }
    Ok(())
}

/// Exact tensor-product translation sweep with fixed profile orientation.
///
/// The profile is anchored at the actual natural-domain start of the path,
/// including for unclamped input. Relative coordinates avoid subtracting two
/// large evaluated world positions. Degrees, knots and U/V domains are retained;
/// output weights are products of independently normalized positive weights.
/// Subnormal weights/products and control hulls whose coordinate spacing exceeds
/// the session solver budget are refused rather than losing significant geometry.
///
/// This constructs geometry only. It does not certify a regular, nonintersecting
/// sheet or solid: the modelling caller must run its complete graph/body checks.
/// In particular, profiles are not transported by a moving frame or path tangent.
pub fn sweep_translation(profile: &BCurve, path: &BCurve, p: &Precision) -> Result<BSurface> {
    check_extent(profile, p)?;
    check_extent(path, p)?;
    check_coordinate_resolution(profile.points(), p)?;
    let profile_weights = normalized_weights(profile)?;
    let path_weights = normalized_weights(path)?;
    let reference = path.points()[0];
    let relative = BCurve::new(
        path.degree(),
        path.knots().to_vec(),
        path.points()
            .iter()
            .map(|&q| {
                let d = q - reference;
                Point3::new(d.x, d.y, d.z)
            })
            .collect(),
        Some(path_weights.clone()),
    )?;
    let anchor = relative.eval(relative.range().lo).point;
    if !anchor.is_finite() {
        return Err(Error::degenerate(
            "translation sweep path anchor is not representable",
        ));
    }
    let count = profile
        .points()
        .len()
        .checked_mul(path.points().len())
        .ok_or_else(|| {
            Error::new(
                ErrorCode::CapacityExceeded,
                "translation sweep control net is too large",
            )
        })?;
    let mut points = Vec::new();
    let mut weights = Vec::new();
    points.try_reserve_exact(count).map_err(|_| {
        Error::new(
            ErrorCode::CapacityExceeded,
            "translation sweep control net allocation failed",
        )
    })?;
    weights.try_reserve_exact(count).map_err(|_| {
        Error::new(
            ErrorCode::CapacityExceeded,
            "translation sweep weight allocation failed",
        )
    })?;
    for (&control, &weight) in profile.points().iter().zip(&profile_weights) {
        for (&path_control, &path_weight) in relative.points().iter().zip(&path_weights) {
            let product = weight * path_weight;
            if !product.is_normal() || product <= 0.0 {
                return Err(Error::degenerate(
                    "translation sweep weight product is not representable",
                ));
            }
            let point = control + (path_control - anchor);
            if !point.is_finite() {
                return Err(Error::invalid_argument(
                    "translation sweep control point is not representable",
                ));
            }
            let mut error = [0.0_f64; 3];
            for (k, component) in error.iter_mut().enumerate() {
                let exact = Interval::point(coordinates(control)[k])
                    + (Interval::point(coordinates(path_control)[k])
                        - Interval::point(coordinates(anchor)[k]));
                let residual = Interval::point(coordinates(point)[k]) - exact;
                *component = residual.lo().abs().max(residual.hi().abs());
            }
            let error = error[0].hypot(error[1]).hypot(error[2]);
            if !error.is_finite() || error > p.solver() {
                return Err(Error::unimplemented(
                    "translation sweep control arithmetic exceeds the session solver budget",
                ));
            }
            points.push(point);
            weights.push(product);
        }
    }
    check_coordinate_resolution(&points, p)?;
    BSurface::new(
        profile.degree(),
        path.degree(),
        profile.knots().to_vec(),
        path.knots().to_vec(),
        profile.points().len(),
        path.points().len(),
        points,
        Some(weights),
    )
}
