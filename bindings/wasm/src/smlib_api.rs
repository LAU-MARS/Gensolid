//! Browser convenience API. Rust callers can use the complete `smlib` re-export.
use super::to_js;
use gensolid::smlib::{Backend, Curve, Surface};
use gensolid_core::Precision;
use gensolid_math::Point3;
use wasm_bindgen::prelude::*;

fn point(values: &[f64]) -> Result<Point3, JsError> {
    let [x, y, z]: [f64; 3] = values
        .try_into()
        .map_err(|_| JsError::new("expected three coordinates"))?;
    Ok(Point3::new(x, y, z))
}
/// A curve/surface arena backed by the packaged SMLib Wasm module.
#[wasm_bindgen]
pub struct SmlibBackend {
    inner: Backend,
}
#[wasm_bindgen]
impl SmlibBackend {
    /// Creates an arena after `initializeSmlib` has completed.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Result<Self, JsError> {
        Ok(Self {
            inner: Backend::new(Precision::DEFAULT).map_err(to_js)?,
        })
    }
    /// Creates an XY circle.
    pub fn circle(&self, center: &[f64], radius: f64) -> Result<SmlibCurve, JsError> {
        Ok(SmlibCurve {
            inner: self.inner.circle(point(center)?, radius).map_err(to_js)?,
        })
    }
    /// Creates a straight curve between two points.
    #[wasm_bindgen(js_name = lineSegment)]
    pub fn line_segment(&self, start: &[f64], end: &[f64]) -> Result<SmlibCurve, JsError> {
        Ok(SmlibCurve {
            inner: self
                .inner
                .line_segment(point(start)?, point(end)?)
                .map_err(to_js)?,
        })
    }
    /// Creates a four-corner surface. Input contains four consecutive XYZ triples.
    #[wasm_bindgen(js_name = cornerSurface)]
    pub fn corner_surface(&self, corners: &[f64]) -> Result<SmlibSurface, JsError> {
        if corners.len() != 12 {
            return Err(JsError::new("expected four XYZ corners"));
        }
        Ok(SmlibSurface {
            inner: self
                .inner
                .corner_surface([
                    point(&corners[0..3])?,
                    point(&corners[3..6])?,
                    point(&corners[6..9])?,
                    point(&corners[9..12])?,
                ])
                .map_err(to_js)?,
        })
    }
}
/// A curve handle. Its session stays alive until the last related handle is freed.
#[wasm_bindgen]
pub struct SmlibCurve {
    inner: Curve,
}
#[wasm_bindgen]
impl SmlibCurve {
    /// Natural parameter interval.
    pub fn domain(&self) -> Result<Vec<f64>, JsError> {
        Ok(self.inner.domain().map_err(to_js)?.to_vec())
    }
    /// Returns position, first derivative and second derivative as nine numbers.
    pub fn evaluate(&self, t: f64) -> Result<Vec<f64>, JsError> {
        let e = self.inner.evaluate(t).map_err(to_js)?;
        Ok(vec![
            e.point.x, e.point.y, e.point.z, e.d1.x, e.d1.y, e.d1.z, e.d2.x, e.d2.y, e.d2.z,
        ])
    }
    /// Returns a reversed copy.
    pub fn reversed(&self) -> Result<SmlibCurve, JsError> {
        Ok(Self {
            inner: self.inner.reversed().map_err(to_js)?,
        })
    }
}
/// A surface handle sharing ownership of its SMLib session.
#[wasm_bindgen]
pub struct SmlibSurface {
    inner: Surface,
}
#[wasm_bindgen]
impl SmlibSurface {
    /// Returns `[u_min, u_max, v_min, v_max]`.
    pub fn domain(&self) -> Result<Vec<f64>, JsError> {
        Ok(self.inner.domain().map_err(to_js)?.to_vec())
    }
    /// Returns XYZ triples for position, du, dv, duu, duv, dvv.
    pub fn evaluate(&self, u: f64, v: f64) -> Result<Vec<f64>, JsError> {
        let e = self.inner.evaluate([u, v]).map_err(to_js)?;
        Ok(vec![
            e.point.x, e.point.y, e.point.z, e.du.x, e.du.y, e.du.z, e.dv.x, e.dv.y, e.dv.z,
            e.duu.x, e.duu.y, e.duu.z, e.duv.x, e.duv.y, e.duv.z, e.dvv.x, e.dvv.y, e.dvv.z,
        ])
    }
}
