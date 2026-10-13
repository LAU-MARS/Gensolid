//! Browser transport. Each call copies spans between the two Wasm memories.
use super::{Output, native_error};
use gensolid_core::{Error, Precision, Result};
use js_sys::{Array, Float64Array, Uint32Array};
use wasm_bindgen::prelude::*;

#[wasm_bindgen(module = "/js/smlib-bridge.js")]
extern "C" {
    #[wasm_bindgen(catch, js_name = initialize)]
    fn initialize_js(module: JsValue) -> std::result::Result<(), JsValue>;
    #[wasm_bindgen(catch, js_name = create)]
    fn create_js(linear: f64, angular: f64, size: f64) -> std::result::Result<Array, JsValue>;
    #[wasm_bindgen(js_name = destroy)]
    pub(super) fn destroy(session: usize);
    #[wasm_bindgen(catch, js_name = call)]
    fn call_js(
        session: usize,
        op: u32,
        ids: &[u32],
        values: &[f64],
        sizes: &[u32],
    ) -> std::result::Result<Array, JsValue>;
}
fn js_error(value: JsValue) -> Error {
    Error::internal(format!("SMLib browser transport: {value:?}"))
}
/// Installs an initialized Emscripten SMLib module before creating a backend.
/// A module cannot be replaced while sessions exist.
pub fn initialize(module: JsValue) -> Result<()> {
    initialize_js(module).map_err(js_error)
}
pub(super) fn create(precision: Precision) -> Result<usize> {
    let result =
        create_js(precision.linear(), precision.angular(), precision.size()).map_err(js_error)?;
    let code = result
        .get(0)
        .as_f64()
        .ok_or_else(|| Error::internal("invalid create status"))? as i32;
    if code != 0 {
        return Err(native_error(code, 0, 0));
    }
    Ok(result
        .get(1)
        .as_f64()
        .ok_or_else(|| Error::internal("invalid session"))? as usize)
}
pub(super) fn call(
    session: usize,
    op: u32,
    ids: &[usize],
    values: &[f64],
    sizes: &[usize],
) -> Result<Output> {
    let ids: Vec<u32> = ids.iter().map(|&v| v as u32).collect();
    let sizes: Vec<u32> = sizes.iter().map(|&v| v as u32).collect();
    let result = call_js(session, op, &ids, values, &sizes).map_err(js_error)?;
    let code = result
        .get(0)
        .as_f64()
        .ok_or_else(|| Error::internal("invalid call status"))? as i32;
    let upstream = result
        .get(1)
        .as_string()
        .ok_or_else(|| Error::internal("invalid upstream status"))?
        .parse::<i64>()
        .map_err(|_| Error::internal("invalid upstream status"))?;
    if code != 0 {
        return Err(native_error(code, upstream, op));
    }
    let values = Float64Array::new(&result.get(2)).to_vec();
    let indices = Uint32Array::new(&result.get(3))
        .to_vec()
        .into_iter()
        .map(|v| v as usize)
        .collect();
    if values.iter().any(|v| !v.is_finite()) {
        return Err(Error::internal("SMLib returned non-finite output"));
    }
    Ok(Output { values, indices })
}
