//! WebAssembly bindings. Method names are camel-cased versions of the PK-style Rust API in
//! `gensolid-api`.

use gensolid_api::{Axis2, Error, Session as ApiSession, SessionOptions, Tag};
use wasm_bindgen::prelude::*;

fn to_js(error: Error) -> JsError {
    JsError::new(&error.to_string())
}

/// A modelling session.
#[wasm_bindgen]
pub struct Session {
    inner: ApiSession,
}

#[wasm_bindgen]
impl Session {
    /// Starts a session with default precision.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            inner: ApiSession::start(SessionOptions::default()),
        }
    }

    /// Creates a solid block centred on the origin; returns the body tag.
    #[wasm_bindgen(js_name = bodyCreateSolidBlock)]
    pub fn body_create_solid_block(&mut self, x: f64, y: f64, z: f64) -> Result<Tag, JsError> {
        self.inner
            .body_create_solid_block(x, y, z, None)
            .map_err(to_js)
    }

    /// Creates a solid block in the given coordinate system; returns the body tag.
    #[wasm_bindgen(js_name = bodyCreateSolidBlockAt)]
    pub fn body_create_solid_block_at(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        location: &[f64],
        axis: &[f64],
        ref_direction: &[f64],
    ) -> Result<Tag, JsError> {
        let triple = |v: &[f64]| -> Result<[f64; 3], JsError> {
            v.try_into()
                .map_err(|_| JsError::new("expected an array of 3 numbers"))
        };
        let basis = Axis2 {
            location: triple(location)?,
            axis: triple(axis)?,
            ref_direction: triple(ref_direction)?,
        };
        self.inner
            .body_create_solid_block(x, y, z, Some(&basis))
            .map_err(to_js)
    }

    /// Face tags of a body.
    #[wasm_bindgen(js_name = bodyAskFaces)]
    pub fn body_ask_faces(&mut self, body: Tag) -> Result<Vec<Tag>, JsError> {
        self.inner.body_ask_faces(body).map_err(to_js)
    }

    /// Edge tags of a body.
    #[wasm_bindgen(js_name = bodyAskEdges)]
    pub fn body_ask_edges(&mut self, body: Tag) -> Result<Vec<Tag>, JsError> {
        self.inner.body_ask_edges(body).map_err(to_js)
    }

    /// Vertex tags of a body.
    #[wasm_bindgen(js_name = bodyAskVertices)]
    pub fn body_ask_vertices(&mut self, body: Tag) -> Result<Vec<Tag>, JsError> {
        self.inner.body_ask_vertices(body).map_err(to_js)
    }

    /// Number of faults found by the checker; 0 for a valid body.
    #[wasm_bindgen(js_name = bodyCheck)]
    pub fn body_check(&mut self, body: Tag) -> Result<usize, JsError> {
        Ok(self.inner.body_check(body).map_err(to_js)?.len())
    }

    /// Position of a vertex as `[x, y, z]`.
    #[wasm_bindgen(js_name = vertexAskPoint)]
    pub fn vertex_ask_point(&self, vertex: Tag) -> Result<Vec<f64>, JsError> {
        Ok(self.inner.vertex_ask_point(vertex).map_err(to_js)?.to_vec())
    }
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}
