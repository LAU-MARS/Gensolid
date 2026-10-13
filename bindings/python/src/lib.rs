//! Python bindings. Method names match the PK-style Rust API in `gensolid-api`.

use gensolid_api::{Axis2, Error, ErrorCode, Session as ApiSession, SessionOptions, Tag};
use pyo3::exceptions::{
    PyKeyError, PyNotImplementedError, PyRuntimeError, PyTypeError, PyValueError,
};
use pyo3::prelude::*;

type PyPoint = (f64, f64, f64);

fn to_py(error: Error) -> PyErr {
    let message = error.to_string();
    match error.code() {
        ErrorCode::NotFound => PyKeyError::new_err(message),
        ErrorCode::WrongEntityClass => PyTypeError::new_err(message),
        ErrorCode::Unimplemented => PyNotImplementedError::new_err(message),
        ErrorCode::InvalidArgument | ErrorCode::Degenerate | ErrorCode::OutOfRange => {
            PyValueError::new_err(message)
        }
        _ => PyRuntimeError::new_err(message),
    }
}

/// A modelling session.
#[pyclass(module = "gensolid")]
struct Session {
    inner: ApiSession,
}

#[pymethods]
impl Session {
    #[new]
    fn new() -> Self {
        Self {
            inner: ApiSession::start(SessionOptions::default()),
        }
    }

    /// Creates a solid block centred on `location`; returns the body tag.
    #[pyo3(signature = (x, y, z, location=[0.0, 0.0, 0.0], axis=[0.0, 0.0, 1.0], ref_direction=[1.0, 0.0, 0.0]))]
    fn body_create_solid_block(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        location: [f64; 3],
        axis: [f64; 3],
        ref_direction: [f64; 3],
    ) -> PyResult<Tag> {
        let basis = Axis2 {
            location,
            axis,
            ref_direction,
        };
        self.inner
            .body_create_solid_block(x, y, z, Some(&basis))
            .map_err(to_py)
    }

    /// Face tags of a body.
    fn body_ask_faces(&mut self, body: Tag) -> PyResult<Vec<Tag>> {
        self.inner.body_ask_faces(body).map_err(to_py)
    }

    /// Edge tags of a body.
    fn body_ask_edges(&mut self, body: Tag) -> PyResult<Vec<Tag>> {
        self.inner.body_ask_edges(body).map_err(to_py)
    }

    /// Vertex tags of a body.
    fn body_ask_vertices(&mut self, body: Tag) -> PyResult<Vec<Tag>> {
        self.inner.body_ask_vertices(body).map_err(to_py)
    }

    /// Faults as `(state, entity_tag, message)`; empty for a valid body.
    fn body_check(&mut self, body: Tag) -> PyResult<Vec<(String, Tag, String)>> {
        let faults = self.inner.body_check(body).map_err(to_py)?;
        Ok(faults
            .into_iter()
            .map(|f| (format!("{:?}", f.state), f.entity, f.message))
            .collect())
    }

    /// Axis-aligned bounding box as `((min_x, min_y, min_z), (max_x, max_y, max_z))`.
    fn body_ask_bounding_box(&self, body: Tag) -> PyResult<(PyPoint, PyPoint)> {
        let ([x0, y0, z0], [x1, y1, z1]) = self.inner.body_ask_bounding_box(body).map_err(to_py)?;
        Ok(((x0, y0, z0), (x1, y1, z1)))
    }

    /// Volume of the solid regions of a body.
    fn body_ask_volume(&self, body: Tag) -> PyResult<f64> {
        self.inner.body_ask_volume(body).map_err(to_py)
    }

    /// Position of a vertex.
    fn vertex_ask_point(&self, vertex: Tag) -> PyResult<(f64, f64, f64)> {
        let [x, y, z] = self.inner.vertex_ask_point(vertex).map_err(to_py)?;
        Ok((x, y, z))
    }
}

/// Gensolid geometric modelling kernel.
#[pymodule]
fn gensolid(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Session>()?;
    Ok(())
}
