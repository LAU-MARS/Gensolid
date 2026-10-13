//! Print a strict validity report for XT files in their original units.
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for argument in std::env::args_os().skip(1) {
        let path = PathBuf::from(argument);
        let bytes = std::fs::read(&path)?;
        let report = gensolid_xt::validity_report_bytes(&bytes);
        println!(
            "{}: {:?}; schema={:?}; bodies={}; {}",
            path.display(),
            report.status,
            report.schema,
            report.decoded_bodies,
            report.message
        );
        if let Ok(imported) = gensolid_xt::read_bytes(&bytes, gensolid_core::Precision::DEFAULT) {
            for body in imported.bodies {
                println!(
                    "  faces={}, edges={}, vertices={}",
                    imported.model.body_faces(body)?.len(),
                    imported.model.body_edges(body)?.len(),
                    imported.model.body_vertices(body)?.len()
                );
            }
        }
    }
    Ok(())
}
