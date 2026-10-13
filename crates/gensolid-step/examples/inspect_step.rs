//! Prints an honest validity report for a STEP physical file.

use gensolid_core::Precision;
use gensolid_step::validity_report;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: inspect_step <path.step>")?;
    let text = std::fs::read_to_string(path)?;
    let report = validity_report(&text, Precision::DEFAULT);
    println!(
        "status={:?} decoded={} checked={} faults={} unchecked={}",
        report.status,
        report.decoded_bodies,
        report.checked_bodies,
        report.faults,
        report.unchecked
    );
    for message in report.messages {
        println!("{message}");
    }
    Ok(())
}
