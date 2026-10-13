//! Parasolid XT (transmit file) support.
//!
//! XT is a first-class format for Gensolid: the data model is designed so that every XT
//! node maps onto a kernel entity without approximation. The format is published by
//! Siemens in the *Parasolid XT Format Reference* for the purpose of writing translators.
//!
//! Status:
//! - [`header`]: text header (`**PART1;` … `**END_OF_HEADER`) and the transmit preamble
//!   (modeller version, schema key). Read and write.
//! - [`tokens`]: tokenizer for the text node stream.
//! - [`import`]: schema-aware text and neutral-binary decoding, exact analytic and open
//!   NURBS carriers, topology mapping, and mandatory kernel validity checking.
//! - [`report`]: distinguishes checked, invalid, malformed, and unsupported files.
//!
//! The public ABC smoke sample exercises the exact Onshape 30.1 embedded schema
//! `SCH_3001168_30100_13006`. Other built-in profiles remain reviewed subsets, and
//! wire/acorn, tolerant fin, periodic spline and compound geometry may be unsupported.
//! Passing the public smoke sample does not establish private benchmark coverage.

pub mod header;
pub mod import;
pub mod report;
mod schema;
pub mod tokens;

pub use header::{TextFile, TextHeader, read_text, write_text_header};
pub use import::{XtImport, read_bytes, read_bytes_with_schema};
/// Exact schema-provider interfaces for callers with an external XT schema catalog.
pub use parasolid_core::{
    InMemorySchemaProvider, SchemaCatalogLimits, SchemaProvider, parse_schema_catalog,
};
pub use report::{
    XtStatus, XtValidityReport, validity_report, validity_report_bytes,
    validity_report_bytes_with_precision,
};
pub use tokens::{Token, TokenKind, tokenize};
