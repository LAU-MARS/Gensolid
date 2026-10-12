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
//! - Schema-driven node decoding, the binary (`.x_b`) encoding, and model writing are not
//!   implemented yet.
//!
//! Only synthetic fixtures written from the format description are tested so far. Real
//! `.x_t` files from Parasolid-based applications must be added to the verifier corpus
//! before this crate is relied on.

pub mod header;
pub mod tokens;

pub use header::{TextFile, TextHeader, read_text, write_text_header};
pub use tokens::{Token, TokenKind, tokenize};
