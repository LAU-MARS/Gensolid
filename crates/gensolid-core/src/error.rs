//! Kernel-wide error type.

use std::fmt;

/// Stable error categories. The numeric values are part of the C API.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(i32)]
#[non_exhaustive]
pub enum ErrorCode {
    /// An argument is out of its documented domain.
    InvalidArgument = 1,
    /// A referenced entity does not exist (never created, or deleted).
    NotFound = 2,
    /// An entity exists but is of the wrong class for the operation.
    WrongEntityClass = 3,
    /// Geometry is degenerate at session precision (zero length, parallel axes, ...).
    Degenerate = 4,
    /// A parameter lies outside the valid range of a curve or surface.
    OutOfRange = 5,
    /// The input is valid but outside what the kernel supports by design.
    Unsupported = 6,
    /// The input is valid and in scope, but the algorithm is not written yet.
    Unimplemented = 7,
    /// An internal container limit was reached.
    CapacityExceeded = 8,
    /// A file or stream could not be parsed.
    Parse = 9,
    /// An internal invariant was violated; always a kernel bug.
    Internal = 10,
}

/// Error returned by kernel operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Error {
    code: ErrorCode,
    message: String,
}

impl Error {
    /// Creates an error with a category and a human-readable message.
    pub fn new(code: ErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    /// Shorthand for [`ErrorCode::InvalidArgument`].
    pub fn invalid_argument(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::InvalidArgument, message)
    }

    /// Shorthand for [`ErrorCode::NotFound`].
    pub fn not_found(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::NotFound, message)
    }

    /// Shorthand for [`ErrorCode::Degenerate`].
    pub fn degenerate(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Degenerate, message)
    }

    /// Shorthand for [`ErrorCode::Unsupported`].
    pub fn unsupported(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Unsupported, message)
    }

    /// Shorthand for [`ErrorCode::Unimplemented`].
    pub fn unimplemented(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Unimplemented, message)
    }

    /// Shorthand for [`ErrorCode::Parse`].
    pub fn parse(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Parse, message)
    }

    /// Shorthand for [`ErrorCode::Internal`].
    pub fn internal(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Internal, message)
    }

    /// Error category.
    pub const fn code(&self) -> ErrorCode {
        self.code
    }

    /// Human-readable description.
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.code, self.message)
    }
}

impl std::error::Error for Error {}

/// Result alias used throughout the kernel.
pub type Result<T, E = Error> = std::result::Result<T, E>;
