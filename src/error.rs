//! Error types for the VCF filter library.

use thiserror::Error;

/// Errors that can occur during VCF filtering operations.
#[derive(Error, Debug)]
pub enum VcfFilterError {
    /// Failed to parse the VCF header.
    #[error("Header parse error: {0}")]
    HeaderParseError(String),

    /// Failed to parse a VCF data row.
    #[error("Row parse error: {0}")]
    RowParseError(String),

    /// Failed to parse a filter expression.
    #[error("Filter parse error: {0}")]
    FilterParseError(String),

    /// Error during filter evaluation.
    #[error("Evaluation error: {0}")]
    EvaluationError(String),

    /// Attempted to access an unknown field.
    #[error("Unknown field: {0}")]
    UnknownField(String),

    /// Invalid array index access.
    #[error("Invalid index {index} for field {field} (length {length})")]
    InvalidIndex {
        field: String,
        index: usize,
        length: usize,
    },

    /// Type mismatch during comparison.
    #[error("Type mismatch: cannot compare {left} with {right}")]
    TypeMismatch { left: String, right: String },
}

/// Result type alias for VCF filter operations.
pub type Result<T> = std::result::Result<T, VcfFilterError>;
