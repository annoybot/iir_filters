//! Errors which originate in this library.
use crate::errors::Error::IllegalArgument;
use std::num::ParseFloatError;
use thiserror::Error as ThisError;

#[derive(ThisError, Debug, PartialEq)]
/// A custom error type for errors which can occur in this library.
pub enum Error {
    /// Errors caused the user providing an illegal argument to a function.
    #[error("Illegal Argument: {0}")]
    IllegalArgument(String),

    /// Errors resulting from logic errors in the code. Not recoverable by the user.
    #[error("Internal error: {0}")]
    InternalError(String),
}

impl From<ParseFloatError> for Error {
    fn from(err: ParseFloatError) -> Self {
        IllegalArgument(err.to_string())
    }
}

