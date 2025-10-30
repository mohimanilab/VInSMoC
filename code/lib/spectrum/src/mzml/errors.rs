use quick_xml::de as xml;
use std::num::{ParseFloatError, ParseIntError};

use thiserror::Error;

use crate::xml::DecodeDataError;

/// An error that can occur when parsing an mzML file.
#[derive(Error, Debug)]
pub enum Error {
    /// An error occurred while using quick-xml to deserialize the contents
    #[error("unable to deserialize XML")]
    Deserialize(#[from] xml::DeError),
    /// Missing start or close mzML tag
    #[error("no <mzML> tag in file")]
    NoMzmlTag,
    /// The index provided in the file doesn't match its position in the spectrum list
    ///
    /// The first index is the provided one, the second is its position in the spectrum list.
    #[error("spectrum index {0} doesn't match position in spectrum list ({1})")]
    MismatchedIndex(usize, usize),
    /// Multiple floating point precision cvParams were provided for a single spectrum
    #[error("duplicate or conflicting float types")]
    MultipleFloat,
    /// Multiple compression type cvParams were provided for a single spectrum
    #[error("duplicate or conflicting compression specifications")]
    MultipleCompression,
    /// Multiple of the same data array were provided for a single spectrum
    #[error("duplicate or conflicting array data types")]
    MultipleArrayDatatype,
    /// A spectrum did not contain the required floating point precision cvParam
    #[error("missing required float type")]
    NoFloat,
    /// A spectrum did not contain the required compression type cvParam
    #[error("missing required compression specification")]
    NoCompression,
    /// A spectrum didn't have an m/z data array
    #[error("spectrum didn't contain an m/z array")]
    NoMzArray,
    /// A spectrum didn't have an intensity data array
    #[error("spectrum didn't contain an intensity array")]
    NoIntensityArray,
    /// A spectrum didn't have a precursor ion element
    #[error("spectrum didn't have a precursor ion")]
    MissingPrecursor,
    /// A spectrum's precursor ion units were not in m/z
    #[error("precursor ion was not in m/z units")]
    InvalidPrecursorUnits,
    /// The retention time units were not in either minutes or seconds
    #[error("invalid retention time units")]
    InvalidRetentionTimeUnits,
    /// A cvParam that requires a value didn't contain one.
    ///
    /// Contained value is the name of the cvParam.
    #[error("cvParam {0} missing required value attribute")]
    NoCvParamValue(String),
    /// Failed to parse floating point value
    #[error("couldn't parse <cvParam> value as floating point")]
    FloatParse(#[from] ParseFloatError),
    /// Failed to parse integer value
    #[error("couldn't parse <cvParam> value as integer")]
    IntParse(#[from] ParseIntError),
    /// Received invalid floating point precision
    #[error("invalid float bit-precision: {0}; expecting 32 or 64")]
    BadFloatPrecision(usize),
    /// An error occurred when decoding binary data array
    #[error("unable to decode data array")]
    BinaryDataArray(#[from] DecodeDataError),
}
