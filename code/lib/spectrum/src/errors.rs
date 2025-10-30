use std::io;

use thiserror::Error;

/// An error that can occur when parsing a spectrum file.
///
/// This provides a unified endpoint combining I/O errors and format-specific parsing errors. See
/// the individual format errors for more information about them.
#[derive(Error, Debug)]
pub enum Error {
    /// Failed to parse an MGF file.
    #[error("failed to parse MGF file")]
    Mgf,
    /// Failed to parse an mzXML file.
    #[error("failed to parse mzXML file")]
    MzXml(#[from] crate::mzxml::Error),
    /// Failed to parse an mzML file.
    #[error("failed to parse mzML file")]
    MzMl(#[from] crate::mzml::Error),
    /// An IO error occurred before parsing began.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// A path given to [`SpectrumCollection::from_path`](crate::SpectrumCollection::from_path)
    /// does not have a known extension.
    ///
    /// Valid extensions are determined by
    /// 1. Retrieving the extension with [`std::path::Path::extension`]
    /// 2. Converting it to lowercase
    /// 3. Checking if it is one of "mgf", "mzxml", or "mzml"
    #[error("unknown spectral format")]
    UnknownFormat,
    /// Missing scan in spectrum.
    #[error("no scan {scan} in spectrum collection")]
    MissingScan {
        /// Attempted scan number
        scan: usize,
    },
}
