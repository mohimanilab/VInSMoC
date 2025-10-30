use std::io;
use thiserror::Error;

use molecule::{adducts::ParseAdductError, BuildMolError};

/// A unified error type for dereplication.
#[derive(Error, Debug)]
pub enum Error {
    /// The chemical database had an invalid SMILES.
    #[error(transparent)]
    SmilesParse(#[from] BuildMolError),
    /// We attempted to use an invalid [`molecule::Adduct`].
    #[error(transparent)]
    AdductParse(#[from] ParseAdductError),
    /// Attempted to construct a [`crate::MsmSmiles`] from a [`crate::Msm`] that didn't contain a
    /// SMILES.
    #[error("tried to create SMILES match without SMILES data")]
    MissingSmiles,
    /// Some sort of IO error occurred.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Parsing the database CSV failed, or writing the output TSV failed.
    #[error(transparent)]
    Csv(#[from] csv::Error),
    /// Spectral parsing failed.
    #[error("couldn't parse spectra")]
    Spectrum(#[from] spectrum::Error),
    /// Missing molecule index in chemical database.
    #[error("no molecule with index {index} in database \"{file}\"")]
    MissingMolIdx {
        /// Attempted molecule index.
        index: usize,
        /// Database CSV file missing the index.
        file: String,
    },
}
