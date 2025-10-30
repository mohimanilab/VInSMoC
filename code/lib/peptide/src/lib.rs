//! This crate handles peptidic molecules. It includes functionality to parse peptides from
//! strings, convert peptidic sequences to general molecules, and ways to represent generalized
//! amino acids.

#![warn(missing_docs)]

mod amino_acid;
pub use amino_acid::{
    backbone, AminoAcid, LargeCluster, NonStdAminoAcid, ParseAminoAcidError,
    ParseStdAminoAcidError, SmallCluster, StdAminoAcid, StdOrNonStdAa, NUM_STANDARD_AMINO_ACIDS,
};
mod peptide;
pub use crate::peptide::Peptide;
mod trans;
pub use trans::{TranslatedCodon, TranslatedPeptide};
pub mod parsers;
