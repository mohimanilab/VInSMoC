//! Abstractions over chemical databases for code de-duplication.
//!
//! The [`chemical_db`] crate exposes three kinds of chemical database:
//! 1. A completely in-memory [`chemical_db::ChemicalDb`]
//! 2. A pre-processed [`chemical_db::pre::PreprocessedDb`]
//! 3. A lazily-loaded [`chemical_db::LazyDb`].
//!
//! Both 2 and 3 have similar access patterns (you can just get a lazy database from a
//! pre-processed one), but the in-memory database is a separate entity. This module exposes traits
//! for generalizing data access from all kinds of chemical database. It also provides iterators
//! for molecule-spectrum match candidates.
//!
//! # [`AsCandidates`]
//! This trait is implemented for types that can be converted into an association between a
//! compound in the database any one or more spectra from a spectral collection. All previously
//! mentioned databases have this feature, but may require different mechanisms to do so. The key
//! entrypoint here is [`AsCandidates::candidates`], which creates a [`Candidates`] iterator over
//! database indices, rows, and all matched scan numbers with their precursor adducts. Since for
//! lazy databases the operation of reading a row is fallible all [`Candidates`] iterator variants
//! return results.
//!
//! This trait is intended to only be implemented for **references**. The `'db` lifetime should be
//! used as the lifetime for the implementing reference.
//!
//! ## Conversion to candidate matches
//! The [`AsCandidates::candidates`] function builds an iterator over [`Row`] and
//! [`MatchedSpectra`] tuples. The [`Row`] stores the database index of a molecule and its original
//! database row, as either a borrowed or owned value. The [`MatchedSpectra`] is a vector of all
//! scan numbers that can explain that molecule's mass. It also includes the precursor adduct that
//! must be used to enable that explanation.

use std::borrow::Cow;

use chemical_db::{CandidateFilter, ChemicalDb, DatabaseRow, LazyDb, SingletonDb};
use molecule::Adduct;
use spectrum::SpectrumCollection;

use crate::{matches::MsmStub, Error};

pub mod iter;
use iter::Candidates;

/// A single row in a chemical database.
///
/// This is a tuple of the 0-based index of the row in the database and the [`DatabaseRow`] itself.
pub type Row<'db> = (usize, Cow<'db, DatabaseRow>);

/// All spectra candidates for a particular molecule or mass.
///
/// This is a list of tuples of adducts and scan numbers. The adduct marks the precursor adduct
/// that makes the scan number match the molecule or mass.
pub type MatchedSpectra<'adduct> = Vec<(&'adduct Adduct, usize)>;

/// A trait for database references that can used to find molecule-spectrum match candidates.
///
/// The constraints on implementers is largely to enforce that this is only implemented on database
/// references. The `'db` lifetime is intended to be the lifetime of the reference this is
/// implemented for.
pub trait AsCandidates<'db>: Sized + Sync + Copy {
    /// Construct candidates for this database.
    ///
    /// It's expected that this trait is only implemented on references. Given a database
    /// reference, the `spectra` we're searching, and a `filter` on candidates this builds the full
    /// iterator over candidates.
    ///
    /// # Errors
    /// In case constructing the iterator is fallible this returns a [`Result`].
    fn candidates<'adduct>(
        self,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Result<Candidates<'db, 'adduct>, Error>;

    /// Construct a [`MsmStub`] for the provided `row`.
    fn stub<'a>(self, row: &'a Row<'db>, spectrum_file: &'a str) -> MsmStub<'a> {
        MsmStub {
            spectrum_file,
            mol_idx: row.0,
            mol_name: row.1.name(),
            smiles: &row.1.smiles,
            inchikey14: row.1.inchikey14(),
            mol_mass: row.1.exact_mass(),
        }
    }
}

impl<'db> AsCandidates<'db> for &'db ChemicalDb {
    fn candidates<'adduct>(
        self,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Result<Candidates<'db, 'adduct>, Error> {
        Ok(Candidates::in_memory(self, spectra, filter))
    }
}

impl<'db> AsCandidates<'db> for &'db SingletonDb {
    fn candidates<'adduct>(
        self,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Result<Candidates<'db, 'adduct>, Error> {
        Ok(Candidates::singleton(self, spectra, filter))
    }
}

impl<'db> AsCandidates<'db> for &'db LazyDb {
    fn candidates<'adduct>(
        self,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Result<Candidates<'db, 'adduct>, Error> {
        Candidates::lazy(self, spectra, filter)
    }
}
