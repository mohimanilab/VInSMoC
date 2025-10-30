//! Iterators for molecule-spectrum match candidates.

use std::{borrow::Cow, collections::BTreeMap, iter::Enumerate, marker::PhantomData};

use chemical_db::{CandidateFilter, ChemicalDb, DatabaseRow, LazyDb, LazyRows, SingletonDb};
use ordered_float::OrderedFloat;
use spectrum::SpectrumCollection;

use crate::Error;

use super::{MatchedSpectra, Row};

/// An iterator over candidate molecule-spectrum matches, grouped together by molecule.
///
/// See [`InMemoryCandidates`] and [`LazyCandidates`] for details on iterators. Note that this
/// yields [`Result`]s for compatibility with [`LazyCandidates`].
pub enum Candidates<'db, 'adduct> {
    /// Candidates built from a [`ChemicalDb`]
    InMemory(InMemoryCandidates<'db, 'adduct>),
    /// Candidates built from a [`SingletonDb`]
    Singleton(SingletonCandidates<'db, 'adduct>),
    /// Candidates built from a [`LazyDb`]
    Lazy(LazyCandidates<'db, 'adduct>),
}

impl<'db, 'adduct> Candidates<'db, 'adduct> {
    /// Construct a new candidate iterator from an in-memory [`ChemicalDb`].
    pub fn in_memory(
        db: &'db ChemicalDb,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Self {
        Self::InMemory(InMemoryCandidates::new(db, spectra, filter))
    }

    /// Construct a new candidate iterator from an in-memory [`SingletonDb`].
    pub fn singleton(
        db: &'db SingletonDb,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Self {
        Self::Singleton(SingletonCandidates::new(db, spectra, filter))
    }

    /// Construct a new candidate iterator from a [`LazyDb`].
    ///
    /// # Errors
    /// Propagates errors from [`LazyCandidates::new`].
    pub fn lazy(
        db: &'db LazyDb,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Result<Self, Error> {
        Ok(Self::Lazy(LazyCandidates::new(db, spectra, filter)?))
    }
}

/// An iterator over candidate molecule-spectrum matches for in-memory [`ChemicalDb`]s, grouped
/// together by molecule.
///
/// This iterator yields [`Row`] and [`MatchedSpectra`] tuples. The [`Row`] characterizes the
/// molecule being matched and the [`MatchedSpectra`] can be used to retrieve all spectra that
/// could explain the mass of the molecule according to a [`CandidateFilter`]. Any [`Row`]s for
/// which there would be no explained spectra are skipped, so typically this iterator will yield
/// fewer values than there are rows in the database. It is guaranteed that all spectra will be
/// found for a molecule when it's first yielded. The [`Cow`] in [`Row`] is guaranteed to always be
/// [`Cow::Borrowed`] with a `'db` lifetime.
pub struct InMemoryCandidates<'db, 'adduct> {
    /// The internal iterator wrapped by this iterator.
    ///
    /// This provides the database row indices and the actual [`DatabaseRow`]s.
    db: Enumerate<<&'db ChemicalDb as IntoIterator>::IntoIter>,
    /// The candidate filter that stores allowed precursor adducts and precursor error tolerance.
    filter: &'adduct CandidateFilter,
    /// The precursors of the spectral dataset converted to masses using
    /// [`InMemoryCandidates::filter`] and grouped by mass.
    ///
    /// The entries for a given mass key are all precursor adduct and scan numbers that can explain
    /// that mass exactly.
    masses: BTreeMap<OrderedFloat<f64>, MatchedSpectra<'adduct>>,
}

impl<'db, 'adduct> InMemoryCandidates<'db, 'adduct> {
    /// Construct a new [`InMemoryCandidates`] iterator for the provided `db` and `spectra`.
    ///
    /// The `filter` is used to determine what precursor adducts are used and convert precursor
    /// m/zs from `spectra` into masses.
    pub fn new(
        db: &'db ChemicalDb,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Self {
        Self {
            db: db.into_iter().enumerate(),
            filter,
            masses: spectra.precursors(&filter.adducts),
        }
    }

    /// Compute the candidate list for the provided `row`.
    fn candidates(&self, row: &DatabaseRow) -> MatchedSpectra<'adduct> {
        let (lower, upper) = self.filter.error_tolerance.bounds(row.exact_mass());
        let lower = OrderedFloat::<f64>::from(lower);
        let upper = OrderedFloat::<f64>::from(upper);
        self.masses
            .range(lower..=upper)
            .flat_map(|(_, v)| v.iter().copied())
            .collect::<Vec<_>>()
    }
}

impl<'db, 'adduct> Iterator for InMemoryCandidates<'db, 'adduct> {
    type Item = (Row<'db>, MatchedSpectra<'adduct>);

    fn next(&mut self) -> Option<Self::Item> {
        let (mut row_idx, mut row) = self.db.next()?;
        let mut candidates = self.candidates(row);

        // this could be done cleaner recursively but it'll be slower
        while candidates.is_empty() {
            let entry = self.db.next()?;
            row_idx = entry.0;
            row = entry.1;
            candidates = self.candidates(row);
        }

        Some(((row_idx, Cow::Borrowed(row)), candidates))
    }
}

/// An iterator over candidate molecule-spectrum matches for [`LazyDb`]s, grouped by molecule.
///
/// This iterator yields [`Row`] and [`MatchedSpectra`] tuples. The [`Row`] characterizes the
/// molecule being matched and the [`MatchedSpectra`] can be used to retrieve all spectra that
/// could explain the mass of the molecule according to a [`CandidateFilter`]. Any [`Row`]s for
/// which there would be no explained spectra are skipped, so typically this iterator will yield
/// fewer values than there are rows in the database. It is guaranteed that all spectra will be
/// found for a molecule when it's first yielded. The [`Cow`] in [`Row`] is guaranteed to always be
/// [`Cow::Owned`].
///
/// This is to [`LazyDb`] as [`InMemoryCandidates`] is to [`ChemicalDb`].
///
/// # Errors
/// The process of lazily reading rows from disks is fallible due to I/O, so this iterator yields
/// [`Result`]s.
pub struct LazyCandidates<'db, 'adduct> {
    /// The internal iterator wrapped by this one.
    ///
    /// This provides the database row indices and the actual [`DatabaseRow`]s.
    db: LazyRows<std::fs::File>,
    /// The candidate filter that stores allowed precursor adducts and precursor error tolerance.
    filter: &'adduct CandidateFilter,
    /// The precursors of the spectral dataset converted to masses using
    /// [`InMemoryCandidates::filter`] and grouped by mass.
    ///
    /// The entries for a given mass key are all precursor adduct and scan numbers that can explain
    /// that mass exactly.
    masses: BTreeMap<OrderedFloat<f64>, MatchedSpectra<'adduct>>,
    /// A placeholder to force the database lifetime into this type.
    marker: PhantomData<&'db ()>,
}

impl<'db, 'adduct> LazyCandidates<'db, 'adduct> {
    /// Construct a new [`LazyCandidates`] iterator for the provided `db` and `spectra`.
    ///
    /// The `filter` is used to determine what precursor adducts are used and convert precursor
    /// m/zs from `spectra` into masses.
    ///
    /// # Errors
    /// This will propagate errors from calls to [`LazyDb::rows`].
    pub fn new(
        db: &'db LazyDb,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Result<Self, Error> {
        Ok(Self {
            db: db.rows()?,
            filter,
            masses: spectra.precursors(&filter.adducts),
            marker: PhantomData,
        })
    }

    /// Compute the candidate list for the provided `row`.
    fn candidates(&self, row: &DatabaseRow) -> MatchedSpectra<'adduct> {
        let (lower, upper) = self.filter.error_tolerance.bounds(row.exact_mass());
        let lower = OrderedFloat::<f64>::from(lower);
        let upper = OrderedFloat::<f64>::from(upper);
        self.masses
            .range(lower..upper)
            .flat_map(|(_, v)| v.iter().copied())
            .collect::<Vec<_>>()
    }
}

impl<'db, 'adduct> Iterator for LazyCandidates<'db, 'adduct> {
    type Item = Result<(Row<'db>, MatchedSpectra<'adduct>), Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let (mut row_idx, row) = self.db.next()?;
        if let Err(e) = row {
            return Some(Err(e.into()));
        }

        let mut row = row.unwrap();
        let mut candidates = self.candidates(&row);
        while candidates.is_empty() {
            let entry = self.db.next()?;
            row_idx = entry.0;
            if let Err(e) = entry.1 {
                return Some(Err(e.into()));
            }
            row = entry.1.unwrap();
            candidates = self.candidates(&row);
        }

        Some(Ok(((row_idx, Cow::Owned(row)), candidates)))
    }
}

/// An iterator over candidate molecule-spectrum matches for [`SingletonDb`]s.
pub struct SingletonCandidates<'db, 'adduct> {
    db: &'db SingletonDb,
    /// The candidate filter that stores allowed precursor adducts and precursor error tolerance.
    filter: &'adduct CandidateFilter,
    /// The precursors of the spectral dataset converted to masses using
    /// [`InMemoryCandidates::filter`] and grouped by mass.
    ///
    /// The entries for a given mass key are all precursor adduct and scan numbers that can explain
    /// that mass exactly.
    masses: BTreeMap<OrderedFloat<f64>, MatchedSpectra<'adduct>>,
    consumed: bool,
}

impl<'db, 'adduct> SingletonCandidates<'db, 'adduct> {
    /// Construct a new [`SingletonCandidates`] iterator for the provided `db` and `spectra`.
    ///
    /// The `filter` is used to determine what precursor adducts are used and convert precursor
    /// m/zs from `spectra` into masses.
    pub fn new(
        db: &'db SingletonDb,
        spectra: &SpectrumCollection,
        filter: &'adduct CandidateFilter,
    ) -> Self {
        Self {
            db,
            filter,
            masses: spectra.precursors(&filter.adducts),
            consumed: false,
        }
    }

    /// Compute the candidate list for the single row in this database.
    pub fn candidates(&self) -> MatchedSpectra<'adduct> {
        let row = &self.db.entry;
        let (lower, upper) = self.filter.error_tolerance.bounds(row.exact_mass());
        let lower = OrderedFloat::<f64>::from(lower);
        let upper = OrderedFloat::<f64>::from(upper);
        self.masses
            .range(lower..=upper)
            .flat_map(|(_, v)| v.iter().copied())
            .collect::<Vec<_>>()
    }
}

impl<'db, 'adduct> Iterator for SingletonCandidates<'db, 'adduct> {
    type Item = (Row<'db>, MatchedSpectra<'adduct>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.consumed {
            return None;
        }

        let (row_idx, row) = (0, self.db.entry.clone());
        let candidates = self.candidates();
        self.consumed = true;
        Some(((row_idx, Cow::Owned(row)), candidates))
    }
}

impl<'db, 'adduct> Iterator for Candidates<'db, 'adduct> {
    type Item = Result<((usize, Cow<'db, DatabaseRow>), MatchedSpectra<'adduct>), Error>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Candidates::InMemory(candidates) => candidates.next().map(Ok),
            Candidates::Lazy(candidates) => candidates.next(),
            Candidates::Singleton(candidates) => candidates.next().map(Ok),
        }
    }
}
