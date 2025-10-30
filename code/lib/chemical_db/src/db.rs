use std::{
    convert::TryFrom,
    fmt::Debug,
    io::{BufReader, Read, Write},
    iter::FromIterator,
    ops::Range,
    path::Path,
};

use paths::FilePath;
use rayon::{
    iter::{
        FromParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelExtend,
        ParallelIterator,
    },
    slice::{ParallelSlice, ParallelSliceMut},
};
use serde::{Deserialize, Serialize};

use crate::DatabaseRow;

/// A whole chemical database.
///
/// This consists of [`DatabaseRow`] for every molecule in underlying chemical database. These
/// databases are canonically stored as CSVs with columns `smiles`, `exact_mass`, `name`, and
/// `inchikey14`. The `inchikey14` column is the first 14 characters of the molecule's InChIKey,
/// which can be used to compare for exactly 2D molecular equality.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ChemicalDb {
    pub(crate) entries: Vec<DatabaseRow>,
}

/// A chemical database with a single entry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SingletonDb {
    /// The molecule entry in the singleton DB.
    pub entry: DatabaseRow,
}

/// A wrapper around a file path to a database CSV intended for [`clap`] use.
#[derive(Debug, Clone, Default, clap::Parser)]
#[clap(help_heading = "Database")]
pub struct DbCsv {
    /// Path to chemical compound database CSV.
    ///
    /// This is expected to be a CSV file with four columns. These columns can appear in any order
    /// but must be named "name", "smiles", "exact_mass", and "inchikey14". The "inchikey14" column
    /// is the first 14 characters of the compound's InChIKey.
    #[clap(name = "DATABASE", short = 'd', long = "database")]
    pub path: FilePath,
}

impl TryFrom<DbCsv> for ChemicalDb {
    type Error = csv::Error;

    fn try_from(csv: DbCsv) -> Result<Self, Self::Error> {
        Self::from_path(csv.path)
    }
}

impl TryFrom<&DbCsv> for ChemicalDb {
    type Error = csv::Error;

    fn try_from(csv: &DbCsv) -> Result<Self, Self::Error> {
        Self::from_path(&csv.path)
    }
}

impl ChemicalDb {
    /// Constructs a chemical database from anything that implements [`Read`].
    ///
    /// Constructed chemical databases are guaranteed to be sorted by ascending exact mass.
    /// The reader is expected to be a CSV as described in [the struct docs](ChemicalDb). If the
    /// input is malformed this function returns the underlying CSV error.
    ///
    /// See also [`ChemicalDb::from_reader_unchecked`].
    ///
    /// # Errors
    /// Returns an error if the reader is an invalid CSV.
    pub fn from_reader<R: Read>(reader: R) -> Result<Self, csv::Error> {
        let mut db = Self::from_reader_unchecked(reader)?;
        db.sort_rows();
        Ok(db)
    }

    /// Constructs a chemical database from anything that implements [`Read`].
    ///
    /// This is identical to [`ChemicalDb::from_reader`] but it sorts the database by mass in
    /// parallel.
    ///
    /// # Errors
    /// Returns an error if the reader is an invalid CSV.
    pub fn par_from_reader<R: Read>(reader: R) -> Result<Self, csv::Error> {
        let mut db = Self::from_reader_unchecked(reader)?;
        db.par_sort_rows();
        Ok(db)
    }

    /// Constructs a chemical database from anything that implements [`Read`].
    ///
    /// Does *not* guarantee the database is sorted by ascending exact mass. This will be faster
    /// than using [`ChemicalDb::from_reader`] but will be incorrect unless the provided CSV
    /// actually is already sorted by ascending exact mass.
    ///
    /// See also [`ChemicalDb::from_reader`].
    ///
    /// # Errors
    /// Returns an error if the reader is an invalid CSV.
    pub fn from_reader_unchecked<R: Read>(reader: R) -> Result<Self, csv::Error> {
        let mut reader = csv::Reader::from_reader(BufReader::new(reader));
        let entries = reader
            .deserialize::<DatabaseRow>()
            .collect::<Result<_, _>>()?;
        Ok(Self { entries })
    }

    /// Constructs a chemical database from a path on the filesystem.
    ///
    /// Constructed chemical databases are guaranteed to be sorted by ascending exact mass.
    /// The path is expected to be a CSV as described in [the struct docs](ChemicalDb). If the
    /// input is malformed this function returns the underlying CSV error.
    ///
    /// See also [`ChemicalDb::from_path_unchecked`].
    ///
    /// # Errors
    /// Returns an error if the underlying path is an invalid CSV.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, csv::Error> {
        let mut db = Self::from_path_unchecked(path)?;
        db.sort_rows();
        Ok(db)
    }

    /// Constructs a chemical database from a path on the filesystem.
    ///
    /// This is identical to [`ChemicalDb::from_path`] but it sorts the database by mass in
    /// parallel.
    ///
    /// # Errors
    /// Returns an error if the underlying path is an invalid CSV.
    pub fn par_from_path(path: impl AsRef<Path>) -> Result<Self, csv::Error> {
        let mut db = Self::from_path_unchecked(path)?;
        db.par_sort_rows();
        Ok(db)
    }

    /// Constructs a chemical database from a path on the filesystem.
    ///
    /// Does *not* guarantee the database is sorted by ascending exact mass. This will be faster
    /// than using [`ChemicalDb::from_path`] but will be incorrect unless the provided CSV
    /// actually is already sorted by ascending exact mass.
    ///
    /// See also [`ChemicalDb::from_path`].
    ///
    /// # Errors
    /// Returns an error if the underlying path is an invalid CSV.
    pub fn from_path_unchecked(path: impl AsRef<Path>) -> Result<Self, csv::Error> {
        let mut reader = csv::Reader::from_path(path)?;
        let entries = reader
            .deserialize::<DatabaseRow>()
            .collect::<Result<_, _>>()?;
        Ok(Self { entries })
    }

    /// Writes this database as as CSV to the provided `writer`.
    ///
    /// # Note
    /// The [`csv`] crate automatically buffers all writers so the provided writer should not be
    /// buffered.
    ///
    /// # Errors
    /// Propagates errors from [`csv`] serialization.
    pub fn write_csv<W: Write>(&self, writer: W) -> csv::Result<()> {
        let mut writer = csv::WriterBuilder::new().from_writer(writer);
        for row in &self.entries {
            writer.serialize(row)?;
        }

        Ok(())
    }

    /// Sorts database by ascending mass.
    ///
    /// This typically will not need to be called, as both [`Self::from_reader`] and
    /// [`Self::from_path`] already ensure that the database is sorted.
    ///
    /// See also [`slice::sort_unstable_by`].
    pub fn sort_rows(&mut self) {
        self.entries
            .sort_unstable_by(|row1, row2| row1.exact_mass.partial_cmp(&row2.exact_mass).unwrap());
    }

    /// Sorts database by ascending mass in parallel.
    ///
    /// This typically will not need to be called, as both [`Self::from_reader`] and
    /// [`Self::from_path`] already ensure that the database is sorted.
    ///
    /// See also [`rayon::slice::ParallelSliceMut::par_sort_unstable_by`]
    pub fn par_sort_rows(&mut self) {
        self.entries.par_sort_unstable_by(|row1, row2| {
            row1.exact_mass.partial_cmp(&row2.exact_mass).unwrap()
        });
    }

    /// Searches database for molecules with exact mass with provided error tolerance of provided
    /// mass.
    ///
    /// If there are no molecules that match these parameters this returns [`None`].
    /// Otherwise, it will give the indices in the sorted database that match the input
    /// parameters.
    // TODO: test binary search with duplicates
    pub fn query(&self, mass: f64, error_tolerance: f64) -> Option<Range<usize>> {
        let lower_bound = mass - error_tolerance;
        let upper_bound = mass + error_tolerance;

        let ind_low = match self
            .entries
            .binary_search_by(|mol| mol.exact_mass.partial_cmp(&lower_bound).unwrap())
        {
            Ok(i) | Err(i) => i,
        };
        if ind_low == self.len() {
            return None;
        }
        let mut ind_high = ind_low;
        while ind_high < self.len() && self.entries[ind_high].exact_mass <= upper_bound {
            ind_high += 1;
        }
        if ind_high == ind_low {
            return None;
        }
        Some(ind_low..ind_high)
    }

    /// Returns the number of molecules in the database.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Checks if the database has no molecules.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Accesses a molecule in the database associated with the given index.
    ///
    /// Returns [`None`] if the given index is out of bounds.
    pub fn get(&self, index: usize) -> Option<&DatabaseRow> {
        self.entries.get(index)
    }

    /// Constructs an iterator over rows of the database
    pub fn iter(&self) -> std::slice::Iter<'_, DatabaseRow> {
        self.entries.iter()
    }
}

impl AsRef<ChemicalDb> for ChemicalDb {
    fn as_ref(&self) -> &ChemicalDb {
        self
    }
}

impl<'a> IntoIterator for &'a ChemicalDb {
    type Item = &'a DatabaseRow;
    type IntoIter = std::slice::Iter<'a, DatabaseRow>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl IntoIterator for ChemicalDb {
    type Item = DatabaseRow;
    type IntoIter = std::vec::IntoIter<DatabaseRow>;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter()
    }
}

impl FromIterator<DatabaseRow> for ChemicalDb {
    /// Construct a chemical database from an iterator over database rows.
    fn from_iter<T: IntoIterator<Item = DatabaseRow>>(iter: T) -> Self {
        Self {
            entries: iter.into_iter().collect(),
        }
    }
}

impl Extend<DatabaseRow> for ChemicalDb {
    fn extend<T: IntoIterator<Item = DatabaseRow>>(&mut self, iter: T) {
        self.entries.extend(iter);
    }
}

impl IntoParallelIterator for ChemicalDb {
    type Iter = <Vec<DatabaseRow> as IntoParallelIterator>::Iter;

    type Item = DatabaseRow;

    fn into_par_iter(self) -> Self::Iter {
        self.entries.into_par_iter()
    }
}

impl<'data> IntoParallelRefIterator<'data> for ChemicalDb {
    type Iter = <Vec<DatabaseRow> as IntoParallelRefIterator<'data>>::Iter;

    type Item = <Vec<DatabaseRow> as IntoParallelRefIterator<'data>>::Item;

    fn par_iter(&'data self) -> Self::Iter {
        self.entries.par_iter()
    }
}

impl FromParallelIterator<DatabaseRow> for ChemicalDb {
    /// Construct a chemical database from a parallel iterator over database rows.
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = DatabaseRow>,
    {
        Self {
            entries: par_iter.into_par_iter().collect(),
        }
    }
}

impl ParallelExtend<DatabaseRow> for ChemicalDb {
    fn par_extend<I>(&mut self, par_iter: I)
    where
        I: IntoParallelIterator<Item = DatabaseRow>,
    {
        self.entries.par_extend(par_iter);
    }
}

impl ParallelSlice<DatabaseRow> for ChemicalDb {
    fn as_parallel_slice(&self) -> &[DatabaseRow] {
        self.entries.as_slice()
    }
}

impl ParallelSliceMut<DatabaseRow> for ChemicalDb {
    fn as_parallel_slice_mut(&mut self) -> &mut [DatabaseRow] {
        self.entries.as_mut_slice()
    }
}
