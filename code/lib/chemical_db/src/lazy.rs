use std::{
    convert::TryFrom,
    fs::File,
    io::{Read, Write},
    path::{Path, PathBuf},
};

use csv::StringRecord;
use memmap2::Mmap;

use crate::{DatabaseRow, DbCsv};

/// A chemical database that is loaded from disk as needed.
///
/// The regular [`crate::ChemicalDb`] reads a full database file from disk and then
/// operates on it in memory. For extremely large databases this is not feasible. This database
/// just wraps a file path on disk and delegates all other responsibilities to [`LazyRows`].
///
/// # Iteration
/// This can be converted to [`LazyRows`], which supports iteration
#[derive(Debug, Clone)]
pub struct LazyDb {
    /// The chemical database file path
    pub path: PathBuf,
}

impl LazyDb {
    /// Construct a new [`LazyDb`] pointing that the given path.
    pub fn new(path: PathBuf) -> Self {
        Self { path }
    }

    /// Count the number of records in this database.
    ///
    /// This actually counts the number of newlines in the file, so this is not a free operation.
    ///
    /// # Errors
    /// Propagates IO errors from opening or memory mapping the file.
    pub fn count(&self) -> std::io::Result<usize> {
        let f = File::open(&self.path)?;
        let mmap = unsafe { Mmap::map(&f)? };
        Ok(memchr::memchr_iter(b'\n', &mmap).count() - 1)
    }

    /// Attempt to construct an iterator over rows in this database.
    ///
    /// Iteration is over tuples of row indices and [`DatabaseRow`]s.
    ///
    /// # Errors
    /// See [`LazyRows::from_path`].
    pub fn rows(&self) -> std::io::Result<LazyRows<File>> {
        self.try_into()
    }

    /// Get the row at `index`.
    ///
    /// This internally uses [`LazyDb::rows`] to get an iterator over rows in the database. If the
    /// provided index doesn't exist in the database and no errors occur this returns `Ok(None)`.
    ///
    /// # Errors
    /// Errors from [`LazyDb::rows`] are converted into [`csv::Error`]s and propagated. Any
    /// deserialization errors when reading the row of interest are propagated.
    pub fn get_row(&self, index: usize) -> csv::Result<Option<DatabaseRow>> {
        self.rows()?.nth(index).map(|(_, row)| row).transpose()
    }
}

impl TryFrom<LazyDb> for LazyRows<File> {
    type Error = std::io::Error;

    fn try_from(value: LazyDb) -> Result<Self, Self::Error> {
        Self::from_path(&value.path)
    }
}

impl TryFrom<&LazyDb> for LazyRows<File> {
    type Error = std::io::Error;

    fn try_from(value: &LazyDb) -> Result<Self, Self::Error> {
        Self::from_path(&value.path)
    }
}

/// An iterator over tuples of row indices and a database row references.
///
/// This iterator is intended to be re-constructed as necessary from a [`LazyDb`], since it cannot
/// be cloned due to its internal reader.
#[derive(Debug)]
pub struct LazyRows<R: Read> {
    reader: csv::Reader<R>,
    record: StringRecord,
    row: usize,
}

impl<R: Read> LazyRows<R> {
    /// Construct a lazy database from a reader.
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader: csv::Reader::from_reader(reader),
            row: 0,
            record: StringRecord::new(),
        }
    }

    /// Construct a lazy database from a file at the provided path.
    ///
    /// # Errors
    /// Propagates any I/O errors.
    pub fn from_path(path: impl AsRef<Path>) -> std::io::Result<LazyRows<File>> {
        let file = File::open(path.as_ref())?;
        Ok(LazyRows {
            reader: csv::Reader::from_reader(file),
            row: 0,
            record: StringRecord::new(),
        })
    }

    /// Writes this database as as CSV to the provided `writer`.
    ///
    /// # Note
    /// The [`csv`] crate automatically buffers all writers so the provided writer should not be
    /// buffered.
    ///
    /// # Errors
    /// Propagates errors from [`csv`] serialization.
    pub fn write_csv<W: Write>(mut self, writer: W) -> csv::Result<()> {
        let mut writer = csv::Writer::from_writer(writer);
        for row in self.reader.deserialize::<DatabaseRow>() {
            writer.serialize(row?)?;
        }
        Ok(())
    }
}

impl<R: Read> Iterator for LazyRows<R> {
    type Item = (usize, csv::Result<DatabaseRow>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record(&mut self.record) {
            Err(e) => Some((self.row, Err(e))),
            Ok(false) => None,
            Ok(true) => {
                let item = Some((
                    self.row,
                    self.record.deserialize(self.reader.headers().ok()),
                ));
                self.row += 1;
                item
            }
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        for _ in 0..n {
            match self.reader.read_record(&mut self.record) {
                Err(e) => return Some((self.row, Err(e))),
                Ok(false) => return None,
                Ok(true) => (),
            };
            self.row += 1;
        }
        self.next()
    }
}

impl TryFrom<DbCsv> for LazyRows<File> {
    type Error = std::io::Error;

    fn try_from(csv: DbCsv) -> Result<Self, Self::Error> {
        Self::from_path(csv.path)
    }
}

impl TryFrom<&DbCsv> for LazyRows<File> {
    type Error = std::io::Error;

    fn try_from(csv: &DbCsv) -> Result<Self, Self::Error> {
        Self::from_path(&csv.path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lazy_rows_nth() {
        let p = "./test_files/test_chem_db.csv".parse().unwrap();
        let mut rows: LazyRows<File> = LazyDb::new(p).try_into().unwrap();

        let (idx, row) = rows.nth(0).unwrap();
        assert_eq!(idx, 0);
        assert_eq!(row.unwrap().inchikey14(), "GWGKNTICBPKKKW");
        let (idx, row) = rows.nth(0).unwrap();
        assert_eq!(idx, 1);
        assert_eq!(row.unwrap().inchikey14(), "GRJSOZDXIUZXEW");
        let (idx, row) = rows.nth(5).unwrap();
        assert_eq!(idx, 7);
        assert_eq!(row.unwrap().inchikey14(), "QURRTAYEASAREY");

        let p = "./test_files/test_chem_db.csv".parse().unwrap();
        let mut rows: LazyRows<File> = LazyDb::new(p).try_into().unwrap();
        let (idx, row) = rows.nth(1).unwrap();
        assert_eq!(idx, 1);
        assert_eq!(row.unwrap().inchikey14(), "GRJSOZDXIUZXEW");
        let (idx, row) = rows.nth(1).unwrap();
        assert_eq!(idx, 3);
        assert_eq!(row.unwrap().inchikey14(), "ZRWKFXOGNFQPMY");

        assert!(rows.nth(200).is_none());
        assert!(rows.nth(0).is_none());
        assert!(rows.next().is_none());
    }
}
