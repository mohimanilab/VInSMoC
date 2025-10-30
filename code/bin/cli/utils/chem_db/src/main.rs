//! Utility methods for working with chemical databases.

use std::fs::File;

use anyhow::Result;
use chemical_db::{CandidateFilter, ChemicalDb, DatabaseRow, DbCsv, LazyDb};
use clap::Parser;
use dereplicate::db::AsCandidates;
use memmap2::Mmap;
use rustc_hash::FxHashSet;
use spectrum::ParsedCollections;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Utility methods for working with chemical databases.
///
/// Chemical databases for our use cases are CSVs with at least four columns: smiles, exact_mass,
/// inchikey14, and name. Additional columns can be present but will not be used. The smiles column
/// is used to represent the molecular structure. The name column is for any sort of human-readable
/// annotation. While the exact_mass and inchikey14 columns could be computed from the smiles
/// column, they are required to avoid extra work that can easily be done ahead of time.
#[derive(Debug, Parser)]
#[clap(name = "chem_db", version = env!("CARGO_PKG_VERSION"))]
struct Opts {
    #[clap(flatten)]
    db: DbCsv,
    #[clap(subcommand)]
    sub: Subcommand,
}

#[derive(Debug, Parser)]
enum Subcommand {
    /// Print the number of molecules in this database.
    Length,
    /// Check validity of all entries in this database.
    IsValid,
    /// Find all entries with mass within --error of the input --mass.
    FindMass {
        /// Error tolerance to report chemical database entries.
        #[clap(long, required = true)]
        error: f64,
        /// Target mass.
        #[clap(long, required = true)]
        mass: f64,
    },
    /// Find all entries whose mass's integer component is --mass.
    FindTruncMass {
        /// Integer component of the target mass.
        #[clap(long, required = true)]
        mass: u64,
    },
    /// Sort this database by ascending mass.
    Sort,
    /// Remove all duplicates from this database based on InChIKey.
    Uniq {
        /// Report unique database base sorted by ascending mass.
        #[clap(long, action)]
        sort: bool,
    },
    /// Computes candidates that match to the provided spectral collections.
    ///
    /// Results are reported as a four-column CSV. The first two columns contain the spectrum file
    /// name and scan. The third column is the 0-based index of the matched molecule in the
    /// database. Finally, the fourth column is the adduct used to match the spectrum to the
    /// candidate molecule.
    Candidates {
        #[clap(flatten)]
        filtering: CandidateFilter,
        #[clap(flatten)]
        spectra: ParsedCollections,
    },
    /// Computes number of molecule-spectrum candidates that would be computed given filtering
    /// parameters.
    CandidateCount {
        #[clap(flatten)]
        filtering: CandidateFilter,
        #[clap(flatten)]
        spectra: ParsedCollections,
    },
}

impl Opts {
    pub fn mmap(&self) -> Result<Mmap> {
        let f = File::open(&self.db.path)?;
        let mmap = unsafe { Mmap::map(&f)? };
        Ok(mmap)
    }
}

fn main() -> Result<()> {
    let opts = Opts::parse();

    match opts.sub {
        Subcommand::Length => {
            let mmap = opts.mmap()?;
            println!("{}", memchr::memchr_iter(b'\n', &mmap).count());
        }
        Subcommand::IsValid => {
            let mmap = opts.mmap()?;
            let mut reader = csv::Reader::from_reader(mmap.as_ref());

            for row in reader.deserialize::<DatabaseRow>() {
                row?;
            }
        }
        Subcommand::FindMass { error, mass } => {
            let mmap = opts.mmap()?;
            let mut reader = csv::Reader::from_reader(mmap.as_ref());

            let stdout = std::io::stdout();
            let stdout = stdout.lock();
            let mut writer = csv::Writer::from_writer(stdout);

            for row in reader.deserialize::<DatabaseRow>() {
                let row = row?;
                if (row.exact_mass() - mass).abs() <= error {
                    writer.serialize(&row)?;
                }
            }
        }
        Subcommand::FindTruncMass { mass } => {
            let mmap = opts.mmap()?;
            let mut reader = csv::Reader::from_reader(mmap.as_ref());

            let stdout = std::io::stdout();
            let stdout = stdout.lock();
            let mut writer = csv::Writer::from_writer(stdout);

            for row in reader.deserialize::<DatabaseRow>() {
                let row = row?;
                if row.exact_mass().trunc() as u64 == mass {
                    writer.serialize(&row)?;
                }
            }
        }
        Subcommand::Sort => {
            let mut db = ChemicalDb::try_from(opts.db)?;
            db.sort_rows();

            let stdout = std::io::stdout();
            let stdout = stdout.lock();
            let mut writer = csv::Writer::from_writer(stdout);

            for row in &db {
                writer.serialize(row)?;
            }
        }
        Subcommand::Uniq { sort } => {
            let db = ChemicalDb::try_from(opts.db)?;
            let mut rows = db
                .into_iter()
                .collect::<FxHashSet<_>>()
                .into_iter()
                .collect::<ChemicalDb>();
            if sort {
                rows.sort_rows();
            }

            let stdout = std::io::stdout();
            let stdout = stdout.lock();
            let mut writer = csv::Writer::from_writer(stdout);

            for row in &rows {
                writer.serialize(row)?;
            }
        }
        Subcommand::Candidates {
            filtering,
            mut spectra,
        } => {
            let db = LazyDb {
                path: opts.db.path.into(),
            };

            println!("spectrum_file,scan,mol_idx,adduct");

            for spectra in &mut spectra {
                let (spectrum_file, collection) = spectra?;
                let spectrum_file = spectrum_file.canonicalize()?.display().to_string();

                for candidate in db.candidates(&collection, &filtering)? {
                    let ((mol_idx, _), matched_spectra) = candidate?;
                    for (adduct, scan) in matched_spectra {
                        println!("{},{},{},{}", &spectrum_file, scan, mol_idx, adduct);
                    }
                }
            }
        }
        Subcommand::CandidateCount {
            filtering,
            mut spectra,
        } => {
            let db = LazyDb {
                path: opts.db.path.into(),
            };

            let mut candidate_count = 0;

            for spectra in &mut spectra {
                let (_, collection) = spectra?;

                for hit in db.candidates(&collection, &filtering)? {
                    let (_, matched_spectra) = hit?;
                    candidate_count += matched_spectra.len();
                }
            }

            println!("{}", candidate_count);
        }
    }

    Ok(())
}
