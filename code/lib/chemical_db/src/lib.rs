//! Contains functionality to parse and query a database of molecules.
//!
//! The main [`ChemicalDb`] struct stores the whole database in memory, but with molecules in the
//! SMILES format. Conversion of database elements to usable [`Mol`](molecule::Mol) structures must
//! be done manually.
//!
//! # Example
//! ```
//! use chemical_db::{ChemicalDb, DatabaseRow};
//!
//! let csv_contents = r#"smiles,exact_mass,name,inchikey14
//! COC1=C2C=COC2=CC=C1C(=O)CC(OC)C1=CC=CC=C1,310.120509056,CCMSLIB00000847598,AADNEQWIZKTMBL
//! COC1C(CO)OC(N2C=NC3=C2N=C[NH]C3=O)C1O,282.096419548,CCMSLIB00000855402,ABXDBVMGRKZFRC
//! "#;
//!
//! // we can only construct databases by parsing a CSV
//! let db = ChemicalDb::from_reader(csv_contents.as_bytes()).unwrap();
//!
//! // and check some important summary statistics
//! assert!(!db.is_empty());
//! assert_eq!(db.len(), 2);
//!
//! // the input is sorted by mass
//! assert_eq!(db.get(0).unwrap().exact_mass(), 282.096419548);
//! assert_eq!(db.get(0).unwrap().name(), "CCMSLIB00000855402");
//! assert_eq!(db.get(0).unwrap().inchikey14(), "ABXDBVMGRKZFRC");
//! assert_eq!(db.get(1).unwrap().exact_mass(), 310.120509056);
//! assert_eq!(db.get(1).unwrap().name(), "CCMSLIB00000847598");
//! assert_eq!(db.get(1).unwrap().inchikey14(), "AADNEQWIZKTMBL");
//!
//! // accessing out of bounds returns None
//! assert_eq!(db.get(3), None);
//!
//! // we can easily query masses in a error-tolerant way
//! assert_eq!(db.query(310.2, 0.1), Some(1..2));
//! assert_eq!(db.query(300.0, 100.0), Some(0..2));
//! assert_eq!(db.query(400.0, 0.1), None);
//! ```

#![warn(
    missing_docs,
    clippy::all,
    clippy::pedantic,
    rust_2018_idioms,
    future_incompatible
)]
#![allow(
    clippy::must_use_candidate,
    clippy::module_name_repetitions,
    clippy::missing_panics_doc,
    clippy::shadow_unrelated,
    clippy::unseparated_literal_suffix
)]
#![cfg_attr(
    test,
    allow(clippy::pedantic, clippy::similar_names, clippy::float_cmp,)
)]

mod db;
mod filter;
mod lazy;
mod row;

pub use db::{ChemicalDb, DbCsv, SingletonDb};
pub use filter::{CandidateFilter, ErrorTol};
pub use lazy::{LazyDb, LazyRows};
pub use row::DatabaseRow;

#[cfg(test)]
mod tests;
