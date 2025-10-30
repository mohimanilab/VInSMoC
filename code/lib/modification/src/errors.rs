use crate::Mod;
use molecule::{BondType, BuildMolError, NodeIndex};
use std::{fmt, io, path::PathBuf};

/// Represents a kind of error that can happen when modifying a molecule.
#[derive(Debug, Eq, PartialEq, Clone)]
pub enum ModificationErrorKind {
    /// The given node did not exist
    NoSuchNode(NodeIndex),
    /// There was no edge between the two nodes
    NoSuchEdge(NodeIndex, NodeIndex),
    /// Attempted to apply a modification with a negative index
    NegativeIndex(i32),
    /// Attempted to use negative node index with no assoicated true node index
    NoIndexMapping(i32),
    /// Attempted to shift hydrogen count to a negative number
    NegativeHydrogen(NodeIndex),
    /// Attempted to connect an atom to itself
    SelfLoop(i32),
}

/// An error that can occur when modifying a molecule.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct ModificationError {
    /// The modification that cause this error
    pub cause: Mod,
    /// The explanation for this error
    pub kind: ModificationErrorKind,
}

impl ModificationError {
    /// Create a new modification error.
    pub fn new(cause: Mod, kind: ModificationErrorKind) -> Self {
        Self { cause, kind }
    }
}

impl fmt::Display for ModificationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.kind {
            ModificationErrorKind::NoSuchNode(x) => write!(
                f,
                "no node with index {} when applying modification {:?}",
                x.index(),
                self.cause
            ),
            ModificationErrorKind::NoSuchEdge(x, y) => write!(
                f,
                "no edge connecting node indices {} and {} when applying modification {:?}",
                x.index(),
                y.index(),
                self.cause
            ),
            ModificationErrorKind::NegativeIndex(x) => {
                write!(f, "attempted to directly use negative node index {}", x)
            }
            ModificationErrorKind::NoIndexMapping(x) => {
                write!(f, "attempted to use negative node index {} that has no associated true node index", x)
            }
            ModificationErrorKind::NegativeHydrogen(x) => {
                write!(f, "attempted to shift hydrogens to a negative hydrogen count for node index {} when applying modification {:?}", x.index(), self.cause)
            }
            ModificationErrorKind::SelfLoop(x) => {
                write!(
                    f,
                    "attempted to connect node {} to itself when applying modification {:?}",
                    x, self.cause
                )
            }
        }
    }
}

impl std::error::Error for ModificationError {}

#[derive(Debug, thiserror::Error)]
/// A general error that encapsulates parsing, serialization, and modification
/// application issues.
pub enum ContextualError {
    /// Failed to read the mods directory
    #[error("unable to read mods directory")]
    DataMods(#[source] io::Error),
    /// Failed to read the mod sequence file in the mod directory
    #[error("error handling mod.txt in {mod_dir}")]
    Mods {
        /// Directory path
        mod_dir: PathBuf,
        /// Source error
        #[source]
        source: Error,
    },
    /// Failed parsing the file containing the target molecule to apply mods on
    #[error("error parsing mod.mol in {mod_dir}")]
    Mol {
        /// Directory path
        mod_dir: PathBuf,
        /// Source error
        #[source]
        source: BuildMolError,
    },
    /// Failed applying the given modification to the target molecule
    #[error("error applying modification on line {line} to mod.mol in {mod_dir}")]
    Application {
        /// Directory path
        mod_dir: PathBuf,
        /// Line number of modification that we failed to apply
        line: usize,
        /// Source error
        #[source]
        source: ModificationError,
    },
    /// Failed applying the given modification to the target molecule after resetting
    /// indices
    #[error("error applying modification on mod.mol in {mod_dir} after resetting its indices")]
    TransApplication {
        /// Directory path
        mod_dir: PathBuf,
        /// Source error
        #[source]
        source: ModificationError,
    },
    /// Failed serializing converted modification in the mod directory
    #[error("failed to serialize converted modification from {mod_dir}")]
    Serde {
        /// Directory path
        mod_dir: PathBuf,
        /// Source error
        #[source]
        source: bincode::Error,
    },
    /// Failed parsing reaction file due to missing arrow `->`
    #[error(".smi file missing -> in reaction")]
    NoArrow,
    /// Failed indexing SMILES
    #[error("invalid indexed smiles")]
    IndexedSmiles(#[from] BuildMolError),
}

/// An error that involves issues in parsing, I/O, and invalid molecule accesses
/// and manipulations.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Issue parsing the given string
    #[error("unable to parse input\n {0}")]
    Parse(String),
    /// Issue in I/O
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Issue applying modification on the given line
    #[error("failed to apply modification on line {line}")]
    Validation {
        /// Line number of modification that we failed to apply due to validation
        /// issue
        line: usize,
        /// Source error
        #[source]
        source: ValidationError,
    },
}

/// An error for invalid indexing of a molecule and invalid bond operations.
#[derive(Debug, thiserror::Error)]
pub enum ValidationError {
    ///
    #[error("missing bonds for Add directive")]
    MissingBonds,
    #[error("duplicate node index {0} in Add directive")]
    DuplicateIndex(i32),
    #[error("attempted to connect hydrogen at {h} to atom {other} with a {bond} bond instead of a single bond")]
    HydrogenBond { h: i32, other: i32, bond: BondType },
    #[error("atom index {0} is invalid where it appears")]
    InvalidIndex(i32),
}
