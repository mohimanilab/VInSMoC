use std::{convert::TryFrom, fmt, str::FromStr};

use constants::amino_acids::NUM_BACKBONE_ATOMS;
use learning::data::CategoricalFeature;
use modification::{Mod, ModSeq};
use molecule::{Atom, BondType, Mol};
use ndarray::{arr1, Array1};
use serde::{Deserialize, Serialize};

use super::AminoAcid;

mod rgroup;

/// One of the 20 standard amino acids.
///
/// See [`AminoAcid`] for general requirements for a struct to be an amino acid.
///
/// # Examples
/// ```
/// use peptide::{AminoAcid, StdAminoAcid};
/// use molecule::Atom;
/// use constants::MASS_WATER_EXACT;
///
/// // we can parse these amino acids from strings
/// let g = "G".parse::<StdAminoAcid>().unwrap();
/// assert_eq!(g, StdAminoAcid::G);
///
/// // glycine's R-group is just a hydrogen
/// assert_eq!(g.r_group_exact_mass(), Atom::H.exact_mass());
/// assert_eq!(g.exact_mass(), g.residue_exact_mass() + MASS_WATER_EXACT);
/// ```
#[derive(Debug, Eq, PartialEq, Hash, Clone, Copy, Serialize, Deserialize)]
pub enum StdAminoAcid {
    /// Arginine
    R,
    /// Histidine
    H,
    /// Lysine
    K,
    /// Aspartic acid
    D,
    /// Glutamic acid
    E,
    /// Serine
    S,
    /// Threonine
    T,
    /// Asparagine
    N,
    /// Glutamine
    Q,
    /// Cysteine
    C,
    /// Glycine
    G,
    /// Proline
    P,
    /// Alanine
    A,
    /// Valine
    V,
    /// Isoleucine
    I,
    /// Leucine
    L,
    /// Methionine
    M,
    /// Phenylalanine
    F,
    /// Tyrosine
    Y,
    /// Tryptophan
    W,
}

impl CategoricalFeature for StdAminoAcid {
    const NUM_CATEGORIES: usize = 20;

    /// These are in alphabetical order to match with the format of HMMER's HMM files.
    #[inline]
    fn index(&self) -> usize {
        match *self {
            Self::A => 0,
            Self::C => 1,
            Self::D => 2,
            Self::E => 3,
            Self::F => 4,
            Self::G => 5,
            Self::H => 6,
            Self::I => 7,
            Self::K => 8,
            Self::L => 9,
            Self::M => 10,
            Self::N => 11,
            Self::P => 12,
            Self::Q => 13,
            Self::R => 14,
            Self::S => 15,
            Self::T => 16,
            Self::V => 17,
            Self::W => 18,
            Self::Y => 19,
        }
    }

    /// These are in alphabetical order to match with the format of HMMER's HMM files.
    #[inline]
    fn from_index(index: usize) -> Option<Self> {
        match index {
            0 => Some(Self::A),
            1 => Some(Self::C),
            2 => Some(Self::D),
            3 => Some(Self::E),
            4 => Some(Self::F),
            5 => Some(Self::G),
            6 => Some(Self::H),
            7 => Some(Self::I),
            8 => Some(Self::K),
            9 => Some(Self::L),
            10 => Some(Self::M),
            11 => Some(Self::N),
            12 => Some(Self::P),
            13 => Some(Self::Q),
            14 => Some(Self::R),
            15 => Some(Self::S),
            16 => Some(Self::T),
            17 => Some(Self::V),
            18 => Some(Self::W),
            19 => Some(Self::Y),
            _ => None,
        }
    }
}

/// Number of standard amino acids.
pub const NUM_STANDARD_AMINO_ACIDS: usize = 20;

impl AminoAcid for StdAminoAcid {
    fn r_group(&self) -> Mol {
        match *self {
            Self::R => rgroup::r(),
            Self::H => rgroup::h(),
            Self::K => rgroup::k(),
            Self::D => rgroup::d(),
            Self::E => rgroup::e(),
            Self::S => rgroup::s(),
            Self::T => rgroup::t(),
            Self::N => rgroup::n(),
            Self::Q => rgroup::q(),
            Self::C => rgroup::c(),
            Self::G => rgroup::g(),
            Self::P => rgroup::p(),
            Self::A => rgroup::a(),
            Self::V => rgroup::v(),
            Self::I => rgroup::i(),
            Self::L => rgroup::l(),
            Self::M => rgroup::m(),
            Self::F => rgroup::f(),
            Self::Y => rgroup::y(),
            Self::W => rgroup::w(),
        }
    }

    fn connection_mod(&self) -> ModSeq {
        // we always connect alpha carbon to first R-group atom
        let mut mods = vec![Mod::AddEdge(0, NUM_BACKBONE_ATOMS as i32, BondType::Single)];

        // special cases
        match self {
            // glycine has an empty R-group, so just change hydrogens on backbone
            Self::G => {
                mods = vec![Mod::HydrogenShift(0, 1)];
            }
            // proline has a ring connection to backbone nitrogen
            // this means we need an extra edge and need to remove the hydrogen on the nitrogen
            Self::P => {
                mods.push(Mod::HydrogenShift(1, -1));
                mods.push(Mod::AddEdge(
                    1,
                    NUM_BACKBONE_ATOMS as i32 + 2,
                    BondType::Single,
                ));
            }
            _ => (),
        }
        ModSeq::new(mods)
    }
}

impl StdAminoAcid {
    /// Returns the standard mass of the amino acid side chain.
    ///
    /// This does not include any part of the backbone.
    pub fn r_group_std_mass(&self) -> f64 {
        match *self {
            Self::G => self.r_group().std_mass() + Atom::H.std_mass(),
            Self::P => self.r_group().std_mass() - Atom::H.std_mass(),
            _ => self.r_group().std_mass(),
        }
    }

    /// Returns the exact mass of the amino acid side chain.
    pub fn r_group_exact_mass(&self) -> f64 {
        match *self {
            Self::G => self.r_group().exact_mass() + Atom::H.exact_mass(),
            Self::P => self.r_group().exact_mass() - Atom::H.exact_mass(),
            _ => self.r_group().exact_mass(),
        }
    }

    /// Gives one hot encoding of an amino acid.
    pub fn one_hot_encoding(&self) -> Array1<u8> {
        let mut ohe_array = arr1(&[0; 20]);
        let idx = self.index();
        ohe_array[idx] = 1;
        ohe_array
    }

    /// Returns the index of the alcohol group in the side chain of the residue, if it exists.
    pub fn sidechain_alcohol_index(&self) -> Option<usize> {
        match *self {
            Self::S => Some(5),
            Self::T => Some(5),
            _ => None,
        }
    }

    /// Returns the index of the thiol group in the side chain of the residue, if it exists.
    pub fn sidechain_thiol_index(&self) -> Option<usize> {
        match *self {
            Self::C => Some(5),
            _ => None,
        }
    }

    /// Returns the optional index of the start of the amino group in the residue.
    ///
    /// This returns the index of the amino nitrogen in the side chain of the residue, if it
    /// exists.
    pub fn sidechain_amino_nitrogen_index(&self) -> Option<usize> {
        match *self {
            Self::K => Some(8),
            _ => None,
        }
    }

    /// Returns the optional index of the start of the carboxyl group in the residue.
    ///
    /// This returns the index of the carboxy carbon in the side chain of the residue, if it
    /// exists.
    pub fn sidechain_carboxy_carbon_index(&self) -> Option<usize> {
        match *self {
            Self::D => Some(5),
            Self::E => Some(6),
            _ => None,
        }
    }
}

/// This error represents an error with parsing a standard amino acid.
///
/// The derived traits enable us to set this error message as a source for downstream errors.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ParseStdAminoAcidError {
    cause: String,
}

impl ParseStdAminoAcidError {
    /// Expose the cause of the error, to use in downstream errors.
    pub fn get_cause(&self) -> String {
        self.cause.clone()
    }
}

impl std::fmt::Display for ParseStdAminoAcidError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Invalid amino acid string {}, must be a valid case-insensitive single- or three-letter amino acid code", self.cause)
    }
}

impl std::error::Error for ParseStdAminoAcidError {}

impl FromStr for StdAminoAcid {
    type Err = ParseStdAminoAcidError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_ref() {
            "R" | "ARG" => Ok(Self::R), // Arginine
            "H" | "HIS" => Ok(Self::H), // Histidine
            "K" | "LYS" => Ok(Self::K), // Lysine
            "D" | "ASP" => Ok(Self::D), // Aspartic acid
            "E" | "GLU" => Ok(Self::E), // Glutamic acid
            "S" | "SER" => Ok(Self::S), // Serine
            "T" | "THR" => Ok(Self::T), // Threonine
            "N" | "ASN" => Ok(Self::N), // Asparagine
            "Q" | "GLN" => Ok(Self::Q), // Glutamine
            "C" | "CYS" => Ok(Self::C), // Cysteine
            "G" | "GLY" => Ok(Self::G), // Glycine
            "P" | "PRO" => Ok(Self::P), // Proline
            "A" | "ALA" => Ok(Self::A), // Alanine
            "V" | "VAL" => Ok(Self::V), // Valine
            "I" | "ILE" => Ok(Self::I), // Isoleucine
            "L" | "LEU" => Ok(Self::L), // Leucine
            "M" | "MET" => Ok(Self::M), // Methionine
            "F" | "PHE" => Ok(Self::F), // Phenylalanine
            "Y" | "TYR" => Ok(Self::Y), // Tyrosine
            "W" | "TRP" => Ok(Self::W), // tryptophan
            _ => Err(Self::Err {
                cause: s.to_string(),
            }), // invalid string
        }
    }
}

impl fmt::Display for StdAminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Self::R => "ARG", // Arginine
            Self::H => "HIS", // Histidine
            Self::K => "LYS", // Lysine
            Self::D => "ASP", // Aspartic acid
            Self::E => "GLU", // Glutamic acid
            Self::S => "SER", // Serine
            Self::T => "THR", // Threonine
            Self::N => "ASN", // Asparagine
            Self::Q => "GLN", // Glutamine
            Self::C => "CYS", // Cysteine
            Self::G => "GLY", // Glycine
            Self::P => "PRO", // Proline
            Self::A => "ALA", // Alanine
            Self::V => "VAL", // Valine
            Self::I => "ILE", // Isoleucine
            Self::L => "LEU", // Leucine
            Self::M => "MET", // Methionine
            Self::F => "PHE", // Phenylalanine
            Self::Y => "TYR", // Tyrosine
            Self::W => "TRP", // tryptophan
        };
        write!(f, "{}", s)
    }
}

impl TryFrom<char> for StdAminoAcid {
    type Error = ParseStdAminoAcidError;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'R' | 'r' => Ok(Self::R), // Arginine
            'H' | 'h' => Ok(Self::H), // Histidine
            'K' | 'k' => Ok(Self::K), // Lysine
            'D' | 'd' => Ok(Self::D), // Aspartic acid
            'E' | 'e' => Ok(Self::E), // Glutamic acid
            'S' | 's' => Ok(Self::S), // Serine
            'T' | 't' => Ok(Self::T), // Threonine
            'N' | 'n' => Ok(Self::N), // Asparagine
            'Q' | 'q' => Ok(Self::Q), // Glutamine
            'C' | 'c' => Ok(Self::C), // Cysteine
            'G' | 'g' => Ok(Self::G), // Glycine
            'P' | 'p' => Ok(Self::P), // Proline
            'A' | 'a' => Ok(Self::A), // Alanine
            'V' | 'v' => Ok(Self::V), // Valine
            'I' | 'i' => Ok(Self::I), // Isoleucine
            'L' | 'l' => Ok(Self::L), // Leucine
            'M' | 'm' => Ok(Self::M), // Methionine
            'F' | 'f' => Ok(Self::F), // Phenylalanine
            'Y' | 'y' => Ok(Self::Y), // Tyrosine
            'W' | 'w' => Ok(Self::W), // tryptophan
            _ => Err(Self::Error {
                cause: value.to_string(),
            }), // invalid string
        }
    }
}

impl From<StdAminoAcid> for char {
    fn from(aa: StdAminoAcid) -> Self {
        match aa {
            StdAminoAcid::R => 'R',
            StdAminoAcid::H => 'H',
            StdAminoAcid::K => 'K',
            StdAminoAcid::D => 'D',
            StdAminoAcid::E => 'E',
            StdAminoAcid::S => 'S',
            StdAminoAcid::T => 'T',
            StdAminoAcid::N => 'N',
            StdAminoAcid::Q => 'Q',
            StdAminoAcid::C => 'C',
            StdAminoAcid::G => 'G',
            StdAminoAcid::P => 'P',
            StdAminoAcid::A => 'A',
            StdAminoAcid::V => 'V',
            StdAminoAcid::I => 'I',
            StdAminoAcid::L => 'L',
            StdAminoAcid::M => 'M',
            StdAminoAcid::F => 'F',
            StdAminoAcid::Y => 'Y',
            StdAminoAcid::W => 'W',
        }
    }
}

#[cfg(test)]
mod test;
