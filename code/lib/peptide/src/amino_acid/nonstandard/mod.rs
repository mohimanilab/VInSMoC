use std::{fmt, str::FromStr};

use constants::amino_acids::NUM_BACKBONE_ATOMS;
use learning::data::CategoricalFeature;
use modification::{Mod, ModSeq};
use molecule::{BondType, Mol};
use serde::{Deserialize, Serialize};

use super::AminoAcid;

mod rgroup;

/// A non-standard amino acid.
///
/// See [`AminoAcid`] for general requirements for a struct to be an amino acid.
///
/// This is a general purpose enum that exists to logically separate these non-standard amino acids
/// from the 20 standard ones. Any new non-standard amino acids can be added as variants to this
/// enum without API breakage. These non-standard amino acids are not assigned a single-letter
/// code, but instead are assigned three-letter codes that match the string of their variant.
/// For example, Citrulline has a three letter code "Cit".
#[derive(Debug, Eq, PartialEq, Hash, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub enum NonStdAminoAcid {
    /// [2-Aminohexanedioic acid](https://pubchem.ncbi.nlm.nih.gov/compound/469)
    Aad,
    /// [L-2-Aminobutyric acid](https://pubchem.ncbi.nlm.nih.gov/compound/80283)
    Abu,
    /// [2-Amino-3-hydroxy-3-(4-hydroxyphenyl)propanoic acid](https://pubchem.ncbi.nlm.nih.gov/compound/13309269)
    Bht,
    /// [3,5-Dihydroxyphenylglycine](https://pubchem.ncbi.nlm.nih.gov/compound/108001)
    Dpg,
    /// [Oxfenicine](https://pubchem.ncbi.nlm.nih.gov/compound/36143)
    Hpg,
    /// [L-Ornithine](https://pubchem.ncbi.nlm.nih.gov/compound/6262)
    Orn,
    /// [L-Phenylglycine](https://pubchem.ncbi.nlm.nih.gov/compound/99291)
    Phg,
    /// [L-Pipecolic acid](https://pubchem.ncbi.nlm.nih.gov/compound/439227)
    Pip,
    /// [Citrulline](https://pubchem.ncbi.nlm.nih.gov/compound/9750)
    Cit,
    /// [L-2,4-Diaminobutyric acid](https://pubchem.ncbi.nlm.nih.gov/compound/134490)
    Dab,
    /// [2-Amino-4-(4-hydroxyphenyl)butanoic Acid](https://pubchem.ncbi.nlm.nih.gov/compound/4153395)
    Hty,
    /// [L-Homoserine](https://pubchem.ncbi.nlm.nih.gov/compound/12647)
    Hse,
    /// [hydroxy Ornithine](https://pubchem.ncbi.nlm.nih.gov/compound/N_5_-Hydroxy-L-ornithine)
    Hrn,
    /// [acetyl hydroxy ornithine](https://pubchem.ncbi.nlm.nih.gov/compound/N5-Acetyl-N5-hydroxy-L-ornithine)
    Han,
    /// [2_3-Diaminopropionic-acid](https://pubchem.ncbi.nlm.nih.gov/compound/2_3-Diaminopropionic-acid)
    Dap,
    /// [Beta-hydroxy enduracidine](https://pubchem.ncbi.nlm.nih.gov/substance/242082644)
    Bnd,
    /// [2-Amino-3,5-dimethyl-4-hexenoic Acid](https://pubchem.ncbi.nlm.nih.gov/compound/2-Amino-3_5-dimethylhex-4-enoic-acid)
    Adh,
    /// [1-(1,1-Dimethylallyl)-L-tryptophan](https://pubchem.ncbi.nlm.nih.gov/compound/22216483)
    Dmw,
    /// [dehydrobutyrine](https://pubchem.ncbi.nlm.nih.gov/compound/6449989)
    Dhb,
    /// [5-chloro tryptophan](https://pubchem.ncbi.nlm.nih.gov/compound/5-Chloro-L-tryptophan)
    Crp,
    /// [2-Amino 5-MethylHexanoic Acid](https://pubchem.ncbi.nlm.nih.gov/compound/220783)
    Mha,
    /// [hydroxy formyl ornithine](https://pubchem.ncbi.nlm.nih.gov/compound/46173094)
    Hfn,
    /// [4-methyl proline](https://pubchem.ncbi.nlm.nih.gov/compound/4-Methyl-proline)
    Mpr,
    /// [beta-alanine](https://pubchem.ncbi.nlm.nih.gov/compound/beta-Alanine)
    Bla,
    /// [beta butyric acid](https://pubchem.ncbi.nlm.nih.gov/compound/10932)
    Bbu,
    /// [2-amino decanoic acid](https://pubchem.ncbi.nlm.nih.gov/compound/307923)
    Ade,
    /// [2-Amino-3-(pyridin-3-yl)propanoic acid](https://pubchem.ncbi.nlm.nih.gov/compound/314100)
    Pyr,
    /// [croto glycine](https://pubchem.ncbi.nlm.nih.gov/compound/168747)
    Cgy,
    /// [o methyl tyrosine](https://pubchem.ncbi.nlm.nih.gov/compound/O-Methyl-L-tyrosine)
    Omy,
    /// [piperazic acid](https://pubchem.ncbi.nlm.nih.gov/compound/1712209)
    Piz,
}

impl CategoricalFeature for NonStdAminoAcid {
    const NUM_CATEGORIES: usize = 30;

    #[inline]
    fn index(&self) -> usize {
        match *self {
            Self::Aad => 0,
            Self::Abu => 1,
            Self::Bht => 2,
            Self::Dpg => 3,
            Self::Hpg => 4,
            Self::Orn => 5,
            Self::Phg => 6,
            Self::Pip => 7,
            Self::Cit => 8,
            Self::Dab => 9,
            Self::Hty => 10,
            Self::Hse => 11,
            Self::Hrn => 12,
            Self::Han => 13,
            Self::Dap => 14,
            Self::Bnd => 15,
            Self::Adh => 16,
            Self::Dmw => 17,
            Self::Dhb => 18,
            Self::Crp => 19,
            Self::Mha => 20,
            Self::Hfn => 21,
            Self::Mpr => 22,
            Self::Bla => 23,
            Self::Bbu => 24,
            Self::Ade => 25,
            Self::Pyr => 26,
            Self::Cgy => 27,
            Self::Omy => 28,
            Self::Piz => 29,
        }
    }

    #[inline]
    fn from_index(index: usize) -> Option<Self> {
        match index {
            0 => Some(Self::Aad),
            1 => Some(Self::Abu),
            2 => Some(Self::Bht),
            3 => Some(Self::Dpg),
            4 => Some(Self::Hpg),
            5 => Some(Self::Orn),
            6 => Some(Self::Phg),
            7 => Some(Self::Pip),
            8 => Some(Self::Cit),
            9 => Some(Self::Dab),
            10 => Some(Self::Hty),
            11 => Some(Self::Hse),
            12 => Some(Self::Hrn),
            13 => Some(Self::Han),
            14 => Some(Self::Dap),
            15 => Some(Self::Bnd),
            16 => Some(Self::Adh),
            17 => Some(Self::Dmw),
            18 => Some(Self::Dhb),
            19 => Some(Self::Crp),
            20 => Some(Self::Mha),
            21 => Some(Self::Hfn),
            22 => Some(Self::Mpr),
            23 => Some(Self::Bla),
            24 => Some(Self::Bbu),
            25 => Some(Self::Ade),
            26 => Some(Self::Pyr),
            27 => Some(Self::Cgy),
            28 => Some(Self::Omy),
            29 => Some(Self::Piz),
            _ => None,
        }
    }
}

impl AminoAcid for NonStdAminoAcid {
    fn r_group(&self) -> Mol {
        match *self {
            Self::Aad => rgroup::aad(),
            Self::Abu => rgroup::abu(),
            Self::Bht => rgroup::bht(),
            Self::Dpg => rgroup::dpg(),
            Self::Hpg => rgroup::hpg(),
            Self::Orn => rgroup::orn(),
            Self::Phg => rgroup::phg(),
            Self::Pip => rgroup::pip(),
            Self::Cit => rgroup::cit(),
            Self::Dab => rgroup::dab(),
            Self::Hty => rgroup::hty(),
            Self::Hse => rgroup::hse(),
            Self::Hrn => rgroup::hrn(),
            Self::Han => rgroup::han(),
            Self::Dap => rgroup::dap(),
            Self::Bnd => rgroup::bnd(),
            Self::Adh => rgroup::adh(),
            Self::Dmw => rgroup::dmw(),
            Self::Dhb => rgroup::dhb(),
            Self::Crp => rgroup::crp(),
            Self::Mha => rgroup::mha(),
            Self::Hfn => rgroup::hfn(),
            Self::Mpr => rgroup::mpr(),
            Self::Bla => rgroup::bla(),
            Self::Bbu => rgroup::bbu(),
            Self::Ade => rgroup::ade(),
            Self::Pyr => rgroup::pyr(),
            Self::Cgy => rgroup::cgy(),
            Self::Omy => rgroup::omy(),
            Self::Piz => rgroup::piz(),
        }
    }

    fn connection_mod(&self) -> ModSeq {
        let mut mods = match *self {
            Self::Dhb => vec![
                Mod::HydrogenShift(0, -1),
                Mod::AddEdge(0, NUM_BACKBONE_ATOMS as i32, BondType::Double),
            ],
            _ => vec![Mod::AddEdge(0, NUM_BACKBONE_ATOMS as i32, BondType::Single)],
        };

        if *self == Self::Pip || *self == Self::Piz {
            mods.push(Mod::HydrogenShift(1, -1));
            mods.push(Mod::AddEdge(
                1,
                NUM_BACKBONE_ATOMS as i32 + 3,
                BondType::Single,
            ));
        }

        if *self == Self::Mpr {
            mods.push(Mod::HydrogenShift(1, -1));
            mods.push(Mod::AddEdge(
                1,
                NUM_BACKBONE_ATOMS as i32 + 2,
                BondType::Single,
            ));
        }

        if *self == Self::Bla {
            mods.push(Mod::RemoveEdge(0, 1));
            mods.push(Mod::AddEdge(1, NUM_BACKBONE_ATOMS as i32, BondType::Single));
            mods.push(Mod::HydrogenShift(0, 1));
        }

        if *self == Self::Bbu {
            mods.push(Mod::RemoveEdge(0, 1));
            mods.push(Mod::AddEdge(1, NUM_BACKBONE_ATOMS as i32, BondType::Single));
            mods.push(Mod::HydrogenShift(0, 1));
        }

        ModSeq::new(mods)
    }
}

impl NonStdAminoAcid {
    /// Returns the index of the amino nitrogen at the N-terminus of the residue.
    pub fn n_term_nitrogen_index(&self) -> usize {
        1
    }

    /// Returns the index of the carboxy carbon at the C-terminus of the residue.
    pub fn c_term_carbon_index(&self) -> usize {
        2
    }

    /// Returns the optional index of the start of the amino group in the residue.
    ///
    /// This returns the index of the start of the amino group in the residue, if it exists.
    pub fn sidechain_amino_nitrogen_index(&self) -> Option<usize> {
        match *self {
            Self::Orn => Some(7),
            Self::Dab => Some(6),
            Self::Dap => Some(5),
            _ => None,
        }
    }

    /// Returns the optional index of the start of the carboxyl group in the residue.
    ///
    /// This returns the index of the carboxy carbon in the side chain of the residue, if it
    /// exists.
    pub fn sidechain_carboxy_carbon_index(&self) -> Option<usize> {
        match *self {
            Self::Aad => Some(8),
            _ => None,
        }
    }
}

/// This error represents an error with parsing a standard amino acid.
///
/// The derived traits enable us to set this error message as a source for downstream errors.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ParseNonStdAminoAcidError {
    cause: String,
}

impl ParseNonStdAminoAcidError {
    /// Expose the cause of the error, to use in downstream errors.
    pub fn get_cause(&self) -> String {
        self.cause.clone()
    }
}

impl std::fmt::Display for ParseNonStdAminoAcidError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Invalid amino acid string {}, must be a valid case-insensitive three-letter amino acid code", self.cause)
    }
}

impl std::error::Error for ParseNonStdAminoAcidError {}

impl FromStr for NonStdAminoAcid {
    type Err = ParseNonStdAminoAcidError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_ref() {
            "AAD" => Ok(Self::Aad),
            "ABU" => Ok(Self::Abu),
            "BHT" => Ok(Self::Bht),
            "DPG" => Ok(Self::Dpg),
            "HPG" => Ok(Self::Hpg),
            "ORN" => Ok(Self::Orn),
            "PHG" => Ok(Self::Phg),
            "PIP" => Ok(Self::Pip),
            "CIT" => Ok(Self::Cit),
            "DAB" => Ok(Self::Dab),
            "HTY" => Ok(Self::Hty),
            "HSE" => Ok(Self::Hse),
            "HRN" => Ok(Self::Hrn),
            "HAN" => Ok(Self::Han),
            "DAP" => Ok(Self::Dap),
            "BND" => Ok(Self::Bnd),
            "ADH" => Ok(Self::Adh),
            "DMW" => Ok(Self::Dmw),
            "DHB" => Ok(Self::Dhb),
            "CRP" => Ok(Self::Crp),
            "MHA" => Ok(Self::Mha),
            "HFN" => Ok(Self::Hfn),
            "MPR" => Ok(Self::Mpr),
            "BLA" => Ok(Self::Bla),
            "BBU" => Ok(Self::Bbu),
            "ADE" => Ok(Self::Ade),
            "PYR" => Ok(Self::Pyr),
            "CGY" => Ok(Self::Cgy),
            "OMY" => Ok(Self::Omy),
            "PIZ" => Ok(Self::Piz),
            _ => Err(Self::Err {
                cause: s.to_string(),
            }), // invalid string
        }
    }
}

impl fmt::Display for NonStdAminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Self::Aad => "AAD",
            Self::Abu => "ABU",
            Self::Bht => "BHT",
            Self::Dpg => "DPG",
            Self::Hpg => "HPG",
            Self::Orn => "ORN",
            Self::Phg => "PHG",
            Self::Pip => "PIP",
            Self::Cit => "CIT",
            Self::Dab => "DAB",
            Self::Hty => "HTY",
            Self::Hse => "HSE",
            Self::Hrn => "HRN",
            Self::Han => "HAN",
            Self::Dap => "DAP",
            Self::Bnd => "BND",
            Self::Adh => "ADH",
            Self::Dmw => "DMW",
            Self::Dhb => "DHB",
            Self::Crp => "CRP",
            Self::Mha => "MHA",
            Self::Hfn => "HFN",
            Self::Mpr => "MPR",
            Self::Bla => "BLA",
            Self::Bbu => "BBU",
            Self::Ade => "ADE",
            Self::Pyr => "PYR",
            Self::Cgy => "CGY",
            Self::Omy => "OMY",
            Self::Piz => "PIZ",
        };
        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod test;
