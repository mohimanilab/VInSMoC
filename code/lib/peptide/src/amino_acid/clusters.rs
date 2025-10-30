use super::{nonstandard::NonStdAminoAcid, standard::StdAminoAcid};
use crate::StdOrNonStdAa;
use learning::data::CategoricalFeature;
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use thiserror::Error;

/// Associated error type for parsing an [`StdOrNonStdAa`] into a [`SmallCluster`] or a
/// [`LargeCluster`].
#[derive(Error, Debug, Clone, Eq, PartialEq, Hash)]
pub enum ParseAminoAcidClusterError {
    /// Indicates that the [`StdOrNonStdAa`] did not belong to any [`SmallCluster`].
    #[error("amino acid {0} does not belong to any small cluster amino acids")]
    NoMatchingSmallCluster(StdOrNonStdAa),
    /// Indicates that the [`StdOrNonStdAa`] did not belong to any [`LargeCluster`].
    #[error("amino acid {0} does not belong to any large cluster of amino acids")]
    NoMatchingLargeCluster(StdOrNonStdAa),
}

/// Small clusters of amino acids, corresponding logical groupings of amino acid.
///
/// Similar amino acids are part of the same group for classification. In cases where we don't have
/// enough training data, it becomes difficult to predict top one amino acids. It is a more
/// amenable task to predict the top one small cluster. For more information, please see the
/// NRPSpredictor2 paper: <https://academic.oup.com/nar/article/39/suppl_2/W362/2506164>
#[derive(Copy, Clone, Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub enum SmallCluster {
    /// Singular cluster of 2-amino-adipic acid
    TwoAminoAdipicAcid,
    /// Cluster corresponding to Dhb, Sal
    DhbSal,
    /// Polar, uncharged (hydroxy-phenyl): Dhpg, Hpg
    PolarUncharged,
    /// Singular cluster of cysteine
    CysSpecific,
    /// Singular cluster of serine
    SerineSpecific,
    /// Singular cluster of threonine
    ThreonineSpecific,
    /// Cluster corresponding to Asp, Asn
    AspAsn,
    /// Cluster corresponding to Orn and hydroxy- Orn
    OrnSpecific,
    /// Aliphatic, branched hydrophobic: Val, Leu, Ile, Abu, Iva
    AliphaticBranchedHydrophobic,
    /// Tiny, hydrophilic, transition to aliphatic: Gly, Ala
    TransitionToAliphatic,
    /// Singular cluster of proline
    ProSpecific,
    /// Polar aromatic ring: Tyr, Bht
    PolarAromaticRing,
    /// Cluster corresponding to Glu, Gln
    GluGln,
    /// Singular cluster of arginine
    ArgSpecific,
    /// Unpolar aromatic ring: Phe, Trp
    UnpolarAromaticRing,
    /// beta amino acids
    BetaAmino,
}

impl TryFrom<StdOrNonStdAa> for SmallCluster {
    type Error = ParseAminoAcidClusterError;

    fn try_from(value: StdOrNonStdAa) -> Result<Self, Self::Error> {
        match value {
            StdOrNonStdAa::NonStd(NonStdAminoAcid::Aad) => Ok(SmallCluster::TwoAminoAdipicAcid),
            StdOrNonStdAa::NonStd(NonStdAminoAcid::Hpg) => Ok(SmallCluster::PolarUncharged),
            StdOrNonStdAa::Std(StdAminoAcid::C) => Ok(SmallCluster::CysSpecific),
            StdOrNonStdAa::Std(StdAminoAcid::S) => Ok(SmallCluster::SerineSpecific),
            StdOrNonStdAa::Std(StdAminoAcid::T) => Ok(SmallCluster::ThreonineSpecific),
            StdOrNonStdAa::Std(StdAminoAcid::D) | StdOrNonStdAa::Std(StdAminoAcid::N) => {
                Ok(SmallCluster::AspAsn)
            }
            StdOrNonStdAa::NonStd(NonStdAminoAcid::Orn)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Hrn) => Ok(SmallCluster::OrnSpecific),
            StdOrNonStdAa::Std(StdAminoAcid::V)
            | StdOrNonStdAa::Std(StdAminoAcid::L)
            | StdOrNonStdAa::Std(StdAminoAcid::I)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Abu) => {
                Ok(SmallCluster::AliphaticBranchedHydrophobic)
            }
            StdOrNonStdAa::Std(StdAminoAcid::G) | StdOrNonStdAa::Std(StdAminoAcid::A) => {
                Ok(SmallCluster::TransitionToAliphatic)
            }
            StdOrNonStdAa::Std(StdAminoAcid::P) | StdOrNonStdAa::NonStd(NonStdAminoAcid::Mpr) => {
                Ok(SmallCluster::ProSpecific)
            }
            StdOrNonStdAa::Std(StdAminoAcid::Y) | StdOrNonStdAa::NonStd(NonStdAminoAcid::Bht) => {
                Ok(SmallCluster::PolarAromaticRing)
            }
            StdOrNonStdAa::Std(StdAminoAcid::E) | StdOrNonStdAa::Std(StdAminoAcid::Q) => {
                Ok(SmallCluster::GluGln)
            }
            StdOrNonStdAa::Std(StdAminoAcid::R) => Ok(SmallCluster::ArgSpecific),
            StdOrNonStdAa::Std(StdAminoAcid::F) | StdOrNonStdAa::Std(StdAminoAcid::W) => {
                Ok(SmallCluster::UnpolarAromaticRing)
            }
            StdOrNonStdAa::NonStd(NonStdAminoAcid::Bla)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Bbu) => Ok(SmallCluster::BetaAmino),
            _ => Err(ParseAminoAcidClusterError::NoMatchingSmallCluster(value)),
        }
    }
}

impl SmallCluster {
    /// Mapping from small cluster to its constituent amino acids.
    pub fn to_constituents(&self) -> &'static [StdOrNonStdAa] {
        match *self {
            Self::TwoAminoAdipicAcid => &[StdOrNonStdAa::NonStd(NonStdAminoAcid::Aad)],
            Self::DhbSal => &[],
            Self::PolarUncharged => &[StdOrNonStdAa::NonStd(NonStdAminoAcid::Hpg)],
            Self::CysSpecific => &[StdOrNonStdAa::Std(StdAminoAcid::C)],
            Self::SerineSpecific => &[StdOrNonStdAa::Std(StdAminoAcid::S)],
            Self::ThreonineSpecific => &[StdOrNonStdAa::Std(StdAminoAcid::T)],
            Self::AspAsn => &[
                StdOrNonStdAa::Std(StdAminoAcid::D),
                StdOrNonStdAa::Std(StdAminoAcid::N),
            ],
            Self::OrnSpecific => &[
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Orn),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Hrn),
            ],
            Self::AliphaticBranchedHydrophobic => &[
                StdOrNonStdAa::Std(StdAminoAcid::V),
                StdOrNonStdAa::Std(StdAminoAcid::L),
                StdOrNonStdAa::Std(StdAminoAcid::I),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Abu),
            ],
            Self::TransitionToAliphatic => &[
                StdOrNonStdAa::Std(StdAminoAcid::G),
                StdOrNonStdAa::Std(StdAminoAcid::A),
            ],
            Self::ProSpecific => &[
                StdOrNonStdAa::Std(StdAminoAcid::P),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Mpr),
            ],
            Self::PolarAromaticRing => &[
                StdOrNonStdAa::Std(StdAminoAcid::Y),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Bht),
            ],
            Self::GluGln => &[
                StdOrNonStdAa::Std(StdAminoAcid::E),
                StdOrNonStdAa::Std(StdAminoAcid::Q),
            ],
            Self::ArgSpecific => &[StdOrNonStdAa::Std(StdAminoAcid::R)],
            Self::UnpolarAromaticRing => &[
                StdOrNonStdAa::Std(StdAminoAcid::F),
                StdOrNonStdAa::Std(StdAminoAcid::W),
            ],
            Self::BetaAmino => &[
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Bla),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Bbu),
            ],
        }
    }
}

impl CategoricalFeature for SmallCluster {
    const NUM_CATEGORIES: usize = 16;

    fn index(&self) -> usize {
        match *self {
            Self::TwoAminoAdipicAcid => 0usize,
            Self::DhbSal => 1usize,
            Self::PolarUncharged => 2usize,
            Self::CysSpecific => 3usize,
            Self::SerineSpecific => 4usize,
            Self::ThreonineSpecific => 5usize,
            Self::AspAsn => 6usize,
            Self::OrnSpecific => 7usize,
            Self::AliphaticBranchedHydrophobic => 8usize,
            Self::TransitionToAliphatic => 9usize,
            Self::ProSpecific => 10usize,
            Self::PolarAromaticRing => 11usize,
            Self::GluGln => 12usize,
            Self::ArgSpecific => 13usize,
            Self::UnpolarAromaticRing => 14usize,
            Self::BetaAmino => 15usize,
        }
    }

    fn from_index(index: usize) -> Option<Self> {
        match index {
            0usize => Some(Self::TwoAminoAdipicAcid),
            1usize => Some(Self::DhbSal),
            2usize => Some(Self::PolarUncharged),
            3usize => Some(Self::CysSpecific),
            4usize => Some(Self::SerineSpecific),
            5usize => Some(Self::ThreonineSpecific),
            6usize => Some(Self::AspAsn),
            7usize => Some(Self::OrnSpecific),
            8usize => Some(Self::AliphaticBranchedHydrophobic),
            9usize => Some(Self::TransitionToAliphatic),
            10usize => Some(Self::ProSpecific),
            11usize => Some(Self::PolarAromaticRing),
            12usize => Some(Self::GluGln),
            13usize => Some(Self::ArgSpecific),
            14usize => Some(Self::UnpolarAromaticRing),
            15usize => Some(Self::BetaAmino),
            _ => None,
        }
    }
}

/// Large clusters of amino acids, corresponding logical groupings of amino acid.
///
/// Similar amino acids are part of the same group for classification. In cases where we don't have
/// enough training data, it becomes difficult to predict top one amino acids. It is a more
/// amenable task to predict the top one large cluster. For more information, please see the
/// NRPSpredictor2 paper: <https://academic.oup.com/nar/article/39/suppl_2/W362/2506164>
#[derive(Copy, Clone, Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub enum LargeCluster {
    /// Hydroxy-benzoic acid derivates: Dhb, Sal
    HydroxybenzoicAcidDerivates,
    /// Polar, uncharged (aliphatic with -SH): Cys
    PolarUncharged,
    /// Aliphatic chain or phenyl group with -OH: Ser, Thr, Dhpg, Hpg
    AliphaticChainOrPhenylGroup,
    /// Aliphatic chain with H-bond donor: Asp, Asn, Glu, Gln, Aad
    AliphaticChainWithHbondDonor,
    /// Apolar, aliphatic: Gly, Ala, Val, Leu, Ile, Abu, Iva
    ApolarAliphatic,
    /// Aromatic side chain: Phe, Trp, Phg, Tyr, Bht
    AromaticSideChain,
    /// Cyclic aliphatic chain (polar NH2 group): Pro, Pip
    CyclicAliphaticChain,
    /// Long positively charged side chain: Orn, Lys, Arg, His
    LongPositivelyChargedSideChain,
    /// Beta amino acids
    BetaAmino,
}

impl TryFrom<StdOrNonStdAa> for LargeCluster {
    type Error = ParseAminoAcidClusterError;

    fn try_from(value: StdOrNonStdAa) -> Result<Self, Self::Error> {
        match value {
            StdOrNonStdAa::Std(StdAminoAcid::C) => Ok(LargeCluster::PolarUncharged),
            StdOrNonStdAa::Std(StdAminoAcid::S)
            | StdOrNonStdAa::Std(StdAminoAcid::T)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Hpg) => {
                Ok(LargeCluster::AliphaticChainOrPhenylGroup)
            }
            StdOrNonStdAa::Std(StdAminoAcid::D)
            | StdOrNonStdAa::Std(StdAminoAcid::N)
            | StdOrNonStdAa::Std(StdAminoAcid::E)
            | StdOrNonStdAa::Std(StdAminoAcid::Q)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Aad) => {
                Ok(LargeCluster::AliphaticChainWithHbondDonor)
            }
            StdOrNonStdAa::Std(StdAminoAcid::G)
            | StdOrNonStdAa::Std(StdAminoAcid::A)
            | StdOrNonStdAa::Std(StdAminoAcid::V)
            | StdOrNonStdAa::Std(StdAminoAcid::L)
            | StdOrNonStdAa::Std(StdAminoAcid::I)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Abu) => Ok(LargeCluster::ApolarAliphatic),
            StdOrNonStdAa::Std(StdAminoAcid::F)
            | StdOrNonStdAa::Std(StdAminoAcid::W)
            | StdOrNonStdAa::Std(StdAminoAcid::Y)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Phg)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Bht) => Ok(LargeCluster::AromaticSideChain),
            StdOrNonStdAa::Std(StdAminoAcid::P)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Pip)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Mpr) => Ok(LargeCluster::CyclicAliphaticChain),
            StdOrNonStdAa::Std(StdAminoAcid::K)
            | StdOrNonStdAa::Std(StdAminoAcid::R)
            | StdOrNonStdAa::Std(StdAminoAcid::H)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Orn) => {
                Ok(LargeCluster::LongPositivelyChargedSideChain)
            }
            StdOrNonStdAa::NonStd(NonStdAminoAcid::Bla)
            | StdOrNonStdAa::NonStd(NonStdAminoAcid::Bbu) => Ok(LargeCluster::BetaAmino),
            _ => Err(ParseAminoAcidClusterError::NoMatchingLargeCluster(value)),
        }
    }
}

impl LargeCluster {
    /// Mapping from large cluster to its constituent amino acids.
    pub fn to_constituents(&self) -> &'static [StdOrNonStdAa] {
        match *self {
            Self::HydroxybenzoicAcidDerivates => &[],
            Self::PolarUncharged => &[StdOrNonStdAa::Std(StdAminoAcid::C)],
            Self::AliphaticChainOrPhenylGroup => &[
                StdOrNonStdAa::Std(StdAminoAcid::S),
                StdOrNonStdAa::Std(StdAminoAcid::T),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Hpg),
            ],
            Self::AliphaticChainWithHbondDonor => &[
                StdOrNonStdAa::Std(StdAminoAcid::D),
                StdOrNonStdAa::Std(StdAminoAcid::N),
                StdOrNonStdAa::Std(StdAminoAcid::E),
                StdOrNonStdAa::Std(StdAminoAcid::Q),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Aad),
            ],
            Self::ApolarAliphatic => &[
                StdOrNonStdAa::Std(StdAminoAcid::G),
                StdOrNonStdAa::Std(StdAminoAcid::A),
                StdOrNonStdAa::Std(StdAminoAcid::V),
                StdOrNonStdAa::Std(StdAminoAcid::L),
                StdOrNonStdAa::Std(StdAminoAcid::I),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Abu),
            ],
            Self::AromaticSideChain => &[
                StdOrNonStdAa::Std(StdAminoAcid::F),
                StdOrNonStdAa::Std(StdAminoAcid::W),
                StdOrNonStdAa::Std(StdAminoAcid::Y),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Phg),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Bht),
            ],
            Self::CyclicAliphaticChain => &[
                StdOrNonStdAa::Std(StdAminoAcid::P),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Pip),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Mpr),
            ],
            Self::LongPositivelyChargedSideChain => &[
                StdOrNonStdAa::Std(StdAminoAcid::K),
                StdOrNonStdAa::Std(StdAminoAcid::R),
                StdOrNonStdAa::Std(StdAminoAcid::H),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Orn),
            ],
            Self::BetaAmino => &[
                StdOrNonStdAa::Std(StdAminoAcid::G),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Bla),
                StdOrNonStdAa::NonStd(NonStdAminoAcid::Bbu),
            ],
        }
    }
}

impl CategoricalFeature for LargeCluster {
    const NUM_CATEGORIES: usize = 9;

    fn index(&self) -> usize {
        match *self {
            Self::HydroxybenzoicAcidDerivates => 0usize,
            Self::PolarUncharged => 1usize,
            Self::AliphaticChainOrPhenylGroup => 2usize,
            Self::AliphaticChainWithHbondDonor => 3usize,
            Self::ApolarAliphatic => 4usize,
            Self::AromaticSideChain => 5usize,
            Self::CyclicAliphaticChain => 6usize,
            Self::LongPositivelyChargedSideChain => 7usize,
            Self::BetaAmino => 8usize,
        }
    }

    fn from_index(index: usize) -> Option<Self> {
        match index {
            0usize => Some(Self::HydroxybenzoicAcidDerivates),
            1usize => Some(Self::PolarUncharged),
            2usize => Some(Self::AliphaticChainOrPhenylGroup),
            3usize => Some(Self::AliphaticChainWithHbondDonor),
            4usize => Some(Self::ApolarAliphatic),
            5usize => Some(Self::AromaticSideChain),
            6usize => Some(Self::CyclicAliphaticChain),
            7usize => Some(Self::LongPositivelyChargedSideChain),
            8usize => Some(Self::BetaAmino),
            _ => None,
        }
    }
}
