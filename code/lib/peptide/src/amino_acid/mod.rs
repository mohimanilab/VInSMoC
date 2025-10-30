use std::{fmt, str::FromStr};

use constants::amino_acids::NUM_BACKBONE_ATOMS;
use learning::data::CategoricalFeature;
use modification::ModSeq;
use molecule::{Atom, BondType, Mol, MolBuilder, MolGraph};
use petgraph::stable_graph::node_index;
use serde::{Deserialize, Serialize};

mod standard;
pub use standard::{ParseStdAminoAcidError, StdAminoAcid, NUM_STANDARD_AMINO_ACIDS};
mod nonstandard;
pub use nonstandard::{NonStdAminoAcid, ParseNonStdAminoAcidError};
mod clusters;
pub use clusters::{LargeCluster, SmallCluster};

/// The residual backbone of an amino acid.
///
/// This does not contain the N-terminal hydrogen or the C-terminal hydroxyl group that would be
/// expected for a standalone amino acid molecule. Instead this is just the residual backbone.
///
/// # Indices
/// The indices of the backbone atoms in the underlying molecular graph are fixed.
///
/// | Atom | Index |
/// | ---- | ----- |
/// | α-Carbon | 0 |
/// | Nitrogen | 1 |
/// | Carbonyl Carbon | 2 |
/// | Carbonyl Oxygen | 3 |
pub fn backbone() -> Mol {
    let mut backbone = MolGraph::with_capacity(4, 3);
    backbone.add_node(Atom::C);
    backbone.add_node(Atom::N);
    backbone.add_node(Atom::C);
    backbone.add_node(Atom::O);

    backbone.add_edge(node_index(0), node_index(1), BondType::Single);
    backbone.add_edge(node_index(0), node_index(2), BondType::Single);
    backbone.add_edge(node_index(2), node_index(3), BondType::Double);

    let hydrogen_hash = [(node_index(0), 1), (node_index(1), 1)]
        .iter()
        .copied()
        .collect();

    MolBuilder::default()
        .graph(backbone)
        .hydrogens(hydrogen_hash)
        .build()
}

/// A general purpose amino acid.
///
/// We define a general amino acid as any molecule that contains the peptidic backbone and branches
/// off of it with some R-group. In order to define a new amino acid we require the implementation
/// of the structure of the R-group and a description of how the R-group connects to the static
/// backbone.
///
/// This trait automatically provides the ability to convert the defined amino acid into a residue
/// molecule ([`AminoAcid::to_residue_mol`]) and calculate residual ([`AminoAcid::residue_std_mass`], [`AminoAcid::residue_exact_mass`])
/// and full amino acid ([`AminoAcid::std_mass`], [`AminoAcid::exact_mass`]) masses.
///
/// # Examples
/// We will implement this trait for a fake amino acid whose R-group is just a chlorine molecule.
/// ```
/// use peptide::AminoAcid;
/// use molecule::{Mol, Atom, BondType};
/// use constants::amino_acids::NUM_BACKBONE_ATOMS;
/// use modification::{Mod, ModSeq};
/// use petgraph::graph::node_index;
///
/// struct ClAa {}
///
/// impl AminoAcid for ClAa {
///     fn r_group(&self) -> Mol {
///         // start with an empty molecule
///         let mut mol = Mol::default();
///         // add in our chlorine atom
///         let cl_idx = mol.graph.add_node(Atom::Cl);
///
///         mol
///     }
///
///     fn connection_mod(&self) -> ModSeq {
///         let mods = vec![
///             Mod::AddEdge(
///                 0,                          // connect alpha carbon
///                 NUM_BACKBONE_ATOMS as i32,  // to the chlorine
///                 BondType::Single            // with a single bond
///             )
///         ];
///         ModSeq::new(mods)
///     }
/// }
/// ```
pub trait AminoAcid {
    /// Returns a molecule representing the side chain of this amino acid.
    ///
    /// To get a full amino acid residue this R-group is combined with the amino acid [`backbone`]
    /// according to the modifications specified in [`AminoAcid::connection_mod`].
    fn r_group(&self) -> Mol;

    /// Returns a modification sequence that should be used to connect the backbone to the R-group
    /// for this amino acid.
    ///
    /// The indices this modification sequence operates on are indices of a graph with 2 connected
    /// components construct by first inserting backbone atoms into a graph and then adding R-group
    /// atoms to the same graph. All of the R-group atom indices will be shifted by
    /// [`constants::amino_acids::NUM_BACKBONE_ATOMS`], which should be taken into account when implementing
    /// this function. If an empty modification sequence is returned the final amino acid will have
    /// 2 connected components.
    fn connection_mod(&self) -> ModSeq;

    /// Returns the residual amino acid molecule.
    ///
    /// This connects a backbone to the R-group associated with this amino acid. The
    /// [`AminoAcid::connection_mod`] function is used to apply the actual connection between the backbone and
    /// the R-group. The backbone will always be the first indices in the resulting residue
    /// molecule, with the 0-index R-group atom appearing in the residue molecule at index 4 (after
    /// the backbone atoms). See also [`backbone`].
    ///
    /// This will panic if the [`Self::connection_mod`] applied to the [`Self::r_group`] and [`backbone`]
    /// fails. See [`modification::ModificationError`] for information about the possible fail
    /// points.
    fn to_residue_mol(&self) -> Mol {
        let mut residue: Mol = backbone();

        // combine backbone and R-group graphs, shift indices
        // this should have 2 connnected components: 1 for the backbone and 1 for the R-group
        let r_group: Mol = self.r_group();
        for node in r_group.graph.node_indices() {
            residue.graph.add_node(r_group.graph[node]);
        }
        for edge in r_group.graph.edge_indices() {
            let (source, target) = r_group.graph.edge_endpoints(edge).unwrap();
            residue.graph.add_edge(
                node_index(source.index() + NUM_BACKBONE_ATOMS),
                node_index(target.index() + NUM_BACKBONE_ATOMS),
                r_group.graph[edge],
            );
        }

        // combine hydrogens field
        residue.extend_hydrogens(
            r_group
                .hydrogens()
                .iter()
                .map(|(node, count)| (node_index(node.index() + NUM_BACKBONE_ATOMS), *count)),
        );

        // combine charges field
        residue.extend_charges(
            r_group
                .charges()
                .iter()
                .map(|(node, charge)| (node_index(node.index() + NUM_BACKBONE_ATOMS), *charge)),
        );

        // use this AA's connection modification to merge 2-cc graph into one residue graph
        // this will panic if you defined the modification illogically
        self.connection_mod()
            .modify(&mut residue)
            .expect("Provided connection modification failed");

        residue.compact();

        residue
    }

    /// Returns the standard mass of the residual amino acid.
    ///
    /// Internally this uses [`AminoAcid::to_residue_mol`], so if you need to convert the amino
    /// acid into a molecule already prefer using [`molecule::Mol::std_mass`].
    fn residue_std_mass(&self) -> f64 {
        self.to_residue_mol().std_mass()
    }

    /// Returns the exact mass of the residual amino acid.
    ///
    /// Internally this uses [`AminoAcid::to_residue_mol`], so if you need to convert the amino
    /// acid into a molecule already prefer using [`molecule::Mol::exact_mass`].
    fn residue_exact_mass(&self) -> f64 {
        self.to_residue_mol().exact_mass()
    }

    /// Returns the standard mass of the whole amino acid
    ///
    /// This is just a sum of the residue mass and the mass of the missing water molecule.
    /// Internally this uses [`AminoAcid::to_residue_mol`], so if you need to convert the amino
    /// acid into a molecule already prefer using [`molecule::Mol::std_mass`] adding the mass of
    /// a water molecule.
    fn std_mass(&self) -> f64 {
        self.to_residue_mol().std_mass() + constants::MASS_WATER
    }

    /// Returns the exact mass of the whole amino acid
    ///
    /// This is just a sum of the residue mass and the mass of the missing water molecule.
    /// Internally this uses [`AminoAcid::to_residue_mol`], so if you need to convert the amino
    /// acid into a molecule already prefer using [`molecule::Mol::exact_mass`] adding the mass of
    /// a water molecule.
    fn exact_mass(&self) -> f64 {
        self.to_residue_mol().exact_mass() + constants::MASS_WATER_EXACT
    }
}

/// An enum to represent either a [`StdAminoAcid`] or [`NonStdAminoAcid`].
#[derive(Debug, Eq, PartialEq, Clone, Copy, Hash)]
pub enum StdOrNonStdAa {
    /// This variant represents a standard amino acid, as defined by the variants of
    /// [`StdAminoAcid`].
    ///
    /// There will always be [`NUM_STANDARD_AMINO_ACIDS`] possible variants.
    Std(StdAminoAcid),
    /// This variant represents a non-standard amino acid, as defined by the variants of
    /// [`NonStdAminoAcid`].
    NonStd(NonStdAminoAcid),
}

impl CategoricalFeature for StdOrNonStdAa {
    const NUM_CATEGORIES: usize = StdAminoAcid::NUM_CATEGORIES + NonStdAminoAcid::NUM_CATEGORIES;

    #[inline]
    fn index(&self) -> usize {
        match *self {
            Self::Std(x) => x.index(),
            Self::NonStd(x) => NUM_STANDARD_AMINO_ACIDS + x.index(),
        }
    }

    #[inline]
    fn from_index(index: usize) -> Option<Self> {
        if index < NUM_STANDARD_AMINO_ACIDS {
            StdAminoAcid::from_index(index).map(Self::Std)
        } else {
            NonStdAminoAcid::from_index(index - NUM_STANDARD_AMINO_ACIDS).map(Self::NonStd)
        }
    }
}

impl AminoAcid for StdOrNonStdAa {
    fn r_group(&self) -> Mol {
        match *self {
            Self::Std(x) => x.r_group(),
            Self::NonStd(x) => x.r_group(),
        }
    }

    fn connection_mod(&self) -> ModSeq {
        match *self {
            Self::Std(x) => x.connection_mod(),
            Self::NonStd(x) => x.connection_mod(),
        }
    }
}

impl StdOrNonStdAa {
    /// Returns a vector of all variants of this enum for use with sequence tags.
    pub fn get_all_aas() -> Vec<Self> {
        let mut aas = Vec::new();

        // Add StdAminoAcid variants
        aas.extend([
            Self::Std(StdAminoAcid::G),
            Self::Std(StdAminoAcid::R),
            Self::Std(StdAminoAcid::H),
            Self::Std(StdAminoAcid::K),
            Self::Std(StdAminoAcid::D),
            Self::Std(StdAminoAcid::E),
            Self::Std(StdAminoAcid::S),
            Self::Std(StdAminoAcid::T),
            Self::Std(StdAminoAcid::N),
            Self::Std(StdAminoAcid::Q),
            Self::Std(StdAminoAcid::C),
            Self::Std(StdAminoAcid::P),
            Self::Std(StdAminoAcid::A),
            Self::Std(StdAminoAcid::V),
            Self::Std(StdAminoAcid::I),
            Self::Std(StdAminoAcid::L),
            Self::Std(StdAminoAcid::M),
            Self::Std(StdAminoAcid::F),
            Self::Std(StdAminoAcid::Y),
            Self::Std(StdAminoAcid::W),
        ]);

        // Add NonStdAminoAcid variants
        aas.extend([
            Self::NonStd(NonStdAminoAcid::Aad),
            Self::NonStd(NonStdAminoAcid::Abu),
            Self::NonStd(NonStdAminoAcid::Bht),
            Self::NonStd(NonStdAminoAcid::Dpg),
            Self::NonStd(NonStdAminoAcid::Hpg),
            Self::NonStd(NonStdAminoAcid::Orn),
            Self::NonStd(NonStdAminoAcid::Phg),
            Self::NonStd(NonStdAminoAcid::Pip),
            Self::NonStd(NonStdAminoAcid::Cit),
            Self::NonStd(NonStdAminoAcid::Dab),
            Self::NonStd(NonStdAminoAcid::Hty),
            Self::NonStd(NonStdAminoAcid::Hse),
            Self::NonStd(NonStdAminoAcid::Hrn),
            Self::NonStd(NonStdAminoAcid::Han),
            Self::NonStd(NonStdAminoAcid::Dap),
            Self::NonStd(NonStdAminoAcid::Bnd),
            Self::NonStd(NonStdAminoAcid::Adh),
            Self::NonStd(NonStdAminoAcid::Dmw),
            Self::NonStd(NonStdAminoAcid::Dhb),
            Self::NonStd(NonStdAminoAcid::Crp),
            Self::NonStd(NonStdAminoAcid::Mha),
            Self::NonStd(NonStdAminoAcid::Hfn),
            Self::NonStd(NonStdAminoAcid::Mpr),
            Self::NonStd(NonStdAminoAcid::Bla),
            Self::NonStd(NonStdAminoAcid::Bbu),
            Self::NonStd(NonStdAminoAcid::Ade),
            Self::NonStd(NonStdAminoAcid::Pyr),
            Self::NonStd(NonStdAminoAcid::Cgy),
            Self::NonStd(NonStdAminoAcid::Omy),
            Self::NonStd(NonStdAminoAcid::Piz),
        ]);

        aas
    }

    /// Gets the underlying standard amino acid, if it exists.
    pub fn as_std(&self) -> Option<StdAminoAcid> {
        if let Self::Std(v) = self {
            Some(*v)
        } else {
            None
        }
    }

    /// Gets the underlying non-standard amino acid, if it exists.
    pub fn as_non_std(&self) -> Option<NonStdAminoAcid> {
        if let Self::NonStd(v) = self {
            Some(*v)
        } else {
            None
        }
    }

    /// Return the [`SmallCluster`] this [`StdOrNonStdAa`] corresponds to, if any.
    pub fn to_small_cluster(self) -> Option<SmallCluster> {
        self.try_into().ok()
    }

    /// Return the [`LargeCluster`] this [`StdOrNonStdAa`] corresponds to, if any.
    pub fn to_large_cluster(self) -> Option<LargeCluster> {
        self.try_into().ok()
    }

    /// Returns the index of the alpha carbon of the residue.
    pub fn alpha_c_idx(&self) -> usize {
        0
    }

    /// Returns the index of the nitrogen at the N-terminus of the residue.
    pub fn n_term_idx(&self) -> usize {
        1
    }

    /// Returns the index of the carbon at the carboxyl group (C-terminus) of the residue.
    pub fn carbonyl_c_idx(&self) -> usize {
        2
    }

    /// Returns the index of the oxygen at the carboxyl group (C-terminus) of the residue.
    pub fn carbonyl_o_idx(&self) -> usize {
        3
    }

    /// Returns the index of the alcohol group in the side chain of the residue, if it exists.
    pub fn sidechain_alcohol_index(&self) -> Option<usize> {
        match *self {
            Self::Std(x) => x.sidechain_alcohol_index(),
            Self::NonStd(_) => None,
        }
    }

    /// Returns the index of the thiol group in the side chain of the residue, if it exists.
    pub fn sidechain_thiol_index(&self) -> Option<usize> {
        match *self {
            Self::Std(x) => x.sidechain_thiol_index(),
            Self::NonStd(_) => None,
        }
    }

    /// Returns the optional index of the start of the amino group in the residue.
    ///
    /// This returns the index of the start of the amino group in the residue, if it exists.
    pub fn sidechain_amino_nitrogen_index(&self) -> Option<usize> {
        match *self {
            Self::Std(x) => x.sidechain_amino_nitrogen_index(),
            Self::NonStd(x) => x.sidechain_amino_nitrogen_index(),
        }
    }

    /// Returns the index of the carbonyl carbon in the side chain of the residue, if it exists.
    pub fn sidechain_carboxy_carbon_index(&self) -> Option<usize> {
        match *self {
            Self::Std(x) => x.sidechain_carboxy_carbon_index(),
            Self::NonStd(x) => x.sidechain_carboxy_carbon_index(),
        }
    }

    /// Returns the index of the carbonyl oxygen in the side chain of the residue, if it exists.
    pub fn sidechain_carboxy_oxygen_index(&self) -> Option<usize> {
        self.sidechain_carboxy_carbon_index().map(|x| x + 1)
    }
}

impl Serialize for StdOrNonStdAa {
    fn serialize<S>(&self, ser: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let s = &self.to_string();
        ser.serialize_str(s)
    }
}

impl<'de> Deserialize<'de> for StdOrNonStdAa {
    fn deserialize<D>(de: D) -> Result<StdOrNonStdAa, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{Error, Unexpected};
        let aa = <&str as Deserialize>::deserialize(de)?;
        let aa: StdOrNonStdAa = aa.parse().map_err(|_| {
            D::Error::invalid_value(Unexpected::Str(aa), &"a standard or nonstandard amino acid")
        })?;
        Ok(aa)
    }
}

impl From<ParseStdAminoAcidError> for ParseAminoAcidError {
    fn from(err: ParseStdAminoAcidError) -> Self {
        ParseAminoAcidError {
            cause: err.get_cause(),
        }
    }
}

impl From<ParseNonStdAminoAcidError> for ParseAminoAcidError {
    fn from(err: ParseNonStdAminoAcidError) -> Self {
        ParseAminoAcidError {
            cause: err.get_cause(),
        }
    }
}

impl From<StdAminoAcid> for StdOrNonStdAa {
    fn from(aa: StdAminoAcid) -> Self {
        Self::Std(aa)
    }
}

impl From<NonStdAminoAcid> for StdOrNonStdAa {
    fn from(aa: NonStdAminoAcid) -> Self {
        Self::NonStd(aa)
    }
}

impl fmt::Display for StdOrNonStdAa {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Std(aa) => write!(f, "{}", aa),
            Self::NonStd(aa) => write!(f, "{}", aa),
        }
    }
}

/// An error which is associated with trying to parse invalid amino acids.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct ParseAminoAcidError {
    cause: String,
}

impl std::fmt::Display for ParseAminoAcidError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Invalid amino acid string {}, must be a valid case-insensitive single- or three-letter amino acid code", self.cause)
    }
}

impl std::error::Error for ParseAminoAcidError {}

impl FromStr for StdOrNonStdAa {
    type Err = ParseAminoAcidError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse::<StdAminoAcid>()
            .map(|x| x.into())
            .or_else(|_| s.parse::<NonStdAminoAcid>().map(|x| x.into()))
            .map_err(|_| ParseAminoAcidError {
                cause: s.to_string(),
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_str() {
        assert_eq!(
            "k".parse::<StdOrNonStdAa>(),
            Ok(StdOrNonStdAa::Std(StdAminoAcid::K))
        );
        assert_eq!(
            "P".parse::<StdOrNonStdAa>(),
            Ok(StdOrNonStdAa::Std(StdAminoAcid::P))
        );
        assert_eq!(
            "x".parse::<StdOrNonStdAa>(),
            Err(ParseAminoAcidError {
                cause: String::from("x")
            })
        );
        assert_eq!(
            "aAd".parse::<StdOrNonStdAa>(),
            Ok(StdOrNonStdAa::NonStd(NonStdAminoAcid::Aad))
        );
        assert_eq!(
            "DAB".parse::<StdOrNonStdAa>(),
            Ok(StdOrNonStdAa::NonStd(NonStdAminoAcid::Dab))
        );
        assert_eq!(
            "hse".parse::<StdOrNonStdAa>(),
            Ok(StdOrNonStdAa::NonStd(NonStdAminoAcid::Hse))
        );
        assert_eq!(
            "zzz".parse::<StdOrNonStdAa>(),
            Err(ParseAminoAcidError {
                cause: String::from("zzz")
            })
        );
    }
}
