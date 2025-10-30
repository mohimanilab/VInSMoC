use indexmap::IndexSet;
use itertools::Itertools;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use regex::Regex;
use rustc_hash::FxHashMap;

use super::{Mol, MolGraph, NodeIndex, ValenceError};
use crate::{
    parsers::{parse_mol_file, parse_smiles, MolFileParseError, SmilesParseError},
    Atom,
};
use std::{fs, io, path::Path};

/// An error that can occur when using [`MolBuilder`].
#[derive(Debug, thiserror::Error)]
pub enum BuildMolError {
    /// Error when parsing a MOL file.
    #[error(transparent)]
    MolFile(#[from] MolFileParseError),
    /// Error when parsing a SMILES.
    #[error(transparent)]
    Smiles(#[from] SmilesParseError),
    /// Error during hydrogen inference.
    #[error(transparent)]
    Valence(#[from] ValenceError),
    /// Error performing I/O.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Indexed SMILES had two of the same index
    #[error("duplicate atom index {0}")]
    DuplicateIndex(usize),
    /// Indexed SMILES is missing index annotations for some atoms
    #[error("indexed SMILES index annotations don't match atom count")]
    MissingIndex,
}

/// Builds a [`Mol`].
///
/// This can be more ergonomic than trying to use standalone parse functions and the various
/// hydrogen inference member methods of [`Mol`]. Right now this supports SMILES and V3000 MDL MOL
/// files.
///
/// # Examples
/// ```
/// # use molecule::MolBuilder;
/// # use std::error::Error;
/// # fn t() -> Result<(), Box<dyn Error>> {
/// // build from SMILES with no changes
/// let mol = MolBuilder::from_smiles("CCC")?.build();
/// // make all the SMILES hydrogens explicit
/// let mol = MolBuilder::from_smiles("CCC")?.explicit_hs().build();
/// // infer all hydrogens then make them explicit
/// let mol = MolBuilder::from_smiles("CCC")?.infer_hs()?.explicit_hs().build();
/// # Ok(())
/// # }
/// # assert!(t().is_ok());
/// ```
#[derive(Debug, Clone, Default)]
#[repr(transparent)]
pub struct MolBuilder(Mol);

impl MolBuilder {
    /// Creates a [`MolBuilder`] from a SMILES string.
    ///
    /// By default a [`Mol`] built using SMILES will have implicit hydrogens correctly inferred for
    /// all atoms not surround in square brackets. Unless the SMILES itself contains explicit
    /// hydrogens then the output molecule will not contain explicit hydrogens. This should not
    /// require any hydrogen inference either, since the SMILES format provides mechanisms for
    /// hydrogen inference when desired.
    ///
    /// See the [SMILES specification](https://daylight.com/dayhtml/doc/theory/theory.smiles.html)
    /// for a detailed breakdown of SMILES features. We completely ignore 3D information including
    /// cis/trans orientations (Section 3.3.2) and chiral centers (Section 3.3.3).
    ///
    /// # Known Issues
    /// 1. We do not support isotopic specification (Section 3.3.1)
    /// 2. Some of our aromatic bond inference is incorrect
    pub fn from_smiles(smiles: &str) -> Result<Self, BuildMolError> {
        Ok(parse_smiles(smiles).map(Self)?)
    }

    /// Creates a [`MolBuilder`] from a V3000 MDL MOL file.
    ///
    /// By default a [`Mol`] built from a MOL file will have explicit hydrogens unless the `HCOUNT`
    /// field is used by the input MOL file. In practice this is uncommon. Therefore, most of the
    /// time this will require hydrogen inference to convert standalone hydrogen nodes to the
    /// implicit format.
    ///
    /// See the [MOL file specification](http://c4.cabrillo.edu/404/ctfile.pdf) for detailed
    /// breakdown of V3000 MOL files. The only extra fields we use are `HCOUNT` and `CHG`, the
    /// rest are just ignored.
    pub fn from_mol_file<P: AsRef<Path>>(mol_file: P) -> Result<Self, BuildMolError> {
        let mol_file = fs::read_to_string(mol_file)?;
        let mol = parse_mol_file(&mol_file)?;
        Ok(Self(mol))
    }

    /// Parse a smiles with additional index markers.
    ///
    /// Indexing is done via a `:X` suffix, where `X` is the desired index of the atom in the final
    /// molecule. Additionally, to distinguish between aromatic bonds to a ring, which may take the
    /// form of `:X`, all atoms in an indexed smiles must be surround by `{}`.
    ///
    /// Note that this does not require a call to [`MolBuilder::build`], since it immediately
    /// returns a [`Mol`]. This is by design to prevent callers from easily performing
    /// modifications that change the indices of atoms in the underlying molecule.
    ///
    /// # Examples
    /// In the simple case we can just add indices to all atoms that appear in the SMILES.
    /// ```
    /// use molecule::{MolBuilder, NodeIndex, Atom};
    ///
    /// let og = MolBuilder::from_smiles(r"C1(=O)CC1").unwrap().build();
    /// let indexed = MolBuilder::from_indexed_smiles(r"{C:0}1(={O:1}){C:2}{C:3}1").unwrap();
    /// assert_eq!(og.graph[NodeIndex::new(0)], Atom::C);
    /// assert_eq!(og.graph[NodeIndex::new(1)], Atom::O);
    /// assert_eq!(indexed.graph[NodeIndex::new(0)], Atom::C);
    /// assert_eq!(indexed.graph[NodeIndex::new(1)], Atom::O);
    /// ```
    /// However we can add any indices we like.
    /// ```
    /// use molecule::{MolBuilder, NodeIndex, Atom};
    ///
    /// let og = MolBuilder::from_smiles(r"C1(=O)CC1").unwrap().build();
    /// let indexed = MolBuilder::from_indexed_smiles(r"{C:1000}1(={O:2}){C:4}{C:2000}1").unwrap();
    /// assert_eq!(og.graph[NodeIndex::new(0)], Atom::C);
    /// assert_eq!(og.graph[NodeIndex::new(1)], Atom::O);
    /// assert_eq!(indexed.graph[NodeIndex::new(1000)], Atom::C);
    /// assert_eq!(indexed.graph[NodeIndex::new(2)], Atom::O);
    /// ```
    /// One thing to note is that all atoms must have indices
    /// ```
    /// use molecule::{MolBuilder, BuildMolError};
    ///
    /// // first carbon has no index
    /// let indexed = MolBuilder::from_indexed_smiles(r"C1(={O:2}){C:4}{C:2000}1");
    /// assert!(indexed.is_err());
    /// assert!(matches!(indexed.unwrap_err(), BuildMolError::MissingIndex));
    /// ```
    /// and no indices can be repeated.
    /// ```
    /// use molecule::{MolBuilder, BuildMolError, NodeIndex, Atom};
    ///
    /// // all indices are 0
    /// let indexed = MolBuilder::from_indexed_smiles(r"{C:0}(={O:0}){C:0}{C:0}1");
    /// assert!(indexed.is_err());
    /// assert!(matches!(indexed.unwrap_err(), BuildMolError::DuplicateIndex(0)));
    ///
    /// // index 1 is repeated
    /// let indexed = MolBuilder::from_indexed_smiles(r"{C:1}(={O:0}){C:1}{C:2}1");
    /// assert!(indexed.is_err());
    /// assert!(matches!(indexed.unwrap_err(), BuildMolError::DuplicateIndex(1)));
    /// ```
    /// When dealing with molecules with charged atoms and/or non-standard valences the hydrogens
    /// must be made their own explicit nodes.
    /// ```
    /// use molecule::{MolBuilder, NodeIndex, Atom};
    ///
    /// // these are exactly the same molecule
    /// let og = MolBuilder::from_smiles(r"C[NH+](O)C").unwrap().build();
    /// let indexed = MolBuilder::from_indexed_smiles(r"{C:0}[{N:1}+]([{H:2}])({O:3}){C:4}").unwrap();
    /// assert_eq!(og.graph[NodeIndex::new(0)], Atom::C);
    /// assert_eq!(og.graph[NodeIndex::new(1)], Atom::N);
    /// assert_eq!(og.charge(&NodeIndex::new(1)), 1);
    /// assert_eq!(og.graph[NodeIndex::new(2)], Atom::H);
    /// assert_eq!(indexed.graph[NodeIndex::new(0)], Atom::C);
    /// assert_eq!(indexed.graph[NodeIndex::new(1)], Atom::N);
    /// assert_eq!(indexed.charge(&NodeIndex::new(1)), 1);
    /// assert_eq!(indexed.graph[NodeIndex::new(2)], Atom::H);
    /// ```
    pub fn from_indexed_smiles(smiles: &str) -> Result<Mol, BuildMolError> {
        let mut annotated_idxs = IndexSet::<usize>::default();

        // this regex captures all bracketed atoms as two groups:
        // 1. the parts preceding the index annotation
        // 2. the index annotation itself
        // the first capture group is used to replace matches
        //
        // .unwrap is safe here since this is a valid regex
        let re = Regex::new(r"\{([^\[:]+):(\d+)}").unwrap();

        // track all of the annotated indices
        // since the SMILES parser guarantees atoms get indices in the same order they appear in
        // the SMILES the mapping simply counts upward to track current atom index
        for cap in re.captures_iter(smiles) {
            // .unwrap() is safe because the regex guarantees a valid usize
            let idx = cap[2].parse::<usize>().unwrap();

            // throw error when a duplicate index is found
            if !annotated_idxs.insert(idx) {
                return Err(BuildMolError::DuplicateIndex(idx));
            }
        }

        // delete the index mapping part
        let smiles = re.replace_all(smiles, "$1");
        let og_mol = Self::from_smiles(&smiles)?.build();
        // terminate if they didn't annotate every atom
        if og_mol.graph.node_count() != annotated_idxs.len() {
            return Err(BuildMolError::MissingIndex);
        }

        let mut indexed_mol = Mol {
            graph: MolGraph::with_capacity(og_mol.graph.node_count(), og_mol.graph.edge_count()),
            ..Mol::default()
        };

        let sorted = annotated_idxs
            .iter()
            .copied()
            .enumerate()
            .sorted_by_key(|(_, new_idx)| *new_idx);
        let mut to_remove =
            Vec::with_capacity(*annotated_idxs.iter().max().unwrap() + 1 - annotated_idxs.len());
        for (old_idx, new_idx) in sorted {
            // hammer indices until we're at the one we want
            // we have to delete placeholders later
            while indexed_mol.graph.node_count() != new_idx {
                to_remove.push(indexed_mol.graph.add_node(Atom::H));
            }

            // copy over node data
            let old_idx = NodeIndex::new(old_idx);
            let new_idx = indexed_mol.graph.add_node(og_mol.graph[old_idx]);
            *indexed_mol.h_entry(new_idx) = og_mol.h_count(&old_idx);
            *indexed_mol.charge_entry(new_idx) = og_mol.charge(&old_idx);
            indexed_mol.set_aromaticity(&new_idx, og_mol.is_aromatic(old_idx));
        }

        // delete placeholder nodes
        for idx in to_remove {
            indexed_mol.graph.remove_node(idx);
        }

        // copy over edge data
        for e in og_mol.graph.edge_references() {
            let source = NodeIndex::new(annotated_idxs[e.source().index()]);
            let target = NodeIndex::new(annotated_idxs[e.target().index()]);
            indexed_mol.graph.add_edge(source, target, *e.weight());
        }

        indexed_mol.compact();
        Ok(indexed_mol)
    }

    /// Use the provided molecular graph for this molecule.
    pub fn graph(mut self, graph: MolGraph) -> Self {
        self.0.graph = graph;
        self
    }

    /// Use the provided hydrogen counts for this molecule.
    pub fn hydrogens(mut self, hydrogens: FxHashMap<NodeIndex, u8>) -> Self {
        self.0.hydrogens = hydrogens;
        self
    }

    /// Use the provided charge counts for this molecule.
    pub fn charges(mut self, charges: FxHashMap<NodeIndex, i8>) -> Self {
        self.0.charges = charges;
        self
    }

    /// Use explicit hydrogens for this molecule.
    ///
    /// See [`Mol::explicit_hs`].
    pub fn explicit_hs(mut self) -> Self {
        self.0.explicit_hs();
        self
    }

    /// Use implicit hydrogens for this molecule.
    ///
    /// See [`Mol::implicit_hs`].
    pub fn implicit_hs(mut self) -> Self {
        self.0.implicit_hs();
        self
    }

    /// Infer hydrogens for this molecule.
    ///
    /// See [`Mol::infer_hs`].
    pub fn infer_hs(mut self) -> Result<Self, BuildMolError> {
        Ok(self.0.infer_hs().map(|_| self)?)
    }

    /// Finalize the molecule.
    pub fn build(self) -> Mol {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use rustc_hash::FxHashSet;

    use crate::{BondType, EdgeIndex};

    use super::*;

    #[test]
    fn indexed_simple() {
        let mol = MolBuilder::from_indexed_smiles("{C:0}{N:1}").unwrap();
        assert_eq!(mol.graph.node_count(), 2);
        assert_eq!(mol.graph.edge_count(), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 3);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 2);
        assert!(mol.charges.is_empty(), "{:?}", mol.charges);
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::N);
        assert_eq!(mol.graph[EdgeIndex::new(0)], BondType::Single);

        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("CN")
                .unwrap()
                .build()
                .to_inchi_key()
        );

        let mol = MolBuilder::from_indexed_smiles("{C:1}{N:0}").unwrap();
        assert_eq!(mol.graph.node_count(), 2);
        assert_eq!(mol.graph.edge_count(), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 2);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 3);
        assert!(mol.charges.is_empty());
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::N);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::C);
        assert_eq!(mol.graph[EdgeIndex::new(0)], BondType::Single);

        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("CN")
                .unwrap()
                .build()
                .to_inchi_key()
        );

        let mol = MolBuilder::from_indexed_smiles("{C:1}{N:1000}").unwrap();
        assert_eq!(mol.graph.node_count(), 2);
        assert_eq!(mol.graph.edge_count(), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(1000)), 2);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 3);
        assert!(mol.charges.is_empty());
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(1000)], Atom::N);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::C);
        assert_eq!(mol.graph[EdgeIndex::new(0)], BondType::Single);

        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("CN")
                .unwrap()
                .build()
                .to_inchi_key()
        );
    }

    #[test]
    fn indexed_ring() {
        let mol = MolBuilder::from_indexed_smiles("{C:0}1{N:1}{O:2}1").unwrap();
        assert_eq!(mol.graph.node_count(), 3);
        assert_eq!(mol.graph.edge_count(), 3);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 2);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(2)), 0);
        assert!(mol.charges.is_empty());
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::N);
        assert_eq!(mol.graph[NodeIndex::new(2)], Atom::O);
        assert_eq!(
            mol.graph
                .edge_weights()
                .copied()
                .collect::<FxHashSet<_>>()
                .into_iter()
                .collect_vec(),
            vec![BondType::Single]
        );
        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("C1NO1")
                .unwrap()
                .build()
                .to_inchi_key()
        );

        let mol = MolBuilder::from_indexed_smiles("{C:2}1{N:0}{O:1}1").unwrap();
        assert_eq!(mol.graph.node_count(), 3);
        assert_eq!(mol.graph.edge_count(), 3);
        assert_eq!(mol.h_count(&NodeIndex::new(2)), 2);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 0);
        assert!(mol.charges.is_empty());
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::N);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::O);
        assert_eq!(mol.graph[NodeIndex::new(2)], Atom::C);
        assert_eq!(
            mol.graph
                .edge_weights()
                .copied()
                .collect::<FxHashSet<_>>()
                .into_iter()
                .collect_vec(),
            vec![BondType::Single]
        );
        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("C1NO1")
                .unwrap()
                .build()
                .to_inchi_key()
        );

        let mol = MolBuilder::from_indexed_smiles("{C:0}1(={O:1}){N:2}{S:3}1").unwrap();
        assert_eq!(mol.graph.node_count(), 4);
        assert_eq!(mol.graph.edge_count(), 4);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(2)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(3)), 0);
        assert!(mol.charges.is_empty());
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::O);
        assert_eq!(mol.graph[NodeIndex::new(2)], Atom::N);
        assert_eq!(mol.graph[NodeIndex::new(3)], Atom::S);
        assert_eq!(mol.graph.edge_weights().copied().counts().len(), 2);
        assert_eq!(
            mol.graph.edge_weights().copied().counts()[&BondType::Single],
            3
        );
        assert_eq!(
            mol.graph.edge_weights().copied().counts()[&BondType::Double],
            1
        );
        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("C1(=O)NS1")
                .unwrap()
                .build()
                .to_inchi_key()
        );
    }

    #[test]
    fn indexed_stereo() {
        let mol = MolBuilder::from_indexed_smiles("{N:0}[{C:1}@@{H:100}]({C:3}){C:4}(={O:5}){O:6}")
            .unwrap();
        assert_eq!(mol.graph.node_count(), 7);
        assert_eq!(mol.graph.edge_count(), 6);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 2);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 0);
        // assert_eq!(mol.h_count(&NodeIndex::new(3)), 3);
        assert_eq!(mol.h_count(&NodeIndex::new(4)), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(5)), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(6)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(100)), 0);
        assert!(mol.charges.is_empty());
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::N);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(3)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(4)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(5)], Atom::O);
        assert_eq!(mol.graph[NodeIndex::new(6)], Atom::O);
        assert_eq!(mol.graph[NodeIndex::new(100)], Atom::H);
        assert_eq!(mol.graph.edge_weights().copied().counts().len(), 2);
        assert_eq!(
            mol.graph.edge_weights().copied().counts()[&BondType::Single],
            5
        );
        assert_eq!(
            mol.graph.edge_weights().copied().counts()[&BondType::Double],
            1
        );
        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("N[CH](C)C(=O)O")
                .unwrap()
                .build()
                .to_inchi_key()
        );
    }

    #[test]
    fn indexed_multi_component() {
        let mol = MolBuilder::from_indexed_smiles("{C:0}.{N:1}").unwrap();
        assert_eq!(mol.graph.node_count(), 2); // hydrogen is implicit
        assert_eq!(mol.graph.edge_count(), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 4);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 3);
        assert!(mol.charges.is_empty());
        assert!(mol.aromatic_atoms.is_empty());
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::N);
        #[cfg(feature = "inchi")]
        assert_eq!(
            mol.to_inchi_key(),
            MolBuilder::from_smiles("C.N")
                .unwrap()
                .build()
                .to_inchi_key()
        );
    }

    #[test]
    fn indexed_aromatic() {
        let mol =
            MolBuilder::from_indexed_smiles("[{Na:0}+].[{O:1}-]{c:2}1{c:3}{c:4}{c:5}{c:6}{c:7}1")
                .unwrap();
        assert_eq!(mol.graph.node_count(), 8); // hydrogen is implicit
        assert_eq!(mol.graph.edge_count(), 7);
        assert_eq!(mol.h_count(&NodeIndex::new(0)), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(1)), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(2)), 0);
        assert_eq!(mol.h_count(&NodeIndex::new(3)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(4)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(5)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(6)), 1);
        assert_eq!(mol.h_count(&NodeIndex::new(7)), 1);
        assert_eq!(mol.charge(&NodeIndex::new(0)), 1);
        assert_eq!(mol.charge(&NodeIndex::new(1)), -1);
        assert_eq!(mol.charge(&NodeIndex::new(2)), 0);
        assert_eq!(mol.charge(&NodeIndex::new(3)), 0);
        assert_eq!(mol.charge(&NodeIndex::new(4)), 0);
        assert_eq!(mol.charge(&NodeIndex::new(5)), 0);
        assert_eq!(mol.charge(&NodeIndex::new(6)), 0);
        assert_eq!(mol.charge(&NodeIndex::new(7)), 0);
        assert_eq!(
            mol.aromatic_atoms
                .iter()
                .copied()
                .map(NodeIndex::index)
                .sorted()
                .collect_vec(),
            vec![2, 3, 4, 5, 6, 7]
        );
        assert_eq!(mol.graph[NodeIndex::new(0)], Atom::Na);
        assert_eq!(mol.graph[NodeIndex::new(1)], Atom::O);
        assert_eq!(mol.graph[NodeIndex::new(2)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(3)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(4)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(5)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(6)], Atom::C);
        assert_eq!(mol.graph[NodeIndex::new(7)], Atom::C);
        assert_eq!(mol.graph.edge_weights().copied().counts().len(), 2);
        assert_eq!(
            mol.graph.edge_weights().copied().counts()[&BondType::Single],
            1
        );
        assert_eq!(
            mol.graph.edge_weights().copied().counts()[&BondType::Aromatic],
            6
        );

        // can't use inchikeys for aromatic bonds
    }
}
