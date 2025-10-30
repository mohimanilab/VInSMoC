use constants::MASS_ELECTRON;
use molecule::{Atom, BondType, Mol, NodeIndex};
use petgraph::{
    graph::node_index,
    visit::{DfsPostOrder, EdgeRef, IntoEdgeReferences},
};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::iter::FromIterator;
use std::num::TryFromIntError;
use std::{
    convert::{TryFrom, TryInto},
    fmt,
};

use crate::{ModificationError, ModificationErrorKind};

mod field;
mod parse;
mod seq;
mod subgraph;

pub use parse::parse_mod;
pub use seq::ModSeq;
pub use subgraph::SubgraphMod;

/// Represents all possible modifications that can be applied to a molecule.
///
/// Each variant includes [`i32`] node indices representing one or more nodes in the molecule that
/// are to be modified. An [`i32`] instead of a [`NodeIndex`](petgraph::graph::NodeIndex) is used
/// to allow for a different convention when adding nodes. The [`Mod::AddNode`] variant
/// accepts a [`molecule::Atom`], a hydrogen count, a charge, and a dummy index. The safest way to
/// use this dummy index is to make it negative, that way it will never clash with existing node
/// indices in the molecular graph.
///
/// See also [`crate::ModSeq`] for a sequence of logically connected modifications applied
/// to the same molecule.
///
/// # Examples
/// ```
/// use modification::{Mod, ModificationError, ModificationErrorKind};
/// use molecule::parsers::parse_smiles;
/// use petgraph::graph::node_index;
///
/// // propane
/// let mut mol = parse_smiles("CCC").unwrap();
///
/// // remove a hydrogen from the last carbon (index 2)
/// let remove_h = Mod::HydrogenShift(2, -1);
/// assert_eq!(mol.h_count(&node_index(2)), 3);
/// assert_eq!(remove_h.modify(&mut mol).unwrap(), None);
/// assert_eq!(mol.h_count(&node_index(2)), 2);
///
/// // trying to remove a node that doesn't exist gives an error
/// let failed_remove = Mod::RemoveNode(10);
/// assert_eq!(
///     failed_remove.modify(&mut mol).unwrap_err(),
///     ModificationError::new(failed_remove, ModificationErrorKind::NoSuchNode(node_index(10)))
/// );
/// ```
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Mod {
    /// Adds a new edge between the two nodes with given bond type
    AddEdge(i32, i32, BondType),
    /// Removes the edge connecting the two nodes
    RemoveEdge(i32, i32),
    /// Adds a new node of with a given atom, hydrogen count, charge, and temporary index
    AddNode(Atom, u8, i8, i32),
    /// Removes the given node
    RemoveNode(i32),
    /// Removes all nodes in the provided node's connected component
    RemoveCC(i32),
    /// Changes the given node's charge
    ChargeShift(i32, i8),
    /// Changes the give node's hydrogen count
    HydrogenShift(i32, i8),
}

impl fmt::Display for Mod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Mod::AddEdge(u, v, bond) => write!(f, "AddEdge {} {} {}", u, v, bond),
            Mod::RemoveEdge(u, v) => write!(f, "RemoveEdge {} {}", u, v),
            Mod::AddNode(atom, h, c, i) => write!(f, "AddNode {} {} {} {}", atom, h, c, i),
            Mod::RemoveNode(u) => write!(f, "RemoveNode {}", u),
            Mod::RemoveCC(u) => write!(f, "RemoveCC {}", u),
            Mod::ChargeShift(u, shift) => write!(f, "ChargeShift {} {}", u, shift),
            Mod::HydrogenShift(u, shift) => write!(f, "HydrogenShift {} {}", u, shift),
        }
    }
}

#[inline(always)]
fn num2idx(num: i32) -> Result<NodeIndex, TryFromIntError> {
    let idx_num = usize::try_from(num)?;
    Ok(NodeIndex::new(idx_num))
}

#[inline(always)]
fn convert_idx(num: i32, mapping: &FxHashMap<i32, i32>) -> Result<i32, TryFromIntError> {
    mapping.get(&num).map(|x| Ok(*x)).unwrap_or(Ok(num))
}

impl Mod {
    /// Checks whether the modification is valid.
    pub fn convert(self, mapping: &FxHashMap<i32, i32>) -> Result<Self, ModificationError> {
        match self {
            Self::AddEdge(node1, node2, bond) => {
                let node1 = convert_idx(node1, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node1))
                })?;
                let node2 = convert_idx(node2, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node2))
                })?;
                if node1 == node2 {
                    return Err(ModificationError::new(
                        self,
                        ModificationErrorKind::SelfLoop(node1),
                    ));
                }
                Ok(Self::AddEdge(node1, node2, bond))
            }
            Self::RemoveEdge(node1, node2) => {
                let node1 = convert_idx(node1, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node1))
                })?;
                let node2 = convert_idx(node2, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node2))
                })?;
                Ok(Self::RemoveEdge(node1, node2))
            }
            Self::AddNode(_, _, _, _) => Ok(self),
            Self::RemoveNode(node) => {
                let node = convert_idx(node, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node))
                })?;
                Ok(Self::RemoveNode(node))
            }
            Self::RemoveCC(node) => {
                let node = convert_idx(node, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node))
                })?;
                Ok(Self::RemoveCC(node))
            }
            Self::ChargeShift(node, shift) => {
                let node = convert_idx(node, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node))
                })?;
                Ok(Self::ChargeShift(node, shift))
            }
            Self::HydrogenShift(node, shift) => {
                let node = convert_idx(node, mapping).map_err(|_| {
                    ModificationError::new(self, ModificationErrorKind::NoIndexMapping(node))
                })?;
                Ok(Self::HydrogenShift(node, shift))
            }
        }
    }

    /// Modify the given molecule.
    ///
    /// If the modification was [`Mod::AddNode`] then this returns the true node index of the added
    /// node. Otherwise the ok variant of the returned [`Result`] can be ignored (it will be [`None`]).
    pub fn modify(&self, mol: &mut Mol) -> Result<Option<NodeIndex>, ModificationError> {
        match *self {
            Self::AddEdge(node1, node2, bond) => {
                // if either node doesn't exist return an error
                let node1 = num2idx(node1).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node1))
                })?;
                mol.graph.node_weight(node1).ok_or_else(|| {
                    ModificationError::new(*self, ModificationErrorKind::NoSuchNode(node1))
                })?;
                let node2 = num2idx(node2).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node2))
                })?;
                mol.graph.node_weight(node2).ok_or_else(|| {
                    ModificationError::new(*self, ModificationErrorKind::NoSuchNode(node2))
                })?;
                mol.graph.update_edge(node1, node2, bond);
            }
            Self::RemoveEdge(node1, node2) => {
                let node1 = num2idx(node1).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node1))
                })?;
                let node2 = num2idx(node2).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node2))
                })?;
                let edge_idx = mol.graph.find_edge(node1, node2).ok_or_else(|| {
                    ModificationError::new(*self, ModificationErrorKind::NoSuchEdge(node1, node2))
                })?;
                mol.graph.remove_edge(edge_idx);
            }
            Self::AddNode(atom, hydrogens, charge, _) => {
                let v = mol.add_node(atom, hydrogens, charge, false);
                return Ok(Some(v));
            }
            Self::RemoveNode(node) => {
                let node = num2idx(node).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node))
                })?;
                mol.remove_node(node).ok_or_else(|| {
                    ModificationError::new(*self, ModificationErrorKind::NoSuchNode(node))
                })?;
            }
            Self::RemoveCC(node) => {
                let node = num2idx(node).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node))
                })?;
                if !mol.graph.contains_node(node) {
                    return Err(ModificationError::new(
                        *self,
                        ModificationErrorKind::NoSuchNode(node),
                    ));
                }

                let mut dfs = DfsPostOrder::new(&mol.graph, node);

                while let Some(node) = dfs.next(&mol.graph) {
                    mol.remove_node(node);
                }
                mol.compact();
            }
            Self::ChargeShift(node, shift) => {
                let node = num2idx(node).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node))
                })?;
                let charge = mol.charge_entry(node);
                *charge += shift;
            }
            Self::HydrogenShift(node, shift) => {
                let node = num2idx(node).map_err(|_| {
                    ModificationError::new(*self, ModificationErrorKind::NegativeIndex(node))
                })?;
                let hcount = mol.h_entry(node);
                if *hcount == 0 && shift < 0 {
                    return Err(ModificationError::new(
                        *self,
                        ModificationErrorKind::NegativeHydrogen(node),
                    ));
                }
                let new_hcount = TryInto::<i8>::try_into(*hcount).unwrap() + shift;
                *hcount = new_hcount.try_into().unwrap();
            }
        };

        Ok(None)
    }

    /// Computes the change in mass that would result from a mod application.
    pub fn mass_shift(&self, mol: &Mol) -> f64 {
        match *self {
            Self::AddNode(atom, hydrogens, charge, _) => {
                atom.exact_mass() + (hydrogens as f64) * Atom::H.exact_mass()
                    - (charge as f64) * MASS_ELECTRON
            }
            Self::RemoveNode(node) => {
                let node = node_index(node as usize);
                -mol.node_exact_mass(node)
            }
            Self::RemoveCC(node) => {
                let node = node_index(node as usize);
                let mut mass = 0.0;
                let mut dfs = DfsPostOrder::new(&mol.graph, node);
                while let Some(node) = dfs.next(&mol.graph) {
                    mass -= mol.node_exact_mass(node);
                }
                mass
            }
            Self::HydrogenShift(_, shift) => (shift as f64) * Atom::H.exact_mass(),
            _ => 0.0,
        }
    }

    /// Adds the given offset to the node indices in the modifications
    pub fn apply_offset(&self, offset: i32) -> Self {
        match *self {
            Self::AddEdge(x, y, bt) => Self::AddEdge(x + offset, y + offset, bt),
            Self::RemoveEdge(x, y) => Self::RemoveEdge(x + offset, y + offset),
            Self::AddNode(_, _, _, _) => *self,
            Self::RemoveNode(x) => Self::RemoveNode(x + offset),
            Self::RemoveCC(x) => Self::RemoveCC(x + offset),
            Self::ChargeShift(x, shift) => Self::ChargeShift(x + offset, shift),
            Self::HydrogenShift(x, shift) => Self::HydrogenShift(x + offset, shift),
        }
    }
}

/// This function finds the reaction necessary to go from reactant to product.
///
/// The reactant and product are both represented as a [`molecule::Mol`], and the
/// output is a [`subgraph::SubgraphMod`]. This function assumes that the atoms in
/// the reactant and product have previously been mapped. This means that if an atom is
/// in both the reactant and product, its node index is the same in both the reactant
/// and product graphs.
///
/// # Examples
/// ```
/// use modification::{reaction, Mod};
/// use molecule::{BondType, MolBuilder, NodeIndex};
///
/// let reactant = MolBuilder::from_smiles("CCCC").unwrap().build();
/// let product = MolBuilder::from_smiles("CC=CC").unwrap().build();
///
/// let changes = reaction(&reactant, &product, "mod_name".into());
/// assert_eq!(
///     changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
///     vec![
///         Mod::HydrogenShift(1, -1),
///         Mod::HydrogenShift(2, -1),
///         Mod::AddEdge(1, 2, BondType::Double),
///     ]
/// );
/// assert_eq!(changes.pattern.graph.node_count(), 2);
/// assert!(changes.pattern.graph.contains_node(NodeIndex::new(1)));
/// assert!(changes.pattern.graph.contains_node(NodeIndex::new(2)));
/// assert_eq!(changes.pattern.h_count(&NodeIndex::new(1)), 2);
/// assert_eq!(changes.pattern.h_count(&NodeIndex::new(2)), 2);
/// assert_eq!(changes.pattern.charge(&NodeIndex::new(1)), 0);
/// assert_eq!(changes.pattern.charge(&NodeIndex::new(2)), 0);
/// assert_eq!(&changes.name, "mod_name");
/// ```
pub fn reaction(reactant: &Mol, product: &Mol, name: String) -> SubgraphMod {
    let mut result = Vec::new();

    for idx in reactant.graph.node_indices() {
        if !product.graph.contains_node(idx) {
            // Record all atoms that are present in reactants but not in product.
            result.push(Mod::RemoveNode(idx.index() as i32));
        } else {
            // If other attributes for the atom differ, change them.
            let h_count = reactant.h_count(&idx);
            let prod_h_count = product.h_count(&idx);
            if h_count != prod_h_count {
                result.push(Mod::HydrogenShift(
                    idx.index() as i32,
                    prod_h_count as i8 - h_count as i8,
                ));
            }
            let charge = reactant.charge(&idx);
            let prod_charge = product.charge(&idx);
            if charge != prod_charge {
                result.push(Mod::ChargeShift(idx.index() as i32, prod_charge - charge));
            }
        }
    }

    // Record all atoms that are present in the product but not reactant and add them.
    for idx in product.graph.node_indices() {
        if !reactant.graph.contains_node(idx) {
            let atom = product.graph[idx];
            let h_count = product.h_count(&idx);
            let charge = product.charge(&idx);
            result.push(Mod::AddNode(atom, h_count, charge, idx.index() as i32));
        }
    }

    // Record bonds that are added or changed in the reaction.
    for edge in product.graph.edge_references() {
        let idx1 = edge.source();
        let idx2 = edge.target();
        let bond_type = *edge.weight();
        match reactant.graph.find_edge(idx1, idx2) {
            Some(x) => {
                let original_bond_type = reactant.graph.edge_weight(x).unwrap();
                if *original_bond_type != bond_type {
                    result.push(Mod::AddEdge(
                        idx1.index() as i32,
                        idx2.index() as i32,
                        bond_type,
                    ));
                }
            }
            None => result.push(Mod::AddEdge(
                idx1.index() as i32,
                idx2.index() as i32,
                bond_type,
            )),
        }
    }

    // Record bonds that are broken in the reaction, but the atoms at both ends are still
    // part of the product.
    for edge in reactant.graph.edge_references() {
        let idx1 = edge.source();
        let idx2 = edge.target();
        if !product.graph.contains_node(idx1) || !product.graph.contains_node(idx2) {
            continue;
        }

        if !product.graph.contains_edge(idx1, idx2) {
            result.push(Mod::RemoveEdge(idx1.index() as i32, idx2.index() as i32))
        }
    }

    let mod_seq = ModSeq::new(result);
    let center_indices = reaction_center(reactant, &mod_seq, 0);
    let pattern = reactant.to_submol(FxHashSet::from_iter(center_indices));

    SubgraphMod::new(pattern, mod_seq, name)
}

/// This function finds the reaction center for a given reaction.
///
/// The reactant is represented as a reference to a [`molecule::Mol`], and the
/// modifications are represented as a reference to a [`seq::ModSeq`]. The
/// output is a vector of [`molecule::NodeIndex`]. The radius is a [`u8`] which
/// specifies how far from changed atoms an atom has to be to be included in the
/// output. This means if the radius is 2, the output includes all atoms which
/// are changed between the reaction and product and all atoms that can be reached
/// from any of the changed atoms by 2 molecular bonds or less. This function
/// assumes that the node indices for each atom are the same in the reactant and
/// product, and if an atom is present in one and not the other then its index is
/// not used when the atom is not present.
///
/// # Examples
/// ```
/// use modification::{reaction_center, Mod, ModSeq};
/// use molecule::{Atom, BondType, Mol, NodeIndex, MolBuilder, MolGraph};
/// use petgraph::graph::node_index;
///
/// let reactant = MolBuilder::from_smiles("CC=CC").unwrap().build();
/// let mods = ModSeq::new(vec![Mod::RemoveEdge(1, 2)]);
///
/// let center = reaction_center(&reactant, &mods, 0);
/// assert_eq!(center, vec![NodeIndex::new(1), NodeIndex::new(2)]);
/// ```
pub fn reaction_center(reactant: &Mol, mods: &ModSeq, radius: u8) -> Vec<NodeIndex> {
    let mut prev_recorded = FxHashSet::default(); // for storing information for bfs
    let mut reaction_atoms = FxHashSet::default(); // for storing output

    for modification in mods.iter() {
        match modification {
            Mod::AddEdge(idx1, idx2, _) => {
                prev_recorded.insert(*idx1 as usize);
                prev_recorded.insert(*idx2 as usize);
                reaction_atoms.insert(*idx1 as usize);
                reaction_atoms.insert(*idx2 as usize);
            }
            Mod::RemoveEdge(idx1, idx2) => {
                prev_recorded.insert(*idx1 as usize);
                prev_recorded.insert(*idx2 as usize);
                reaction_atoms.insert(*idx1 as usize);
                reaction_atoms.insert(*idx2 as usize);
            }
            Mod::RemoveNode(idx) => {
                prev_recorded.insert(*idx as usize);
                reaction_atoms.insert(*idx as usize);
            }
            Mod::ChargeShift(idx, _) => {
                prev_recorded.insert(*idx as usize);
                reaction_atoms.insert(*idx as usize);
            }
            Mod::HydrogenShift(idx, _) => {
                prev_recorded.insert(*idx as usize);
                reaction_atoms.insert(*idx as usize);
            }
            _ => (),
        }
    }

    // use BFS to find the atoms within the specified radius of the reacting atoms
    for _ in 0..radius {
        let mut curr_recorded = FxHashSet::default();
        for idx in prev_recorded.iter() {
            for neighbor in reactant.graph.neighbors(NodeIndex::new(*idx as usize)) {
                if !reaction_atoms.contains(&neighbor.index()) {
                    curr_recorded.insert(neighbor.index());
                    reaction_atoms.insert(neighbor.index());
                }
            }
        }
        prev_recorded = curr_recorded;
    }

    reaction_atoms.iter().map(|x| NodeIndex::new(*x)).collect()
}

#[cfg(test)]
mod tests;
