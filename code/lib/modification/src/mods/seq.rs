use super::{Mod, ModificationError};
use molecule::Mol;
use petgraph::visit::DfsPostOrder;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::{fmt, iter::FromIterator};

/// A sequence of modifications to be applied.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct ModSeq {
    /// The vector of modifications
    mods: Vec<Mod>,
}

impl fmt::Display for ModSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for curr_mod in &self.mods {
            writeln!(f, "{}", curr_mod)?;
        }
        Ok(())
    }
}

impl ModSeq {
    /// Creates a new ModSeq from a vector of modifications.
    pub fn new(mods: Vec<Mod>) -> Self {
        Self { mods }
    }

    /// Applies the modifications to the given target molecule.
    ///
    /// First applies all AddNode modifications and then applies all remaining
    /// modifications in order.
    pub fn modify(&self, mol: &mut Mol) -> Result<(), ModificationError> {
        let mut mapping = FxHashMap::<i32, i32>::default();

        // build the index mapping from the add node modifications first
        for modification in self.mods.iter() {
            if let Mod::AddNode(_, _, _, idx) = modification {
                // unwrap is safe here since the AddNode modify always results the Some variant
                let new_node_idx = modification.modify(mol)?.unwrap();
                // this will cause problems if the usize somehow overflows an i32
                // at that point we have bigger problems (a 2 billion node molecule?)
                mapping.insert(*idx, new_node_idx.index() as i32);
            }
        }

        // actually apply all the modifications in sequence, skipping the add node ones since we
        // already did those in the previous loop
        for modification in self.iter() {
            if matches!(modification, Mod::AddNode(_, _, _, _)) {
                continue;
            }
            modification.convert(&mapping)?.modify(mol)?;
        }

        mol.compact();

        Ok(())
    }

    /// Computes the set of atom indices in the modified molecule that were changed in some way.
    ///
    /// Any added atoms are automatically considered part of the modified site. Any atoms at the
    /// endpoints of a modified bond, whether it was removed or added, are also considered part of
    /// the modified site. Any heavy atoms whose hydrogen count or charge were changed are also
    /// considered part of the modified site. Finally, the neighbors of any removed atoms are also
    /// considered part of the modified site.
    ///
    /// # Errors
    /// Propagates errors from [`ModSeq::modify`].
    pub fn modified_site(&self, mol: &Mol) -> Result<FxHashSet<i32>, ModificationError> {
        let mut mol = mol.clone();

        let mut mapping = FxHashMap::<i32, i32>::default();
        let mut modified_atoms = FxHashSet::default();

        // build the index mapping from the add node modifications first
        for modification in self.mods.iter() {
            if let Mod::AddNode(_, _, _, idx) = modification {
                // unwrap is safe here since the AddNode modify always results the Some variant
                let new_node_idx = modification.modify(&mut mol)?.unwrap();
                // this will cause problems if the usize somehow overflows an i32
                // at that point we have bigger problems (a 2 billion node molecule?)
                mapping.insert(*idx, new_node_idx.index() as i32);
                modified_atoms.insert(new_node_idx.index() as i32);
            }
        }

        // actually apply all the modifications in sequence, skipping the add node ones since we
        // already did those in the previous loop
        for modification in self.iter() {
            if matches!(modification, Mod::AddNode(_, _, _, _)) {
                continue;
            }
            match modification.convert(&mapping)? {
                Mod::AddEdge(u, v, _) | Mod::RemoveEdge(u, v) => {
                    modified_atoms.insert(u);
                    modified_atoms.insert(v);
                }
                Mod::ChargeShift(u, _) | Mod::HydrogenShift(u, _) => {
                    modified_atoms.insert(u);
                }
                Mod::RemoveNode(u) => {
                    modified_atoms.remove(&u);
                    let node_index = molecule::NodeIndex::new(u as usize);
                    let neighbors = mol.graph.neighbors(node_index).map(|u| u.index() as i32);
                    modified_atoms.extend(neighbors);
                }
                Mod::RemoveCC(u) => {
                    let node = molecule::NodeIndex::new(u as usize);

                    let mut dfs = DfsPostOrder::new(&mol.graph, node);
                    let mut removed_atoms = FxHashSet::default();
                    while let Some(node) = dfs.next(&mol.graph) {
                        removed_atoms.insert(node);
                        modified_atoms.remove(&(node.index() as i32));
                    }

                    modified_atoms.extend(
                        removed_atoms
                            .iter()
                            .flat_map(|&node| mol.graph.neighbors(node))
                            .collect::<FxHashSet<_>>()
                            .difference(&removed_atoms)
                            .map(|u| u.index() as i32),
                    );
                }
                Mod::AddNode(_, _, _, _) => (),
            }
        }

        Ok(modified_atoms)
    }

    /// Computes the net mass change shift would occur from applying the mod sequence
    /// to the target molecule.assert_eq!
    pub fn mass_shift(&self, mol: &Mol) -> f64 {
        self.mods.iter().map(|modif| modif.mass_shift(mol)).sum()
    }

    /// Returns an iterator for the mod sequence.
    pub fn iter(&self) -> impl Iterator<Item = &Mod> {
        self.mods.iter()
    }

    /// Apply the given node index offset to all the modifications in the sequence.
    pub fn apply_offset(self, offset: i32) -> Self {
        self.iter().map(|x| x.apply_offset(offset)).collect()
    }

    /// Perform a validity check on each of the modifications in the sequence.
    pub fn convert(self, mapping: &FxHashMap<i32, i32>) -> Result<Self, ModificationError> {
        let mut mods = Vec::new();
        for modif in self.mods {
            mods.push(modif.convert(mapping)?);
        }
        Ok(Self { mods })
    }
}

impl IntoIterator for ModSeq {
    type Item = Mod;

    type IntoIter = <Vec<Mod> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.mods.into_iter()
    }
}

impl FromIterator<Mod> for ModSeq {
    fn from_iter<T: IntoIterator<Item = Mod>>(iter: T) -> Self {
        Self {
            mods: iter.into_iter().collect(),
        }
    }
}
