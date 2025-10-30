use std::{
    collections::{HashMap, HashSet},
    convert::TryInto,
};

use crate::Mod;
use molecule::{Atom, BondType, Mol, NodeIndex};
use petgraph::graph::node_index;

use crate::errors::ValidationError;

pub type Atoms = Vec<(Atom, i32)>;
pub type Bonds = Vec<(i32, i32, BondType)>;

#[derive(Debug)]
pub enum Field {
    Disconnect(i32, i32),
    Remove(i32),
    Connect(i32, i32, BondType),
    Charge(i32, i8),
    Add {
        start: i32,
        end: i32,
        bond: BondType,
        atoms: Atoms,
        bonds: Option<Bonds>,
    },
}

impl From<(i32, i32)> for Field {
    fn from((start, end): (i32, i32)) -> Self {
        Self::Disconnect(start, end)
    }
}

impl From<i32> for Field {
    fn from(idx: i32) -> Self {
        Self::Remove(idx)
    }
}

impl From<(i32, i32, BondType)> for Field {
    fn from((start, end, bond): (i32, i32, BondType)) -> Self {
        Self::Connect(start, end, bond)
    }
}

impl From<(i32, i32, BondType, Atoms, Option<Bonds>)> for Field {
    fn from((start, end, bond, atoms, bonds): (i32, i32, BondType, Atoms, Option<Bonds>)) -> Self {
        Self::Add {
            start,
            end,
            bond,
            atoms,
            bonds,
        }
    }
}

impl From<(i32, i8)> for Field {
    fn from((idx, charge): (i32, i8)) -> Self {
        Self::Charge(idx, charge)
    }
}

fn node(idx: i32) -> NodeIndex {
    node_index(idx.try_into().unwrap())
}

impl Field {
    pub fn apply(&self, mol: &Mol, new_idx: &mut i32) -> Result<Option<Vec<Mod>>, ValidationError> {
        match self {
            Field::Disconnect(idx1, idx2) => {
                let (idx1, idx2) = (*idx1 - 1, *idx2 - 1);
                let atom1 = mol
                    .graph
                    .node_weight(node(idx1))
                    .ok_or(ValidationError::InvalidIndex(idx1))?;
                let atom2 = mol
                    .graph
                    .node_weight(node(idx2))
                    .ok_or(ValidationError::InvalidIndex(idx2))?;

                Ok(Some(vec![if atom1 == &Atom::H {
                    Mod::HydrogenShift(idx2, -1)
                } else if atom2 == &Atom::H {
                    Mod::HydrogenShift(idx1, -1)
                } else {
                    Mod::RemoveEdge(idx1, idx2)
                }]))
            }
            Field::Charge(idx, charge) => {
                let idx = *idx - 1;
                Ok(Some(vec![Mod::ChargeShift(idx, *charge)]))
            }
            Field::Remove(idx) => {
                let idx = *idx - 1;
                // if it's a hydrogen the disconnect removes it
                Ok((mol.graph[node(idx)] != Atom::H).then(|| vec![Mod::RemoveCC(idx)]))
            }
            Field::Connect(idx1, idx2, bond) => {
                Ok(Some(vec![Mod::AddEdge(*idx1 - 1, *idx2 - 1, *bond)]))
            }
            Field::Add {
                start,
                end,
                bond,
                atoms,
                bonds,
            } => {
                if (atoms.len() > 1) && bonds.is_none() {
                    return Err(ValidationError::MissingBonds);
                }
                if atoms.len() == 1 {
                    return Ok(Some(self.add_single_atom(new_idx)));
                }

                // check that indices only used once
                let mut indices = atoms.iter().map(|(_, idx)| *idx).collect::<Vec<_>>();
                indices.sort_unstable();
                for (prev, curr) in indices.iter().zip(indices.iter().skip(1)) {
                    if *prev == *curr {
                        return Err(ValidationError::DuplicateIndex(*prev));
                    }
                }

                // store the indices of newly added hydrogens
                let hydrogens = atoms
                    .iter()
                    .filter_map(|(atom, idx)| (atom == &Atom::H).then(|| *idx))
                    .collect::<HashSet<i32>>();

                // translate from old modification Add directive indices to new, negative indices
                let translation = atoms
                    .iter()
                    // we skip hydrogens here since they won't be explicit nodes
                    .filter_map(|(atom, idx)| (atom != &Atom::H).then(|| *idx))
                    // now map these old indices to the new negative ones
                    .enumerate()
                    .map(|(offset, old_idx)| (old_idx, *new_idx - offset as i32))
                    // and turn into a hash map
                    .collect::<HashMap<i32, i32>>();

                // update running negative index allowance
                *new_idx -= translation.len() as i32;

                // map from new negative indices to atoms to add
                let mut atoms = atoms
                    .iter()
                    // we skip hydrogens here since they won't be explicit nodes
                    // map from negative index to (atom, h_count) tuples
                    .filter_map(|(atom, idx)| {
                        (atom != &Atom::H).then(|| (translation[idx], (*atom, 0u8)))
                    })
                    .collect::<HashMap<i32, (Atom, u8)>>();

                // add edges
                let bonds = bonds.as_ref().unwrap();
                let mut bond_mods = Vec::<Mod>::with_capacity(bonds.len());
                for (start_idx, end_idx, bt) in bonds.iter().copied() {
                    // increment hydrogen count if one of the atoms in the bond is a hydrogen
                    if hydrogens.contains(&start_idx) {
                        let to_add = atoms
                            .get_mut(
                                translation
                                    .get(&end_idx)
                                    .ok_or(ValidationError::InvalidIndex(end_idx))?,
                            )
                            .unwrap();
                        if bt != BondType::Single {
                            return Err(ValidationError::HydrogenBond {
                                h: start_idx,
                                other: end_idx,
                                bond: bt,
                            });
                        }
                        to_add.1 += 1;
                        continue;
                    }
                    if hydrogens.contains(&end_idx) {
                        let to_add = atoms
                            .get_mut(
                                translation
                                    .get(&start_idx)
                                    .ok_or(ValidationError::InvalidIndex(start_idx))?,
                            )
                            .unwrap();
                        if bt != BondType::Single {
                            return Err(ValidationError::HydrogenBond {
                                h: end_idx,
                                other: start_idx,
                                bond: bt,
                            });
                        }
                        to_add.1 += 1;
                        continue;
                    }

                    // make sure that the indices we have are valid
                    for idx in [start_idx, end_idx] {
                        if !mol
                            .graph
                            .node_indices()
                            .any(|x| x.index() == (idx - 1) as usize)
                            && !translation.contains_key(&idx)
                        {
                            eprintln!("helo {:?}", self);
                            return Err(ValidationError::InvalidIndex(idx));
                        }
                    }

                    let start_idx = *translation.get(&start_idx).unwrap_or(&(start_idx - 1));
                    let end_idx = *translation.get(&end_idx).unwrap_or(&(end_idx - 1));
                    bond_mods.push(Mod::AddEdge(start_idx, end_idx, bt));
                }

                let mut mods = atoms
                    .into_iter()
                    .map(|(idx, (atom, h_count))| Mod::AddNode(atom, h_count, 0, idx))
                    .collect::<Vec<Mod>>();
                // add first edge
                mods.push(Mod::AddEdge(*start - 1, translation[end], *bond));
                // then all the edges
                mods.extend(bond_mods);

                Ok(Some(mods))
            }
        }
    }

    fn add_single_atom(&self, new_idx: &mut i32) -> Vec<Mod> {
        match self {
            Field::Add {
                start,
                bond,
                atoms,
                bonds,
                ..
            } => {
                let start = *start - 1;
                assert_eq!(atoms.len(), 1);
                let (atom, idx) = atoms[0];
                if atom == Atom::H {
                    // in this case we just add a hydrogen to the start node
                    return vec![Mod::HydrogenShift(start, 1)];
                }

                // here we create the node first
                // then we connect it to the start
                let mut mods = vec![
                    Mod::AddNode(atom, 0, 0, *new_idx),
                    Mod::AddEdge(start, *new_idx, *bond),
                ];

                if let Some(bonds) = bonds {
                    for (mut start, mut end, bt) in bonds {
                        if start == idx {
                            start = *new_idx;
                        } else {
                            start -= 1;
                        }
                        if end == idx {
                            end = *new_idx;
                        } else {
                            end -= 1;
                        }
                        mods.push(Mod::AddEdge(start, end, *bt))
                    }
                }

                // make sure we don't use this dummy index again
                *new_idx -= 1;

                mods
            }
            _ => unreachable!(),
        }
    }
}
