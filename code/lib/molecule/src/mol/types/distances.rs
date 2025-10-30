use std::{collections::VecDeque, num::TryFromIntError};

use super::{Mol, NodeIndex};
use itertools::Itertools;
use learning::Distance;
use ordered_float::OrderedFloat;
use rustc_hash::{FxHashMap, FxHashSet};
use std::convert::AsRef;

pub(crate) type Of64 = OrderedFloat<f64>;

/// A wrapper struct to denote a pair of node indices, such that the first node index is always
/// smaller than the second one.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct AtomPair(NodeIndex, NodeIndex);

impl AtomPair {
    fn new(node1: NodeIndex, node2: NodeIndex) -> Self {
        if node1.index() < node2.index() {
            Self(node1, node2)
        } else {
            Self(node2, node1)
        }
    }
}

/// A structure that is used to derive information about pairwise distances between atoms within a
/// given molecule.
///
/// This functions by performing a one-time expensive operation upon construction to calculate
/// pairwise distances between all atoms in `molecule`.
#[derive(Clone, Debug)]
pub struct DistanceInfo<'a> {
    /// A reference to the molecule being described.
    pub molecule: &'a Mol,
    /// A map which is built during initialization which stores the pairwise distances between each
    /// pair of atoms within a molecule.
    distance_map: FxHashMap<AtomPair, u16>,
}

impl<'a> DistanceInfo<'a> {
    /// Runs BFS on each node in the molecule to compute all pairwise distances between atom
    /// indices, and stores the result.
    pub fn new(molecule: &'a Mol) -> Self {
        let distance_map = Self::all_pairwise_distances(molecule);
        Self {
            molecule,
            distance_map,
        }
    }

    /// Calculate the symmetric set distance between two sets of atoms within the same molecule.
    ///
    /// This is calculated as the average of the asymmetric set distance of `set1` to `set2` and
    /// the asymmetric set distance of `set2` to `set1`. For more information about the asymmetric
    /// set distance, see [`Self::asymmetric_set_distance`].
    ///
    /// # Panics
    /// If `set1` or `set2` are empty then the distance is not meaningful, and therefore this
    /// function panics. This function also panics if a `NodeIndex` is provided that does not
    /// belong in the molecule that `self` represents.
    pub fn set_distance(&self, set1: &FxHashSet<NodeIndex>, set2: &FxHashSet<NodeIndex>) -> Of64 {
        let dist_fwd = self.asymmetric_set_distance(set1, set2);
        let dist_rev = self.asymmetric_set_distance(set2, set1);
        Of64::from((dist_fwd + dist_rev) / 2.)
    }

    /// Calculate the distance of `source` from a set of atoms `set`.
    ///
    /// This distance is defined to be the shortest distance from `source` to any point within
    /// `set`. Shortest distance in a molecular graph is defined as the number of edges (or bonds)
    /// that must be traversed to get from one atom to another.
    ///
    /// # Panics
    /// If `set` is empty then the distance is not meaningful, and therefore this function panics.
    /// This function also panics if a `NodeIndex` is provided that does not belong in the molecule
    /// that `self` represents.
    pub fn distance_from_set(&self, source: NodeIndex, set: &FxHashSet<NodeIndex>) -> u16 {
        set.iter()
            .map(|pt2| {
                *self
                    .distance_map
                    .get(&AtomPair::new(source, *pt2))
                    .expect("all of the node indices should be valid")
            })
            .min()
            .expect("the provided set of node indices should not be empty")
    }

    /// Calculate the asymmetric distance of `set1` from `set2`, where `set1` and `set2` refer to
    /// possibly overlapping sets of `NodeIndex`'s within the molecule that `self` refers to.
    ///
    /// This set distance is defined as the average distance of atoms from `set1` from `set2`. For
    /// a definition on what the distance between an atom from `set1` and `set2` would mean, see
    /// [`Self::distance_from_set`].
    ///
    /// # Panics
    /// If `set2` is empty then the distance is not meaningful, and therefore this function panics.
    /// This function also panics if a `NodeIndex` is provided that does not belong in the molecule
    /// that `self` represents.
    pub fn asymmetric_set_distance(
        &self,
        set1: &FxHashSet<NodeIndex>,
        set2: &FxHashSet<NodeIndex>,
    ) -> f64 {
        let cum_dist = set1
            .iter()
            .map(|pt1| self.distance_from_set(*pt1, set2))
            .sum::<u16>() as f64;

        cum_dist / (set1.len() as f64)
    }

    /// A function to calculate all pairwise distances between the atoms in a molecule.
    ///
    /// Using Seidel's algorithm requires some crazy matrix algebra and also 5 (|V| x |V|) matrix
    /// allocations during each recursive step to achieve a O(V^(2.3) log V) runtime. Simply
    /// running BFS on a single source takes O(|V| + |E|), so running it on all sources takes
    /// O(|V|^2 |E|*|V|). However in molecular graphs, |E| is equal to |V| + (# cycles) - 1,
    /// which is usually not significantly higher than |V|, so the amortized runtime is comparable
    /// if not better in most cases.
    fn all_pairwise_distances(molecule: &Mol) -> FxHashMap<AtomPair, u16> {
        let num_atoms = molecule.graph.node_count();
        let mut distance_map = FxHashMap::with_capacity_and_hasher(
            num_atoms * num_atoms,
            std::default::Default::default(),
        );

        // get this memory once since each hashmap is the exact same size
        let mut distances_from_source = FxHashMap::default();
        distances_from_source.reserve(num_atoms);

        // keep track of node and distance
        let mut bfs_queue = VecDeque::<(NodeIndex, u16)>::with_capacity(num_atoms);

        for source_node in molecule.graph.node_indices() {
            Self::all_distances_from_source(
                molecule,
                source_node,
                &mut distances_from_source,
                &mut bfs_queue,
            );

            for (neigh_node, dist) in &distances_from_source {
                // use the `AtomPair` struct to save half the memory
                distance_map.insert(AtomPair::new(source_node, *neigh_node), *dist);
            }

            distances_from_source.clear();
            bfs_queue.clear();
        }

        distance_map
    }

    /// A function which calculates pairwise distances between all node `source` and all other
    /// nodes in `molecule`.
    ///
    /// This function visits each node once and all of its outgoing edges once, so the time
    /// complexity is O(|V| + |E|), where |V| is the number of atoms in the molecule and |E| is the
    /// number of bonds in the molecule.
    fn all_distances_from_source(
        molecule: &Mol,
        source: NodeIndex,
        distances_from_source: &mut FxHashMap<NodeIndex, u16>,
        bfs_queue: &mut VecDeque<(NodeIndex, u16)>,
    ) {
        bfs_queue.push_back((source, 0));
        while let Some((curr_node, dist_to_source)) = bfs_queue.pop_front() {
            distances_from_source.insert(curr_node, dist_to_source);
            // add the neighbors to the queue
            molecule.graph.neighbors(curr_node).for_each(|neigh_idx| {
                // the shortest distance is calculated the first time we encounter a node in a BFS,
                // so we don't consider it if we see it again
                if !distances_from_source.contains_key(&neigh_idx) {
                    bfs_queue.push_back((neigh_idx, dist_to_source + 1));
                }
            })
        }
    }
}

/// A structure to store a submolecule, defined as a subset of the atoms within a molecule, along
/// with a reference to a [`DistanceInfo`] which was calculated for that molecule.
///
/// This can be used to define a subset of the atom indices of the molecule. The distance from a
/// [`DistanceSubMol`] to another [`DistanceSubMol`] can then be calculated on the internal atom
/// indices. For more information, see [`DistanceInfo::set_distance`].
#[derive(Clone, Debug)]
pub struct DistanceSubMol<'a> {
    /// A reference to the distance info for the molecule that this submol refers to
    distance_info: &'a DistanceInfo<'a>,
    /// The atom indices of the molecule that this submol refers to
    mol_indices: FxHashSet<NodeIndex>,
}

impl std::cmp::PartialEq for DistanceSubMol<'_> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.distance_info, other.distance_info)
            && self.mol_indices == other.mol_indices
    }
}

impl Distance for DistanceSubMol<'_> {
    fn distance(&self, other: &Self) -> Of64 {
        assert!(std::ptr::eq(self.distance_info, other.distance_info));
        self.distance_info
            .set_distance(&self.mol_indices, &other.mol_indices)
    }
}

impl<'a> AsRef<DistanceSubMol<'a>> for DistanceSubMol<'a> {
    #[inline]
    fn as_ref(&self) -> &DistanceSubMol<'a> {
        self
    }
}

impl<'a> DistanceSubMol<'a> {
    /// Construct a new instance of [`Self`] from a subgraph isomorphism.
    ///
    /// This function fails if `index_map` does not refer to valid u16 integers.
    pub fn from_distance_info_and_index_map(
        distance_info: &'a DistanceInfo<'a>,
        index_map: &FxHashMap<i32, i32>,
    ) -> Result<Self, TryFromIntError> {
        let mol_indices: FxHashSet<NodeIndex> = index_map
            .values()
            .copied()
            .map(|i| u16::try_from(i).map(NodeIndex::from))
            .collect::<Result<_, _>>()?;

        Ok(Self {
            distance_info,
            mol_indices,
        })
    }

    /// Obtain a list of [`Self`] given the distance info for a molecule and a list of subgraph
    /// isomorphisms corresponding to a specific modification.
    ///
    /// Note that `subgraph_isomorphisms` must live at least as long as `distance_info`.
    pub fn from_distance_info_and_subgraph_isomorphisms(
        distance_info: &'a DistanceInfo<'a>,
        subgraph_isomorphisms: &[FxHashMap<NodeIndex, NodeIndex>],
    ) -> Vec<Self> {
        subgraph_isomorphisms
            .iter()
            .map(|index_conversions| {
                // convert this into a set of node indices within the molecule
                let mol_indices = index_conversions
                    .values()
                    .copied()
                    .collect::<FxHashSet<_>>();

                Self {
                    distance_info,
                    mol_indices,
                }
            })
            .collect::<Vec<_>>()
    }

    /// Choose a representative modification within the cluster.
    ///
    /// This is chosen as the modification which has the least distance to all other modifications
    /// within the cluster.
    pub fn representative<R>(cluster: &[R]) -> Option<R>
    where
        R: Copy + AsRef<DistanceSubMol<'a>>,
    {
        cluster
            .iter()
            .enumerate()
            .map(|(idx, submol)| {
                cluster
                    .iter()
                    .enumerate()
                    .filter_map(|(idx2, submol2)| {
                        // only calculate distance to other mods since distance to self is always 0
                        (idx != idx2).then(|| submol.as_ref().distance(submol2.as_ref()))
                    })
                    .sum::<Of64>()
            })
            .position_min()
            .map(|centroid_idx| cluster[centroid_idx])
    }

    /// Get the atom indices that this submol refers to.
    pub fn mol_indices(&self) -> &FxHashSet<NodeIndex> {
        &self.mol_indices
    }

    /// Get the molecule that this submol is with respect to.
    ///
    /// A submol only makes sense in the context of a molecule. This function returns a reference
    /// to the molecule that was used to define `self`.
    pub fn mol(&self) -> &Mol {
        self.distance_info.molecule
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, MolBuilder};
    use approx::assert_relative_eq;

    #[test]
    fn test_bfs() {
        let mol = MolBuilder::from_smiles("S(C1C(C(C(C(O1)Br)Cl)F)N)P")
            .unwrap()
            .build();
        let sulfur_idx = NodeIndex::new(0);
        let num_atoms = mol.graph.node_count();

        let mut distances_from_sulfur =
            FxHashMap::with_capacity_and_hasher(num_atoms, std::default::Default::default());
        let mut queue = VecDeque::with_capacity(num_atoms);
        DistanceInfo::all_distances_from_source(
            &mol,
            sulfur_idx,
            &mut distances_from_sulfur,
            &mut queue,
        );

        for (node, dist) in distances_from_sulfur {
            let atom_type = mol.graph.node_weight(node).unwrap();
            match *atom_type {
                Atom::S => assert_eq!(dist, 0),
                Atom::P => assert_eq!(dist, 1),
                Atom::O => assert_eq!(dist, 2),
                Atom::N => assert_eq!(dist, 3),
                Atom::Br => assert_eq!(dist, 4),
                Atom::F => assert_eq!(dist, 4),
                Atom::Cl => assert_eq!(dist, 5),
                _ => continue,
            }
        }
    }

    #[test]
    fn test_distances() {
        // first try with a linear molecule
        let mol = MolBuilder::from_smiles("O=CC(CS)NC(=O)C(CS)N")
            .unwrap()
            .build();
        // first cysteine residue, including relevant part of backbone
        let set1 = vec![
            NodeIndex::new(0),
            NodeIndex::new(1),
            NodeIndex::new(2),
            NodeIndex::new(3),
            NodeIndex::new(4),
            NodeIndex::new(5),
        ]
        .into_iter()
        .collect::<FxHashSet<_>>();
        // second cysteine residue, including relevant part of backbone
        let set2 = vec![
            NodeIndex::new(6),
            NodeIndex::new(7),
            NodeIndex::new(8),
            NodeIndex::new(9),
            NodeIndex::new(10),
            NodeIndex::new(11),
        ]
        .into_iter()
        .collect::<FxHashSet<_>>();

        let dinfo = DistanceInfo::new(&mol);
        // make sure asymmetric distances are correct
        assert_relative_eq!(dinfo.asymmetric_set_distance(&set1, &set2), 17. / 6.);
        assert_relative_eq!(dinfo.asymmetric_set_distance(&set2, &set1), 15. / 6.);

        // make sure symmetric set distance is same regardless of order
        assert_relative_eq!(dinfo.set_distance(&set1, &set2).into_inner(), 32. / 12.);
        assert_relative_eq!(dinfo.set_distance(&set2, &set1).into_inner(), 32. / 12.);

        // make sure the set is distance 0 from itself
        assert_relative_eq!(dinfo.set_distance(&set1, &set1).into_inner(), 0.0);
        assert_relative_eq!(dinfo.set_distance(&set2, &set2).into_inner(), 0.0);

        // try with a cyclic molecule
        let mol = MolBuilder::from_smiles("S(C1C(C(C(C(O1)Br)Cl)F)N)P")
            .unwrap()
            .build();

        // [S, Br] -> [3, 3]
        let set1 = vec![NodeIndex::new(0), NodeIndex::new(7)]
            .into_iter()
            .collect::<FxHashSet<_>>();
        // [Cl, F, N] -> [3, 4, 3]
        let set2 = vec![NodeIndex::new(10), NodeIndex::new(8), NodeIndex::new(9)]
            .into_iter()
            .collect::<FxHashSet<_>>();

        let dinfo = DistanceInfo::new(&mol);
        // make sure asymmetric distances are correct
        assert_relative_eq!(dinfo.asymmetric_set_distance(&set1, &set2), 3.);
        assert_relative_eq!(dinfo.asymmetric_set_distance(&set2, &set1), 10. / 3.);

        // make sure symmetric set distance is same regardless of order
        assert_relative_eq!(dinfo.set_distance(&set1, &set2).into_inner(), 19. / 6.);
        assert_relative_eq!(dinfo.set_distance(&set2, &set1).into_inner(), 19. / 6.);

        // make sure the set is distance 0 from itself
        assert_relative_eq!(dinfo.set_distance(&set1, &set1).into_inner(), 0.0);
        assert_relative_eq!(dinfo.set_distance(&set2, &set2).into_inner(), 0.0);
    }

    #[test]
    fn test_representative_point() {
        let mol = MolBuilder::from_smiles("NC(CS)C(=O)NC(CS)C(=O)NC(CS)C(=O)O")
            .unwrap()
            .build();
        let distance_info = DistanceInfo::new(&mol);

        let first_cysteine = vec![
            NodeIndex::new(1),
            NodeIndex::new(2),
            NodeIndex::new(3),
            NodeIndex::new(4),
            NodeIndex::new(5),
            NodeIndex::new(6),
        ]
        .into_iter()
        .collect::<FxHashSet<_>>();

        let first_cysteine = DistanceSubMol {
            distance_info: &distance_info,
            mol_indices: first_cysteine,
        };

        let second_cysteine = vec![
            NodeIndex::new(7),
            NodeIndex::new(8),
            NodeIndex::new(9),
            NodeIndex::new(10),
            NodeIndex::new(11),
            NodeIndex::new(12),
        ]
        .into_iter()
        .collect::<FxHashSet<_>>();

        let second_cysteine = DistanceSubMol {
            distance_info: &distance_info,
            mol_indices: second_cysteine,
        };

        let third_cysteine = vec![
            NodeIndex::new(13),
            NodeIndex::new(14),
            NodeIndex::new(15),
            NodeIndex::new(16),
            NodeIndex::new(17),
            NodeIndex::new(18),
        ]
        .into_iter()
        .collect::<FxHashSet<_>>();

        let third_cysteine = DistanceSubMol {
            distance_info: &distance_info,
            mol_indices: third_cysteine,
        };

        let residues = &[&first_cysteine, &second_cysteine, &third_cysteine];
        let centroid = DistanceSubMol::representative(residues);
        assert_eq!(Some(&second_cysteine), centroid);
    }

    #[test]
    fn test_representative_single_element_cluster() {
        let mol = MolBuilder::from_smiles("NC(CS)C(=O)O").unwrap().build();
        let distance_info = DistanceInfo::new(&mol);

        let first_cysteine = vec![
            NodeIndex::new(1),
            NodeIndex::new(2),
            NodeIndex::new(3),
            NodeIndex::new(4),
            NodeIndex::new(5),
            NodeIndex::new(6),
        ]
        .into_iter()
        .collect::<FxHashSet<_>>();

        let first_cysteine = DistanceSubMol {
            distance_info: &distance_info,
            mol_indices: first_cysteine,
        };

        let representative = DistanceSubMol::representative(&[&first_cysteine]);
        assert_eq!(representative, Some(&first_cysteine));
    }
}
