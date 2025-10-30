//! Contains tooling for cycle detection and tracking.
//!
//! Maintaining cycle information for molecule graphs can be of importance
//! when working with aromatic compounds that might require kekulization.
//! Here we define a number of structures and methods to generate and query
//! said information in a semi-efficient manner (we cannot do this in polynomial
//! time and space since finding all cycles is NP-complete--for proof, note that
//! we could reduce finding a Hamiltonian cycle to generating all cycles--but we
//! attempted to use the most efficient algorithms available). Specifically, we use
//! [Paton's Algorithm](https://dl.acm.org/doi/10.1145/363219.363232) to generate
//! a cycle basis and then use [Gibb's Algorithm](https://dspace.mit.edu/bitstream/handle/1721.1/68106/FTL_R_1982_07.pdf)
//! to reconstruct all cycles from the basis. For those uninterested in the references,
//! the basic idea is as follows:
//! 1. Generate a spanning tree, T, for the graph.
//! 2. For every edge, e not in T, DFS for the singleton cycle in T U {e}. This yields one basis cycle.
//! 3. Once we have all the basis cycles, we combine them in every logical manner to get all the cycles.
//!
//! Once the cycles are generated, they can be printed as vertex sequences, and we provide methods for
//! checking if two vertices in a molecule graph belong to the same cycle.
use itertools::Itertools;
use petgraph::visit::{depth_first_search, DfsEvent, EdgeRef};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;

use super::{EdgeIndex, Mol, NodeIndex};

/// Used to gather and maintain information on cycles in a given molecule.
#[derive(Debug, Clone)]
pub struct CycleInfo<'a> {
    /// A list of the cycles in the molecule given as edge lists.
    ///
    /// Each edge list is sorted for ease of comparison and querying.
    /// We sort using the default `EdgeIndex` comparator, which basically
    /// means that the edges are sorted by their arbitrary index assignment.
    pub(crate) cycles: Vec<Vec<EdgeIndex>>,
    /// A mapping from nodes to indices of the cycle list.
    pub(crate) cycle_memberships: FxHashMap<NodeIndex, Vec<usize>>,
    /// A reference to the molecule being described.
    pub molecule: &'a Mol,
}

impl<'a> CycleInfo<'a> {
    /// Compute and return the cycle info for the given molecule.
    ///
    /// # Example
    /// ```
    /// use molecule::{CycleInfo, MolBuilder};
    /// # use std::error::Error;
    /// # fn t() -> Result<(), Box<dyn Error>> {
    /// // Get the cycle info for naphthalene
    /// let mol = MolBuilder::from_smiles("C1=CC=C2C=CC=CC2=C1")?.build();
    /// assert_eq!(CycleInfo::new(&mol).cycles().len(), 3);
    /// # Ok(())
    /// # }
    /// # assert!(t().is_ok());
    /// ```
    pub fn new(molecule: &'a Mol) -> Self {
        let mut out = Self {
            molecule,
            cycles: Vec::new(),
            cycle_memberships: FxHashMap::default(),
        };

        out.find_cycles();
        out
    }

    /// Getter function to access the cycle list. Returns a reference to a `Vec<Vec<EdgeIndex>>`, where each
    /// element is a list of the edges in the cycle.
    pub fn cycles(&self) -> &Vec<Vec<EdgeIndex>> {
        &self.cycles
    }

    /// Getter function to access the cycle membership mapping.
    ///
    /// When given a node index, the resultant map can return a vector of indices, each of which
    /// corresponds to a cycle contained in the result of calling `cycles()`. These
    /// indices can also be passed to `get_vertices` to reconstruct the cycles as node sequences.
    pub fn cycle_memberships(&self) -> &FxHashMap<NodeIndex, Vec<usize>> {
        &self.cycle_memberships
    }

    /// Checks to see if `u` and `v` belong to the same cycle.
    ///
    /// # Example
    /// ```
    /// use molecule::{CycleInfo, MolBuilder, NodeIndex};
    /// # use std::error::Error;
    /// # fn t() -> Result<(), Box<dyn Error>> {
    /// // In benzene, every atom is cocyclic with every other
    /// let mol = MolBuilder::from_smiles("C1=CC=CC=C1")?.build();
    /// assert!(CycleInfo::new(&mol).cocyclic(NodeIndex::from(0), NodeIndex::from(1)));
    ///
    /// // In butane, no atoms are cocyclic
    /// let mol = MolBuilder::from_smiles("CCCC")?.build();
    /// assert!(!CycleInfo::new(&mol).cocyclic(NodeIndex::from(0), NodeIndex::from(1)));
    /// # Ok(())
    /// # }
    /// # assert!(t().is_ok());
    /// ```
    pub fn cocyclic(&self, u: NodeIndex, v: NodeIndex) -> bool {
        let cycle_list1 = self.cycle_memberships.get(&u);
        let cycle_list2 = self.cycle_memberships.get(&v);

        if cycle_list1.is_none() || cycle_list2.is_none() {
            return false;
        }

        let cycle_list1 = cycle_list1.unwrap();
        let cycle_list2 = cycle_list2.unwrap();

        for ind in cycle_list1 {
            if cycle_list2[..].binary_search(ind).is_ok() {
                return true;
            }
        }

        // No shared cycles found
        false
    }

    /// Given an index, get the cycle at that index and return its vertex list representation.
    ///
    /// Returns [`None`] if the index is out of bounds for the cycle list. We do not impose a canonical
    /// ordering on the vertex cycle (the first vertex and direction are arbitrary). NOTE: This is NOT the
    /// same as just indexing into the cycle list. The cycle list stores cycles as lists of edges. However,
    /// it is often more useful to get the cycle as an ordered list of vertices.
    ///
    /// # Example
    /// ```
    /// use molecule::{CycleInfo, MolBuilder, NodeIndex};
    /// # use std::error::Error;
    /// # fn t() -> Result<(), Box<dyn Error>> {
    /// // In benzene, there should be just one cycle, and it should include all vertices
    /// let mol = MolBuilder::from_smiles("C1=CC=CC=C1")?.build();
    /// let info = CycleInfo::new(&mol);
    ///
    /// // NOTE: In this case we happen to get the cycle 0->5->4->3->2->1->0, but this is somewhat
    /// // arbitrary. What we know is that the cycle should include all indices and that the indices
    /// // should be adjacent.
    /// assert_eq!(info.get_vertices(0), Some(vec![0, 5, 4, 3, 2, 1].into_iter().map(NodeIndex::from).collect()));
    /// assert_eq!(info.get_vertices(1), None);
    /// # Ok(())
    /// # }
    /// # assert!(t().is_ok());
    /// ```
    pub fn get_vertices(&self, cycle_index: usize) -> Option<Vec<NodeIndex>> {
        if cycle_index >= self.cycles.len() {
            return None;
        }

        let e_cycle = &self.cycles[cycle_index];

        // A valid cycle must have at least 3 edges
        debug_assert!(e_cycle.len() >= 3);

        // Since we should be the only ones constructing cycles, assume all indices exist
        debug_assert!(self.molecule.graph.edge_endpoints(e_cycle[0]).is_some());

        let (start, _) = self.molecule.graph.edge_endpoints(e_cycle[0]).unwrap();

        // We will now DFS the graph using only edges from e_cycle. This will reconstruct the cycle in vertex form
        let mut v_cycle = vec![];
        let mut curr_node = start;
        let mut used = FxHashSet::default();
        let mut found: bool;

        loop {
            found = false;
            for edge_ref in self.molecule.graph.edges(curr_node) {
                let edge = edge_ref.id();
                let next = edge_ref.target();

                if e_cycle.binary_search(&edge).is_ok() && !used.contains(&edge) {
                    v_cycle.push(curr_node);
                    used.insert(edge);
                    curr_node = next;
                    found = true;
                    break;
                }
            }

            // This should never happen in a cycle that satisfies our assumptions
            assert!(found);

            if curr_node == start {
                return Some(v_cycle);
            }
        }
    }

    /// Populate our cycle list for the [`Mol`] passed to our uninitialized [`CycleInfo`].
    fn find_cycles(&mut self) {
        let basis = self.cycle_basis();

        self.cycles = CycleInfo::cycle_basis_to_all_cycles(basis);

        // Now compute cycle memberships
        let mut memberships = FxHashMap::<NodeIndex, Vec<usize>>::default();
        for cycle_index in 0..self.cycles.len() {
            for vertex in self.get_vertices(cycle_index).unwrap() {
                memberships
                    .entry(vertex)
                    .or_insert_with(Vec::new)
                    .push(cycle_index);
            }
        }

        self.cycle_memberships = memberships;
    }

    /// Use Paton's algorithm to compute a cycle basis for our associated molecule graph.
    ///
    /// (K. Paton, An algorithm for finding a fundamental set of cycles for an undirected linear graph, Comm. ACM 12 (1969), pp. 514-518.)
    fn cycle_basis(&self) -> Vec<FxHashSet<EdgeIndex>> {
        let mut spanning_tree = self.spanning_tree();
        let mut basis = Vec::new();

        // For every edge, e, not in the spanning tree, T, we find the unique cycle in T U {e}
        for edge in self.molecule.graph.edge_indices() {
            if spanning_tree.contains(&edge) {
                continue;
            }

            spanning_tree.insert(edge);

            // Choose an arbitrary starting vertex for our traversal
            let (v, _) = self.molecule.graph.edge_endpoints(edge).unwrap();

            // Get the cycle
            let cycle_path = self.dfs_single_cycle(v, &spanning_tree);

            // Now convert the cycle vecdeque into an edge hashset and add it to our basis
            let cycle = cycle_path.into_iter().collect::<FxHashSet<_>>();
            basis.push(cycle);

            // Get rid of the non-spanning tree edge for the next iteration
            spanning_tree.remove(&edge);
        }

        basis
    }

    /// DFS the associated molecule graph to compute an arbitrary spanning tree.
    fn spanning_tree(&self) -> FxHashSet<EdgeIndex> {
        let mut tree_edges = FxHashSet::default();
        let starts = self.molecule.graph.node_indices().next();

        depth_first_search(&self.molecule.graph, starts, |event| {
            if let DfsEvent::TreeEdge(u, v) = event {
                tree_edges.insert(self.molecule.graph.find_edge(u, v).unwrap());
            }
        });

        tree_edges
    }

    /// DFS using only the specified edges to find the unique cycle in a graph composed of a spanning tree plus one edge.
    fn dfs_single_cycle(
        &self,
        source: NodeIndex,
        edge_set: &FxHashSet<EdgeIndex>,
    ) -> VecDeque<EdgeIndex> {
        // To do this iteratively, we simulate a callstack by pushing both visit events and exit
        // events for when we're done with a node
        enum Event {
            // If the first arm is (u, v, e), we are visiting u from v via edge e
            Visit(NodeIndex, Option<NodeIndex>, Option<EdgeIndex>),
            // This arm tells us that we have finished a node, so we should pop from the path
            Finish,
        }

        let mut events = vec![Event::Visit(source, None, None)];
        let mut visited = FxHashSet::<NodeIndex>::default();
        let mut curr_path = VecDeque::<EdgeIndex>::new();
        let mut done = false;

        while !events.is_empty() && !done {
            // Retrieve a node. Visit it. Push its neighbors
            let event = events.pop().unwrap();

            match event {
                Event::Visit(curr_node, prev, e) => {
                    // If we're visiting for the first time, mark as visited, push an exit
                    // event, push neighbors, and add edge to the path

                    // Mark visited
                    visited.insert(curr_node);

                    // Remember to exit
                    events.push(Event::Finish);

                    // Add edge to the path
                    if let Some(incoming_edge) = e {
                        curr_path.push_back(incoming_edge);
                    }

                    // Handle all the neighbors
                    for edge_ref in self.molecule.graph.edges(curr_node) {
                        let edge = edge_ref.id();
                        let next = edge_ref.target();

                        // Restrict available edges to those in the subgraph
                        if !edge_set.contains(&edge) {
                            continue;
                        }

                        if !visited.contains(&next) {
                            // The successor is unvisited. Create an event to visit it
                            events.push(Event::Visit(next, Some(curr_node), Some(edge)));
                        } else if prev.is_some() && next != prev.unwrap() {
                            // If we are revisiting a node that is not our predecessor,
                            // then we found a cycle!
                            curr_path.push_back(edge);
                            done = true;
                            break;
                        }
                    }
                }
                Event::Finish => {
                    // If we're leaving, pop the edge that brought us here
                    curr_path.pop_back();
                }
            }
        }

        // Now we check if the edge stack is a perfect loop or a "lollipop" (where
        // a lollipop is a cycle connected to a path)
        let last = curr_path[curr_path.len() - 1];
        if self.edge_intersection(last, curr_path[0]).is_none() {
            // The last edge is not coincident with the first, so we trim the lollipop
            while self.edge_intersection(last, curr_path[0]).is_none() {
                curr_path.pop_front();
            }

            // Now pop the last one (three edges met at the final node)
            curr_path.pop_front();
        }

        curr_path
    }

    // Find the vertex that two edges are coincident upon. Return None if one can't be found
    fn edge_intersection(&self, e1: EdgeIndex, e2: EdgeIndex) -> Option<NodeIndex> {
        if self.molecule.graph.edge_endpoints(e1).is_none()
            || self.molecule.graph.edge_endpoints(e2).is_none()
        {
            return None;
        }

        let (v1, v2) = self.molecule.graph.edge_endpoints(e1).unwrap();
        let (v3, v4) = self.molecule.graph.edge_endpoints(e2).unwrap();

        if v1 == v3 || v1 == v4 {
            Some(v1)
        } else if v2 == v3 || v2 == v4 {
            Some(v2)
        } else {
            None
        }
    }

    /// Use Gibb's Algorithm (page 14 in the linked paper) to reconstruct all cycles from the cycle basis.
    ///
    /// https://dspace.mit.edu/bitstream/handle/1721.1/68106/FTL_R_1982_07.pdf
    fn cycle_basis_to_all_cycles(basis: Vec<FxHashSet<EdgeIndex>>) -> Vec<Vec<EdgeIndex>> {
        if basis.is_empty() {
            return vec![];
        }

        // The basis cycle products that produce single cycles
        let mut good_cycles: Vec<FxHashSet<EdgeIndex>> = vec![];

        // All cycle products, including unions of two disjoint cycles
        let mut cycle_products: Vec<FxHashSet<EdgeIndex>> = vec![];

        // The new, good cycles found in the most recent iteration
        let mut new_good_candidates: Vec<FxHashSet<EdgeIndex>> = vec![];

        // The candidates that are not single cycles but are kept to make other cycle products later
        let mut new_bad_candidates: Vec<FxHashSet<EdgeIndex>> = vec![];

        good_cycles.push(basis[0].clone());
        cycle_products.push(basis[0].clone());

        for base_cycle in basis.iter().skip(1) {
            for cycle in cycle_products.iter() {
                let intersection = cycle.intersection(base_cycle);
                let xor = cycle
                    .symmetric_difference(base_cycle)
                    .into_iter()
                    .cloned()
                    .collect::<FxHashSet<_>>();

                // If intersection is non-empty
                if intersection.into_iter().next().is_some() {
                    new_good_candidates.push(xor);
                } else {
                    new_bad_candidates.push(xor);
                }
            }

            let mut supersets: Vec<usize> = vec![];
            for (ind, v) in new_good_candidates.iter().enumerate() {
                for u in new_good_candidates.iter() {
                    if u != v && u.is_subset(v) {
                        supersets.push(ind);
                    }
                }
            }

            supersets.dedup();

            // Retain the cycles that don't appear in supersets. Insert the rest into R*
            for &ind in supersets.iter().rev() {
                new_bad_candidates.push(new_good_candidates.remove(ind));
            }

            good_cycles.extend_from_slice(&new_good_candidates);
            good_cycles.push(base_cycle.clone());

            cycle_products.extend_from_slice(&new_good_candidates);
            cycle_products.extend_from_slice(&new_bad_candidates);
            cycle_products.push(base_cycle.clone());

            new_good_candidates = vec![];
            new_bad_candidates = vec![];
        }

        // Now output all the cycles in canonical form
        good_cycles
            .iter()
            .map(|e_set| {
                let mut e_vec = e_set.iter().cloned().collect_vec();

                // Sort the edges in the cycle to impose a canonical ordering
                e_vec.sort_unstable();
                e_vec
            })
            .collect_vec()
    }
}

#[cfg(all(test, debug_assertions))]
/// Simple brute force implementations for testing
impl<'a> CycleInfo<'a> {
    /// Compute cycle info via brute force for debugging purposes.
    pub fn get_cycle_info_slow(molecule: &'a Mol) -> Self {
        let mut out = Self {
            cycles: Vec::new(),
            cycle_memberships: FxHashMap::default(),
            molecule,
        };

        out.find_cycles_slow();
        out
    }

    /// Populate our cycle list using a brute force algorithm.
    ///
    /// Only for debugging purposes.
    fn find_cycles_slow(&mut self) {
        let mut e_path = vec![];
        let mut v_path = vec![];

        // Collect nodes so that there's no confusion about (im)mutably borrowing self.
        let nodes = self.molecule.graph.node_indices().collect::<Vec<_>>();
        for node in nodes {
            let mut visit_set = FxHashSet::default();
            self.explore_cycles(node, node, &mut visit_set, &mut e_path, &mut v_path);
        }
    }

    /// Explore all possible paths starting from `source_node`. If we find a cycle, add it to the cycle list.
    fn explore_cycles(
        &mut self,
        source_node: NodeIndex,
        curr_node: NodeIndex,
        visited: &mut FxHashSet<NodeIndex>,
        e_path: &mut Vec<EdgeIndex>,
        v_path: &mut Vec<NodeIndex>,
    ) {
        visited.insert(curr_node);
        v_path.push(curr_node);

        // Collect everything before the loop to avoid issues with mutably borrowing self while immutably referenced
        let edges = self
            .molecule
            .graph
            .edges(curr_node)
            .map(|edge_ref| (edge_ref.id(), edge_ref.target()))
            .collect::<Vec<_>>();

        for (edge, next) in edges {
            if !visited.contains(&next) {
                e_path.push(edge);
                self.explore_cycles(source_node, next, visited, e_path, v_path);
                e_path.pop();
            } else if next == source_node && e_path.len() > 1 {
                e_path.push(edge);
                self.insert_cycle(e_path, v_path);
                e_path.pop();
            }
        }

        v_path.pop();
        visited.remove(&curr_node);
    }

    /// Insert a cycle into the cycle list.
    ///
    /// `e_path` and `v_path` both represent the cycle, but `e_path` holds its edge sequence whereas
    /// `v_path` stores a vertex sequence.
    fn insert_cycle(&mut self, e_path: &[EdgeIndex], v_path: &[NodeIndex]) {
        let mut cycle = e_path.to_vec();

        // Sort the cycle edges to impose a canonical ordering
        cycle.sort_unstable();

        // Check if this cycle is new
        for existing_cycle in &self.cycles {
            if existing_cycle.eq(&cycle) {
                return;
            }
        }

        // This cycle didn't match any, so it must be new
        self.cycles.push(cycle.to_vec());
        let new_index = self.cycles.len() - 1;

        // Note cycle membership for every vertex
        for vertex in v_path {
            self.cycle_memberships
                .entry(*vertex)
                .or_insert_with(Vec::new)
                .push(new_index);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    use crate::{parsers::parse_smiles, MolBuilder};

    use super::*;

    // For testing purposes, checks if two cycle info structs for the same molecule contain the same cycles
    fn cycle_info_matches(c1: CycleInfo, c2: CycleInfo) -> bool {
        let c_list_1 = c1.cycles();
        let c_list_2 = c2.cycles();

        // To get away with only matching all cycles in c1, we need to verify that the lists have
        // the same sizes and no duplicates. For now, we will only test for the former since the
        // latter should not be possible with either implementation.

        if c_list_1.len() != c_list_2.len() {
            return false;
        }

        let mut found;
        // For each cycle in c1, we should be able to find a matching one in c2
        for cycle_1 in c_list_1 {
            found = false;
            for cycle_2 in c_list_2 {
                if cycle_1.eq(cycle_2) {
                    found = true;
                    break;
                }
            }

            if !found {
                return false;
            }
        }

        true
    }

    #[test]
    fn count_cycles() {
        // https://pubchem.ncbi.nlm.nih.gov/compound/Cubane#section=Canonical-SMILES
        let mol = parse_smiles("C12C3C4C1C5C2C3C45").unwrap();
        // https://mathworld.wolfram.com/HypercubeGraph.html
        assert_eq!(CycleInfo::new(&mol).cycles.len(), 28);
        assert_eq!(CycleInfo::get_cycle_info_slow(&mol).cycles.len(), 28);

        let mol = parse_smiles("CCCC").unwrap();
        assert_eq!(CycleInfo::new(&mol).cycles.len(), 0);
        assert_eq!(CycleInfo::get_cycle_info_slow(&mol).cycles.len(), 0);

        // https://pubchem.ncbi.nlm.nih.gov/compound/benzene#section=InChI
        let mol = parse_smiles("C1=CC=CC=C1").unwrap();
        assert_eq!(CycleInfo::new(&mol).cycles.len(), 1);
        assert_eq!(CycleInfo::get_cycle_info_slow(&mol).cycles.len(), 1);

        // https://pubchem.ncbi.nlm.nih.gov/compound/naphthalene
        let mol = parse_smiles("C1=CC=C2C=CC=CC2=C1").unwrap();
        assert_eq!(CycleInfo::new(&mol).cycles.len(), 3);
        assert_eq!(CycleInfo::get_cycle_info_slow(&mol).cycles.len(), 3);
    }

    #[test]
    fn full_test_against_brute_force() {
        let smi_file = BufReader::new(File::open("test_files/test_compounds.smi").unwrap());
        let lines = smi_file.lines().collect::<Vec<_>>();

        for line in lines.iter() {
            let smiles = line.as_ref().unwrap();
            let mol = MolBuilder::from_smiles(smiles).unwrap().build();
            let c1 = CycleInfo::new(&mol);
            let c2 = CycleInfo::get_cycle_info_slow(&mol);

            assert!(cycle_info_matches(c1, c2));
        }
    }
}
