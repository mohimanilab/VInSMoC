use std::collections::VecDeque;

use petgraph::{graph::DiGraph, visit::EdgeRef, EdgeDirection};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::{Mol, MolGraph, NodeIndex};

/// Graph isomorphism implementations
impl Mol {
    /// Implementation of the RI-DS Subgraph Isomorphism Algorithm
    ///
    /// Outputs a FxHashMap which pertains to a isomorphic matching between some subgraph of a target graph and a
    /// pattern graph. This algorithm, as outlined in its corresponding paper, makes use of the above
    /// preprocessing function to process vertices in an order that reduces the search space in a brute force
    /// isomorphism finder.
    ///
    /// # Citations
    /// \[1\]: [Bonnici, V., Giugno, R., Pulvirenti, A. et al. A subgraph isomorphism algorithm and
    /// its application to biochemical data. BMC Bioinformatics 14, S13 (2013).](https://doi.org/10.1186/1471-2105-14-S7-S13)
    pub fn ri_subgraph_isomorphism(&self, pattern: &Mol) -> Vec<FxHashMap<NodeIndex, NodeIndex>> {
        let mu = pattern.ri_pattern_preprocess();
        let mut parent: Vec<Option<NodeIndex>> = Vec::new();

        // Generate parent vector
        for node_ind in 0..mu.len() {
            for pv_ind in 0..node_ind {
                if pattern.graph.find_edge(mu[pv_ind], mu[node_ind]).is_some() {
                    parent.push(Some(mu[pv_ind]));
                    break;
                }
            }
            if parent.len() < (node_ind + 1) {
                parent.push(None);
            }
        }
        let mut search_tree =
            DiGraph::<(NodeIndex, NodeIndex, Vec<(NodeIndex, NodeIndex)>), (), u16>::default();
        // Finding initial "roots" of the search tree
        for init_node in self.graph.node_indices() {
            if self.h_count(&init_node) >= pattern.h_count(&mu[0])
                && self.graph[init_node] == pattern.graph[mu[0]]
            {
                search_tree.add_node((mu[0], init_node, vec![(mu[0], init_node)]));
            }
        }

        // Iterating over remaining pattern vertices to match
        for mu_ind in 1..mu.len() {
            let curr_mu: NodeIndex = mu[mu_ind];
            // This contains the set of possible search tree parents
            let mut possible_pars: Vec<NodeIndex> = Vec::new();
            for node in search_tree.node_indices() {
                if search_tree[node].0 == mu[mu_ind - 1] {
                    possible_pars.push(node);
                }
            }

            // Iterating over possible self nodes to match with curr_mu
            for par in &possible_pars {
                // This FxHashMap represents the matchings we've made so far (self -> pattern), used for checking conditions 1/4
                let mut cur_mapping: FxHashMap<usize, usize> = FxHashMap::default();
                for matches in &search_tree[*par].2 {
                    cur_mapping.insert(matches.0.index(), matches.1.index());
                }

                // Set of possible self nodes to match with curr_mu
                let possible_self_matches: Vec<_> = if let Some(par_idx) = parent[mu_ind] {
                    self.graph
                        .neighbors(NodeIndex::new(cur_mapping[&par_idx.index()]))
                        .collect()
                } else {
                    self.graph.node_indices().collect()
                };
                for self_node in possible_self_matches {
                    let mut match_found = true;
                    // Condition 1: Make sure self_node hasn't been matched
                    if search_tree[*par].2.iter().any(|(_, x)| *x == self_node) {
                        continue;
                    }

                    // Condition 2: Node labels should be the same
                    // this includes a check on hydrogen compatibility
                    if self.graph[self_node] != pattern.graph[curr_mu]
                        || self.h_count(&self_node) < pattern.h_count(&curr_mu)
                    {
                        continue;
                    }

                    // Condition 3: Degree Condition
                    if self.graph.neighbors(self_node).count()
                        < pattern.graph.neighbors(curr_mu).count()
                    {
                        continue;
                    }

                    // Condition 4: Topological Conditions, Compatibility of Edge Labels
                    for mu_neighbor_edge in pattern.graph.edges(curr_mu) {
                        let mu_neighbor_node = mu_neighbor_edge.target();
                        // do we have a mapping for this pattern graph neighbor node in our current
                        // index mapping path?
                        if let Some((_, self_neighbor_idx)) = search_tree[*par]
                            .2
                            .iter()
                            .find(|(x, _)| *x == mu_neighbor_node)
                        {
                            // does an edge in the original graph exist if we applying the
                            // candidate mapping?
                            if let Some(self_neighbor_edge) =
                                self.graph.find_edge(*self_neighbor_idx, self_node)
                            {
                                // if the label is wrong this can't be a match
                                if self.graph[self_neighbor_edge] != *mu_neighbor_edge.weight() {
                                    match_found = false;
                                    break;
                                }
                                continue;
                            }
                            // no corresponding edge even though we mapped this neighbor, this can't be a match
                            match_found = false;
                            break;
                        } // this neighbor isn't mapped yet, carry on
                    }

                    if match_found {
                        let mut new_path = search_tree[*par].2.clone();
                        new_path.push((curr_mu, self_node));
                        let new_node = search_tree.add_node((curr_mu, self_node, new_path));
                        search_tree.add_edge(*par, new_node, ());
                    }
                }
            }
        }

        // Adding all final matches to iso_matches
        let mut iso_matches: Vec<FxHashMap<NodeIndex, NodeIndex>> = Vec::new();
        for sink in search_tree.externals(EdgeDirection::Outgoing) {
            if search_tree[sink].2.len() == mu.len() {
                let mapping = search_tree[sink]
                    .2
                    .clone()
                    .into_iter()
                    .collect::<FxHashMap<NodeIndex, NodeIndex>>();
                iso_matches.push(mapping);
            }
        }
        iso_matches
    }

    /// Implementation of the GreedyConstraintFirst Function for RI-DS
    ///
    /// Outputs a greedy ordering of the vertices in an input graph for use in the RI-DS
    /// subgraph isomorphism algorithm. Specifically, vertices are iterated through and
    /// inserted into an "order" vector on the basis of a scoring function which involves the
    /// number of visited/unvisited neighbors of each node at a particular step. This serves
    /// to reduce the search space of the actual subgraph isomorphism algorithm.
    ///
    /// # Citations
    /// \[1\]: [Bonnici, V., Giugno, R., Pulvirenti, A. et al. A subgraph isomorphism algorithm and
    /// its application to biochemical data. BMC Bioinformatics 14, S13 (2013).](https://doi.org/10.1186/1471-2105-14-S7-S13)
    pub fn ri_pattern_preprocess(&self) -> Vec<NodeIndex> {
        let mut ord = Vec::new();
        let mut remaining = self.graph.node_indices().collect::<Vec<_>>();

        if self.graph.node_count() == 0 {
            return ord;
        }

        let mu_0 = self
            .graph
            .node_indices()
            .max_by_key(|node| self.graph.neighbors(*node).count())
            .unwrap(); // assumes there is at least one node in molecule

        ord.push(mu_0);
        remaining.retain(|&x| x != mu_0);

        while !remaining.is_empty() {
            let mut cur_best = (0, 0, 0);
            let mut next_vertex = NodeIndex::new(0);
            for node in remaining.iter() {
                let neighbors = self.graph.neighbors(*node).collect::<FxHashSet<_>>();
                let vis_size = ord.iter().filter(|x| neighbors.contains(x)).count();
                let mut neig_size: usize = 0;
                let mut unv_size: usize = 0;

                for prev_node in ord.iter() {
                    if self.graph.find_edge(*prev_node, *node).is_none() {
                        for unvis_node in remaining.iter() {
                            if self.graph.find_edge(*node, *unvis_node).is_some()
                                && self.graph.find_edge(*prev_node, *unvis_node).is_some()
                            {
                                neig_size += 1;
                                break;
                            }
                        }
                    }
                }

                for neighbor in self.graph.neighbors(*node) {
                    if !ord
                        .iter()
                        .any(|vertex| self.graph.find_edge(neighbor, *vertex).is_some())
                        && !ord.iter().any(|vertex| *vertex == neighbor)
                    {
                        unv_size += 1;
                    }
                }

                let score: (usize, usize, usize) = (vis_size, neig_size, unv_size);
                if score >= cur_best {
                    next_vertex = *node;
                    cur_best = score;
                }
            }
            ord.push(next_vertex);
            remaining.retain(|&x| x != next_vertex);
        }
        ord
    }

    /// Implementation of the process outlined in Lemma 5 for the Eppstein Planar Subgraph Isomorphism ALgorithm
    ///
    /// Outputs a partition of the vertices of a given planar graph, paired with resultant induced subgraph.
    /// Specifically, this groups the vertices of an input graph by their distance from an arbitrary vertex,
    /// which we consider as the vertex with index 0 here. This aids in efficiently computing the number of
    /// isomorphs on subgraphs with lower tree-width using dynamic programming, as outlined in the respective paper.
    ///
    /// # Citations
    /// \[1\]: [Eppstein, D. (2002). Subgraph isomorphism in planar graphs and related problems. In Graph
    /// Algorithms and Applications I (pp. 283-309).](https://arxiv.org/abs/cs/9911003v1)
    pub fn eppstein_partition(&self, w: usize) -> Vec<(Vec<usize>, MolGraph)> {
        let node_count = self.graph.node_count();

        // This maps each distance i to a vector of nodes that are at distance i from 0
        let mut sets = vec![Vec::<usize>::new(); node_count];

        // BFS to get distances
        let mut dist = vec![0usize; node_count];
        let mut vis = vec![false; node_count];
        let mut bfs_queue = VecDeque::new();
        bfs_queue.push_back(NodeIndex::new(0));
        vis[0] = true;
        while !bfs_queue.is_empty() {
            let cur_node = bfs_queue.pop_front().unwrap();
            sets[dist[cur_node.index()]].push(cur_node.index());
            for next_node in self.graph.neighbors(cur_node) {
                if !vis[next_node.index()] {
                    dist[next_node.index()] = dist[cur_node.index()] + 1;
                    vis[next_node.index()] = true;
                    bfs_queue.push_back(next_node);
                }
            }
        }

        // This variable contains our final partition of the vertices, pairing each set with its respective subgraph
        let mut partition: Vec<(Vec<usize>, MolGraph)> = Vec::new();
        for (dist_ind, cur_set) in sets.iter().enumerate().take(node_count) {
            if cur_set.is_empty() {
                continue;
            }
            // Create Induced Subgraph
            let induced_subgraph = self.graph.filter_map(
                |nodeind, node| {
                    if (dist_ind..(w + dist_ind)).contains(&dist[nodeind.index()]) {
                        Some(*node)
                    } else {
                        None
                    }
                },
                |_, edge| Some(*edge),
            );
            partition.push((cur_set.clone(), induced_subgraph));
        }
        partition
    }
}

#[cfg(test)]
mod tests {
    use petgraph::{graph::node_index, stable_graph::StableGraph as Graph};
    use rustc_hash::{FxHashMap, FxHashSet};

    use crate::{parsers::parse_smiles, Atom, BondType::*, Mol, NodeIndex};

    #[test]
    fn ri_pattern_preprocess() {
        // Currently just testing ri_pattern_preprocess with the same case as outlined above. Pending further testing once the main method is completed

        let mut backbone = Graph::with_capacity(5, 5);
        backbone.add_node(Atom::C);
        backbone.add_node(Atom::N);
        backbone.add_node(Atom::C);
        backbone.add_node(Atom::O);
        backbone.add_node(Atom::O);

        backbone.add_edge(node_index(0), node_index(1), Single);
        backbone.add_edge(node_index(0), node_index(2), Single);
        backbone.add_edge(node_index(2), node_index(3), Single);
        backbone.add_edge(node_index(2), node_index(4), Double);
        let hydrogen_hash = [
            (node_index(0), 1),
            (node_index(1), 2),
            (node_index(2), 0),
            (node_index(3), 1),
            (node_index(4), 0),
        ]
        .iter()
        .cloned()
        .collect();

        let charge_hash = (0..backbone.node_count())
            .map(|x| (node_index(x), 0))
            .collect();

        let mol = Mol {
            graph: backbone,
            hydrogens: hydrogen_hash,
            charges: charge_hash,
            aromatic_atoms: FxHashSet::default(),
        };

        let init_ord = mol.ri_pattern_preprocess();
        assert_eq!(init_ord, [2.into(), 0.into(), 4.into(), 3.into(), 1.into()]);
    }

    #[test]
    fn ri_subgraph_isomorphism() {
        let mut target = parse_smiles("C(=O)NC#N").unwrap();
        let mut query = parse_smiles("[C][N]").unwrap();
        target.implicit_hs();
        query.implicit_hs();
        // should map carbon in query (idx 0) to carbon at idx 0 in target and carbon at idx 3 in target, since both are connected with a single bond to a nitrogen
        let mapping1 = vec![
            (node_index(0), node_index(0)),
            (node_index(1), node_index(2)),
        ]
        .into_iter()
        .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        let mapping2 = vec![
            (node_index(0), node_index(3)),
            (node_index(1), node_index(2)),
        ]
        .into_iter()
        .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        let solution = target.ri_subgraph_isomorphism(&query);
        assert_eq!(solution.len(), 2);
        assert!(solution.contains(&mapping1));
        assert!(solution.contains(&mapping2));

        let mut target = parse_smiles("C(=O)NC#N").unwrap();
        let mut query = parse_smiles("[C]=[N]").unwrap();
        target.implicit_hs();
        query.implicit_hs();
        // should have no mapping due to bond type
        let solution = target.ri_subgraph_isomorphism(&query);
        assert_eq!(solution.len(), 0);

        let mut target = parse_smiles("C(=O)NC#N").unwrap();
        let mut query = parse_smiles("[CH][N]").unwrap();
        target.implicit_hs();
        query.implicit_hs();
        // should have just one mapping due to hydrogen constraint
        let solution = target.ri_subgraph_isomorphism(&query);
        let mapping = vec![(0.into(), 0.into()), (2.into(), 2.into())]
            .into_iter()
            .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        assert_eq!(solution.len(), 1);
        assert_eq!(solution, vec![mapping]);

        // multiple connected components
        let mut target = parse_smiles("C(=O)NC#N").unwrap();
        let mut query = parse_smiles("[CH].[N]").unwrap();
        target.implicit_hs();
        query.implicit_hs();
        let solution = target.ri_subgraph_isomorphism(&query);
        let mapping1 = vec![(0.into(), 0.into()), (2.into(), 2.into())]
            .into_iter()
            .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        let mapping2 = vec![(0.into(), 0.into()), (2.into(), 4.into())]
            .into_iter()
            .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        assert_eq!(solution.len(), 2);
        assert!(solution.contains(&mapping1));
        assert!(solution.contains(&mapping2));

        let mut target = parse_smiles("C(=O)NC#N").unwrap();
        let mut query = parse_smiles("[C].[N]").unwrap();
        target.implicit_hs();
        query.implicit_hs();
        let solution = target.ri_subgraph_isomorphism(&query);
        let mapping1 = vec![(0.into(), 0.into()), (1.into(), 2.into())]
            .into_iter()
            .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        let mapping2 = vec![(0.into(), 0.into()), (1.into(), 4.into())]
            .into_iter()
            .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        let mapping3 = vec![(0.into(), 3.into()), (1.into(), 2.into())]
            .into_iter()
            .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        let mapping4 = vec![(0.into(), 3.into()), (1.into(), 4.into())]
            .into_iter()
            .collect::<FxHashMap<NodeIndex, NodeIndex>>();
        assert_eq!(solution.len(), 4);
        assert!(solution.contains(&mapping1));
        assert!(solution.contains(&mapping2));
        assert!(solution.contains(&mapping3));
        assert!(solution.contains(&mapping4));
    }

    #[test]
    fn eppstein_partition() {
        let mut backbone = Graph::with_capacity(5, 5);
        backbone.add_node(Atom::C);
        backbone.add_node(Atom::N);
        backbone.add_node(Atom::C);
        backbone.add_node(Atom::O);
        backbone.add_node(Atom::O);

        backbone.add_edge(node_index(0), node_index(1), Single);
        backbone.add_edge(node_index(0), node_index(2), Single);
        backbone.add_edge(node_index(2), node_index(3), Single);
        backbone.add_edge(node_index(2), node_index(4), Double);
        let hydrogen_hash = [
            (node_index(0), 1),
            (node_index(1), 2),
            (node_index(2), 0),
            (node_index(3), 1),
            (node_index(4), 0),
        ]
        .iter()
        .cloned()
        .collect();

        let charge_hash = (0..backbone.node_count())
            .map(|x| (node_index(x), 0))
            .collect();

        let mol = Mol {
            graph: backbone,
            hydrogens: hydrogen_hash,
            charges: charge_hash,
            aromatic_atoms: FxHashSet::default(),
        };

        let part = mol.eppstein_partition(2);
        assert_eq!(part[0].0, vec![0]);
        assert_eq!(part[1].0, vec![2, 1]);
        assert_eq!(part[2].0, vec![4, 3]);
        assert_eq!(part.len(), 3);
    }

    #[test]
    fn ri_subgraph_isomorphism_pubchem() {
        // https://pubchem.ncbi.nlm.nih.gov/compound/6993429
        let target = parse_smiles("CC(C)(C)OC(=O)NC(CC1=CN(C=N1)CC2=CC=CC=C2)C(=O)O").unwrap();

        // benzene
        let pattern = parse_smiles("[C]=1[C]=[C][C]=[C][C]1").unwrap();
        let mut solution = target.ri_subgraph_isomorphism(&pattern);
        // 6 mappings because of cyclic rotations and reversals

        // rotations
        assert_eq!(solution.len(), 6);
        let mapping = vec![
            (1.into(), 16.into()),
            (2.into(), 17.into()),
            (3.into(), 18.into()),
            (4.into(), 19.into()),
            (5.into(), 20.into()),
            (0.into(), 21.into()),
        ]
        .into_iter()
        .collect::<FxHashMap<_, _>>();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());

        let mapping = vec![
            (1.into(), 18.into()),
            (2.into(), 19.into()),
            (3.into(), 20.into()),
            (4.into(), 21.into()),
            (5.into(), 16.into()),
            (0.into(), 17.into()),
        ]
        .into_iter()
        .collect::<FxHashMap<_, _>>();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());

        let mapping = vec![
            (1.into(), 20.into()),
            (2.into(), 21.into()),
            (3.into(), 16.into()),
            (4.into(), 17.into()),
            (5.into(), 18.into()),
            (0.into(), 19.into()),
        ]
        .into_iter()
        .collect::<FxHashMap<_, _>>();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());

        // reversals
        let mapping = vec![
            (0.into(), 16.into()),
            (1.into(), 21.into()),
            (2.into(), 20.into()),
            (3.into(), 19.into()),
            (4.into(), 18.into()),
            (5.into(), 17.into()),
        ]
        .into_iter()
        .collect::<FxHashMap<_, _>>();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());

        let mapping = vec![
            (0.into(), 18.into()),
            (1.into(), 17.into()),
            (2.into(), 16.into()),
            (3.into(), 21.into()),
            (4.into(), 20.into()),
            (5.into(), 19.into()),
        ]
        .into_iter()
        .collect::<FxHashMap<_, _>>();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());

        let mapping = vec![
            (0.into(), 20.into()),
            (1.into(), 19.into()),
            (2.into(), 18.into()),
            (3.into(), 17.into()),
            (4.into(), 16.into()),
            (5.into(), 21.into()),
        ]
        .into_iter()
        .collect::<FxHashMap<_, _>>();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());

        // should have no more solutions now
        assert_eq!(solution, vec![]);

        // double nitrogen ring
        let pattern = parse_smiles("[N]1[C]=[N][C]=[C]1").unwrap();
        let mut solution = target.ri_subgraph_isomorphism(&pattern);
        // this ring doesn't have rotational symmetry, so we should only find one match
        assert_eq!(solution.len(), 1);

        let mapping = vec![
            (0.into(), 12.into()),
            (1.into(), 13.into()),
            (2.into(), 14.into()),
            (3.into(), 10.into()),
            (4.into(), 11.into()),
        ]
        .into_iter()
        .collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        // should have no more solutions now
        assert_eq!(solution, vec![]);

        // all methyl groups
        let mut pattern = parse_smiles("[CH3]").unwrap();
        pattern.implicit_hs();
        let mut solution = target.ri_subgraph_isomorphism(&pattern);
        // we have an isopropyl structure at one end of the molecule
        assert_eq!(solution.len(), 3);
        let mapping = vec![(0.into(), 0.into())].into_iter().collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        let mapping = vec![(0.into(), 2.into())].into_iter().collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        let mapping = vec![(0.into(), 3.into())].into_iter().collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        // should have no more solutions now
        assert_eq!(solution, vec![]);

        // any nitrogen
        let pattern = parse_smiles("[N]").unwrap();
        let mut solution = target.ri_subgraph_isomorphism(&pattern);
        assert_eq!(solution.len(), 3);
        let mapping = vec![(0.into(), 7.into())].into_iter().collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        let mapping = vec![(0.into(), 12.into())].into_iter().collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        let mapping = vec![(0.into(), 14.into())].into_iter().collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        // should have no more solutions now
        assert_eq!(solution, vec![]);

        // just the nitrogen that has a hydrogen
        let mut pattern = parse_smiles("[NH]").unwrap();
        pattern.implicit_hs();
        let mut solution = target.ri_subgraph_isomorphism(&pattern);
        // now there's only one option
        assert_eq!(solution.len(), 1);
        let mapping = vec![(0.into(), 7.into())].into_iter().collect();
        let idx = solution.iter().position(|x| *x == mapping);
        assert!(idx.is_some());
        solution.remove(idx.unwrap());
        // should have no more solutions now
        assert_eq!(solution, vec![]);

        // the two carobyxl-ish groups at the same time
        let pattern = parse_smiles("O=[C][O].O=[C][O]").unwrap();
        let mut solution = target.ri_subgraph_isomorphism(&pattern);
        assert_eq!(solution.len(), 2);
        dbg!(&solution);
        let mapping1 = vec![
            (3.into(), 6.into()),
            (4.into(), 5.into()),
            (5.into(), 4.into()),
            (0.into(), 23.into()),
            (1.into(), 22.into()),
            (2.into(), 24.into()),
        ]
        .into_iter()
        .collect();
        let mapping2 = vec![
            (3.into(), 23.into()),
            (4.into(), 22.into()),
            (5.into(), 24.into()),
            (0.into(), 6.into()),
            (1.into(), 5.into()),
            (2.into(), 4.into()),
        ]
        .into_iter()
        .collect();
        let idx1 = solution.iter().position(|x| *x == mapping1);
        let idx2 = solution.iter().position(|x| *x == mapping2);
        assert!(idx1.is_some());
        assert!(idx2.is_some());
        solution.remove(idx.unwrap());

        let target = parse_smiles("CCCC").unwrap();
        let query = parse_smiles("[C].[O]").unwrap();
        let solution = target.ri_subgraph_isomorphism(&query);
        assert_eq!(solution.len(), 0);

        // test for overlapping multi-CCs
        // if multi-CCs wrong then you'll probably get 4 instead of 2
        let target = parse_smiles("CC").unwrap();
        let query = parse_smiles("[C].[C]").unwrap();
        let solution = target.ri_subgraph_isomorphism(&query);
        assert_eq!(solution.len(), 2);
    }
}
