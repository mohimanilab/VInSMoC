use std::{
    convert::{TryFrom, TryInto},
    error::Error,
    fmt,
    io::Write,
};

use crate::{
    atoms::{Atom, OrganicAtom},
    bonds::BondType,
};

use constants::MASS_ELECTRON;
use petgraph::{
    algo::{dijkstra, tarjan_scc},
    graph::{EdgeIndex as PgEdgeIndex, NodeIndex as PgNodeIndex},
    stable_graph::StableUnGraph,
    visit::{EdgeRef, IntoEdgeReferences},
};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

use rand::Rng;

mod builder;
mod cycles;
mod distances;
mod fpt;
mod iso;

pub use builder::{BuildMolError, MolBuilder};
pub use cycles::CycleInfo;
pub use distances::{DistanceInfo, DistanceSubMol};

#[cfg(feature = "inchi")]
mod inchi;

#[cfg(feature = "draw")]
mod draw;
#[cfg(feature = "draw")]
pub use draw::Point;

pub type NodeIndex = PgNodeIndex<u16>;
pub type EdgeIndex = PgEdgeIndex<u16>;
pub type MolGraph = StableUnGraph<Atom, BondType, u16>;

/// Indicates that an atom has an impossible valence.
#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct ValenceError {
    atom_idx: usize,
    atom_type: OrganicAtom,
}

impl ValenceError {
    pub fn new(atom_idx: usize, atom_type: OrganicAtom) -> Self {
        Self {
            atom_idx,
            atom_type,
        }
    }
}

impl fmt::Display for ValenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid valence for atom {} with index {}",
            self.atom_type, self.atom_idx
        )
    }
}

impl Error for ValenceError {}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct Mol {
    /// The graph representation of the molecule's explicit bonds
    pub graph: MolGraph,
    /// The implicit hydrogens in the molecule by node
    pub(crate) hydrogens: FxHashMap<NodeIndex, u8>,
    /// The charge of each atom in the molecule by node
    pub(crate) charges: FxHashMap<NodeIndex, i8>,
    /// A set of atom nodes that are aromatic
    pub(crate) aromatic_atoms: FxHashSet<NodeIndex>,
}

impl Mol {
    /// The total bond weight of the provided atom index.
    ///
    /// This bond weight is determined by summing over all coincident bonds, assigning weights to
    /// bonds according to their multiplicity. Aromatic bonds are given a bond weight of 1, since
    /// this simplifies the calculation of unfulfilled valence.
    ///
    /// This does not include any implicit hydrogens.
    pub fn total_bonds(&self, atom_idx: NodeIndex) -> u8 {
        self.graph
            .edges(atom_idx)
            .map(|x| match x.weight() {
                BondType::Single => 1,
                BondType::Double => 2,
                BondType::Triple => 3,
                BondType::Aromatic => 1,
            })
            .sum()
    }

    /// Checks the aromaticity of the provided atom index.
    ///
    /// The aromaticity is determined by checking the aromaticity of all the incident bonds
    /// and if any of them are aromatic, the atom is also deemed aromatic. This is based on
    /// how bond aromaticity is determined while parsing.
    ///
    /// NOTE: The assumption that aromatic bonds imply aromatic atoms at both endpoints
    /// is incorrect and will be fixed in the future. This function is a placeholder
    /// in the meantime.
    pub fn is_aromatic(&self, atom_idx: NodeIndex) -> bool {
        self.graph
            .edges(atom_idx)
            .any(|x| matches!(x.weight(), BondType::Aromatic))
    }

    /// The allocated valence of the provided atom index.
    pub fn allocated_valence(&self, atom_idx: NodeIndex) -> i8 {
        (self.total_bonds(atom_idx) + self.h_count(&atom_idx)) as i8 - self.charge(&atom_idx)
    }

    /// The allowed valence of the provided atom index.
    ///
    /// The allowed valence of the provided atom is extracted if it is an organic atom.
    /// Otherwise, since allowed valence isn't tracked for inorganic atoms, the [`None`] option is returned.
    pub fn allowed_valence(&self, atom_idx: NodeIndex) -> Option<u8> {
        let atom: Result<OrganicAtom, _> = self.graph[atom_idx].try_into();

        atom.ok()
            .and_then(|atom| atom.allowed_valence(self.allocated_valence(atom_idx)))
    }

    /// The standard mass of this molecule.
    pub fn std_mass(&self) -> f64 {
        let non_hydrogens: f64 = self
            .graph
            .node_indices()
            .map(|x| self.graph[x].std_mass())
            .sum();
        let h_mass: f64 = (Atom::H).std_mass();
        let hydrogens: u32 = self.hydrogens.values().map(|x: &u8| *x as u32).sum();

        non_hydrogens + (hydrogens as f64) * h_mass
    }

    /// The exact mass of this molecule.
    pub fn exact_mass(&self) -> f64 {
        let non_hydrogens: f64 = self
            .graph
            .node_indices()
            .map(|x| self.graph[x].exact_mass())
            .sum();
        let h_mass: f64 = (Atom::H).exact_mass();
        let hydrogens: u32 = self.hydrogens.values().map(|x: &u8| *x as u32).sum();

        non_hydrogens + (hydrogens as f64) * h_mass
    }

    /// The number of explicit hydrogens on a particular atom
    ///
    /// The explicit hydrogens are found by iterating through neighbors and checking
    /// if they are the hydrogen atom.
    pub fn explicit_h_count(&self, node: NodeIndex) -> u8 {
        self.graph
            .neighbors(node)
            .filter(|atom_idx| matches!(self.graph[*atom_idx], Atom::H))
            .count() as u8
    }

    /// The number of implicit hydrogens on a particular atom.
    pub fn h_count(&self, node: &NodeIndex) -> u8 {
        *self.hydrogens.get(node).unwrap_or(&0)
    }

    /// A mutable reference to the hydrogen count of a particular atom.
    ///
    /// Modifying this count through the mutable reference changes the underlying hydrogen account.
    pub fn h_entry(&mut self, node: NodeIndex) -> &mut u8 {
        self.hydrogens.entry(node).or_insert(0)
    }

    /// The charge of a particular atom
    pub fn charge(&self, node: &NodeIndex) -> i8 {
        *self.charges.get(node).unwrap_or(&0)
    }

    /// A mutable reference to the charge of a particular atom.
    ///
    /// Modifying this charge through the mutable reference changes the underlying charge value.
    pub fn charge_entry(&mut self, node: NodeIndex) -> &mut i8 {
        self.charges.entry(node).or_insert(0)
    }

    /// The aromaticity of a particular atom
    pub fn aromaticity(&self, node: &NodeIndex) -> bool {
        self.aromatic_atoms.contains(node)
    }

    /// Set the aromaticity of a particular atom
    pub fn set_aromaticity(&mut self, node: &NodeIndex, aromaticity: bool) {
        if aromaticity {
            self.aromatic_atoms.insert(*node);
        } else {
            self.aromatic_atoms.remove(node);
        }
    }

    /// Removes an atom from the molecule.
    ///
    /// If the atom node didn't exist, then will return [`None`].
    /// Otherwise will return the atom, hydrogen count, charge, and aromaticity of the removed node.
    pub fn remove_node(&mut self, node: NodeIndex) -> Option<(Atom, u8, i8, bool)> {
        let atom = self.graph.remove_node(node)?;
        let hcount = self.hydrogens.remove(&node).unwrap_or(0);
        let charge = self.charges.remove(&node).unwrap_or(0);
        let aromatic = self.aromatic_atoms.remove(&node);

        Some((atom, hcount, charge, aromatic))
    }

    /// Adds an atom to the molecule.
    ///
    /// Returns the index of the added atom.
    pub fn add_node(&mut self, atom: Atom, hcount: u8, charge: i8, aromatic: bool) -> NodeIndex {
        let node_idx = self.graph.add_node(atom);
        if hcount != 0 {
            self.hydrogens.insert(node_idx, hcount);
        }
        if charge != 0 {
            self.charges.insert(node_idx, charge);
        }
        if aromatic {
            self.aromatic_atoms.insert(node_idx);
        }

        node_idx
    }

    /// Clean up storage by removing 0 hydrogen counts and charges.
    ///
    /// This operates in place and does not affect the molecular graph.
    pub fn compact(&mut self) {
        let zero_keys = self
            .hydrogens
            .iter()
            .filter_map(|(node, count)| if *count != 0 { None } else { Some(*node) })
            .collect::<Vec<_>>();
        for key in zero_keys {
            self.hydrogens.remove(&key);
        }
        self.hydrogens.shrink_to_fit();
        let zero_keys = self
            .charges
            .iter()
            .filter_map(|(node, charge)| if *charge != 0 { None } else { Some(*node) })
            .collect::<Vec<_>>();
        for key in zero_keys {
            self.charges.remove(&key);
        }
        self.charges.shrink_to_fit();
        self.aromatic_atoms.shrink_to_fit();
    }

    pub fn extend_charges<I: IntoIterator<Item = (NodeIndex, i8)>>(&mut self, charges: I) {
        self.charges.extend(charges);
    }

    pub fn extend_hydrogens<I: IntoIterator<Item = (NodeIndex, u8)>>(&mut self, hydrogens: I) {
        self.hydrogens.extend(hydrogens);
    }

    pub fn extend_aromatic_atoms<I: IntoIterator<Item = NodeIndex>>(&mut self, aromatic_atoms: I) {
        self.aromatic_atoms.extend(aromatic_atoms);
    }

    pub fn hydrogens(&self) -> &FxHashMap<NodeIndex, u8> {
        &self.hydrogens
    }

    pub fn charges(&self) -> &FxHashMap<NodeIndex, i8> {
        &self.charges
    }

    pub fn aromatic_atoms(&self) -> &FxHashSet<NodeIndex> {
        &self.aromatic_atoms
    }

    /// The exact mass of the provided node, including all implicit hydrogens.
    ///
    /// This also accounts for any charges on the provided node by using [`MASS_ELECTRON`].
    pub fn node_exact_mass(&self, node: NodeIndex) -> f64 {
        self.graph[node].exact_mass() + (self.h_count(&node) as f64) * Atom::H.exact_mass()
            - (self.charge(&node) as f64) * MASS_ELECTRON
    }

    /// Retrieve the total formal charge of this molecule.
    pub fn formal_charge(&self) -> i8 {
        self.charges.values().sum()
    }

    /// Represents all hydrogens in the molecule as nodes in the molecular graph.
    ///
    /// By default all hydrogens are stored as extra pieces of information and each node in the
    /// underlying molecular graph is a heavy atom. This function adds all the hydrogens as
    /// explicit nodes instead of implicit information attached to heavy atoms. Note, this will
    /// cause the Mol struct to use a lot more memory and will make most graph operations a lot
    /// slower.
    ///
    /// This operation happens in place.
    pub fn explicit_hs(&mut self) {
        // add in hydrogen nodes
        for (heavy_atom_idx, h_count) in self.hydrogens.iter() {
            for _ in 0..*h_count {
                let h_node_idx = self.graph.add_node(Atom::H);
                self.graph
                    .add_edge(*heavy_atom_idx, h_node_idx, BondType::Single);
            }
        }

        self.hydrogens = FxHashMap::default();
    }

    /// Represents all hydrogens in the molecule as entries in a hydrogen count map.
    ///
    /// This removes all nodes that are hydrogen atoms and instead stores them as extra pieces of
    /// information attached to the heavy atoms. This function does not change the semantic value
    /// of the molecule, it just modifies how the same molecule is represented. Most algorithms
    /// that operate on molecules assume that hydrogens are implicit. This ends up being more
    /// efficient generally since this greatly reduces the number of nodes in the molecular graph.
    ///
    /// This operation happens in place.
    ///
    /// # Panics
    /// If the molecule existed of only hydrogens this function will panic with an informative
    /// message.
    pub fn implicit_hs(&mut self) {
        let mut to_remove = Vec::new();
        for node_idx in self.graph.node_indices() {
            if self.graph[node_idx] != Atom::H {
                continue;
            }

            // add to hydrogen count of node neighboring this hydrogen
            for neigh_idx in self.graph.neighbors(node_idx) {
                *self.hydrogens.entry(neigh_idx).or_insert(0) += 1;
            }

            // delete the hydrogen node
            // since we use a stable graph, this doesn't change any indices
            to_remove.push(node_idx);
        }

        if to_remove.len() == self.graph.node_count() {
            panic!("attempted to remove hydrogens in a molecule consisting of only hydrogens");
        }

        for node_idx in to_remove {
            self.graph.remove_node(node_idx);
            self.hydrogens.remove(&node_idx);
            self.charges.remove(&node_idx);
        }
    }

    /// Adds hydrogens to heavy atoms to complete their valence.
    ///
    /// For all organic heavy atoms in this molecule this adds hydrogens to complete their valence.
    /// This takes atomic charge into account. This function will not remove extra hydrogens, and
    /// if it believes it needs to it will return a [`ValenceError`].
    /// See [`OrganicAtom::remaining_valence`](crate::OrganicAtom::remaining_valence) for more
    /// information about valence calculation.
    ///
    /// This calls [`Mol::implicit_hs()`] internally.
    /// This operation happens in place.
    pub fn infer_hs(&mut self) -> Result<(), ValenceError> {
        // TODO: remove when we have a kekulization function
        assert!(!self
            .graph
            .edge_indices()
            .any(|x| self.graph[x] == BondType::Aromatic));

        self.implicit_hs();

        for node_idx in self.graph.node_indices() {
            let atom: Result<OrganicAtom, _> = self.graph[node_idx].try_into();

            // keep going if this atom isn't organic
            if let Ok(atom) = atom {
                let allocated_valence = self.allocated_valence(node_idx);
                *self.hydrogens.entry(node_idx).or_insert(0) += atom
                    .remaining_valence(allocated_valence)
                    .ok_or(ValenceError {
                        atom_idx: node_idx.index(),
                        atom_type: atom,
                    })?;
            }
        }

        self.compact();

        Ok(())
    }

    /// Checks if all organic heavy atoms have their valence properly filled.
    ///
    /// For all organic heavy atoms in this molecule this checks if hydrogens complete their valence.
    /// This takes atomic charge into account. This function will not remove extra hydrogens, and
    /// if it believes it needs to it will return a [`ValenceError`]. Also, if any any atoms have
    /// unallocated valence this will return a [`ValenceError`].
    ///
    /// See [`OrganicAtom::remaining_valence`](crate::OrganicAtom::remaining_valence) for more
    /// information about valence calculation.
    ///
    /// This calls [`Mol::implicit_hs()`] internally.
    /// This operation happens in place.
    pub fn check_valence(&mut self) -> Result<(), ValenceError> {
        assert!(!self
            .graph
            .edge_indices()
            .any(|x| self.graph[x] == BondType::Aromatic));

        self.implicit_hs();

        for node_idx in self.graph.node_indices() {
            let atom: Result<OrganicAtom, _> = self.graph[node_idx].try_into();

            // keep going if this atom isn't organic
            if let Ok(atom) = atom {
                let allocated_valence = self.allocated_valence(node_idx);
                let valence = atom
                    .remaining_valence(allocated_valence)
                    .ok_or(ValenceError {
                        atom_idx: node_idx.index(),
                        atom_type: atom,
                    })?;
                if valence != 0 {
                    return Err(ValenceError {
                        atom_idx: node_idx.index(),
                        atom_type: atom,
                    });
                }
            }
        }

        Ok(())
    }

    /// Extracts a sub-molecule containing only the provided nodes.
    ///
    /// All edges whose endpoints are not both in the provided node set will be removed. This does
    /// not change any indices on its own. However, unless the provided node indices are contiguous
    /// the resulting molecule will have non-contiguous node indices. To correct for this you can
    /// use [`reset_indices`](Mol::reset_indices) on the output of this function.
    pub fn to_submol(&self, nodes: FxHashSet<NodeIndex>) -> Self {
        let mut graph = self.graph.clone();
        graph.retain_nodes(|_, x| nodes.contains(&x));
        let mut hydrogens = self.hydrogens.clone();
        hydrogens.retain(|x, _| nodes.contains(x));
        let mut charges = self.charges.clone();
        charges.retain(|x, _| nodes.contains(x));
        let mut aromatic_atoms = self.aromatic_atoms.clone();
        aromatic_atoms.retain(|x| nodes.contains(x));

        Self {
            graph,
            hydrogens,
            charges,
            aromatic_atoms,
        }
    }

    /// Makes indices in this molecule go from 0 to the number of atoms.
    ///
    /// Having contiguous node indices is necessary for many applications. This function takes
    /// ownership of the molecule, sets all of its node indices to range from 0 to the number of
    /// atoms, and returns the updated molecule and a map that can convert from the new contiguous
    /// indices to the original indices.
    pub fn reset_indices(self) -> (Self, FxHashMap<NodeIndex, NodeIndex>) {
        let mut index_map = FxHashMap::default();
        let mut graph = MolGraph::with_capacity(self.graph.node_count(), self.graph.edge_count());
        let mut hydrogens = FxHashMap::default();
        let mut charges = FxHashMap::default();
        let mut aromatic_atoms = FxHashSet::default();
        for old_idx in self.graph.node_indices() {
            let new_idx = graph.add_node(self.graph[old_idx]);
            if self.h_count(&old_idx) != 0 {
                hydrogens.insert(new_idx, self.hydrogens[&old_idx]);
            }
            if self.charge(&old_idx) != 0 {
                charges.insert(new_idx, self.charges[&old_idx]);
            }
            if self.aromaticity(&old_idx) {
                aromatic_atoms.insert(new_idx);
            }
            index_map.insert(old_idx, new_idx);
        }

        for edge in self.graph.edge_references() {
            graph.add_edge(
                index_map[&edge.source()],
                index_map[&edge.target()],
                *edge.weight(),
            );
        }

        let mol = Self {
            graph,
            hydrogens,
            charges,
            aromatic_atoms,
        };

        let index_map = mol
            .graph
            .node_indices()
            .zip(self.graph.node_indices())
            .collect();
        (mol, index_map)
    }

    /// Converts this into a V3000 MDL MOL file.
    ///
    /// This accepts anything that implements the [`std::io::Write`] trait. You can pass it a file
    /// or a string easily this way. The most common usage will probably be as follows.
    /// ```rust
    /// # use molecule::{Mol, MolBuilder};
    /// # use std::error::Error;
    /// # fn t() -> Result<(), Box<dyn Error>> {
    /// let mol = MolBuilder::from_smiles("C")?.build();
    /// let mut mol_file = std::fs::File::create("/tmp/mymol.mol")?;
    /// mol.to_mol_file(&mut mol_file)?;
    /// # Ok(())
    /// # }
    /// # assert!(t().is_ok());
    /// ```
    ///
    /// See [the specification](http://c4.cabrillo.edu/404/ctfile.pdf) for more information about
    /// the format. For our molecules every explicit node in the molecular graph gets its own atom
    /// line. If there are any implicit hydrogens for an atom that is encoded in the HCOUNT field
    /// of that atom. Any non-zero charges are encoded in the CHG field.
    pub fn to_mol_file<W: Write>(&self, out: &mut W) -> std::io::Result<()> {
        // writes out the 3-line header
        writeln!(out)?;
        writeln!(out, "Molecule {}", env!("CARGO_PKG_VERSION"))?;
        writeln!(out)?;

        // some magic line that no one ever looks at
        writeln!(out, "  0  0  0     0  0            999 V3000")?;

        let pre = "M  V30 ";

        // start actual MOL and write out counts
        writeln!(out, "{}BEGIN CTAB", pre)?;
        writeln!(
            out,
            "{}COUNTS {} {} 0 0 0",
            pre,
            self.graph.node_count(),
            self.graph.edge_count()
        )?;

        // atoms block
        writeln!(out, "{}BEGIN ATOM", pre)?;
        let mut new_idx_map = FxHashMap::<NodeIndex, usize>::default();
        for (mol_idx, node_idx) in self.graph.node_indices().enumerate() {
            // for some reason MOL files are always 1-indexed
            let mol_idx = mol_idx + 1;
            // also, we need these new indices when outputting bonds
            new_idx_map.insert(node_idx, mol_idx);

            write!(out, "{}{} {} 0 0 0 0", pre, mol_idx, self.graph[node_idx])?;

            // only write out CHG and HCOUNT if they're nonzero
            if self.charge(&node_idx) != 0 {
                write!(out, " CHG={}", self.charge(&node_idx))?;
            }
            if self.h_count(&node_idx) != 0 {
                write!(out, " HCOUNT={}", self.h_count(&node_idx))?;
            }
            writeln!(out)?;
        }
        writeln!(out, "{}END ATOM", pre)?;

        // bonds block
        writeln!(out, "{}BEGIN BOND", pre)?;
        for (mol_idx, edge_idx) in self.graph.edge_indices().enumerate() {
            let mol_idx = mol_idx + 1;
            let bt = match self.graph[edge_idx] {
                BondType::Single => 1,
                BondType::Double => 2,
                BondType::Triple => 3,
                BondType::Aromatic => 4,
            };

            let (start, end) = self.graph.edge_endpoints(edge_idx).unwrap();
            let start = new_idx_map[&start];
            let end = new_idx_map[&end];

            writeln!(out, "{}{} {} {} {}", pre, mol_idx, bt, start, end)?;
        }
        writeln!(out, "{}END BOND", pre)?;

        // finish it up
        writeln!(out, "{}END CTAB", pre)?;
        writeln!(out, "M  END")?;

        Ok(())
    }

    // Adds an atom with charge and ring closure/ring bond data if applicable to
    // the provided SMILES string. Uses adding_atom1 as a special flag for adding
    // atom1 rather than atom2. It is used when the very first atom is visited.
    fn add_atom_to_smiles(
        &self,
        smiles: &mut String,
        atom1: NodeIndex,
        atom2: Option<NodeIndex>,
        ring_tracker: &mut FxHashMap<NodeIndex, Vec<(u8, BondType)>>,
    ) {
        let node_idx = match atom2 {
            None => atom1,
            Some(atom) => atom,
        };

        let atom = self.graph[node_idx];

        // smiles bracketing rules:
        // - any charge
        // - inorganic (or hydrogen)
        // - un-inferrable hydrogen count
        let organic = OrganicAtom::try_from(atom);
        let charge = self.charges.get(&node_idx);
        let mut bracketed = organic.is_err() || charge.is_some();

        // check for non-standard hydrogen count
        if let Ok(oatom) = organic {
            bracketed = bracketed || oatom == OrganicAtom::H;
            let allocated = self.allocated_valence(node_idx);
            match oatom.remaining_valence(allocated) {
                Some(0) => (),
                None | Some(_) => bracketed = true,
            }
        }

        if bracketed {
            smiles.push('[');
        }

        smiles.push_str(atom.symbol());

        if bracketed {
            // when an atom is bracketed it must have a hydrogen count
            let h_count = self.h_count(&node_idx);
            if h_count != 0 {
                smiles.push('H');
                smiles.push_str(&h_count.to_string());
            }
            if let Some(&num) = charge {
                if num.abs() > 1 {
                    smiles.push_str(&format!("{:+}", num));
                } else {
                    match num.signum() {
                        1 => smiles.push('+'),
                        -1 => smiles.push('-'),
                        0 => (),
                        _ => unreachable!(),
                    }
                }
            }
            smiles.push(']');
        }

        // Adding ring data
        if let Some(ring_bonds) = ring_tracker.remove(&node_idx) {
            for (ring_idx, bond_type) in ring_bonds {
                match bond_type {
                    BondType::Double => smiles.push('='),
                    BondType::Triple => smiles.push('#'),
                    _ => (),
                }
                if ring_idx > 9 {
                    smiles.push('%');
                }
                smiles.push_str(&ring_idx.to_string());
            }
        }
    }

    /// Convert a molecule into SMILES.
    ///
    /// This builds a DFS-based spanning tree of the molecular graph starting at an arbitrary atom.
    /// Rings are tracked by cross or forward edges during the initial DFS traversal. We then
    /// traverse the spanning tree and add atoms to the SMILES string as we visit them.
    ///
    /// Note this follows whatever structure the molecule happens to be in. If hydrogens are
    /// implicit then it will generate a SMILES where hydrogens are implicit. Otherwise it will
    /// add an entry for every hydrogen. There are use cases for this, but it's likely that you
    /// will want to call [`Mol::implicit_hs`] before running this function.
    pub fn to_smiles(&self) -> String {
        // Identify isolated components, to be separated by '.'
        let isolated_components = tarjan_scc(&self.graph);
        let mut smiles = String::with_capacity(self.graph.node_count());
        let mut first_component = true;

        for component in isolated_components {
            let start = component[0];
            let mut component_smiles = String::with_capacity(component.len());
            let mut dfs_tree = FxHashMap::<NodeIndex, Vec<(NodeIndex, BondType)>>::default();
            dfs_tree.reserve(component.len());
            // add the start node in case this component is just that node
            dfs_tree.insert(start, vec![]);
            let mut ring_tracker: FxHashMap<NodeIndex, Vec<(u8, BondType)>> = FxHashMap::default();
            let mut ring_closure_id = 1;

            let mut visited_e = FxHashSet::<EdgeIndex>::default();
            visited_e.reserve(self.graph.edge_count());
            let mut stack = self.graph.edges(start).collect::<Vec<_>>();

            while let Some(e) = stack.pop() {
                let u = e.source();
                let v = e.target();
                let bond = *e.weight();

                if !visited_e.insert(e.id()) {
                    continue;
                }

                if dfs_tree.contains_key(&v) {
                    // check if this is a back edge
                    // if not add to the ring tracker
                    ring_tracker
                        .entry(u)
                        .or_insert_with(Vec::new)
                        .push((ring_closure_id, bond));
                    ring_tracker
                        .entry(v)
                        .or_insert_with(Vec::new)
                        .push((ring_closure_id, bond));
                    ring_closure_id += 1;
                    continue;
                }

                dfs_tree.entry(u).or_insert_with(Vec::new).push((v, bond));
                dfs_tree.entry(v).or_insert_with(Vec::new);

                // add neighbors
                for e in self.graph.edges(v) {
                    stack.push(e);
                }
            }

            // Build mapping of atoms to the number of their neighbors from the DFS tree
            let mut num_neighbors = dfs_tree
                .iter()
                .map(|(idx, neighbors)| (*idx, neighbors.len()))
                .collect::<FxHashMap<NodeIndex, usize>>();
            // Initialize the first atom to have 0 neighbors, in case there are no edges
            // to successfully produce dfs_tree
            num_neighbors.entry(start).or_insert(0);
            let num_neighbors_original = num_neighbors.clone();

            // Add the first atom since SMILES strings must begin with an atom
            self.add_atom_to_smiles(&mut component_smiles, start, None, &mut ring_tracker);
            let mut num_open_paren = 0;

            // This case accounts for multiple branches out of the first atom
            // since otherwise, the branch check is after the SMILES string is
            // first extended by a bond in the DFS below.
            // This access will be safe since the number of neighbors for the start node
            // has been initialized and updated above
            if num_neighbors[&start] > 1 {
                component_smiles.push('(');
                num_open_paren += 1;
            }

            // Traverse DFS tree to construct the SMILES string.
            // Extend the current string with a bond and atom (and ring number,
            // charge, branch parentheses if applicable)
            let mut visited = FxHashSet::<NodeIndex>::default();
            visited.reserve(dfs_tree.len());
            let mut stack = dfs_tree[&start]
                .iter()
                .copied()
                .map(|(v, bond)| (start, v, bond))
                .collect::<Vec<_>>();

            while let Some((u, v, weight)) = stack.pop() {
                if visited.contains(&v) {
                    continue;
                }

                let num_neighbors_u = num_neighbors[&u];
                let num_neighbors_v = num_neighbors[&v];

                // Compare values in the hashmap, if different AND there are
                // remaining branches to explore out of u then we are on a new branch
                if (num_neighbors_u > 1) && (num_neighbors_original[&u] > num_neighbors_u) {
                    component_smiles.push('(');
                    num_open_paren += 1;
                }

                // Add bond to SMILES string
                match weight {
                    BondType::Double => component_smiles.push('='),
                    BondType::Triple => component_smiles.push('#'),
                    _ => (),
                }

                self.add_atom_to_smiles(&mut component_smiles, u, Some(v), &mut ring_tracker);

                // If we reached the end of a branch and have open parenthesis
                // then we need to close the branch
                if (num_neighbors_v == 0) && (num_open_paren > 0) {
                    component_smiles.push(')');
                    num_open_paren -= 1;
                }

                // Begin a new branch leaving v
                if num_neighbors_v >= 2 {
                    component_smiles.push('(');
                    num_open_paren += 1;
                }

                // update number of remaining branches
                *num_neighbors.get_mut(&u).unwrap() -= 1;

                visited.insert(v);

                for (next, bond) in &dfs_tree[&v] {
                    stack.push((v, *next, *bond));
                }
            }

            debug_assert_eq!(num_open_paren, 0);

            if first_component {
                first_component = false;
            } else {
                smiles.push('.');
            }
            smiles.push_str(&component_smiles);
        }
        smiles
    }

    /// Checks if switching two edges would result in an invalid graph. e.g loop, edges with more
    /// than triple bond or not any changes in the graph.
    ///
    /// This is used for edge_switching fn
    fn invalid_bond(&self, edge_idx_1: EdgeIndex, edge_idx_2: EdgeIndex) -> bool {
        if edge_idx_1 == edge_idx_2 {
            return true;
        }
        let (start1, end1) = self.graph.edge_endpoints(edge_idx_1).unwrap();
        let (start2, end2) = self.graph.edge_endpoints(edge_idx_2).unwrap();

        if start1 == start2 || end1 == end2 {
            return true;
        }

        if start1 == end2 || start2 == end1 {
            return true;
        }

        // `.unwrap()` safe because we check for None before using it.
        let index = self.graph.find_edge(start1, end2);
        if index != None && self.graph[index.unwrap()] == BondType::Triple {
            return true;
        }
        let index = self.graph.find_edge(start2, end1);
        if index != None && self.graph[index.unwrap()] == BondType::Triple {
            return true;
        }
        false
    }

    /// Increment the number of bonds on an edge
    fn increment_bond(&mut self, node_idx_1: NodeIndex, node_idx_2: NodeIndex) {
        let index = self.graph.find_edge(node_idx_1, node_idx_2);
        if index == None {
            self.graph
                .add_edge(node_idx_1, node_idx_2, BondType::Single);
        } else {
            match self.graph[index.unwrap()] {
                BondType::Single => {
                    self.graph[index.unwrap()] = BondType::Double;
                }
                BondType::Double => {
                    self.graph[index.unwrap()] = BondType::Triple;
                }
                BondType::Triple => {
                    unreachable!("attempted to increment a triple bond");
                }
                BondType::Aromatic => {
                    unreachable!("attempted to decrement an aromatic bond, molecule must be kekulized before calling Mol::remove_bond");
                }
            }
        };
    }

    /// Decrement the number of bonds on an edge
    /// if the bond is single already, the edge is removed
    fn decrement_bond(&mut self, index: EdgeIndex) {
        match self.graph[index] {
            BondType::Single => {
                self.graph.remove_edge(index);
            }
            BondType::Double => {
                self.graph[index] = BondType::Single;
            }
            BondType::Triple => {
                self.graph[index] = BondType::Double;
            }
            BondType::Aromatic => {
                unreachable!("attempted to decrement an aromatic bond, molecule must be kekulized before calling Mol::remove_bond");
            }
        };
    }

    /// Selects a random edge in the graph
    fn random_edge(&mut self) -> EdgeIndex {
        // `.unwrap()` safe because the edge index is called by `.nth()` on `StableGraph::edge_indices()`
        // we therefore know that it's a valid edge indice
        return self
            .graph
            .edge_indices()
            .nth(rand::thread_rng().gen_range(0..self.graph.edge_count()))
            .unwrap();
    }

    /// Mutates the graph of a molecule by switching the end-points of two randomly selected edges
    /// while keeping the graph connected.
    ///
    /// If the the result is not connected anymore, we try again and keep doing it until it is connected
    pub fn edge_switch(&mut self, distance_limit: u32) {
        for _ in 0..3 {
            let mut i1 = self.random_edge();
            let mut i2 = self.random_edge();

            while self.invalid_bond(i1, i2) {
                i1 = self.random_edge();
                i2 = self.random_edge();
            }

            let (start1, end1) = self.graph.edge_endpoints(i1).unwrap();
            let (start2, end2) = self.graph.edge_endpoints(i2).unwrap();

            let dis1 = dijkstra(&self.graph, start1, Some(end2), |_| 1);

            self.increment_bond(start1, end2);
            self.increment_bond(start2, end1);
            self.decrement_bond(i1);
            self.decrement_bond(i2);

            let dis2 = dijkstra(&self.graph, start1, Some(end1), |_| 1);

            if tarjan_scc(&self.graph).len() == 1
                && match dis1.get(&end2) {
                    Some(&number) => number <= distance_limit,
                    _ => false,
                }
                && match dis2.get(&end1) {
                    Some(&number) => number <= distance_limit,
                    _ => false,
                }
            {
                break;
            } else {
                // if the number of connected components is not 1 anymore we revert what we did and try again
                // `.unwrap()` safe because both we added the same edge in the previous step
                // therefore know that all edge indices are valid
                self.decrement_bond(self.graph.find_edge(start1, end2).unwrap());
                self.decrement_bond(self.graph.find_edge(start2, end1).unwrap());
                self.increment_bond(start1, end1);
                self.increment_bond(start2, end2);
            }
        }
    }

    /// Checks if the provided `bond` is a peptidic bond.
    ///
    /// A peptide bond in this case is defined as a nitrogen to carbon single bond where the carbon
    /// additionally has a double-bonded oxygen and the nitrogen is bound to another carbon. See
    /// the bond marked with a `*` below for the selected peptide bonds.
    /// ```txt
    /// C - N * C = O
    /// ```
    ///
    /// # Panics
    /// If the bond index doesn't exist in the molecule this function panics.
    pub fn is_peptidic(&self, bond: EdgeIndex) -> bool {
        let (source, target) = self.graph.edge_endpoints(bond).expect("invalid bond index");
        let atom1 = self.graph[source];

        // if the endpoints of this edge are N and C pull out relevant indices
        let (carbonyl_idx, n_idx) = match atom1 {
            Atom::C => {
                let atom2 = self.graph[target];
                if atom2 != Atom::N {
                    return false;
                }
                (source, target)
            }
            Atom::N => {
                let atom2 = self.graph[target];
                if atom2 != Atom::C {
                    return false;
                }
                (target, source)
            }
            _ => return false,
        };

        // check if carbonyl carbon really is carbonyl
        self.graph
            .edges(carbonyl_idx)
            .any(|e| e.weight() == &BondType::Double && self.graph[e.target()] == Atom::O)
        &&
        // check for alpha carbon
        self.graph
            .edges(n_idx)
            .filter(|e| e.target() != carbonyl_idx)
            .any(|e| e.weight() == &BondType::Single && self.graph[e.target()] == Atom::C)
    }

    /// Finds all bond indices that correspond to generalized peptide bonds.
    ///
    /// Generalized peptide bonds are taken from the original Dereplicator method. See
    /// Supplementary Figure 6 from the paper [here](https://doi.org/10.1038/nchembio.2219).
    /// The returned iterator may contain duplicate edge indices.
    pub fn generalized_peptide_bonds(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        // These SMILES are set up so the target bond is between nodes 0 and 1
        let subgraphs = vec![
            // peptide bonds SF6d
            MolBuilder::from_smiles("[N](C(=O)[C])[C][C]")
                .unwrap()
                .build(),
            // Dha/Dhb bonds SF6e
            MolBuilder::from_smiles("[N](C(=O)[C])[C]=[C]")
                .unwrap()
                .build(),
            // oxazole/methyl-oxazole bonds SF6f
            MolBuilder::from_smiles("[C]C1O[C]=[C]N=1").unwrap().build(),
            // oxazoline/methyl-oxazoline bonds SF6g
            MolBuilder::from_smiles("[C]C1O[C][C]N=1").unwrap().build(),
            // thiazole bonds SF6h
            MolBuilder::from_smiles("[C]C1=N[C]=[C]S1").unwrap().build(),
        ];

        // can unwrap unless ri_subgraph_isomorphism buggy
        subgraphs
            .into_iter()
            .flat_map(|pattern| self.ri_subgraph_isomorphism(&pattern).into_iter())
            .map(|hit| {
                self.graph
                    .find_edge(hit[&0.into()], hit[&1.into()])
                    .unwrap()
            })
    }

    /// Get all bond indices that correspond to peptide bonds in this molecule.
    ///
    /// See [`Mol::is_peptidic`] for the definition of peptide bonds.
    pub fn peptide_bonds(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.graph
            .edge_indices()
            .filter(|&bond| self.is_peptidic(bond))
    }

    /// Selects a uniformly random sample of connected graph with the same degree sequence.
    /// step is the number of mutations
    pub fn uniform_selection(&mut self, step: u32) {
        for _ in 0..step {
            self.edge_switch(u32::MAX);
        }
    }
}
