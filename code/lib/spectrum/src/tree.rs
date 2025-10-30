use std::{collections::BTreeMap, path::PathBuf};

use balloon::of64;
use molecule::Adduct;
use petgraph::{graphmap::DiGraphMap, Direction::Incoming};

use crate::{Ms, Spectrum};

/// [`MultistageTree`] is a directed graph map whose nodes are indices in a vector of
/// multistage ingested spectra and edges exist from precursor spectra to their product spectra.
///
/// Only spectra with precursor-product relationships will have nodes in the tree.
#[derive(Default, Debug)]
pub struct MultistageTree {
    graph: DiGraphMap<usize, ()>,
}

impl MultistageTree {
    /// Constructs a [`MultistageTree`] from the `spectral_data` and `metadata` produced by
    /// a struct upon spectral ingestion.
    pub fn from_ingested_spectra(
        spectral_data: &[Spectrum],
        metadata: &[(PathBuf, usize)],
    ) -> Self {
        // if metadata is empty, the data was likely not parsed from ParsedCollections, meaning
        // they are likely from the same spectrum collection
        if metadata.is_empty() {
            return Self::from_spectra_in_collection(spectral_data);
        }
        // Each tuple in `metadata` corresponds to an ingested spectrum collection, whose last
        // spectrum was ingested into `spectral_data` one index below that specified by the
        // tuple's second field.
        let mut coll_start = 0;
        let mut tree = Self::default();
        for (_, coll_end) in metadata {
            let spectrum_collection = &spectral_data[coll_start..*coll_end];
            // append the tree for this spectrum collection to the current tree, offsetting
            // added indices by the collection's start index
            tree.append(
                &Self::from_spectra_in_collection(spectrum_collection),
                coll_start,
            );
            coll_start = *coll_end;
        }
        tree
    }

    /// Constructs a [`MultistageTree`] from a slice of spectra corresponding to a
    /// complete spectrum collection.
    pub fn from_spectra_in_collection(spectra: &[Spectrum]) -> Self {
        Self {
            graph: spectra
                .iter()
                .enumerate()
                .filter_map(|(product_idx, product_spectrum)| {
                    // zipping here only retains (precursor, product) index pairs for product
                    // indices with:
                    // 1. a defined precursor scan
                    // 2. ...that is found in the current collection
                    product_spectrum
                        .precursor_scan()
                        .and_then(|precursor_scan| {
                            spectra
                                .binary_search_by(|spectrum| spectrum.scan.cmp(&precursor_scan))
                                .ok()
                        })
                        .zip(Some(product_idx))
                })
                .collect::<DiGraphMap<usize, ()>>(),
        }
    }

    /// Append a [`MultistageTree`] onto the current [`MultistageTree`], offsetting the appended
    /// tree's indices by the provided `offset`.
    ///
    /// The provided `offset` should ensure:
    /// 1. the correspondence of spectral indices with spectral data ingested from more than
    ///    one collection.
    /// 2. that each appended tree is disconnected from any previously added trees.
    fn append(&mut self, other: &MultistageTree, offset: usize) {
        self.graph.extend(
            other
                .graph
                .all_edges()
                .map(|(precursor, product, ())| (precursor + offset, product + offset, ())),
        );
    }

    /// Returns an iterator over all directed precursor-product edges found in the
    /// `[MultistageTree]`.
    pub fn all_edges(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.graph
            .all_edges()
            .map(|(precursor_idx, product_idx, _)| (precursor_idx, product_idx))
    }

    /// Returns an iterator over indices of each component's root spectrum for this [`MultistageTree`].
    ///
    /// While not enforced, these are expected to have an MS level of 1.
    pub fn component_roots(&self) -> impl Iterator<Item = usize> + '_ {
        self.graph
            .nodes()
            .filter(|&node| self.graph.edges_directed(node, Incoming).next().is_none())
    }

    /// Returns an iterator over indices of all spectra which are direct children of the
    /// component roots of this [`MultistageTree`]. In multistage spectra, this corresponds to
    /// spectra whose precursor mass is that of an unfragmented precursor molecule, hence the
    /// usage of the term `precursor`.
    ///
    /// While not enforced, these are expected to have an MS level of 2.
    pub fn precursor_roots(&self) -> impl Iterator<Item = usize> + '_ {
        self.component_roots()
            .flat_map(|root_idx| self.graph.neighbors(root_idx))
    }

    /// Returns an iterator over indices of all product spectra for the precursor spectrum
    /// specified by `precursor_idx`.
    pub fn product_idxs(&self, precursor_idx: usize) -> impl Iterator<Item = usize> + '_ {
        self.graph.neighbors(precursor_idx)
    }

    /// Computes a map from precursor m/zs to scans. See [`MultistageTree::precursor_roots`]
    /// for more information about this usage of the term `precursor`.
    pub fn precursor_map<'a>(
        &self,
        adducts: &'a [Adduct],
        spectra: &[Spectrum],
    ) -> BTreeMap<of64, Vec<(&'a Adduct, usize)>> {
        let mut precursor_map = BTreeMap::default();
        for spectrum_idx in self.precursor_roots() {
            for adduct in adducts {
                let mass = adduct.to_mass(spectra[spectrum_idx].pepmass());
                precursor_map
                    .entry(of64::try_from(mass).unwrap())
                    .or_insert_with(Vec::new)
                    .push((adduct, spectrum_idx));
            }
        }
        precursor_map
    }

    /// Return an iterator over indices of all members of the subgraph rooted at the given
    /// `subtree_root_idx`, including the subtree root.
    pub fn descendants(&self, subtree_root_idx: usize) -> impl Iterator<Item = usize> + '_ {
        SubtreeMembers::new(self, subtree_root_idx)
    }

    /// Return an iterator over all paths in the subgraph rooted at the given `subtree_root_idx`,
    /// including the subtree root. Paths are obtained via depth-first-search starting at the
    /// subtree root.
    pub fn paths(&self, subtree_root_idx: usize) -> impl Iterator<Item = Vec<usize>> + '_ {
        SubtreePaths::new(self, subtree_root_idx)
    }

    /// Returns the number of spectra in the subtree rooted at the given `subtree_root_idx`,
    /// including the subtree root.
    pub fn subtree_size(&self, subtree_root_idx: usize) -> usize {
        self.descendants(subtree_root_idx).count()
    }
}

/// A data structure for efficient computation of members of [`MultistageTree`] subtrees.
struct SubtreeMembers<'tree> {
    tree: &'tree MultistageTree,
    q: Vec<usize>,
}

impl SubtreeMembers<'_> {
    /// Construct a [`SubtreeMembers`] by initializing a queue containing the provided
    /// `subtree_root_idx`.
    fn new(tree: &MultistageTree, subtree_root_idx: usize) -> SubtreeMembers<'_> {
        SubtreeMembers {
            tree,
            q: vec![subtree_root_idx],
        }
    }
}

impl Iterator for SubtreeMembers<'_> {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        self.q.pop().map(|precursor_idx| {
            self.q.extend(self.tree.product_idxs(precursor_idx));
            precursor_idx
        })
    }
}

/// A data structure for computation of paths from [`MultistageTree`] subtrees.
struct SubtreePaths<'tree> {
    tree: &'tree MultistageTree,
    q: Vec<(Vec<usize>, usize)>,
}

impl SubtreePaths<'_> {
    /// Construct a [`SubtreePaths`] by initializing a queue containing the provided
    /// `subtree_root_idx`.
    fn new(tree: &MultistageTree, subtree_root_idx: usize) -> SubtreePaths<'_> {
        SubtreePaths {
            tree,
            q: vec![(Vec::new(), subtree_root_idx)],
        }
    }
}

impl Iterator for SubtreePaths<'_> {
    type Item = Vec<usize>;
    fn next(&mut self) -> Option<Self::Item> {
        self.q.pop().map(|(mut path_before, precursor_idx)| {
            path_before.push(precursor_idx);
            self.q.extend(
                self.tree
                    .product_idxs(precursor_idx)
                    .map(|product_idx| (path_before.clone(), product_idx)),
            );
            path_before
        })
    }
}
