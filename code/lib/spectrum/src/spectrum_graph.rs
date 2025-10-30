//! This file defines a spectrum graph.
//!
//! This is used to then compute sequence tags of amino acids in peptides.

use petgraph::{
    graph::{DiGraph, NodeIndex},
    EdgeDirection,
};
use rustc_hash::FxHashSet;

use crate::Peak;

/// A structure to represent a spectrum graph.
///
/// A node in this data structure represents a peak in a mass spectrum, and edges between peaks
/// represent a mass shift whose mass matches some amino acid.
pub struct SpectrumGraph<'a, P> {
    pub(crate) graph: DiGraph<Peak, String>,
    pub(crate) peaks: &'a P,
}

impl<'a, P> SpectrumGraph<'a, P> {
    /// This function returns the peaks which were used to construct the spectrum graph `self`.
    pub fn peaks(&self) -> &'a P {
        self.peaks
    }

    /// This function returns the number of peaks in the spectrum graph `self`.
    pub fn size(&self) -> usize {
        self.graph.node_count()
    }

    /// This function returns all paths of length `length` in the spectrum graph `self`.
    ///
    /// The paths are represented as a vector of strings, where each string is the name of the
    /// amino acid corresponding to the peak in the path. The mass peaks corresponding to these
    /// mass differences are also returned.
    pub fn sequence_tags(&self, length: usize) -> Vec<(Vec<String>, Vec<f64>)> {
        if self.graph.node_count() == 0 {
            return vec![];
        }

        let mut unvisited = self.graph.node_indices().collect::<FxHashSet<_>>();

        let mut all_paths_aa = vec![];
        let mut all_paths_masses = vec![];

        while let Some(&start) = unvisited.iter().next() {
            let mut curr_names = vec![];
            let mut curr_masses = vec![];
            self.dfs(
                start,
                length,
                &mut curr_names,
                &mut all_paths_aa,
                &mut curr_masses,
                &mut all_paths_masses,
                &mut unvisited,
            );
        }

        all_paths_aa
            .into_iter()
            .zip(all_paths_masses.into_iter())
            .collect()
    }

    /// Runs depth first search on the spectrum graph, and appends any new paths with the specified
    /// `length`.
    #[allow(clippy::too_many_arguments)]
    fn dfs(
        &self,
        source: NodeIndex,
        length: usize,
        curr_names: &mut Vec<String>,
        paths_aa: &mut Vec<Vec<String>>,
        curr_masses: &mut Vec<f64>,
        paths_masses: &mut Vec<Vec<f64>>,
        unvisited: &mut FxHashSet<NodeIndex>,
    ) {
        // add the source node to the current path
        curr_masses.push(self.graph[source].mz);
        unvisited.remove(&source);

        // if the current path is longer than the specified length, add the last `length` elements
        // of it to `paths`
        if curr_masses.len() > length {
            // append the last `length` elements to paths
            let start = curr_masses.len() - length - 1;
            paths_aa.push(curr_names[start..].to_vec());
            paths_masses.push(curr_masses[start..].to_vec());
        }

        // deal with the children of this node
        for child in self
            .graph
            .neighbors_directed(source, EdgeDirection::Outgoing)
        {
            // add the name of the monomer with this edge mass
            let edge_index = self.graph.find_edge(source, child).unwrap();
            if let Some(edge_weight) = self.graph.edge_weight(edge_index) {
                curr_names.push(edge_weight.clone());
            }

            // get the edge value from parent to child
            self.dfs(
                child,
                length,
                curr_names,
                paths_aa,
                curr_masses,
                paths_masses,
                unvisited,
            );

            curr_names.pop();
        }

        curr_masses.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Peaks;

    #[derive(Clone)]
    struct MySpectrum {
        peaks: Vec<Peak>,
    }

    impl std::ops::Deref for MySpectrum {
        type Target = Vec<Peak>;
        fn deref(&self) -> &Self::Target {
            &self.peaks
        }
    }

    impl std::ops::DerefMut for MySpectrum {
        fn deref_mut(&mut self) -> &mut Self::Target {
            &mut self.peaks
        }
    }

    impl From<MySpectrum> for Vec<Peak> {
        fn from(spectrum: MySpectrum) -> Self {
            spectrum.peaks
        }
    }

    #[test]
    fn test_dfs() {
        let peak1 = Peak::new(100.0000, 1.0);
        let peak2 = Peak::new(200.0000, 1.0);
        let peak3 = Peak::new(250.0000, 1.0);
        let peak4 = Peak::new(300.0000, 1.0);
        let peak5 = Peak::new(330.0000, 1.0);

        let my_spectrum = MySpectrum {
            peaks: vec![peak1, peak2, peak3, peak4, peak5],
        };
        let masses = vec![
            balloon::of64::try_from(50.0000).unwrap(),
            balloon::of64::try_from(100.0000).unwrap(),
            balloon::of64::try_from(30.0000).unwrap(),
        ];

        let names = vec!["50".to_owned(), "100".to_owned(), "30".to_owned()];

        let spectrum_graph = my_spectrum
            .spectrum_graph(&masses, &names, 331.0, 1e-4)
            .unwrap();

        let paths = spectrum_graph.sequence_tags(4);
        assert_eq!(paths.len(), 1);

        let paths = spectrum_graph.sequence_tags(3);
        assert_eq!(paths.len(), 3);

        let paths = spectrum_graph.sequence_tags(2);
        assert_eq!(paths.len(), 5);
    }
}
