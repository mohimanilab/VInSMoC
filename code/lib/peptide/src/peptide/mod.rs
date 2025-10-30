use std::{fmt, iter::FromIterator, ops::Range};

use constants::{
    amino_acids::{CARBONYL_C_IDX, NUM_BACKBONE_ATOMS, N_IDX},
    MASS_WATER, MASS_WATER_EXACT,
};
use itertools::Itertools;
use molecule::{Atom, BondType, Mol};
use petgraph::graph::node_index;
use serde::{Deserialize, Serialize};

use crate::{amino_acid::StdOrNonStdAa, AminoAcid, NonStdAminoAcid, StdAminoAcid};

/// A linear peptide.
///
/// This represents a linear peptide as a sequence of amino acids.
#[derive(Debug, PartialEq, Eq, Hash, Clone, Default, Serialize, Deserialize)]
pub struct Peptide {
    /// The sequence of amino acids that make up this peptide (from N to C terminal)
    sequence: Vec<StdOrNonStdAa>,
}

impl FromIterator<StdAminoAcid> for Peptide {
    fn from_iter<I: IntoIterator<Item = StdAminoAcid>>(iter: I) -> Self {
        Self {
            sequence: iter.into_iter().map(|x| x.into()).collect(),
        }
    }
}

impl FromIterator<NonStdAminoAcid> for Peptide {
    fn from_iter<I: IntoIterator<Item = NonStdAminoAcid>>(iter: I) -> Self {
        Self {
            sequence: iter.into_iter().map(|x| x.into()).collect(),
        }
    }
}

impl FromIterator<StdOrNonStdAa> for Peptide {
    fn from_iter<I: IntoIterator<Item = StdOrNonStdAa>>(iter: I) -> Self {
        Self {
            sequence: iter.into_iter().collect(),
        }
    }
}

impl Peptide {
    /// Creates new peptide
    ///
    /// Constructs a peptide using a vector of Amino acids for the sequence
    pub fn new(seq: Vec<StdOrNonStdAa>) -> Peptide {
        Peptide { sequence: seq }
    }

    /// Total standard mass of peptide.
    ///
    /// This is the sum of the standard residual masses of each constituent amino acid and the mass
    /// of a single water molecule.
    pub fn std_mass(&self) -> f64 {
        self.sequence
            .iter()
            .map(|x| x.residue_std_mass())
            .sum::<f64>()
            + MASS_WATER
    }

    /// Total exact mass of peptide.
    ///
    /// This is the sum of the exact residual masses of each constituent amino acid and the mass
    /// of a single water molecule.
    pub fn exact_mass(&self) -> f64 {
        self.sequence
            .iter()
            .map(|x| x.residue_exact_mass())
            .sum::<f64>()
            + MASS_WATER_EXACT
    }

    /// Conversion into a molecule.
    ///
    /// Connects all amino acid residues with peptide bonds. The first amino acid has a hydrogen
    /// added to its backbone nitrogen. The last amino acid has a hydroxyl group added to its
    /// backbone α-carbon.
    pub fn to_mol(&self) -> Mol {
        // empty list of amino acids
        if self.sequence.is_empty() {
            return Mol::default();
        }

        let mut pep_mol: Mol = Mol::default();
        let mut prev_node_count: usize = 0;
        for i in 0..self.sequence.len() {
            let curr_node_count: usize = pep_mol.graph.node_count();
            let curr_amino_acid = self.sequence[i].to_residue_mol();

            // add nodes from new amino acid to graph
            for node in curr_amino_acid.graph.node_indices() {
                pep_mol.graph.add_node(curr_amino_acid.graph[node]);
            }

            // last amino acid has extra oxygen
            if i == self.sequence.len() - 1 {
                pep_mol.graph.add_node(Atom::O);
            }

            // add edges from new amino acid to graph
            for edge in curr_amino_acid.graph.edge_indices() {
                let (src, end) = curr_amino_acid.graph.edge_endpoints(edge).unwrap();
                let src = src.index() + curr_node_count;
                let end = end.index() + curr_node_count;

                pep_mol.graph.add_edge(
                    node_index(src),
                    node_index(end),
                    curr_amino_acid.graph[edge],
                );
            }

            // add extra edge from C terminus to added O on last amino acid
            if i == self.sequence.len() - 1 {
                let src = curr_node_count + CARBONYL_C_IDX;
                let end = pep_mol.graph.node_count() - 1;
                pep_mol
                    .graph
                    .add_edge(node_index(src), node_index(end), BondType::Single);
            }

            // Connect amino acids
            if i > 0 {
                pep_mol.graph.add_edge(
                    node_index(prev_node_count + 2),
                    node_index(curr_node_count + 1),
                    BondType::Single,
                );
            }

            // create hydrogen hash map
            for node_idx in curr_amino_acid.hydrogens().keys() {
                let new_idx = node_idx.index();

                let h_count = curr_amino_acid.h_count(node_idx);
                if h_count != 0 {
                    *pep_mol.h_entry(node_index(new_idx + curr_node_count)) =
                        curr_amino_acid.h_count(node_idx);
                }
            }

            // Add O's hydrogens on last amino acid
            if i == self.sequence.len() - 1 {
                *pep_mol.h_entry(node_index(pep_mol.graph.node_count() - 1)) = 1;
            }

            prev_node_count = curr_node_count;
        }

        // account for extra hydrogen on nitrogen of first amino acid
        *pep_mol.h_entry(node_index(N_IDX)) += 1;

        pep_mol
    }

    /// Construct a molecule from the cyclic version of this peptide.
    ///
    /// This does a head-to-tail cycle connecting the N-terminal nitrogen with the C-terminal
    /// carbonyl carbon. The returned molecule will have a single macrocycle.
    pub fn to_cyclic(&self) -> Mol {
        // this uses a lot of the info from [`peptide::amino_acid::backbone`]
        // review that for how backbone indices work

        let mut mol = self.to_mol();
        // since r-groups always added after backone, N-terminal nitrogen always has index 1
        let n_terminal_n_idx = node_index(1);
        // get_offset gives us starting point for last amino acid
        // the carbonyl carbon always index 2 of a given AA
        // .unwrap() safe because `self.len() - 1` always in bounds
        let c_terminal_co = node_index(self.get_offset(self.len() - 1).unwrap() + 2);

        // delete OH from C-term COOH
        // this OH is always the very last atom
        mol.remove_node(node_index(mol.graph.node_count() - 1));
        // delete H from N-terminal nitrogen
        *mol.h_entry(n_terminal_n_idx) -= 1;

        // create the new peptide bond
        mol.graph
            .add_edge(n_terminal_n_idx, c_terminal_co, BondType::Single);

        mol
    }

    /// Returns the NodeIndex offset for an amino acid in the Peptide molecule
    ///
    /// This function gives the index of the ith amino acid's alpha carbon. It
    /// iterates through the amino acid chain and adds the number of atoms
    /// of each amino acid until we get to the ith amino acid. This function
    /// can be used to find the offset of specific atoms in any amino acid by
    /// adding the offset from the amino acid's alpha carbon.
    ///
    /// Will return [`None`] if the provided index is out of bounds.
    pub fn get_offset(&self, ind: usize) -> Option<usize> {
        if ind >= self.sequence.len() {
            return None;
        }

        let mut curr_node_count: usize = 0;

        let mut i: usize = 0;
        while i < ind {
            curr_node_count += NUM_BACKBONE_ATOMS;

            let aa = &self.sequence[i];
            curr_node_count += aa.r_group().graph.node_count();

            i += 1;
        }

        Some(curr_node_count)
    }

    /// Get the number of amino acids in this peptide
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if the peptide contains no amino acids
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Returns all NodeIndex offsets for amino acids that matches the input
    ///
    /// This function iterates through the list of amino acids and identifies
    /// the indices that match with the input. It returns the node indices of
    /// these amino acids in order. If the input is `None` it matches with
    /// every amino acid in the range. `start_idx` and `end_idx` can be negative:
    /// they will refer to the index moving from the back of the peptide.
    pub fn find<T: Into<StdOrNonStdAa>>(
        &self,
        input: Option<T>,
        start_idx: i32,
        end_idx: i32,
    ) -> Vec<usize> {
        let mut aa_list: Vec<usize> = Vec::new();

        let start: usize = if start_idx >= 0 {
            start_idx as usize
        } else {
            (start_idx + (self.sequence.len() as i32)) as usize
        };

        let end: usize = if end_idx > 0 {
            end_idx as usize
        } else {
            (end_idx + (self.sequence.len() as i32)) as usize
        };

        if input.is_none() {
            for i in start..end {
                aa_list.push(self.get_offset(i).unwrap());
            }
            return aa_list;
        }

        let aa: StdOrNonStdAa = input.unwrap().into();
        for i in start..end {
            let aa_candidate = self.sequence[i];
            if aa == aa_candidate {
                aa_list.push(self.get_offset(i).unwrap());
            }
        }

        aa_list
    }

    /// Creates an iterator over constituent amino acids.
    pub fn iter(&self) -> impl Iterator<Item = &StdOrNonStdAa> {
        self.sequence.iter()
    }

    /// Construct a new peptide that has amino acids from the provided `indices`.
    ///
    /// # Panics
    /// If the provided range has out of bounds indices.
    pub fn sub_peptide(&self, indices: Range<usize>) -> Self {
        Self::new(self.sequence[indices].to_vec())
    }

    /// Computes all index ranges of the peptide that have length in `length_range`.
    ///
    /// To avoid length filtering run with `0..self.len()+1`.
    pub fn sub_peptide_ranges(
        &self,
        length_range: Range<usize>,
    ) -> impl Iterator<Item = Range<usize>> + '_ {
        let (min_len, max_len) = (length_range.start, length_range.end);
        (0..self.len()).flat_map(move |start| {
            let min_end = start + min_len;
            let max_end = if min_end > self.len() {
                // effectively returning an empty iterator
                0
            } else {
                std::cmp::min(start + max_len, self.len() + 1)
            };
            (min_end..max_end).map(move |end| start..end)
        })
    }

    /// Computes all contiguous substrings of this peptide with lengths in `length_range`.
    ///
    /// This is just a wrapper around [`Peptide::sub_peptide_ranges`].
    pub fn sub_peptides(&self, length_range: Range<usize>) -> impl Iterator<Item = Peptide> + '_ {
        self.sub_peptide_ranges(length_range)
            .map(move |range| self.sub_peptide(range))
    }

    /// Shuffle the order of amino acids in this peptide.
    #[cfg(feature = "rand")]
    pub fn shuffle<R: rand::Rng + ?Sized>(&mut self, rng: &mut R) {
        use rand::seq::SliceRandom;
        self.sequence.shuffle(rng);
    }

    /// Predicts half-open index ranges of cyclic lactone cores that have length in
    /// `length_range`.
    ///
    /// Cores of RiPP class cyclic lactone are subpeptides of the parent ORF whose end
    /// index (exclusive) is the start index of a negatively charged region on the c-terminal
    /// of the ORF.
    ///
    /// This method first predicts the n-terminal start idx of the negatively charged region by
    /// scoring each occurrence of `[StdAminoAcid::D]` or `[StdAminoAcid::A]` based on the amino
    /// acid composition of the region from its index to the ORF's c-terminal. If this region's
    /// amino acids are at least 1/3 negatively charged, its start idx is considered a
    /// potential core end idx. Since core ranges are half-open intervals, this first negatively
    /// charged residue will not be included in the generated cores from this method.
    ///
    /// Core ranges with lengths in `length_range` are then exhaustively generated from each
    /// predicted core end idx.
    ///
    /// ```txt
    /// Peptide:            MENIFNLFIKFFTTILEFIGTVAGDSVCASYFDEPEVPEELTKLYE
    /// Negatively-charged region:                          --_-__--_____-
    ///                                                     ^            ^
    ///                                                 start            end
    /// Generated cores:                           GDSVCASYF
    ///                                           AGDSVCASYF
    ///                                       ...VAGDSVCASYF
    /// ```
    pub fn cyclic_lactone_ranges(
        &self,
        length_range: Range<usize>,
    ) -> impl Iterator<Item = Range<usize>> + '_ {
        let (min_len, max_len) = (length_range.start, length_range.end);
        self.sequence
            .iter()
            .positions(|aa| {
                aa.as_std() == Some(StdAminoAcid::D) || aa.as_std() == Some(StdAminoAcid::E)
            })
            .rev()
            .enumerate()
            .filter(|(num_negative, idx)| {
                2 * (num_negative + 1) >= self.len() - idx - (num_negative + 1)
            })
            .flat_map(move |(_, core_end_idx)| {
                (core_end_idx.saturating_sub(max_len)..core_end_idx.saturating_sub(min_len))
                    .into_iter()
                    .map(move |core_start_idx| core_start_idx..core_end_idx)
            })
    }
}

impl IntoIterator for Peptide {
    type Item = StdOrNonStdAa;

    type IntoIter = std::vec::IntoIter<StdOrNonStdAa>;

    fn into_iter(self) -> Self::IntoIter {
        self.sequence.into_iter()
    }
}

impl fmt::Display for Peptide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.iter().map(|x| x.to_string()).collect::<String>()
        )
    }
}

#[cfg(test)]
mod test;
