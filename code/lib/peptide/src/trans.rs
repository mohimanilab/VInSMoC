//! Peptides and amino acids that were translated from DNA.

use std::{fmt, iter::FromIterator, ops::Range, str::FromStr};

use crate::{amino_acid::ParseStdAminoAcidError, Peptide, StdAminoAcid};

/// A translated codon, allowing for stop codon representation.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum TranslatedCodon {
    /// A translated non-stop codon
    Aa(StdAminoAcid),
    /// The stop codon
    Stop,
    /// The codon contained an ambiguous nucleotide
    X,
}

impl TranslatedCodon {
    /// Get the underlying amino acid.
    ///
    /// If this is a translated stop codon or `X` returns [`None`].
    pub fn get_amino_acid(self) -> Option<StdAminoAcid> {
        match self {
            Self::Aa(aa) => Some(aa),
            _ => None,
        }
    }
}

impl From<TranslatedCodon> for char {
    fn from(codon: TranslatedCodon) -> Self {
        match codon {
            TranslatedCodon::Aa(aa) => aa.into(),
            TranslatedCodon::Stop => '*',
            TranslatedCodon::X => 'X',
        }
    }
}

impl FromStr for TranslatedCodon {
    type Err = ParseStdAminoAcidError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "*" => Ok(TranslatedCodon::Stop),
            "X" | "x" => Ok(TranslatedCodon::X),
            x => x.parse::<StdAminoAcid>().map(Self::Aa),
        }
    }
}

impl fmt::Display for TranslatedCodon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", char::from(*self))
    }
}

impl From<StdAminoAcid> for TranslatedCodon {
    fn from(aa: StdAminoAcid) -> Self {
        Self::Aa(aa)
    }
}

/// A [`Peptide`](crate::Peptide) that was translated from a DNA sequence.
///
/// It's possible that such a peptide contains intermediate stop codons and therefore would not in
/// reality exist as a peptide. This is primarily used as an intermediate representation.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct TranslatedPeptide {
    seq: Vec<TranslatedCodon>,
}

impl TranslatedPeptide {
    /// Construct a new [`TranslatedPeptide`].
    pub fn new(seq: Vec<TranslatedCodon>) -> Self {
        Self { seq }
    }

    /// Get a representation of the underlying translated codon sequence.
    pub fn seq(&self) -> &Vec<TranslatedCodon> {
        &self.seq
    }

    /// Get the number of [`TranslatedCodon`]s in this translated peptide.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Check if there are any [`TranslatedCodon`]s in this translated peptide.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Get an iterator over the translated codons.
    pub fn iter(&self) -> std::slice::Iter<TranslatedCodon> {
        self.seq.iter()
    }

    /// Get a mutable iterator over the translated codons.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<TranslatedCodon> {
        self.seq.iter_mut()
    }

    /// Get a single open reading frame corresponding to the given `indices` in this translated
    /// peptide.
    ///
    /// If any of the [`TranslatedCodon`]s at the given `indices` are [`TranslatedCodon::Stop`] or
    /// [`TranslatedCodon::X`] then this returns [`None`]. If used with
    /// [`TranslatedPeptide::orf_ranges`] then this is guaranteed to not return [`None`].
    pub fn get_orf(&self, indices: Range<usize>) -> Option<Peptide> {
        let sub_seq = &self.seq[indices];
        sub_seq
            .iter()
            .all(|aa| aa.get_amino_acid().is_some())
            .then(|| {
                sub_seq
                    .iter()
                    .map(|aa| aa.get_amino_acid().unwrap())
                    .collect()
            })
    }

    /// Computes an iterator over index ranges corresponding to open reading frames.
    ///
    /// An open reading frame (ORF) by our definition is a contiguous region of a
    /// [`TranslatedPeptide`] that begins with a start codon and ends right before either
    /// [`TranslatedCodon::Stop`] or [`TranslatedCodon::X`].
    ///
    /// We also support filtering of ORFs according to their length. The input `length_range`
    /// guarantees that all [`Range`]s in the output iterator fall within `length_range`. If you
    /// simply want all ORFs you can pass in `0..self.len()+1` for this range.
    pub fn orf_ranges(
        &self,
        length_range: Range<usize>,
    ) -> impl Iterator<Item = Range<usize>> + '_ {
        // this is to allow easy moving of the input range
        let (min_len, max_len) = (length_range.start, length_range.end);

        // get all the points where one of our ORFs could end
        // we don't allow ORFs to contain stop or X translated codons
        let end_idxs =
            (0..self.seq.len()).filter(move |idx| self.seq[*idx].get_amino_acid().is_none());

        // iterate over the (curr end index, prev end index) tuples
        // the previous end index for first end index is 0
        let end_idxs = end_idxs.clone().zip(std::iter::once(0).chain(end_idxs));

        end_idxs.flat_map(move |(end, prev_end)| {
            // valid start indices can't be before the last barrier AA (stop or X), since ORFs
            // cannot contain either of these values. Also, they must start at a translated start
            // codon (Aa(M))
            // this also guarantees that start..end is at least min_len
            let start_idxs = (prev_end..end.saturating_sub(min_len) + 1)
                .filter(move |idx| matches!(self.seq[*idx], TranslatedCodon::Aa(StdAminoAcid::M)));

            start_idxs.filter_map(move |start| {
                debug_assert!(end - start >= min_len);

                // we still might find a valid start later
                if end - start >= max_len {
                    return None;
                }

                debug_assert!(self.seq[start..end]
                    .iter()
                    .all(|aa| matches!(aa, TranslatedCodon::Aa(_))));
                Some(start..end)
            })
        })
    }

    /// Computes an iterator over all open reading frames in this [`TranslatedPeptide`].
    ///
    /// This is a wrapper around [`TranslatedPeptide::orf_ranges`], so all documentation from there
    /// applies here.
    pub fn orfs(&self, length_range: Range<usize>) -> impl Iterator<Item = Peptide> + '_ {
        self.orf_ranges(length_range)
            .map(move |range| self.get_orf(range).unwrap())
    }
}

impl<AA: Into<TranslatedCodon>> FromIterator<AA> for TranslatedPeptide {
    fn from_iter<T: IntoIterator<Item = AA>>(iter: T) -> Self {
        Self {
            seq: iter.into_iter().map(|aa| aa.into()).collect(),
        }
    }
}

impl fmt::Display for TranslatedPeptide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for codon in &self.seq {
            write!(f, "{}", codon)?;
        }
        Ok(())
    }
}

impl IntoIterator for TranslatedPeptide {
    type Item = TranslatedCodon;

    type IntoIter = std::vec::IntoIter<TranslatedCodon>;

    fn into_iter(self) -> Self::IntoIter {
        self.seq.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn get_amino_acid() {
        assert_eq!(
            TranslatedCodon::Aa(StdAminoAcid::A).get_amino_acid(),
            Some(StdAminoAcid::A)
        );
        assert_eq!(
            TranslatedCodon::Aa(StdAminoAcid::Q).get_amino_acid(),
            Some(StdAminoAcid::Q)
        );
        assert_eq!(TranslatedCodon::X.get_amino_acid(), None);
        assert_eq!(TranslatedCodon::Stop.get_amino_acid(), None);
    }

    #[test]
    fn to_char() {
        assert_eq!(char::from(TranslatedCodon::Aa(StdAminoAcid::A)), 'A');
        assert_eq!(char::from(TranslatedCodon::Aa(StdAminoAcid::Q)), 'Q');
        assert_eq!(char::from(TranslatedCodon::X), 'X');
        assert_eq!(char::from(TranslatedCodon::Stop), '*');
    }

    #[test]
    fn from_str() {
        assert_eq!(
            "A".parse::<TranslatedCodon>().unwrap(),
            TranslatedCodon::Aa(StdAminoAcid::A)
        );
        assert_eq!(
            "a".parse::<TranslatedCodon>().unwrap(),
            TranslatedCodon::Aa(StdAminoAcid::A)
        );
        assert_eq!(
            "Q".parse::<TranslatedCodon>().unwrap(),
            TranslatedCodon::Aa(StdAminoAcid::Q)
        );
        assert_eq!(
            "q".parse::<TranslatedCodon>().unwrap(),
            TranslatedCodon::Aa(StdAminoAcid::Q)
        );
        assert_eq!("X".parse::<TranslatedCodon>().unwrap(), TranslatedCodon::X);
        assert_eq!("x".parse::<TranslatedCodon>().unwrap(), TranslatedCodon::X);
        assert_eq!(
            "*".parse::<TranslatedCodon>().unwrap(),
            TranslatedCodon::Stop
        );
        assert!("*toolong".parse::<TranslatedCodon>().is_err());
        assert!(".".parse::<TranslatedCodon>().is_err());
        assert!("foo".parse::<TranslatedCodon>().is_err());
    }

    #[test]
    fn display() {
        assert_eq!(
            TranslatedCodon::Aa(StdAminoAcid::A).to_string().as_str(),
            "A"
        );
        assert_eq!(
            TranslatedCodon::Aa(StdAminoAcid::Q).to_string().as_str(),
            "Q"
        );
        assert_eq!(TranslatedCodon::X.to_string().as_str(), "X");
        assert_eq!(TranslatedCodon::Stop.to_string().as_str(), "*");

        assert_eq!(
            TranslatedPeptide::new(vec![
                TranslatedCodon::Aa(StdAminoAcid::A),
                TranslatedCodon::Aa(StdAminoAcid::Q),
                TranslatedCodon::X,
                TranslatedCodon::Stop,
            ])
            .to_string()
            .as_str(),
            "AQX*"
        );
        assert_eq!(
            TranslatedPeptide::new(vec![
                TranslatedCodon::Aa(StdAminoAcid::A),
                TranslatedCodon::Stop,
                TranslatedCodon::X,
                TranslatedCodon::Aa(StdAminoAcid::S),
                TranslatedCodon::Aa(StdAminoAcid::S),
                TranslatedCodon::Aa(StdAminoAcid::S),
                TranslatedCodon::Aa(StdAminoAcid::R),
                TranslatedCodon::Stop,
                TranslatedCodon::X,
                TranslatedCodon::X,
            ])
            .to_string()
            .as_str(),
            "A*XSSSR*XX"
        );
    }

    #[test]
    fn iter() {
        let pep = TranslatedPeptide::new(vec![
            TranslatedCodon::Aa(StdAminoAcid::A),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::R),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::X,
        ]);

        assert_eq!(&pep.clone().into_iter().collect::<Vec<_>>(), pep.seq());
        assert_eq!(
            &pep.clone().into_iter().collect::<TranslatedPeptide>(),
            &pep
        );
        assert_eq!(
            &pep.clone().iter().copied().collect::<TranslatedPeptide>(),
            &pep
        );
        assert_eq!(
            &pep.clone()
                .iter_mut()
                .map(|&mut aa| aa)
                .collect::<TranslatedPeptide>(),
            &pep
        );
    }

    #[test]
    fn orfs() {
        let pep = TranslatedPeptide::new(vec![
            TranslatedCodon::Aa(StdAminoAcid::A),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::R),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::X,
        ]);
        assert_eq!(pep.orf_ranges(0..pep.len()).next(), None);
        assert_eq!(pep.orfs(0..pep.len()).next(), None);

        let pep = TranslatedPeptide::new(vec![
            TranslatedCodon::Aa(StdAminoAcid::A),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::R),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::X,
        ]);

        let mut ranges = pep.orf_ranges(0..pep.len());
        assert_eq!(ranges.next(), Some(3..7));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(0..pep.len());
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::R
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(orfs.next(), None);
        assert_eq!(
            pep.get_orf(3..7),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::R
                ]
                .into_iter()
                .collect()
            )
        );

        // length range should include start length
        let mut ranges = pep.orf_ranges(4..pep.len());
        assert_eq!(ranges.next(), Some(3..7));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(4..pep.len());
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::R
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(orfs.next(), None);

        // but not end length
        let mut ranges = pep.orf_ranges(0..4);
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(0..4);
        assert_eq!(orfs.next(), None);

        let pep = TranslatedPeptide::new(vec![
            TranslatedCodon::Aa(StdAminoAcid::A),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Aa(StdAminoAcid::R),
            TranslatedCodon::Stop,
            TranslatedCodon::X,
            TranslatedCodon::X,
        ]);

        // orfs can overlap if lengths allow it
        let mut ranges = pep.orf_ranges(0..pep.len());
        assert_eq!(ranges.next(), Some(3..7));
        assert_eq!(ranges.next(), Some(5..7));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(0..pep.len());
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::M,
                    StdAminoAcid::R
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(
            orfs.next(),
            Some(vec![StdAminoAcid::M, StdAminoAcid::R].into_iter().collect())
        );
        assert_eq!(
            pep.get_orf(3..7),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::M,
                    StdAminoAcid::R
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(
            pep.get_orf(5..7),
            Some(vec![StdAminoAcid::M, StdAminoAcid::R].into_iter().collect())
        );

        // make sure it works when it starts with M
        let pep = TranslatedPeptide::new(vec![
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Aa(StdAminoAcid::S),
            TranslatedCodon::Stop,
            TranslatedCodon::Aa(StdAminoAcid::R),
        ]);

        let mut ranges = pep.orf_ranges(0..pep.len());
        assert_eq!(ranges.next(), Some(0..4));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(0..pep.len());
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::S
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(orfs.next(), None);
        assert_eq!(
            pep.get_orf(0..4),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::S
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(pep.get_orf(0..5), None);

        let mut ranges = pep.orf_ranges(2..pep.len());
        assert_eq!(ranges.next(), Some(0..4));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(2..pep.len());
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::S
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(orfs.next(), None);

        let mut ranges = pep.orf_ranges(4..pep.len());
        assert_eq!(ranges.next(), Some(0..4));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(4..pep.len());
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::S
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(orfs.next(), None);

        let mut ranges = pep.orf_ranges(4..5);
        assert_eq!(ranges.next(), Some(0..4));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(4..5);
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::S,
                    StdAminoAcid::S,
                    StdAminoAcid::S
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(orfs.next(), None);

        // make sure it works when ends with Stop or X
        let pep = TranslatedPeptide::new(vec![
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Aa(StdAminoAcid::M),
            TranslatedCodon::Stop,
        ]);

        let mut ranges = pep.orf_ranges(0..pep.len());
        assert_eq!(ranges.next(), Some(0..4));
        assert_eq!(ranges.next(), Some(1..4));
        assert_eq!(ranges.next(), Some(2..4));
        assert_eq!(ranges.next(), Some(3..4));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(0..pep.len());
        assert_eq!(
            orfs.next(),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::M,
                    StdAminoAcid::M,
                    StdAminoAcid::M
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(
            orfs.next(),
            Some(
                vec![StdAminoAcid::M, StdAminoAcid::M, StdAminoAcid::M]
                    .into_iter()
                    .collect()
            )
        );
        assert_eq!(
            orfs.next(),
            Some(vec![StdAminoAcid::M, StdAminoAcid::M].into_iter().collect())
        );
        assert_eq!(
            orfs.next(),
            Some(vec![StdAminoAcid::M].into_iter().collect())
        );
        assert_eq!(orfs.next(), None);
        assert_eq!(
            pep.get_orf(0..4),
            Some(
                vec![
                    StdAminoAcid::M,
                    StdAminoAcid::M,
                    StdAminoAcid::M,
                    StdAminoAcid::M
                ]
                .into_iter()
                .collect()
            )
        );
        assert_eq!(
            pep.get_orf(1..4),
            Some(
                vec![StdAminoAcid::M, StdAminoAcid::M, StdAminoAcid::M]
                    .into_iter()
                    .collect()
            )
        );
        assert_eq!(
            pep.get_orf(2..4),
            Some(vec![StdAminoAcid::M, StdAminoAcid::M].into_iter().collect())
        );
        assert_eq!(
            pep.get_orf(3..4),
            Some(vec![StdAminoAcid::M].into_iter().collect())
        );
        assert_eq!(pep.get_orf(4..5), None);

        let mut ranges = pep.orf_ranges(2..4);
        assert_eq!(ranges.next(), Some(1..4));
        assert_eq!(ranges.next(), Some(2..4));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(2..4);
        assert_eq!(
            orfs.next(),
            Some(
                vec![StdAminoAcid::M, StdAminoAcid::M, StdAminoAcid::M]
                    .into_iter()
                    .collect()
            )
        );
        assert_eq!(
            orfs.next(),
            Some(vec![StdAminoAcid::M, StdAminoAcid::M].into_iter().collect())
        );
        assert_eq!(orfs.next(), None);

        let mut ranges = pep.orf_ranges(2..3);
        assert_eq!(ranges.next(), Some(2..4));
        assert_eq!(ranges.next(), None);
        let mut orfs = pep.orfs(2..3);
        assert_eq!(
            orfs.next(),
            Some(vec![StdAminoAcid::M, StdAminoAcid::M].into_iter().collect())
        );
        assert_eq!(orfs.next(), None);
    }
}
