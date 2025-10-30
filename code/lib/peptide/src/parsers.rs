//! Parsers to extract peptides from strings and FASTA files.

use std::iter::FromIterator;
use std::path::Path;

use crate::amino_acid::StdOrNonStdAa;
use crate::{Peptide, StdAminoAcid};

use ::parsers::{
    fasta::{parse_fasta_path, FastaParseError},
    identity_parser,
};

use nom::{
    bytes::complete::take,
    combinator::{all_consuming, map, map_res},
    error::VerboseError,
    multi::many1,
    IResult,
};

fn parse_one_letter_amino_acid(input: &str) -> IResult<&str, StdAminoAcid, VerboseError<&str>> {
    map_res(take(1usize), |s: &str| s.parse::<StdAminoAcid>())(input)
}

fn parse_three_letter_amino_acid(input: &str) -> IResult<&str, StdOrNonStdAa, VerboseError<&str>> {
    map_res(take(3usize), |s: &str| s.parse::<StdOrNonStdAa>())(input)
}

/// Converts an amino acid string into a peptide.
///
/// This amino acid string should consist of single-letter amino acid codes with no whitespace in
/// between consecutive amino acids. It is assumed the amino acid string is given from N to C
/// terminus.
pub fn parse_one_letter_peptide(input: &str) -> IResult<&str, Peptide, VerboseError<&str>> {
    all_consuming(map(
        many1(parse_one_letter_amino_acid),
        |sequence: Vec<StdAminoAcid>| sequence.into_iter().collect::<Peptide>(),
    ))(input)
}

/// Converts an amino acid string into a peptide.
///
/// This amino acid string should consist of three-letter amino acid codes with no whitespace in
/// between consecutive amino acids. It is assumed the amino acid string is given from N to C
/// terminus. Nonstandard amino acids are allowed.
pub fn parse_three_letter_peptide(input: &str) -> IResult<&str, Peptide, VerboseError<&str>> {
    all_consuming(map(
        many1(parse_three_letter_amino_acid),
        |sequence: Vec<StdOrNonStdAa>| sequence.into_iter().collect::<Peptide>(),
    ))(input)
}

/// Parses a peptide sequence from a FASTA file assuming one-letter amino acids.
///
/// The output is a vector of tuples with the format (id, rest of header, parsed [`Peptide`]).
pub fn parse_one_letter_peptide_fasta(
    path: impl AsRef<Path>,
) -> Result<Vec<(String, String, Peptide)>, FastaParseError> {
    let res = parse_fasta_path(path, identity_parser, parse_one_letter_peptide)?;

    let mut out = Vec::<(String, String, Peptide)>::new();
    for (id, comment, pep) in res {
        out.push((id, String::from_iter(comment), pep));
    }

    Ok(out)
}

/// Parses a peptide sequence from a FASTA file assuming three-letter amino acids.
///
/// The output is a vector of tuples with the format (id, rest of header, parsed [`Peptide`]).
pub fn parse_three_letter_peptide_fasta(
    path: impl AsRef<Path>,
) -> Result<Vec<(String, String, Peptide)>, FastaParseError> {
    let res = parse_fasta_path(path, identity_parser, parse_three_letter_peptide)?;

    let mut out = Vec::<(String, String, Peptide)>::new();
    for (id, comment, pep) in res {
        out.push((id, String::from_iter(comment), pep));
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino_acid::StdOrNonStdAa;

    #[test]
    fn parse_peptide_works() {
        assert_eq!(
            parse_one_letter_peptide("F").unwrap().1,
            std::iter::once(StdOrNonStdAa::from(StdAminoAcid::F)).collect::<Peptide>()
        );

        assert_eq!(
            parse_one_letter_peptide("RAVSTE").unwrap().1,
            vec![
                StdOrNonStdAa::from(StdAminoAcid::R),
                StdOrNonStdAa::from(StdAminoAcid::A),
                StdOrNonStdAa::from(StdAminoAcid::V),
                StdOrNonStdAa::from(StdAminoAcid::S),
                StdOrNonStdAa::from(StdAminoAcid::T),
                StdOrNonStdAa::from(StdAminoAcid::E)
            ]
            .into_iter()
            .collect::<Peptide>()
        );

        assert!(parse_one_letter_peptide("DKHX").is_err());
        assert!(parse_one_letter_peptide("R A V S T E").is_err());
    }

    #[test]
    fn parse_fasta_one_letter() {
        let path = "test_files/one_letter_pep.fasta";
        let res = parse_one_letter_peptide_fasta(path);

        // Should be ok and parse one entry
        assert!(res.is_ok());
        let res = res.unwrap();
        assert!(res.len() == 1);

        // Side test: 3-letter parser should fail on this sequence
        let res2 = parse_three_letter_peptide_fasta(path);
        assert!(res2.is_err());

        let expected_pep = vec![
            StdOrNonStdAa::from(StdAminoAcid::R),
            StdOrNonStdAa::from(StdAminoAcid::A),
            StdOrNonStdAa::from(StdAminoAcid::V),
            StdOrNonStdAa::from(StdAminoAcid::S),
            StdOrNonStdAa::from(StdAminoAcid::T),
            StdOrNonStdAa::from(StdAminoAcid::E),
        ]
        .into_iter()
        .collect::<Peptide>();

        let (id, comment, observed_pep) = res[0].clone();

        assert_eq!(observed_pep, expected_pep);
        assert_eq!(id, "id");
        assert_eq!(comment, "comment");
    }

    #[test]
    fn parse_fasta_three_letter() {
        let path = "test_files/three_letter_pep.fasta";
        let res = parse_three_letter_peptide_fasta(path);

        // Should be ok and parse one entry
        assert!(res.is_ok());
        let res = res.unwrap();
        assert!(res.len() == 1);

        // Side test: 1-letter parser should fail on this sequence
        let res2 = parse_one_letter_peptide_fasta(path);
        assert!(res2.is_err());

        let expected_pep = vec![
            StdOrNonStdAa::from(StdAminoAcid::R),
            StdOrNonStdAa::from(StdAminoAcid::A),
            StdOrNonStdAa::from(StdAminoAcid::V),
            StdOrNonStdAa::from(StdAminoAcid::S),
            StdOrNonStdAa::from(StdAminoAcid::T),
            StdOrNonStdAa::from(StdAminoAcid::E),
        ]
        .into_iter()
        .collect::<Peptide>();

        let (id, comment, observed_pep) = res[0].clone();

        assert_eq!(observed_pep, expected_pep);
        assert_eq!(id, "id");
        assert_eq!(comment, "comment");
    }

    #[test]
    fn invalid_whitespace() {
        let path = "test_files/bad1.fasta";
        let res = parse_three_letter_peptide_fasta(path);
        let res2 = parse_one_letter_peptide_fasta(path);

        // Should fail due to the spacing in the sequence
        assert!(res.is_err());
        assert!(res2.is_err());
    }

    #[test]
    fn invalid_aminos() {
        let path = "test_files/bad2.fasta";
        let res = parse_three_letter_peptide_fasta(path);
        let res2 = parse_one_letter_peptide_fasta(path);

        // Should fail due to nonexistant amino acids in the sequence
        assert!(res.is_err());
        assert!(res2.is_err());
    }

    #[test]
    fn invalid_fasta() {
        let path = "test_files/bad3.fasta";
        let res = parse_three_letter_peptide_fasta(path);
        let res2 = parse_one_letter_peptide_fasta(path);

        // Should fail due to invalid fasta format
        assert!(res.is_err());
        assert!(res2.is_err());
    }
}
