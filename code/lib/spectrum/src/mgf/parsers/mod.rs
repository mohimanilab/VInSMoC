// The parser provided here is largely based on the C++ parser from libmgf
// here: https://github.com/kirchnerlab/libmgf/blob/master/src/Parser.ypp
use super::{GlobalMgfKey, Mgf, Spectrum};
use crate::Peak;
use std::convert::TryFrom;

use nom::{
    bytes::complete::tag,
    character::complete::{char, line_ending, multispace0, not_line_ending, space0, space1},
    combinator::{all_consuming, map},
    error::{Error, ErrorKind},
    multi::{many0, many1},
    number::complete::double,
    sequence::{pair, preceded, separated_pair, terminated},
    Err as NomErr, IResult,
};

#[cfg(test)]
mod tests;

mod params;
use params::{global_params, local_params};

/// An alternative to [`line_ending`] that allows any number of preceding whitespace.
fn line_end(input: &str) -> IResult<&str, ()> {
    map(pair(space0, line_ending), |_| ())(input)
}

fn comment(input: &str) -> IResult<&str, &str> {
    preceded(
        multispace0,
        preceded(char('#'), terminated(not_line_ending, line_ending)),
    )(input)
}

// Parse peaks that represent ions:
// Peaks are represented as two floating point numbers on the same line,
// seperated by a single space.
fn ion(input: &str) -> IResult<&str, Peak> {
    terminated(
        map(
            separated_pair(double, space1, double),
            |(m_z, intensity)| Peak::new(m_z, intensity),
        ),
        line_end,
    )(input)
}

// Repeats the above function to get a vector of Peaks for a given spectra.
fn ions(input: &str) -> IResult<&str, Vec<Peak>> {
    many0(ion)(input)
}

fn begin_ions(input: &str) -> IResult<&str, &str> {
    terminated(tag("BEGIN IONS"), line_end)(input)
}

fn end_ions(input: &str) -> IResult<&str, &str> {
    terminated(tag("END IONS"), line_end)(input)
}

// One block consists of a Spectrum. Each spectrum in an MGF file is represented as
// a list of local variables followed by a set of peaks. These two are wrapped between
// BEGIN IONS...END IONS. Combining the data stored inside, we can get the result for
// an individual Spectrum.
fn block(input: &str) -> IResult<&str, Spectrum> {
    let (input, _) = begin_ions(input)?;
    let (input, local_vec) = local_params(input)?;
    let (input, peak_vec) = ions(input)?;
    let (input, _) = end_ions(input)?;

    match Spectrum::try_from((local_vec, peak_vec)) {
        Ok(spec) => Ok((input, spec)),
        Err(_) => Err(NomErr::Error(Error::new(input, ErrorKind::Tag))),
    }
}

// This parser collects all the Spectrum queries provided in the MGF file.
fn blocks(input: &str) -> IResult<&str, Vec<Spectrum>> {
    many1(terminated(block, many0(line_end)))(input)
}

fn mgf_file(input: &str) -> IResult<&str, (Vec<GlobalMgfKey>, Vec<Spectrum>)> {
    let (input, attrs) = global_params(input)?;
    let (input, _) = many0(line_end)(input)?;
    let (input, spectra) = all_consuming(blocks)(input)?;
    Ok((input, (attrs, spectra)))
}

/// This parser puts together the global parameters and the set of queries
/// that consist of the MGF file.
///
/// An MGF file consists of the global parameters at the beginning of the file, followed by data from queries that provide
/// information for different Spectra, with each query corresponding to one Spectrum. This function extracts all the data
/// for the MGF file and produces it in an MgfFile type.
pub fn parse_mgf(input: &str) -> Result<Mgf, &str> {
    match mgf_file(input) {
        Ok((_, (attrs, spectra))) => Ok(Mgf { spectra, attrs }),
        Err(x) => {
            eprintln!("{:?}", x);
            Err("Failed to parse MGF")
        }
    }
}
