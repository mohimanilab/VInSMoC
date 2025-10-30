use super::{comment, line_end};
use crate::{
    mgf::{GlobalMgfKey, LocalMgfKey},
    ChargeVec,
};
use nom::{
    branch::alt,
    bytes::complete::{tag, take_while_m_n},
    character::complete::{alpha1, alphanumeric1, digit1, not_line_ending, space0, space1},
    combinator::{cut, map, map_res, opt, recognize},
    multi::{many0, separated_list1},
    number::complete::double,
    sequence::{delimited, pair, preceded, terminated, tuple},
    IResult,
};

#[cfg(test)]
mod tests;

// The functions in this module parse a distinct parameter, which is in the form
// <NAME>=<Value>. A comprehensive list of names can be found at
// https://www.matrixscience.com/help/data_file_help.html .

// Creates parser functions for simple string/numeric local/global keys
macro_rules! key_parser {
    (string local => $func_name:ident, $tag:literal, $variant:expr) => {
        key_parser!(string crate::mgf::LocalMgfKey => $func_name, $tag, $variant);
    };
    (string global => $func_name:ident, $tag:literal, $variant:expr) => {
        key_parser!(string crate::mgf::GlobalMgfKey => $func_name, $tag, $variant);
    };
    (num local $ft:ty => $func_name:ident, $tag:literal, $variant:expr) => {
        key_parser!(num crate::mgf::LocalMgfKey,$ft => $func_name, $tag, $variant);
    };
    (num global $ft:ty => $func_name:ident, $tag:literal, $variant:expr) => {
        key_parser!(num crate::mgf::GlobalMgfKey,$ft => $func_name, $tag, $variant);
    };
    (string $keytype:path => $func_name:ident, $tag:literal, $variant:expr) => {
        fn $func_name(input: &str) -> IResult<&str, $keytype> {
            use $keytype::*;
            preceded(
                tag($tag),
                terminated(map(not_line_ending, |s: &str| $variant(s.to_string())), line_end),
            )(input)
        }
    };
    (num $keytype:path,$ft:ty => $func_name:ident, $tag:literal, $variant:expr) => {
        fn $func_name(input: &str) -> IResult<&str, $keytype> {
            use $keytype::*;
            preceded(
                tag($tag),
                cut(terminated(map_res(not_line_ending, |s: &str| {
                    match s.parse::<$ft>() {
                        Ok(x) => Ok($variant(x)),
                        Err(x) => Err(x)
                    }
                }), line_end)),
            )(input)
        }
    }
}

// Local keys
key_parser!(string local => title_local, "TITLE=", Title);
key_parser!(string local => comp_local, "COMP=", Comp);
key_parser!(string local => instrument_local, "INSTRUMENT=", Instrument);
key_parser!(string local => it_mods_local, "IT_MODS=", ItMods);
key_parser!(string local => tolu_local, "TOLU=", Tolu);
key_parser!(string local => seq_local, "SEQ=", Seq);
key_parser!(string local => tag_local, "TAG=", Tag);
key_parser!(string local => etag_local, "ETAG=", Etag);
key_parser!(num local f64 => tol_local, "TOL=", Tol);

// Global keys
key_parser!(string global => cle_global, "CLE=", Cle);
key_parser!(string global => com_global, "COM=", Com);
key_parser!(string global => db_global, "DB=", Db);
key_parser!(string global => format_global, "FORMAT=", Format);
key_parser!(string global => instrument_global, "INSTRUMENT=", Instrument);
key_parser!(string global => it_mods_global, "IT_MODS=", ItMods);
key_parser!(string global => itolu_global, "ITOLU=", Itolu);
key_parser!(string global => mass_global, "MASS=", Mass);
key_parser!(string global => mods_global, "MODS=", Mods);
key_parser!(string global => quantitation_global, "QUANTITATION=", Quantitation);
key_parser!(string global => report_global, "REPORT=", Report);
key_parser!(string global => rep_type_global, "REPTYPE=", RepType);
key_parser!(string global => search_global, "SEARCH=", Search);
key_parser!(string global => taxonomy_global, "TAXONOMY=", Taxonomy);
key_parser!(string global => tolu_global, "TOLU=", Tolu);
key_parser!(string global => user_global, "USER=", User);
key_parser!(string global => user_email_global, "USEREMAIL=", UserEmail);
key_parser!(string global => username_global, "USERNAME=", UserName);
key_parser!(num global i32 => pfa_global, "PFA=", Pfa); // TODO: make u8
key_parser!(num global f64 => itol_global, "ITOL=", Itol);
key_parser!(num global f64 => pep_isotope_error_global, "PEP_ISOTOPE_ERROR=", PepIsotopeError);
key_parser!(num global f64 => precursor_global, "PRECURSOR=", Precursor);
key_parser!(num global f64 => seg_global, "SEG=", Seg);
key_parser!(num global f64 => tol_global, "TOL=", Tol);
key_parser!(string global => frames_global, "FRAMES=", Frames);

// Special case for PEPMASS (due to possible intensity)
fn peptide_mass_local(input: &str) -> IResult<&str, LocalMgfKey> {
    map(
        delimited(
            tag("PEPMASS="),
            terminated(double, opt(preceded(space1, double))),
            line_end,
        ),
        LocalMgfKey::PepMass,
    )(input)
}

// Special case for CHARGE= keys

fn charge_sign(input: &str) -> IResult<&str, bool> {
    map(take_while_m_n(0, 1, |c| "+-".contains(c)), |c: &str| {
        // if no sign specified assume positive mode
        c.is_empty() || c == "+"
    })(input)
}

fn charge(input: &str) -> IResult<&str, ChargeVec> {
    let (input, charge_vec) = preceded(
        tag("CHARGE="),
        cut(separated_list1(
            // allow multiple charges separated by , or and
            // whitespace around delimiters are discarded
            delimited(space0, alt((tag(","), tag("and"))), space0),
            alt((
                map_res(
                    tuple((digit1, charge_sign)),
                    |(charge, sign): (&str, bool)| {
                        let charge: i8 = match charge.parse::<i8>() {
                            Ok(x) => x,
                            Err(x) => return Err(x),
                        };
                        Ok(match sign {
                            true => charge,
                            false => -charge,
                        })
                    },
                ),
                map_res(
                    tuple((charge_sign, digit1)),
                    |(sign, charge): (bool, &str)| {
                        let charge: i8 = match charge.parse::<i8>() {
                            Ok(x) => x,
                            Err(x) => return Err(x),
                        };
                        Ok(match sign {
                            true => charge,
                            false => -charge,
                        })
                    },
                ),
            )),
        )),
    )(input)?;
    let charge_vec = ChargeVec::from_vec(charge_vec);
    Ok((input, charge_vec))
}

fn charge_local(input: &str) -> IResult<&str, LocalMgfKey> {
    terminated(map(charge, LocalMgfKey::Charge), line_end)(input)
}

fn charge_global(input: &str) -> IResult<&str, GlobalMgfKey> {
    terminated(map(charge, GlobalMgfKey::Charge), line_end)(input)
}

// Special cases for range-based keys (RTINSECONDS and SCANS)
fn rt_in_seconds_local(input: &str) -> IResult<&str, LocalMgfKey> {
    preceded(
        tag("RTINSECONDS="),
        terminated(
            map(
                tuple((double, opt(preceded(tag("-"), double)))),
                |(start, end): (f64, Option<f64>)| {
                    let end = end.unwrap_or(start + f64::EPSILON);
                    LocalMgfKey::RtInSeconds(start..end)
                },
            ),
            line_end,
        ),
    )(input)
}

fn scans_local(input: &str) -> IResult<&str, LocalMgfKey> {
    preceded(
        tag("SCANS="),
        terminated(
            map_res(
                tuple((digit1, opt(preceded(tag("-"), digit1)))),
                |(start, end): (&str, Option<&str>)| {
                    let start = match start.parse::<u32>() {
                        Ok(x) => x,
                        Err(x) => return Err(x),
                    };

                    let end = match end {
                        None => start,
                        Some(x) => match x.parse::<u32>() {
                            Ok(y) => y,
                            Err(y) => return Err(y),
                        },
                    } + 1;
                    Ok(LocalMgfKey::Scans(start..end))
                },
            ),
            line_end,
        ),
    )(input)
}

// Special cases for boolean keys (DECOY and ERRORTOLERANT)
fn decoy_global(input: &str) -> IResult<&str, GlobalMgfKey> {
    preceded(
        tag("DECOY="),
        terminated(
            map_res(not_line_ending, |s: &str| {
                match s.parse::<u8>() {
                    Ok(x) => match x {
                        0 => Ok(GlobalMgfKey::Decoy(false)),
                        _ => Ok(GlobalMgfKey::Decoy(true)), // allow anything non-zero to be true
                    },
                    Err(x) => Err(x),
                }
            }),
            line_end,
        ),
    )(input)
}

fn error_tolerant_global(input: &str) -> IResult<&str, GlobalMgfKey> {
    preceded(
        tag("ERRORTOLERANT="),
        terminated(
            map_res(not_line_ending, |s: &str| {
                match s.parse::<u8>() {
                    Ok(x) => match x {
                        0 => Ok(GlobalMgfKey::ErrorTolerant(false)),
                        _ => Ok(GlobalMgfKey::ErrorTolerant(true)), // allow anything non-zero to be true
                    },
                    Err(x) => Err(x),
                }
            }),
            line_end,
        ),
    )(input)
}

// take from the nom recipes module
fn identifier(input: &str) -> IResult<&str, &str> {
    recognize(pair(
        alt((alpha1, tag("_"))),
        many0(alt((alphanumeric1, tag("_")))),
    ))(input)
}

// Cleans up any weird unsupported attributes in an otherwise correct mgf file
fn unsupported_attribute(input: &str) -> IResult<&str, ()> {
    let (input, _) = tuple((identifier, tag("="), not_line_ending, line_end))(input)?;

    Ok((input, ()))
}

fn local_param(input: &str) -> IResult<&str, LocalMgfKey> {
    alt((
        peptide_mass_local,
        charge_local,
        title_local,
        comp_local,
        instrument_local,
        it_mods_local,
        rt_in_seconds_local,
        scans_local,
        tolu_local,
        seq_local,
        tag_local,
        etag_local,
        tol_local,
    ))(input)
}

fn global_param(input: &str) -> IResult<&str, GlobalMgfKey> {
    alt((
        cle_global,
        com_global,
        db_global,
        format_global,
        instrument_global,
        it_mods_global,
        itolu_global,
        mass_global,
        mods_global,
        mods_global,
        quantitation_global,
        report_global,
        rep_type_global,
        search_global,
        taxonomy_global,
        tolu_global,
        user_global,
        user_email_global,
        username_global,
        decoy_global,
        alt((
            error_tolerant_global,
            pfa_global,
            itol_global,
            pep_isotope_error_global,
            precursor_global,
            seg_global,
            tol_global,
            charge_global,
            frames_global,
        )),
    ))(input)
}

enum KeyOrComment {
    Unsupported,
    Comment,
    Local(LocalMgfKey),
    Global(GlobalMgfKey),
}

// This parser finds the full list of local parameters in order, until
// we reach the peaks of the spectra in the MGF file.
pub fn local_params(input: &str) -> IResult<&str, Vec<LocalMgfKey>> {
    // many0_comments(local_param)(input)
    let (input, key_or_comment) = many0(alt((
        map(local_param, KeyOrComment::Local),
        map(comment, |_| KeyOrComment::Comment),
        map(unsupported_attribute, |_| KeyOrComment::Unsupported),
    )))(input)?;

    Ok((
        input,
        key_or_comment
            .into_iter()
            .filter_map(|x| match x {
                KeyOrComment::Local(key) => Some(key),
                KeyOrComment::Comment | KeyOrComment::Unsupported => None,
                _ => unreachable!(),
            })
            .collect(),
    ))
}

// This parser finds the full list of global parameters in order, until
// we reach the queries for spectra in the MGF file.
pub fn global_params(input: &str) -> IResult<&str, Vec<GlobalMgfKey>> {
    let (input, key_or_comment) = many0(alt((
        map(global_param, KeyOrComment::Global),
        map(comment, |_| KeyOrComment::Comment),
        map(unsupported_attribute, |_| KeyOrComment::Unsupported),
    )))(input)?;

    Ok((
        input,
        key_or_comment
            .into_iter()
            .filter_map(|x| match x {
                KeyOrComment::Global(key) => Some(key),
                KeyOrComment::Comment | KeyOrComment::Unsupported => None,
                _ => unreachable!(),
            })
            .collect(),
    ))
}
