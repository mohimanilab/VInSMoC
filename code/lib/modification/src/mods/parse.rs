use std::{convert::TryFrom, path::Path};

use molecule::{Atom, BondType};
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alpha1, digit1, line_ending, space1},
    combinator::{all_consuming, into, map_res, opt},
    error::{context, convert_error, VerboseError},
    multi::separated_list1,
    sequence::{preceded, separated_pair, terminated},
    IResult,
};

use crate::{
    errors::Error,
    mods::field::{Atoms, Bonds, Field},
};

// start node, end node, bond type for that, newly added atoms, optionally newly added bonds
type AddTuple = (i32, i32, BondType, Atoms, Option<Bonds>);

fn index(input: &str) -> IResult<&str, i32, VerboseError<&str>> {
    map_res(digit1, |idx: &str| idx.parse::<i32>())(input)
}

fn prespace_index(input: &str) -> IResult<&str, i32, VerboseError<&str>> {
    preceded(space1, index)(input)
}

fn bond(input: &str) -> IResult<&str, BondType, VerboseError<&str>> {
    map_res(digit1, |bond: &str| match bond {
        "1" => Ok(BondType::Single),
        "2" => Ok(BondType::Double),
        "3" => Ok(BondType::Triple),
        _ => Err("invalid bond type"),
    })(input)
}

fn disconnect(input: &str) -> IResult<&str, (i32, i32), VerboseError<&str>> {
    let (input, _) = tag("Disconnect")(input)?;
    let (input, idx1) = context("disconnect first index", prespace_index)(input)?;
    let (input, idx2) = context("disconnect second index", prespace_index)(input)?;

    Ok((input, (idx1, idx2)))
}

fn remove(input: &str) -> IResult<&str, i32, VerboseError<&str>> {
    let (input, _) = tag("Remove")(input)?;
    context("remove index", prespace_index)(input)
}

fn connect(input: &str) -> IResult<&str, (i32, i32, BondType), VerboseError<&str>> {
    let (input, _) = tag("Connect")(input)?;

    let (input, idx1) = context("connect first index", prespace_index)(input)?;
    let (input, idx2) = context("connect second index", prespace_index)(input)?;
    let (input, bond) = context("connect bond type", preceded(space1, bond))(input)?;

    Ok((input, (idx1, idx2, bond)))
}

fn atom(input: &str) -> IResult<&str, Atom, VerboseError<&str>> {
    map_res(alpha1, |atom: &str| atom.parse::<Atom>())(input)
}

fn add_atom(input: &str) -> IResult<&str, (Atom, i32), VerboseError<&str>> {
    separated_pair(atom, tag("_"), index)(input)
}

fn add_atoms(input: &str) -> IResult<&str, Atoms, VerboseError<&str>> {
    separated_list1(tag(";"), add_atom)(input)
}

fn add_bond(input: &str) -> IResult<&str, (i32, i32, BondType), VerboseError<&str>> {
    let (input, idx1) = terminated(index, tag(":"))(input)?;
    let (input, (idx2, bond)) = separated_pair(index, tag(":"), bond)(input)?;

    Ok((input, (idx1, idx2, bond)))
}

fn add_bonds(input: &str) -> IResult<&str, Bonds, VerboseError<&str>> {
    separated_list1(tag(";"), add_bond)(input)
}

fn add_last_field(input: &str) -> IResult<&str, (Atoms, Option<Bonds>), VerboseError<&str>> {
    let (input, atoms) = preceded(space1, add_atoms)(input)?;
    let (input, bonds) = opt(preceded(tag(","), add_bonds))(input)?;
    Ok((input, (atoms, bonds)))
}

fn add(input: &str) -> IResult<&str, AddTuple, VerboseError<&str>> {
    let (input, _) = tag("Add")(input)?;

    let (input, idx1) = context("add first index", prespace_index)(input)?;
    let (input, idx2) = context("add second index", prespace_index)(input)?;
    let (input, bond) = context("add bond type", preceded(space1, bond))(input)?;
    let (input, (atoms, bonds)) = context("add final field", add_last_field)(input)?;

    Ok((input, (idx1, idx2, bond, atoms, bonds)))
}

fn charge(input: &str) -> IResult<&str, (i32, i8), VerboseError<&str>> {
    let (input, _) = tag("Charge")(input)?;
    let (input, idx1) = context("charge index", prespace_index)(input)?;
    let (input, charge) = map_res(context("charge charge", prespace_index), i8::try_from)(input)?;
    Ok((input, (idx1, charge)))
}

fn mods(input: &str) -> IResult<&str, Vec<Field>, VerboseError<&str>> {
    all_consuming(separated_list1(
        line_ending,
        alt((
            into(disconnect),
            into(remove),
            into(connect),
            into(add),
            into(charge),
        )),
    ))(input.trim())
}

/// Parses the modification file into a vector of disconnect, removal, connect,
/// add, and charge modifications.
pub fn parse_mod(path: impl AsRef<Path>) -> Result<Vec<Field>, Error> {
    let data = std::fs::read_to_string(path.as_ref())?;
    match mods(&data) {
        Err(nom::Err::Error(e) | nom::Err::Failure(e)) => {
            dbg!(&e);
            Err(Error::Parse(convert_error(data.as_str(), e)))
        }
        Ok((_, fields)) => Ok(fields),
        Err(nom::Err::Incomplete(_)) => unreachable!("we did not use incomplete parsers"),
    }
}
