use crate::Atom;

use nom::{
    branch::alt,
    bytes::complete::take_while_m_n,
    character::complete::digit1,
    combinator::{map_res, peek},
    error::VerboseError,
    multi::many1_count,
    sequence::preceded,
    IResult,
};
use std::convert::TryFrom;

mod smiles;
pub use self::smiles::{parse_indexed_smiles, parse_smiles, SmilesParseError};

mod formula_adduct;
pub use formula_adduct::{parse_adduct, parse_formula};

mod mol_file;
pub use mol_file::{parse_mol_file, MolFileParseError};

fn charge_sign(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    take_while_m_n(1, 1, |c| "+-".contains(c))(input)
}

/// Parses a non-neutral charge from a string.
/// First, if the charge string is an empty string we return an error.
/// Next, if the charge is not empty it must begin with its sign character, either + or -.
/// After the sign character there can be a number, denoting the absolute value of the
/// charge. Alternatively, there can be multiple of the same sign character, in which case the
/// absolute charge is the number of times that sign character appears. For example, ++ and +2 both
/// denote the absolute charge 2 with a positive sign.
fn forced_charge(input: &str) -> IResult<&str, i8, VerboseError<&str>> {
    // first find out the charge sign character
    // if there is no such charcter, our charge must be empty and therefore we should return an error
    let (input, charge_char) = peek(charge_sign)(input)?;

    // match exactly one of our determined sign character
    let sign_matcher = take_while_m_n(1, 1, |c| charge_char.starts_with(c));

    // Try, in order
    //  1. To parse an absolute charge number
    //  2. To count the number of the charge character
    // and parse the resulting absolute charge into a usize.
    // Afterwards convert the usize into an i8, propagating error up on failure.
    // This gives us our absolute charge in a signed integer.
    let (input, absolute_charge) = map_res(
        alt((
            map_res(preceded(&sign_matcher, digit1), |s: &str| {
                s.parse::<usize>()
            }),
            many1_count(&sign_matcher),
        )),
        i8::try_from,
    )(input)?;

    // Add the sign of the charge to the absolute charge
    match charge_char {
        "+" => Ok((input, absolute_charge)),
        "-" => Ok((input, -absolute_charge)),
        _ => unreachable!(),
    }
}

fn remaining_hydrogens(atom: &Atom, num_bonds: u8, charge: Option<i8>) -> Result<u8, &'static str> {
    let num_bonds_charge = match charge {
        Some(x) => (num_bonds as i8) - x,
        None => (num_bonds as i8),
    };

    let allowed_valence: i8 = match atom {
        Atom::B => 3,
        Atom::C => 4,
        Atom::N => {
            if num_bonds_charge <= 3 {
                3
            } else if num_bonds_charge <= 5 {
                5
            } else {
                0
            }
        }
        Atom::O => 2,
        Atom::P => {
            if num_bonds_charge <= 3 {
                3
            } else if num_bonds_charge <= 5 {
                5
            } else {
                0
            }
        }
        Atom::S => {
            if num_bonds_charge <= 2 {
                2
            } else if num_bonds_charge <= 4 {
                4
            } else if num_bonds_charge <= 6 {
                6
            } else {
                0
            }
        }
        Atom::F => 1,
        Atom::Cl => 1,
        Atom::Br => 1,
        Atom::I => 1,
        _ => unreachable!(),
    };

    match allowed_valence.checked_sub(num_bonds_charge) {
        Some(x) => Ok(x as u8),
        None => Err("Too many bonds on atom"),
    }
}
