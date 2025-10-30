use crate::{Adduct, Atom, ChemFormula};
use nom::{
    branch::alt,
    bytes::complete::{tag, take_while_m_n},
    character::complete::digit0,
    combinator::{all_consuming, map, map_res},
    error::VerboseError,
    multi::{many0, many1},
    sequence::{delimited, pair, tuple},
    Err as NomErr, IResult,
};
use rustc_hash::FxHashMap;

use super::forced_charge;
use crate::adducts::{CoefficientSymbol, Symbol};

fn single_atom(input: &str) -> IResult<&str, Atom, VerboseError<&str>> {
    alt((
        map_res(
            take_while_m_n(2, 2, |c: char| c.is_ascii_alphabetic()),
            |s: &str| s.parse::<Atom>(),
        ),
        map_res(
            take_while_m_n(1, 1, |c: char| c.is_ascii_alphabetic()),
            |s: &str| s.parse::<Atom>(),
        ),
    ))(input)
}

fn atom_and_count(input: &str) -> IResult<&str, (Atom, u32), VerboseError<&str>> {
    tuple((
        single_atom,
        map_res(digit0, |num: &str| match num.is_empty() {
            true => Ok(1),
            false => num.parse::<u32>(),
        }),
    ))(input)
}

fn single_formula(input: &str) -> IResult<&str, ChemFormula, VerboseError<&str>> {
    map(many1(atom_and_count), |mut x| {
        let mut cnts: FxHashMap<Atom, u32> = FxHashMap::default();
        x.drain(0..)
            .for_each(|(atom, cnt)| *cnts.entry(atom).or_insert(0) += cnt);
        ChemFormula { cnts }
    })(input)
}

pub fn parse_formula(input: &str) -> Result<ChemFormula, VerboseError<&str>> {
    match all_consuming(single_formula)(input) {
        Ok((_, x)) => Ok(x),
        Err(NomErr::Error(e)) | Err(NomErr::Failure(e)) => Err(e),
        Err(_) => unreachable!(),
    }
}

fn charge_sign_verbose(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    take_while_m_n(1, 1, |c| "+-".contains(c))(input)
}

fn expanded_formula(input: &str) -> IResult<&str, CoefficientSymbol, VerboseError<&str>> {
    map(
        tuple((charge_sign_verbose, digit0, single_formula)),
        |(sign_str, num, formula): (&str, &str, ChemFormula)| {
            let sign = if sign_str == "+" { 1 } else { -1 };

            let coefficient = match num.parse::<i8>() {
                Ok(x) => x * sign,
                _ => sign,
            };
            CoefficientSymbol::new(Symbol::Formula(formula), coefficient)
        },
    )(input)
}

fn init_molecule(input: &str) -> IResult<&str, CoefficientSymbol, VerboseError<&str>> {
    map(pair(digit0, tag("M")), |(num, _): (&str, &str)| {
        let coefficient: i8 = match num.parse() {
            Ok(x) => x,
            _ => 1,
        };
        CoefficientSymbol::new(Symbol::M, coefficient)
    })(input)
}

fn many_expanded_formula(input: &str) -> IResult<&str, Vec<CoefficientSymbol>, VerboseError<&str>> {
    map(
        delimited(
            tag("["),
            pair(init_molecule, many0(expanded_formula)),
            tag("]"),
        ),
        |(init_m, mut added_formulae): (CoefficientSymbol, Vec<CoefficientSymbol>)| {
            added_formulae.insert(0, init_m);
            added_formulae
        },
    )(input)
}

pub fn parse_adduct(input: &str) -> IResult<&str, Adduct, VerboseError<&str>> {
    map(
        all_consuming(pair(many_expanded_formula, forced_charge)),
        |(mols, charge)| Adduct::new(mols, charge),
    )(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_atom_works() {
        assert_eq!(single_atom("H"), Ok(("", Atom::H)));
        assert_eq!(single_atom("HCa"), Ok(("Ca", Atom::H)));
        assert_eq!(single_atom("C2Ca"), Ok(("2Ca", Atom::C)));
        assert_eq!(single_atom("Canonsense"), Ok(("nonsense", Atom::Ca)));
        assert_eq!(single_atom("PtPP"), Ok(("PP", Atom::Pt)));
    }

    #[test]
    fn atom_and_count_works() {
        assert_eq!(atom_and_count("H2"), Ok(("", (Atom::H, 2))));
        assert_eq!(atom_and_count("C"), Ok(("", (Atom::C, 1))));
        assert_eq!(atom_and_count("CH2"), Ok(("H2", (Atom::C, 1))));
        assert_eq!(atom_and_count("CCCCCH2"), Ok(("CCCCH2", (Atom::C, 1))));
        assert!(atom_and_count("2").is_err());
    }

    #[test]
    fn single_formula_works() {
        let cnts: FxHashMap<Atom, u32> = vec![(Atom::C, 2), (Atom::H, 2)].into_iter().collect();
        assert_eq!(single_formula("C2HH").unwrap().1.cnts, cnts);

        let cnts: FxHashMap<Atom, u32> = vec![(Atom::Na, 1), (Atom::C, 12)].into_iter().collect();
        assert_eq!(single_formula("NaCCC10").unwrap().1.cnts, cnts);

        assert!(single_formula("aklsjhfd").is_err());
    }

    #[test]
    fn parse_formula_works() {
        assert_eq!(
            parse_formula("C2HH"),
            Ok(ChemFormula {
                cnts: vec![(Atom::C, 2), (Atom::H, 2)].into_iter().collect()
            })
        );
        assert_eq!(
            parse_formula("NaCCC10"),
            Ok(ChemFormula {
                cnts: vec![(Atom::Na, 1), (Atom::C, 12),].into_iter().collect()
            })
        );

        assert!(parse_formula("NaCCC10A").is_err())
    }

    #[test]
    fn many_expanded_formula_works() {
        let mols: Vec<CoefficientSymbol> = vec![
            CoefficientSymbol::new(Symbol::M, 3),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::N, 1), (Atom::H, 4)].into_iter().collect(),
                }),
                1,
            ),
        ];
        assert_eq!(many_expanded_formula("[3M+NH4]").unwrap().1, mols);

        let mols: Vec<CoefficientSymbol> = vec![
            CoefficientSymbol::new(Symbol::M, 1),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::Na, 1)].into_iter().collect(),
                }),
                2,
            ),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::C, 2), (Atom::H, 3), (Atom::N, 1)]
                        .into_iter()
                        .collect(),
                }),
                -1,
            ),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::H, 1)].into_iter().collect(),
                }),
                1,
            ),
        ];
        assert_eq!(many_expanded_formula("[M+2Na-CH3CN+H]").unwrap().1, mols);

        assert!(many_expanded_formula("[M+abcde]").is_err());

        assert!(many_expanded_formula("[M+Na+]").is_err());
        assert!(many_expanded_formula("[+Na]").is_err());
    }

    #[test]
    fn parse_adduct_works() {
        let mols: Vec<CoefficientSymbol> = vec![
            CoefficientSymbol::new(Symbol::M, 3),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::N, 1), (Atom::H, 4)].into_iter().collect(),
                }),
                1,
            ),
        ];
        let adduct: Adduct = Adduct::new(mols, 1);
        assert_eq!(parse_adduct("[3M+NH4]+").unwrap().1, adduct);

        let mols: Vec<CoefficientSymbol> = vec![
            CoefficientSymbol::new(Symbol::M, 1),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::Na, 1)].into_iter().collect(),
                }),
                2,
            ),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::C, 2), (Atom::H, 3), (Atom::N, 1)]
                        .into_iter()
                        .collect(),
                }),
                -1,
            ),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::H, 1)].into_iter().collect(),
                }),
                1,
            ),
        ];
        let adduct: Adduct = Adduct::new(mols, 2);
        assert_eq!(parse_adduct("[M+2Na-CH3CN+H]+2").unwrap().1, adduct);

        let mols: Vec<CoefficientSymbol> = vec![
            CoefficientSymbol::new(Symbol::M, 1),
            CoefficientSymbol::new(
                Symbol::Formula(ChemFormula {
                    cnts: vec![(Atom::Cl, 1)].into_iter().collect(),
                }),
                2,
            ),
        ];
        let adduct: Adduct = Adduct::new(mols, -2);
        assert_eq!(parse_adduct("[M+2Cl]--").unwrap().1, adduct);

        assert!(parse_adduct("[3M+NH4]").is_err()); // nonzero charge required

        // check M+ adduct
        assert_eq!(
            parse_adduct("[M]+").unwrap().1,
            Adduct::new(vec![CoefficientSymbol::new(Symbol::M, 1)], 1)
        );
    }
}
