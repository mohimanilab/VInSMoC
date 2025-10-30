//! Adducts for conversion of neutral masses to ionic masses and back again.
//!
//! For a discussion on the definition of an adduct see [`Adduct`]. In addition to the core
//! [`Adduct`] this module contains an [`Symbol`] and a [`CoefficientSymbol`] which build up
//! the equation internal to an adduct.

#![warn(missing_docs)]

use crate::Atom;
use crate::{parsers::parse_adduct, ChemFormula};
use nom::Err as NomErr;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{fmt, str::FromStr};

/// A symbol that can appear in an adduct.
///
/// An adduct is made up of various molecular formulas along with coefficients for each of them.
/// These formulas also can take the special form of `M`, which represents the original target
/// molecule.
#[derive(Debug, Eq, PartialEq, Clone, Hash, Serialize, Deserialize)]
pub enum Symbol {
    /// A variable that represents the original molecule
    M,
    /// Any chemical formula that is given a coefficient in an adduct.
    Formula(ChemFormula),
}

impl Symbol {
    /// Construct an adduct symbol that consists of a single atom.
    pub fn from_atom(atom: Atom) -> Self {
        let mut formula = ChemFormula {
            cnts: FxHashMap::default(),
        };
        formula.cnts.insert(atom, 1);
        Self::Formula(formula)
    }

    /// Calculates the exact mass of this symbol, returning [`0.0`] for [`Symbol::M`].
    pub fn exact_mass(&self) -> f64 {
        match self {
            Self::M => 0.0,
            Self::Formula(formula) => formula.exact_mass(),
        }
    }

    /// Calculates the standard mass of this symbol, returning [`0.0`] for [`Symbol::M`].
    pub fn std_mass(&self) -> f64 {
        match self {
            Self::M => 0.0,
            Self::Formula(formula) => formula.std_mass(),
        }
    }
}

/// The combination of an [`Symbol`] and its coefficient.
#[derive(Debug, Eq, PartialEq, Clone, Hash, Serialize, Deserialize)]
pub struct CoefficientSymbol {
    symbol: Symbol,
    coefficient: i8,
}

impl CoefficientSymbol {
    /// Construct a new [`CoefficientSymbol`].
    pub fn new(symbol: Symbol, coefficient: i8) -> Self {
        Self {
            symbol,
            coefficient,
        }
    }

    /// Get a reference to the expanded formula's symbol.
    pub fn symbol(&self) -> &Symbol {
        &self.symbol
    }

    /// Get the expanded formula's coefficient.
    pub fn coefficient(&self) -> i8 {
        self.coefficient
    }

    /// Calculate the exact mass of this [`CoefficientSymbol`], incorporating the coefficient into
    /// the returned value.
    pub fn exact_mass(&self) -> f64 {
        self.coefficient as f64 * self.symbol.exact_mass()
    }

    /// Calculate the standard mass of this [`CoefficientSymbol`], incorporating the coefficient
    /// into the returned value.
    pub fn std_mass(&self) -> f64 {
        self.coefficient as f64 * self.symbol.std_mass()
    }
}

/// The recipe for an ion formed from some neutral molecule.
///
/// Adducts associate some coefficient with the neutral molecule mass (often denoted `M`) and an
/// equation of added and subtracted moieties, also including coefficients. They are often used to
/// represent ionization of some neutral structure by addition/subtraction of atoms.
///
/// Our `Adduct` is a representation of this equation that can be used to convert a neutral
/// compounds mass to its ionized mass following the adduct equation.
///
/// # Construction
/// It's recommended to use the [`std::str::FromStr`] implementation to parse `Adduct`s from
/// strings. An `Adduct` string consists of an equation containing various molecular formulae with
/// different coefficients. This equation is wrapped with square brackets and suffixed with a
/// charge. Our constructors do not validate the provided charge, instead electing to accept it as
/// is.
/// ```
/// # use molecule::{Adduct, Atom, adducts::{Symbol, CoefficientSymbol}};
/// use approx::assert_abs_diff_eq;
///
/// let m_h = "[M+H]+".parse::<Adduct>().unwrap();
/// assert_eq!(m_h.charge(), 1);
/// assert_eq!(m_h.formulae(), &[
///     CoefficientSymbol::new(Symbol::M, 1),
///     CoefficientSymbol::new(Symbol::from_atom(Atom::H), 1)
/// ]);
/// assert_abs_diff_eq!(
///     m_h.to_ion(Atom::H.exact_mass()),
///     Atom::H.exact_mass() * 2.0 - constants::MASS_ELECTRON,
///     epsilon = 0.00001
/// );
///
/// let m_h = "[M-H]-".parse::<Adduct>().unwrap();
/// assert_eq!(m_h.charge(), -1);
/// assert_eq!(m_h.formulae(), &[
///     CoefficientSymbol::new(Symbol::M, 1),
///     CoefficientSymbol::new(Symbol::from_atom(Atom::H), -1)
/// ]);
/// assert_abs_diff_eq!(
///     m_h.to_ion(Atom::C.exact_mass()),
///     Atom::C.exact_mass() - Atom::H.exact_mass() + constants::MASS_ELECTRON,
///     epsilon = 0.00001
/// );
///
/// let m_h = "[M+H]++".parse::<Adduct>().unwrap();
/// assert_eq!(m_h.charge(), 2);
/// assert_eq!(m_h.formulae(), &[
///     CoefficientSymbol::new(Symbol::M, 1),
///     CoefficientSymbol::new(Symbol::from_atom(Atom::H), 1)
/// ]);
/// assert_abs_diff_eq!(
///     m_h.to_ion(Atom::C.exact_mass()),
///     Atom::C.exact_mass() + Atom::H.exact_mass() - 2.0 * constants::MASS_ELECTRON,
///     epsilon = 0.00001
/// );
/// ```
#[derive(Debug, PartialEq, Clone, Default, Serialize, Deserialize)]
pub struct Adduct {
    formulae: Vec<CoefficientSymbol>,
    charge: i8,
    exact_mass: f64,
}

impl Adduct {
    /// Construct a new adduct.
    ///
    /// This function takes the coefficient-paired adduct symbols as `formulae` and the intended
    /// `charge` of the adduct. Note that this does no checks on the provided parameters. Users
    /// should make sure that the charge is consistent with the adduct's formulae and that the
    /// provided formulae include a [`Symbol::M`] somewhere.
    pub fn new(formulae: Vec<CoefficientSymbol>, charge: i8) -> Self {
        let exact_mass = formulae
            .iter()
            .map(CoefficientSymbol::exact_mass)
            .sum::<f64>()
            - (charge as f64) * constants::MASS_ELECTRON;
        Self {
            formulae,
            charge,
            exact_mass,
        }
    }

    /// Get a reference to the adduct's formuale.
    pub fn formulae(&self) -> &[CoefficientSymbol] {
        &self.formulae
    }

    /// Get the adduct's charge.
    pub fn charge(&self) -> i8 {
        self.charge
    }

    /// Compute the standard mass of the formulae in the adduct.
    ///
    /// The [`Symbol::M`] does not contribute to this mass computation whatsoever.
    pub fn std_mass(&self) -> f64 {
        self.formulae
            .iter()
            .map(CoefficientSymbol::std_mass)
            .sum::<f64>()
            - (self.charge as f64) * constants::MASS_ELECTRON
    }

    /// Compute the exact mass of the formulae in this adduct.
    ///
    /// The [`Symbol::M`] does not contribute to this mass computation whatsoever. Also, this
    /// value is cached during construction of the adduct, so this is a simple `O(1)` getter.
    pub fn exact_mass(&self) -> f64 {
        self.exact_mass
    }

    /// Convert the provided neutral mass into an ion mass using this adduct.
    ///
    /// This uses exact masses throughout.
    #[inline]
    pub fn to_ion(&self, mass: f64) -> f64 {
        self.exact_mass() + mass * (self.formulae[0].coefficient as f64)
    }

    /// Convert the provided neutral mass into a mass-to-charge ratio.
    ///
    /// This uses exact masses throughout.
    #[inline]
    pub fn to_mz(&self, mass: f64) -> f64 {
        self.to_ion(mass) / (self.charge.abs() as f64)
    }

    /// Convert the provided mass-to-charge ratio to a neutral mass.
    ///
    /// This uses exact masses throughout.
    #[inline]
    pub fn to_mass(&self, m_z: f64) -> f64 {
        (((self.charge.abs() as f64) * m_z) - self.exact_mass())
            / (self.formulae[0].coefficient as f64)
    }
}

impl fmt::Display for CoefficientSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.symbol {
            Symbol::M => {
                debug_assert!(self.coefficient > 0);
                if self.coefficient != 1 {
                    write!(f, "{}", self.coefficient)?;
                }
                write!(f, "M")?;
            }
            Symbol::Formula(ref x) => {
                if self.coefficient != 1 {
                    write!(f, "{:+}", self.coefficient)?;
                } else {
                    write!(f, "+")?;
                }

                write!(f, "{}", x)?;
            }
        };

        Ok(())
    }
}

impl fmt::Display for Adduct {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let m = self
            .formulae
            .iter()
            .find(|x| x.symbol == Symbol::M)
            .expect("Impossible adduct with no M symbol found");
        write!(f, "[")?;
        write!(f, "{}", m)?;

        for mol in self.formulae.iter() {
            if mol.symbol == Symbol::M {
                continue;
            }
            write!(f, "{}", mol)?;
        }
        write!(f, "]")?;

        let charge_char = match self.charge > 0 {
            true => '+',
            false => '-',
        };

        for _ in 0..self.charge.abs() {
            write!(f, "{}", charge_char)?;
        }

        Ok(())
    }
}

/// An error that occurred when parsing an adduct from a string.
#[derive(Debug, thiserror::Error)]
#[error("unable to parse adduct: {0}")]
pub struct ParseAdductError(String);

impl<'a> FromStr for Adduct {
    type Err = ParseAdductError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match parse_adduct(s) {
            Ok((_, x)) => Ok(x),
            Err(NomErr::Error(e)) | Err(NomErr::Failure(e)) => {
                Err(ParseAdductError(nom::error::convert_error(s, e)))
            }
            Err(_) => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn to_mz_works() {
        assert_relative_eq!(
            "[M+3H]+3".parse::<Adduct>().unwrap().to_mz(853.33089),
            285.45090645209, // from rdkit
            max_relative = 0.000001
        );
        assert_relative_eq!(
            "[M+Na-2H]-".parse::<Adduct>().unwrap().to_mz(853.33089),
            874.305556,
            max_relative = 0.000001
        );
        assert_relative_eq!(
            "[2M+K]+".parse::<Adduct>().unwrap().to_mz(853.33089),
            1745.624938,
            max_relative = 0.000001
        );
    }

    #[test]
    fn to_mass_works() {
        assert_relative_eq!(
            "[M+3H]+3".parse::<Adduct>().unwrap().to_mass(285.451425032),
            853.33089,
            max_relative = 0.0001
        );
        assert_relative_eq!(
            "[M+Na-2H]-".parse::<Adduct>().unwrap().to_mass(874.305556),
            853.33089,
            max_relative = 0.0001
        );
        assert_relative_eq!(
            "[2M+K]+".parse::<Adduct>().unwrap().to_mass(1745.624938),
            853.33089,
            max_relative = 0.0001
        );
    }

    #[test]
    fn display_expandedformula() {
        let m = CoefficientSymbol {
            symbol: Symbol::M,
            coefficient: 2,
        };
        assert_eq!(m.to_string(), "2M");

        let f = CoefficientSymbol {
            symbol: Symbol::Formula("Na2".parse::<ChemFormula>().unwrap()),
            coefficient: 1,
        };
        assert_eq!(f.to_string(), "+Na2");
        let f = CoefficientSymbol {
            symbol: Symbol::Formula("COOH".parse::<ChemFormula>().unwrap()),
            coefficient: -2,
        };
        assert_eq!(f.to_string(), "-2HCO2");
    }

    #[test]
    fn display_adduct() {
        assert_eq!(
            "[M+3H]+3".parse::<Adduct>().unwrap().to_string(),
            "[M+3H]+++",
        );
        assert_eq!(
            "[M+Na-2H]-".parse::<Adduct>().unwrap().to_string(),
            "[M+Na-2H]-",
        );
        assert_eq!("[2M+K]+".parse::<Adduct>().unwrap().to_string(), "[2M+K]+");
    }
}
