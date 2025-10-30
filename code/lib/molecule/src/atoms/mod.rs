use error::Error;
use serde::{Deserialize, Serialize};
use std::{convert::TryFrom, error, fmt, hash::Hash, str::FromStr};

mod organic_atom;
pub use organic_atom::OrganicAtom;

mod inorganic_atom;
pub use inorganic_atom::InorganicAtom;

mod atomic_pattern;
pub use atomic_pattern::AtomicPattern;

use constants::{PERIODIC_TABLE, STANDARD_ATOMIC_WEIGHTS};

#[rustfmt::skip]
#[derive(Debug, Eq, PartialEq, Clone, Copy, Serialize, Deserialize, Hash, PartialOrd, Ord, schemars::JsonSchema)]
pub enum Atom {
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca,
    Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr,
    Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd,
    Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg,
    Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm,
    Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
}

impl Atom {
    /// Attempt to construct an atom from an atomic number.
    ///
    /// # Errors
    /// If the provided atomic number is not in the valid range (`1..=118`) this will return an
    /// error.
    pub fn try_from_atomic_num(atomic_num: usize) -> Result<Self, InvalidAtomicNum> {
        atomic_num.try_into()
    }

    pub fn atomic_num(&self) -> usize {
        *self as usize + 1
    }

    pub fn std_mass(&self) -> f64 {
        STANDARD_ATOMIC_WEIGHTS[*self as usize]
    }

    pub fn exact_mass(&self) -> f64 {
        PERIODIC_TABLE[*self as usize][0].mass()
    }

    /// Convert this atom into its symbol.
    pub fn symbol(&self) -> &'static str {
        self.into()
    }
}

impl From<OrganicAtom> for Atom {
    #[rustfmt::skip]
    fn from(atom: OrganicAtom) -> Self {
        use OrganicAtom::*;
        match atom {
            H => Atom::H, B => Atom::B, C => Atom::C,
            N => Atom::N, O => Atom::O, P => Atom::P,
            S => Atom::S, F => Atom::F, Cl => Atom::Cl,
            Br => Atom::Br, I => Atom::I,
        }
    }
}

#[derive(Debug)]
pub struct OrganicAtomConversionError {
    atom: Atom,
}

impl fmt::Display for OrganicAtomConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Atom {} is not an organic atom", self.atom)
    }
}

impl Error for OrganicAtomConversionError {}

impl TryFrom<Atom> for OrganicAtom {
    type Error = OrganicAtomConversionError;

    fn try_from(value: Atom) -> Result<Self, Self::Error> {
        use OrganicAtom::*;
        match value {
            Atom::H => Ok(H),
            Atom::B => Ok(B),
            Atom::C => Ok(C),
            Atom::N => Ok(N),
            Atom::O => Ok(O),
            Atom::P => Ok(P),
            Atom::S => Ok(S),
            Atom::F => Ok(F),
            Atom::Cl => Ok(Cl),
            Atom::Br => Ok(Br),
            Atom::I => Ok(I),
            _ => Err(OrganicAtomConversionError { atom: value }),
        }
    }
}

impl From<InorganicAtom> for Atom {
    #[rustfmt::skip]
    fn from(atom: InorganicAtom) -> Self {
        use InorganicAtom::*;
        match atom {
            He => Atom::He, Li => Atom::Li, Be => Atom::Be, Ne => Atom::Ne, Na => Atom::Na,
            Mg => Atom::Mg, Al => Atom::Al, Si => Atom::Si, Ar => Atom::Ar, K => Atom::K,
            Ca => Atom::Ca, Sc => Atom::Sc, Ti => Atom::Ti, V => Atom::V, Cr => Atom::Cr,
            Mn => Atom::Mn, Fe => Atom::Fe, Co => Atom::Co, Ni => Atom::Ni, Cu => Atom::Cu,
            Zn => Atom::Zn, Ga => Atom::Ga, Ge => Atom::Ge, As => Atom::As, Se => Atom::Se,
            Kr => Atom::Kr, Rb => Atom::Rb, Sr => Atom::Sr, Y => Atom::Y, Zr => Atom::Zr,
            Nb => Atom::Nb, Mo => Atom::Mo, Tc => Atom::Tc, Ru => Atom::Ru, Rh => Atom::Rh,
            Pd => Atom::Pd, Ag => Atom::Ag, Cd => Atom::Cd, In => Atom::In, Sn => Atom::Sn,
            Sb => Atom::Sb, Te => Atom::Te, Xe => Atom::Xe, Cs => Atom::Cs, Ba => Atom::Ba,
            La => Atom::La, Ce => Atom::Ce, Pr => Atom::Pr, Nd => Atom::Nd, Pm => Atom::Pm,
            Sm => Atom::Sm, Eu => Atom::Eu, Gd => Atom::Gd, Tb => Atom::Tb, Dy => Atom::Dy,
            Ho => Atom::Ho, Er => Atom::Er, Tm => Atom::Tm, Yb => Atom::Yb, Lu => Atom::Lu,
            Hf => Atom::Hf, Ta => Atom::Ta, W => Atom::W, Re => Atom::Re, Os => Atom::Os,
            Ir => Atom::Ir, Pt => Atom::Pt, Au => Atom::Au, Hg => Atom::Hg, Tl => Atom::Tl,
            Pb => Atom::Pb, Bi => Atom::Bi, Po => Atom::Po, At => Atom::At, Rn => Atom::Rn,
            Fr => Atom::Fr, Ra => Atom::Ra, Ac => Atom::Ac, Th => Atom::Th, Pa => Atom::Pa,
            U => Atom::U, Np => Atom::Np, Pu => Atom::Pu, Am => Atom::Am, Cm => Atom::Cm,
            Bk => Atom::Bk, Cf => Atom::Cf, Es => Atom::Es, Fm => Atom::Fm, Md => Atom::Md,
            No => Atom::No, Lr => Atom::Lr, Rf => Atom::Rf, Db => Atom::Db, Sg => Atom::Sg,
            Bh => Atom::Bh, Hs => Atom::Hs, Mt => Atom::Mt, Ds => Atom::Ds, Rg => Atom::Rg,
            Cn => Atom::Cn, Nh => Atom::Nh, Fl => Atom::Fl, Mc => Atom::Mc, Lv => Atom::Lv,
            Ts => Atom::Ts, Og => Atom::Og,
        }
    }
}

#[derive(Debug)]
pub struct ParseAtomError;

impl fmt::Display for ParseAtomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ParseAtomError")
    }
}

impl Error for ParseAtomError {}

impl FromStr for Atom {
    type Err = ParseAtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let atom = s.parse::<OrganicAtom>();
        let atom = match atom {
            Ok(x) => Ok(x.into()),
            Err(_) => s.parse::<InorganicAtom>().map(|a| a.into()),
        };
        atom.map_err(|_| ParseAtomError {})
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

impl From<Atom> for &'static str {
    #[rustfmt::skip]
    fn from(atom: Atom) -> Self {
        use Atom::*;
        match atom {
            H => "H", He => "He", Li => "Li", Be => "Be", B => "B",
            C => "C", N => "N", O => "O", F => "F", Ne => "Ne",
            Na => "Na", Mg => "Mg", Al => "Al", Si => "Si", P => "P",
            S => "S", Cl => "Cl", Ar => "Ar", K => "K", Ca => "Ca",
            Sc => "Sc", Ti => "Ti", V => "V", Cr => "Cr", Mn => "Mn",
            Fe => "Fe", Co => "Co", Ni => "Ni", Cu => "Cu", Zn => "Zn",
            Ga => "Ga", Ge => "Ge", As => "As", Se => "Se", Br => "Br",
            Kr => "Kr", Rb => "Rb", Sr => "Sr", Y => "Y", Zr => "Zr",
            Nb => "Nb", Mo => "Mo", Tc => "Tc", Ru => "Ru", Rh => "Rh",
            Pd => "Pd", Ag => "Ag", Cd => "Cd", In => "In", Sn => "Sn",
            Sb => "Sb", Te => "Te", I => "I", Xe => "Xe", Cs => "Cs",
            Ba => "Ba", La => "La", Ce => "Ce", Pr => "Pr", Nd => "Nd",
            Pm => "Pm", Sm => "Sm", Eu => "Eu", Gd => "Gd", Tb => "Tb",
            Dy => "Dy", Ho => "Ho", Er => "Er", Tm => "Tm", Yb => "Yb",
            Lu => "Lu", Hf => "Hf", Ta => "Ta", W => "W", Re => "Re",
            Os => "Os", Ir => "Ir", Pt => "Pt", Au => "Au", Hg => "Hg",
            Tl => "Tl", Pb => "Pb", Bi => "Bi", Po => "Po", At => "At",
            Rn => "Rn", Fr => "Fr", Ra => "Ra", Ac => "Ac", Th => "Th",
            Pa => "Pa", U => "U", Np => "Np", Pu => "Pu", Am => "Am",
            Cm => "Cm", Bk => "Bk", Cf => "Cf", Es => "Es", Fm => "Fm",
            Md => "Md", No => "No", Lr => "Lr", Rf => "Rf", Db => "Db",
            Sg => "Sg", Bh => "Bh", Hs => "Hs", Mt => "Mt", Ds => "Ds",
            Rg => "Rg", Cn => "Cn", Nh => "Nh", Fl => "Fl", Mc => "Mc",
            Lv => "Lv", Ts => "Ts", Og => "Og",
        }
    }
}

impl From<&Atom> for &'static str {
    fn from(atom: &Atom) -> Self {
        Self::from(*atom)
    }
}

#[derive(Debug, thiserror::Error)]
#[error("{0} is not a valid atomic number")]
pub struct InvalidAtomicNum(usize);

impl TryFrom<usize> for Atom {
    type Error = InvalidAtomicNum;

    #[rustfmt::skip]
    fn try_from(atomic_num: usize) -> Result<Self, Self::Error> {
        match atomic_num {
            1 => Ok(Atom::H), 2 => Ok(Atom::He), 3 => Ok(Atom::Li), 4 => Ok(Atom::Be),
            5 => Ok(Atom::B), 6 => Ok(Atom::C), 7 => Ok(Atom::N), 8 => Ok(Atom::O),
            9 => Ok(Atom::F), 10 => Ok(Atom::Ne), 11 => Ok(Atom::Na), 12 => Ok(Atom::Mg),
            13 => Ok(Atom::Al), 14 => Ok(Atom::Si), 15 => Ok(Atom::P), 16 => Ok(Atom::S),
            17 => Ok(Atom::Cl), 18 => Ok(Atom::Ar), 19 => Ok(Atom::K), 20 => Ok(Atom::Ca),
            21 => Ok(Atom::Sc), 22 => Ok(Atom::Ti), 23 => Ok(Atom::V), 24 => Ok(Atom::Cr),
            25 => Ok(Atom::Mn), 26 => Ok(Atom::Fe), 27 => Ok(Atom::Co), 28 => Ok(Atom::Ni),
            29 => Ok(Atom::Cu), 30 => Ok(Atom::Zn), 31 => Ok(Atom::Ga), 32 => Ok(Atom::Ge),
            33 => Ok(Atom::As), 34 => Ok(Atom::Se), 35 => Ok(Atom::Br), 36 => Ok(Atom::Kr),
            37 => Ok(Atom::Rb), 38 => Ok(Atom::Sr), 39 => Ok(Atom::Y), 40 => Ok(Atom::Zr),
            41 => Ok(Atom::Nb), 42 => Ok(Atom::Mo), 43 => Ok(Atom::Tc), 44 => Ok(Atom::Ru),
            45 => Ok(Atom::Rh), 46 => Ok(Atom::Pd), 47 => Ok(Atom::Ag), 48 => Ok(Atom::Cd),
            49 => Ok(Atom::In), 50 => Ok(Atom::Sn), 51 => Ok(Atom::Sb), 52 => Ok(Atom::Te),
            53 => Ok(Atom::I), 54 => Ok(Atom::Xe), 55 => Ok(Atom::Cs), 56 => Ok(Atom::Ba),
            57 => Ok(Atom::La), 58 => Ok(Atom::Ce), 59 => Ok(Atom::Pr), 60 => Ok(Atom::Nd),
            61 => Ok(Atom::Pm), 62 => Ok(Atom::Sm), 63 => Ok(Atom::Eu), 64 => Ok(Atom::Gd),
            65 => Ok(Atom::Tb), 66 => Ok(Atom::Dy), 67 => Ok(Atom::Ho), 68 => Ok(Atom::Er),
            69 => Ok(Atom::Tm), 70 => Ok(Atom::Yb), 71 => Ok(Atom::Lu), 72 => Ok(Atom::Hf),
            73 => Ok(Atom::Ta), 74 => Ok(Atom::W), 75 => Ok(Atom::Re), 76 => Ok(Atom::Os),
            77 => Ok(Atom::Ir), 78 => Ok(Atom::Pt), 79 => Ok(Atom::Au), 80 => Ok(Atom::Hg),
            81 => Ok(Atom::Tl), 82 => Ok(Atom::Pb), 83 => Ok(Atom::Bi), 84 => Ok(Atom::Po),
            85 => Ok(Atom::At), 86 => Ok(Atom::Rn), 87 => Ok(Atom::Fr), 88 => Ok(Atom::Ra),
            89 => Ok(Atom::Ac), 90 => Ok(Atom::Th), 91 => Ok(Atom::Pa), 92 => Ok(Atom::U),
            93 => Ok(Atom::Np), 94 => Ok(Atom::Pu), 95 => Ok(Atom::Am), 96 => Ok(Atom::Cm),
            97 => Ok(Atom::Bk), 98 => Ok(Atom::Cf), 99 => Ok(Atom::Es), 100 => Ok(Atom::Fm),
            101 => Ok(Atom::Md), 102 => Ok(Atom::No), 103 => Ok(Atom::Lr), 104 => Ok(Atom::Rf),
            105 => Ok(Atom::Db), 106 => Ok(Atom::Sg), 107 => Ok(Atom::Bh), 108 => Ok(Atom::Hs),
            109 => Ok(Atom::Mt), 110 => Ok(Atom::Ds), 111 => Ok(Atom::Rg), 112 => Ok(Atom::Cn),
            113 => Ok(Atom::Nh), 114 => Ok(Atom::Fl), 115 => Ok(Atom::Mc), 116 => Ok(Atom::Lv),
            117 => Ok(Atom::Ts), 118 => Ok(Atom::Og),
            _ => Err(InvalidAtomicNum(atomic_num)),
        }
    }
}

#[cfg(test)]
mod tests;
