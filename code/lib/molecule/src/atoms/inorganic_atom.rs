use std::error::Error;
use std::fmt;
use std::str::FromStr;

#[rustfmt::skip]
#[derive(Debug, Eq, PartialEq)]
pub enum InorganicAtom {
    He, Li, Be, Ne, Na, Mg, Al, Si, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu,
    Zn, Ga, Ge, As, Se, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn,
    Sb, Te, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
    Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa,
    U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg,
    Cn, Nh, Fl, Mc, Lv, Ts, Og,
}

#[derive(Debug)]
pub struct ParseInorganicAtomError;

impl fmt::Display for ParseInorganicAtomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ParseInorganicAtomError")
    }
}

impl Error for ParseInorganicAtomError {}

impl FromStr for InorganicAtom {
    type Err = ParseInorganicAtomError;

    #[rustfmt::skip]
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use InorganicAtom::*;
        match s {
            // valid string
            "He" => Ok(He), "Li" => Ok(Li), "Be" => Ok(Be), "Ne" => Ok(Ne), "Na" => Ok(Na),
            "Mg" => Ok(Mg), "Al" => Ok(Al), "Si" => Ok(Si), "Ar" => Ok(Ar), "K" => Ok(K),
            "Ca" => Ok(Ca), "Sc" => Ok(Sc), "Ti" => Ok(Ti), "V" => Ok(V), "Cr" => Ok(Cr),
            "Mn" => Ok(Mn), "Fe" => Ok(Fe), "Co" => Ok(Co), "Ni" => Ok(Ni), "Cu" => Ok(Cu),
            "Zn" => Ok(Zn), "Ga" => Ok(Ga), "Ge" => Ok(Ge), "As" => Ok(As), "Se" => Ok(Se),
            "Kr" => Ok(Kr), "Rb" => Ok(Rb), "Sr" => Ok(Sr), "Y" => Ok(Y), "Zr" => Ok(Zr),
            "Nb" => Ok(Nb), "Mo" => Ok(Mo), "Tc" => Ok(Tc), "Ru" => Ok(Ru), "Rh" => Ok(Rh),
            "Pd" => Ok(Pd), "Ag" => Ok(Ag), "Cd" => Ok(Cd), "In" => Ok(In), "Sn" => Ok(Sn),
            "Sb" => Ok(Sb), "Te" => Ok(Te), "Xe" => Ok(Xe), "Cs" => Ok(Cs), "Ba" => Ok(Ba),
            "La" => Ok(La), "Ce" => Ok(Ce), "Pr" => Ok(Pr), "Nd" => Ok(Nd), "Pm" => Ok(Pm),
            "Sm" => Ok(Sm), "Eu" => Ok(Eu), "Gd" => Ok(Gd), "Tb" => Ok(Tb), "Dy" => Ok(Dy),
            "Ho" => Ok(Ho), "Er" => Ok(Er), "Tm" => Ok(Tm), "Yb" => Ok(Yb), "Lu" => Ok(Lu),
            "Hf" => Ok(Hf), "Ta" => Ok(Ta), "W" => Ok(W), "Re" => Ok(Re), "Os" => Ok(Os),
            "Ir" => Ok(Ir), "Pt" => Ok(Pt), "Au" => Ok(Au), "Hg" => Ok(Hg), "Tl" => Ok(Tl),
            "Pb" => Ok(Pb), "Bi" => Ok(Bi), "Po" => Ok(Po), "At" => Ok(At), "Rn" => Ok(Rn),
            "Fr" => Ok(Fr), "Ra" => Ok(Ra), "Ac" => Ok(Ac), "Th" => Ok(Th), "Pa" => Ok(Pa),
            "U" => Ok(U), "Np" => Ok(Np), "Pu" => Ok(Pu), "Am" => Ok(Am), "Cm" => Ok(Cm),
            "Bk" => Ok(Bk), "Cf" => Ok(Cf), "Es" => Ok(Es), "Fm" => Ok(Fm), "Md" => Ok(Md),
            "No" => Ok(No), "Lr" => Ok(Lr), "Rf" => Ok(Rf), "Db" => Ok(Db), "Sg" => Ok(Sg),
            "Bh" => Ok(Bh), "Hs" => Ok(Hs), "Mt" => Ok(Mt), "Ds" => Ok(Ds), "Rg" => Ok(Rg),
            "Cn" => Ok(Cn), "Nh" => Ok(Nh), "Fl" => Ok(Fl), "Mc" => Ok(Mc), "Lv" => Ok(Lv),
            "Ts" => Ok(Ts), "Og" => Ok(Og),
            // invalid string
            _ => Err(Self::Err {}),
        }
    }
}
