use std::{error::Error, fmt, str::FromStr};

/// An organic atom.
///
/// Organic atoms are a subset of atoms that appear in organic compounds. These are logically
/// differentiated from the catch-all [Atom](crate::Atom) type by their ability to take
/// well-defined valences. This particular organic set is chosen to match the organic atom subset
/// that have their valences auto-inferred by the [SMILES specification](https://daylight.com/dayhtml/doc/theory/theory.smiles.html).
#[derive(Debug, Eq, PartialEq, Clone, Copy, Hash)]
pub enum OrganicAtom {
    H,
    B,
    C,
    N,
    O,
    P,
    S,
    F,
    Cl,
    Br,
    I,
}

impl OrganicAtom {
    /// Returns the allowed valence of the provided organic atom given its allocated valence.
    /// Taken from the [SMILES specification](https://daylight.com/dayhtml/doc/theory/theory.smiles.html), Section 3.2.1.
    /// If the allocated valence exceeds the allowed valence for the atom then this returns [`None`].
    ///
    /// | Atom | Allowed Valence |
    /// | ---- | --------------- |
    /// | B  | 3 |
    /// | C  | 4 |
    /// | O  | 2 |
    /// | N, P | 3, 5 |
    /// | S  | 2, 4, 6 |
    /// | F, Cl, Br, I, H | 1 |
    ///
    /// # Multiple Valences
    /// For atoms with multiple possible allowed valences this uses the lowest possible one based on the allocated valence.
    pub fn allowed_valence(&self, allocated_valence: i8) -> Option<u8> {
        match *self {
            Self::B => Some(3),
            Self::C => Some(4),
            Self::N | Self::P => match allocated_valence {
                x if x <= 3 => Some(3),
                x if x <= 5 => Some(5),
                _ => None,
            },
            Self::O => Some(2),
            Self::S => match allocated_valence {
                x if x <= 2 => Some(2),
                x if x <= 4 => Some(4),
                x if x <= 6 => Some(6),
                _ => None,
            },
            Self::H | Self::F | Self::Cl | Self::Br | Self::I => Some(1),
        }
    }

    /// Calculates unallocated valence.
    ///
    /// The provided allocated valence is taken as is. This allocated valence should include all
    /// bond multiplicities, existing hydrogens, and atomic charges. If the allocated valence
    /// exceeds the allowed valence for this atom then this returns [`None`]. It's possible for the
    /// allocated valence to be negative, since this can be addressed with the addition of
    /// hydrogens. However, the remaining valence must be positive otherwise this returns [`None`].
    ///
    /// ```
    /// # use molecule::OrganicAtom;
    /// assert_eq!(OrganicAtom::S.remaining_valence(-1), Some(3));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 0), Some(2));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 1), Some(1));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 2), Some(0));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 3), Some(1));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 4), Some(0));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 5), Some(1));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 6), Some(0));
    /// assert_eq!(OrganicAtom::S.remaining_valence( 7), None);
    /// ```
    pub fn remaining_valence(&self, allocated_valence: i8) -> Option<u8> {
        self.allowed_valence(allocated_valence)
            .and_then(
                |allowed_valence| match (allowed_valence as i8) - allocated_valence {
                    x if x < 0 => None,
                    x => Some(x as u8),
                },
            )
    }
}

#[derive(Debug)]
pub struct ParseOrganicAtomError;

impl fmt::Display for ParseOrganicAtomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ParseOrganicAtomError")
    }
}

impl Error for ParseOrganicAtomError {}

impl FromStr for OrganicAtom {
    type Err = ParseOrganicAtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(OrganicAtom::H),
            "B" => Ok(OrganicAtom::B),
            "C" => Ok(OrganicAtom::C),
            "N" => Ok(OrganicAtom::N),
            "O" => Ok(OrganicAtom::O),
            "P" => Ok(OrganicAtom::P),
            "S" => Ok(OrganicAtom::S),
            "F" => Ok(OrganicAtom::F),
            "Cl" => Ok(OrganicAtom::Cl),
            "Br" => Ok(OrganicAtom::Br),
            "I" => Ok(OrganicAtom::I),
            _ => Err(Self::Err {}),
        }
    }
}

impl fmt::Display for OrganicAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn remaining_valence() {
        // sulfur tested in doc tests
        assert_eq!(OrganicAtom::N.remaining_valence(-1), Some(4));
        assert_eq!(OrganicAtom::N.remaining_valence(2), Some(1));
        assert_eq!(OrganicAtom::N.remaining_valence(3), Some(0));
        assert_eq!(OrganicAtom::N.remaining_valence(5), Some(0));
        assert_eq!(OrganicAtom::N.remaining_valence(6), None);
        assert_eq!(OrganicAtom::O.remaining_valence(1), Some(1));
        assert_eq!(OrganicAtom::C.remaining_valence(1), Some(3));
        assert_eq!(OrganicAtom::C.remaining_valence(3), Some(1));
        assert_eq!(OrganicAtom::C.remaining_valence(10), None);
        assert_eq!(OrganicAtom::B.remaining_valence(2), Some(1));

        // atoms with the same valence
        for i in -5..10 {
            assert_eq!(
                OrganicAtom::N.remaining_valence(i),
                OrganicAtom::P.remaining_valence(i)
            );
            assert_eq!(
                OrganicAtom::H.remaining_valence(i),
                OrganicAtom::F.remaining_valence(i)
            );
            assert_eq!(
                OrganicAtom::F.remaining_valence(i),
                OrganicAtom::Cl.remaining_valence(i)
            );
            assert_eq!(
                OrganicAtom::Cl.remaining_valence(i),
                OrganicAtom::Br.remaining_valence(i)
            );
            assert_eq!(
                OrganicAtom::Br.remaining_valence(i),
                OrganicAtom::I.remaining_valence(i)
            );
        }
    }
}
