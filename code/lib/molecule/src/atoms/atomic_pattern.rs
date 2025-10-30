use super::Atom;
use crate::{Mol, NodeIndex};
use error::Error;
use serde::{Deserialize, Serialize};
use std::{convert::TryFrom, error, fmt, hash::Hash};

/// A placeholder for an atom.  
///
/// Atomic patterns are wildcards that act as placeholders for atoms with particular
/// characteristics. These are differentiated from the catch-all [`Atom`](crate::Atom)
/// type by their ability to take partial information about atoms, including whether
/// they are aromatic or aliphatic, their degree, ring membership and more. This
/// particular set of patterns was chosen to match the atomic primitives as defined
/// by the [SMARTS specification](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html).
#[derive(Debug, Eq, PartialEq, Clone, Copy, Serialize, Deserialize, Hash, PartialOrd, Ord)]
pub enum AtomicPattern {
    /// Represents any arbitrary atom.  
    ///
    /// Variant used when parsing `*`.
    Wildcard,
    /// Represents any arbitrary atom in a ring.  
    ///
    /// Variant used when parsing `r` or `R` followed by no integer.
    RingWildcard,
    /// Defines an arbitrary aromatic atom.  
    ///
    /// Variant used when parsing `a`.
    Aromatic,
    /// Defines an arbitrary aliphatic atom.    
    ///
    /// Variant used when parsing `A`.
    Aliphatic,
    /// Represents an atom with the provided degree or explicit bonds (defaults to 1 if not given).  
    ///
    /// Variant used when parsing from `D<n>` where `<n>` is the degree.
    ExplicitBondCount(u8),
    /// Represents an atom that is attached to the provided number of hydrogens (defaults to aliphatic hydrogen atom if not given).  
    ///
    /// Variant used when parsing `H<n>` where `<n>` is the total number of attached hydrogens.
    HCount(u8),
    /// Represents an atom that is implicitly attached to the provided number of hydrogens (defaults to 1 if not given).  
    ///
    /// Variant used when parsing `h<n>` where `<n>` is the number of implicit hydrogens attached.
    ImplicitHCount(u8),
    /// Represents an atom that is contained in the provided number of [SSSR rings](https://www.daylight.com/dayhtml/doc/theory/theory.mol.html).  
    ///
    /// Variant used when parsing `R<n>` where `<n>` is the number of rings.
    RingMembership(u8),
    /// Represents an atom whose smallest surrounding ring is an [SSSR ring](https://www.daylight.com/dayhtml/doc/theory/theory.mol.html) of the provided size.  
    ///
    /// Variant used when parsing `r<n>` where `<n>` is the size of smallest ring the atom is a part of.
    RingSize(u8),
    /// Represents an atom that has the provided valence (defaults to 1 if not given).  
    ///
    /// Variant used when parsing `v<n>` where `<n>` is the total bond order.
    Valence(u8),
    /// Represents an atom that has the provided total number of bonds, both explicit and implicit (defaults to 1 if not given).  
    ///
    /// Variant used when parsing `X<n>` where `<n>` is the total number of connections.
    TotalBondCount(u8),
    /// Represents an atom that is connected to the provided number of rings (defaults to 1 if not given).  
    ///
    /// Variant used when parsing `x<n>` where `<n>` is the number of connected rings.
    RingConnectivity(u8),
    /// Represents an atom with the provided signed charge (defaults to +1/-1 if not given).  
    ///
    /// Variant used when parsing `+<n>` or `-<n>` where `<n>` is the charge.
    Charge(i8),
    /// Represents an atom with the provided atomic number (no default).  
    ///
    /// Variant used when parsing `#<n>` where `<n>` is the atomic number.
    AtomicNumber(u8),
    /// Defines an atom as specified by the [`Atom`](crate::Atom) argument.  
    ///
    /// Variant used when parsing any aromatic or aliphatic element abbreviation.
    SpecificAtom(Atom),
}

impl AtomicPattern {
    /// Converts an `Atom` into its `AtomicPattern` equivalent.
    pub fn from_atom(atom: Atom) -> Self {
        Self::SpecificAtom(atom)
    }

    /// Checks if the given atomic pattern matches the atom's characteristics
    /// at the index in the provided molecule.
    pub fn matches_atom(&self, mol: &Mol, idx: NodeIndex) -> bool {
        use AtomicPattern::*;
        match self {
            Wildcard => true,
            Aromatic => mol.is_aromatic(idx),
            Aliphatic => !mol.is_aromatic(idx),
            ExplicitBondCount(n) => (mol.graph.edges(idx).count() as u8) == *n,
            HCount(n) => (mol.explicit_h_count(idx) + mol.h_count(&idx)) == *n,
            ImplicitHCount(n) => mol.h_count(&idx) == *n,
            Valence(n) => mol.allowed_valence(idx).map(|v| v == *n).unwrap_or(false),
            TotalBondCount(n) => ((mol.graph.edges(idx).count() as u8) + mol.h_count(&idx)) == *n,
            Charge(n) => mol.charge(&idx) == *n,
            AtomicNumber(n) => (mol.graph[idx].atomic_num() as u8) == *n,
            SpecificAtom(atom) => mol.graph[idx] == *atom,
            RingMembership(_) | RingSize(_) | RingConnectivity(_) | RingWildcard => {
                unimplemented!("Ring handling not implemented")
            }
        }
    }
}

#[derive(Debug)]
pub struct AtomicPatternConversionError {
    atomic_pattern: AtomicPattern,
}

impl Error for AtomicPatternConversionError {}

impl fmt::Display for AtomicPatternConversionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "AtomicPattern {:?} is not a well-defined atom",
            self.atomic_pattern
        )
    }
}

impl TryFrom<AtomicPattern> for Atom {
    type Error = AtomicPatternConversionError;

    fn try_from(atomic_pattern: AtomicPattern) -> Result<Self, Self::Error> {
        match atomic_pattern {
            AtomicPattern::SpecificAtom(atom) => Ok(atom),
            _ => Err(Self::Error { atomic_pattern }),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{parsers::parse_smiles, Atom::*, AtomicPattern::*, NodeIndex};
    #[test]
    fn primitive_matching() {
        let mut mol = parse_smiles(&"C").unwrap();
        let carbon = NodeIndex::new(0);
        assert!(Wildcard.matches_atom(&mol, carbon));
        assert!(Aliphatic.matches_atom(&mol, carbon));
        assert!(ExplicitBondCount(0).matches_atom(&mol, carbon));
        assert!(HCount(4).matches_atom(&mol, carbon));
        assert!(ImplicitHCount(4).matches_atom(&mol, carbon));
        assert!(Valence(4).matches_atom(&mol, carbon));
        assert!(TotalBondCount(4).matches_atom(&mol, carbon));
        assert!(Charge(0).matches_atom(&mol, carbon));
        assert!(AtomicNumber(6).matches_atom(&mol, carbon));
        assert!(SpecificAtom(C).matches_atom(&mol, carbon));

        mol = parse_smiles(&"[Na+].[O-]c1ccccc1").unwrap();
        let sodium = NodeIndex::new(0);
        let oxygen = NodeIndex::new(1);
        let carbon = NodeIndex::new(2);
        assert!(Aliphatic.matches_atom(&mol, sodium));
        assert!(Aliphatic.matches_atom(&mol, oxygen));
        assert!(Aromatic.matches_atom(&mol, carbon));
        assert!(ExplicitBondCount(0).matches_atom(&mol, sodium));
        assert!(ExplicitBondCount(1).matches_atom(&mol, oxygen));
        assert!(ExplicitBondCount(3).matches_atom(&mol, carbon));
        assert!(HCount(0).matches_atom(&mol, sodium));
        assert!(HCount(0).matches_atom(&mol, oxygen));
        assert!(HCount(0).matches_atom(&mol, carbon));
        assert!(ImplicitHCount(0).matches_atom(&mol, sodium));
        assert!(ImplicitHCount(0).matches_atom(&mol, oxygen));
        assert!(ImplicitHCount(0).matches_atom(&mol, carbon));
        assert!(!Valence(0).matches_atom(&mol, sodium));
        assert!(Valence(2).matches_atom(&mol, oxygen));
        assert!(Valence(4).matches_atom(&mol, carbon));
        assert!(TotalBondCount(0).matches_atom(&mol, sodium));
        assert!(TotalBondCount(1).matches_atom(&mol, oxygen));
        assert!(TotalBondCount(3).matches_atom(&mol, carbon));
        assert!(Charge(1).matches_atom(&mol, sodium));
        assert!(Charge(-1).matches_atom(&mol, oxygen));
        assert!(Charge(0).matches_atom(&mol, carbon));
        assert!(AtomicNumber(11).matches_atom(&mol, sodium));
        assert!(AtomicNumber(8).matches_atom(&mol, oxygen));
        assert!(AtomicNumber(6).matches_atom(&mol, carbon));
        assert!(SpecificAtom(Na).matches_atom(&mol, sodium));
        assert!(SpecificAtom(O).matches_atom(&mol, oxygen));
        assert!(SpecificAtom(C).matches_atom(&mol, carbon));

        // Check for explicit hydrogens
        // mol = parse_smiles(&"[CH2]").unwrap();
        // let carbon = NodeIndex::new(0);
        // assert!(SpecificAtom(C).matches_atom(&mol, carbon));
        // assert!(ImplicitHCount(1).matches_atom(&mol, carbon));
        // assert!(HCount(3).matches_atom(&mol, carbon));
        // assert!(ExplicitBondCount(2).matches_atom(&mol, carbon));
        // assert!(TotalBondCount(3).matches_atom(&mol, carbon));

        // TODO: Add more targeted tests for each variant after issue resolved.
    }
}
