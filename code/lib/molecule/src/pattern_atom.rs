use crate::{atoms::AtomicPattern, Mol, NodeIndex};
use serde::{Deserialize, Serialize};
use std::hash::Hash;

/// A logical operator.  
///
/// Logical operators in the [SMARTS specification](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
/// are used to establish relationships between primitives and form expressions.
/// The operators have the following precedence:  
/// 1. Not operator: `!`
/// 2. High priority And operator: `&`
/// 3. Or operator: `,`
/// 4. Low priority And operator: `;`
#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize, Hash)]
pub enum LogicalOp {
    /// Represents the Not operator.  
    ///
    /// Variant used when parsing `!`.
    Not(Box<PatternAtom>),
    /// Represents the high and low priority And operator.  
    ///
    /// Variant used when parsing `&` or `;`.
    And(Box<PatternAtom>, Box<PatternAtom>),
    /// Represents the Or operator.  
    ///
    /// Variant used when parsing `,`.
    Or(Box<PatternAtom>, Box<PatternAtom>),
}

/// A placeholder for an atom with potentially multiple characteristics.  
///
/// A pattern atom can consist of a simple [`AtomicPattern`]
/// primitive or be an expression to convey multiple aspects of the same atom. This
/// data structure supports recursive and non-recursive [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) strings.
#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize, Hash)]
pub enum PatternAtom {
    /// A primitive is used to represent a single characteristic of an atom
    Primitive(AtomicPattern),
    /// An expression consists of one or more primitives joined by operators
    /// to convey many aspects of the same atom.
    Expression(LogicalOp),
}

impl PatternAtom {
    /// Checks if the given [`PatternAtom`] instance matches an [`Atom`](crate::atoms::Atom) in context
    /// of a molecule. Assumes a semantically valid [`PatternAtom`] instance.
    pub fn matches_atom(&self, mol: &Mol, idx: NodeIndex) -> bool {
        use crate::pattern_atom::LogicalOp::*;
        match self {
            Self::Primitive(atomic_pattern) => {
                AtomicPattern::matches_atom(atomic_pattern, mol, idx)
            }
            Self::Expression(logical_op) => match logical_op {
                Not(pattern_atom) => !(pattern_atom.matches_atom(mol, idx)),
                And(pattern_atom1, pattern_atom2) => {
                    (pattern_atom1.matches_atom(mol, idx)) && (pattern_atom2.matches_atom(mol, idx))
                }
                Or(pattern_atom1, pattern_atom2) => {
                    (pattern_atom1.matches_atom(mol, idx)) || (pattern_atom2.matches_atom(mol, idx))
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        parsers::parse_smiles, Atom::*, AtomicPattern::*, LogicalOp::*, NodeIndex, PatternAtom::*,
    };
    #[test]
    fn matching() {
        let mut mol = parse_smiles(&"C").unwrap();
        let carbon = NodeIndex::new(0);
        assert!(Primitive(Wildcard).matches_atom(&mol, carbon));

        let mut pattern_atom = Expression(Or(
            Box::new(Expression(And(
                Box::new(Primitive(Aliphatic)),
                Box::new(Primitive(Valence(4))),
            ))),
            Box::new(Primitive(AtomicNumber(12))),
        ));
        assert!(pattern_atom.matches_atom(&mol, carbon));
        pattern_atom = Expression(Not(Box::new(Primitive(TotalBondCount(2)))));
        assert!(pattern_atom.matches_atom(&mol, carbon));
        pattern_atom = Expression(Or(
            Box::new(Primitive(Charge(0))),
            Box::new(Primitive(Charge(1))),
        ));
        assert!(pattern_atom.matches_atom(&mol, carbon));

        mol = parse_smiles(&"[Na+].[O-]c1ccccc1").unwrap();
        let sodium = NodeIndex::new(0);
        pattern_atom = Expression(And(
            Box::new(Expression(Or(
                Box::new(Primitive(Aliphatic)),
                Box::new(Primitive(ExplicitBondCount(1))),
            ))),
            Box::new(Primitive(Charge(1))),
        ));
        assert!(pattern_atom.matches_atom(&mol, sodium));
        pattern_atom = Expression(Or(
            Box::new(Expression(And(
                Box::new(Primitive(AtomicNumber(11))),
                Box::new(Primitive(TotalBondCount(1))),
            ))),
            Box::new(Primitive(SpecificAtom(H))),
        ));
        assert!(!pattern_atom.matches_atom(&mol, sodium));
        pattern_atom = Expression(Or(
            Box::new(Expression(Or(
                Box::new(Primitive(SpecificAtom(C))),
                Box::new(Primitive(SpecificAtom(O))),
            ))),
            Box::new(Primitive(SpecificAtom(Na))),
        ));
        assert!(pattern_atom.matches_atom(&mol, sodium));

        mol = parse_smiles("OC(=O)C(Br)(Cl)N").unwrap();
        let oxygen = NodeIndex::new(0);
        let carbon = NodeIndex::new(1);
        let bromine = NodeIndex::new(4);
        let chlorine = NodeIndex::new(5);
        let nitrogen = NodeIndex::new(6);
        assert!(Primitive(SpecificAtom(O)).matches_atom(&mol, oxygen));
        assert!(Primitive(SpecificAtom(C)).matches_atom(&mol, carbon));
        assert!(Primitive(SpecificAtom(Br)).matches_atom(&mol, bromine));
        assert!(Primitive(SpecificAtom(Cl)).matches_atom(&mol, chlorine));
        assert!(Primitive(SpecificAtom(N)).matches_atom(&mol, nitrogen));

        pattern_atom = Expression(And(
            Box::new(Expression(And(
                Box::new(Primitive(AtomicNumber(35))),
                Box::new(Primitive(ExplicitBondCount(1))),
            ))),
            Box::new(Expression(And(
                Box::new(Primitive(HCount(0))),
                Box::new(Primitive(TotalBondCount(1))),
            ))),
        ));
        assert!(pattern_atom.matches_atom(&mol, bromine));
    }
}
