use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::error::Error;
use std::fmt;
use std::hash::Hash;
use std::str::FromStr;

/// A placeholder for a bond.  
///
/// Bond patterns are wildcards that act as placeholders for bonds with particular
/// characteristics. [`BondPattern`]s support all variants of the [`BondType`](crate::bonds::BondType)
/// type and accept `Wildcard` and `Ring` bonds as well. This particular set of
/// patterns was chosen to match the bond primitives as defined by the
/// [SMARTS specification](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html).
#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize, Hash, Copy)]
pub enum BondPattern {
    /// Represents a single bond.  
    ///
    /// Variant used when parsing `-`.
    Single,
    /// Represents a double bond.  
    ///
    /// Variant used when parsing `=`.
    Double,
    /// Represents a double bond.  
    ///
    /// Variant used when parsing `#`.
    Triple,
    /// Represents a double bond.  
    ///
    /// Variant used when parsing `:`.
    Aromatic,
    /// Represents a double bond.  
    ///
    /// Variant used when parsing `~`.
    Wildcard,
    /// Represents a double bond.  
    ///
    /// Variant used when parsing `@`.
    Ring,
}

#[derive(Debug)]
pub struct ParseBondPatternError<T: fmt::Debug + fmt::Display> {
    cause: T,
}

impl<T: fmt::Debug + fmt::Display> fmt::Display for ParseBondPatternError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "can't parse BondPattern from {}", self.cause)
    }
}

impl<T: fmt::Debug + fmt::Display> Error for ParseBondPatternError<T> {}

impl FromStr for BondPattern {
    type Err = ParseBondPatternError<String>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "-" => Ok(Self::Single),
            "=" => Ok(Self::Double),
            "#" => Ok(Self::Triple),
            ":" => Ok(Self::Aromatic),
            "~" => Ok(Self::Wildcard),
            "@" => Ok(Self::Ring),
            _ => Err(Self::Err {
                cause: s.to_string(),
            }),
        }
    }
}

impl fmt::Display for BondPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Self::Single => "-",
            Self::Double => "=",
            Self::Triple => "#",
            Self::Aromatic => ":",
            Self::Wildcard => "~",
            Self::Ring => "@",
        };
        write!(f, "{}", s)
    }
}

impl TryFrom<u8> for BondPattern {
    type Error = ParseBondPatternError<u8>;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(Self::Single),
            2 => Ok(Self::Double),
            3 => Ok(Self::Triple),
            _ => Err(Self::Error { cause: value }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        assert_eq!(BondPattern::Single.to_string(), String::from("-"));
        assert_eq!(BondPattern::Double.to_string(), String::from("="));
        assert_eq!(BondPattern::Triple.to_string(), String::from("#"));
        assert_eq!(BondPattern::Aromatic.to_string(), String::from(":"));
        assert_eq!(BondPattern::Wildcard.to_string(), String::from("~"));
        assert_eq!(BondPattern::Ring.to_string(), String::from("@"));
    }

    #[test]
    fn from_str() {
        assert_eq!("-".parse::<BondPattern>().unwrap(), BondPattern::Single);
        assert_eq!("=".parse::<BondPattern>().unwrap(), BondPattern::Double);
        assert_eq!("#".parse::<BondPattern>().unwrap(), BondPattern::Triple);
        assert_eq!(":".parse::<BondPattern>().unwrap(), BondPattern::Aromatic);
        assert_eq!("~".parse::<BondPattern>().unwrap(), BondPattern::Wildcard);
        assert_eq!("@".parse::<BondPattern>().unwrap(), BondPattern::Ring);
        assert!(" ".parse::<BondPattern>().is_err());
        assert!(".".parse::<BondPattern>().is_err());
        assert!("ABCDAD".parse::<BondPattern>().is_err());
        assert!("-=#:".parse::<BondPattern>().is_err());
    }

    #[test]
    fn try_from_u8() {
        assert_eq!(BondPattern::try_from(1u8).unwrap(), BondPattern::Single);
        assert_eq!(BondPattern::try_from(2u8).unwrap(), BondPattern::Double);
        assert_eq!(BondPattern::try_from(3u8).unwrap(), BondPattern::Triple);
        assert!(BondPattern::try_from(0u8).is_err());
        assert!(BondPattern::try_from(4u8).is_err());
        assert!(BondPattern::try_from(8u8).is_err());
    }
}
