use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::error::Error;
use std::fmt;
use std::hash::Hash;
use std::str::FromStr;

/// A molecular bond type.
#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize, Hash, Copy, schemars::JsonSchema)]
#[repr(u8)]
pub enum BondType {
    Single = 1,
    Double,
    Triple,
    Aromatic,
}

#[derive(Debug)]
pub struct ParseBondTypeError<T: fmt::Debug + fmt::Display> {
    cause: T,
}

impl<T: fmt::Debug + fmt::Display> fmt::Display for ParseBondTypeError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "can't parse BondType from {}", self.cause)
    }
}

impl<T: fmt::Debug + fmt::Display> Error for ParseBondTypeError<T> {}

impl FromStr for BondType {
    type Err = ParseBondTypeError<String>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" | "-" => Ok(BondType::Single),
            "=" => Ok(BondType::Double),
            "#" => Ok(BondType::Triple),
            ":" => Ok(BondType::Aromatic),
            _ => Err(Self::Err {
                cause: s.to_string(),
            }),
        }
    }
}

impl fmt::Display for BondType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            BondType::Single => "-",
            BondType::Double => "=",
            BondType::Triple => "#",
            BondType::Aromatic => ":",
        };
        write!(f, "{}", s)
    }
}

impl TryFrom<u8> for BondType {
    type Error = ParseBondTypeError<u8>;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(BondType::Single),
            2 => Ok(BondType::Double),
            3 => Ok(BondType::Triple),
            _ => Err(Self::Error { cause: value }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        assert_eq!(BondType::Single.to_string(), String::from("-"));
        assert_eq!(BondType::Double.to_string(), String::from("="));
        assert_eq!(BondType::Triple.to_string(), String::from("#"));
        assert_eq!(BondType::Aromatic.to_string(), String::from(":"));
    }

    #[test]
    fn from_str() {
        assert_eq!("".parse::<BondType>().unwrap(), BondType::Single);
        assert_eq!("-".parse::<BondType>().unwrap(), BondType::Single);
        assert_eq!("=".parse::<BondType>().unwrap(), BondType::Double);
        assert_eq!("#".parse::<BondType>().unwrap(), BondType::Triple);
        assert_eq!(":".parse::<BondType>().unwrap(), BondType::Aromatic);
        assert!(" ".parse::<BondType>().is_err());
        assert!(".".parse::<BondType>().is_err());
        assert!("ABCDAD".parse::<BondType>().is_err());
        assert!("-=#:".parse::<BondType>().is_err());
    }

    #[test]
    fn try_from_u8() {
        assert_eq!(BondType::try_from(1u8).unwrap(), BondType::Single);
        assert_eq!(BondType::try_from(2u8).unwrap(), BondType::Double);
        assert_eq!(BondType::try_from(3u8).unwrap(), BondType::Triple);
        assert!(BondType::try_from(0u8).is_err());
        assert!(BondType::try_from(4u8).is_err());
        assert!(BondType::try_from(8u8).is_err());
    }
}
