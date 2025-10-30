//! Approximately equal floating point collections.
//!
//! The core of this module is the [`Tol`] value for both absolute and relative error tolerances.
//! In addition to this tolerance we provide a [`FuzzyMap`] and [`FuzzySet`] for storing
//! collections of approximately equal floating point values (as determined by some [`Tol.`]).

use std::ops::RangeInclusive;

use crate::{of64, ordered::ArchivedOf64};
use rkyv::{out_field, Archive, Deserialize, Serialize};

pub mod map;
pub mod set;

pub use map::FuzzyMap;
pub use set::FuzzySet;

/// The fuzzy comparison error tolerance.
///
/// Fuzzy floating point comparison is used to test approximate equality between two floating point
/// values. In order to test this equality there needs to be some definition of the tolerance for
/// two values to be considered equal. We support two kinds of tolerance: absolute and relative.
/// Each simply wraps a floating point value.
#[derive(Archive, Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Tol {
    /// An absolute tolerance.
    ///
    /// Two values are considered equal using an absolute tolerance if their absolute difference is
    /// smaller than the absolute tolerance.
    Abs(of64),
    /// A relative tolerance.
    ///
    /// Two values are considered equal using a relative tolerance if their absolute difference is
    /// smaller than one of their values scaled down by the relative error.
    Rel(of64),
}

impl Tol {
    /// Checks if two floating point values are equal using this tolerance.
    ///
    /// ```
    /// use balloon::{of64, fuzzy::Tol};
    ///
    /// let f = |x: f64| -> of64 { x.try_into().unwrap() };
    ///
    /// // to be equal f2 must be within ±1 of f1
    /// let tol = Tol::Abs(f(1.0));
    /// assert!(tol.are_equal(f(0.0), f(1.0)));
    /// assert!(tol.are_equal(f(0.0), f(0.5)));
    /// assert!(!tol.are_equal(f(-1.0), f(0.5)));
    /// assert!(!tol.are_equal(f(-1.0), f(0.1)));
    ///
    /// // to be equal f2 must be within ±f1/2
    /// let tol = Tol::Rel(f(0.5));
    /// assert!(tol.are_equal(f(1.0), f(0.75)));
    /// assert!(tol.are_equal(f(1.0), f(0.5)));
    /// assert!(tol.are_equal(f(1.0), f(1.5)));
    /// assert!(tol.are_equal(f(200.0), f(150.0)));
    /// assert!(!tol.are_equal(f(200.0), f(50.0)));
    /// assert!(!tol.are_equal(f(1.0), f(1.6)));
    /// ```
    pub fn are_equal(&self, f1: of64, f2: of64) -> bool {
        match self {
            Self::Abs(tol) => &(f1 - f2).abs() <= tol,
            Self::Rel(tol) => (f1 - f2).abs() <= tol * f1,
        }
    }
    /// Checks if two floating point values are equal using this tolerance.
    ///
    /// ```
    /// use balloon::{of64, fuzzy::Tol};
    ///
    /// let f = |x: f64| -> of64 { x.try_into().unwrap() };
    ///
    /// // to be equal f2 must be within ±1 of f1
    /// let tol = Tol::Abs(f(1.0));
    /// assert!(tol.are_f64_equal(0.0, 1.0));
    /// assert!(tol.are_f64_equal(0.0, 0.5));
    /// assert!(!tol.are_f64_equal(-1.0, 0.5));
    /// assert!(!tol.are_f64_equal(-1.0, 0.1));
    ///
    /// // to be equal f2 must be within ±f1/2
    /// let tol = Tol::Rel(f(0.5));
    /// assert!(tol.are_f64_equal(1.0, 0.75));
    /// assert!(tol.are_f64_equal(1.0, 0.5));
    /// assert!(tol.are_f64_equal(1.0, 1.5));
    /// assert!(tol.are_f64_equal(200.0, 150.0));
    /// assert!(!tol.are_f64_equal(200.0, 50.0));
    /// assert!(!tol.are_f64_equal(1.0, 1.6));
    /// ```
    pub fn are_f64_equal(&self, f1: f64, f2: f64) -> bool {
        match self {
            Self::Abs(tol) => (f1 - f2).abs() <= tol.into_inner(),
            Self::Rel(tol) => (f1 - f2).abs() <= tol.into_inner() * f1,
        }
    }

    /// Convert this tolerance to an absolute tolerance.
    ///
    /// The parameter `f1` is unused if this is a [`Tol::Abs`]. Otherwise, it's used to scaled the
    /// relative error tolerance to retrieve an absolute one.
    pub fn to_absolute(&self, f1: of64) -> Self {
        match self {
            Self::Abs(_) => *self,
            Self::Rel(tol) => Self::Abs(f1 * tol),
        }
    }

    /// Get the inner tolerance.
    pub fn tolerance(&self) -> of64 {
        match self {
            Self::Abs(tol) | Self::Rel(tol) => *tol,
        }
    }

    /// Retrieve the range of values that would be equal to `f1` using this tolerance.
    pub fn range(&self, f1: of64) -> RangeInclusive<of64> {
        let abs = self.to_absolute(f1).tolerance();
        RangeInclusive::new(f1 - abs, f1 + abs)
    }
}

impl ArchivedTol {
    /// Convert this tolerance to an absolute tolerance.
    ///
    /// The parameter `f1` is unused if this is a [`ArchivedTol::Abs`]. Otherwise, it's used to scaled the
    /// relative error tolerance to retrieve an absolute one.
    pub fn to_absolute(&self, f1: ArchivedOf64) -> Self {
        match self {
            Self::Abs(tol) => Self::Abs(*tol),
            Self::Rel(tol) => Self::Abs(f1 * tol),
        }
    }

    /// Get the inner tolerance.
    pub fn tolerance(&self) -> ArchivedOf64 {
        match self {
            Self::Abs(tol) | Self::Rel(tol) => *tol,
        }
    }

    /// Checks if two floating point values are equal using this tolerance.
    ///
    /// ```
    /// use balloon::{of64, ArchivedOf64, fuzzy::ArchivedTol};
    ///
    /// let f = |x: f64| {
    ///    let of = of64::try_from(x).expect("NaN values cannot be ordered");
    ///    ArchivedOf64::from(of)
    /// };
    ///
    /// // to be equal f2 must be within ±1 of f1
    /// let tol = ArchivedTol::Abs(f(1.0));
    /// assert!(tol.are_equal(f(0.0), f(1.0)));
    /// assert!(tol.are_equal(f(0.0), f(0.5)));
    /// assert!(!tol.are_equal(f(-1.0), f(0.5)));
    /// assert!(!tol.are_equal(f(-1.0), f(0.1)));
    ///
    /// // to be equal f2 must be within ±f1/2
    /// let tol = ArchivedTol::Rel(f(0.5));
    /// assert!(tol.are_equal(f(1.0), f(0.75)));
    /// assert!(tol.are_equal(f(1.0), f(0.5)));
    /// assert!(tol.are_equal(f(1.0), f(1.5)));
    /// assert!(tol.are_equal(f(200.0), f(150.0)));
    /// assert!(!tol.are_equal(f(200.0), f(50.0)));
    /// assert!(!tol.are_equal(f(1.0), f(1.6)));
    /// ```
    pub fn are_equal(&self, f1: ArchivedOf64, f2: ArchivedOf64) -> bool {
        match self {
            Self::Abs(tol) => (f1 - f2).abs() <= *tol,
            Self::Rel(tol) => (f1 - f2).abs() <= tol * f1,
        }
    }
    /// Checks if two floating point values (of type f64) are equal using this tolerance.
    ///
    /// ```
    /// use balloon::{of64, ArchivedOf64, fuzzy::ArchivedTol};
    ///
    /// let f = |x: f64| {
    ///    let of = of64::try_from(x).expect("NaN values cannot be ordered");
    ///    ArchivedOf64::from(of)
    /// };
    ///
    /// // to be equal f2 must be within ±1 of f1
    /// let tol = ArchivedTol::Abs(f(1.0));
    /// assert!(tol.are_f64_equal(0.0, 1.0));
    /// assert!(tol.are_f64_equal(0.0, 0.5));
    /// assert!(!tol.are_f64_equal(-1.0, 0.5));
    /// assert!(!tol.are_f64_equal(-1.0, 0.1));
    ///
    /// // to be equal f2 must be within ±f1/2
    /// let tol = ArchivedTol::Rel(f(0.5));
    /// assert!(tol.are_f64_equal(1.0, 0.75));
    /// assert!(tol.are_f64_equal(1.0, 0.5));
    /// assert!(tol.are_f64_equal(1.0, 1.5));
    /// assert!(tol.are_f64_equal(200.0, 150.0));
    /// assert!(!tol.are_f64_equal(200.0, 50.0));
    /// assert!(!tol.are_f64_equal(1.0, 1.6));
    /// ```
    pub fn are_f64_equal(&self, f1: f64, f2: f64) -> bool {
        match self {
            Self::Abs(tol) => (f1 - f2).abs() <= tol.inner(),
            Self::Rel(tol) => (f1 - f2).abs() <= tol.inner() * f1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_are_equal() {
        let f = |x: f64| {
            let of = of64::try_from(x).expect("NaN values cannot be ordered");
            ArchivedOf64::from(of)
        };

        let tol = ArchivedTol::Abs(f(1.0));
        assert!(tol.are_equal(f(0.0), f(1.0)));
        assert!(tol.are_equal(f(0.0), f(0.5)));
        assert!(!tol.are_equal(f(-1.0), f(0.5)));
        assert!(!tol.are_equal(f(-1.0), f(0.1)));

        // to be equal f2 must be within ±f1/2
        let tol = ArchivedTol::Rel(f(0.5));
        assert!(tol.are_equal(f(1.0), f(0.75)));
        assert!(tol.are_equal(f(1.0), f(0.5)));
        assert!(tol.are_equal(f(1.0), f(1.5)));
        assert!(tol.are_equal(f(200.0), f(150.0)));
        assert!(!tol.are_equal(f(200.0), f(50.0)));
        assert!(!tol.are_equal(f(1.0), f(1.6)));
    }

    #[test]
    fn test_tolerance() {
        let tol = Tol::Abs(of64::try_from(2.0).unwrap());
        assert_eq!(tol.tolerance().into_inner(), 2.0);
    }
}
