#![warn(missing_docs)]

//! A crate for utilities around floating point types.
//!
//! # Ordered Floats
//! This crate provides [`of64`] type, which is a floating point type that cannot be constructed
//! as a [`f64::NAN`] and therefore safely implements [`Ord`], [`Eq`],  and [`Hash`]. Values are
//! checked for NaN upon construction, or optionally unchecked via an `unsafe` constructor. Because
//! of this comparisons are completely unchecked so they don't have performance penalties.
//!
//! There is a known potential issue with soundness though. Since arithmetic operations cannot be
//! unsafe we currently do not enforce that [`of64`]s are not combined in ways that can result in a
//! NaN value, even if none of the operands were NaN.
//!
//! # Fuzzy Floats
//! This crate also provides collections that can store approximately equal floating point values,
//! either with relative or absolute error tolerances. These collections are should have a machine
//! word plus a byte of overhead over [`std::collections::BTreeMap`] and
//! [`std::collections::BTreeSet`]. These are best used for unique collections of floating points
//! that should be approximately equal and either do not have known bounds or have impractically
//! large bounds.
//!
//! # Histograms
//! For collections of bounded floating point types that are approximately equal we provide
//! histograms in [`hist`].

/// A macro to implement iterator traits on structs that wrap an `inner` iterator.
macro_rules! inheriting_iter {
    ($it:ty, $item:ty) => {
        impl Iterator for $it {
            type Item = $item;

            fn next(&mut self) -> Option<Self::Item> {
                self.inner.next()
            }

            fn size_hint(&self) -> (usize, Option<usize>) {
                self.inner.size_hint()
            }
        }
    };
    ($it:ty, $item:ty, $($life:tt),+) => {
        impl<$($life),+> Iterator for $it {
            type Item = $item;

            fn next(&mut self) -> Option<Self::Item> {
                self.inner.next()
            }

            fn size_hint(&self) -> (usize, Option<usize>) {
                self.inner.size_hint()
            }
        }
    };
    (@const $it:ty, $item:ty, $($life:tt),+, $cg:tt, $cgt:ty) => {
        impl<$($life),+, const $cg: $cgt> Iterator for $it {
            type Item = $item;

            fn next(&mut self) -> Option<Self::Item> {
                self.inner.next()
            }

            fn size_hint(&self) -> (usize, Option<usize>) {
                self.inner.size_hint()
            }
        }
    };
    (@exact $it:ty) => {
        impl ExactSizeIterator for $it {
            fn len(&self) -> usize {
                self.inner.len()
            }
        }
    };
    (@double $it:ty) => {
        impl DoubleEndedIterator for $it {
            fn next_back(&mut self) -> Option<Self::Item> {
                self.inner.next_back()
            }
        }
    };
    (@fused $it:ty) => {
        impl std::iter::FusedIterator for $it {}
    };
    (@all $it:ty) => {
        inheriting_iter!{ @exact $it }
        inheriting_iter!{ @double $it }
        inheriting_iter!{ @fused $it }
    };
    (@exact $it:ty, $($life:tt),+) => {
        impl<$($life),+> ExactSizeIterator for $it {
            fn len(&self) -> usize {
                self.inner.len()
            }
        }
    };
    (@double $it:ty, $($life:tt),+) => {
        impl<$($life),+> DoubleEndedIterator for $it {
            fn next_back(&mut self) -> Option<Self::Item> {
                self.inner.next_back()
            }
        }
    };
    (@fused $it:ty, $($life:tt),+) => {
        impl<$($life),+> std::iter::FusedIterator for $it {}
    };
    (@all $it:ty, $($life:tt),+) => {
        inheriting_iter!{ @exact $it, $($life),+ }
        inheriting_iter!{ @double $it, $($life),+ }
        inheriting_iter!{ @fused $it, $($life),+ }
    }
}

pub mod fuzzy;
pub mod hist;
mod ordered;

pub use ordered::{of64, ArchivedOf64, NanError};
