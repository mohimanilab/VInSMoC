//! An alternative to [`std::collections::BTreeSet`] that compares keys for approximate equality

use crate::of64;

use super::{FuzzyMap, Tol};

/// A set of fuzzy floats.
///
/// When querying this set the [`Tol`] given during constructed is used to check for equality.
pub struct FuzzySet {
    map: FuzzyMap<()>,
}

impl FuzzySet {
    /// Construct an empty set with the provided error tolerance.
    pub fn new(tol: Tol) -> Self {
        Self {
            map: FuzzyMap::new(tol),
        }
    }

    /// Construct a set from the provided iterator and error tolerance.
    pub fn from_iter<I: IntoIterator<Item = of64>>(iter: I, tol: Tol) -> Self {
        let mut out = Self::new(tol);
        out.extend(iter);
        out
    }

    /// Get the internal error tolerance.
    pub fn tol(&self) -> Tol {
        self.map.tol()
    }

    /// Remove all values from this set.
    ///
    /// See [`FuzzyMap::clear`] for more details.
    pub fn clear(&mut self) {
        self.map.clear();
    }

    /// Get the value equal to `value` in this set.
    ///
    /// Since this set does fuzzy comparison it's possible that more than one value is equal to the
    /// query. In this case the returned value will be the smallest value in the set equal to the
    /// query. If no equal values are found this returns [`None`].
    pub fn get(&self, value: of64) -> Option<&of64> {
        self.map.get_key_value(value).map(|(k, _)| k)
    }

    /// Get all values equal to `value` in this set.
    ///
    /// If only a single value is desired use [`FuzzyMap::get`].
    pub fn get_all(&self, value: of64) -> impl Iterator<Item = &of64> {
        self.map.get_all_keys(value)
    }

    /// Check if there this set contains `value`.
    pub fn contains(&self, value: of64) -> bool {
        self.map.contains_key(value)
    }

    /// Add `value` to this set.
    ///
    /// If this set did not have an equal element, this returns [`true`]. If there was an equal
    /// element this returns [`false`] and does not update the value.
    pub fn insert(&mut self, value: of64) -> bool {
        self.map.insert(value, ()).is_none()
    }

    /// Get an iterator over references to values in this set.
    pub fn iter(&self) -> Iter<'_> {
        Iter {
            inner: self.map.keys(),
        }
    }
}

impl Extend<of64> for FuzzySet {
    fn extend<T: IntoIterator<Item = of64>>(&mut self, iter: T) {
        for i in iter {
            self.insert(i);
        }
    }
}

/// An iterator over references to values in a [`FuzzySet`].
///
/// Constructed using [`FuzzySet::iter`].
#[derive(Debug, Clone)]
pub struct Iter<'a> {
    inner: super::map::Keys<'a, ()>,
}
inheriting_iter! { Iter<'a>, &'a of64, 'a }
inheriting_iter! { @all Iter<'_> }

/// An iterator over owned values in a [`FuzzySet`].
///
/// Constructed using the [`IntoIterator`] implementation on a [`FuzzySet`].
#[derive(Debug)]
pub struct IntoIter {
    inner: super::map::IntoKeys<()>,
}
inheriting_iter! { IntoIter, of64 }
inheriting_iter! { @all IntoIter }

impl IntoIterator for FuzzySet {
    type Item = of64;
    type IntoIter = IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            inner: self.map.into_keys(),
        }
    }
}
