//! An alternative to [`std::collections::BTreeMap`] that compares keys for approximate equality
//! using a [`Tol`].

use std::{collections::BTreeMap, ops::Index};

use super::Tol;
use crate::{fuzzy::ArchivedTol, ordered::of64, ordered::ArchivedOf64};
use rkyv::{
    out_field,
    ser::{ScratchSpace, Serializer},
    vec::{ArchivedVec, VecResolver},
    Archive, Deserialize, Fallible, Serialize,
};

/// A map keyed by fuzzy floats.
///
/// When querying this map the [`Tol`] used during its construction is used to check for key
/// equality. Internally this is a [`BTreeMap`] keyed with [`of64`]s.
#[derive(Debug, Clone)]
pub struct FuzzyMap<T> {
    map: BTreeMap<of64, T>,
    tol: Tol,
}

impl<T> FuzzyMap<T> {
    /// Construct an empty map with the provided error tolerance.
    pub fn new(tol: Tol) -> Self {
        Self {
            map: BTreeMap::new(),
            tol,
        }
    }

    /// Construct a map from the provided iterator and error tolerance.
    pub fn from_iter<I: IntoIterator<Item = (of64, T)>>(iter: I, tol: Tol) -> Self {
        let mut out = Self::new(tol);
        out.extend(iter);
        out
    }

    /// Get the internal error tolerance.
    pub fn tol(&self) -> Tol {
        self.tol
    }

    /// Remove all keys and values from this map.
    ///
    /// See [`BTreeMap::clear`] for more details.
    pub fn clear(&mut self) {
        self.map.clear();
    }

    /// Get the value associated with this `key`, if any.
    ///
    /// Since this map does fuzzy comparison for its keys it's possible that more than one key in
    /// the map is equal to the query key. In this case the value associated with the smallest such
    /// key is returned. If no equal keys in the map exist this returns [`None`].
    pub fn get(&self, key: of64) -> Option<&T> {
        self.map.range(self.tol.range(key)).next().map(|(_, v)| v)
    }

    /// Get all values associated with with this `key`.
    ///
    /// If only a single value is desired use [`FuzzyMap::get`].
    pub fn get_all(&self, key: of64) -> impl Iterator<Item = &T> {
        self.map.range(self.tol.range(key)).map(|(_, v)| v)
    }

    /// Get all keys that match with the current `key`.
    ///
    /// If only a single value is desired use [`FuzzyMap::get`].
    pub fn get_all_keys(&self, key: of64) -> impl Iterator<Item = &of64> {
        self.map.range(self.tol.range(key)).map(|(k, _)| k)
    }

    /// Get the key and value associated with this `key`, if any.
    ///
    /// Since this map does fuzzy comparison it's possible that there are multiple hits to the
    /// provided `key`. In those cases the key and value associated with the smallest of the equal
    /// keys in the map is returned. If no equal keys in the map exist this returns [`None`].
    pub fn get_key_value(&self, key: of64) -> Option<(&of64, &T)> {
        self.map.range(self.tol.range(key)).next()
    }

    /// Get all key-value pairs that match with the current `key`.
    ///
    /// If only a single key, value is desired use [`FuzzyMap::get_key_value`].
    pub fn get_all_key_value(&self, key: of64) -> impl Iterator<Item = (&of64, &T)> {
        self.map.range(self.tol.range(key))
    }

    /// Check if this map contains a key equal to `key`.
    pub fn contains_key(&self, key: of64) -> bool {
        self.map.range(self.tol.range(key)).next().is_some()
    }

    /// Get a mutable reference to the value associated with `key`, if any exist.
    ///
    /// Since this map does fuzzy comparison it's possible that there are multiple hits to the
    /// provided `key`. In those cases the mutable value reference returned is the one associated
    /// with the smallest of the equal keys. IF no equal keys in the map exist this returns [`None`].
    pub fn get_mut(&mut self, key: of64) -> Option<&mut T> {
        self.map
            .range_mut(self.tol.range(key))
            .next()
            .map(|(_, v)| v)
    }

    /// Get key and mutable reference to the value associated with `key`, if any exists.
    ///
    /// Since this map does fuzzy comparison it's possible that there are multiple hits to the
    /// provided `key`. In those cases the key, value pair returned is the one associated
    /// with the smallest of the equal keys. IF no equal keys in the map exist this returns [`None`].
    pub fn get_key_value_mut(&mut self, key: of64) -> Option<(&of64, &mut T)> {
        self.map.range_mut(self.tol.range(key)).next()
    }

    /// Add a `key` with `value` to this map.
    ///
    /// If a value with a key fuzzily equal to `key` already exists in this map then the old value
    /// will be overwritten and returned. However, the key is not overwritten. If no such key
    /// already exists this inserts the provided `key` and `value`, returning [`None`].
    pub fn insert(&mut self, key: of64, mut value: T) -> Option<T> {
        if let Some(old_value) = self.get_mut(key) {
            std::mem::swap(old_value, &mut value);
            Some(value)
        } else {
            self.map.insert(key, value);
            None
        }
    }

    /// Convert into an iterator over owned keys.
    pub fn into_keys(self) -> IntoKeys<T> {
        IntoKeys {
            inner: self.map.into_keys(),
        }
    }

    /// Convert into an iterator over owned values.
    pub fn into_values(self) -> IntoValues<T> {
        IntoValues {
            inner: self.map.into_values(),
        }
    }

    /// Get a by-reference iterator over keys and values.
    pub fn iter(&self) -> Iter<'_, T> {
        Iter {
            inner: self.map.iter(),
        }
    }

    /// Get a by-reference iterator over keys.
    pub fn keys(&self) -> Keys<'_, T> {
        Keys {
            inner: self.map.keys(),
        }
    }

    /// Get a by-reference iterator over keys.
    pub fn values(&self) -> Values<'_, T> {
        Values {
            inner: self.map.values(),
        }
    }

    /// Get the number of key-value pairs in this map.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Check if there are any key-value pairs in this map.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<T> Extend<(of64, T)> for FuzzyMap<T> {
    fn extend<I: IntoIterator<Item = (of64, T)>>(&mut self, iter: I) {
        for (k, v) in iter {
            self.insert(k, v);
        }
    }
}

impl<T> Index<of64> for FuzzyMap<T> {
    type Output = T;

    fn index(&self, index: of64) -> &Self::Output {
        self.get(index).expect("could not find key")
    }
}

/// An iterator over references to key-value pairs in a [`FuzzyMap`].
///
/// Constructed using [`FuzzyMap::iter`].
#[derive(Debug, Clone)]
pub struct Iter<'a, T> {
    inner: std::collections::btree_map::Iter<'a, of64, T>,
}
inheriting_iter! { Iter<'a, T>, (&'a of64, &'a T), 'a, T }
inheriting_iter! { @all Iter<'a, T>, 'a, T }

/// An iterator over references to keys in a [`FuzzyMap`].
///
/// Constructed using [`FuzzyMap::keys`].
#[derive(Debug, Clone)]
pub struct Keys<'a, T> {
    inner: std::collections::btree_map::Keys<'a, of64, T>,
}
inheriting_iter! { Keys<'a, T>, &'a of64, 'a, T }
inheriting_iter! { @all Keys<'a, T>, 'a, T }

/// An iterator over references to values in a [`FuzzyMap`].
///
/// Constructed using [`FuzzyMap::values`].
#[derive(Debug, Clone)]
pub struct Values<'a, T> {
    inner: std::collections::btree_map::Values<'a, of64, T>,
}
inheriting_iter! { Values<'a, T>, &'a T, 'a, T }
inheriting_iter! { @all Values<'a, T>, 'a, T }

/// An iterator over owned key-value pairs in a [`FuzzyMap`].
///
/// Constructed using the [`IntoIterator`] implementation on a [`FuzzyMap`].
#[derive(Debug)]
pub struct IntoIter<T> {
    inner: std::collections::btree_map::IntoIter<of64, T>,
}
inheriting_iter! { IntoIter<T>, (of64, T), T }
inheriting_iter! { @all IntoIter<T>, T }

/// An iterator over owned keys in a [`FuzzyMap`].
///
/// Constructed using [`FuzzyMap::into_keys`].
#[derive(Debug)]
pub struct IntoKeys<T> {
    inner: std::collections::btree_map::IntoKeys<of64, T>,
}
inheriting_iter! { IntoKeys<T>, of64, T }
inheriting_iter! { @all IntoKeys<T>, T }

/// An iterator over owned values in a [`FuzzyMap`].
///
/// Constructed using [`FuzzyMap::into_values`].
#[derive(Debug)]
pub struct IntoValues<T> {
    inner: std::collections::btree_map::IntoValues<of64, T>,
}
inheriting_iter! { IntoValues<T>, T, T }
inheriting_iter! { @all IntoValues<T>, T }

impl<T> IntoIterator for FuzzyMap<T> {
    type Item = (of64, T);
    type IntoIter = IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIter {
            inner: self.map.into_iter(),
        }
    }
}

/// [`ArchivedFuzzyMap`] is a type that implements the [`Archive`] trait for [`FuzzyMap`].
///
/// [`FuzzyMap`] is a wrapper over [`BTreeMap`] with the key type set to float.
/// Our implementation of the Archived type does not maintain the tree structure and instead
/// saves sorted vectors of keys and their corresponding values using the default iterators from [`BTreeMap`].
/// These saved vectors can then support the fast range queries we are interesed in using binary search.
pub struct ArchivedFuzzyMap<T: Archive> {
    tol: ArchivedTol,
    keys: ArchivedVec<ArchivedOf64>,
    values: ArchivedVec<<T as Archive>::Archived>,
}

impl<T: Archive> Archive for FuzzyMap<T> {
    type Archived = ArchivedFuzzyMap<T>;
    type Resolver = (<Tol as Archive>::Resolver, VecResolver, VecResolver);

    unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
        let (fp, fo) = out_field!(out.tol);
        self.tol.resolve(pos + fp, resolver.0, fo);
        let (fp, fo) = out_field!(out.keys);
        ArchivedVec::resolve_from_len(self.len(), pos + fp, resolver.1, fo);
        let (fp, fo) = out_field!(out.values);
        ArchivedVec::resolve_from_len(self.len(), pos + fp, resolver.2, fo);
    }
}

impl<S, T> Serialize<S> for FuzzyMap<T>
where
    T: Archive + Serialize<S> + Clone,
    S: ScratchSpace + Serializer + ?Sized,
{
    fn serialize(&self, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        let tol_resolver = self.tol.serialize(serializer)?;
        let keys: Vec<of64> = self.keys().into_iter().cloned().collect::<Vec<of64>>();
        let keys_resolver = keys.serialize(serializer)?;
        let values: Vec<T> = self.values().into_iter().cloned().collect::<Vec<T>>();
        let values_resolver = values.serialize(serializer)?;
        Ok((tol_resolver, keys_resolver, values_resolver))
    }
}

impl<T, D> Deserialize<FuzzyMap<T>, D> for ArchivedFuzzyMap<T>
where
    T: Archive,
    T::Archived: Deserialize<T, D>,
    D: Fallible + ?Sized,
{
    fn deserialize(&self, deserializer: &mut D) -> Result<FuzzyMap<T>, D::Error> {
        let mut result: FuzzyMap<T> = FuzzyMap::new(self.tol.deserialize(deserializer)?);
        for (key, value) in self.keys.iter().zip(self.values.iter()) {
            result.insert(
                key.deserialize(deserializer)?,
                value.deserialize(deserializer)?,
            );
        }
        Ok(result)
    }
}

impl<T: Archive> ArchivedFuzzyMap<T> {
    /// This function returns the tolerance used for the internal [`FuzzyMap`]
    pub fn tol(&self) -> &ArchivedTol {
        &self.tol
    }

    /// This function fetches the first matching key from the [`ArchivedFuzzyMap`]
    ///
    /// This function uses binary search on the saved slice of keys.
    pub fn get(&self, key: ArchivedOf64) -> Option<&<T as Archive>::Archived> {
        let abs = self.tol.to_absolute(key).tolerance();
        let key_lo = key - abs;
        let lo = self.keys.partition_point(|x| x < &key_lo);
        if lo < self.keys.len() && self.tol.are_equal(self.keys[lo], key) {
            Some(&self.values[lo])
        } else {
            None
        }
    }

    /// This function fetches all matching keys from the [`ArchivedFuzzyMap`]
    ///
    /// This function uses binary search on the saved slice of keys.
    pub fn get_all_keys(&self, key: ArchivedOf64) -> impl Iterator<Item = &ArchivedOf64> {
        let abs = self.tol.to_absolute(key).tolerance();
        let key_lo = key - abs;
        let key_hi = key + abs;
        let lo = self.keys.partition_point(|x| x < &key_lo);
        let hi = self.keys.partition_point(|x| x <= &key_hi);
        self.keys[lo..hi].iter()
    }

    /// This function fetches all matching key-value pairs from the [`ArchivedFuzzyMap`]
    ///
    /// This function uses binary search on the saved slice of keys.
    pub fn get_all_key_value(
        &self,
        key: ArchivedOf64,
    ) -> impl Iterator<Item = (&ArchivedOf64, &<T as Archive>::Archived)> {
        let abs = self.tol.to_absolute(key).tolerance();
        let key_lo = key - abs;
        let key_hi = key + abs;
        let lo = self.keys.partition_point(|x| x < &key_lo);
        let hi: usize = self.keys.partition_point(|x| x <= &key_hi);
        (lo..hi).map(|i| (&self.keys[i], &self.values[i]))
    }

    /// This function fetches the first matching key-value pair from the [`ArchivedFuzzyMap`]
    ///
    /// This function uses binary search on the saved slice of keys.
    pub fn get_key_value(
        &self,
        key: ArchivedOf64,
    ) -> Option<(&ArchivedOf64, &<T as Archive>::Archived)> {
        let abs = self.tol.to_absolute(key).tolerance();
        let key_lo = key - abs;
        let lo = self.keys.partition_point(|x| x < &key_lo);
        if lo < self.keys.len() && self.tol.are_equal(self.keys[lo], key) {
            Some((&self.keys[lo], &self.values[lo]))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get() {
        let f = |x: f64| of64::try_from(x).unwrap();
        let af = |x: f64| ArchivedOf64::from(of64::try_from(x).unwrap());
        let tol = Tol::Abs(of64::try_from(2.0).unwrap());
        let mut fmap: FuzzyMap<usize> = FuzzyMap::new(tol);
        fmap.insert(f(2.0), 10);
        fmap.insert(f(2.5), 12);
        fmap.insert(f(3.0), 15);
        fmap.insert(f(4.1), 17);
        let bytes = rkyv::to_bytes::<_, 256>(&fmap).unwrap();
        let archived = unsafe { rkyv::archived_root::<FuzzyMap<usize>>(&bytes) };
        assert_eq!(archived.get(af(1.8)).unwrap(), &15);
        assert_eq!(archived.get(af(2.2)).unwrap(), &15);
        assert_eq!(archived.get(af(4.9)).unwrap(), &17);
        assert_eq!(archived.get(af(6.0)).unwrap(), &17);
        assert_eq!(archived.get(af(6.3)), None);
        assert_eq!(archived.get(af(-0.3)), None);

        assert_eq!(archived.get_key_value(af(2.7)).unwrap(), (&af(2.0), &15));
        assert_eq!(archived.get_key_value(af(-0.5)), None);
        assert_eq!(archived.get_key_value(af(0.1)).unwrap(), (&af(2.0), &15));
        assert_eq!(archived.get_key_value(af(6.0)).unwrap(), (&af(4.1), &17));
        assert_eq!(archived.get_key_value(af(6.2)), None);
    }

    #[test]
    fn test_get_all_keys() {
        let f = |x: f64| of64::try_from(x).unwrap();
        let af = |x: f64| ArchivedOf64::from(of64::try_from(x).unwrap());
        let tol = Tol::Abs(of64::try_from(2.0).unwrap());
        let mut fmap: FuzzyMap<usize> = FuzzyMap::new(tol);
        fmap.insert(f(0.4), 10);
        fmap.insert(f(2.5), 12);
        fmap.insert(f(4.6), 15);
        let bytes = rkyv::to_bytes::<_, 256>(&fmap).unwrap();
        let archived = unsafe { rkyv::archived_root::<FuzzyMap<usize>>(&bytes) };
        assert_eq!(
            archived
                .get_all_keys(af(-2.0))
                .collect::<Vec<&ArchivedOf64>>(),
            Vec::<&ArchivedOf64>::new()
        );
        assert_eq!(
            archived
                .get_all_keys(af(-1.5))
                .collect::<Vec<&ArchivedOf64>>(),
            vec![&af(0.4)]
        );
        assert_eq!(
            archived
                .get_all_keys(af(1.8))
                .collect::<Vec<&ArchivedOf64>>(),
            vec![&af(0.4), &af(2.5)]
        );
        assert_eq!(
            archived
                .get_all_keys(af(2.7))
                .collect::<Vec<&ArchivedOf64>>(),
            vec![&af(2.5), &af(4.6)]
        );
        assert_eq!(
            archived
                .get_all_keys(af(4.55))
                .collect::<Vec<&ArchivedOf64>>(),
            vec![&af(4.6)]
        );
        assert_eq!(
            archived
                .get_all_keys(af(6.6))
                .collect::<Vec<&ArchivedOf64>>(),
            vec![&af(4.6)]
        );
        assert_eq!(
            archived
                .get_all_keys(af(6.7))
                .collect::<Vec<&ArchivedOf64>>(),
            Vec::<&ArchivedOf64>::new()
        );
        assert_eq!(
            archived
                .get_all_key_value(af(2.7))
                .collect::<Vec<(&ArchivedOf64, &u32)>>(),
            vec![(&af(2.5), &12), (&af(4.6), &15)]
        );
    }
}
