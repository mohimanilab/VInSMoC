//! Drop-in replacements for [`indexmap::IndexMap`] and a [`rustc_hash::FxHasher`] version of
//! [`indexmap::IndexMap`] that serialize as vectors of tuples.
//!
//! This exposes [`IndexTupleMap`] and [`FxIndexTupleMap`] that are suitable for JSON serialization
//! with no extra steps, even when key types are non-trivial. Both types re-implement all methods
//! for their corresponding map types. [`IndexTupleMap`] re-implements all traits implemented on
//! [`indexmap::IndexMap`] except for [`indexmap::map::MutableKeys`] and
//! [`serde::de::IntoDeserializer`].

pub use indexmap::map::{
    Drain, Entry, IntoIter, IntoKeys, IntoValues, Iter, IterMut, Keys, Values, ValuesMut,
};

use std::{
    cmp::Ordering,
    collections::hash_map::RandomState,
    fmt,
    hash::{BuildHasher, BuildHasherDefault, Hash},
    ops::{Index, IndexMut, RangeBounds},
};

use indexmap::{Equivalent, IndexMap};
use rustc_hash::FxHasher;

/// A drop-in replacement for [`indexmap::IndexMap`] that serializes as a vector of tuples with
/// [`serde`].
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(bound(serialize = "K: serde::Serialize, V: serde::Serialize"))]
#[serde(bound(
    deserialize = "K: serde::Deserialize<'de> + Hash + Eq, V: serde::Deserialize<'de>, S: BuildHasher + Default"
))]
#[serde(transparent)]
pub struct IndexTupleMap<K, V, S = RandomState> {
    #[serde(with = "crate::serde")]
    inner: IndexMap<K, V, S>,
}

/// A drop-in replacment for [`indexmap::IndexMap`] with a [`rustc_hash::FxHasher`] that serializes
/// as a vector of tuples with [`serde`].
#[allow(clippy::module_name_repetitions)]
pub type FxIndexTupleMap<K, V> = IndexTupleMap<K, V, BuildHasherDefault<FxHasher>>;

macro_rules! inherit {
    (ctor $name:ident, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        #[must_use]
        pub fn $name($($v:$t),*) -> Self {
            Self {
                inner: IndexMap::$name($($v),*)
            }
        }
    };
    ($name:ident, void, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name(&self, $($v:$t),*) {
            self.inner.$name($($v),*)
        }
    };
    (mut $name:ident, void, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name(&mut self, $($v:$t),*) {
            self.inner.$name($($v),*)
        }
    };
    (owned $name:ident, void, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name(self, $($v:$t),*) {
            self.inner.$name($($v),*)
        }
    };
    ($name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name(&self, $($v:$t),*) -> $ret {
            self.inner.$name($($v),*)
        }
    };
    (mut $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name(&mut self, $($v:$t),*) -> $ret {
            self.inner.$name($($v),*)
        }
    };
    (owned $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name(self, $($v:$t),*) -> $ret {
            self.inner.$name($($v),*)
        }
    };
    (get $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name<Q>(&self, $($v:$t),*) -> $ret
        where
            Q: Hash + Equivalent<K> + ?Sized,
        {
            self.inner.$name($($v),*)
        }
    };
    (mut get $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name<Q>(&mut self, $($v:$t),*) -> $ret
        where
            Q: Hash + Equivalent<K> + ?Sized,
        {
            self.inner.$name($($v),*)
        }
    };
    (sort $name:ident, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name<F>(&mut self, $($v:$t),*)
        where
            F: FnMut(&K, &V, &K, &V) -> Ordering,
        {
            self.inner.$name($($v),*)
        }
    };
    (owned sort $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`IndexMap::", stringify!($name), "`].")]
        pub fn $name<F>(self, $($v:$t),*) -> $ret
        where
            F: FnMut(&K, &V, &K, &V) -> Ordering,
        {
            self.inner.$name($($v),*)
        }
    };
}

impl<K, V> IndexTupleMap<K, V> {
    inherit! { ctor new, }
    inherit! { ctor with_capacity, n: usize }
}

impl<K, V, S> IndexTupleMap<K, V, S> {
    inherit! { ctor with_hasher, hash_builder: S }
    inherit! { ctor with_capacity_and_hasher, capacity: usize, hasher: S }
    inherit! { capacity, usize, }
    inherit! { hasher, &S, }
    inherit! { len, usize, }
    inherit! { is_empty, bool, }
    inherit! { iter, Iter<'_, K, V>, }
    inherit! { mut iter_mut, IterMut<'_, K, V>, }
    inherit! { keys, Keys<'_, K, V>, }
    inherit! { owned into_keys, IntoKeys<K, V>, }
    inherit! { values, Values<'_, K, V>, }
    inherit! { mut values_mut, ValuesMut<'_, K, V>, }
    inherit! { owned into_values, IntoValues<K, V>, }
    inherit! { mut clear, void, }
    inherit! { mut truncate, void, len: usize}

    /// See [`IndexMap::drain`].
    pub fn drain<R>(&mut self, range: R) -> Drain<'_, K, V>
    where
        R: RangeBounds<usize>,
    {
        self.inner.drain(range)
    }

    /// See [`IndexMap::split_off`].
    pub fn split_off(&mut self, at: usize) -> Self
    where
        S: Clone,
    {
        Self {
            inner: self.inner.split_off(at),
        }
    }

    inherit! { get_index, Option<(&K, &V)>, index: usize }
    inherit! { mut get_index_mut, Option<(&mut K, &mut V)>, index: usize }
    inherit! { first, Option<(&K, &V)>, }
    inherit! { mut first_mut, Option<(&K, &mut V)>, }
    inherit! { last, Option<(&K, &V)>, }
    inherit! { mut last_mut, Option<(&K, &mut V)>, }
    inherit! { mut swap_remove_index, Option<(K, V)>, index: usize }
    inherit! { mut shift_remove_index, Option<(K, V)>, index: usize }
    inherit! { mut move_index, void, from: usize, to: usize }
    inherit! { mut swap_indices, void, a: usize, b: usize }
}

impl<K, V, S> IndexTupleMap<K, V, S>
where
    K: Hash + Eq,
    S: BuildHasher,
{
    inherit! { mut reserve, void, additional: usize }
    inherit! { mut shrink_to_fit, void, }
    inherit! { mut shrink_to, void, min_capacity: usize }
    inherit! { mut insert, Option<V>, key: K, value: V }
    inherit! { mut insert_full, (usize, Option<V>), key: K, value: V }
    inherit! { mut entry, Entry<'_, K, V>, key: K }
    inherit! { get contains_key, bool, key: &Q }
    inherit! { get get, Option<&V>, key: &Q }
    inherit! { get get_key_value, Option<(&K, &V)>, key: &Q }
    inherit! { get get_full, Option<(usize, &K, &V)>, key: &Q }
    inherit! { get get_index_of, Option<usize>, key: &Q }
    inherit! { mut get get_mut, Option<&mut V>, key: &Q }
    inherit! { mut get get_full_mut, Option<(usize, &K, &mut V)>, key: &Q }
    inherit! { mut get remove, Option<V>, key: &Q }
    inherit! { mut get remove_entry, Option<(K, V)>, key: &Q }
    inherit! { mut get swap_remove, Option<V>, key: &Q }
    inherit! { mut get swap_remove_entry, Option<(K, V)>, key: &Q }
    inherit! { mut get swap_remove_full, Option<(usize, K, V)>, key: &Q }
    inherit! { mut get shift_remove, Option<V>, key: &Q }
    inherit! { mut get shift_remove_entry, Option<(K, V)>, key: &Q }
    inherit! { mut get shift_remove_full, Option<(usize, K, V)>, key: &Q }
    inherit! { mut pop, Option<(K, V)>, }

    /// See [`IndexMap::retain`].
    pub fn retain<F>(&mut self, keep: F)
    where
        F: FnMut(&K, &mut V) -> bool,
    {
        self.inner.retain(keep);
    }

    inherit! { sort sort_by, cmp: F }
    inherit! { owned sort sorted_by, IntoIter<K, V>, cmp: F }
    inherit! { sort sort_unstable_by, cmp: F }
    inherit! { owned sort sorted_unstable_by, IntoIter<K, V>, cmp: F }
    inherit! { mut reverse, void, }
}

impl<K, V, S> IndexTupleMap<K, V, S>
where
    K: Hash + Eq + Ord,
    S: BuildHasher,
{
    inherit! { mut sort_keys, void, }
    inherit! { mut sort_unstable_keys, void, }
}

// trait implementations

// skipping arbitrary

impl<K, V, S> fmt::Debug for IndexTupleMap<K, V, S>
where
    K: fmt::Debug,
    V: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.inner.fmt(f)
    }
}

impl<K, V, S: Default> Default for IndexTupleMap<K, V, S> {
    fn default() -> Self {
        Self {
            inner: IndexMap::default(),
        }
    }
}

// serde::Deserialize implemented as derive macro

impl<'a, K, V, S> Extend<(&'a K, &'a V)> for IndexTupleMap<K, V, S>
where
    K: Hash + Eq + Copy + 'a,
    V: Copy + 'a,
    S: BuildHasher,
{
    fn extend<T: IntoIterator<Item = (&'a K, &'a V)>>(&mut self, iter: T) {
        self.inner.extend(iter);
    }
}

impl<K, V, S> Extend<(K, V)> for IndexTupleMap<K, V, S>
where
    K: Hash + Eq,
    S: BuildHasher,
{
    fn extend<T: IntoIterator<Item = (K, V)>>(&mut self, iter: T) {
        self.inner.extend(iter);
    }
}

impl<K, V, const N: usize> From<[(K, V); N]> for IndexTupleMap<K, V, RandomState>
where
    K: Hash + Eq,
{
    fn from(arr: [(K, V); N]) -> Self {
        Self { inner: arr.into() }
    }
}

impl<K, V, S> FromIterator<(K, V)> for IndexTupleMap<K, V, S>
where
    K: Hash + Eq,
    S: BuildHasher + Default,
{
    fn from_iter<T: IntoIterator<Item = (K, V)>>(iter: T) -> Self {
        Self {
            inner: iter.into_iter().collect(),
        }
    }
}

impl<K, V, Q, S> Index<&Q> for IndexTupleMap<K, V, S>
where
    Q: Hash + Equivalent<K> + ?Sized,
    K: Hash + Eq,
    S: BuildHasher,
{
    type Output = V;

    fn index(&self, index: &Q) -> &Self::Output {
        self.inner.index(index)
    }
}

impl<K, V, S> Index<usize> for IndexTupleMap<K, V, S> {
    type Output = V;

    fn index(&self, index: usize) -> &Self::Output {
        self.inner.index(index)
    }
}

impl<K, V, Q, S> IndexMut<&Q> for IndexTupleMap<K, V, S>
where
    Q: Hash + Equivalent<K> + ?Sized,
    K: Hash + Eq,
    S: BuildHasher,
{
    fn index_mut(&mut self, index: &Q) -> &mut Self::Output {
        self.inner.index_mut(index)
    }
}

impl<K, V, S> IndexMut<usize> for IndexTupleMap<K, V, S> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.inner.index_mut(index)
    }
}

// IntoDeserializer is skipped

impl<'a, K, V, S> IntoIterator for &'a IndexTupleMap<K, V, S> {
    type Item = (&'a K, &'a V);

    type IntoIter = Iter<'a, K, V>;

    fn into_iter(self) -> Self::IntoIter {
        (&self.inner).into_iter()
    }
}

impl<'a, K, V, S> IntoIterator for &'a mut IndexTupleMap<K, V, S> {
    type Item = (&'a K, &'a mut V);

    type IntoIter = IterMut<'a, K, V>;

    fn into_iter(self) -> Self::IntoIter {
        (&mut self.inner).into_iter()
    }
}

impl<K, V, S> IntoIterator for IndexTupleMap<K, V, S> {
    type Item = (K, V);

    type IntoIter = IntoIter<K, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter()
    }
}

// have to skip mutablekeys

impl<K, V1, S1, V2, S2> PartialEq<IndexTupleMap<K, V2, S2>> for IndexTupleMap<K, V1, S1>
where
    K: Hash + Eq,
    V1: PartialEq<V2>,
    S1: BuildHasher,
    S2: BuildHasher,
{
    fn eq(&self, other: &IndexTupleMap<K, V2, S2>) -> bool {
        self.inner.eq(&other.inner)
    }
}

impl<K, V, S> Eq for IndexTupleMap<K, V, S>
where
    K: Eq + Hash,
    V: Eq,
    S: BuildHasher,
{
}

impl<K, V, S> From<IndexMap<K, V, S>> for IndexTupleMap<K, V, S> {
    fn from(inner: IndexMap<K, V, S>) -> Self {
        Self { inner }
    }
}

impl<K, V, S> From<IndexTupleMap<K, V, S>> for IndexMap<K, V, S> {
    fn from(outer: IndexTupleMap<K, V, S>) -> Self {
        outer.inner
    }
}

impl<K, V, S> AsRef<IndexMap<K, V, S>> for IndexTupleMap<K, V, S> {
    fn as_ref(&self) -> &IndexMap<K, V, S> {
        &self.inner
    }
}

#[cfg(test)]
mod tests {
    // for the macros
    #![allow(unused_mut)]

    use super::*;

    type FxIndexMap<K, V> = IndexMap<K, V, BuildHasherDefault<FxHasher>>;

    macro_rules! method_test {
        (ctor $name:ident, $($params:expr),*) => {
            let tuple_map = IndexTupleMap::<u64, u64>::$name($($params),*);
            let map = IndexMap::<u64, u64>::$name($($params),*);
            assert_eq!(tuple_map.inner, map);
        };
        (other ctor $name:ident, $($params:expr),*) => {
            let tuple_map = IndexTupleMap::<u64, u64, _>::$name($($params),*);
            let map = IndexMap::<u64, u64, _>::$name($($params),*);
            assert_eq!(tuple_map.inner, map);
        };
        ($name:ident, $($params:expr),*) => {
            let mut tuple_map = IndexTupleMap::<u64, u64>::default();
            let mut map = IndexMap::<u64, u64>::default();
            assert_eq!(tuple_map.$name($($params),*), map.$name($($params),*));
            assert_eq!(tuple_map.inner, map);

            let mut tuple_map = IndexTupleMap::<u64, u64>::default();
            tuple_map.insert(1,1);
            let mut map = IndexMap::<u64, u64>::default();
            map.insert(1,1);
            assert_eq!(tuple_map.$name($($params),*), map.$name($($params),*));
            assert_eq!(tuple_map.inner, map);


            let mut tuple_map = IndexTupleMap::<u64, u64>::default();
            tuple_map.insert(1,1);
            tuple_map.insert(2,3);
            tuple_map.insert(4,5);
            tuple_map.insert(6,7);
            let mut map = IndexMap::<u64, u64>::default();
            map.insert(1,1);
            map.insert(2,3);
            map.insert(4,5);
            map.insert(6,7);
            assert_eq!(tuple_map.$name($($params),*), map.$name($($params),*));
            assert_eq!(tuple_map.inner, map);

            let mut tuple_map = FxIndexTupleMap::<String, i8>::default();
            tuple_map.insert("hello".into(), 5);
            tuple_map.insert("bye".into(), 5);
            tuple_map.insert("hello".into(), 3);
            let mut map = FxIndexMap::<String, i8>::default();
            map.insert("hello".into(), 5);
            map.insert("bye".into(), 5);
            map.insert("hello".into(), 3);
            assert_eq!(tuple_map.$name($($params),*), map.$name($($params),*));
            assert_eq!(tuple_map.inner, map);

        };
    }

    #[test]
    fn new() {
        method_test! { ctor new,  }
    }

    #[test]
    fn with_capacity() {
        method_test! { ctor with_capacity, 10 }
        method_test! { ctor with_capacity, 20 }
        method_test! { ctor with_capacity, 0 }
    }

    #[test]
    fn with_hasher() {
        method_test! { other ctor with_hasher, BuildHasherDefault::<FxHasher>::default() }
    }

    #[test]
    fn with_capacity_and_hasher() {
        method_test! { other ctor with_capacity_and_hasher, 10, BuildHasherDefault::<FxHasher>::default() }
        method_test! { other ctor with_capacity_and_hasher, 20, BuildHasherDefault::<FxHasher>::default() }
        method_test! { other ctor with_capacity_and_hasher, 0, BuildHasherDefault::<FxHasher>::default() }
    }

    #[test]
    fn capacity() {
        method_test! { capacity, }
        let tuple_map = IndexTupleMap::<i8, i8>::with_capacity(0);
        assert_eq!(tuple_map.capacity(), 0);
        let tuple_map = IndexTupleMap::<i8, i8>::with_capacity(100);
        assert!(tuple_map.capacity() >= 100);
        let tuple_map = IndexTupleMap::<i8, i8, _>::with_capacity_and_hasher(
            0,
            BuildHasherDefault::<FxHasher>::default(),
        );
        assert_eq!(tuple_map.capacity(), 0);
        let tuple_map = IndexTupleMap::<i8, i8, _>::with_capacity_and_hasher(
            10,
            BuildHasherDefault::<FxHasher>::default(),
        );
        assert!(tuple_map.capacity() >= 10);
    }

    // hasher doesn't impl eq always, skipping tests

    // tests keys, values, values_mut, iter, iter_mut, into_keys, into_values, into_iter
    #[test]
    fn iterators() {
        // keys
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.keys().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.keys().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![&1, &3, &5]);
        assert_eq!(tuple_iter, hash_iter);

        // values
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.values().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.values().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![&2, &4, &6]);
        assert_eq!(tuple_iter, hash_iter);

        // iter
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.iter().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.iter().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![(&1, &2), (&3, &4), (&5, &6)]);
        assert_eq!(tuple_iter, hash_iter);

        // values_mut
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.values_mut().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.values_mut().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![&mut 2, &mut 4, &mut 6]);
        assert_eq!(tuple_iter, hash_iter);

        // iter_mut
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.iter_mut().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.iter_mut().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![(&1, &mut 2), (&3, &mut 4), (&5, &mut 6)]);
        assert_eq!(tuple_iter, hash_iter);

        // into_keys
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.into_keys().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.into_keys().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![1, 3, 5]);
        assert_eq!(tuple_iter, hash_iter);

        // into_values
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.into_values().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.into_values().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![2, 4, 6]);
        assert_eq!(tuple_iter, hash_iter);

        // into_iter
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        let mut tuple_iter = tuple_map.into_iter().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.into_iter().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![(1, 2), (3, 4), (5, 6)]);
        assert_eq!(tuple_iter, hash_iter);
    }

    #[test]
    fn len() {
        method_test! { len, }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        assert_eq!(map.len(), 3);
        map.insert(7, 8);
        assert_eq!(map.len(), 4);
        map.insert(7, 9);
        assert_eq!(map.len(), 4);
        map.remove(&1);
        assert_eq!(map.len(), 3);
    }

    #[test]
    fn is_empty() {
        method_test! { is_empty, }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        assert!(!map.is_empty());
        map.remove(&1);
        assert!(!map.is_empty());
        map.remove(&3);
        assert!(!map.is_empty());
        map.remove(&5);
        assert!(map.is_empty());
    }

    #[test]
    fn clear() {
        method_test! { clear, }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        map.clear();
        assert_eq!(map, IndexTupleMap::<i32, i32>::default());
    }

    #[test]
    fn truncate() {
        method_test! { truncate, 10 }
        method_test! { truncate, 2 }
        method_test! { truncate, 1 }
    }

    #[test]
    fn drain() {
        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        let mut tuple_drained = tuple_map.drain(..).collect::<Vec<_>>();
        tuple_drained.sort_by(|(k1, v1), (k2, v2)| k1.cmp(k2).then(v1.cmp(v2)));
        let mut hash_drained = hash_map.drain(..).collect::<Vec<_>>();
        hash_drained.sort_by(|(k1, v1), (k2, v2)| k1.cmp(k2).then(v1.cmp(v2)));
        assert_eq!(tuple_drained, hash_drained);
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.len(), 0);
    }

    #[test]
    fn split_off() {
        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();

        assert_eq!(tuple_map.split_off(1).inner, hash_map.split_off(1));
        assert_eq!(tuple_map.inner, hash_map);
    }

    #[test]
    fn get_index() {
        method_test! { get_index, 0 }
        method_test! { get_index, 1 }
        method_test! { get_index, 2 }
        method_test! { get_index, 100 }
    }

    #[test]
    fn get_index_mut() {
        method_test! { get_index_mut, 0 }
        method_test! { get_index_mut, 1 }
        method_test! { get_index_mut, 2 }
        method_test! { get_index_mut, 100 }
    }

    #[test]
    fn first() {
        method_test! { first, }
    }

    #[test]
    fn first_mut() {
        method_test! { first, }
    }

    #[test]
    fn last() {
        method_test! { last, }
    }

    #[test]
    fn last_mut() {
        method_test! { last, }
    }

    #[test]
    fn swap_remove_index() {
        method_test! { swap_remove_index, 0 }
    }

    #[test]
    fn shift_remove_index() {
        method_test! { swap_remove_index, 0 }
    }

    #[test]
    fn move_index() {
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&map.inner, &hash_map);
        map.move_index(0, 1);
        hash_map.move_index(0, 1);
        assert_eq!(&map.inner, &hash_map);
    }

    #[test]
    fn swap_indices() {
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&map.inner, &hash_map);
        map.swap_indices(0, 1);
        hash_map.swap_indices(0, 1);
        assert_eq!(&map.inner, &hash_map);
    }

    #[test]
    fn reserve() {
        method_test! { reserve, 10 }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        map.reserve(20);
        assert!(map.capacity() >= 23);
        map.reserve(10);
        assert!(map.capacity() >= 23);
        map.reserve(30);
        assert!(map.capacity() >= 33);
        let mut map = IndexTupleMap::<i32, i32>::new();
        assert_eq!(map.capacity(), 0);
        map.reserve(10);
        assert!(map.capacity() >= 10);
    }

    #[test]
    fn shrink_to_fit() {
        method_test! { shrink_to_fit, }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        map.reserve(100);
        assert!(map.capacity() >= 103);
        map.shrink_to_fit();
        assert!(map.capacity() >= 3);
        let mut map = IndexTupleMap::<u64, u64>::new();
        map.insert(10, 20);
        map.insert(20, 20);
        map.reserve(100);
        assert!(map.capacity() >= 102);
        map.shrink_to_fit();
        assert!(map.capacity() >= 2);
    }

    #[test]
    fn shrink_to() {
        method_test! { shrink_to, 5 }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        map.reserve(100);
        assert!(map.capacity() >= 103);
        map.shrink_to(50);
        assert!(map.capacity() >= 50);
        map.shrink_to(0);
        assert!(map.capacity() >= 3);
        let mut map = IndexTupleMap::<i32, i32>::new();
        map.insert(1, 1);
        map.reserve(100);
        assert!(map.capacity() >= 100);
        map.shrink_to(50);
        assert!(map.capacity() >= 50);
        map.shrink_to(0);
        assert!(map.capacity() >= 1);
    }

    #[test]
    fn insert_remove() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.remove(&1), hash_map.remove(&1));
        assert_eq!(tuple_map.remove(&3), hash_map.remove(&3));
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn insert_full() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert_full(10, 20), hash_map.insert_full(10, 20));
        assert_eq!(tuple_map.insert_full(1, 3), hash_map.insert_full(1, 3));
        assert_eq!(tuple_map.insert_full(2, 2), hash_map.insert_full(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn entry() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.entry(1).key(), hash_map.entry(1).key());
        assert_eq!(
            tuple_map.entry(1).or_default(),
            hash_map.entry(1).or_default()
        );
        assert_eq!(tuple_map.entry(3).key(), hash_map.entry(3).key());
        assert_eq!(
            tuple_map.entry(3).or_default(),
            hash_map.entry(3).or_default()
        );
        assert_eq!(tuple_map.entry(4).key(), hash_map.entry(4).key());
        assert_eq!(
            tuple_map.entry(3).or_insert(20),
            hash_map.entry(3).or_insert(20)
        );
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn contains_key() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.contains_key(&1), hash_map.contains_key(&1));
        assert_eq!(tuple_map.contains_key(&2), hash_map.contains_key(&2));
        assert_eq!(tuple_map.contains_key(&3), hash_map.contains_key(&3));
        assert_eq!(tuple_map.contains_key(&4), hash_map.contains_key(&4));
        assert_eq!(tuple_map.contains_key(&5), hash_map.contains_key(&5));

        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.contains_key("a"), hash_map.contains_key("a"));
        assert_eq!(tuple_map.contains_key("b"), hash_map.contains_key("b"));
        assert_eq!(tuple_map.contains_key("c"), hash_map.contains_key("c"));
        assert_eq!(tuple_map.contains_key("d"), hash_map.contains_key("d"));
        assert_eq!(tuple_map.contains_key("e"), hash_map.contains_key("e"));
    }

    #[test]
    fn get() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get(&1), hash_map.get(&1));
        assert_eq!(tuple_map.get(&2), hash_map.get(&2));
        assert_eq!(tuple_map.get(&3), hash_map.get(&3));
        assert_eq!(tuple_map.get(&4), hash_map.get(&4));
        assert_eq!(tuple_map.get(&5), hash_map.get(&5));

        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get("a"), hash_map.get("a"));
        assert_eq!(tuple_map.get("b"), hash_map.get("b"));
        assert_eq!(tuple_map.get("c"), hash_map.get("c"));
        assert_eq!(tuple_map.get("d"), hash_map.get("d"));
        assert_eq!(tuple_map.get("e"), hash_map.get("e"));
    }

    #[test]
    fn get_key_value() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_key_value(&1), hash_map.get_key_value(&1));
        assert_eq!(tuple_map.get_key_value(&2), hash_map.get_key_value(&2));
        assert_eq!(tuple_map.get_key_value(&3), hash_map.get_key_value(&3));
        assert_eq!(tuple_map.get_key_value(&4), hash_map.get_key_value(&4));
        assert_eq!(tuple_map.get_key_value(&5), hash_map.get_key_value(&5));

        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_key_value("a"), hash_map.get_key_value("a"));
        assert_eq!(tuple_map.get_key_value("b"), hash_map.get_key_value("b"));
        assert_eq!(tuple_map.get_key_value("c"), hash_map.get_key_value("c"));
        assert_eq!(tuple_map.get_key_value("d"), hash_map.get_key_value("d"));
        assert_eq!(tuple_map.get_key_value("e"), hash_map.get_key_value("e"));
    }

    #[test]
    fn get_full() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_full(&1), hash_map.get_full(&1));
        assert_eq!(tuple_map.get_full(&2), hash_map.get_full(&2));
        assert_eq!(tuple_map.get_full(&3), hash_map.get_full(&3));
        assert_eq!(tuple_map.get_full(&4), hash_map.get_full(&4));
        assert_eq!(tuple_map.get_full(&5), hash_map.get_full(&5));

        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_full("a"), hash_map.get_full("a"));
        assert_eq!(tuple_map.get_full("b"), hash_map.get_full("b"));
        assert_eq!(tuple_map.get_full("c"), hash_map.get_full("c"));
        assert_eq!(tuple_map.get_full("d"), hash_map.get_full("d"));
        assert_eq!(tuple_map.get_full("e"), hash_map.get_full("e"));
    }

    #[test]
    fn get_index_of() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_index_of(&1), hash_map.get_index_of(&1));
        assert_eq!(tuple_map.get_index_of(&2), hash_map.get_index_of(&2));
        assert_eq!(tuple_map.get_index_of(&3), hash_map.get_index_of(&3));
        assert_eq!(tuple_map.get_index_of(&4), hash_map.get_index_of(&4));
        assert_eq!(tuple_map.get_index_of(&5), hash_map.get_index_of(&5));

        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_index_of("a"), hash_map.get_index_of("a"));
        assert_eq!(tuple_map.get_index_of("b"), hash_map.get_index_of("b"));
        assert_eq!(tuple_map.get_index_of("c"), hash_map.get_index_of("c"));
        assert_eq!(tuple_map.get_index_of("d"), hash_map.get_index_of("d"));
        assert_eq!(tuple_map.get_index_of("e"), hash_map.get_index_of("e"));
    }

    #[test]
    fn get_mut() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_mut(&1), hash_map.get_mut(&1));
        assert_eq!(tuple_map.get_mut(&2), hash_map.get_mut(&2));
        assert_eq!(tuple_map.get_mut(&3), hash_map.get_mut(&3));
        assert_eq!(tuple_map.get_mut(&4), hash_map.get_mut(&4));
        assert_eq!(tuple_map.get_mut(&5), hash_map.get_mut(&5));
        *tuple_map.get_mut(&1).unwrap() = 5;
        *hash_map.get_mut(&1).unwrap() = 5;
        assert_eq!(tuple_map.get(&1), hash_map.get(&1));
        assert_eq!(&tuple_map.inner, &hash_map);

        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_mut("a"), hash_map.get_mut("a"));
        assert_eq!(tuple_map.get_mut("b"), hash_map.get_mut("b"));
        assert_eq!(tuple_map.get_mut("c"), hash_map.get_mut("c"));
        assert_eq!(tuple_map.get_mut("d"), hash_map.get_mut("d"));
        assert_eq!(tuple_map.get_mut("e"), hash_map.get_mut("e"));
        *tuple_map.get_mut("b").unwrap() = 5;
        *hash_map.get_mut("b").unwrap() = 5;
        assert_eq!(tuple_map.get("b"), hash_map.get("b"));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn get_full_mut() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_full_mut(&1), hash_map.get_full_mut(&1));
        assert_eq!(tuple_map.get_full_mut(&2), hash_map.get_full_mut(&2));
        assert_eq!(tuple_map.get_full_mut(&3), hash_map.get_full_mut(&3));
        assert_eq!(tuple_map.get_full_mut(&4), hash_map.get_full_mut(&4));
        assert_eq!(tuple_map.get_full_mut(&5), hash_map.get_full_mut(&5));
        *tuple_map.get_full_mut(&1).unwrap().2 = 5;
        *hash_map.get_full_mut(&1).unwrap().2 = 5;
        assert_eq!(tuple_map.get(&1), hash_map.get(&1));
        assert_eq!(&tuple_map.inner, &hash_map);

        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_full_mut("a"), hash_map.get_full_mut("a"));
        assert_eq!(tuple_map.get_full_mut("b"), hash_map.get_full_mut("b"));
        assert_eq!(tuple_map.get_full_mut("c"), hash_map.get_full_mut("c"));
        assert_eq!(tuple_map.get_full_mut("d"), hash_map.get_full_mut("d"));
        assert_eq!(tuple_map.get_full_mut("e"), hash_map.get_full_mut("e"));
        *tuple_map.get_full_mut("b").unwrap().2 = 5;
        *hash_map.get_full_mut("b").unwrap().2 = 5;
        assert_eq!(tuple_map.get("b"), hash_map.get("b"));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn remove_entry() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.remove(&1), hash_map.remove(&1));
        assert_eq!(tuple_map.remove(&3), hash_map.remove(&3));
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn swap_remove() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.swap_remove(&1), hash_map.swap_remove(&1));
        assert_eq!(tuple_map.swap_remove(&3), hash_map.swap_remove(&3));
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn swap_remove_entry() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(
            tuple_map.swap_remove_entry(&1),
            hash_map.swap_remove_entry(&1)
        );
        assert_eq!(
            tuple_map.swap_remove_entry(&3),
            hash_map.swap_remove_entry(&3)
        );
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn swap_remove_full() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(
            tuple_map.swap_remove_full(&1),
            hash_map.swap_remove_full(&1)
        );
        assert_eq!(
            tuple_map.swap_remove_full(&3),
            hash_map.swap_remove_full(&3)
        );
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn shift_remove() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.shift_remove(&1), hash_map.shift_remove(&1));
        assert_eq!(tuple_map.shift_remove(&3), hash_map.shift_remove(&3));
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn shift_remove_entry() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(
            tuple_map.shift_remove_entry(&1),
            hash_map.shift_remove_entry(&1)
        );
        assert_eq!(
            tuple_map.shift_remove_entry(&3),
            hash_map.shift_remove_entry(&3)
        );
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn shift_remove_full() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.insert(10, 20), hash_map.insert(10, 20));
        assert_eq!(tuple_map.insert(1, 3), hash_map.insert(1, 3));
        assert_eq!(tuple_map.insert(2, 2), hash_map.insert(2, 2));
        assert_eq!(tuple_map.get(&10), Some(&20));
        assert_eq!(tuple_map.get(&1), Some(&3));
        assert_eq!(tuple_map.get(&2), Some(&2));
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(
            tuple_map.shift_remove_full(&1),
            hash_map.shift_remove_full(&1)
        );
        assert_eq!(
            tuple_map.shift_remove_full(&3),
            hash_map.shift_remove_full(&3)
        );
        assert!(!tuple_map.contains_key(&1));
        assert!(!tuple_map.contains_key(&3));
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn pop() {
        method_test! { pop, }
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.pop(), hash_map.pop());
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.pop(), hash_map.pop());
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.pop(), hash_map.pop());
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.pop(), hash_map.pop());
        assert_eq!(&tuple_map.inner, &hash_map);
        assert!(tuple_map.is_empty());
    }

    #[test]
    fn retain() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        tuple_map.retain(|&k, &mut v| (k + v) % 2 == 0);
        hash_map.retain(|&k, &mut v| (k + v) % 2 == 0);
        assert!(tuple_map.is_empty());
        assert_eq!(&tuple_map.inner, &hash_map);

        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        tuple_map.retain(|&k, &mut v| (k + v) % 2 == 1);
        hash_map.retain(|&k, &mut v| (k + v) % 2 == 1);
        assert_eq!(tuple_map.len(), 3);
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    // tests sort_by, sorted_by, sort_unstable_by, sorted_unstable_by
    #[test]
    fn sorting() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();

        tuple_map.sort_by(|_, v1, _, v2| v1.cmp(v2).reverse());
        hash_map.sort_by(|_, v1, _, v2| v1.cmp(v2).reverse());
        assert_eq!(tuple_map.inner, hash_map);

        tuple_map.sort_unstable_by(|_, v1, _, v2| v1.cmp(v2));
        hash_map.sort_unstable_by(|_, v1, _, v2| v1.cmp(v2));
        assert_eq!(tuple_map.inner, hash_map);

        assert_eq!(
            tuple_map
                .clone()
                .sorted_by(|_, v1, _, v2| v1.cmp(v2).reverse())
                .collect::<Vec<_>>(),
            hash_map
                .clone()
                .sorted_by(|_, v1, _, v2| v1.cmp(v2).reverse())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            tuple_map
                .sorted_unstable_by(|_, v1, _, v2| v1.cmp(v2).reverse())
                .collect::<Vec<_>>(),
            hash_map
                .sorted_unstable_by(|_, v1, _, v2| v1.cmp(v2).reverse())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn sort_keys() {
        let mut tuple_map = [(3, 2), (2, 4), (1, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(3, 2), (2, 4), (1, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();

        assert_eq!(&tuple_map.inner, &hash_map);
        tuple_map.sort_keys();
        hash_map.sort_keys();
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn sort_unstable_keys() {
        let mut tuple_map = [(3, 2), (2, 4), (1, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(3, 2), (2, 4), (1, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();

        assert_eq!(&tuple_map.inner, &hash_map);
        tuple_map.sort_unstable_keys();
        hash_map.sort_unstable_keys();
        assert_eq!(&tuple_map.inner, &hash_map);
    }

    #[test]
    fn reverse() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexTupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<IndexMap<_, _>>();
        tuple_map.reverse();
        hash_map.reverse();
        assert_eq!(tuple_map.inner, hash_map);
    }
}
