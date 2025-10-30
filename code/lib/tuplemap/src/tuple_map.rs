//! Drop-in replacements for [`std::collections::HashMap`] and [`rustc_hash::FxHashMap`] that
//! serialize as vectors of tuples.
//!
//! This exposes [`TupleMap`] and [`FxTupleMap`] that are suitable for JSON serialization with no
//! extra steps, even when key types are non-trivial. Both types re-implement all methods for their
//! corresponding map types. [`TupleMap`] re-implements all traits implemented on
//! [`std::collections::HashMap`].

pub use std::collections::hash_map::{
    Drain, Entry, IntoIter, IntoKeys, IntoValues, Iter, IterMut, Keys, Values, ValuesMut,
};
use std::{
    borrow::Borrow,
    collections::{hash_map::RandomState, HashMap, TryReserveError},
    fmt,
    hash::{BuildHasher, BuildHasherDefault, Hash},
    ops::Index,
    panic::UnwindSafe,
};

use rustc_hash::FxHasher;

/// A drop-in replacement for [`std::collections::HashMap`] that serializes as a vector of tuples
/// with [`serde`].
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(bound(serialize = "K: serde::Serialize, V: serde::Serialize"))]
#[serde(bound(
    deserialize = "K: serde::Deserialize<'de> + Hash + Eq, V: serde::Deserialize<'de>, S: BuildHasher + Default"
))]
#[serde(transparent)]
pub struct TupleMap<K, V, S = RandomState> {
    #[serde(with = "crate::serde")]
    inner: HashMap<K, V, S>,
}

/// A drop-in replacment for [`rustc_hash::FxHashMap`] that serializes as a vector of tuples with
/// [`serde`].
#[allow(clippy::module_name_repetitions)]
pub type FxTupleMap<K, V> = TupleMap<K, V, BuildHasherDefault<FxHasher>>;

macro_rules! inherit {
    (ctor $name:ident, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        #[must_use]
        pub fn $name($($v:$t),*) -> Self {
            Self {
                inner: HashMap::$name($($v),*)
            }
        }
    };
    ($name:ident, void, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name(&self, $($v:$t),*) {
            self.inner.$name($($v),*)
        }
    };
    (mut $name:ident, void, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name(&mut self, $($v:$t),*) {
            self.inner.$name($($v),*)
        }
    };
    (owned $name:ident, void, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name(self, $($v:$t),*) {
            self.inner.$name($($v),*)
        }
    };
    ($name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name(&self, $($v:$t),*) -> $ret {
            self.inner.$name($($v),*)
        }
    };
    (mut $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name(&mut self, $($v:$t),*) -> $ret {
            self.inner.$name($($v),*)
        }
    };
    (owned $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name(self, $($v:$t),*) -> $ret {
            self.inner.$name($($v),*)
        }
    };
    (get $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name<Q>(&self, $($v:$t),*) -> $ret
        where
            K: Borrow<Q>,
            Q: Hash + Eq + ?Sized
        {
            self.inner.$name($($v),*)
        }
    };
    (mut get $name:ident, $ret:ty, $($v:ident: $t:ty),*) => {
        #[doc = concat!("See [`HashMap::", stringify!($name), "`].")]
        pub fn $name<Q>(&mut self, $($v:$t),*) -> $ret
        where
            K: Borrow<Q>,
            Q: Hash + Eq + ?Sized
        {
            self.inner.$name($($v),*)
        }
    };
}

impl<K, V> TupleMap<K, V, RandomState> {
    inherit! { ctor new, }
    inherit! { ctor with_capacity, capacity: usize }
}

impl<K, V, S> TupleMap<K, V, S> {
    inherit! { ctor with_hasher, hash_builder: S }
    inherit! { ctor with_capacity_and_hasher, capacity: usize, hasher: S }
    inherit! { capacity, usize, }
    inherit! { keys, Keys<'_, K, V>, }
    inherit! { values, Values<'_, K, V>, }
    inherit! { mut values_mut, ValuesMut<'_, K, V>, }
    inherit! { iter, Iter<'_,K,V>, }
    inherit! { mut iter_mut, IterMut<'_, K, V>, }
    inherit! { len, usize, }
    inherit! { is_empty, bool, }
    inherit! { mut drain, Drain<'_, K, V>, }

    inherit! { mut clear, void, }
    inherit! { hasher, &S, }
}

impl<K, V, S> TupleMap<K, V, S>
where
    K: Eq + Hash,
    S: BuildHasher,
{
    inherit! { mut reserve, void, additional: usize }
    inherit! { mut try_reserve, Result<(), TryReserveError>, additional: usize }
    inherit! { mut shrink_to_fit, void, }
    inherit! { mut shrink_to, void, min_capacity: usize }
    inherit! { mut entry, Entry<'_, K, V>, key: K }
    inherit! { get get, Option<&V>, k: &Q}
    inherit! { get get_key_value, Option<(&K, &V)>, k: &Q}
    inherit! { get contains_key, bool, k: &Q }
    inherit! { mut get get_mut, Option<&mut V>, k: &Q }
    inherit! { mut insert, Option<V>, k: K, v: V }
    inherit! { mut get remove, Option<V>, k: &Q }
    inherit! { mut get remove_entry, Option<(K, V)>, k: &Q }

    // in 1.59.0 the below methods moved into the other impl block, but we're using 1.58.1

    // special case for retain, no need to deal with generics in macro for one method

    /// See [`HashMap::retain`].
    pub fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&K, &mut V) -> bool,
    {
        self.inner.retain(f);
    }

    inherit! { owned into_keys, IntoKeys<K, V>, }
    inherit! { owned into_values, IntoValues<K, V>, }
}

// trait implementations from std lib redone below

impl<K, V, S> fmt::Debug for TupleMap<K, V, S>
where
    K: fmt::Debug,
    V: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.inner.fmt(f)
    }
}

impl<K, V, S: Default> Default for TupleMap<K, V, S> {
    fn default() -> Self {
        Self {
            inner: HashMap::default(),
        }
    }
}

impl<'a, K, V, S> Extend<(&'a K, &'a V)> for TupleMap<K, V, S>
where
    K: Eq + Hash + Copy + 'a,
    V: Copy + 'a,
    S: BuildHasher,
{
    fn extend<T: IntoIterator<Item = (&'a K, &'a V)>>(&mut self, iter: T) {
        self.inner.extend(iter);
    }
}

impl<K, V, S> Extend<(K, V)> for TupleMap<K, V, S>
where
    K: Eq + Hash,
    S: BuildHasher,
{
    fn extend<T: IntoIterator<Item = (K, V)>>(&mut self, iter: T) {
        self.inner.extend(iter);
    }
}

impl<K, V, const N: usize> From<[(K, V); N]> for TupleMap<K, V, RandomState>
where
    K: Eq + Hash,
{
    fn from(arr: [(K, V); N]) -> Self {
        Self { inner: arr.into() }
    }
}

impl<K, V, S> FromIterator<(K, V)> for TupleMap<K, V, S>
where
    K: Eq + Hash,
    S: BuildHasher + Default,
{
    fn from_iter<T: IntoIterator<Item = (K, V)>>(iter: T) -> Self {
        Self {
            inner: iter.into_iter().collect(),
        }
    }
}

impl<K, Q, V, S> Index<&Q> for TupleMap<K, V, S>
where
    K: Eq + Hash + Borrow<Q>,
    Q: Eq + Hash + ?Sized,
    S: BuildHasher,
{
    type Output = V;

    fn index(&self, index: &Q) -> &Self::Output {
        self.inner.index(index)
    }
}

impl<'a, K, V, S> IntoIterator for &'a TupleMap<K, V, S> {
    type Item = (&'a K, &'a V);

    type IntoIter = Iter<'a, K, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter()
    }
}

impl<'a, K, V, S> IntoIterator for &'a mut TupleMap<K, V, S> {
    type Item = (&'a K, &'a mut V);

    type IntoIter = IterMut<'a, K, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter_mut()
    }
}

impl<K, V, S> IntoIterator for TupleMap<K, V, S> {
    type Item = (K, V);

    type IntoIter = IntoIter<K, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter()
    }
}

impl<K, V, S> PartialEq<TupleMap<K, V, S>> for TupleMap<K, V, S>
where
    K: Eq + Hash,
    V: PartialEq,
    S: BuildHasher,
{
    fn eq(&self, other: &TupleMap<K, V, S>) -> bool {
        self.inner == other.inner
    }
}

impl<K, V, S> Eq for TupleMap<K, V, S>
where
    K: Eq + Hash,
    V: Eq,
    S: BuildHasher,
{
}

impl<K, V, S> UnwindSafe for TupleMap<K, V, S>
where
    K: UnwindSafe,
    V: UnwindSafe,
    S: UnwindSafe,
{
}

impl<K, V, S> From<HashMap<K, V, S>> for TupleMap<K, V, S> {
    fn from(inner: HashMap<K, V, S>) -> Self {
        Self { inner }
    }
}

impl<K, V, S> From<TupleMap<K, V, S>> for HashMap<K, V, S> {
    fn from(outer: TupleMap<K, V, S>) -> Self {
        outer.inner
    }
}

impl<K, V, S> AsRef<HashMap<K, V, S>> for TupleMap<K, V, S> {
    fn as_ref(&self) -> &HashMap<K, V, S> {
        &self.inner
    }
}

#[cfg(test)]
mod tests {
    // for the macros
    #![allow(unused_mut)]

    use super::*;

    use rustc_hash::FxHashMap;

    macro_rules! method_test {
        (ctor $name:ident, $($params:expr),*) => {
            let tuple_map = TupleMap::<u64, u64>::$name($($params),*);
            let map = HashMap::<u64, u64>::$name($($params),*);
            assert_eq!(tuple_map.inner, map);
        };
        (other ctor $name:ident, $($params:expr),*) => {
            let tuple_map = TupleMap::<u64, u64, _>::$name($($params),*);
            let map = HashMap::<u64, u64, _>::$name($($params),*);
            assert_eq!(tuple_map.inner, map);
        };
        ($name:ident, $($params:expr),*) => {
            let mut tuple_map = TupleMap::<u64, u64>::default();
            let mut map = HashMap::<u64, u64>::default();
            assert_eq!(tuple_map.$name($($params),*), map.$name($($params),*));
            assert_eq!(tuple_map.inner, map);

            let mut tuple_map = TupleMap::<u64, u64>::default();
            tuple_map.insert(1,1);
            let mut map = HashMap::<u64, u64>::default();
            map.insert(1,1);
            assert_eq!(tuple_map.$name($($params),*), map.$name($($params),*));
            assert_eq!(tuple_map.inner, map);


            let mut tuple_map = TupleMap::<u64, u64>::default();
            tuple_map.insert(1,1);
            tuple_map.insert(2,3);
            tuple_map.insert(4,5);
            tuple_map.insert(6,7);
            let mut map = HashMap::<u64, u64>::default();
            map.insert(1,1);
            map.insert(2,3);
            map.insert(4,5);
            map.insert(6,7);
            assert_eq!(tuple_map.$name($($params),*), map.$name($($params),*));
            assert_eq!(tuple_map.inner, map);

            let mut tuple_map = FxTupleMap::<String, i8>::default();
            tuple_map.insert("hello".into(), 5);
            tuple_map.insert("bye".into(), 5);
            tuple_map.insert("hello".into(), 3);
            let mut map = FxHashMap::<String, i8>::default();
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
        let tuple_map = TupleMap::<i8, i8>::with_capacity(0);
        assert_eq!(tuple_map.capacity(), 0);
        let tuple_map = TupleMap::<i8, i8>::with_capacity(100);
        assert!(tuple_map.capacity() >= 100);
        let tuple_map = TupleMap::<i8, i8, _>::with_capacity_and_hasher(
            0,
            BuildHasherDefault::<FxHasher>::default(),
        );
        assert_eq!(tuple_map.capacity(), 0);
        let tuple_map = TupleMap::<i8, i8, _>::with_capacity_and_hasher(
            10,
            BuildHasherDefault::<FxHasher>::default(),
        );
        assert!(tuple_map.capacity() >= 10);
    }

    // tests keys, values, values_mut, iter, iter_mut, into_keys, into_values, into_iter
    #[test]
    fn iterators() {
        // keys
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let mut tuple_iter = tuple_map.keys().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.keys().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![&1, &3, &5]);
        assert_eq!(tuple_iter, hash_iter);

        // values
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let mut tuple_iter = tuple_map.values().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.values().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![&2, &4, &6]);
        assert_eq!(tuple_iter, hash_iter);

        // iter
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let mut tuple_iter = tuple_map.iter().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.iter().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![(&1, &2), (&3, &4), (&5, &6)]);
        assert_eq!(tuple_iter, hash_iter);

        // values_mut
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let mut tuple_iter = tuple_map.values_mut().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.values_mut().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![&mut 2, &mut 4, &mut 6]);
        assert_eq!(tuple_iter, hash_iter);

        // iter_mut
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let mut tuple_iter = tuple_map.iter_mut().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.iter_mut().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![(&1, &mut 2), (&3, &mut 4), (&5, &mut 6)]);
        assert_eq!(tuple_iter, hash_iter);

        // into_keys
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let mut tuple_iter = tuple_map.into_keys().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.into_keys().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![1, 3, 5]);
        assert_eq!(tuple_iter, hash_iter);

        // into_values
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let mut tuple_iter = tuple_map.into_values().collect::<Vec<_>>();
        tuple_iter.sort();
        let mut hash_iter = hash_map.into_values().collect::<Vec<_>>();
        hash_iter.sort();
        assert_eq!(tuple_iter, vec![2, 4, 6]);
        assert_eq!(tuple_iter, hash_iter);

        // into_iter
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
            .collect::<TupleMap<_, _>>();
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
            .collect::<TupleMap<_, _>>();
        assert!(!map.is_empty());
        map.remove(&1);
        assert!(!map.is_empty());
        map.remove(&3);
        assert!(!map.is_empty());
        map.remove(&5);
        assert!(map.is_empty());
    }

    #[test]
    fn drain() {
        let mut tuple_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<TupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        let mut tuple_drained = tuple_map.drain().collect::<Vec<_>>();
        tuple_drained.sort_by(|(k1, v1), (k2, v2)| k1.cmp(k2).then(v1.cmp(v2)));
        let mut hash_drained = hash_map.drain().collect::<Vec<_>>();
        hash_drained.sort_by(|(k1, v1), (k2, v2)| k1.cmp(k2).then(v1.cmp(v2)));
        assert_eq!(tuple_drained, hash_drained);
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.len(), 0);
    }

    #[test]
    fn clear() {
        method_test! { clear, }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        map.clear();
        assert_eq!(map, TupleMap::<_, _>::default());
    }

    // hasher doesn't impl eq always, skipping tests

    #[test]
    fn reserve() {
        method_test! { reserve, 10 }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        map.reserve(20);
        assert!(map.capacity() >= 23);
        map.reserve(10);
        assert!(map.capacity() >= 23);
        map.reserve(30);
        assert!(map.capacity() >= 33);
        let mut map = TupleMap::<i32, i32>::new();
        assert_eq!(map.capacity(), 0);
        map.reserve(10);
        assert!(map.capacity() >= 10);
    }

    #[test]
    fn try_reserve() {
        method_test! { try_reserve, 10 }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        map.try_reserve(20).unwrap();
        assert!(map.capacity() >= 23);
        map.try_reserve(10).unwrap();
        assert!(map.capacity() >= 23);
        map.try_reserve(30).unwrap();
        assert!(map.capacity() >= 33);
        let mut map = TupleMap::<i32, i32>::new();
        assert_eq!(map.capacity(), 0);
        map.try_reserve(10).unwrap();
        assert!(map.capacity() >= 10);
        assert!(map.try_reserve(usize::MAX).is_err());
    }

    #[test]
    fn shrink_to_fit() {
        method_test! { shrink_to_fit, }
        let mut map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        map.reserve(100);
        assert!(map.capacity() >= 103);
        map.shrink_to_fit();
        assert!(map.capacity() >= 3);
        let mut map = TupleMap::<u64, u64>::new();
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
            .collect::<TupleMap<_, _>>();
        map.reserve(100);
        assert!(map.capacity() >= 103);
        map.shrink_to(50);
        assert!(map.capacity() >= 50);
        map.shrink_to(0);
        assert!(map.capacity() >= 3);
        let mut map = TupleMap::<i32, i32>::new();
        map.insert(1, 1);
        map.reserve(100);
        assert!(map.capacity() >= 100);
        map.shrink_to(50);
        assert!(map.capacity() >= 50);
        map.shrink_to(0);
        assert!(map.capacity() >= 1);
    }

    #[test]
    fn entry() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
    fn get() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
        .collect::<TupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
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
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
        .collect::<TupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.get_key_value("a"), hash_map.get_key_value("a"));
        assert_eq!(tuple_map.get_key_value("b"), hash_map.get_key_value("b"));
        assert_eq!(tuple_map.get_key_value("c"), hash_map.get_key_value("c"));
        assert_eq!(tuple_map.get_key_value("d"), hash_map.get_key_value("d"));
        assert_eq!(tuple_map.get_key_value("e"), hash_map.get_key_value("e"));
    }

    #[test]
    fn contains_key() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
        .collect::<TupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        assert_eq!(&tuple_map.inner, &hash_map);
        assert_eq!(tuple_map.contains_key("a"), hash_map.contains_key("a"));
        assert_eq!(tuple_map.contains_key("b"), hash_map.contains_key("b"));
        assert_eq!(tuple_map.contains_key("c"), hash_map.contains_key("c"));
        assert_eq!(tuple_map.contains_key("d"), hash_map.contains_key("d"));
        assert_eq!(tuple_map.contains_key("e"), hash_map.contains_key("e"));
    }

    #[test]
    fn get_mut() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
        .collect::<TupleMap<_, _>>();
        let mut hash_map = [
            ("a".to_owned(), 1),
            ("b".to_owned(), 2),
            ("c".to_owned(), 3),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
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
    fn insert_remove() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
    fn remove_entry() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
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
    fn retain() {
        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        tuple_map.retain(|&k, &mut v| (k + v) % 2 == 0);
        hash_map.retain(|&k, &mut v| (k + v) % 2 == 0);
        assert!(tuple_map.is_empty());
        assert_eq!(&tuple_map.inner, &hash_map);

        let mut tuple_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<TupleMap<_, _>>();
        let mut hash_map = [(1, 2), (3, 4), (5, 6)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        tuple_map.retain(|&k, &mut v| (k + v) % 2 == 1);
        hash_map.retain(|&k, &mut v| (k + v) % 2 == 1);
        assert_eq!(tuple_map.len(), 3);
        assert_eq!(&tuple_map.inner, &hash_map);
    }
}
