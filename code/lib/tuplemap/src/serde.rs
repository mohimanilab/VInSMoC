//! [`serde`] serialization and deserialization for maps that represents them as sequences of
//! tuples.
//!
//! These (de)serialization routines are intended for use with formats like JSON that require keys
//! of map types to be non-composite types. Instead of using some work-around intermediate type
//! like a [`String`], maps are serialized as sequences of key-value tuples. This module can be
//! used as part of the derive macros of [`serde`].
//!
//! # Example
//! Below we implement a custom counter type, which wraps a hash map from a generic type `T` to a
//! [`usize`] counter. We would like to be able use JSON (de)serialization with this type, so below
//! we show an example of the necessarily pieces of the [`serde`] derive macros to serialize this
//! type as a sequence of tuples. Note that we need implement [`FromIterator`] for this type to be
//! able to deserialize it.
//!
//! ```
//! use std::{collections::HashMap, hash::Hash, cmp::{PartialEq, Eq}};
//!
//! #[derive(Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
//! #[serde(bound(serialize = "T: serde::Serialize"))]
//! #[serde(bound(deserialize = "T: serde::Deserialize<'de> + Hash + Eq"))]
//! struct Counter<T: Eq + Hash> {
//!     #[serde(with = "tuplemap::serde")]
//!     inner: HashMap<T, usize>
//! }
//!
//! impl<T: Eq + Hash> Counter<T> {
//!     fn increment(&mut self, key: T) {
//!         *self.inner.entry(key).or_default() += 1;
//!     }
//! }
//!
//! impl<T: Eq + Hash> FromIterator<(T, usize)> for Counter<T> {
//!     fn from_iter<I: IntoIterator<Item=(T, usize)>>(iter: I) -> Self {
//!         Self { inner: iter.into_iter().collect() }
//!     }
//! }
//!
//! #[derive(Debug, Hash, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
//! struct ComplexStruct {
//!     a: u32,
//!     b: u32,
//!     c: u32,
//! }
//!
//! // we can serialize to/from JSON with a complex key
//! let mut counter = Counter::<ComplexStruct> { inner: HashMap::new() };
//! counter.increment(ComplexStruct { a: 1, b: 1, c: 1 });
//! counter.increment(ComplexStruct { a: 1, b: 1, c: 1 });
//! counter.increment(ComplexStruct { a: 2, b: 1, c: 1 });
//! let json_string = serde_json::to_string(&counter).unwrap();
//! let rebuilt = serde_json::from_str(&json_string).unwrap();
//! assert_eq!(counter, rebuilt);
//!
//! // we can't if we just used a hashmap
//! let mut simple_counter = HashMap::new();
//! simple_counter.insert(ComplexStruct { a: 1, b: 1, c: 1 }, 2);
//! simple_counter.insert(ComplexStruct { a: 2, b: 1, c: 1 }, 1);
//! assert_eq!(&counter.inner, &simple_counter);
//! assert!(serde_json::to_string(&simple_counter).is_err());
//! ```

use serde::{
    de::Deserializer,
    ser::{SerializeSeq, Serializer},
    Deserialize, Serialize,
};

/// Serialization routine for maps.
///
/// In formats like JSON keys must be serialized into a trivial (as opposed to composite) type.
/// What this ends up meaning for practical use is that a map structure whose keys are structs
/// needs to undergo some sort of intermediate step if they are to be used with JSON formats.
/// This routine serializes a map type as a sequence of key-value tuples, making it suitable for
/// use with formats like JSON.
///
/// The input generic `T` accepts maps as any type that can be converted into an owned exactly
/// sized iterator over tuples whose two elements implement [`serde::Serialize`]. A notable type
/// that will work with this routine that isn't a map is a vector of tuples.
///
/// This should be used together with [`deserialize`].
///
/// # Errors
/// Propagates [`serde`] errors from calls to [`serde::Serializer::serialize_seq`] and
/// [`serde::ser::SerializeSeq::serialize_element`].
pub fn serialize<'a, S, T, K, V>(map: &'a T, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    &'a T: IntoIterator<Item = (&'a K, &'a V)>,
    <&'a T as IntoIterator>::IntoIter: ExactSizeIterator,
    K: Serialize + 'a,
    V: Serialize + 'a,
{
    let it = map.into_iter();
    let mut seq = serializer.serialize_seq(Some(it.len()))?;
    for key_val_tuple in it {
        seq.serialize_element(&key_val_tuple)?;
    }
    seq.end()
}

/// Deserialization routine for maps.
///
/// In formats like JSON keys must be serialized into a trivial (as opposite to composite) type.
/// What this ends up meaning for practical use is that a map structure whose keys are structs
/// needs to undergo some sort of intermediate stepp if they are to be used with JSON formats. This
/// routine deserializes map types from a sequence of key-value tuples, making it suitable for use
/// with formats with JSON.
///
/// Thbei nput generic `T` only needs to be constructable from an iterator of owned key-value
/// tuples, where these tuples implement [`serde::Deserialize`]. A notable type that will work with
/// this routine that isn't a map is a vector of tuples.
///
/// This routine allocates an intermediate vector of tuples which is then converted into the target
/// type `T` via [`FromIterator`].
///
/// This should be used together with [`serialize`].
///
/// # Errors
/// Propagates [`serde`] errors from deserialization of a vector of tuples.
pub fn deserialize<'de, D, T, K, V>(deserializer: D) -> Result<T, D::Error>
where
    D: Deserializer<'de>,
    T: FromIterator<(K, V)>,
    (K, V): Deserialize<'de>,
{
    let tup_vec = Vec::<(K, V)>::deserialize(deserializer)?;
    Ok(tup_vec.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use crate::{IndexTupleMap, TupleMap};

    #[test]
    fn tuplemap_ser_de() {
        let mut map = TupleMap::<u64, Vec<u32>>::default();
        map.insert(10, vec![1, 2, 3]);
        map.insert(20, vec![6, 7, 8]);
        map.insert(30, vec![29]);
        map.insert(12, vec![900, 20]);
        map.insert(0, vec![20]);

        assert_eq!(
            serde_json::from_str::<TupleMap<u64, Vec<u32>>>(&serde_json::to_string(&map).unwrap())
                .unwrap(),
            map
        );
    }

    #[test]
    fn indextuplemap_ser_de() {
        let mut map = IndexTupleMap::<u64, Vec<u32>>::default();
        map.insert(10, vec![1, 2, 3]);
        map.insert(20, vec![6, 7, 8]);
        map.insert(30, vec![29]);
        map.insert(12, vec![900, 20]);
        map.insert(0, vec![20]);

        assert_eq!(
            serde_json::from_str::<IndexTupleMap<u64, Vec<u32>>>(
                &serde_json::to_string(&map).unwrap()
            )
            .unwrap(),
            map
        );
    }
}
