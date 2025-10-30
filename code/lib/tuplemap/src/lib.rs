//! Map types serialized as vectors of tuples.
//!
//! Using [`std::collections::HashMap`] with [`serde`] has a pitfall. If we're considering JSON,
//! the data model the keys in a map must be simple (as opposed to composite) types. If you have a
//! hashable struct as the key of a map, you cannot directly (de)serialize it to JSON. One
//! workaround is to manually implement the relevant [`serde`] traits and use some simple (often
//! string) representation of the complex key type. This can be fine, but it can get tedious. This
//! crate provides map types that serialize as vectors of tuples, making them suitable drop-in
//! replacements that work with JSON serialization.
//!
//! This crate provides a few map types that all serve as drop-in replacements to their
//! counterparts.
//!
//! | `tuplemap` type | Corresponding type |
//! | --------------- | ------------------ |
//! | [`TupleMap`]    | [`std::collections::HashMap`] |
//! | [`FxTupleMap`]  | [`rustc_hash::FxHashMap`] |
//! | [`IndexTupleMap`] | [`indexmap::IndexMap`] |
//! | [`FxIndexTupleMap`] | [`indexmap::IndexMap`] with [`rustc_hash`] hasher |
//!
//! # Examples
//! ```
//! // a complex key type we'll be using as an example
//! #[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Deserialize, serde::Serialize)]
//! struct ComplexKey {
//!     a: i32,
//!     b: String,
//!     c: bool,
//! }
//!
//! // if we build a TupleMap, we can serialize to JSON
//! use tuplemap::TupleMap;
//!
//! let mut tuple_map = TupleMap::<ComplexKey, i32>::new();
//! tuple_map.insert(ComplexKey { a: 1, b: "hello".into(), c: false }, 32);
//! tuple_map.insert(ComplexKey { a: 2, b: "goodbye".into(), c: true }, 66);
//! let json_str = serde_json::to_string(&tuple_map).unwrap();
//! let rebuilt = serde_json::from_str(&json_str).unwrap();
//! assert_eq!(&tuple_map, &rebuilt);
//!
//! // if we do the same thing for HashMap, it fails
//! use std::collections::HashMap;
//!
//! let mut hash_map = HashMap::<ComplexKey, i32>::new();
//! hash_map.insert(ComplexKey { a: 1, b: "hello".into(), c: false }, 32);
//! hash_map.insert(ComplexKey { a: 2, b: "goodbye".into(), c: true }, 66);
//! assert!(serde_json::to_string(&hash_map).is_err());
//! ```

#![warn(
    missing_docs,
    clippy::all,
    clippy::pedantic,
    rust_2018_idioms,
    future_incompatible
)]

pub mod index_tuple_map;
pub mod serde;
pub mod tuple_map;

pub use index_tuple_map::{FxIndexTupleMap, IndexTupleMap};
pub use tuple_map::{FxTupleMap, TupleMap};
