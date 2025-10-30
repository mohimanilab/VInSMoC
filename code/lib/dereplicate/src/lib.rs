#![warn(clippy::all, clippy::pedantic, rust_2018_idioms, future_incompatible)]
#![allow(
    clippy::must_use_candidate,
    clippy::module_name_repetitions,
    clippy::missing_panics_doc,
    clippy::shadow_unrelated,
    clippy::unseparated_literal_suffix
)]
#![cfg_attr(
    test,
    allow(clippy::pedantic, clippy::similar_names, clippy::float_cmp,)
)]

mod errors;
pub use errors::Error;
mod matches;
pub use matches::{DereplicationRun, Msm, MsmPValue, MsmSmiles, MsmStub};

pub mod db;
