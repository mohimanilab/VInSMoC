use pyo3::prelude::*;
use pyo3::types::PyList;

use rustc_hash::FxHashMap;

use crate::{BondType, Bond};
use crate::atoms;

#[pyclass]
#[derive(Debug, Eq, PartialEq)]
pub struct Mol {
    #[pyo3(get, set)]
    pub atoms: Vec<atoms::Atom>,
    #[pyo3(get, set)]
    pub bonds: Vec<Bond>,
    #[pyo3(get, set)]
    pub hydrogens: Vec<u8>,
    #[pyo3(get, set)]
    pub charges: Vec<i8>,
    #[pyo3(get, set)]
    pub connectivity: FxHashMap<usize, FxHashMap<usize, usize>>,
}

#[pymethods]
impl Mol {
    pub fn total_bonds(&self, atom_idx: usize) -> u8 {
        self.bonds
            .iter()
            .filter_map(|x| {
                if x.start_atom == atom_idx || x.end_atom == atom_idx {
                    Some(match x.bond_type {
                        BondType::Single => 1,
                        BondType::Double => 2,
                        BondType::Triple => 3,
                        BondType::Aromatic => 1,
                    })
                } else {
                    None
                }
            })
            .sum()
    }
}
