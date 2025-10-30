use crate::atoms::Atom;
use crate::bond::Bond;
use crate::bond_type::BondType;

use pyo3::prelude::*;
use rustc_hash::FxHashMap;

#[pyclass]
#[derive(Debug, Eq, PartialEq)]
pub struct Mol {
    #[pyo3(get, set)]
    pub graph: UnGraph<Atom, BondType>,
    #[pyo3(get, set)]
    pub hydrogens: FxHashMap<NodeIndex<u32>, u8>,
    #[pyo3(get, set)]
    pub charges: FxHashMap<NodeIndex<u32>, i8>,
}

#[pymethods]
impl Mol {
    pub fn total_bonds(&self, atom_idx: NodeIndex<u32>) -> u8 {
        self.graph
            .edges(atom_idx)
            .map(|x| match x.weight() {
                BondType::Single => 1,
                BondType::Double => 2,
                BondType::Triple => 3,
                BondType::Aromatic => 1,
            })
            .sum()
    }
}
