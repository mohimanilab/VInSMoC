use pyo3::prelude::*;
use pyo3::types::PyString;
use pyo3::exceptions::PyValueError;

use crate::bond_type::BondType;

pub const SINGLE: &'static str = "SINGLE";
pub const DOUBLE: &'static str = "DOUBLE";
pub const TRIPLE: &'static str = "TRIPLE";
pub const AROMATIC: &'static str = "AROMATIC";

#[pyclass]
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Bond {
    #[pyo3(get, set)]
    pub start_atom: usize,
    #[pyo3(get, set)]
    pub end_atom: usize,
    #[pyo3(get, set)]
    pub bond_type: BondType,
}

#[pymethods]
impl Bond {
    #[new]
    pub fn new(start_atom: usize, end_atom: usize, bond_type: BondType) -> Self {
        Bond {
            start_atom,
            end_atom,
            bond_type,
        }
    }
}

impl IntoPy<PyObject> for BondType {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            BondType::Single => SINGLE.into_py(py),
            BondType::Double => DOUBLE.into_py(py),
            BondType::Triple => TRIPLE.into_py(py),
            BondType::Aromatic => AROMATIC.into_py(py),
        }
    }
}

impl<'source> FromPyObject<'source> for BondType {
    fn extract(obj: &'source PyAny) -> PyResult<Self> {
        let str_val = <PyString as PyTryFrom>::try_from(obj)?.to_str()?;
        match str_val {
            SINGLE => Ok(BondType::Single),
            DOUBLE => Ok(BondType::Double),
            TRIPLE => Ok(BondType::Triple),
            AROMATIC => Ok(BondType::Aromatic),
            _ => Err(PyValueError::new_err("Invalid bond type string"))
        }
    }
}

impl std::fmt::Display for Bond {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Bond({} {} {:?})",
            self.start_atom, self.end_atom, self.bond_type
        )
    }
}
