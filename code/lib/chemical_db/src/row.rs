use std::{convert::TryFrom, hash};

use molecule::{BuildMolError, Mol, MolBuilder};
use serde::{Deserialize, Serialize};

/// Radius to use for [`Mol::ecfp_hash`] when computing [`DatabaseRow::inchikey14`].
const ECFP_RADIUS: u8 = 8;

/// A candidate molecule in a chemical database.
///
/// These entries cache the mass and InChIKey for each molecule, but only store the SMILES, not the
/// full molecular graph.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct DatabaseRow {
    /// The structure of this compound
    pub smiles: String,
    pub(crate) exact_mass: f64,
    pub(crate) name: String,
    #[serde(with = "crate::row::inchi_serde")]
    pub(crate) inchikey14: [u8; 14],
}

/// [`serde`] implementations for [u8; 14]
mod inchi_serde {
    use std::convert::TryInto;

    use serde::{de::Unexpected, Deserialize, Deserializer, Serializer};

    pub fn serialize<S>(bytes: &[u8; 14], ser: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::Error;
        let s = std::str::from_utf8(bytes).map_err(S::Error::custom)?;
        ser.serialize_str(s)
    }

    pub fn deserialize<'de, D>(de: D) -> Result<[u8; 14], D::Error>
    where
        D: Deserializer<'de>,
    {
        use serde::de::Error;
        let inchikey = <&str as Deserialize>::deserialize(de)?;
        let inchikey: [u8; 14] = inchikey.as_bytes().try_into().map_err(|_| {
            D::Error::invalid_value(Unexpected::Str(inchikey), &"a 14 character ASCII string")
        })?;
        Ok(inchikey)
    }
}

impl DatabaseRow {
    /// Construct a new row from its raw pieces.
    // FIXME: make the params in new the same order as struct decl
    #[allow(clippy::inconsistent_struct_constructor)]
    pub fn new(name: String, smiles: String, exact_mass: f64, inchikey14: [u8; 14]) -> Self {
        Self {
            name,
            smiles,
            exact_mass,
            inchikey14,
        }
    }

    /// Construct a [`DatabaseRow`] from a named molecule.
    ///
    /// # Errors
    /// Propagates errors from [`Mol::to_inchi_key`].
    pub fn from_named_mol(name: String, mol: &Mol) -> Self {
        Self {
            name,
            smiles: mol.to_smiles(),
            inchikey14: mol.ecfp_hash(ECFP_RADIUS),
            exact_mass: mol.exact_mass(),
        }
    }

    /// Construct a [`DatabaseRow`] from a molecule.
    ///
    /// The name of the compound in the resulting [`DatabaseRow`] will be the same as its InChIKey.
    /// This is identical to the [`std::convert::From`] implementation.
    pub fn from_mol(mol: &Mol) -> Self {
        let hash = mol.ecfp_hash(ECFP_RADIUS);
        debug_assert_eq!(hash.len(), 14);
        let name = String::from_utf8(hash.to_vec()).unwrap();

        Self {
            name,
            smiles: mol.to_smiles(),
            // we already checked the len, .unwrap() is safe
            inchikey14: hash,
            exact_mass: mol.exact_mass(),
        }
    }

    /// Retrieves the exact mass of this compound.
    pub fn exact_mass(&self) -> f64 {
        self.exact_mass
    }

    /// Retrieves the name of this compound.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Retrieves 14-character InChIKey prefix of this compound.
    pub fn inchikey14(&self) -> &str {
        std::str::from_utf8(&self.inchikey14).unwrap()
    }

    /// Retrieves the raw byte buffer of the InChIKey prefix of this compound.
    pub fn inchikey14_bytes(&self) -> &[u8; 14] {
        &self.inchikey14
    }

    /// Converts this compound into a [`Mol`](molecule::Mol).
    ///
    /// # Errors
    /// If the SMILES is invalid.
    pub fn mol(&self) -> Result<Mol, BuildMolError> {
        Mol::try_from(self)
    }
}

impl PartialEq for DatabaseRow {
    fn eq(&self, other: &Self) -> bool {
        self.inchikey14 == other.inchikey14
    }
}

impl Eq for DatabaseRow {}

impl hash::Hash for DatabaseRow {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        self.inchikey14.hash(state);
    }
}

impl TryFrom<&DatabaseRow> for Mol {
    type Error = BuildMolError;

    fn try_from(db_row: &DatabaseRow) -> Result<Self, Self::Error> {
        Ok(MolBuilder::from_smiles(&db_row.smiles)?
            .implicit_hs()
            .build())
    }
}

impl From<Mol> for DatabaseRow {
    /// See [`DatabaseRow::from_mol`] for details.
    fn from(value: Mol) -> Self {
        Self::from_mol(&value)
    }
}

impl From<&Mol> for DatabaseRow {
    /// See [`DatabaseRow::from_mol`] for details.
    fn from(value: &Mol) -> Self {
        Self::from_mol(value)
    }
}

impl From<(String, &Mol)> for DatabaseRow {
    /// See [`DatabaseRow::from_named_mol`] for details.
    fn from((name, mol): (String, &Mol)) -> Self {
        Self::from_named_mol(name, mol)
    }
}

impl From<(String, Mol)> for DatabaseRow {
    /// See [`DatabaseRow::from_named_mol`] for details.
    fn from((name, mol): (String, Mol)) -> Self {
        Self::from_named_mol(name, &mol)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_hash_eq() {
        fn calculate_hash<T: hash::Hash>(t: &T) -> u64 {
            use std::hash::Hasher;
            let mut s = std::collections::hash_map::DefaultHasher::new();
            t.hash(&mut s);
            s.finish()
        }

        assert_eq!(
            DatabaseRow {
                smiles: String::from("CCC"),
                exact_mass: 1.0,
                name: String::from("name"),
                inchikey14: *b"AAAAAAAAAAAAAA",
            },
            DatabaseRow {
                smiles: String::from("CCCCCC"),
                exact_mass: 10.0,
                name: String::from("othername"),
                inchikey14: *b"AAAAAAAAAAAAAA",
            },
        );
        assert_ne!(
            DatabaseRow {
                smiles: String::from("CCC"),
                exact_mass: 1.0,
                name: String::from("name"),
                inchikey14: *b"AAAAAAAAAAAAAA",
            },
            DatabaseRow {
                smiles: String::from("CCC"),
                exact_mass: 1.0,
                name: String::from("name"),
                inchikey14: *b"AAAAAAAABAAAAA",
            },
        );

        assert_eq!(
            calculate_hash(&DatabaseRow {
                smiles: String::from("CCC"),
                exact_mass: 1.0,
                name: String::from("name"),
                inchikey14: *b"AAAAAAAAAAAAAA",
            }),
            calculate_hash(&DatabaseRow {
                smiles: String::from("CCCCCC"),
                exact_mass: 10.0,
                name: String::from("othername"),
                inchikey14: *b"AAAAAAAAAAAAAA",
            }),
        );
        assert_ne!(
            calculate_hash(&DatabaseRow {
                smiles: String::from("CCC"),
                exact_mass: 1.0,
                name: String::from("name"),
                inchikey14: *b"AAAAAAAAAAAAAA",
            }),
            calculate_hash(&DatabaseRow {
                smiles: String::from("CCC"),
                exact_mass: 1.0,
                name: String::from("name"),
                inchikey14: *b"AAAAAAAABAAAAA",
            }),
        );
    }

    #[test]
    fn serde_bincode() {
        let row = DatabaseRow {
            smiles: String::from("CCC"),
            exact_mass: 1.0,
            name: String::from("name"),
            inchikey14: *b"AAAAAAAABAAAAA",
        };
        let bytes = bincode::serialize(&row).unwrap();
        let deser_row: DatabaseRow = bincode::deserialize(&bytes).unwrap();
        assert_eq!(deser_row, row);

        let row = DatabaseRow {
            smiles: String::from("CCC"),
            exact_mass: 1.0,
            name: String::from("name"),
            inchikey14: *b"AAAAAAAABAAAAA",
        };
        let bytes = bincode::serialize(&vec![row.clone()]).unwrap();
        let deser_row: Vec<DatabaseRow> = bincode::deserialize(&bytes).unwrap();
        assert_eq!(deser_row, vec![row]);
    }
}
