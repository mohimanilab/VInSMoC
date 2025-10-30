use crate::{mods::parse, ModSeq};
use molecule::{Mol, MolBuilder, NodeIndex};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// A modification specified by a subgraph pattern.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubgraphMod {
    /// The pattern molecule
    pub pattern: Mol,
    /// The modifications applied
    pub mod_seq: ModSeq,
    /// The name of this modification
    pub name: String,
}

impl SubgraphMod {
    /// Create a new modification from a pattern graph and a modification sequence.
    pub fn new(pattern: Mol, mod_seq: ModSeq, name: String) -> Self {
        Self {
            pattern,
            mod_seq,
            name,
        }
    }

    /// Generate a new [`SubgraphMod`] from a mol file and a mod file and a name.
    ///
    /// This function clones the input molecule twice. The first clone is to make sure that the
    /// input molecule has implicit hydrogens. The second clone is to keep track of the
    /// modifications with a `checker` molecule, to which the modifications are sequentially
    /// applied. Any potential errors are propagated.
    pub fn from_mol_and_mod_file(
        mol_file: impl AsRef<Path>,
        mod_file: impl AsRef<Path>,
        name: &str,
    ) -> anyhow::Result<Self> {
        // set up a checker molecule to sequentially apply modifications and make sure they are
        // valid
        let mut mol = MolBuilder::from_mol_file(mol_file)?.build();
        let mut checker = mol.clone();
        checker.implicit_hs();

        let mod_file = PathBuf::from(mod_file.as_ref());

        let fields = parse::parse_mod(&mod_file)?;
        let mut mods = Vec::with_capacity(fields.len());
        let mut new_idx = -1;

        for (line, field) in fields.into_iter().enumerate() {
            // people expect 1-based line numbers
            let line = line + 1;

            // get the mods that correspond to this line
            let new_mods = field
                .apply(&mol, &mut new_idx)
                .map_err(|source| crate::Error::Validation { line, source })
                .map_err(|source| crate::ContextualError::Mods {
                    mod_dir: mod_file.clone(),
                    source,
                })?;

            // if there are any, then check that we can apply them to the Mol
            if let Some(new_mods) = new_mods {
                ModSeq::new(new_mods.clone())
                    .modify(&mut checker)
                    .map_err(|source| crate::ContextualError::Application {
                        mod_dir: mod_file.clone(),
                        line,
                        source,
                    })?;
                mods.extend(new_mods);
            }
        }
        mol.implicit_hs();

        let (mol, index_map) = mol.reset_indices();
        let index_map = index_map
            .into_iter()
            .map(|(k, v)| (v.index() as i32, k.index() as i32))
            .collect();

        // and build the final modification sequence
        let mods = mods
            .into_iter()
            .map(|x| x.convert(&index_map).unwrap())
            .collect::<ModSeq>();

        // check it still works
        mods.modify(&mut (mol.clone())).map_err(|source| {
            crate::ContextualError::TransApplication {
                mod_dir: mod_file.clone(),
                source,
            }
        })?;

        Ok(Self::new(mol, mods, String::from(name)))
    }

    /// Find all locations of the input molecule that this modification maps to.
    ///
    /// Results are returned as a vector of index mappings that convert internal pattern indices
    /// into the target molecule's indices.
    ///
    /// See also [`molecule::Mol::ri_subgraph_isomorphism`].
    pub fn conversions(&self, mol: &Mol) -> Vec<FxHashMap<NodeIndex, NodeIndex>> {
        mol.ri_subgraph_isomorphism(&self.pattern)
    }
}
