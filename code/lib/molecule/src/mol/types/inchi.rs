use std::{
    convert::TryInto,
    ffi::{CStr, CString},
    os::raw::c_char,
};

use libinchi_sys::{
    inchi_Atom, inchi_Input, inchi_Output, tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE,
    tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE, tagINCHIBondType_INCHI_BOND_TYPE_NONE,
    tagINCHIBondType_INCHI_BOND_TYPE_SINGLE, tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE,
    tagINCHIRadical_INCHI_RADICAL_NONE, FreeStdINCHI, GetStdINCHI, GetStdINCHIKeyFromStdINCHI,
    AT_NUM, NO_ATOM, S_CHAR,
};
use petgraph::visit::EdgeRef;
use rustc_hash::FxHashMap;

use crate::{BondType, Mol, NodeIndex};

impl Mol {
    /// Converts a molecule to its standard InChI string.
    ///
    /// No stereochemical or 3D information is included when construction the InChI string. The
    /// returned result uses the libinchi error code as its `Err` variant. In the case of an error
    /// this function will print the InChI logs to stderr. This is behavior that will be improved
    /// in the future, possibly in a separate `inchi-sys` crate.
    pub fn to_inchi(&self) -> Result<String, i32> {
        let indices = self
            .graph
            .node_indices()
            .enumerate()
            .map(|(x, y)| (y, x.try_into().unwrap()))
            .collect::<FxHashMap<NodeIndex, i16>>();
        let mut inchi_atoms = vec![None; self.graph.node_count()];

        for (&old_idx, &new_idx) in indices.iter() {
            let mut neighbor: [AT_NUM; 20usize] = [NO_ATOM as i16; 20];
            let mut bond_type: [S_CHAR; 20usize] =
                [tagINCHIBondType_INCHI_BOND_TYPE_NONE as i8; 20usize];
            let mut num_bonds: AT_NUM = 0;
            for (idx, neigh) in self.graph.edges(old_idx).enumerate() {
                neighbor[idx] = indices[&neigh.target()];
                bond_type[idx] = match neigh.weight() {
                    BondType::Single => tagINCHIBondType_INCHI_BOND_TYPE_SINGLE,
                    BondType::Double => tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE,
                    BondType::Triple => tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE,
                    BondType::Aromatic => {
                        unimplemented!("Aromatic bond types fall apart for InCHIs")
                    }
                } as i8;
                num_bonds += 1;
            }

            let mut elname: [c_char; 6] = [0; 6];
            for (idx, atom_char) in format!("{:?}", self.graph[old_idx]).chars().enumerate() {
                elname[idx] = atom_char as i8;
            }

            inchi_atoms[new_idx as usize] = Some(inchi_Atom {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                neighbor,
                bond_type,
                bond_stereo: [tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE as i8; 20],
                elname,
                num_bonds,
                num_iso_H: [self.h_count(&old_idx).try_into().unwrap(), 0, 0, 0],
                isotopic_mass: 0,
                radical: tagINCHIRadical_INCHI_RADICAL_NONE as i8,
                charge: self.charge(&old_idx),
            });
        }

        let mut inchi_atoms = inchi_atoms
            .into_iter()
            .map(|x| x.unwrap())
            .collect::<Vec<_>>();
        let inchi_atoms = inchi_atoms.as_mut_ptr();
        let mut inchi_input = inchi_Input {
            atom: inchi_atoms,
            stereo0D: std::ptr::null_mut(),
            szOptions: std::ptr::null_mut(),
            num_atoms: self.graph.node_count() as AT_NUM,
            num_stereo0D: 0,
        };
        let inchi_input_ptr: *mut _ = &mut inchi_input;
        let mut inchi_output = inchi_Output {
            szInChI: std::ptr::null_mut(),
            szAuxInfo: std::ptr::null_mut(),
            szMessage: std::ptr::null_mut(),
            szLog: std::ptr::null_mut(),
        };
        let inchi_output_ptr: *mut _ = &mut inchi_output;

        let ret_code = unsafe { GetStdINCHI(inchi_input_ptr, inchi_output_ptr) };

        if ret_code != 0 && ret_code != 1 {
            unsafe {
                eprintln!("##### MESSAGE");
                eprintln!(
                    "{}",
                    CStr::from_ptr(inchi_output.szMessage).to_str().unwrap()
                );
                eprintln!("##### LOG");
                eprintln!("{}", CStr::from_ptr(inchi_output.szLog).to_str().unwrap());
                FreeStdINCHI(inchi_output_ptr);
            };
            return Err(ret_code);
        }

        let output = String::from(
            unsafe { CStr::from_ptr(inchi_output.szInChI) }
                .to_str()
                .unwrap(),
        ); //.into_string().unwrap();
        unsafe {
            FreeStdINCHI(inchi_output_ptr);
        };
        Ok(output)
    }

    /// Converts a molecule into its 14-character InChIKey prefix.
    ///
    /// This can be used to check molecules for uniqueness. The `Err` variant of the returned
    /// result tracks the error code from the underlying libinchi implementation. This is a very
    /// shaky, ugly implementation right now and can panic if anything goes wrong in conversion.
    pub fn to_inchi_key(&self) -> Result<String, i32> {
        let inchi = CString::new(self.to_inchi()?)
            .expect("Couldn't create CString")
            .into_raw();
        let inchikey = CString::new(vec![1u8; 28])
            .expect("Couldn't create CString")
            .into_raw();
        let retcode = unsafe { GetStdINCHIKeyFromStdINCHI(inchi, inchikey) };
        // take ownership of inchi again, this is just to drop it
        // TODO: there's probably a better way to do this
        let _inchi = unsafe { CString::from_raw(inchi) };
        let result = unsafe { CString::from_raw(inchikey) };

        if retcode != 0 {
            return Err(retcode);
        }
        let mut result = result.into_string().unwrap();
        result.truncate(14);
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use petgraph::graph::node_index;

    use crate::MolBuilder;

    #[test]
    fn to_inchi() {
        let m = MolBuilder::from_smiles("C").unwrap().build();
        let inchi = m.to_inchi().unwrap();
        assert_eq!(inchi, "InChI=1S/CH4/h1H4");

        let m = MolBuilder::from_smiles("CNOC").unwrap().build();
        let inchi = m.to_inchi().unwrap();
        assert_eq!(inchi, "InChI=1S/C2H7NO/c1-3-4-2/h3H,1-2H3");

        let mut m = MolBuilder::from_smiles("CC").unwrap().build();
        m.graph.remove_node(node_index(0));
        m.hydrogens.remove(&node_index(0));
        m.charges.remove(&node_index(0));
        *m.hydrogens.get_mut(&node_index(1)).unwrap() += 1;
        let inchi = m.to_inchi().unwrap();
        assert_eq!(inchi, "InChI=1S/CH4/h1H4");

        let mut m = MolBuilder::from_smiles("CC").unwrap().build();
        m.graph.remove_node(node_index(1));
        m.hydrogens.remove(&node_index(1));
        m.charges.remove(&node_index(1));
        *m.hydrogens.get_mut(&node_index(0)).unwrap() += 1;
        let inchi = m.to_inchi().unwrap();
        assert_eq!(inchi, "InChI=1S/CH4/h1H4");
    }

    #[test]
    fn to_inchi_key() {
        let m = MolBuilder::from_smiles("C").unwrap().build();
        let inchi = m.to_inchi_key().unwrap();
        assert_eq!(inchi, "VNWKTOKETHGBQD");

        let m = MolBuilder::from_smiles("CNOC").unwrap().build();
        let inchi = m.to_inchi_key().unwrap();
        assert_eq!(inchi, "KRKPYFLIYNGWTE");

        let m = MolBuilder::from_smiles("CONC").unwrap().build();
        let inchi = m.to_inchi_key().unwrap();
        assert_eq!(inchi, "KRKPYFLIYNGWTE");

        let m = MolBuilder::from_smiles("C=CC(C)(C)OC(C)C1NC(=O)C2CSC(=N2)C(CC(C)C)NC(=O)C(C(C)OC(C)(C)C=C)NC(=O)C2CCCN2C(=O)C(C(C)C)NC1=O").unwrap().build();
        let inchi = m.to_inchi_key().unwrap();
        assert_eq!(inchi, "PXBQBGADJRIZRZ");
    }
}
