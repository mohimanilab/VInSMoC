use super::*;
use crate::parsers::smiles::parse_smiles;

use petgraph::{algo::is_isomorphic_matching, Graph};
use std::fs;

macro_rules! mol_file_smiles_equality_helper {
    ($smiles:literal, $mol_file_path:literal) => {
        let mut mol_smiles: Mol = parse_smiles($smiles).unwrap();
        mol_smiles.explicit_hs();
        let content = fs::read_to_string($mol_file_path).expect("Unable to read file");
        assert!(parse_mol_file(&content).is_ok());
        let mol_file_parsed = parse_mol_file(&content).unwrap();

        assert!(is_isomorphic_matching::<
            &Graph<_, _, _, _>,
            &Graph<_, _, _, _>,
            _,
            _,
        >(
            &mol_smiles.graph.into(),
            &mol_file_parsed.graph.into(),
            |x, y| x == y,
            |x, y| x == y
        ));
        assert_eq!(mol_smiles.charges, mol_file_parsed.charges, "charges");
        assert_eq!(mol_smiles.hydrogens, mol_file_parsed.hydrogens, "hydrogens");
    };
}

#[test]
fn parse_mol_file_works() {
    // all molfiles have hydrogens explictly added by RDKit
    // https://pubchem.ncbi.nlm.nih.gov/compound/2244
    mol_file_smiles_equality_helper!("CC(=O)OC1=CC=CC=C1C(=O)O", "test_files/aspirin.mol");

    // https://pubchem.ncbi.nlm.nih.gov/compound/241
    mol_file_smiles_equality_helper!("C1=CC=CC=C1", "test_files/benzene.mol");

    // https://pubchem.ncbi.nlm.nih.gov/compound/687061
    mol_file_smiles_equality_helper!("C1C(OC2=CC=CC=C2O1)C(=O)O", "test_files/ex1.mol");

    // https://pubchem.ncbi.nlm.nih.gov/compound/240
    mol_file_smiles_equality_helper!("C1=CC=C(C=C1)C=O", "test_files/benzaldehyde.mol");

    // https://pubchem.ncbi.nlm.nih.gov/compound/9549299
    mol_file_smiles_equality_helper!(
        "C1CC1C(=O)NC2=CC=CC(=C2)NC3=NC=NC(=C3)NC4=CC=CC(=C4)C(F)(F)F",
        "test_files/egfr.mol"
    );

    // https://pubchem.ncbi.nlm.nih.gov/compound/79842
    mol_file_smiles_equality_helper!("CCCCN(CCCC)CCC1=CC=CC=C1", "test_files/ex2.mol");

    // https://pubchem.ncbi.nlm.nih.gov/compound/389
    mol_file_smiles_equality_helper!("C(CC(C(=O)O)N)CN", "test_files/ex3.mol");

    // https://pubchem.ncbi.nlm.nih.gov/compound/305
    mol_file_smiles_equality_helper!("C[N+](C)(C)CCO", "test_files/choline.mol");

    // https://pubchem.ncbi.nlm.nih.gov/compound/223
    mol_file_smiles_equality_helper!("[NH4+]", "test_files/ammonium.mol");

    // just chlorine
    // FIXME: this fails and reveals a pretty horrifying state of the MOL file parser
    // mol_file_smiles_equality_helper!("[Cl-]", "test_files/negative_charge.mol");
}
