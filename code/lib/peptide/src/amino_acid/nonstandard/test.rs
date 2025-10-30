use super::*;
use molecule::parsers::parse_smiles;
use petgraph::{algo::is_isomorphic_matching, Graph};

macro_rules! check_smiles {
    ($aa:expr, $smiles:literal) => {
        let mut mol = parse_smiles($smiles).unwrap();
        mol.implicit_hs();
        let mol_graph = mol.graph.map(
            |node_idx, node| (node, mol.h_count(&node_idx)),
            |_, edge| edge,
        );
        let aa = $aa.to_residue_mol();
        let aa_graph = aa.graph.map(
            |node_idx, node| (node, aa.h_count(&node_idx)),
            |_, edge| edge,
        );
        assert!(is_isomorphic_matching::<
            &Graph<_, _, _, _>,
            &Graph<_, _, _, _>,
            _,
            _,
        >(
            &mol_graph.into(),
            &aa_graph.into(),
            |x, y| x == y,
            |x, y| x == y
        ));
    };
}

#[test]
fn check_full_molecule_w_smiles() {
    check_smiles!(NonStdAminoAcid::Aad, "C(CC([C](=O))[NH])CC(=O)O");
    check_smiles!(NonStdAminoAcid::Abu, "CCC([NH])[C](=O)");
    check_smiles!(NonStdAminoAcid::Bht, "[NH]C([C](=O))C(O)C1=CC=C(O)C=C1");
    check_smiles!(NonStdAminoAcid::Dpg, "[NH]C([C](=O))C1=CC(O)=CC(O)=C1");
    check_smiles!(NonStdAminoAcid::Hpg, "[NH]C([C](=O))C1C=CC(O)=CC=1");
    check_smiles!(NonStdAminoAcid::Orn, "NCCC[C@H]([NH])[C](=O)");
    check_smiles!(NonStdAminoAcid::Phg, "[NH]C([C](=O))C1=CC=CC=C1");
    check_smiles!(NonStdAminoAcid::Pip, "O=[C]C1CCCC[N]1");
    check_smiles!(NonStdAminoAcid::Cit, "NC(=O)NCCC[C@H]([NH])[C](=O)");
    check_smiles!(NonStdAminoAcid::Dab, "NCC[C@H]([NH])[C](=O)");
    check_smiles!(NonStdAminoAcid::Hty, "[NH]C(CCC1=CC=C(O)C=C1)[C](=O)");
    check_smiles!(NonStdAminoAcid::Hse, "[NH]C(CCO)[C](=O)");
    check_smiles!(NonStdAminoAcid::Hrn, "ONCCC[C@H]([NH])[C](=O)");
    check_smiles!(NonStdAminoAcid::Han, "ON(C(=O)C)CCC[C@H]([NH])[C](=O)");
    check_smiles!(NonStdAminoAcid::Adh, "CC(C)=CC(C)C([C]=O)[NH1]");
    check_smiles!(
        NonStdAminoAcid::Dmw,
        "C=CC(C)(C)N1C=C(CC([NH])[C]=O)C2=CC=CC=C21"
    );
    check_smiles!(NonStdAminoAcid::Dhb, "CC=C([NH])[C](=O)");
    check_smiles!(
        NonStdAminoAcid::Crp,
        "C1(Cl)=CC=C2C(=C1)C(=CN2)CC([C](=O))[NH]"
    );
    check_smiles!(NonStdAminoAcid::Mha, "CC(C)CCC([C](=O))[NH]");
    check_smiles!(NonStdAminoAcid::Hfn, "ON(C(=O))CCC[C@H]([NH])[C](=O)");
    check_smiles!(NonStdAminoAcid::Mpr, "O=[C]C1CC(C)C[N]1");
    check_smiles!(NonStdAminoAcid::Bla, "C([NH])C[C](=O)");
    check_smiles!(NonStdAminoAcid::Bbu, "C([NH])(C)C[C](=O)");
}
