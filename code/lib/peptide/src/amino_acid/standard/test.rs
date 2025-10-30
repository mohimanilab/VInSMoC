use crate::amino_acid::{backbone, AminoAcid, StdAminoAcid};

use approx::assert_abs_diff_eq;
use constants::{
    amino_acids::{MASS_PEPTIDIC_BACKBONE, MASS_PEPTIDIC_BACKBONE_EXACT},
    MASS_WATER, MASS_WATER_EXACT,
};
use molecule::BondType::{Double, Single};
use molecule::{parsers, Atom, Mol, MolGraph, NodeIndex};

use petgraph::{algo::is_isomorphic_matching, stable_graph::node_index, Graph};
use rustc_hash::FxHashMap;

#[test]
fn check_fromstr() {
    assert_eq!(StdAminoAcid::R, "R".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::H, "H".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::K, "K".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::D, "D".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::E, "E".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::S, "S".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::T, "T".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::N, "N".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::Q, "Q".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::C, "C".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::G, "G".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::P, "P".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::A, "A".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::V, "V".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::I, "I".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::L, "L".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::M, "M".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::F, "F".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::Y, "Y".parse::<StdAminoAcid>().unwrap());
    assert_eq!(StdAminoAcid::W, "W".parse::<StdAminoAcid>().unwrap());
}

#[test]
fn check_aa_std_mass() {
    // Compare calculated mass with Pubchem standard mass, allow 0.01 measurement error
    let r_mass = StdAminoAcid::R.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 174.2, epsilon = 0.01);

    let r_mass = StdAminoAcid::H.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 155.15, epsilon = 0.01);

    let r_mass = StdAminoAcid::K.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 146.19, epsilon = 0.01);

    let r_mass = StdAminoAcid::D.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 133.1, epsilon = 0.01);

    let r_mass = StdAminoAcid::E.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 147.13, epsilon = 0.01);

    let r_mass = StdAminoAcid::S.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 105.09, epsilon = 0.01);

    let r_mass = StdAminoAcid::T.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 119.12, epsilon = 0.01);

    let r_mass = StdAminoAcid::N.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 132.12, epsilon = 0.01);

    let r_mass = StdAminoAcid::Q.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 146.14, epsilon = 0.01);

    let r_mass = StdAminoAcid::C.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 121.16, epsilon = 0.01);

    let r_mass = StdAminoAcid::G.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 75.07, epsilon = 0.01);

    let r_mass = StdAminoAcid::P.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 115.13, epsilon = 0.01);

    let r_mass = StdAminoAcid::A.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 89.09, epsilon = 0.01);

    let r_mass = StdAminoAcid::V.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 117.15, epsilon = 0.01);

    let r_mass = StdAminoAcid::I.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 131.17, epsilon = 0.01);

    let r_mass = StdAminoAcid::L.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 131.17, epsilon = 0.01);

    let r_mass = StdAminoAcid::M.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 149.21, epsilon = 0.01);

    let r_mass = StdAminoAcid::F.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 165.19, epsilon = 0.01);

    let r_mass = StdAminoAcid::Y.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 181.19, epsilon = 0.01);

    let r_mass = StdAminoAcid::W.r_group_std_mass() + MASS_WATER + backbone().std_mass();
    assert_abs_diff_eq!(r_mass, 204.22, epsilon = 0.01);
}

#[test]
fn check_aa_exact_mass() {
    // Compare calculated mass with Pubchem mass, allow 0.1~0.2 measurement error
    let r_mass = StdAminoAcid::R.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 174.111676, epsilon = 0.00001);

    let r_mass = StdAminoAcid::H.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 155.069477, epsilon = 0.00001);

    let r_mass = StdAminoAcid::K.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 146.105528, epsilon = 0.00001);

    let r_mass = StdAminoAcid::D.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 133.037508, epsilon = 0.00001);

    let r_mass = StdAminoAcid::E.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 147.053158, epsilon = 0.00001);

    let r_mass = StdAminoAcid::S.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 105.042593, epsilon = 0.00001);

    let r_mass = StdAminoAcid::T.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 119.058243, epsilon = 0.00001);

    let r_mass = StdAminoAcid::N.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 132.053492, epsilon = 0.00001);

    let r_mass = StdAminoAcid::Q.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 146.069142, epsilon = 0.00001);

    let r_mass = StdAminoAcid::C.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 121.01975, epsilon = 0.00001);

    let r_mass = StdAminoAcid::G.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 75.032028, epsilon = 0.00001);

    let r_mass = StdAminoAcid::P.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 115.063329, epsilon = 0.00001);

    let r_mass = StdAminoAcid::A.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 89.047678, epsilon = 0.00001);

    let r_mass = StdAminoAcid::V.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 117.078979, epsilon = 0.00001);

    let r_mass = StdAminoAcid::I.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 131.094629, epsilon = 0.00001);

    let r_mass = StdAminoAcid::L.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 131.094629, epsilon = 0.00001);

    let r_mass = StdAminoAcid::M.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 149.05105, epsilon = 0.00001);

    let r_mass = StdAminoAcid::F.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 165.078979, epsilon = 0.00001);

    let r_mass = StdAminoAcid::Y.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 181.073893, epsilon = 0.00001);

    let r_mass = StdAminoAcid::W.r_group_exact_mass() + MASS_WATER_EXACT + backbone().exact_mass();
    assert_abs_diff_eq!(r_mass, 204.089878, epsilon = 0.00001);
}

#[test]
// Check if the mol is correct. If yes, should have same weight as peptide mass constant.
fn check_backbone_mass() {
    let backbone = backbone();
    println!("{}", backbone.std_mass());
    assert_abs_diff_eq!(
        backbone.std_mass(),
        MASS_PEPTIDIC_BACKBONE,
        epsilon = 0.0001
    );
    assert_abs_diff_eq!(
        backbone.exact_mass(),
        MASS_PEPTIDIC_BACKBONE_EXACT,
        epsilon = 0.0001
    );
}

#[test]
fn check_to_residue_mol() {
    // test on lysine
    let mut lysine = MolGraph::with_capacity(10, 9);
    lysine.add_node(Atom::C);
    lysine.add_node(Atom::N);
    lysine.add_node(Atom::C);
    lysine.add_node(Atom::O);

    lysine.add_node(Atom::C);
    lysine.add_node(Atom::C);
    lysine.add_node(Atom::C);
    lysine.add_node(Atom::C);
    lysine.add_node(Atom::N);

    lysine.add_edge(node_index(0), node_index(1), Single);
    lysine.add_edge(node_index(0), node_index(2), Single);
    lysine.add_edge(node_index(2), node_index(3), Double);

    lysine.add_edge(node_index(0), node_index(4), Single);
    lysine.add_edge(node_index(4), node_index(5), Single);
    lysine.add_edge(node_index(5), node_index(6), Single);
    lysine.add_edge(node_index(6), node_index(7), Single);
    lysine.add_edge(node_index(7), node_index(8), Single);

    let hydrogen_hash: FxHashMap<NodeIndex, u8> = [
        (node_index(0), 1),
        (node_index(1), 1),
        (node_index(4), 2),
        (node_index(5), 2),
        (node_index(6), 2),
        (node_index(7), 2),
        (node_index(8), 2),
    ]
    .iter()
    .cloned()
    .collect();

    let lysine_predicted: Mol = StdAminoAcid::K.to_residue_mol();

    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &lysine.into(),
        &lysine_predicted.clone().graph.into(),
        |x, y| x == y,
        |x, y| x == y
    ));
    assert_eq!(hydrogen_hash, lysine_predicted.hydrogens().clone());

    // special case 1: glycine
    let mut glycine = MolGraph::with_capacity(10, 9);
    glycine.add_node(Atom::C);
    glycine.add_node(Atom::N);
    glycine.add_node(Atom::C);
    glycine.add_node(Atom::O);

    glycine.add_edge(node_index(0), node_index(1), Single);
    glycine.add_edge(node_index(0), node_index(2), Single);
    glycine.add_edge(node_index(2), node_index(3), Double);

    let hydrogen_hash: FxHashMap<NodeIndex, u8> = [(node_index(0), 2), (node_index(1), 1)]
        .iter()
        .cloned()
        .collect();

    let glycine_predicted: Mol = StdAminoAcid::G.to_residue_mol();

    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &glycine.into(),
        &glycine_predicted.clone().graph.into(),
        |x, y| x == y,
        |x, y| x == y
    ));
    assert_eq!(hydrogen_hash, glycine_predicted.hydrogens().clone());

    // special case 1: Proline
    let mut proline = MolGraph::with_capacity(10, 9);
    proline.add_node(Atom::C);
    proline.add_node(Atom::N);
    proline.add_node(Atom::C);
    proline.add_node(Atom::O);

    proline.add_node(Atom::C);
    proline.add_node(Atom::C);
    proline.add_node(Atom::C);

    proline.add_edge(node_index(0), node_index(1), Single);
    proline.add_edge(node_index(0), node_index(2), Single);
    proline.add_edge(node_index(2), node_index(3), Double);

    proline.add_edge(node_index(0), node_index(4), Single);
    proline.add_edge(node_index(4), node_index(5), Single);
    proline.add_edge(node_index(5), node_index(6), Single);
    proline.add_edge(node_index(6), node_index(1), Single);

    let hydrogen_hash: FxHashMap<NodeIndex, u8> = [
        (node_index(0), 1),
        (node_index(4), 2),
        (node_index(5), 2),
        (node_index(6), 2),
    ]
    .iter()
    .cloned()
    .collect();

    let proline_predicted: Mol = StdAminoAcid::P.to_residue_mol();

    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &proline.into(),
        &proline_predicted.clone().graph.into(),
        |x, y| x == y,
        |x, y| x == y
    ));
    assert_eq!(hydrogen_hash, proline_predicted.hydrogens().clone());
}

#[test]
fn check_full_molecule_w_smiles() {
    let smiles = "CC([C]=O)[NH]";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::A.to_residue_mol();
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

    let smiles = "N=C(N)NCCCC([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::R.to_residue_mol();
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

    let smiles = "NC(=O)CC([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::N.to_residue_mol();
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

    let smiles = "[NH]C([C]=O)CC(=O)O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::D.to_residue_mol();
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

    let smiles = "[NH]C([C]=O)CS";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::C.to_residue_mol();
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

    let smiles = "NC(=O)CCC([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::Q.to_residue_mol();
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

    let smiles = "[NH]C([C]=O)CCC(=O)O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::E.to_residue_mol();
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

    let smiles = "[NH]C[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::G.to_residue_mol();
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

    let smiles = "C1=C(NC=N1)CC([C](=O))[NH]";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::H.to_residue_mol();
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

    let smiles = "CCC(C)C([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::I.to_residue_mol();
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

    let smiles = "CC(C)CC([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::L.to_residue_mol();
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

    let smiles = "NCCCCC([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::K.to_residue_mol();
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

    let smiles = "CSCCC([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::M.to_residue_mol();
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

    let smiles = "C1=CC=C(C=C1)CC([C](=O))[NH]";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::F.to_residue_mol();
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

    let smiles = "O=[C]C1CCC[N]1";
    let mol = parsers::parse_smiles(smiles).unwrap();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::P.to_residue_mol();
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

    let smiles = "[NH]C([C]=O)CO";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::S.to_residue_mol();
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

    let smiles = "CC(O)C([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::T.to_residue_mol();
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

    let smiles = "C1=CC=C2C(=C1)C(=CN2)CC([C](=O))[NH]";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::W.to_residue_mol();
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

    let smiles = "C1=CC(=CC=C1CC([C](=O))[NH])O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::Y.to_residue_mol();
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

    let smiles = "CC(C)C([NH])[C]=O";
    let mut mol = parsers::parse_smiles(smiles).unwrap();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |node_idx, node| (node, mol.h_count(&node_idx)),
        |_, edge| edge,
    );
    let aa = StdAminoAcid::V.to_residue_mol();
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
}

#[test]
// Check that one hot encoding works properly
fn check_ohe() {
    let aa = StdAminoAcid::A;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 0 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::C;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 1 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::D;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 2 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::E;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 3 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::F;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 4 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::G;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 5 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::H;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 6 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::I;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 7 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::K;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 8 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::L;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 9 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::M;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 10 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::N;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 11 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::P;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 12 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::Q;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 13 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::R;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 14 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::S;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 15 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::T;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 16 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::V;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 17 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::W;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 18 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }

    let aa = StdAminoAcid::Y;
    let ohe = aa.one_hot_encoding();
    assert!(ohe.shape() == &[20]);
    for i in 0..20 {
        if i == 19 {
            assert!(ohe[i] == 1);
        } else {
            assert!(ohe[i] == 0);
        }
    }
}
