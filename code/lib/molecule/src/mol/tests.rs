use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::MolBuilder;
use crate::{
    atoms::Atom,
    bonds::BondType::{Double, Single},
    parsers::{parse_mol_file, parse_smiles},
    BondType, Mol, NodeIndex, OrganicAtom,
};

use approx::assert_abs_diff_eq;

use itertools::Itertools;
use petgraph::{
    algo::{is_isomorphic_matching, tarjan_scc},
    stable_graph::{node_index, StableGraph, StableUnGraph},
    Graph,
};
use std::collections::HashSet;

use super::types::ValenceError;

#[test]
fn mol_mass_works() {
    // iodine needs a lower precision, since we don't include the error digits
    assert_abs_diff_eq!(
        parse_smiles("CI").unwrap().exact_mass(),
        141.927948096,
        epsilon = 0.00001
    );

    let f = File::open("test_files/test_chem_db.csv").unwrap();
    let br = BufReader::new(f);

    // compare to RDKit
    for (i, line) in br.lines().enumerate() {
        if i > 0 {
            let l = line.unwrap();
            let data: Vec<&str> = l.split(",").collect();
            let smiles: &str = data[0];
            let reference_mass: f64 = data[1].parse().unwrap();
            let mol: Mol = parse_smiles(smiles).unwrap();

            // these are due to differences in precision with RDKit in Iodine precision
            // they use one more digit than we do, but the current (21 Aug 2018) version of CIAAW
            // reports that digit as error, so we don't use it. Because of this we can safely
            // ignore the following indices
            if i == 94 || i == 70 {
                continue;
            }
            assert_abs_diff_eq!(mol.exact_mass(), reference_mass, epsilon = 0.000001);
        }
    }

    // compare to pubchem
    // https://pubchem.ncbi.nlm.nih.gov/compound/75791
    let mol = parse_smiles("C1=CC=C(C=C1)C(C(=O)O)C(=O)O").unwrap();
    assert_abs_diff_eq!(mol.exact_mass(), 180.042259, epsilon = 0.000001);
    assert_eq!(mol.formal_charge(), 0);

    // https://pubchem.ncbi.nlm.nih.gov/compound/420
    let mol = parse_smiles("CC(C)CCCC(C)C1CCC2C1(CCC3C2=CCC4C3(CCC(C4)O)C)C").unwrap();
    assert_abs_diff_eq!(mol.exact_mass(), 386.354866, epsilon = 0.000001);
    assert_eq!(mol.formal_charge(), 0);

    // https://pubchem.ncbi.nlm.nih.gov/compound/6934
    let mol = parse_smiles("CC1=CC(=C(C=C1)N)S(=O)(=O)O").unwrap();
    assert_abs_diff_eq!(mol.exact_mass(), 187.030314, epsilon = 0.000001);
    assert_eq!(mol.formal_charge(), 0);

    // https://pubchem.ncbi.nlm.nih.gov/compound/12345
    let mol = parse_smiles("CC(=O)OCOC(=O)C").unwrap();
    assert_abs_diff_eq!(mol.exact_mass(), 132.042259, epsilon = 0.000001);
    assert_eq!(mol.formal_charge(), 0);

    // https://pubchem.ncbi.nlm.nih.gov/compound/23409
    let mol = parse_smiles("C1=CC=C(C(=C1)C=C(C#N)C#N)O").unwrap();
    assert_abs_diff_eq!(mol.exact_mass(), 170.048013, epsilon = 0.000001);

    // https://pubchem.ncbi.nlm.nih.gov/compound/10000
    let mol = parse_smiles("C(CCl)C(F)(F)F").unwrap();
    assert_abs_diff_eq!(mol.exact_mass(), 131.995362, epsilon = 0.000001);

    // https://pubchem.ncbi.nlm.nih.gov/compound/76823
    let mol = parse_smiles("CCCCCCCCCCCCCCCCSCCCCCCCCCCCCCCCC").unwrap();
    assert_abs_diff_eq!(mol.exact_mass(), 482.488523, epsilon = 0.000001);

    // Ammonium
    // https://pubchem.ncbi.nlm.nih.gov/compound/Ammonium-ion
    let mol = parse_smiles("[NH4+]").unwrap();
    assert_eq!(mol.formal_charge(), 1);

    // Sulfate
    // https://pubchem.ncbi.nlm.nih.gov/compound/Sulfate
    let mol = parse_smiles("[O-]S(=O)(=O)[O-]").unwrap();
    assert_eq!(mol.formal_charge(), -2);

    // Phosphate
    // https://pubchem.ncbi.nlm.nih.gov/compound/Phosphate
    let mol = parse_smiles("[O-]P(=O)([O-])[O-]").unwrap();
    assert_eq!(mol.formal_charge(), -3);
}

#[test]
fn explict_hs() {
    // normal molecule
    let mut mol = parse_smiles("CC").unwrap();
    assert_eq!(mol.graph.node_count(), 2);
    assert_eq!(
        mol.graph
            .node_indices()
            .filter(|x| mol.graph[*x] == Atom::H)
            .count(),
        0
    );
    assert_eq!(
        mol.hydrogens.values().copied().collect::<HashSet<u8>>(),
        [3].iter().copied().collect::<HashSet<u8>>()
    );
    let mut ref_graph = StableUnGraph::<Atom, BondType, u16>::with_capacity(8, 7);
    ref_graph.add_node(Atom::C);
    ref_graph.add_node(Atom::C);
    ref_graph.add_edge(node_index(0), node_index(1), BondType::Single);
    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &mol.graph.clone().into(),
        &ref_graph.into(),
        |x, y| x == y,
        |x, y| x == y,
    ));

    // now add in the explicit Hs
    mol.explicit_hs();
    assert_eq!(mol.graph.node_count(), 8);
    assert_eq!(
        mol.graph
            .node_indices()
            .filter(|x| mol.graph[*x] == Atom::H)
            .count(),
        6
    );
    assert!(mol.hydrogens.is_empty());
    let mut ref_graph = StableUnGraph::<Atom, BondType, u16>::with_capacity(8, 7);
    ref_graph.add_node(Atom::C);
    ref_graph.add_node(Atom::C);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_edge(node_index(0), node_index(1), BondType::Single);
    ref_graph.add_edge(node_index(0), node_index(2), BondType::Single);
    ref_graph.add_edge(node_index(0), node_index(3), BondType::Single);
    ref_graph.add_edge(node_index(0), node_index(4), BondType::Single);
    ref_graph.add_edge(node_index(1), node_index(5), BondType::Single);
    ref_graph.add_edge(node_index(1), node_index(6), BondType::Single);
    ref_graph.add_edge(node_index(1), node_index(7), BondType::Single);
    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &mol.graph.clone().into(),
        &ref_graph.into(),
        |x, y| x == y,
        |x, y| x == y,
    ));

    // double calling should be a no-op
    mol.explicit_hs();
    assert_eq!(mol.graph.node_count(), 8);
    assert_eq!(
        mol.graph
            .node_indices()
            .filter(|x| mol.graph[*x] == Atom::H)
            .count(),
        6
    );
    assert!(mol.hydrogens.is_empty());
    let mut ref_graph = StableUnGraph::<Atom, BondType, u16>::with_capacity(8, 7);
    ref_graph.add_node(Atom::C);
    ref_graph.add_node(Atom::C);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_node(Atom::H);
    ref_graph.add_edge(node_index(0), node_index(1), BondType::Single);
    ref_graph.add_edge(node_index(0), node_index(2), BondType::Single);
    ref_graph.add_edge(node_index(0), node_index(3), BondType::Single);
    ref_graph.add_edge(node_index(0), node_index(4), BondType::Single);
    ref_graph.add_edge(node_index(1), node_index(5), BondType::Single);
    ref_graph.add_edge(node_index(1), node_index(6), BondType::Single);
    ref_graph.add_edge(node_index(1), node_index(7), BondType::Single);
    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &mol.graph.into(),
        &ref_graph.into(),
        |x, y| x == y,
        |x, y| x == y,
    ));
}

macro_rules! round_trip_hs_smiles {
    ($smiles:literal) => {
        dbg!($smiles);
        let mut ref_mol = parse_smiles($smiles).unwrap();
        ref_mol.implicit_hs();
        let ref_graph = ref_mol.graph.map(
            |idx, atom| (atom, ref_mol.h_count(&idx), ref_mol.charge(&idx)),
            |_, bond| bond,
        );
        let mut mol = ref_mol.clone();
        mol.explicit_hs();
        mol.implicit_hs();
        let mol_graph = mol.graph.map(
            |idx, atom| (atom, ref_mol.h_count(&idx), ref_mol.charge(&idx)),
            |_, bond| bond,
        );
        assert!(is_isomorphic_matching::<
            &Graph<_, _, _, _>,
            &Graph<_, _, _, _>,
            _,
            _,
        >(
            &ref_graph.into(),
            &mol_graph.into(),
            |x, y| x == y,
            |x, y| x == y,
        ));
    };
}

#[test]
fn implicit_hs() {
    round_trip_hs_smiles!("CC");
    round_trip_hs_smiles!("CC1=C(C=CC(=C1)Cl)OCC(=O)OCC(C)C");
    round_trip_hs_smiles!("CN1CCC2=CC=CC=C2C3=C(C1)C=CC(=C3O)OC");
    round_trip_hs_smiles!("C[N+](C)(C1=CC=C(C=C1)N=NC2=CC=CC=C2)[O-]");
    round_trip_hs_smiles!("C1=CC=C2C(=C1)C3=CC=CC=C3C(=O)C2=O");
    round_trip_hs_smiles!("CC[NH+](CC)CCSC1=C2CCCCCC2=[NH+]C3=C1C=C(C=C3)Cl.[Cl-].[Cl-]");
    round_trip_hs_smiles!("CCOC(=O)NC(=O)CC#N");
    round_trip_hs_smiles!("[NH3+]O.[NH3+]O.[O-]S(=O)(=O)[O-]");
    round_trip_hs_smiles!("CN1CCN(CC1)CC2=CNC3=C2C=C(C=C3)[N+](=O)[O-]");
    round_trip_hs_smiles!("CN1CC2CCCC(C1)N2CCOC(=O)CCCCCCCC(=O)OCCN3C4CCCC3CN(C4)C");
    round_trip_hs_smiles!("CC(C)C(=O)NNC(C)C");
    round_trip_hs_smiles!("C1CN=C(N1)CC2=C3C(=CC=C2)C4=CC=CC=C4S3");
    round_trip_hs_smiles!("C1C(=O)N=C(S1)N");
    round_trip_hs_smiles!("C1=C(NC(=O)NC1=O)C(=O)[O-].[Na+]");
    round_trip_hs_smiles!(
        "C[N+]1(C2CCCC1CN(C2)CC3=CC=CC=C3)CCCCCCC[N+]4(C5CCCC4CN(C5)CC6=CC=CC=C6)C"
    );
    round_trip_hs_smiles!("CC(=O)N=C(N)S");
    round_trip_hs_smiles!("C1=CC(=CN=C1)CCCO");
    round_trip_hs_smiles!("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CCCCC3");
    round_trip_hs_smiles!("C1CCC(CC1)C2(OCC(O2)C3CCCCN3)C4CCCCC4");

    // make sure double-calling is a no-op
    let ref_mol = parse_smiles("CCCC").unwrap();
    let ref_graph = ref_mol.graph.map(
        |idx, atom| (atom, ref_mol.h_count(&idx), ref_mol.charge(&idx)),
        |_, bond| bond,
    );
    let mut mol = ref_mol.clone();
    mol.explicit_hs();
    mol.implicit_hs();
    mol.implicit_hs();
    let mol_graph = mol.graph.map(
        |idx, atom| (atom, ref_mol.h_count(&idx), ref_mol.charge(&idx)),
        |_, bond| bond,
    );
    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &ref_graph.into(),
        &mol_graph.into(),
        |x, y| x == y,
        |x, y| x == y,
    ));
}

macro_rules! infer_mol_file {
    ($smiles:literal, $mol_file_path:literal) => {
        let mut mol_smiles: Mol = parse_smiles($smiles).unwrap();
        mol_smiles.implicit_hs();
        let content = std::fs::read_to_string($mol_file_path).expect("Unable to read file");
        assert!(parse_mol_file(&content).is_ok());
        let mut mol_file_parsed = parse_mol_file(&content).unwrap();
        mol_file_parsed.infer_hs().unwrap();

        let mol_smiles_graph = mol_smiles.graph.map(
            |idx, atom| (atom, mol_smiles.h_count(&idx), mol_smiles.charge(&idx)),
            |_, bond| bond,
        );

        let mol_file_graph = mol_file_parsed.graph.map(
            |idx, atom| {
                (
                    atom,
                    mol_file_parsed.h_count(&idx),
                    mol_file_parsed.charge(&idx),
                )
            },
            |_, bond| bond,
        );

        dbg!(&mol_smiles_graph);
        dbg!(&mol_file_graph);

        assert!(is_isomorphic_matching::<
            &Graph<_, _, _, _>,
            &Graph<_, _, _, _>,
            _,
            _,
        >(
            &mol_smiles_graph.into(),
            &mol_file_graph.into(),
            |x, y| x == y,
            |x, y| x == y
        ));
        assert!(!mol_file_parsed
            .graph
            .node_indices()
            .any(|x| mol_file_parsed.graph[x] == Atom::H));
    };
}

#[test]
fn mol_file_infer_hs() {
    // this is basically the old MOL file parser test
    // https://pubchem.ncbi.nlm.nih.gov/compound/2244
    infer_mol_file!(
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "test_files/aspirin.mol.implicit"
    );
    // https://pubchem.ncbi.nlm.nih.gov/compound/241
    infer_mol_file!("C1=CC=CC=C1", "test_files/benzene.mol.implicit");
    // https://pubchem.ncbi.nlm.nih.gov/compound/687061
    infer_mol_file!("C1C(OC2=CC=CC=C2O1)C(=O)O", "test_files/ex1.mol.implicit");
    // https://pubchem.ncbi.nlm.nih.gov/compound/240
    infer_mol_file!("C1=CC=C(C=C1)C=O", "test_files/benzaldehyde.mol.implicit");
    // https://pubchem.ncbi.nlm.nih.gov/compound/9549299
    infer_mol_file!(
        "C1CC1C(=O)NC2=CC=CC(=C2)NC3=NC=NC(=C3)NC4=CC=CC(=C4)C(F)(F)F",
        "test_files/egfr.mol.implicit"
    );
    // https://pubchem.ncbi.nlm.nih.gov/compound/79842
    infer_mol_file!("CCCCN(CCCC)CCC1=CC=CC=C1", "test_files/ex2.mol.implicit");
    // https://pubchem.ncbi.nlm.nih.gov/compound/389
    infer_mol_file!("C(CC(C(=O)O)N)CN", "test_files/ex3.mol.implicit");
    // https://pubchem.ncbi.nlm.nih.gov/compound/305
    infer_mol_file!("C[N+](C)(C)CCO", "test_files/choline.mol.implicit");
    // https://pubchem.ncbi.nlm.nih.gov/compound/223
    infer_mol_file!("[NH4+]", "test_files/ammonium.mol.implicit");
}

#[test]
fn infer_hs() {
    // exceeded valence
    let mut mol = parse_smiles("[CH5]").unwrap();
    assert_eq!(mol.infer_hs(), Err(ValenceError::new(0, OrganicAtom::C)));

    // skip inorganic
    let mol = parse_smiles("[Au]").unwrap();
    let mut inferred = mol.clone();
    inferred.infer_hs().unwrap();
    assert_eq!(mol.hydrogens, inferred.hydrogens);

    // 4 valence nitrogen
    let mut mol = parse_smiles("[N]([C])([C])([C])([C])").unwrap();
    mol.infer_hs().unwrap();
    assert_eq!(mol.hydrogens[&node_index(0)], 1);
    assert_eq!(mol.hydrogens[&node_index(1)], 3);
    assert_eq!(mol.hydrogens[&node_index(2)], 3);
    assert_eq!(mol.hydrogens[&node_index(3)], 3);
    assert_eq!(mol.hydrogens[&node_index(4)], 3);

    // check it also works on explicit versions
    // 4 valence nitrogen
    let mut mol = parse_smiles("[N]([C])([C])([C])([C])").unwrap();
    mol.explicit_hs();
    mol.infer_hs().unwrap();
    assert_eq!(mol.hydrogens[&node_index(0)], 1);
    assert_eq!(mol.hydrogens[&node_index(1)], 3);
    assert_eq!(mol.hydrogens[&node_index(2)], 3);
    assert_eq!(mol.hydrogens[&node_index(3)], 3);
    assert_eq!(mol.hydrogens[&node_index(4)], 3);
}

// Using amino acid backbone to check edge switching function
#[test]
fn edge_switching_test() {
    // decalaring the backbone structure
    let mut backbone = StableGraph::with_capacity(5, 5);
    backbone.add_node(Atom::C);
    backbone.add_node(Atom::N);
    backbone.add_node(Atom::C);
    backbone.add_node(Atom::O);
    backbone.add_node(Atom::O);

    backbone.add_edge(node_index(0), node_index(1), Single);
    backbone.add_edge(node_index(0), node_index(2), Single);
    backbone.add_edge(node_index(2), node_index(3), Single);
    backbone.add_edge(node_index(2), node_index(4), Double);

    let hydrogen_hash = [
        (node_index(0), 1),
        (node_index(1), 2),
        (node_index(2), 0),
        (node_index(3), 1),
        (node_index(4), 0),
    ]
    .iter()
    .cloned()
    .collect();

    let charge_hash = (0..backbone.node_count())
        .map(|x| (node_index(x), 0))
        .collect();

    let mut mol = Mol {
        graph: backbone,
        hydrogens: hydrogen_hash,
        charges: charge_hash,
        aromatic_atoms: HashSet::default(),
    };

    for _ in 0..10 {
        let prev_graph = mol.graph.clone();
        mol.edge_switch(u32::MAX);
        assert_eq!(tarjan_scc(&mol.graph).len(), 1);
        assert_eq!(
            is_isomorphic_matching::<&Graph<_, _, _, _>, &Graph<_, _, _, _>, _, _>(
                &prev_graph.into(),
                &mol.graph.clone().into(),
                |x, y| x == y,
                |x, y| x == y,
            ),
            false
        );
    }
}

macro_rules! mol_file_roundtrip {
    ($molfile:literal, $checkstr:literal) => {
        let og_str = std::fs::read_to_string($molfile).unwrap();
        let og_mol = match parse_mol_file(&og_str) {
            Ok(x) => x,
            Err(x) => panic!("{}", x),
        };
        let og_str = Itertools::intersperse(og_str.lines().skip(3), "\n").collect::<String>();
        let mut s = Vec::<u8>::new();
        og_mol.to_mol_file(&mut s).expect("Failed to mol file");
        let s = String::from_utf8(s).unwrap();

        let new_mol = match parse_mol_file(&s) {
            Ok(x) => x,
            Err(x) => panic!("{}", x),
        };
        assert!(is_isomorphic_matching::<
            &Graph<_, _, _, _>,
            &Graph<_, _, _, _>,
            _,
            _,
        >(
            &og_mol
                .graph
                .map(
                    |idx, x| (x, og_mol.h_count(&idx), og_mol.charge(&idx)),
                    |_, x| x
                )
                .into(),
            &new_mol
                .graph
                .map(
                    |idx, x| (x, new_mol.h_count(&idx), new_mol.charge(&idx)),
                    |_, x| x
                )
                .into(),
            |x, y| x == y,
            |x, y| x == y,
        ));

        if $checkstr {
            let s = Itertools::intersperse(s.lines().skip(3), "\n").collect::<String>();
            println!("{}", &og_str);
            println!("{}", &s);
            assert_eq!(og_str, s);
        }
    };
}

#[test]
fn to_mol_file() {
    mol_file_roundtrip!("test_files/acetone.mol", false);
    mol_file_roundtrip!("test_files/acetone.mol.implicit", true);
    mol_file_roundtrip!("test_files/ammonium.mol", false);
    mol_file_roundtrip!("test_files/ammonium.mol.implicit", true);
    mol_file_roundtrip!("test_files/aspirin.mol", false);
    mol_file_roundtrip!("test_files/aspirin.mol.implicit", true);
    mol_file_roundtrip!("test_files/benzaldehyde.mol", false);
    mol_file_roundtrip!("test_files/benzaldehyde.mol.implicit", true);
    mol_file_roundtrip!("test_files/benzene.mol", false);
    mol_file_roundtrip!("test_files/benzene.mol.implicit", true);
    mol_file_roundtrip!("test_files/choline.mol", false);
    mol_file_roundtrip!("test_files/choline.mol.implicit", false);
    mol_file_roundtrip!("test_files/egfr.mol", false);
    mol_file_roundtrip!("test_files/egfr.mol.implicit", true);
    mol_file_roundtrip!("test_files/ex1.mol", false);
    mol_file_roundtrip!("test_files/ex1.mol.implicit", true);
    mol_file_roundtrip!("test_files/ex2.mol", false);
    mol_file_roundtrip!("test_files/ex2.mol.implicit", true);
    mol_file_roundtrip!("test_files/ex3.mol", false);
    mol_file_roundtrip!("test_files/ex3.mol.implicit", true);
    // FIXME: this fails because our mol parser is bad
    // mol_file_roundtrip!("test_files/negative_charge.mol", true);
}

macro_rules! round_trip_smiles {
    ($smiles:literal) => {
        let mut mol1 = parse_smiles($smiles).unwrap();
        println!("Original: {}", &($smiles));
        let smiles = mol1.to_smiles();
        println!("From Graph: {}\n", &smiles);
        let mut mol2 = parse_smiles(&smiles).unwrap();

        assert!(mol1.infer_hs().is_ok());
        assert!(mol2.infer_hs().is_ok());

        let node_mapper1 = |node_idx: NodeIndex, node: &Atom| {
            (
                node.atomic_num(),
                mol1.charge(&node_idx),
                mol1.h_count(&node_idx),
                mol1.aromaticity(&node_idx),
            )
        };
        let node_mapper2 = |node_idx: NodeIndex, node: &Atom| {
            (
                node.atomic_num(),
                mol2.charge(&node_idx),
                mol2.h_count(&node_idx),
                mol2.aromaticity(&node_idx),
            )
        };

        let edge_mapper = |_, edge| edge;
        let matches_all = is_isomorphic_matching::<&Graph<_, _, _, _>, &Graph<_, _, _, _>, _, _>(
            &mol1.graph.map(node_mapper1, edge_mapper).into(),
            &mol2.graph.map(node_mapper2, edge_mapper).into(),
            |x, y| x == y,
            |x, y| x == y,
        );

        assert!(matches_all);
    };
}

#[test]
fn to_smiles() {
    round_trip_smiles!("CC(=O)OCOC(=O)C");
    round_trip_smiles!("C([N+2](CCO)(CCO)CCO)C[N+](CCO)(CCO)CCO");
    round_trip_smiles!("C1=CC=CC=C1");
    round_trip_smiles!("O1CCCCC1N1CCCCC1");
    round_trip_smiles!("CC1=CC(CCC1)Br");
    round_trip_smiles!("C1=CC=C(C(=C1)C=C(C#N)C#N)O");
    round_trip_smiles!("C1(CCCCC1)P(C1CCCCC1)C1CCCCC1");
    round_trip_smiles!("C(CCl)C(F)(F)F");
    round_trip_smiles!("CCCCCCCCCCCCCCCCSCCCCCCCCCCCCCCCC");
    round_trip_smiles!("C(CCl)SC(F)(F)F");
    round_trip_smiles!("CSC(C(F)Cl)(F)F");
    round_trip_smiles!("C(C(C(F)(F)F)(Br)Br)Cl");
    round_trip_smiles!("C(Cl)(Cl)(Cl)Cl");
    round_trip_smiles!("C1CC2CC1CC2Cl");
    round_trip_smiles!("CCC(C)C(C(=O)NC(CC(=O)NC)C(=O)NC(C(C)C)C(=O)NC(C(C(=O)NC)O)C(=O)NC(C)C(=O)NC(CC(=O)NC)C(=O)NC(C(C)C)C(=O)NC(CO)C(=O)NC(C(C)C)C(=O)NC(CC(=O)N)C(=O)NC(C(=O)NC(CC(=O)N)C(=O)NC(CCC(=O)N)C(=O)NC(C(C)O)C(=O)NC(C(C)O)C(=O)O)C(C)(C)CS(=O)C)NC(=O)C(CC(=O)NC)NC(=O)CNC(=O)C(C(C)(C)O)NC(=O)C(C(C)(C)C)NC(=O)C(C(C(=O)NC)O)NC(=O)C(C(C)CC)NC(=O)C(CC(=O)NC)NC(=O)CNC(=O)CNC(=O)C(C)NC(=O)C(C(C)(C)O)NC(=O)C(C(C)CC(=O)N)NC(=O)C(CC(=O)NC)NC(=O)C(C(C)(C)C)NC(=O)CNC(=O)C(C)NC(=O)CNC(=O)C(C(C)(C)O)NC(=O)C(CC(=O)NC)NC(=O)C(C)NC(=O)C(C(C)(C)C)NC(=O)C(C)NC(=O)CNC(=O)C(C)NC(=O)C(C(C)(C)C)NC(=O)C(C(C)(C)C)NC(=O)C(C)NC(=O)C(C(C)(C)C)NC(=O)C(C(C)(C)C)NC(=O)C(C(C)(C)C)NC(=O)CNC(=O)C(C(C)(C)CC)NC(=O)CNC(=O)C(=O)CCC(C)(C)C");
    round_trip_smiles!("CC1C2=NC(=CS2)C(=O)NC(=C)C(=O)NC(=C)C3=NC(=C(O3)C)C4=C(C=CC(=N4)C(=O)NC(=C)C(=O)NC(=C)C(=O)N)C5=NC(=CS5)C(=O)NC(C(=O)NCC6=NC(=C(O6)C)C(=O)NC(C7=NC(=CS7)C(=O)N1)C)C(C)C");
    round_trip_smiles!("CC([C]=O)[NH]");
    round_trip_smiles!("C[N+](C)(C)CCO.[OH-]");
    round_trip_smiles!("O=C(C(CC1=CC=CC=C1)NC(C2N=C(C3=COC(C4=COC(C5=COC(C6=COC(C(C(CC)C)NC(C(C(CC)C)NC(C7=C(C)OC(C8=C(C)OC(C9=CSC(C%10=C(C)OC(C%11=CSC(C(CCCNC(N)=N)N(C)C)=N%11)=N%10)=N9)=N8)=N7)=O)=O)=N6)=N5)=N4)=N3)OC2C)=O)O");
}

#[test]
#[cfg(feature = "inchi")]
fn to_smiles_inchi() {
    use crate::MolBuilder;

    macro_rules! smiles_inchi {
        ($inchi:literal, $smiles:literal) => {{
            // implicit
            let mol = MolBuilder::from_smiles($smiles)
                .unwrap()
                .implicit_hs()
                .build();
            let new_smiles = mol.to_smiles();
            let round_trip = MolBuilder::from_smiles(&new_smiles)
                .unwrap_or_else(|e| {
                    panic!(
                        "failed to parse generated smiles {} from {}. Error {}",
                        new_smiles, $smiles, e
                    )
                })
                .implicit_hs()
                .build();
            assert_eq!(
                $inchi,
                mol.to_inchi_key().unwrap(),
                "failed getting {} from {}",
                $inchi,
                $smiles
            );
            assert_eq!(
                $inchi,
                round_trip.to_inchi_key().unwrap(),
                "failed round tripping {} from {}. Generated {}",
                $inchi,
                $smiles,
                new_smiles
            );

            // unimplicit
            let mol = MolBuilder::from_smiles($smiles).unwrap().build();
            let new_smiles = mol.to_smiles();
            let round_trip = MolBuilder::from_smiles(&new_smiles)
                .unwrap_or_else(|e| {
                    panic!(
                        "failed to parse generated smiles {} from {}. Error {}",
                        new_smiles, $smiles, e
                    )
                })
                .build();
            assert_eq!(
                $inchi,
                mol.to_inchi_key().unwrap(),
                "failed getting {} from {}",
                $inchi,
                $smiles
            );
            assert_eq!(
                $inchi,
                round_trip.to_inchi_key().unwrap(),
                "failed round tripping {} from {}. Generated {}",
                $inchi,
                $smiles,
                new_smiles
            );
        }};
    }

    // these are all from pubchem
    smiles_inchi!(
        "DXOVFKCXTQMPEM",
        "C1=CC(=C(C=C1[As](=O)(O)[O-])NCS(=O)[O-])O.[Na+].[Na+]"
    );
    smiles_inchi!(
        "QJJXOEFWXSQISU",
        "CN(C)CCN1C(=O)CC(SC2=CC=CC=C21)C3=CC=CC=C3"
    );
    smiles_inchi!("WROWUHCAAFNMJL", "CNC(=O)OC1=CC=CC=C1[N+](=O)[O-]");
    smiles_inchi!(
        "BBFQZRXNYIEMAW",
        "COC1=CC=CC2=C3C(=C(C=C21)[N+](=O)[O-])C(=CC4=C3OCO4)C(=O)O"
    );
    smiles_inchi!(
        "ZLVMAMIPILWYHQ",
        "CCOC(=O)OC1=C(C=C(C=C1)CCNC(=O)C(CCSC)NC(=O)C)OC(=O)OCC"
    );
    smiles_inchi!("FHIDNBAQOFJWCA", "C1=C(C(=O)NC(=O)N1C2C(C(C(O2)CO)O)O)F");
    smiles_inchi!("HFPZCAJZSCWRBC", "CC1=CC=C(C=C1)C(C)C");
    smiles_inchi!(
        "KDFKJOFJHSVROC",
        "COC1=C(C2=C(CC3C4=CC(=C(C=C4CCN3C2)OC)O)C=C1)OC"
    );
    smiles_inchi!("WKNFADCGOAHBPG", "CN1CCC(CC1)C2=CNC3=C2C=C(C=C3)O");
    smiles_inchi!(
        "IHDMJEPWSOHAAL",
        "C[NH+](C)CCOC(=O)C(C1=CC=CC=C1)C2=CC=CC=C2.[Cl-]"
    );
    smiles_inchi!( "OIMACDRJUANHTJ",  "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)OCC4C(C(C(O4)N5C=NC6=C(N=CN=C65)N)O)O)O)O)N");
    smiles_inchi!("KKIMDKMETPPURN", "C1CN(CCN1)C2=CC=CC(=C2)C(F)(F)F");
    smiles_inchi!(
        "KTHVBAZBLKXIHZ",
        "CCN1CCCC(C1)OC(=O)C(C2=CC=CC=C2)C3=CC=CC=C3"
    );
    smiles_inchi!("DKPFZGUDAPQIHT", "CCCCOC(=O)C");
    smiles_inchi!(
        "QDLHCMPXEPAAMD",
        "CC(=O)OC1CC2(C(CCC2=O)C3=C1C4(C(OC(=O)C5=COC(=C54)C3=O)COC)C)C"
    );
    smiles_inchi!("ZEZJPIDPVXJEME", "C1=CNC(=CC1=O)O");
    smiles_inchi!("NKAAEMMYHLFEFN", "C(C(C(=O)[O-])O)(C(=O)O)O.[Na+]");
    smiles_inchi!("UVRVNMDNGVIBGM", "CC1CCC2=C1C=C(C=C2)C");
    smiles_inchi!("IISBACLAFKSPIT", "CC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O");
    smiles_inchi!("BPLQKQKXWHCZSS", "COC1=CC(=CC(=C1OC)OC)CC=C");
    smiles_inchi!("DVSPHWCZXKPJEQ", "C[N+](C)(C)CCOC");
    smiles_inchi!("CQMYCPZZIPXILQ", "CCOC(=O)CCCCCCCC(=O)OCC");
    smiles_inchi!("QZWNXXINFABALM", "C1C2CC3CC1CC(C2)C3N");
    smiles_inchi!("KTDILAZSFSAPMD", "C1=CC=C(C=C1)NCCCl");
    smiles_inchi!("QWLIHSCHURTGJX", "C[Se]C1=CC=CC=C1C(=O)O");
    smiles_inchi!("OJOWICOBYCXEKR", "CC=C1CC2CC1C=C2");
    smiles_inchi!("XPPSKHVGJRBAAP", "CN(C)C(=O)C=CC1=CC=C(O1)[N+](=O)[O-]");
    smiles_inchi!("ZCOGQSHZVSZAHH", "CN(C)C(=O)N1CC1");
    smiles_inchi!(
        "MUMTUZWXHQLLMC",
        "C[N+](C)(C)CC#CC(C1=CC=CC=C1)(C2=CC=CC=N2)O"
    );
    smiles_inchi!("LXCFILQKKLGQFO", "COC(=O)C1=CC=C(C=C1)O");

    // these are manual
    smiles_inchi!("HZVOZRGWRWCICA", "[CH2]");
    smiles_inchi!("AHKZTVQIVOEVFO", "[O-2]");
    smiles_inchi!("SKALCVOFYPVXLA", "O=C(C(CC1=CC=CC=C1)NC(C2N=C(C3=COC(C4=COC(C5=COC(C6=COC(C(C(CC)C)NC(C(C(CC)C)NC(C7=C(C)OC(C8=C(C)OC(C9=CSC(C%10=C(C)OC(C%11=CSC(C(CCCNC(N)=N)N(C)C)=N%11)=N%10)=N9)=N8)=N7)=O)=O)=N6)=N5)=N4)=N3)OC2C)=O)O")
}

#[test]
fn peptide_bonds() {
    // AHACG
    let mol =
        MolBuilder::from_smiles("CC(C(NC(CC=1NC=NC=1)C(NC(C)C(NC(CS)C(NCC(O)=O)=O)=O)=O)=O)N")
            .unwrap()
            .build();
    assert_eq!(
        mol.peptide_bonds().count(),
        4,
        "AHACG wrong number of peptide bonds"
    );
    assert_eq!(
        mol.generalized_peptide_bonds()
            .collect::<HashSet<_>>()
            .len(),
        4,
        "AHACG wronger number of generalized peptide bonds"
    );

    // GAHAPEP
    let mol = MolBuilder::from_smiles("CC(C(NC(CC=3NC=NC=3)C(NC(C)C(N2C(CCC2)C(NC(CCC(=O)O)C(N1C(CCC1)C(O)=O)=O)=O)=O)=O)=O)NC(=O)CN").unwrap().build();
    assert_eq!(
        mol.peptide_bonds().count(),
        6,
        "GAHAPEP wrong number of peptide bonds"
    );
    assert_eq!(
        mol.generalized_peptide_bonds()
            .collect::<HashSet<_>>()
            .len(),
        6,
        "GAHAPEP wronger number of generalized peptide bonds"
    );

    // PPPPP
    let mol = MolBuilder::from_smiles(
        "N5C(CCC5)C(N4C(CCC4)C(N3C(CCC3)C(N2C(CCC2)C(N1C(CCC1)C(O)=O)=O)=O)=O)=O",
    )
    .unwrap()
    .build();
    assert_eq!(
        mol.peptide_bonds().count(),
        4,
        "PPPPP wrong number of peptide bonds"
    );
    assert_eq!(
        mol.generalized_peptide_bonds()
            .collect::<HashSet<_>>()
            .len(),
        4,
        "PPPPP wronger number of generalized peptide bonds"
    );

    // plantazolicin https://pubchem.ncbi.nlm.nih.gov/compound/Plantazolicin
    let mol = MolBuilder::from_smiles("CCC(C)C(C1=NC(=CO1)C2=NC(=CO2)C3=NC(=CO3)C4=NC(=CO4)C5=NC(C(O5)C)C(=O)NC(CC6=CC=CC=C6)C(=O)O)NC(=O)C(C(C)CC)NC(=O)C7=C(OC(=N7)C8=C(OC(=N8)C9=CSC(=N9)C1=C(OC(=N1)C1=CSC(=N1)C(CCCN=C(N)N)N(C)C)C)C)C").unwrap().build();
    assert_eq!(
        mol.peptide_bonds().count(),
        3,
        "plantazolicin wrong number of peptide bonds"
    );
    assert_eq!(
        mol.generalized_peptide_bonds()
            .collect::<HashSet<_>>()
            .len(),
        13,
        "plantazolicin wronger number of generalized peptide bonds"
    );

    // SapB https://pubchem.ncbi.nlm.nih.gov/compound/44227772
    let mol = MolBuilder::from_smiles("CCC(C)C1C(=O)NC(C(=O)NC(C(=O)NC(CSCC(C(=O)NC(C(=O)NC(C(=O)NC(=C)C(=O)N1)CC(C)C)CO)NC(=O)C(CC(=O)O)NC(=O)CNC(=O)C2CSCC(C(=O)NC(C(=O)NC(C(=O)NC(=C)C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N2)CC(C)C)CC(C)C)CC(C)C)C)CCCNC(=N)N)NC(=O)CNC(=O)C(C(C)O)N)C(=O)NC(CC(=O)N)C(=O)O)C(C)O)C(C)O").unwrap().build();
    assert_eq!(
        mol.peptide_bonds().count(),
        20,
        "SapB wrong number of peptide bonds"
    );
    assert_eq!(
        mol.generalized_peptide_bonds()
            .collect::<HashSet<_>>()
            .len(),
        20,
        "SapB wronger number of generalized peptide bonds"
    );

    // actinomycin X2
    let mol = MolBuilder::from_smiles("CC1C(C(=O)NC(C(=O)N2CCCC2C(=O)N(CC(=O)N(C(C(=O)O1)C(C)C)C)C)C(C)C)NC(=O)C3=C4C(=C(C=C3)C)OC5=C(C(=O)C(=C(C5=N4)C(=O)NC6C(OC(=O)C(N(C(=O)CN(C(=O)C7CC(=O)CN7C(=O)C(NC6=O)C(C)C)C)C)C(C)C)C)N)C").unwrap().build();
    assert_eq!(
        mol.peptide_bonds().count(),
        10,
        "actinomycin X2 wrong number of peptide bonds"
    );
    assert_eq!(
        mol.generalized_peptide_bonds()
            .collect::<HashSet<_>>()
            .len(),
        10,
        "actinomycin X2 wronger number of generalized peptide bonds"
    );
}
