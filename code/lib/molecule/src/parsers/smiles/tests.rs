use super::*;

#[derive(Eq, PartialEq, Clone, Debug)]
struct Bond {
    pub start_atom: usize,
    pub end_atom: usize,
    pub bond_type: BondType,
}

impl Bond {
    fn new(start_atom: usize, end_atom: usize, bond_type: BondType) -> Self {
        Bond {
            start_atom,
            end_atom,
            bond_type,
        }
    }
}

use petgraph::{graph::edge_index, prelude::EdgeIndex};

// We don't test more granular parsers than this
#[test]
fn organic_atom_works() {
    assert_eq!(organic_atom("CC"), Ok(("C", (OrganicAtom::C, false))));
    assert_eq!(organic_atom("Cl"), Ok(("", (OrganicAtom::Cl, false))));
    assert_eq!(organic_atom("Bl"), Ok(("l", (OrganicAtom::B, false))));
    assert_eq!(organic_atom("ClC"), Ok(("C", (OrganicAtom::Cl, false))));
    assert_eq!(organic_atom("c"), Ok(("", (OrganicAtom::C, true))));
    assert_eq!(organic_atom("n"), Ok(("", (OrganicAtom::N, true))));
    assert_eq!(organic_atom("o"), Ok(("", (OrganicAtom::O, true))));
    assert_eq!(organic_atom("s"), Ok(("", (OrganicAtom::S, true))));
    assert!(organic_atom("A").is_err());
    assert!(organic_atom("b").is_err());
    assert!(organic_atom("i").is_err());
    assert!(organic_atom("Y").is_err());
    assert!(organic_atom("H").is_err());
}

#[test]
fn inorganic_atom_works() {
    assert_eq!(inorganic_atom("Au"), Ok(("", InorganicAtom::Au)));
    assert_eq!(inorganic_atom("YYb"), Ok(("Yb", InorganicAtom::Y)));
    assert_eq!(inorganic_atom("SeCCCC"), Ok(("CCCC", InorganicAtom::Se)));
    assert!(inorganic_atom("H").is_err());
}

#[test]
fn hydrogen_count_works() {
    assert_eq!(hydrogen_count("H3"), Ok(("", 3u8)));
    assert_eq!(hydrogen_count("H"), Ok(("", 1)));
    assert_eq!(hydrogen_count(""), Ok(("", 0)));
    assert_eq!(hydrogen_count("H2"), Ok(("", 2)));
}

#[test]
fn charge_works() {
    assert_eq!(charge("+"), Ok(("", 1i8)));
    assert_eq!(charge("+++"), Ok(("", 3i8)));
    assert_eq!(charge("-"), Ok(("", -1i8)));
    assert_eq!(charge("--"), Ok(("", -2i8)));
    assert_eq!(charge("+22"), Ok(("", 22i8)));
    assert_eq!(charge("-16"), Ok(("", -16i8)));
    assert_eq!(charge(""), Ok(("", 0i8)));
}

#[test]
fn chirality_works() {
    assert_eq!(chirality("@"), Ok(("", ())));
    assert_eq!(chirality("@@"), Ok(("", ())));
    assert_eq!(chirality("@SP1"), Ok(("", ())));
    assert_eq!(chirality("@AL1"), Ok(("", ())));
    assert_eq!(chirality("@OH20"), Ok(("", ())));
}

#[test]
fn bracketed_atom_works() {
    assert_eq!(bracketed_atom("[C]"), Ok(("", ((Atom::C, false), 0, 0))));
    assert_eq!(bracketed_atom("[C@@]"), Ok(("", ((Atom::C, false), 0, 0))));
    assert_eq!(
        bracketed_atom("[C]trash"),
        Ok(("trash", ((Atom::C, false), 0, 0)))
    );
    assert_eq!(
        bracketed_atom("[BH2++]"),
        Ok(("", ((Atom::B, false), 2, 2)))
    );
    assert_eq!(
        bracketed_atom("[B@AL2H2++]"),
        Ok(("", ((Atom::B, false), 2, 2)))
    );
    assert_eq!(bracketed_atom("[B-]"), Ok(("", ((Atom::B, false), 0, -1))));
    assert_eq!(bracketed_atom("[F+2]"), Ok(("", ((Atom::F, false), 0, 2))));
    assert_eq!(bracketed_atom("[OH]"), Ok(("", ((Atom::O, false), 1, 0))));
    assert_eq!(
        bracketed_atom("[XeH1++]"),
        Ok(("", ((Atom::Xe, false), 1, 2)))
    );
    assert_eq!(bracketed_atom("[Se]"), Ok(("", ((Atom::Se, false), 0, 0))));
    assert_eq!(bracketed_atom("[H]"), Ok(("", ((Atom::H, false), 0, 0))));
}

#[test]
fn atom_works() {
    assert_eq!(
        atom("C"),
        Ok(("", SmilesSyntaxNode::Atom(Atom::C, None, 0, false)))
    );
    assert_eq!(
        atom("NNC"),
        Ok(("NC", SmilesSyntaxNode::Atom(Atom::N, None, 0, false)))
    );
    assert_eq!(
        atom("[OH-]NNC"),
        Ok(("NNC", SmilesSyntaxNode::Atom(Atom::O, Some(1), -1, false)))
    );
    assert_eq!(
        atom("oNNC"),
        Ok(("NNC", SmilesSyntaxNode::Atom(Atom::O, None, 0, true)))
    );
    assert_eq!(
        atom("[AuH++]CBB"),
        Ok(("CBB", SmilesSyntaxNode::Atom(Atom::Au, Some(1), 2, false)))
    );
}

#[test]
fn bond_type_works() {
    assert_eq!(bond_type("-"), Ok(("", BondType::Single)));
    assert_eq!(bond_type("-=-"), Ok(("=-", BondType::Single)));
    assert_eq!(bond_type("="), Ok(("", BondType::Double)));
    assert_eq!(bond_type("#="), Ok(("=", BondType::Triple)));
    assert_eq!(bond_type(":="), Ok(("=", BondType::Aromatic)));
    assert!(bond_type("CCC").is_err());
}

#[test]
fn bond_works() {
    assert_eq!(
        bond("-"),
        Ok(("", SmilesSyntaxNode::Bond(BondType::Single)))
    );
    assert_eq!(
        bond("-=-"),
        Ok(("=-", SmilesSyntaxNode::Bond(BondType::Single)))
    );
    assert_eq!(
        bond("="),
        Ok(("", SmilesSyntaxNode::Bond(BondType::Double)))
    );
    assert_eq!(
        bond("#="),
        Ok(("=", SmilesSyntaxNode::Bond(BondType::Triple)))
    );
    assert_eq!(
        bond(":="),
        Ok(("=", SmilesSyntaxNode::Bond(BondType::Aromatic)))
    );
}

#[test]
fn ring_closure_digit_works() {
    assert_eq!(ring_closure_digit("1"), Ok(("", 1)));
    assert_eq!(ring_closure_digit("12"), Ok(("2", 1)));
    assert_eq!(ring_closure_digit("%12"), Ok(("", 12)));
    assert_eq!(ring_closure_digit("%121"), Ok(("1", 12)));
    assert_eq!(ring_closure_digit("1%221"), Ok(("%221", 1)));
    assert_eq!(ring_closure_digit("%221"), Ok(("1", 22)));
}

#[test]
fn ring_closure_works() {
    assert_eq!(
        ring_closure("1"),
        Ok(("", SmilesSyntaxNode::RingClosure(1)))
    );
    assert_eq!(
        ring_closure("12"),
        Ok(("2", SmilesSyntaxNode::RingClosure(1)))
    );
    assert_eq!(
        ring_closure("%12"),
        Ok(("", SmilesSyntaxNode::RingClosure(12)))
    );
}

#[test]
fn branched_atom_works() {
    assert_eq!(
        branched_atom("C"),
        Ok(("", vec![SmilesSyntaxNode::Atom(Atom::C, None, 0, false)]))
    );
    assert_eq!(
        branched_atom("CCC"),
        Ok(("CC", vec![SmilesSyntaxNode::Atom(Atom::C, None, 0, false)]))
    );
    assert_eq!(
        branched_atom("C1CC1"),
        Ok((
            "CC1",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(1)
            ]
        ))
    );
    assert_eq!(
        branched_atom("C123"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::RingClosure(2),
                SmilesSyntaxNode::RingClosure(3),
            ]
        ))
    );
    assert_eq!(
        branched_atom("C(CC)"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
            ]
        ))
    );
    assert_eq!(
        branched_atom("C(CC)(NCC)"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::N, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
            ]
        ))
    );
    assert_eq!(
        branched_atom("C(CC(NCC))"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::N, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::BranchEnd,
            ]
        ))
    );

    assert_eq!(
        branched_atom("[Se]"),
        Ok((
            "",
            vec![SmilesSyntaxNode::Atom(Atom::Se, Some(0), 0, false)]
        ))
    );
}

#[test]
fn atom_chain_works() {
    assert_eq!(
        atom_chain("CCC"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
            ]
        ))
    );

    assert_eq!(
        atom_chain("CC(O)C1(CC)C"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
            ]
        ))
    );

    assert_eq!(
        atom_chain("CC(O)C-1(CC)C"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Bond(BondType::Single),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
            ]
        ))
    );

    assert_eq!(
        atom_chain("C=C(=O)C1(CC)C"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
            ]
        ))
    );

    // aspirin from pubchem https://pubchem.ncbi.nlm.nih.gov/compound/2244
    assert_eq!(
        atom_chain("CC(=O)OC1=CC=CC=C1C(=O)O"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
            ]
        ))
    );

    // another from pubchem https://pubchem.ncbi.nlm.nih.gov/compound/687061
    assert_eq!(
        atom_chain("C1C(OC2=CC=CC=C2O1)C(=O)O"),
        Ok((
            "",
            vec![
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(2),
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::RingClosure(2),
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::RingClosure(1),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::C, None, 0, false),
                SmilesSyntaxNode::BranchBegin,
                SmilesSyntaxNode::Bond(BondType::Double),
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
                SmilesSyntaxNode::BranchEnd,
                SmilesSyntaxNode::Atom(Atom::O, None, 0, false),
            ]
        ))
    );

    assert!(atom_chain("=").is_err());
    let c = atom_chain("C[Se]CC=NNc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]");
    assert!(c.is_ok());
    let (leftover, _syntax_tree) = c.unwrap();
    assert_eq!(leftover, "");
}

fn mol_equality_helper(
    mol: Mol,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    charges: Vec<i8>,
    hydrogens: Vec<u8>,
    connect: FxHashMap<usize, FxHashMap<usize, usize>>,
) {
    assert_eq!(
        mol.graph
            .node_indices()
            .map(|x| mol.graph[x])
            .collect::<Vec<Atom>>(),
        atoms
    );
    assert_eq!(
        mol.graph
            .edge_indices()
            .map(|x| mol.graph[x])
            .collect::<Vec<BondType>>(),
        bonds
            .iter()
            .cloned()
            .map(|x| x.bond_type)
            .collect::<Vec<BondType>>()
    );

    let mut found: FxHashSet<EdgeIndex<u16>> = FxHashSet::default();
    for (node, neighbor_map) in connect {
        for (neighbor, bond_idx) in neighbor_map {
            let node_idx = NodeIndex::new(node);
            let neighbor_idx = NodeIndex::new(neighbor);
            let edge_idx = mol.graph.find_edge(node_idx, neighbor_idx);
            assert!(edge_idx.is_some());
            let edge_idx = edge_idx.unwrap();
            assert_eq!(bonds[bond_idx].bond_type, mol.graph[edge_idx]);
            found.insert(edge_idx);
        }
    }

    assert_eq!(found, mol.graph.edge_indices().collect::<FxHashSet<_>>());

    assert_eq!(
        (0..mol.graph.node_count())
            .map(|x| mol.charge(&NodeIndex::new(x)))
            .collect::<Vec<_>>(),
        charges
    );
    assert_eq!(
        (0..mol.graph.node_count())
            .map(|x| mol.h_count(&NodeIndex::new(x)))
            .collect::<Vec<_>>(),
        hydrogens
    );
}

#[test]
fn parse_smiles_works() {
    let mut connect: FxHashMap<usize, FxHashMap<usize, usize>> = FxHashMap::default();
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(1, 0);
    let atoms: Vec<Atom> = vec![Atom::C, Atom::C];
    let bonds: Vec<Bond> = vec![Bond::new(0, 1, BondType::Single)];
    let charges: Vec<i8> = vec![0, 0];
    let hydrogens: Vec<u8> = vec![3, 3];

    let mol: Mol = parse_smiles("CC").unwrap();
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    let mut connect: FxHashMap<usize, FxHashMap<usize, usize>> = FxHashMap::default();
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(1, 0);
    let inner = connect.entry(1).or_insert(FxHashMap::default());
    inner.insert(2, 1);
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(2, 2);
    let atoms = vec![Atom::C, Atom::C, Atom::C];
    let bonds = vec![
        Bond::new(0, 1, BondType::Single),
        Bond::new(1, 2, BondType::Single),
        Bond::new(0, 2, BondType::Single),
    ];
    let charges = vec![0, 0, 0];
    let hydrogens = vec![2, 2, 2];

    let mol: Mol = parse_smiles("C1CC1").unwrap();
    eprintln!("C1CC1");
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    let mut connect: FxHashMap<usize, FxHashMap<usize, usize>> = FxHashMap::default();
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(1, 0);
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(2, 1);
    let inner = connect.entry(2).or_insert(FxHashMap::default());
    inner.insert(3, 2);
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(3, 3);
    let atoms = vec![Atom::C, Atom::O, Atom::C, Atom::C];
    let bonds = vec![
        Bond::new(0, 1, BondType::Double),
        Bond::new(0, 2, BondType::Single),
        Bond::new(2, 3, BondType::Single),
        Bond::new(0, 3, BondType::Single),
    ];
    let charges = vec![0, 0, 0, 0];
    let hydrogens = vec![0, 0, 2, 2];

    let mol: Mol = parse_smiles("C1(=O)CC1").unwrap();
    eprintln!("C1(=O)CC1");
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    let mut connect: FxHashMap<usize, FxHashMap<usize, usize>> = FxHashMap::default();
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(1, 0);
    inner.insert(3, 2);
    inner.insert(4, 4);
    let inner = connect.entry(1).or_insert(FxHashMap::default());
    inner.insert(2, 1);
    let inner = connect.entry(3).or_insert(FxHashMap::default());
    inner.insert(4, 3);
    let mol = parse_smiles("C1(=C(=O))CC1");
    assert!(mol.is_ok());
    let mol = mol.unwrap();
    let atoms = vec![Atom::C, Atom::C, Atom::O, Atom::C, Atom::C];
    let bonds = vec![
        Bond::new(0, 1, BondType::Double),
        Bond::new(1, 2, BondType::Double),
        Bond::new(0, 3, BondType::Single),
        Bond::new(3, 4, BondType::Single),
        Bond::new(0, 4, BondType::Single),
    ];
    let charges = vec![0, 0, 0, 0, 0];
    let hydrogens = vec![0, 0, 0, 2, 2];
    eprintln!("C1(=C(=O))CC1");
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    let mut connect: FxHashMap<usize, FxHashMap<usize, usize>> = FxHashMap::default();
    let inner = connect.entry(0).or_insert(FxHashMap::default());
    inner.insert(1, 0);
    inner.insert(9, 10);
    let inner = connect.entry(1).or_insert(FxHashMap::default());
    inner.insert(2, 1);
    inner.insert(10, 11);
    let inner = connect.entry(2).or_insert(FxHashMap::default());
    inner.insert(3, 2);
    let _inner = connect
        .entry(3)
        .or_insert(FxHashMap::default())
        .insert(4, 3);
    let _inner = connect
        .entry(4)
        .or_insert(FxHashMap::default())
        .insert(5, 4);
    let _inner = connect
        .entry(5)
        .or_insert(FxHashMap::default())
        .insert(6, 5);
    let _inner = connect
        .entry(6)
        .or_insert(FxHashMap::default())
        .insert(7, 6);
    let _inner = connect
        .entry(7)
        .or_insert(FxHashMap::default())
        .insert(8, 7);
    let _inner = connect
        .entry(3)
        .or_insert(FxHashMap::default())
        .insert(8, 8);
    let _inner = connect
        .entry(8)
        .or_insert(FxHashMap::default())
        .insert(9, 9);
    let _inner = connect
        .entry(10)
        .or_insert(FxHashMap::default())
        .insert(11, 12);
    let _inner = connect
        .entry(10)
        .or_insert(FxHashMap::default())
        .insert(12, 13);
    let mol = parse_smiles("C1C(OC2=CC=CC=C2O1)C(=O)O");
    assert!(mol.is_ok());
    let mol = mol.unwrap();
    let atoms = vec![
        Atom::C,
        Atom::C,
        Atom::O,
        Atom::C,
        Atom::C,
        Atom::C,
        Atom::C,
        Atom::C,
        Atom::C,
        Atom::O,
        Atom::C,
        Atom::O,
        Atom::O,
    ];

    let bonds = vec![
        Bond::new(0, 1, BondType::Single),
        Bond::new(1, 2, BondType::Single),
        Bond::new(2, 3, BondType::Single),
        Bond::new(3, 4, BondType::Double),
        Bond::new(4, 5, BondType::Single),
        Bond::new(5, 6, BondType::Double),
        Bond::new(6, 7, BondType::Single),
        Bond::new(7, 8, BondType::Double),
        Bond::new(3, 8, BondType::Single),
        Bond::new(8, 9, BondType::Single),
        Bond::new(0, 9, BondType::Single),
        Bond::new(1, 10, BondType::Single),
        Bond::new(10, 11, BondType::Double),
        Bond::new(10, 12, BondType::Single),
    ];
    let charges = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let hydrogens = vec![2, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1];
    eprintln!("C1C(OC2=CC=CC=C2O1)C(=O)O");
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    // This is for the very simple method of ignoring double bond orientation
    let mol = parse_smiles("N[C@@H](C/)C\\(=O/)O/////");
    assert!(mol.is_ok());
    let mol = mol.unwrap();
    let atoms = vec![
        Atom::N,
        Atom::C,
        Atom::H,
        Atom::C,
        Atom::C,
        Atom::O,
        Atom::O,
    ];
    let bonds = vec![
        Bond::new(0, 1, BondType::Single),
        Bond::new(1, 2, BondType::Single),
        Bond::new(1, 3, BondType::Single),
        Bond::new(1, 4, BondType::Single),
        Bond::new(4, 5, BondType::Double),
        Bond::new(4, 6, BondType::Single),
    ];
    let mut connect: FxHashMap<usize, FxHashMap<usize, usize>> = FxHashMap::default();
    connect
        .entry(0)
        .or_insert(FxHashMap::default())
        .insert(1, 0);
    connect
        .entry(1)
        .or_insert(FxHashMap::default())
        .insert(2, 1);
    connect
        .entry(1)
        .or_insert(FxHashMap::default())
        .insert(3, 2);
    connect
        .entry(1)
        .or_insert(FxHashMap::default())
        .insert(4, 3);
    connect
        .entry(4)
        .or_insert(FxHashMap::default())
        .insert(5, 4);
    connect
        .entry(4)
        .or_insert(FxHashMap::default())
        .insert(6, 5);

    let charges = vec![0, 0, 0, 0, 0, 0, 0];
    let hydrogens = vec![2, 0, 0, 3, 0, 0, 1];
    eprintln!("N[C@@H](C/)C\\(=O/)O/////");
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    let mol = parse_smiles("C.C");
    assert!(mol.is_ok());
    let mol = mol.unwrap();
    let atoms = vec![Atom::C, Atom::C];
    let bonds = Vec::new();
    let connect = FxHashMap::default();
    let hydrogens = vec![4, 4];
    let charges = vec![0, 0];
    eprintln!("C.C");
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    let mol = parse_smiles("[Na+].[O-]c1ccccc1");
    assert!(mol.is_ok());
    let mol = mol.unwrap();
    let atoms = vec![
        Atom::Na,
        Atom::O,
        Atom::C,
        Atom::C,
        Atom::C,
        Atom::C,
        Atom::C,
        Atom::C,
    ];
    let bonds = vec![
        Bond::new(1, 2, BondType::Single),
        Bond::new(2, 3, BondType::Aromatic),
        Bond::new(3, 4, BondType::Aromatic),
        Bond::new(4, 5, BondType::Aromatic),
        Bond::new(5, 6, BondType::Aromatic),
        Bond::new(6, 7, BondType::Aromatic),
        Bond::new(2, 7, BondType::Aromatic),
    ];

    let hydrogens = vec![0, 0, 0, 1, 1, 1, 1, 1];
    let charges = vec![1, -1, 0, 0, 0, 0, 0, 0];
    let mut connect = FxHashMap::default();
    connect
        .entry(1)
        .or_insert(FxHashMap::default())
        .insert(2, 0);
    connect
        .entry(2)
        .or_insert(FxHashMap::default())
        .insert(3, 1);
    connect
        .entry(3)
        .or_insert(FxHashMap::default())
        .insert(4, 2);
    connect
        .entry(4)
        .or_insert(FxHashMap::default())
        .insert(5, 3);
    connect
        .entry(5)
        .or_insert(FxHashMap::default())
        .insert(6, 4);
    connect
        .entry(6)
        .or_insert(FxHashMap::default())
        .insert(7, 5);
    connect
        .entry(2)
        .or_insert(FxHashMap::default())
        .insert(7, 6);

    eprintln!("[Na+].[O-]c1ccccc1");
    mol_equality_helper(mol, atoms, bonds, charges, hydrogens, connect);

    assert!(parse_smiles("CC(=O)c1c(C)cc2c(c1O)C(=O)c1c(OC3OC(CO)C(O)C(O)C3O)cc(O)cc1C2C1c2cc(O)cc(OC3OC(CO)C(O)C(O)C3O)c2C(=O)c2c1cc(C)c(C(C)=O)c2O").is_ok());
}

#[test]
fn explicit_h_parsing() {
    let mol = parse_smiles("N[C@@H](C)C(=O)O");
    assert!(mol.is_ok());
    let mol = mol.unwrap();
    let atoms = vec![
        Atom::N,
        Atom::C,
        Atom::H,
        Atom::C,
        Atom::C,
        Atom::O,
        Atom::O,
    ];

    let bonds = vec![
        Bond::new(0, 1, BondType::Single),
        Bond::new(1, 2, BondType::Single),
        Bond::new(1, 3, BondType::Single),
        Bond::new(1, 4, BondType::Single),
        Bond::new(4, 5, BondType::Double),
        Bond::new(4, 6, BondType::Single),
    ];
    let mut connect: FxHashMap<usize, FxHashMap<usize, usize>> = FxHashMap::default();
    let _inner = connect
        .entry(0)
        .or_insert(FxHashMap::default())
        .insert(1, 0);
    let _inner = connect
        .entry(1)
        .or_insert(FxHashMap::default())
        .insert(2, 1);
    let _inner = connect
        .entry(1)
        .or_insert(FxHashMap::default())
        .insert(3, 2);
    let _inner = connect
        .entry(1)
        .or_insert(FxHashMap::default())
        .insert(4, 3);
    let _inner = connect
        .entry(4)
        .or_insert(FxHashMap::default())
        .insert(5, 4);
    let _inner = connect
        .entry(4)
        .or_insert(FxHashMap::default())
        .insert(6, 5);

    assert_eq!(
        mol.graph
            .node_indices()
            .map(|x| mol.graph[x])
            .collect::<Vec<Atom>>(),
        atoms
    );
    assert_eq!(
        mol.graph
            .edge_indices()
            .map(|x| mol.graph[x])
            .collect::<Vec<BondType>>(),
        bonds
            .iter()
            .cloned()
            .map(|x| x.bond_type)
            .collect::<Vec<BondType>>()
    );
    assert_eq!(
        *(mol.hydrogens().get(&NodeIndex::new(1)).unwrap_or(&0u8)),
        0u8
    );
    assert_eq!(*(mol.hydrogens().get(&NodeIndex::new(3)).unwrap()), 3u8);
    assert_eq!(*(mol.hydrogens().get(&NodeIndex::new(4)).unwrap()), 0u8);
    assert_eq!(*(mol.hydrogens().get(&NodeIndex::new(6)).unwrap()), 1u8);

    for (node, neighbor_map) in connect {
        for (neighbor, bond_idx) in neighbor_map {
            let node_idx = NodeIndex::new(node);
            let neighbor_idx = NodeIndex::new(neighbor);
            assert!(mol.graph.contains_edge(node_idx, neighbor_idx));
            assert_eq!(
                bonds[bond_idx].bond_type,
                mol.graph
                    [EdgeIndex::new(mol.graph.find_edge(node_idx, neighbor_idx).unwrap().index())]
            )
        }
    }
}

#[test]
fn parse_smiles_ring_bonds() {
    let m = parse_smiles("C=1C=CC=CC1").unwrap();
    assert_eq!(m.graph[edge_index(5)], BondType::Double);
    let m = parse_smiles("C=1C=CC=CC=1").unwrap();
    assert_eq!(m.graph[edge_index(5)], BondType::Double);
    let m = parse_smiles("C1C=CC=CC=1").unwrap();
    assert_eq!(m.graph[edge_index(5)], BondType::Double);

    assert_eq!(
        parse_smiles("C#1C=CC=CC=1").unwrap_err(),
        SmilesParseError::InconsistentRingBonds(1)
    );
}

#[test]
#[cfg(feature = "inchi")]
fn gresimycin() {
    // gresimycin with the '\' and '@' marks gives strange results
    let og_smiles = r#"CC(C(N/C(C(N1C(CCC1)C(NC(C)C(NC(C(NC(C(NC(C(NC(C(NC(C(NC(C(NC(C(NC(C(NC(C(N/C(C(NC(C(NC(C(NC(C(NC2C(C)C)=O)CC(C)C)=O)CS/C=C\NC2=O)=O)[C@H](C)CC)=O)=C/C)=O)CO)=O)[H])=O)CCC(N)=O)=O)[C@H](C)CC)=O)C(C)C)=O)CC3=CC=CC=C3)=O)CCC(N)=O)=O)C)=O)C(C)C)=O)=O)=O)=C\C)=O)N(C)C"#;
    let mut mol = parse_smiles(og_smiles).unwrap();
    assert_eq!(mol.graph.node_count(), 132);
    assert_eq!(mol.graph.edge_count(), 134);
    assert!(mol
        .graph
        .edge_indices()
        .map(|e| mol.graph.edge_endpoints(e).unwrap())
        .all(|ns| ns.0.index() < mol.graph.node_count() && ns.1.index() < mol.graph.node_count()));
    assert_eq!(mol.to_inchi_key().unwrap().as_str(), "OIQYGWMODDNMOB");

    mol.implicit_hs();
    assert_eq!(mol.graph.node_count(), 129);
    assert_eq!(mol.graph.edge_count(), 131);
    // inverse of this was failing before but it actually isn't a valid assumption
    // when we delete nodes their indices don't necessarily disappear
    // therefore we expect there to be node indices that are larger than the number of nodes
    //
    // the confusion arose because the .molecule.graph.node_holes index is not in the
    // .molecule.graph.nodes list. Instead this must be taken into account and all indices after
    // those in node_holes must be bumped manually.
    // this was a bad decision for serde representation by the petgraph team
    assert!(!mol
        .graph
        .edge_indices()
        .map(|e| mol.graph.edge_endpoints(e).unwrap())
        .all(|ns| ns.0.index() < mol.graph.node_count() && ns.1.index() < mol.graph.node_count()));
    assert_eq!(mol.to_inchi_key().unwrap().as_str(), "OIQYGWMODDNMOB");
}

#[test]
fn failure_590() {
    // this is to test that we can successfully parse
    // '[C-]#[N+][C@@H]([C@@]1(C)C=C)[C@]([C@@](C[C@H]1Cl)([H])C(C)2C)([H])C3=C(C(C)(C)C=C)NC4=C3C2=CC=C4'
    // from issue 590

    // mini version around the error
    // this is ok
    let smiles = r#"C(CCCl)([H])C(C)C"#;
    assert!(
        parse_smiles(smiles).is_ok(),
        "mini\n{}",
        parse_smiles(smiles).unwrap_err()
    );

    // slightly larger version
    // this is also ok
    let smiles = r#"C1(C)C=CC(C(CC1Cl))"#;
    assert!(
        parse_smiles(smiles).is_ok(),
        "bigger mini\n{}",
        parse_smiles(smiles).unwrap_err()
    );

    // this failed
    let smiles = r#"C(C)C=CC(C(CCCl)C(C)2C)CCC2"#;
    assert!(
        parse_smiles(smiles).is_ok(),
        "ring bond post ()\n{}",
        parse_smiles(smiles).unwrap_err()
    );

    // but this passed
    let smiles = r#"C(C)C=CC(C(CCCl)C2C)CCC2"#;
    assert!(
        parse_smiles(smiles).is_ok(),
        "ring bond no ()\n{}",
        parse_smiles(smiles).unwrap_err()
    );

    // simplifying the problem
    let smiles = r#"C(N)1SCCCO1"#;
    let parsed = parse_smiles(smiles);
    assert!(parsed.is_ok());
    let mol = parsed.unwrap();
    assert_eq!(mol.graph.node_count(), 7);
    assert_eq!(mol.graph.edge_count(), 7);
    assert_eq!(mol.graph[NodeIndex::new(0)], Atom::C);
    assert_eq!(mol.graph[NodeIndex::new(1)], Atom::N);
    assert_eq!(mol.graph[NodeIndex::new(2)], Atom::S);
    assert_eq!(mol.graph[NodeIndex::new(6)], Atom::O);
    assert_eq!(
        mol.graph
            .neighbors(NodeIndex::new(0))
            .collect::<FxHashSet<_>>(),
        [NodeIndex::new(1), NodeIndex::new(6), NodeIndex::new(2)]
            .into_iter()
            .collect::<FxHashSet<_>>()
    );

    // are the [] the problem?
    // no, still failed
    let smiles =
        r#"[C-]#[N+][CH](C1(C)C=C)C(C(CC1Cl)([H])C(C)2C)([H])C3=C(C(C)(C)C=C)NC4=C3C2=CC=C4"#;
    assert!(
        parse_smiles(smiles).is_ok(),
        "less [] {}",
        parse_smiles(smiles).unwrap_err()
    );

    // are the stereomarkers the problem?
    // no, still failed
    let smiles = r#"[C-]#[N+][CH]([C]1(C)C=C)[C]([C](C[CH]1Cl)([H])C(C)2C)([H])C3=C(C(C)(C)C=C)NC4=C3C2=CC=C4"#;
    assert!(
        parse_smiles(smiles).is_ok(),
        "no stereomarkers {}",
        parse_smiles(smiles).unwrap_err()
    );

    // original
    let smiles = r#"[C-]#[N+][C@@H]([C@@]1(C)C=C)[C@]([C@@](C[C@H]1Cl)([H])C(C)2C)([H])C3=C(C(C)(C)C=C)NC4=C3C2=CC=C4"#;
    assert!(
        parse_smiles(smiles).is_ok(),
        "{}",
        parse_smiles(smiles).unwrap_err()
    );

    #[cfg(feature = "inchi")]
    assert_eq!(
        parse_smiles(smiles).unwrap().to_inchi_key().unwrap(),
        "GHYIJWADNLIVDB"
    );
}
