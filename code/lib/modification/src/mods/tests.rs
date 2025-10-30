use super::*;
use molecule::BondType::*;
use molecule::{Atom, MolBuilder, MolGraph};

use petgraph::{algo::is_isomorphic_matching, prelude::EdgeIndex, stable_graph::node_index, Graph};

fn test_with_amino_acid_h(modif: Mod, r_group: Mol) {
    let mut amino = MolGraph::with_capacity(6, 6);
    amino.add_node(Atom::C);
    amino.add_node(Atom::C);
    amino.add_node(Atom::C);
    amino.add_node(Atom::N);
    amino.add_node(Atom::C);
    amino.add_node(Atom::N);

    amino.add_edge(node_index(0), node_index(1), Single);
    amino.add_edge(node_index(1), node_index(2), Aromatic);
    amino.add_edge(node_index(2), node_index(3), Aromatic);
    amino.add_edge(node_index(3), node_index(4), Aromatic);
    amino.add_edge(node_index(4), node_index(5), Aromatic);
    amino.add_edge(node_index(5), node_index(1), Aromatic);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 0),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(4), 1),
        (node_index(5), 0),
    ]
    .iter()
    .cloned()
    .collect();

    let mut amino = MolBuilder::default()
        .graph(amino)
        .hydrogens(hydrogen_hash)
        .build();

    assert!(
        f64::abs(modif.mass_shift(&amino) - (r_group.exact_mass() - amino.exact_mass())) < 0.000001
    );

    assert!(modif.modify(&mut amino).is_ok());

    assert_eq!(amino.hydrogens().clone(), r_group.hydrogens().clone());
    assert_eq!(amino.charges().clone(), r_group.charges().clone());
    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &amino.graph.into(),
        &r_group.graph.into(),
        |x, y| x == y,
        |x, y| x == y
    ));
}

#[test]
fn modify_works() {
    let remove_node = Mod::RemoveNode(5);
    let add_edge = Mod::AddEdge(0, 3, BondType::Single);
    let h_shift = Mod::HydrogenShift(4, -1);

    // test remove_node
    let mut rgroup = MolGraph::with_capacity(6, 6);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Aromatic);
    rgroup.add_edge(node_index(2), node_index(3), Aromatic);
    rgroup.add_edge(node_index(3), node_index(4), Aromatic);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 0),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(4), 1),
    ]
    .iter()
    .cloned()
    .collect();

    let h_node_remove = MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build();

    test_with_amino_acid_h(remove_node, h_node_remove);

    // test add_edge
    let mut rgroup = MolGraph::with_capacity(6, 6);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Aromatic);
    rgroup.add_edge(node_index(2), node_index(3), Aromatic);
    rgroup.add_edge(node_index(3), node_index(4), Aromatic);
    rgroup.add_edge(node_index(4), node_index(5), Aromatic);
    rgroup.add_edge(node_index(5), node_index(1), Aromatic);
    rgroup.add_edge(node_index(0), node_index(3), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 0),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(4), 1),
        (node_index(5), 0),
    ]
    .iter()
    .cloned()
    .collect();

    let h_edge_add = MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build();

    test_with_amino_acid_h(add_edge, h_edge_add);

    // test h_shift
    let mut rgroup = MolGraph::with_capacity(6, 6);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Aromatic);
    rgroup.add_edge(node_index(2), node_index(3), Aromatic);
    rgroup.add_edge(node_index(3), node_index(4), Aromatic);
    rgroup.add_edge(node_index(4), node_index(5), Aromatic);
    rgroup.add_edge(node_index(5), node_index(1), Aromatic);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 0),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(4), 0),
        (node_index(5), 0),
    ]
    .iter()
    .cloned()
    .collect();
    let h_h_shift = MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build();

    test_with_amino_acid_h(h_shift, h_h_shift);
}

#[test]
fn reaction_remove_edge() {
    // test removing an edge
    let reactant = MolBuilder::from_smiles("CCCC").unwrap().build();
    let mut product = MolBuilder::from_smiles("CCCC").unwrap().build();
    product.graph.remove_edge(EdgeIndex::new(1));

    let changes = reaction(&reactant, &product, "name".into());
    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![Mod::RemoveEdge(1, 2)]
    );
    assert_eq!(changes.pattern.graph.node_count(), 2);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(1)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(2)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(1)), 2);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(2)), 2);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(1)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(2)), 0);
    assert_eq!(&changes.name, "name");
}

#[test]
fn reaction_add_edge() {
    // test adding an edge and shifting hydrogen count
    let reactant = MolBuilder::from_smiles("CCCCCC").unwrap().build();
    let product = MolBuilder::from_smiles("C1CCCCC1").unwrap().build();

    let changes = reaction(&reactant, &product, "another_name".into());

    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![
            Mod::HydrogenShift(0, -1),
            Mod::HydrogenShift(5, -1),
            Mod::AddEdge(0, 5, Single),
        ]
    );
    assert_eq!(changes.pattern.graph.node_count(), 2);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(0)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(5)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(0)), 3);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(5)), 3);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(0)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(5)), 0);
    assert_eq!(&changes.name, "another_name");
}

#[test]
fn reaction_remove_node() {
    // test removing a node and shifting hydrogen count
    let reactant = MolBuilder::from_smiles("CCCC").unwrap().build();
    let product = MolBuilder::from_smiles("CC").unwrap().build();

    let changes = reaction(&reactant, &product, "remove_node".into());
    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![
            Mod::HydrogenShift(1, 1),
            Mod::RemoveNode(2),
            Mod::RemoveNode(3),
        ]
    );
    assert_eq!(changes.pattern.graph.node_count(), 3);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(1)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(2)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(3)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(1)), 2);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(2)), 2);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(3)), 3);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(1)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(2)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(3)), 0);
    assert_eq!(&changes.name, "remove_node");
}

#[test]
fn reaction_change_charge() {
    // test changing a node's hydrogen count and charge
    let reactant = MolBuilder::from_smiles("CCCC").unwrap().build();
    let product = MolBuilder::from_smiles("CCC[C+]").unwrap().build();

    let changes = reaction(&reactant, &product, "change charge".into());

    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![Mod::HydrogenShift(3, -3), Mod::ChargeShift(3, 1),]
    );
    assert_eq!(changes.pattern.graph.node_count(), 1);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(3)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(3)), 3);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(3)), 0);
    assert_eq!(&changes.name, "change charge");
}

#[test]
fn reaction_add_node() {
    // test adding a node
    let mut reactant = MolBuilder::from_smiles("CCC").unwrap().build();
    let product = MolBuilder::from_smiles("CCCC").unwrap().build();

    let changes = reaction(&reactant, &product, "add node".into());

    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![
            Mod::HydrogenShift(2, -1),
            Mod::AddNode(Atom::C, 3, 0, 3),
            Mod::AddEdge(2, 3, Single)
        ]
    );
    assert_eq!(changes.pattern.graph.node_count(), 1);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(2)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(2)), 3);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(2)), 0);
    assert_eq!(&changes.name, "add node");

    let result = changes.mod_seq.modify(&mut reactant);
    assert!(result.is_ok());
    assert!(reactant
        .graph
        .contains_edge(NodeIndex::new(2), NodeIndex::new(3)));
}

#[test]
fn reaction_complex() {
    // test with complex smiles
    let reactant = MolBuilder::from_smiles("COC1=C(C=C(C=C1)C2CC3=C(C(=CC=C3)O)C(=O)O2)O")
        .unwrap()
        .build();
    let mut product = reactant.clone();
    let mods = ModSeq::new(vec![
        Mod::HydrogenShift(0, -1),
        Mod::HydrogenShift(20, -1),
        Mod::AddEdge(0, 20, Single),
    ]);
    let result = mods.modify(&mut product);
    assert!(result.unwrap() == ());
    let changes = reaction(&reactant, &product, "complex".into());

    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![
            Mod::HydrogenShift(0, -1),
            Mod::HydrogenShift(20, -1),
            Mod::AddEdge(0, 20, Single)
        ]
    );
    assert_eq!(changes.pattern.graph.node_count(), 2);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(0)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(20)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(0)), 3);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(20)), 1);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(0)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(20)), 0);
    assert_eq!(&changes.name, "complex");

    let mut product = reactant.clone();
    let mods = ModSeq::new(vec![
        Mod::RemoveEdge(6, 7),
        Mod::RemoveEdge(2, 7),
        Mod::AddEdge(2, 6, Single),
    ]);
    let result = mods.modify(&mut product);
    assert!(result.unwrap() == ());
    let changes = reaction(&reactant, &product, "complex2".into());

    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![
            Mod::AddEdge(2, 6, Single),
            Mod::RemoveEdge(6, 7),
            Mod::RemoveEdge(2, 7)
        ]
    );
    assert_eq!(changes.pattern.graph.node_count(), 3);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(2)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(6)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(7)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(2)), 0);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(6)), 1);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(7)), 1);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(2)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(6)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(7)), 0);
    assert_eq!(&changes.name, "complex2");

    let reactant = MolBuilder::from_smiles("CNCC(C1=CC(=C(C=C1)O)O)O")
        .unwrap()
        .build();
    let mut product = reactant.clone();
    let mods = ModSeq::new(vec![
        Mod::AddEdge(4, 5, Single),
        Mod::AddEdge(4, 12, Single),
        Mod::ChargeShift(5, 1),
    ]);
    let result = mods.modify(&mut product);
    assert!(result.unwrap() == ());
    let changes = reaction(&reactant, &product, "complex3".into());

    assert_eq!(
        changes.mod_seq.iter().map(|x| *x).collect::<Vec<Mod>>(),
        vec![
            Mod::ChargeShift(5, 1),
            Mod::AddEdge(4, 5, Single),
            Mod::AddEdge(4, 12, Single)
        ]
    );
    assert_eq!(changes.pattern.graph.node_count(), 3);
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(4)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(5)));
    assert!(changes.pattern.graph.contains_node(NodeIndex::new(12)));
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(4)), 0);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(5)), 1);
    assert_eq!(changes.pattern.h_count(&NodeIndex::new(12)), 1);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(4)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(5)), 0);
    assert_eq!(changes.pattern.charge(&NodeIndex::new(12)), 0);
    assert_eq!(&changes.name, "complex3");
}

#[test]
fn reaction_center_add_edge() {
    // simple test
    let mut reactant = MolBuilder::from_smiles("CCCC").unwrap().build();
    reactant.graph.remove_edge(EdgeIndex::new(1));

    let center = reaction_center(&reactant, &ModSeq::new(vec![Mod::AddEdge(1, 2, Single)]), 0);
    let correct_center = vec![NodeIndex::new(1), NodeIndex::new(2)];
    let mut found = 0;
    for idx in center {
        for val in &correct_center {
            if *val == idx {
                found += 1;
                break;
            }
        }
    }
    assert!(found == correct_center.len());
}

#[test]
fn reaction_center_radius() {
    // check that radius works
    let mut reactant = MolBuilder::from_smiles("CCCCCC").unwrap().build();
    reactant.graph.remove_edge(EdgeIndex::new(1));

    let center = reaction_center(&reactant, &ModSeq::new(vec![Mod::AddEdge(1, 2, Single)]), 1);
    let correct_center = vec![
        NodeIndex::new(0),
        NodeIndex::new(1),
        NodeIndex::new(2),
        NodeIndex::new(3),
    ];
    let mut found = 0;
    for idx in center {
        for val in &correct_center {
            if *val == idx {
                found += 1;
                break;
            }
        }
    }
    assert!(found == correct_center.len());

    let center = reaction_center(&reactant, &ModSeq::new(vec![Mod::AddEdge(1, 2, Single)]), 2);
    let correct_center = vec![
        NodeIndex::new(0),
        NodeIndex::new(1),
        NodeIndex::new(2),
        NodeIndex::new(3),
        NodeIndex::new(4),
    ];
    let mut found = 0;
    for idx in center {
        for val in &correct_center {
            if *val == idx {
                found += 1;
                break;
            }
        }
    }
    assert!(found == correct_center.len());
}

#[test]
fn reaction_center_h_charge() {
    // check that differences in hydrogen and charge are noted as
    // being part of a reaction
    let reactant = MolBuilder::from_smiles("CCCC").unwrap().build();

    let center = reaction_center(
        &reactant,
        &ModSeq::new(vec![
            Mod::HydrogenShift(0, -1),
            Mod::ChargeShift(1, -1),
            Mod::ChargeShift(2, -1),
        ]),
        0,
    );
    let correct_center = vec![NodeIndex::new(0), NodeIndex::new(1), NodeIndex::new(2)];
    let mut found = 0;
    for idx in center {
        for val in &correct_center {
            if *val == idx {
                found += 1;
                break;
            }
        }
    }
    assert!(found == correct_center.len());
}

#[test]
fn reaction_center_remove_edge() {
    // check that edge removal is seen as part of a reaction correctly
    let reactant = MolBuilder::from_smiles("CCCC").unwrap().build();

    let center = reaction_center(&reactant, &ModSeq::new(vec![Mod::RemoveEdge(1, 2)]), 0);
    let correct_center = vec![NodeIndex::new(1), NodeIndex::new(2)];
    let mut found = 0;
    for idx in center {
        for val in &correct_center {
            if *val == idx {
                found += 1;
                break;
            }
        }
    }
    assert!(found == correct_center.len());
}

#[test]
fn reaction_center_change_bond() {
    // check a change in bond type is seen as part of a reaction correctly
    let reactant = MolBuilder::from_smiles("CC=CC").unwrap().build();

    let center = reaction_center(&reactant, &ModSeq::new(vec![Mod::AddEdge(1, 2, Single)]), 0);
    let correct_center = vec![NodeIndex::new(1), NodeIndex::new(2)];
    let mut found = 0;
    for idx in center {
        for val in &correct_center {
            if *val == idx {
                found += 1;
                break;
            }
        }
    }
    assert!(found == correct_center.len());
}

#[test]
fn display_mod() {
    assert_eq!(
        Mod::AddEdge(1, 2, BondType::Single).to_string().as_str(),
        "AddEdge 1 2 -"
    );
    assert_eq!(
        Mod::AddEdge(1, 2, BondType::Double).to_string().as_str(),
        "AddEdge 1 2 ="
    );
    assert_eq!(Mod::RemoveEdge(1, 2).to_string().as_str(), "RemoveEdge 1 2");
    assert_eq!(
        Mod::AddNode(Atom::H, 1, 0, 1).to_string().as_str(),
        "AddNode H 1 0 1"
    );
    assert_eq!(
        Mod::AddNode(Atom::C, 1, -2, 5).to_string().as_str(),
        "AddNode C 1 -2 5"
    );
    assert_eq!(Mod::RemoveNode(1).to_string().as_str(), "RemoveNode 1");
    assert_eq!(Mod::RemoveNode(2).to_string().as_str(), "RemoveNode 2");
    assert_eq!(Mod::RemoveCC(1).to_string().as_str(), "RemoveCC 1");
    assert_eq!(Mod::RemoveCC(2).to_string().as_str(), "RemoveCC 2");
    assert_eq!(
        Mod::ChargeShift(2, 1).to_string().as_str(),
        "ChargeShift 2 1"
    );
    assert_eq!(
        Mod::ChargeShift(100, -1).to_string().as_str(),
        "ChargeShift 100 -1"
    );
    assert_eq!(
        Mod::HydrogenShift(2, 1).to_string().as_str(),
        "HydrogenShift 2 1"
    );
    assert_eq!(
        Mod::HydrogenShift(100, -1).to_string().as_str(),
        "HydrogenShift 100 -1"
    );
}
