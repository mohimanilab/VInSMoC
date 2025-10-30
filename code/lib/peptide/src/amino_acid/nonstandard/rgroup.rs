use molecule::BondType::{Double, Single};
use molecule::{Atom, Mol, MolBuilder, MolGraph};
// use petgraph::graph::{node_index, Graph};
use petgraph::stable_graph::node_index;

pub fn aad() -> Mol {
    let mut rgroup = MolGraph::with_capacity(6, 5);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Double);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(4), 1),
    ]
    .iter()
    .cloned()
    .collect();

    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn abu() -> Mol {
    let mut rgroup = MolGraph::with_capacity(2, 1);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 3)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn bht() -> Mol {
    let mut rgroup = MolGraph::with_capacity(9, 9);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(0), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Double);
    rgroup.add_edge(node_index(2), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Single);
    rgroup.add_edge(node_index(4), node_index(6), Double);
    rgroup.add_edge(node_index(5), node_index(7), Double);
    rgroup.add_edge(node_index(6), node_index(7), Single);
    rgroup.add_edge(node_index(7), node_index(8), Single);

    let hydrogen_hash = [
        (node_index(0), 1),
        (node_index(1), 1),
        (node_index(3), 1),
        (node_index(4), 1),
        (node_index(5), 1),
        (node_index(6), 1),
        (node_index(8), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn dpg() -> Mol {
    let mut rgroup = MolGraph::with_capacity(8, 8);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Double);
    rgroup.add_edge(node_index(0), node_index(2), Single);
    rgroup.add_edge(node_index(1), node_index(3), Single);
    rgroup.add_edge(node_index(2), node_index(4), Double);
    rgroup.add_edge(node_index(3), node_index(5), Double);
    rgroup.add_edge(node_index(4), node_index(5), Single);
    rgroup.add_edge(node_index(3), node_index(6), Single);
    rgroup.add_edge(node_index(4), node_index(7), Single);

    let hydrogen_hash = [
        (node_index(1), 1),
        (node_index(2), 1),
        (node_index(5), 1),
        (node_index(6), 1),
        (node_index(7), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn hpg() -> Mol {
    let mut rgroup = MolGraph::with_capacity(7, 7);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(0), node_index(2), Double);
    rgroup.add_edge(node_index(1), node_index(3), Double);
    rgroup.add_edge(node_index(2), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Single);
    rgroup.add_edge(node_index(4), node_index(5), Double);
    rgroup.add_edge(node_index(5), node_index(6), Single);

    let hydrogen_hash = [
        (node_index(1), 1),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(4), 1),
        (node_index(6), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn orn() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(3), 2),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn phg() -> Mol {
    let mut rgroup = MolGraph::with_capacity(6, 6);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Double);
    rgroup.add_edge(node_index(0), node_index(2), Single);
    rgroup.add_edge(node_index(1), node_index(3), Single);
    rgroup.add_edge(node_index(2), node_index(4), Double);
    rgroup.add_edge(node_index(3), node_index(5), Double);
    rgroup.add_edge(node_index(4), node_index(5), Single);

    let hydrogen_hash = [
        (node_index(1), 1),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(4), 1),
        (node_index(5), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn pip() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(3), 2),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn cit() -> Mol {
    let mut rgroup = MolGraph::with_capacity(7, 6);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(4), node_index(5), Single);
    rgroup.add_edge(node_index(4), node_index(6), Double);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(3), 1),
        (node_index(5), 2),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn dab() -> Mol {
    let mut rgroup = MolGraph::with_capacity(3, 2);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 2), (node_index(2), 2)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn hty() -> Mol {
    let mut rgroup = MolGraph::with_capacity(9, 9);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Double);
    rgroup.add_edge(node_index(2), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Single);
    rgroup.add_edge(node_index(4), node_index(6), Double);
    rgroup.add_edge(node_index(5), node_index(7), Double);
    rgroup.add_edge(node_index(6), node_index(7), Single);
    rgroup.add_edge(node_index(7), node_index(8), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(3), 1),
        (node_index(4), 1),
        (node_index(5), 1),
        (node_index(6), 1),
        (node_index(8), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn hse() -> Mol {
    let mut rgroup = MolGraph::with_capacity(3, 2);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 2), (node_index(2), 1)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn hrn() -> Mol {
    let mut rgroup = MolGraph::with_capacity(5, 4);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(3), 1),
        (node_index(4), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn han() -> Mol {
    let mut rgroup = MolGraph::with_capacity(8, 7);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Single);
    rgroup.add_edge(node_index(5), node_index(6), Single);
    rgroup.add_edge(node_index(5), node_index(7), Double);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(4), 1),
        (node_index(6), 3),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn dap() -> Mol {
    let mut rgroup = MolGraph::with_capacity(3, 2);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 2)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn bnd() -> Mol {
    let mut rgroup = MolGraph::with_capacity(8, 8);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(0), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(2), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Single);
    rgroup.add_edge(node_index(4), node_index(6), Single);
    rgroup.add_edge(node_index(5), node_index(6), Single);
    rgroup.add_edge(node_index(6), node_index(7), Double);

    let hydrogen_hash = [
        (node_index(0), 1),
        (node_index(1), 1),
        (node_index(2), 1),
        (node_index(3), 2),
        (node_index(4), 1),
        (node_index(5), 1),
        (node_index(7), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn adh() -> Mol {
    let mut rgroup = MolGraph::with_capacity(6, 5);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(0), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Double);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Single);

    let hydrogen_hash = [
        (node_index(0), 1),
        (node_index(1), 3),
        (node_index(2), 1),
        (node_index(4), 3),
        (node_index(5), 3),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn dmw() -> Mol {
    let mut rgroup = MolGraph::with_capacity(15, 16);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Double);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(4), node_index(5), Double);
    rgroup.add_edge(node_index(4), node_index(9), Single);
    rgroup.add_edge(node_index(5), node_index(6), Single);
    rgroup.add_edge(node_index(6), node_index(7), Double);
    rgroup.add_edge(node_index(7), node_index(8), Single);
    rgroup.add_edge(node_index(8), node_index(9), Double);
    rgroup.add_edge(node_index(9), node_index(1), Single);
    rgroup.add_edge(node_index(3), node_index(10), Single);
    rgroup.add_edge(node_index(10), node_index(11), Single);
    rgroup.add_edge(node_index(10), node_index(12), Single);
    rgroup.add_edge(node_index(10), node_index(13), Single);
    rgroup.add_edge(node_index(13), node_index(14), Double);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(2), 1),
        (node_index(5), 1),
        (node_index(6), 1),
        (node_index(7), 1),
        (node_index(8), 1),
        (node_index(11), 3),
        (node_index(12), 3),
        (node_index(13), 1),
        (node_index(14), 2),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn dhb() -> Mol {
    let mut rgroup = MolGraph::with_capacity(2, 1);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);

    let hydrogen_hash = [(node_index(0), 1), (node_index(1), 3)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn crp() -> Mol {
    let mut rgroup = MolGraph::with_capacity(11, 11);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::Cl);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Double);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(4), node_index(5), Double);
    rgroup.add_edge(node_index(4), node_index(9), Single);
    rgroup.add_edge(node_index(5), node_index(6), Single);
    rgroup.add_edge(node_index(6), node_index(7), Double);
    rgroup.add_edge(node_index(7), node_index(8), Single);
    rgroup.add_edge(node_index(8), node_index(9), Double);
    rgroup.add_edge(node_index(9), node_index(1), Single);
    rgroup.add_edge(node_index(7), node_index(10), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(5), 1),
        (node_index(6), 1),
        (node_index(8), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn mha() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(2), node_index(4), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 1),
        (node_index(3), 3),
        (node_index(4), 3),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn hfn() -> Mol {
    let mut rgroup = MolGraph::with_capacity(8, 7);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(3), node_index(5), Single);
    rgroup.add_edge(node_index(5), node_index(6), Double);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(4), 1),
        (node_index(5), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn mpr() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(1), node_index(3), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 1),
        (node_index(2), 2),
        (node_index(3), 3),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn bla() -> Mol {
    let mut rgroup = MolGraph::with_capacity(1, 0);
    rgroup.add_node(Atom::C);

    // leave one hydrogen out since this is a beta amino acid
    let hydrogen_hash = [(node_index(0), 2)].iter().cloned().collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn bbu() -> Mol {
    let mut rgroup = MolGraph::with_capacity(2, 1);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);

    // add an extra hydrogen for the upcoming beta shift
    let hydrogen_hash = [(node_index(0), 1), (node_index(1), 3)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn ade() -> Mol {
    let mut rgroup = MolGraph::with_capacity(8, 7);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(4), node_index(5), Single);
    rgroup.add_edge(node_index(5), node_index(6), Single);
    rgroup.add_edge(node_index(6), node_index(7), Single);

    // add an extra hydrogen for the upcoming beta shift
    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(3), 2),
        (node_index(4), 2),
        (node_index(5), 2),
        (node_index(6), 2),
        (node_index(7), 3),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn pyr() -> Mol {
    let mut rgroup = MolGraph::with_capacity(7, 7);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Double);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(4), node_index(5), Double);
    rgroup.add_edge(node_index(5), node_index(6), Single);
    rgroup.add_edge(node_index(6), node_index(1), Double);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(4), 1),
        (node_index(6), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn cgy() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Double);
    rgroup.add_edge(node_index(2), node_index(3), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 1),
        (node_index(2), 1),
        (node_index(3), 3),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn omy() -> Mol {
    let mut rgroup = MolGraph::with_capacity(9, 9);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Double);
    rgroup.add_edge(node_index(3), node_index(4), Single);
    rgroup.add_edge(node_index(4), node_index(5), Single);
    rgroup.add_edge(node_index(4), node_index(6), Double);
    rgroup.add_edge(node_index(6), node_index(7), Single);
    rgroup.add_edge(node_index(7), node_index(1), Double);
    rgroup.add_edge(node_index(5), node_index(8), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(6), 1),
        (node_index(7), 1),
        (node_index(8), 3),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn piz() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(3), 1),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}
