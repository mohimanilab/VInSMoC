use molecule::{
    Atom,
    BondType::{Double, Single},
    Mol, MolBuilder, MolGraph,
};
use petgraph::graph::node_index;

pub fn r() -> Mol {
    let mut rgroup = MolGraph::with_capacity(7, 6);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::N);

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

pub fn h() -> Mol {
    let mut rgroup = MolGraph::with_capacity(6, 6);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Double);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Double);
    rgroup.add_edge(node_index(4), node_index(5), Single);
    rgroup.add_edge(node_index(5), node_index(1), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(2), 1),
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

pub fn k() -> Mol {
    let mut rgroup = MolGraph::with_capacity(5, 4);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(3), node_index(4), Single);

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(1), 2),
        (node_index(2), 2),
        (node_index(3), 2),
        (node_index(4), 2),
    ]
    .iter()
    .cloned()
    .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn d() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(1), node_index(3), Double);

    let hydrogen_hash = [(node_index(0), 2), (node_index(2), 1)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn e() -> Mol {
    let mut rgroup = MolGraph::with_capacity(5, 4);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);
    rgroup.add_edge(node_index(2), node_index(4), Double);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 2), (node_index(3), 1)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn s() -> Mol {
    let mut rgroup = MolGraph::with_capacity(2, 1);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);

    rgroup.add_edge(node_index(0), node_index(1), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 1)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn t() -> Mol {
    let mut rgroup = MolGraph::with_capacity(3, 2);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(0), node_index(2), Single);

    let hydrogen_hash = [(node_index(0), 1), (node_index(1), 1), (node_index(2), 3)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn n() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Double);
    rgroup.add_edge(node_index(1), node_index(3), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(3), 2)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn q() -> Mol {
    let mut rgroup = MolGraph::with_capacity(5, 4);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
    rgroup.add_node(Atom::N);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Double);
    rgroup.add_edge(node_index(2), node_index(4), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 2), (node_index(4), 2)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn c() -> Mol {
    let mut rgroup = MolGraph::with_capacity(2, 1);
    rgroup.add_node(Atom::C);

    rgroup.add_node(Atom::S);

    rgroup.add_edge(node_index(0), node_index(1), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 1)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn g() -> Mol {
    Mol::default()
}

pub fn p() -> Mol {
    let mut rgroup = MolGraph::with_capacity(3, 2);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

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

pub fn a() -> Mol {
    let mut rgroup = MolGraph::with_capacity(1, 0);
    rgroup.add_node(Atom::C);

    let hydrogen_hash = [(node_index(0), 3)].iter().cloned().collect();

    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn v() -> Mol {
    let mut rgroup = MolGraph::with_capacity(3, 2);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(0), node_index(2), Single);

    let hydrogen_hash = [(node_index(0), 1), (node_index(1), 3), (node_index(2), 3)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn i() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(0), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);

    let hydrogen_hash = [
        (node_index(0), 1),
        (node_index(1), 3),
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

pub fn l() -> Mol {
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
        (node_index(2), 3),
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

pub fn m() -> Mol {
    let mut rgroup = MolGraph::with_capacity(4, 3);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::S);
    rgroup.add_node(Atom::C);

    rgroup.add_edge(node_index(0), node_index(1), Single);
    rgroup.add_edge(node_index(1), node_index(2), Single);
    rgroup.add_edge(node_index(2), node_index(3), Single);

    let hydrogen_hash = [(node_index(0), 2), (node_index(1), 2), (node_index(3), 3)]
        .iter()
        .cloned()
        .collect();
    MolBuilder::default()
        .graph(rgroup)
        .hydrogens(hydrogen_hash)
        .build()
}

pub fn f() -> Mol {
    let mut rgroup = MolGraph::with_capacity(7, 7);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
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
        (node_index(5), 1),
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

pub fn y() -> Mol {
    let mut rgroup = MolGraph::with_capacity(8, 8);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::C);
    rgroup.add_node(Atom::O);
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

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(2), 1),
        (node_index(3), 1),
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

pub fn w() -> Mol {
    let mut rgroup = MolGraph::with_capacity(10, 10);
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

    let hydrogen_hash = [
        (node_index(0), 2),
        (node_index(2), 1),
        (node_index(3), 1),
        (node_index(5), 1),
        (node_index(6), 1),
        (node_index(7), 1),
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
