use molecule::MolBuilder;
use petgraph::{algo::is_isomorphic_matching, Graph};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let smiles1 = args.next().unwrap();
    let smiles2 = args.next().unwrap();

    let mol1 = MolBuilder::from_smiles(&smiles1)?.implicit_hs().build();
    let mol2 = MolBuilder::from_smiles(&smiles2)?.implicit_hs().build();

    println!(
        "{}",
        is_isomorphic_matching(
            &Graph::from(mol1.graph.map(|idx, x| (x, mol1.h_count(&idx)), |_, x| x)),
            &Graph::from(mol2.graph.map(|idx, x| (x, mol2.h_count(&idx)), |_, x| x)),
            |x, y| x == y,
            |x, y| x == y,
        )
    );

    Ok(())
}
