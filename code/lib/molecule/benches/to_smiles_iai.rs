use iai::black_box;

use molecule::MolBuilder;

fn simple() {
    let smiles = black_box("CCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC");
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.to_smiles();
}

fn complex() {
    let smiles = black_box("CC(=O)OC1CC(OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CC(O)C5(C)C(C7=CC(=O)OC7)C(O)CC65O)C4)OC3C)OC2C)OC(C)C1OC1OC(CO)C(O)C(O)C1O");
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.to_smiles();
}

iai::main!(simple, complex);
