use iai::black_box;

use molecule::MolBuilder;

fn simple_r4() {
    let smiles = "CCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.ecfp(black_box(4));
}

fn simple_r6() {
    let smiles = "CCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.ecfp(black_box(6));
}

fn simple_r8() {
    let smiles = "CCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.ecfp(black_box(8));
}

fn complex_r4() {
    let smiles = "CC(=O)OC1CC(OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CC(O)C5(C)C(C7=CC(=O)OC7)C(O)CC65O)C4)OC3C)OC2C)OC(C)C1OC1OC(CO)C(O)C(O)C1O";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.ecfp(black_box(4));
}

fn complex_r6() {
    let smiles = "CC(=O)OC1CC(OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CC(O)C5(C)C(C7=CC(=O)OC7)C(O)CC65O)C4)OC3C)OC2C)OC(C)C1OC1OC(CO)C(O)C(O)C1O";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.ecfp(black_box(6));
}

fn complex_r8() {
    let smiles = "CC(=O)OC1CC(OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CC(O)C5(C)C(C7=CC(=O)OC7)C(O)CC65O)C4)OC3C)OC2C)OC(C)C1OC1OC(CO)C(O)C(O)C1O";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    mol.ecfp(black_box(8));
}

iai::main!(simple_r4, simple_r6, simple_r8, complex_r4, complex_r6, complex_r8);
