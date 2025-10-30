use criterion::{black_box, criterion_group, criterion_main, Criterion};

use molecule::MolBuilder;

fn simple(c: &mut Criterion) {
    let smiles = black_box("CCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC");
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    c.bench_function("simple", |b| b.iter(|| mol.to_smiles()));
}

fn one_ring(c: &mut Criterion) {
    let smiles = black_box("CCC(C)CCCCCCC1CC(O)=NC(CCC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(O)=NC(C(C)C)C(O)=NC(CC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(=O)O1");
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    c.bench_function("one_ring", |b| b.iter(|| mol.to_smiles()));
}

fn multi_ring(c: &mut Criterion) {
    let smiles = black_box("CC(=O)OC1CC(OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CC(O)C5(C)C(C7=CC(=O)OC7)C(O)CC65O)C4)OC3C)OC2C)OC(C)C1OC1OC(CO)C(O)C(O)C1O");
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    c.bench_function("multi_ring", |b| b.iter(|| mol.to_smiles()));
}

criterion_group!(to_smiles, simple, one_ring, multi_ring);
criterion_main!(to_smiles);
