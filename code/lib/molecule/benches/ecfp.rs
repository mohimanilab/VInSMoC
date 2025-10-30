use std::time::Duration;

use criterion::{
    black_box, criterion_group, criterion_main, measurement::WallTime, BenchmarkId, Criterion,
};

use molecule::{Mol, MolBuilder};

fn bench_ecfp(group: &mut criterion::BenchmarkGroup<WallTime>, mol: &Mol) {
    for radius in [4, 6, 8] {
        let id = BenchmarkId::new("ecfp", radius);
        group.bench_with_input(id, &radius, |b, &radius| {
            b.iter(|| mol.ecfp(black_box(radius)))
        });
    }
}

fn bench_hash(group: &mut criterion::BenchmarkGroup<WallTime>, mol: &Mol) {
    for radius in [4, 6, 8] {
        let id = BenchmarkId::new("hash", radius);
        group.bench_with_input(id, &radius, |b, &radius| {
            b.iter(|| mol.ecfp_hash(black_box(radius)))
        });
    }
}

fn simple(c: &mut Criterion) {
    let mut group = c.benchmark_group("simple");
    group
        .warm_up_time(Duration::from_millis(250))
        .measurement_time(Duration::from_secs(2))
        .noise_threshold(0.15);
    let smiles = "CCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    bench_ecfp(&mut group, &mol);
    bench_hash(&mut group, &mol);
}

fn one_ring(c: &mut Criterion) {
    let mut group = c.benchmark_group("one_ring");
    group
        .warm_up_time(Duration::from_millis(250))
        .measurement_time(Duration::from_secs(2))
        .noise_threshold(0.15);
    let smiles = "CCC(C)CCCCCCC1CC(O)=NC(CCC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(O)=NC(C(C)C)C(O)=NC(CC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(=O)O1";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    bench_ecfp(&mut group, &mol);
    bench_hash(&mut group, &mol);
}

fn multi_ring(c: &mut Criterion) {
    let mut group = c.benchmark_group("multi_ring");
    group
        .warm_up_time(Duration::from_millis(250))
        .measurement_time(Duration::from_secs(2))
        .noise_threshold(0.15);
    let smiles = "CC(=O)OC1CC(OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CC(O)C5(C)C(C7=CC(=O)OC7)C(O)CC65O)C4)OC3C)OC2C)OC(C)C1OC1OC(CO)C(O)C(O)C1O";
    let mol = MolBuilder::from_smiles(smiles).unwrap().build();
    bench_ecfp(&mut group, &mol);
    bench_hash(&mut group, &mol);
}

criterion_group!(ecfp, simple, one_ring, multi_ring);
criterion_main!(ecfp);
