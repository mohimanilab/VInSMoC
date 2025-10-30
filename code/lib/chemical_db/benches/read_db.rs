use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use chemical_db::ChemicalDb;

fn read(c: &mut Criterion) {
    let mut group = c.benchmark_group("read_chem_db");

    group.bench_function(BenchmarkId::new("unchecked", "sorted"), |b| {
        b.iter_with_large_drop(|| {
            ChemicalDb::from_path_unchecked("test_files/test_chem_db.sorted.csv").unwrap();
        })
    });
    group.bench_function(BenchmarkId::new("unchecked", "unsorted"), |b| {
        b.iter_with_large_drop(|| {
            ChemicalDb::from_path_unchecked("test_files/test_chem_db.csv").unwrap();
        })
    });
    group.bench_function(BenchmarkId::new("unchecked", "duplicated"), |b| {
        b.iter_with_large_drop(|| {
            ChemicalDb::from_path_unchecked("test_files/test_chem_db.dup.csv").unwrap();
        })
    });

    group.bench_function(BenchmarkId::new("checked", "sorted"), |b| {
        b.iter_with_large_drop(|| {
            ChemicalDb::from_path("test_files/test_chem_db.sorted.csv").unwrap();
        })
    });
    group.bench_function(BenchmarkId::new("checked", "unsorted"), |b| {
        b.iter_with_large_drop(|| {
            ChemicalDb::from_path("test_files/test_chem_db.csv").unwrap();
        })
    });
    group.bench_function(BenchmarkId::new("checked", "duplicated"), |b| {
        b.iter_with_large_drop(|| {
            ChemicalDb::from_path("test_files/test_chem_db.dup.csv").unwrap();
        })
    });
}

criterion_group!(reading, read);
criterion_main!(reading);
