use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};

use chemical_db::ChemicalDb;

fn small(c: &mut Criterion) {
    let mut group = c.benchmark_group("sort_small_chem_db");

    // serial impls
    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
            },
            BatchSize::SmallInput,
        )
    });
    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "unsorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
            },
            BatchSize::SmallInput,
        )
    });

    // parallel impls
    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("parallel", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
            },
            BatchSize::SmallInput,
        )
    });
    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.csv").unwrap();
    group.bench_function(BenchmarkId::new("parallel", "unsorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
            },
            BatchSize::SmallInput,
        )
    });
}

fn large(c: &mut Criterion) {
    let mut group = c.benchmark_group("sort_large_chem_db");

    // serial impls
    let db = ChemicalDb::from_path_unchecked("test_files/large_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
            },
            BatchSize::SmallInput,
        )
    });
    let db = ChemicalDb::from_path_unchecked("test_files/large_db.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "unsorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
            },
            BatchSize::SmallInput,
        )
    });

    // parallel impls
    let db = ChemicalDb::from_path_unchecked("test_files/large_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("parallel", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
            },
            BatchSize::SmallInput,
        )
    });
    let db = ChemicalDb::from_path_unchecked("test_files/large_db.csv").unwrap();
    group.bench_function(BenchmarkId::new("parallel", "unsorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
            },
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(sorting, small, large);
criterion_main!(sorting);
