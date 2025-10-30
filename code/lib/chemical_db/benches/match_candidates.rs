use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};

use chemical_db::{CandidateFilter, ChemicalDb};
use itertools::Itertools;
use molecule::Adduct;
use ordered_float::OrderedFloat;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use rustc_hash::FxHashMap;
use spectrum::{Ms, SpectrumCollection};

fn old_match_candidates<'a, D: AsRef<ChemicalDb>>(
    mgf: &SpectrumCollection,
    chem_db: D,
    filter: &'a CandidateFilter,
) -> FxHashMap<usize, Vec<(&'a Adduct, usize)>> {
    let chem_db = chem_db.as_ref();

    mgf.spectra
        .iter()
        .flat_map(|spectrum| {
            filter
                .adducts
                .iter()
                .filter(move |adduct| {
                    !filter.respect_charges || spectrum.charges().contains(&adduct.charge())
                })
                .flat_map(move |adduct| {
                    let mass = adduct.to_mass(spectrum.pepmass());
                    match chem_db.query(mass, filter.error_tolerance) {
                        None => vec![],
                        Some(chem_db_indices) => chem_db_indices
                            .map(|chem_db_ind| (chem_db_ind, adduct, spectrum.scan))
                            .collect(),
                    }
                    .into_iter()
                })
        })
        .sorted_by(|x, y| x.0.cmp(&y.0))
        .group_by(|(chem_db_ind, _, _)| *chem_db_ind)
        .into_iter()
        .map(|(chem_db_ind, group)| {
            (
                chem_db_ind,
                group
                    .map(|(_, adduct, scan)| (adduct, scan))
                    .collect::<Vec<(&Adduct, usize)>>(),
            )
        })
        .collect()
}

fn match_candidates<'a, D: AsRef<ChemicalDb>>(
    mgf: &SpectrumCollection,
    chem_db: &D,
    filter: &'a CandidateFilter,
) -> Vec<(usize, Vec<(&'a Adduct, usize)>)> {
    let masses = mgf.precursors(&filter.adducts);
    chem_db
        .as_ref()
        .iter()
        .enumerate()
        .filter_map(|(row_idx, row)| {
            let lower = OrderedFloat::<f64>::from(row.exact_mass() - filter.error_tolerance);
            let upper = OrderedFloat::<f64>::from(row.exact_mass() + filter.error_tolerance);
            let candidates = masses
                .range(lower..upper)
                .flat_map(|(_, v)| v.iter().copied())
                .collect::<Vec<_>>();
            (!candidates.is_empty()).then(|| (row_idx, candidates))
        })
        .collect()
}

fn par_match_candidates<'a, D: AsRef<ChemicalDb>>(
    mgf: &SpectrumCollection,
    chem_db: &D,
    filter: &'a CandidateFilter,
) -> Vec<(usize, Vec<(&'a Adduct, usize)>)> {
    let masses = mgf.precursors(&filter.adducts);
    chem_db
        .as_ref()
        .par_iter()
        .enumerate()
        .filter_map(|(row_idx, row)| {
            let lower = OrderedFloat::<f64>::from(row.exact_mass() - filter.error_tolerance);
            let upper = OrderedFloat::<f64>::from(row.exact_mass() + filter.error_tolerance);
            let candidates = masses
                .range(lower..upper)
                .flat_map(|(_, v)| v.iter().copied())
                .collect::<Vec<_>>();
            (!candidates.is_empty()).then(|| (row_idx, candidates))
        })
        .collect()
}

fn small_old(c: &mut Criterion) {
    let mut group = c.benchmark_group("old_match_small");

    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.csv").unwrap();
    let spectra = SpectrumCollection::from_path("test_files/mgf_example.mgf").unwrap();
    let filter = CandidateFilter::default();
    group.bench_function(BenchmarkId::new("serial", "unsorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "unsorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });

    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "sorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "sorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });

    group.bench_function(BenchmarkId::new("serial", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
}

fn small_new(c: &mut Criterion) {
    let mut group = c.benchmark_group("new_match_small");

    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.csv").unwrap();
    let spectra = SpectrumCollection::from_path("test_files/mgf_example.mgf").unwrap();
    let filter = CandidateFilter::default();
    group.bench_function(BenchmarkId::new("serial", "unsorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });

    let db = ChemicalDb::from_path_unchecked("test_files/test_chem_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
}

fn large_old(c: &mut Criterion) {
    let mut group = c.benchmark_group("old_match_large");

    let db = ChemicalDb::from_path_unchecked("test_files/med_db.csv").unwrap();
    let spectra = SpectrumCollection::from_path("test_files/large.mgf").unwrap();
    let filter = CandidateFilter::default();
    group.bench_function(BenchmarkId::new("serial", "unsorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "unsorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });

    let db = ChemicalDb::from_path_unchecked("test_files/med_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "sorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "sorted_resort"), |b| {
        b.iter_batched(
            || db.clone(),
            |mut db| {
                db.par_sort_rows();
                old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter))
            },
            BatchSize::SmallInput,
        )
    });

    group.bench_function(BenchmarkId::new("serial", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| old_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
}

fn large_new(c: &mut Criterion) {
    let mut group = c.benchmark_group("new_match_large");

    let db = ChemicalDb::from_path_unchecked("test_files/med_db.csv").unwrap();
    let spectra = SpectrumCollection::from_path("test_files/large.mgf").unwrap();
    let filter = CandidateFilter::default();
    group.bench_function(BenchmarkId::new("serial", "unsorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "unsorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| par_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });

    let db = ChemicalDb::from_path_unchecked("test_files/med_db.sorted.csv").unwrap();
    group.bench_function(BenchmarkId::new("serial", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function(BenchmarkId::new("parallel", "sorted"), |b| {
        b.iter_batched(
            || db.clone(),
            |db| par_match_candidates(black_box(&spectra), black_box(&db), black_box(&filter)),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(matching, small_old, small_new, large_old, large_new);
criterion_main!(matching);
