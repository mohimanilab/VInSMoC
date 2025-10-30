use std::str::FromStr;

use criterion::{criterion_group, criterion_main, Criterion};
use paths::{FilePath, NewDirPath};
use spectrum::{index::disk::DiskIndex, index::Ingest, ParsedCollections, SpectrumPreprocessing};
use tempfile::tempdir;

fn bench_ingest(c: &mut Criterion) {
    let tmp_dir = tempdir().unwrap();
    let tmp_dir_path = tmp_dir.path().to_path_buf();
    let mut index = DiskIndex::new(
        NewDirPath::try_from(tmp_dir_path.clone()).unwrap(),
        2.0,
        0.02,
    )
    .unwrap();

    c.bench_function("ingest", |b| {
        b.iter_with_large_drop(|| {
            let colls = ParsedCollections::new(
                vec![FilePath::from_str("test_files/B3326_R9.mzXML.multiple_charged.mgf").unwrap()],
                SpectrumPreprocessing::default(),
            );
            index.ingest_spectra(colls).unwrap();
        });
    });
    tmp_dir.close().unwrap();
}

criterion_group!(ingesting, bench_ingest);
criterion_main!(ingesting);
