use std::str::FromStr;

use criterion::{criterion_group, criterion_main, Criterion};
use paths::{FilePath, NewDirPath};
use spectrum::{
    index::disk::DiskIndex, index::Ingest, index::Score, mgf::Mgf, MsFormat, ParsedCollections,
    SpectrumCollection, SpectrumPreprocessing,
};
use tempfile::tempdir;

fn bench_score(c: &mut Criterion) {
    let tmp_dir = tempdir().unwrap();
    let tmp_dir_path = tmp_dir.path().to_path_buf().join("test");
    let mut index = DiskIndex::new(
        NewDirPath::try_from(tmp_dir_path.clone()).unwrap(),
        2.0,
        0.02,
    )
    .unwrap();

    let query_mgf = Mgf::from_path("test_files/B3326_R9.mzXML.multiple_charged.mgf").unwrap();
    let query_spectra = SpectrumCollection::from(query_mgf);

    let colls = ParsedCollections::new(
        vec![FilePath::from_str("test_files/B3326_R9.mzXML.multiple_charged.mgf").unwrap()],
        SpectrumPreprocessing::default(),
    );
    index.ingest_spectra(colls).unwrap();

    c.bench_function("score_vec", |b| {
        b.iter_with_large_drop(|| {
            index
                .score(query_spectra.spectra.iter().next().unwrap())
                .unwrap();
        });
    });
    tmp_dir.close().unwrap();
}

criterion_group!(scoring, bench_score);
criterion_main!(scoring);
