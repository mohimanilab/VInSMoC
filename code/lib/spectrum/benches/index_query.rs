use std::str::FromStr;

use criterion::{criterion_group, criterion_main, Criterion};
use paths::{FilePath, NewDirPath};
use spectrum::{
    index::disk::DiskIndex, index::Ingest, ParsedCollections, Peak, SpectrumPreprocessing,
};
use tempfile::tempdir;

fn bench_query(c: &mut Criterion) {
    let tmp_dir = tempdir().unwrap();
    let tmp_dir_path = tmp_dir.path().to_path_buf().join("test");
    let mut index = DiskIndex::new(
        NewDirPath::try_from(tmp_dir_path.clone()).unwrap(),
        2.0,
        0.02,
    )
    .unwrap();
    let colls = ParsedCollections::new(
        vec![FilePath::from_str("test_files/B3326_R9.mzXML.multiple_charged.mgf").unwrap()],
        SpectrumPreprocessing::default(),
    );
    index.ingest_spectra(colls).unwrap();
    let queries = [
        // picked from the mgf file used to create the index
        (324.18, Peak::from((154.721, 0.5))),
        (589.25, Peak::from((89.966, 0.5))),
        (409.21, Peak::from((113.07, 0.5))),
        (692.27, Peak::from((232.22, 0.5))),
        (301.14, Peak::from((109.88, 0.5))),
        // picked from other mgf files not in the index
        (297.17, Peak::from((155.11, 0.5))),
        (508.33, Peak::from((91.06, 0.5))),
        (463.30, Peak::from((85.10, 0.5))),
        (590.42, Peak::from((169.08, 0.5))),
        (485.30, Peak::from((226.28, 0.5))),
    ];

    c.bench_function("query_10", |b| {
        b.iter_with_large_drop(|| {
            let mut _s = 0;
            for q in queries.iter() {
                _s += index
                    .query(q.1.mz, q.0, None)
                    .unwrap()
                    .map(|s| s.spectrum_id)
                    .sum::<u32>();
            }
        });
    });
    tmp_dir.close().unwrap();
}

criterion_group!(querying, bench_query);
criterion_main!(querying);
