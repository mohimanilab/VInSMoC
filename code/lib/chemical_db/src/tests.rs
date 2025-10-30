use std::{fs::File, ops::Range, path::Path};

use clap::Parser;

use super::*;

#[test]
fn new_and_query() {
    let f = File::open("./test_files/test_chem_db.csv").unwrap();
    let chem_db: ChemicalDb = ChemicalDb::from_reader(f).unwrap();

    //checks that the database is sorted in ascending order by exact_mass
    let mut largest: f64 = 0.0;
    for i in 0..100 {
        let entry = chem_db.get(i).unwrap();
        let mass = entry.exact_mass();
        let name = entry.name();
        println!("{}, {}: {}", i, name, mass);
        assert!(mass >= largest);
        largest = mass;
    }

    //Checks that query works correctly
    let r1: Option<Range<usize>> = chem_db.query(1000.0, 100.0);
    assert_eq!(r1, Some(88..91));

    let r2: Option<Range<usize>> = chem_db.query(500.0, 200.0);
    assert_eq!(r2, Some(27..74));

    let r3: Option<Range<usize>> = chem_db.query(650.0, 10.0);
    assert_eq!(r3, Some(68..70));

    let r4: Option<Range<usize>> = chem_db.query(420.0, 5.0);
    assert_eq!(r4, None);

    let r5: Option<Range<usize>> = chem_db.query(1700.0, 10.0);
    assert_eq!(r5, None);

    let r6: Option<Range<usize>> = chem_db.query(10.0, 10.0);
    assert_eq!(r6, None);

    let r7: Option<Range<usize>> = chem_db.query(1700.0, 50.0);
    assert_eq!(r7, Some(99..100));

    let r8: Option<Range<usize>> = chem_db.query(1700.0, 600.0);
    assert_eq!(r8, Some(91..100));

    // check if it can get smalles and largest compounds by mass
    assert_eq!(chem_db.query(151.1, 0.1), Some(0..1));
    assert_eq!(chem_db.query(1665.93, 0.02), Some(99..100));
}

#[test]
fn get_methods() {
    let f = File::open("./test_files/test_chem_db.csv").unwrap();
    let chem_db: ChemicalDb = ChemicalDb::from_reader(f).unwrap();

    //check that bounds checking works for invalid indices
    assert!(chem_db.get(100).is_none());
    assert!(chem_db.get(200).is_none());

    //check that get methods all work
    let entry = chem_db.get(20).unwrap();
    let exact_mass = entry.exact_mass();
    let inchikey14 = entry.inchikey14();
    let name = entry.name();

    assert_eq!(exact_mass, 271.21474379);
    assert_eq!(name, "test_mol42");
    assert_eq!(inchikey14, "NVXHLYRVFXNJGK");

    //check that get_mol correctly gives back the corresponding molecule structure for the
    //database row.
    assert!(entry.mol().is_ok());
}

#[test]
fn into_iterator() {
    let f = File::open("./test_files/test_chem_db.csv").unwrap();
    let chem_db = ChemicalDb::from_reader(f).unwrap();
    let ref_db = chem_db.clone();

    for (idx, row) in chem_db.into_iter().enumerate() {
        assert_eq!(ref_db.get(idx).unwrap(), &row);
    }
}

#[test]
fn dbcsv_opts() {
    let args = "fakebin -d src/lib.rs".split_whitespace();
    let opts = DbCsv::try_parse_from(args).unwrap();
    assert_eq!(&opts.path, &Path::new("src/lib.rs"));

    let args = "fakebin -d fake".split_whitespace();
    let opts = DbCsv::try_parse_from(args);
    assert!(opts.is_err());

    // path required, no default
    let args = "fakebin".split_whitespace();
    let opts = DbCsv::try_parse_from(args);
    assert!(opts.is_err());
}
