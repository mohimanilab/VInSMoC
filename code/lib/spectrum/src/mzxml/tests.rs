use approx::assert_abs_diff_eq;

use super::MzXml;
use crate::{mgf::Mgf, traits::MsFormat};

#[test]
fn mzxml_mgf_round_trip() {
    let mzxml = MzXml::from_path("test_files/mzxml/P434_E7_P3-6I_crude.mzXML");
    assert!(mzxml.is_ok(), "{:?}", mzxml);
    let mzxml = mzxml.unwrap();
    let mzxml = crate::SpectrumCollection::from(mzxml);
    let mgf = Mgf::from_path("test_files/mzxml/P434_E7_P3-6I_crude.mgf").unwrap();
    let mgf = crate::SpectrumCollection::from(mgf);

    assert_eq!(mgf.spectra.len(), mzxml.spectra.len());
    for (mgf_spectrum, mzxml_spectrum) in mgf.spectra.into_iter().zip(mzxml.spectra.into_iter()) {
        assert_abs_diff_eq!(mgf_spectrum, mzxml_spectrum, epsilon = 0.00001);
    }
}

#[test]
fn nested_scans() {
    let mzxml = MzXml::from_path("test_files/mzxml/nested.mzXML");
    assert!(mzxml.is_ok(), "{:?}", mzxml);
    let mzxml = mzxml.unwrap();
    let mzxml = crate::SpectrumCollection::from(mzxml);
    // doesn't exist
    assert!(mzxml.get_scan(0).is_none());
    // is MS1
    assert!(mzxml.get_scan(1).is_none());

    // all good
    for scan in 2..=6 {
        assert!(mzxml.get_scan(scan).is_some());
        // check that indexing is correct
        assert_eq!(mzxml.spectra[scan - 2].scan, scan);
    }

    // is MS1
    assert!(mzxml.get_scan(7).is_none());

    // all good
    for scan in 8..=12 {
        assert!(mzxml.get_scan(scan).is_some());
        // check that indexing is correct
        assert_eq!(mzxml.spectra[scan - 3].scan, scan);
    }

    assert_eq!(mzxml.get_scan(6).unwrap().peaks.len(), 286);
    assert_eq!(mzxml.get_scan(9).unwrap().peaks.len(), 274);
}

#[test]
fn missing_content_type() {
    let mzxml = MzXml::from_path("test_files/mzxml/CNP082_MS_M.mzXML");
    assert!(mzxml.is_ok(), "{:?}", mzxml);
}
