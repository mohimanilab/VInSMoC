use std::path::PathBuf;

use crate::{mgf::Mgf, mzml::Mzml, mzxml::MzXml, traits::MsFormat};

use super::*;
use approx::{assert_abs_diff_eq, AbsDiffEq};
use smallvec::SmallVec;
use spectrum::SpectrumMultistageMeta;

impl AbsDiffEq for Spectrum {
    type Epsilon = <f64 as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> <f64 as AbsDiffEq>::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: <f64 as AbsDiffEq>::Epsilon) -> bool {
        if !f64::abs_diff_eq(&self.pepmass, &other.pepmass, epsilon) {
            return false;
        }

        if self.charges.ne(&other.charges) {
            return false;
        }

        for (ind, peak) in self.peaks.iter().enumerate() {
            if !(f64::abs_diff_eq(&peak.mz, &other.peaks[ind].mz, epsilon)
                && f64::abs_diff_eq(&peak.intensity, &other.peaks[ind].intensity, epsilon))
            {
                return false;
            }
        }

        if self.retention_time.is_none() && other.retention_time.is_none() {
            return true;
        } else if self.retention_time.is_some() && other.retention_time.is_some() {
            return f64::abs_diff_eq(
                &self.retention_time.unwrap(),
                &other.retention_time.unwrap(),
                epsilon,
            );
        }

        false
    }
}

#[test]
fn test_log_rank() {
    let spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 0.0,
                intensity: 1.0,
            },
            Peak {
                mz: 0.0,
                intensity: 3.9,
            },
            Peak {
                mz: 0.0,
                intensity: 1.8,
            },
            Peak {
                mz: 0.0,
                intensity: 6.9,
            },
            Peak {
                mz: 0.0,
                intensity: 4.2,
            },
            Peak {
                mz: 0.0,
                intensity: 2.5,
            },
        ],
        pepmass: 0.0,
        charges: ChargeVec::new(),
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        scan: 0,
        denormalize_factor: None,
    };
    assert_eq!(vec![3, 2, 3, 1, 2, 3], spec.log_ranks(7));
    assert_eq!(vec![2, 2, 2, 1, 2, 2], spec.log_ranks(2));
}

#[test]
fn filter_peaks() {
    let mut spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 9.0,
                intensity: 0.1,
            },
            Peak {
                mz: 1.0,
                intensity: 2.0,
            },
            Peak {
                mz: 21.0,
                intensity: 2.0,
            },
            Peak {
                mz: 20.0,
                intensity: 302.0,
            },
            Peak {
                mz: 7.0,
                intensity: 2.2,
            },
        ],
        pepmass: 0.0,
        charges: ChargeVec::new(),
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        scan: 0,
        denormalize_factor: None,
    };

    let cmp_spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 7.0_f32 as f64,
                intensity: 2.2_f32 as f64,
            },
            Peak {
                mz: 9.0_f32 as f64,
                intensity: 0.1_f32 as f64,
            },
            Peak {
                mz: 20.0_f32 as f64,
                intensity: 302.0_f32 as f64,
            },
        ],
        pepmass: 0.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    spec.filter_peaks(7.0, 1);
    assert_abs_diff_eq!(spec, cmp_spec, epsilon = 1e-5);
}

#[test]
fn merge_peaks() {
    let mut spec1 = Spectrum {
        peaks: vec![
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 3.0_f32 as f64,
            },
            Peak {
                mz: 4.0_f32 as f64,
                intensity: 5.0_f32 as f64,
            },
            Peak {
                mz: 6.0_f32 as f64,
                intensity: 7.0_f32 as f64,
            },
        ],
        pepmass: 1.0_f32 as f64,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    let cmp_spec = Spectrum {
        peaks: vec![Peak {
            mz: 6.0_f32 as f64,
            intensity: 12.0_f32 as f64,
        }],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    spec1.merge_peaks(3.0, 1);
    assert_abs_diff_eq!(spec1, cmp_spec, epsilon = 1e-5);

    // test min peaks limit
    let mut spec2 = Spectrum {
        peaks: vec![
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 3.0_f32 as f64,
            },
            Peak {
                mz: 4.0_f32 as f64,
                intensity: 5.0_f32 as f64,
            },
            Peak {
                mz: 6.0_f32 as f64,
                intensity: 7.0_f32 as f64,
            },
        ],
        pepmass: 1.0_f32 as f64,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    spec2.merge_peaks(3.0, 2);

    let cmp_spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 3.0_f32 as f64,
            },
            Peak {
                mz: 4.0_f32 as f64,
                intensity: 5.0_f32 as f64,
            },
            Peak {
                mz: 6.0_f32 as f64,
                intensity: 7.0_f32 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    assert_abs_diff_eq!(spec2, cmp_spec, epsilon = 1e-5);

    let mut spec3 = Spectrum {
        peaks: vec![
            Peak {
                mz: 9.0_f32 as f64,
                intensity: 0.1_f32 as f64,
            },
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 2.0_f32 as f64,
            },
            Peak {
                mz: 21.0_f32 as f64,
                intensity: 2.0_f32 as f64,
            },
            Peak {
                mz: 20.0_f32 as f64,
                intensity: 302.0_f32 as f64,
            },
            Peak {
                mz: 7.0_f32 as f64,
                intensity: 2.2_f32 as f64,
            },
        ],
        pepmass: 0.0_f32 as f64,
        charges: ChargeVec::new(),
        multistage_meta: SpectrumMultistageMeta::default(),
        retention_time: None,
        scan: 0,
        denormalize_factor: None,
    };

    let cmp_spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 2.0_f32 as f64,
            },
            Peak {
                mz: 7.0_f32 as f64,
                intensity: 2.3_f32 as f64,
            },
            Peak {
                mz: 20.0_f32 as f64,
                intensity: 304.0_f32 as f64,
            },
        ],
        pepmass: 0.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    spec3.merge_peaks(2.1, 1);
    assert_abs_diff_eq!(spec3, cmp_spec, epsilon = 1e-5);
}

#[test]
fn filter_top_n_spectrum_peaks() {
    let mut spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 2.0_f32 as f64,
            },
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 3.0_f32 as f64,
            },
            Peak {
                mz: 4.0_f32 as f64,
                intensity: 5.0_f32 as f64,
            },
            Peak {
                mz: 6.0_f32 as f64,
                intensity: 7.0_f32 as f64,
            },
            Peak {
                mz: 7.0_f32 as f64,
                intensity: 2.3_f32 as f64,
            },
            Peak {
                mz: 20.0_f32 as f64,
                intensity: 304.0_f32 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let cmp_spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 4.0_f32 as f64,
                intensity: 5.0_f32 as f64,
            },
            Peak {
                mz: 6.0_f32 as f64,
                intensity: 7.0_f32 as f64,
            },
            Peak {
                mz: 20.0_f32 as f64,
                intensity: 304.0_f32 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    spec.filter_top_n_spectrum_peaks(3);
    assert_abs_diff_eq!(spec, cmp_spec, epsilon = 1e-5);
}

#[test]
fn normalize_peaks() {
    let mut spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 3.0_f32 as f64,
            },
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 4.0_f32 as f64,
            },
            Peak {
                mz: 3.0_f32 as f64,
                intensity: 9.0_f32 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    let res = Spectrum {
        peaks: vec![
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 0.433012_f32 as f64,
            },
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 0.500000_f32 as f64,
            },
            Peak {
                mz: 3.0_f32 as f64,
                intensity: 0.750000_f32 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let mut denormalize_factor = spec.normalize_peaks();
    assert_eq!(denormalize_factor, Some(4.0));
    assert_abs_diff_eq!(spec, res, epsilon = 1e-4);

    spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 0.0_f32 as f64,
            },
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 0.0_f32 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    let res = Spectrum {
        peaks: vec![
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 0.0 as f64,
            },
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 0.0 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    denormalize_factor = spec.normalize_peaks();
    assert_eq!(denormalize_factor, Some(0.0));
    assert_abs_diff_eq!(spec, res, epsilon = 1e-7);
}

#[test]
fn denormalize_peaks() {
    let initial_spec = Spectrum {
        peaks: vec![
            Peak {
                mz: 1.0_f32 as f64,
                intensity: 3.0_f32 as f64,
            },
            Peak {
                mz: 2.0_f32 as f64,
                intensity: 4.0_f32 as f64,
            },
            Peak {
                mz: 3.0_f32 as f64,
                intensity: 9.0_f32 as f64,
            },
        ],
        pepmass: 1.0,
        retention_time: None,
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let mut round_trip_spec = initial_spec.clone();

    let denormalize_factor = round_trip_spec.normalize_peaks();
    assert_eq!(denormalize_factor, Some(4.0));
    // make sure that denormalizing based on the returned denormalize_factor returns the
    // original spectrum
    round_trip_spec.denormalize_peaks(denormalize_factor.unwrap());
    assert_abs_diff_eq!(initial_spec, round_trip_spec, epsilon = 1e-4);
}

#[test]
fn round_trip_all() {
    let mgf = SpectrumCollection::from(
        Mgf::from_path("test_files/round_trip/tiny1.msconvert.mgf").unwrap(),
    );
    let mzml = SpectrumCollection::from(
        Mzml::from_path("test_files/round_trip/tiny1.msconvert.mzML").unwrap(),
    );
    let mzxml2 = SpectrumCollection::from(
        MzXml::from_path("test_files/round_trip/tiny1.mzXML2.0.mzXML").unwrap(),
    );
    let mzxml3 = SpectrumCollection::from(
        MzXml::from_path("test_files/round_trip/tiny1.msconvert.mzXML").unwrap(),
    );

    assert_eq!(mgf.spectra.len(), mzxml2.spectra.len());
    for (mgf_s, xml_s) in mgf.spectra.iter().zip(mzxml2.spectra.into_iter()) {
        assert_abs_diff_eq!(mgf_s, &xml_s, epsilon = 0.0001);
    }
    assert_eq!(mgf.spectra.len(), mzxml3.spectra.len());
    for (mgf_s, xml_s) in mgf.spectra.iter().zip(mzxml3.spectra.into_iter()) {
        assert_abs_diff_eq!(mgf_s, &xml_s, epsilon = 0.0001);
    }
    assert_eq!(mgf.spectra.len(), mzml.spectra.len());
    for (mgf_s, xml_s) in mgf.spectra.iter().zip(mzml.spectra.into_iter()) {
        assert_abs_diff_eq!(mgf_s, &xml_s, epsilon = 0.0001);
    }
}

#[test]
fn test_preprocessing_singleton0() {
    let mgf = SpectrumCollection::from(
        Mgf::from_path("test_files/singletons/singleton_lib0.mgf").unwrap(),
    );

    let mut target_spectra = Spectrum {
        peaks: vec![
            Peak {
                mz: 91.0541_f32 as f64,
                intensity: 17.0_f32 as f64,
            },
            Peak {
                mz: 91.2815_f32 as f64,
                intensity: 32.0_f32 as f64,
            },
            Peak {
                mz: 95.0505_f32 as f64,
                intensity: 55.0_f32 as f64,
            },
            Peak {
                mz: 119.049_f32 as f64,
                intensity: 237.0_f32 as f64,
            },
            Peak {
                mz: 123.044_f32 as f64,
                intensity: 530.0_f32 as f64,
            },
            Peak {
                mz: 136.076_f32 as f64,
                intensity: 1013.0_f32 as f64,
            },
            Peak {
                mz: 137.086_f32 as f64,
                intensity: 99.0_f32 as f64,
            },
            Peak {
                mz: 147.044_f32 as f64,
                intensity: 126.0_f32 as f64,
            },
            Peak {
                mz: 159.947_f32 as f64,
                intensity: 33.0_f32 as f64,
            },
            Peak {
                mz: 165.019_f32 as f64,
                intensity: 35.0_f32 as f64,
            },
            Peak {
                mz: 165.054_f32 as f64,
                intensity: 1060.0_f32 as f64,
            },
            Peak {
                mz: 182.081_f32 as f64,
                intensity: 633.0_f32 as f64,
            },
            Peak {
                mz: 183.078_f32 as f64,
                intensity: 52.0_f32 as f64,
            },
        ],
        pepmass: 182.0807695,
        retention_time: Some(30.774),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    let normalized_target = Spectrum {
        peaks: vec![
            Peak {
                mz: 91.0541_f32 as f64,
                intensity: 0.0658371_f32 as f64,
            },
            Peak {
                mz: 91.2815_f32 as f64,
                intensity: 0.0903278_f32 as f64,
            },
            Peak {
                mz: 95.0505_f32 as f64,
                intensity: 0.118421_f32 as f64,
            },
            Peak {
                mz: 119.049_f32 as f64,
                intensity: 0.245822_f32 as f64,
            },
            Peak {
                mz: 123.044_f32 as f64,
                intensity: 0.367607_f32 as f64,
            },
            Peak {
                mz: 136.076_f32 as f64,
                intensity: 0.508219_f32 as f64,
            },
            Peak {
                mz: 137.086_f32 as f64,
                intensity: 0.158878_f32 as f64,
            },
            Peak {
                mz: 147.044_f32 as f64,
                intensity: 0.179239_f32 as f64,
            },
            Peak {
                mz: 159.947_f32 as f64,
                intensity: 0.0917283_f32 as f64,
            },
            Peak {
                mz: 165.019_f32 as f64,
                intensity: 0.094467_f32 as f64,
            },
            Peak {
                mz: 165.054_f32 as f64,
                intensity: 0.519875_f32 as f64,
            },
            Peak {
                mz: 182.081_f32 as f64,
                intensity: 0.401743_f32 as f64,
            },
            Peak {
                mz: 183.078_f32 as f64,
                intensity: 0.115146_f32 as f64,
            },
        ],
        pepmass: 182.0807695,
        retention_time: Some(30.774),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let mut mgf_s = mgf.spectra[0].clone();
    mgf_s.filter_peaks_masst(50.0, 5);
    assert_abs_diff_eq!(mgf_s, target_spectra, epsilon = 1e-2);
    target_spectra.normalize_peaks();
    assert_abs_diff_eq!(target_spectra, normalized_target, epsilon = 1e-4);
}

#[test]
fn test_preprocessing_singleton1() {
    let mgf = SpectrumCollection::from(
        Mgf::from_path("test_files/singletons/singleton_lib1.mgf").unwrap(),
    );
    let mut target_spectra = Spectrum {
        peaks: vec![
            Peak {
                mz: 75.5112_f32 as f64,
                intensity: 38.0_f32 as f64,
            },
            Peak {
                mz: 79.0537_f32 as f64,
                intensity: 25.0_f32 as f64,
            },
            Peak {
                mz: 84.9611_f32 as f64,
                intensity: 112.0_f32 as f64,
            },
            Peak {
                mz: 93.069_f32 as f64,
                intensity: 24.0_f32 as f64,
            },
            Peak {
                mz: 103.054_f32 as f64,
                intensity: 95.0_f32 as f64,
            },
            Peak {
                mz: 120.081_f32 as f64,
                intensity: 3263.0_f32 as f64,
            },
            Peak {
                mz: 121.085_f32 as f64,
                intensity: 228.0_f32 as f64,
            },
            Peak {
                mz: 131.05_f32 as f64,
                intensity: 102.0_f32 as f64,
            },
            Peak {
                mz: 149.059_f32 as f64,
                intensity: 105.0_f32 as f64,
            },
            Peak {
                mz: 161.415_f32 as f64,
                intensity: 35.0_f32 as f64,
            },
            Peak {
                mz: 166.087_f32 as f64,
                intensity: 1216.0_f32 as f64,
            },
            Peak {
                mz: 167.078_f32 as f64,
                intensity: 25.0_f32 as f64,
            },
            Peak {
                mz: 168.435_f32 as f64,
                intensity: 18.0_f32 as f64,
            },
        ],
        pepmass: 166.08576565,
        retention_time: Some(52.251),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    let normalized_target = Spectrum {
        peaks: vec![
            Peak {
                mz: 75.5112_f32 as f64,
                intensity: 0.0847868_f32 as f64,
            },
            Peak {
                mz: 79.0537_f32 as f64,
                intensity: 0.0687712_f32 as f64,
            },
            Peak {
                mz: 84.9611_f32 as f64,
                intensity: 0.145561_f32 as f64,
            },
            Peak {
                mz: 93.069_f32 as f64,
                intensity: 0.0673817_f32 as f64,
            },
            Peak {
                mz: 103.054_f32 as f64,
                intensity: 0.13406_f32 as f64,
            },
            Peak {
                mz: 120.081_f32 as f64,
                intensity: 0.785679_f32 as f64,
            },
            Peak {
                mz: 121.085_f32 as f64,
                intensity: 0.207684_f32 as f64,
            },
            Peak {
                mz: 131.05_f32 as f64,
                intensity: 0.138911_f32 as f64,
            },
            Peak {
                mz: 149.059_f32 as f64,
                intensity: 0.140939_f32 as f64,
            },
            Peak {
                mz: 161.415_f32 as f64,
                intensity: 0.0813711_f32 as f64,
            },
            Peak {
                mz: 166.087_f32 as f64,
                intensity: 0.479627_f32 as f64,
            },
            Peak {
                mz: 167.078_f32 as f64,
                intensity: 0.0687712_f32 as f64,
            },
            Peak {
                mz: 168.435_f32 as f64,
                intensity: 0.0583543_f32 as f64,
            },
        ],
        pepmass: 166.08576565,
        retention_time: Some(52.251),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let mut mgf_s = mgf.spectra[0].clone();
    mgf_s.filter_peaks_masst(50.0, 5);
    assert_abs_diff_eq!(mgf_s, target_spectra, epsilon = 1e-2);
    target_spectra.normalize_peaks();
    assert_abs_diff_eq!(target_spectra, normalized_target, epsilon = 1e-4);
}

#[test]
fn test_preprocessing_singleton2() {
    let mgf = SpectrumCollection::from(
        Mgf::from_path("test_files/singletons/singleton_lib2.mgf").unwrap(),
    );
    let mut target_spectra = Spectrum {
        peaks: vec![
            Peak {
                mz: 84.9594_f32 as f64,
                intensity: 82.0_f32 as f64,
            },
            Peak {
                mz: 103.057_f32 as f64,
                intensity: 147.0_f32 as f64,
            },
            Peak {
                mz: 120.081_f32 as f64,
                intensity: 6200.0_f32 as f64,
            },
            Peak {
                mz: 121.085_f32 as f64,
                intensity: 288.0_f32 as f64,
            },
            Peak {
                mz: 131.05_f32 as f64,
                intensity: 143.0_f32 as f64,
            },
            Peak {
                mz: 149.06_f32 as f64,
                intensity: 231.0_f32 as f64,
            },
            Peak {
                mz: 162.015_f32 as f64,
                intensity: 28.0_f32 as f64,
            },
            Peak {
                mz: 165.11_f32 as f64,
                intensity: 30.0_f32 as f64,
            },
            Peak {
                mz: 166.086_f32 as f64,
                intensity: 2437.0_f32 as f64,
            },
            Peak {
                mz: 166.981_f32 as f64,
                intensity: 30.0_f32 as f64,
            },
            Peak {
                mz: 167.089_f32 as f64,
                intensity: 124.0_f32 as f64,
            },
        ],
        pepmass: 166.08574209,
        retention_time: Some(54.807),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let normalized_target = Spectrum {
        peaks: vec![
            Peak {
                mz: 84.9594_f32 as f64,
                intensity: 0.0917545_f32 as f64,
            },
            Peak {
                mz: 103.057_f32 as f64,
                intensity: 0.122851_f32 as f64,
            },
            Peak {
                mz: 120.081_f32 as f64,
                intensity: 0.797841_f32 as f64,
            },
            Peak {
                mz: 121.085_f32 as f64,
                intensity: 0.171956_f32 as f64,
            },
            Peak {
                mz: 131.05_f32 as f64,
                intensity: 0.121168_f32 as f64,
            },
            Peak {
                mz: 149.06_f32 as f64,
                intensity: 0.154002_f32 as f64,
            },
            Peak {
                mz: 162.015_f32 as f64,
                intensity: 0.0536166_f32 as f64,
            },
            Peak {
                mz: 165.11_f32 as f64,
                intensity: 0.0554985_f32 as f64,
            },
            Peak {
                mz: 166.086_f32 as f64,
                intensity: 0.500205_f32 as f64,
            },
            Peak {
                mz: 166.981_f32 as f64,
                intensity: 0.0554985_f32 as f64,
            },
            Peak {
                mz: 167.089_f32 as f64,
                intensity: 0.112832_f32 as f64,
            },
        ],
        pepmass: 166.08574209,
        retention_time: Some(54.807),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let mut mgf_s = mgf.spectra[0].clone();
    mgf_s.filter_peaks_masst(50.0, 5);
    assert_abs_diff_eq!(mgf_s, target_spectra, epsilon = 1e-2);
    target_spectra.normalize_peaks();
    assert_abs_diff_eq!(target_spectra, normalized_target, epsilon = 1e-4);
}

#[test]
fn test_preprocessing_singleton3() {
    let mgf = SpectrumCollection::from(
        Mgf::from_path("test_files/singletons/singleton_lib3.mgf").unwrap(),
    );
    let mut target_spectra = Spectrum {
        peaks: vec![
            Peak {
                mz: 118.065_f32 as f64,
                intensity: 160.0_f32 as f64,
            },
            Peak {
                mz: 136.071_f32 as f64,
                intensity: 31.0_f32 as f64,
            },
            Peak {
                mz: 143.071_f32 as f64,
                intensity: 30.0_f32 as f64,
            },
            Peak {
                mz: 144.081_f32 as f64,
                intensity: 50.0_f32 as f64,
            },
            Peak {
                mz: 146.06_f32 as f64,
                intensity: 176.0_f32 as f64,
            },
            Peak {
                mz: 159.09_f32 as f64,
                intensity: 61.0_f32 as f64,
            },
            Peak {
                mz: 164.087_f32 as f64,
                intensity: 33.0_f32 as f64,
            },
            Peak {
                mz: 170.061_f32 as f64,
                intensity: 46.0_f32 as f64,
            },
            Peak {
                mz: 188.07_f32 as f64,
                intensity: 745.0_f32 as f64,
            },
            Peak {
                mz: 189.074_f32 as f64,
                intensity: 107.0_f32 as f64,
            },
            Peak {
                mz: 200.134_f32 as f64,
                intensity: 29.0_f32 as f64,
            },
            Peak {
                mz: 204.104_f32 as f64,
                intensity: 32.0_f32 as f64,
            },
            Peak {
                mz: 205.097_f32 as f64,
                intensity: 119.0_f32 as f64,
            },
            Peak {
                mz: 206.099_f32 as f64,
                intensity: 22.0_f32 as f64,
            },
            Peak {
                mz: 222.933_f32 as f64,
                intensity: 41.0_f32 as f64,
            },
        ],
        pepmass: 205.09703008,
        retention_time: Some(120.105),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let normalized_target = Spectrum {
        peaks: vec![
            Peak {
                mz: 118.065_f32 as f64,
                intensity: 0.308423_f32 as f64,
            },
            Peak {
                mz: 136.071_f32 as f64,
                intensity: 0.135759_f32 as f64,
            },
            Peak {
                mz: 143.071_f32 as f64,
                intensity: 0.133551_f32 as f64,
            },
            Peak {
                mz: 144.081_f32 as f64,
                intensity: 0.172414_f32 as f64,
            },
            Peak {
                mz: 146.06_f32 as f64,
                intensity: 0.323477_f32 as f64,
            },
            Peak {
                mz: 159.09_f32 as f64,
                intensity: 0.190437_f32 as f64,
            },
            Peak {
                mz: 164.087_f32 as f64,
                intensity: 0.14007_f32 as f64,
            },
            Peak {
                mz: 170.061_f32 as f64,
                intensity: 0.165374_f32 as f64,
            },
            Peak {
                mz: 188.07_f32 as f64,
                intensity: 0.665526_f32 as f64,
            },
            Peak {
                mz: 189.074_f32 as f64,
                intensity: 0.25222_f32 as f64,
            },
            Peak {
                mz: 200.134_f32 as f64,
                intensity: 0.131306_f32 as f64,
            },
            Peak {
                mz: 204.104_f32 as f64,
                intensity: 0.137931_f32 as f64,
            },
            Peak {
                mz: 205.097_f32 as f64,
                intensity: 0.265987_f32 as f64,
            },
            Peak {
                mz: 206.099_f32 as f64,
                intensity: 0.114366_f32 as f64,
            },
            Peak {
                mz: 222.933_f32 as f64,
                intensity: 0.156127_f32 as f64,
            },
        ],
        pepmass: 205.09703008,
        retention_time: Some(120.105),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let mut mgf_s = mgf.spectra[0].clone();
    mgf_s.filter_peaks_masst(50.0, 5);
    assert_abs_diff_eq!(mgf_s, target_spectra, epsilon = 1e-2);
    target_spectra.normalize_peaks();
    assert_abs_diff_eq!(target_spectra, normalized_target, epsilon = 1e-4);
}

#[test]
fn test_preprocessing_singleton4() {
    let mgf = SpectrumCollection::from(
        Mgf::from_path("test_files/singletons/singleton_lib4.mgf").unwrap(),
    );
    let mut target_spectra = Spectrum {
        peaks: vec![
            Peak {
                mz: 118.064_f32 as f64,
                intensity: 145.0_f32 as f64,
            },
            Peak {
                mz: 130.062_f32 as f64,
                intensity: 72.0_f32 as f64,
            },
            Peak {
                mz: 132.083_f32 as f64,
                intensity: 95.0_f32 as f64,
            },
            Peak {
                mz: 144.08_f32 as f64,
                intensity: 122.0_f32 as f64,
            },
            Peak {
                mz: 146.059_f32 as f64,
                intensity: 646.0_f32 as f64,
            },
            Peak {
                mz: 154.302_f32 as f64,
                intensity: 49.0_f32 as f64,
            },
            Peak {
                mz: 159.091_f32 as f64,
                intensity: 117.0_f32 as f64,
            },
            Peak {
                mz: 170.059_f32 as f64,
                intensity: 58.0_f32 as f64,
            },
            Peak {
                mz: 188.069_f32 as f64,
                intensity: 2116.0_f32 as f64,
            },
            Peak {
                mz: 189.072_f32 as f64,
                intensity: 118.0_f32 as f64,
            },
            Peak {
                mz: 204.007_f32 as f64,
                intensity: 22.0_f32 as f64,
            },
            Peak {
                mz: 204.923_f32 as f64,
                intensity: 19.0_f32 as f64,
            },
            Peak {
                mz: 205.095_f32 as f64,
                intensity: 201.0_f32 as f64,
            },
            Peak {
                mz: 206.104_f32 as f64,
                intensity: 54.0_f32 as f64,
            },
            Peak {
                mz: 206.132_f32 as f64,
                intensity: 18.0_f32 as f64,
            },
        ],
        pepmass: 205.09632527,
        retention_time: Some(120.996),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let normalized_target = Spectrum {
        peaks: vec![
            Peak {
                mz: 118.064_f32 as f64,
                intensity: 0.194017_f32 as f64,
            },
            Peak {
                mz: 130.062_f32 as f64,
                intensity: 0.136717_f32 as f64,
            },
            Peak {
                mz: 132.083_f32 as f64,
                intensity: 0.157043_f32 as f64,
            },
            Peak {
                mz: 144.08_f32 as f64,
                intensity: 0.177966_f32 as f64,
            },
            Peak {
                mz: 146.059_f32 as f64,
                intensity: 0.409518_f32 as f64,
            },
            Peak {
                mz: 154.302_f32 as f64,
                intensity: 0.112786_f32 as f64,
            },
            Peak {
                mz: 159.091_f32 as f64,
                intensity: 0.174281_f32 as f64,
            },
            Peak {
                mz: 170.059_f32 as f64,
                intensity: 0.122707_f32 as f64,
            },
            Peak {
                mz: 188.069_f32 as f64,
                intensity: 0.741165_f32 as f64,
            },
            Peak {
                mz: 189.072_f32 as f64,
                intensity: 0.175024_f32 as f64,
            },
            Peak {
                mz: 204.007_f32 as f64,
                intensity: 0.0755733_f32 as f64,
            },
            Peak {
                mz: 204.923_f32 as f64,
                intensity: 0.0702318_f32 as f64,
            },
            Peak {
                mz: 205.095_f32 as f64,
                intensity: 0.228431_f32 as f64,
            },
            Peak {
                mz: 206.104_f32 as f64,
                intensity: 0.118401_f32 as f64,
            },
            Peak {
                mz: 206.132_f32 as f64,
                intensity: 0.0683586_f32 as f64,
            },
        ],
        pepmass: 205.09632527,
        retention_time: Some(120.996),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };
    let mut mgf_s = mgf.spectra[0].clone();
    mgf_s.filter_peaks_masst(50.0, 5);
    assert_abs_diff_eq!(mgf_s, target_spectra, epsilon = 1e-2);
    target_spectra.normalize_peaks();
    assert_abs_diff_eq!(target_spectra, normalized_target, epsilon = 1e-4);
}

#[test]
fn write_mgf_newlines() {
    let spectrum = Spectrum {
        peaks: vec![Peak {
            mz: 1.0,
            intensity: 2.0,
        }],
        pepmass: 3.0,
        retention_time: Some(4.0),
        multistage_meta: SpectrumMultistageMeta::default(),
        charges: ChargeVec::new(),
        scan: 0,
        denormalize_factor: None,
    };

    let mut buffer = Vec::new();
    spectrum.write_mgf(&mut buffer).unwrap();
    let mgf = String::from_utf8(buffer).unwrap();
    let answer = r#"BEGIN IONS
PEPMASS=3
SCANS=0
RTINSECONDS=4
1 2
END IONS
"#;
    assert_eq!(mgf, answer);
}

fn dummy_msn(ms_level: usize, scan: usize, precursor_scan: Option<usize>) -> Spectrum {
    Spectrum::new(
        Vec::new(),
        0.0,
        SmallVec::new(),
        None,
        SpectrumMultistageMeta {
            ms_level: Some(ms_level),
            precursor_intensity: None,
            precursor_scan,
        },
        scan,
    )
}

#[test]
fn multistage_tree_precursor_idx_search() {
    // first spectrum collection
    let precursor_1 = dummy_msn(1, 1, None);
    let product_1 = dummy_msn(2, 3, Some(1));
    let product_2 = dummy_msn(2, 5, Some(1));

    // second spectrum collection
    let precursor_3 = dummy_msn(1, 3, None);
    let product_3 = dummy_msn(2, 4, Some(3));
    let product_4 = dummy_msn(2, 5, Some(3));

    let spectral_data = vec![
        precursor_1,
        product_1,
        product_2,
        precursor_3,
        product_3,
        product_4,
    ];
    let metadata = vec![(PathBuf::default(), 3), (PathBuf::default(), 6)];

    let tree = MultistageTree::from_ingested_spectra(&spectral_data, &metadata);
    assert_eq!(tree.product_idxs(0).collect::<Vec<_>>(), vec![1, 2]);
    assert_eq!(tree.product_idxs(3).collect::<Vec<_>>(), vec![4, 5]);
    assert_eq!(tree.all_edges().count(), 4);
}
