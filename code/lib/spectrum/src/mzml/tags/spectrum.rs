use std::ops::{Deref, DerefMut};

use crate::{mzml::Error, spectrum::SpectrumMultistageMeta, ChargeVec, Peak};
use serde::Deserialize;

use super::{
    bindata::BinDataList,
    cvparam::{CvParam, CvParams},
    precursor::PrecursorList,
    scan::ScanList,
};

/// A single spectrum from an mzML file.
///
/// This can only be created as part of an [`Mzml`](super::Mzml).
#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Spectrum {
    index: usize,
    id: Option<String>,
    #[serde(rename = "cvParam", default)]
    cv_params: Vec<CvParam>,
    binary_data_array_list: Option<BinDataList>,
    precursor_list: Option<PrecursorList>,
    scan_list: Option<ScanList>,

    // These will never be parsed, but they are populated
    // during validation for easy future access
    #[serde(skip)]
    peaks: Vec<Peak>,
    #[serde(skip)]
    pepmass: f64,
    #[serde(skip)]
    precursor_intensity: Option<f64>,
    #[serde(skip)]
    retention_time: Option<f64>,
    #[serde(skip)]
    charges: ChargeVec,

    // These are used to determine parent-child relationships in tree construction.
    #[serde(skip)]
    precursor_id: Option<String>,
    #[serde(skip)]
    precursor_idx: Option<usize>,
    #[serde(skip)]
    ms_level: Option<usize>,
}

impl CvParams for Spectrum {
    fn params(&self) -> &Vec<CvParam> {
        &self.cv_params
    }
}

const MS_LEVEL_SPECIFIER: &str = "ms level";

impl Spectrum {
    /// Parse the contents of this spectrum.
    ///
    /// This converts the raw spectral data into peak vectors, the raw precursor data into a
    /// precursor m/z and charge, and the raw retention time into a floating point.
    ///
    /// This will happily leave fields blank, so it should only be called as part of the mzML
    /// parsing process, where pre-conditions are validated.
    ///
    /// # Errors
    /// See [`Error`] for a superset possible error states.
    pub(crate) fn parse(&mut self, multistage: bool) -> Result<(), Error> {
        self.ms_level = self
            .param_with_name(MS_LEVEL_SPECIFIER)
            .map(|x| x.value::<usize>())
            .transpose()?
            .flatten();
        if let Some(bin_data) = self.binary_data_array_list.as_ref() {
            self.peaks = bin_data.parse()?;
        }

        if let Some(precursor_list) = self.precursor_list.as_ref() {
            let (precursor, charges, precursor_id, precursor_intensity) =
                precursor_list.parse(multistage)?;
            self.pepmass = precursor;
            self.charges = charges;
            self.precursor_id = precursor_id;
            self.precursor_intensity = precursor_intensity;
        }

        if let Some(scan_list) = self.scan_list.as_ref() {
            self.retention_time = scan_list.parse()?;
        }

        Ok(())
    }

    /// Checks if this spectrum is a tandem mass spectrum.
    ///
    /// We use the presence of a child precursorList element has an indicated that this spectrum is
    /// a tandem mass spectrum.
    pub fn is_tandem(&self) -> bool {
        self.precursor_list.is_some()
    }

    /// Get the spectrum's index.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Get the unique string identifier for this spectrum.
    pub fn id(&self) -> Option<String> {
        self.id.to_owned()
    }

    /// Get the unique string identifier for this spectrum's precursor spectrum.
    pub fn precursor_id(&self) -> Option<String> {
        self.precursor_id.to_owned()
    }

    /// Set the spectrum index corresponding to this spectrum's precursor spectrum.
    pub fn set_precursor_index(&mut self, precursor_idx: Option<usize>) {
        self.precursor_idx = precursor_idx;
    }
}

impl Deref for Spectrum {
    type Target = Vec<Peak>;

    fn deref(&self) -> &Self::Target {
        &self.peaks
    }
}

impl DerefMut for Spectrum {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.peaks
    }
}

impl From<Spectrum> for Vec<Peak> {
    fn from(spectrum: Spectrum) -> Self {
        spectrum.peaks
    }
}

impl crate::Ms for Spectrum {
    fn pepmass(&self) -> f64 {
        self.pepmass
    }

    fn charges(&self) -> &crate::ChargeVec {
        &self.charges
    }

    fn retention_time(&self) -> Option<f64> {
        self.retention_time
    }

    fn ms_level(&self) -> Option<usize> {
        self.ms_level
    }

    fn precursor_scan(&self) -> Option<usize> {
        self.precursor_idx
    }

    fn precursor_intensity(&self) -> Option<f64> {
        self.precursor_intensity
    }

    fn get_scan(&self) -> Option<usize> {
        Some(self.index)
    }
}

impl From<Spectrum> for crate::Spectrum {
    fn from(spectrum: Spectrum) -> Self {
        Self {
            peaks: spectrum.peaks,
            pepmass: spectrum.pepmass,
            multistage_meta: SpectrumMultistageMeta {
                ms_level: spectrum.ms_level,
                precursor_intensity: spectrum.precursor_intensity,
                precursor_scan: spectrum.precursor_idx,
            },
            charges: spectrum.charges,
            retention_time: spectrum.retention_time,
            scan: spectrum.index,
            denormalize_factor: None,
        }
    }
}
