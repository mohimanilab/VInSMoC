use std::str::FromStr;

use rustc_hash::FxHashMap;
use serde::Deserialize;

use crate::{mzml::Error, Ms, MsFormat};

use super::{run::Run, Spectrum};

/// A fully parsed mzML file.
///
/// The parsed version will only contain spectra with msLevel 2, meaning we discard all spectra
/// that aren't tandem mass spectra.
///
/// This should be constructed using its [`FromStr`] or [`MsFormat`] implementations.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Mzml {
    run: Run,
}

impl Mzml {
    /// Retrieve spectrum at given index.
    ///
    /// Returns [`None`] if the index is out of bounds.
    pub fn get(&self, spectrum_index: usize) -> Option<&super::Spectrum> {
        self.run.spectrum_list.spectra.get(spectrum_index)
    }

    /// Retrieve spectrum by scan number.
    ///
    /// Returns [`None`] if the scan number is not found.
    pub fn get_by_scan(&self, scan: usize) -> Option<&super::Spectrum> {
        self.run
            .spectrum_list
            .spectra
            .iter()
            .find(|spectrum| spectrum.get_scan() == Some(scan))
    }

    /// Get the number of spectra in this mzML file.
    pub fn len(&self) -> usize {
        self.run.spectrum_list.spectra.len()
    }

    /// Check if mzML contained any spectra
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Checks that the parsed mzml object is valid.
    ///
    /// Current checks: counts are correct; mz_array and intensity array exist with valid cvParams;
    ///                 indices are correct; needed values exist and can be parsed successfully.
    fn validate_and_parse(mut self, multistage: bool) -> Result<Self, Error> {
        let spectra = self.run.spectrum_list.spectra;

        self.run.spectrum_list.spectra = spectra
            .into_iter()
            .enumerate()
            .filter_map(|(ind, spec)| {
                if !multistage && !spec.is_tandem() {
                    return None;
                }
                Some((ind, spec))
            })
            .map(|(ind, mut spec)| -> Result<Spectrum, Error> {
                // Check the spectrum position in the spectrum list
                if spec.index() != ind {
                    return Err(Error::MismatchedIndex(spec.index(), ind));
                }

                spec.parse(multistage)?;

                Ok(spec)
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(self)
    }
}

impl FromStr for Mzml {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // get the actual string contents we want
        // we want to start with <mzML> and ...
        let start_idx = s.find("<mzML").ok_or(Error::NoMzmlTag)?;
        let (_, s) = s.split_at(start_idx);

        // ... end with </mzML>
        let end_tag = "</mzML>";
        let end_idx = s.find(end_tag).ok_or(Error::NoMzmlTag)? + end_tag.len();
        let end_idx = (end_idx + 1).min(s.len() - 1);
        let (s, _) = s.split_at(end_idx);

        let mzml = quick_xml::de::from_str::<Self>(s)?;
        // pass false because we don't want MS1
        mzml.validate_and_parse(false)
    }
}

impl MsFormat for Mzml {}

impl From<Mzml> for crate::SpectrumCollection {
    fn from(mzml: Mzml) -> Self {
        let mut spectra = mzml
            .run
            .spectrum_list
            .spectra
            .into_iter()
            .map(crate::Spectrum::from)
            .collect::<Vec<_>>();
        spectra.sort_by_key(|spectrum| spectrum.scan);

        Self { spectra }
    }
}

/// A fully parsed mzML file appropriate for multistage processing.
///
/// The parsed version will contain all spectra regardless of msLevel.
///
/// This should be constructed using its [`FromStr`] or [`MsFormat`] implementations.
#[derive(Debug)]
pub struct MzmlMultistage {
    spectra: Vec<Spectrum>,
}

impl MzmlMultistage {
    /// Retrieve spectrum at given index.
    ///
    /// Returns [`None`] if the index is out of bounds.
    pub fn get(&self, spectrum_index: usize) -> Option<&super::Spectrum> {
        self.spectra.get(spectrum_index)
    }

    /// Get the number of spectra in this mzML file.
    pub fn len(&self) -> usize {
        self.spectra.len()
    }

    /// Check if mzML contained any spectra
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl FromStr for MzmlMultistage {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // get the actual string contents we want
        // we want to start with <mzML> and ...
        let start_idx = s.find("<mzML").ok_or(Error::NoMzmlTag)?;
        let (_, s) = s.split_at(start_idx);

        // ... end with </mzML>
        let end_tag = "</mzML>";
        let end_idx = s.find(end_tag).ok_or(Error::NoMzmlTag)? + end_tag.len();
        let end_idx = (end_idx + 1).min(s.len() - 1);
        let (s, _) = s.split_at(end_idx);

        let mzml = quick_xml::de::from_str::<Mzml>(s)?;
        // pass true because we want all MS levels present
        mzml.validate_and_parse(true).map(|mzml| mzml.into())
    }
}

impl From<Mzml> for MzmlMultistage {
    fn from(mzml: Mzml) -> Self {
        let mut spectra = mzml.run.spectrum_list.spectra;
        // assumes that precursor scans always occur prior to product scans in a given spectrum
        // file
        let mut precursor_id_to_idx = <FxHashMap<String, usize>>::default();
        for spectrum in spectra.iter_mut() {
            if let Some(id) = spectrum.id() {
                precursor_id_to_idx.insert(id, spectrum.index());
            }
            if let Some(precursor_id) = spectrum.precursor_id() {
                spectrum.set_precursor_index(precursor_id_to_idx.get(&precursor_id).cloned())
            }
        }
        Self { spectra }
    }
}

impl From<MzmlMultistage> for crate::SpectrumCollection {
    fn from(mzml_multistage: MzmlMultistage) -> Self {
        let mut spectra = mzml_multistage
            .spectra
            .into_iter()
            .map(crate::Spectrum::from)
            .collect::<Vec<_>>();
        spectra.sort_by_key(|spectrum| spectrum.scan);

        Self { spectra }
    }
}

impl MsFormat for MzmlMultistage {}
