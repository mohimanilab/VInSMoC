use std::{collections::BTreeMap, fs::File, io, path::Path};

use molecule::Adduct;
use ordered_float::OrderedFloat;

use crate::{
    mgf::Mgf,
    mzml::{Mzml, MzmlMultistage},
    mzxml::{MzXml, MzXmlMultistage},
    traits::MsFormat,
    Error, Spectrum, SpectrumPreprocessing,
};

/// A general purpose collection of tandem mass spectra, not tied to any specific format.
#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct SpectrumCollection {
    /// The constituent spectra for this collection
    pub spectra: Vec<Spectrum>,
}

impl SpectrumCollection {
    /// Parse a spectrum collection from an MGF file.
    pub fn from_mgf_file(f: File) -> Result<Self, Error> {
        let f = Mgf::from_file(f)?;
        Ok(f.into())
    }

    /// Parse a spectrum collection from an mzMXL file.
    pub fn from_mzxml_file(f: File) -> Result<Self, Error> {
        let f = MzXml::from_file(f)?;
        Ok(f.into())
    }

    /// Parse a multistage spectrum collection from an mzMXL file.
    pub fn from_mzxml_file_multistage(f: File) -> Result<Self, Error> {
        let f = MzXmlMultistage::from_file(f)?;
        Ok(f.into())
    }

    /// Parse a spectrum collection from an mzML file.
    pub fn from_mzml_file(f: File) -> Result<Self, Error> {
        let f = Mzml::from_file(f)?;
        Ok(f.into())
    }

    /// Parse a multistage spectrum collection from an mzML file.
    pub fn from_mzml_file_multistage(f: File) -> Result<Self, Error> {
        let f = MzmlMultistage::from_file(f)?;
        Ok(f.into())
    }

    /// Pre-process this collection.
    ///
    /// For every scan in this collection it performs the pre-processing steps described in
    /// [`SpectrumPreprocessing`]. Additionally, the spectra in this collection are sorted by
    /// ascending scan.
    pub fn preprocess(&mut self, params: &SpectrumPreprocessing) {
        let mut tmp = Vec::new();
        std::mem::swap(&mut tmp, &mut self.spectra);
        self.spectra = tmp
            .into_iter()
            .filter_map(|spectrum| params.preprocess(spectrum))
            .collect();
        self.spectra.sort_unstable_by_key(|spectrum| spectrum.scan);
    }

    /// Parse a spectrum collection from a path.
    ///
    /// The spectrum format is selected based on the extension of the file. See
    /// [`Error::UnknownFormat`] for information about the matching process.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, Error> {
        let file = File::open(path.as_ref())?;

        let ext = path
            .as_ref()
            .extension()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "spectrum file had no extension",
                )
            })?
            .to_ascii_lowercase();

        let ext = ext.to_str().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "file extension wasn't valid unicode",
            )
        })?;

        let spectrum_collection = match ext {
            "mgf" => SpectrumCollection::from_mgf_file(file),
            "mzxml" => SpectrumCollection::from_mzxml_file(file),
            "mzml" => SpectrumCollection::from_mzml_file(file),
            _ => return Err(Error::UnknownFormat),
        }?;

        Ok(spectrum_collection)
    }

    /// Parse a multistage spectrum collection from a path.
    ///
    /// The spectrum format is selected based on the extension of the file. See
    /// [`Error::UnknownFormat`] for information about the matching process.
    pub fn from_path_multistage(path: impl AsRef<Path>) -> Result<Self, Error> {
        let file = File::open(path.as_ref())?;

        let ext = path
            .as_ref()
            .extension()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "spectrum file had no extension",
                )
            })?
            .to_ascii_lowercase();

        let ext = ext.to_str().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "file extension wasn't valid unicode",
            )
        })?;

        let spectrum_collection = match ext {
            // mgf doesn't have any MS level information, so we can just keep its original
            // parsing
            "mgf" => SpectrumCollection::from_mgf_file(file),
            // mzml and mzxml however need to be multistage-parsed, not discarding any MS levels
            "mzxml" => SpectrumCollection::from_mzxml_file_multistage(file),
            "mzml" => SpectrumCollection::from_mzml_file_multistage(file),
            _ => return Err(Error::UnknownFormat),
        }?;

        Ok(spectrum_collection)
    }

    /// Retrieves the spectrum with the specified scan number from this collection.
    ///
    /// If no such spectrum exists then this returns [`None`].
    /// Note, this runs in `O(log n)` and requires that the spectra are sorted by ascending scan
    /// number.
    pub fn get_scan(&self, scan: usize) -> Option<&Spectrum> {
        self.spectra
            .binary_search_by_key(&scan, |spectrum| spectrum.scan)
            .ok()
            .map(|idx| &self.spectra[idx])
    }

    /// Retrieves the spectrum with the specified scan number from this collection.
    ///
    /// This is identical to [`SpectrumCollection::get_scan`] except that it returns an error
    /// instead of [`None`] when the scan isn't found.
    ///
    /// # Errors
    /// If `scan` doesn't exist in this collection then this will always return
    /// [`Error::MissingScan`](crate::Error::MissingScan).
    pub fn try_get_scan(&self, scan: usize) -> Result<&Spectrum, crate::Error> {
        self.get_scan(scan)
            .ok_or(crate::Error::MissingScan { scan })
    }

    /// Convert this whole collection into a single scan number.
    ///
    /// If no such spectrum exists this returns [`None`]. The same runtime complexity and
    /// pre-conditions from [`SpectrumCollection::get_scan`] apply here as well.
    pub fn into_scan(self, scan: usize) -> Option<Spectrum> {
        let idx = self
            .spectra
            .binary_search_by_key(&scan, |spectrum| spectrum.scan)
            .ok()?;

        self.spectra.into_iter().nth(idx)
    }

    /// Computes a map from spectrum precursor m/zs to scans.
    ///
    /// This map will be recomputed every time this is called, so as much as possible try to re-use
    /// this structure instead of re-calling this function.
    pub fn precursors<'a>(
        &self,
        adducts: &'a [Adduct],
    ) -> BTreeMap<OrderedFloat<f64>, Vec<(&'a Adduct, usize)>> {
        let mut map = BTreeMap::default();
        for adduct in adducts {
            for spectrum in &self.spectra {
                let mass = adduct.to_mass(spectrum.pepmass);
                map.entry(mass.into())
                    .or_insert_with(Vec::new)
                    .push((adduct, spectrum.scan));
            }
        }
        map
    }
}
