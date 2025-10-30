use std::{
    ops::{Deref, DerefMut},
    str::FromStr,
};

use crate::{
    spectrum::SpectrumMultistageMeta,
    traits::MsFormat,
    xml::{decode_data_array, Endianness},
    ChargeVec, Peak,
};
use roxmltree::{Document, Node};

mod errors;
pub use errors::Error;

/// A single scan from an mzXML file.
///
/// This can only be constructed from within a [`MzXml`].
#[derive(Debug, Clone)]
pub struct Spectrum {
    scan_num: usize,
    precursor: f64,
    precursor_intensity: Option<f64>,
    precursor_scan: Option<usize>,
    ms_level: Option<usize>,
    charges: ChargeVec,
    retention_time: Option<f64>,
    peaks: Vec<Peak>,
}

impl Spectrum {
    /// Retrieve the precursor mz
    pub fn precursor(&self) -> f64 {
        self.precursor
    }

    /// Retrieve the retention time
    pub fn retention_time(&self) -> Option<f64> {
        self.retention_time
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

impl From<Spectrum> for crate::Spectrum {
    fn from(spectrum: Spectrum) -> Self {
        Self {
            peaks: spectrum.peaks,
            charges: spectrum.charges,
            pepmass: spectrum.precursor,
            retention_time: spectrum.retention_time.map(|x| x as f64),
            multistage_meta: SpectrumMultistageMeta {
                ms_level: spectrum.ms_level,
                precursor_intensity: spectrum.precursor_intensity,
                precursor_scan: spectrum.precursor_scan,
            },
            scan: spectrum.scan_num,
            denormalize_factor: None,
        }
    }
}

/// A fully parsed mzXML file.
///
/// The parsed version will only contain spectra with msLevel 2, meaning we discard all spectra
/// that aren't tandem mass spectra.
///
/// This should be constructed using its [`FromStr`] or [`MsFormat`] implementations.
#[derive(Debug)]
pub struct MzXml {
    /// The spectrum that make up this mzXML file
    pub spectra: Vec<Spectrum>,
}

impl MzXml {
    /// Returns the spectrum with a given scan number.
    ///
    /// Returns [`None`] if the scan number is not found.
    pub fn get_by_scan(&self, scan: usize) -> Option<&Spectrum> {
        self.spectra
            .iter()
            .find(|spectrum| spectrum.scan_num == scan)
    }
}

impl From<MzXml> for crate::SpectrumCollection {
    fn from(mzxml: MzXml) -> Self {
        let mut spectra = mzxml
            .spectra
            .into_iter()
            .map(crate::Spectrum::from)
            .collect::<Vec<_>>();
        spectra.sort_by_key(|spectrum| spectrum.scan);

        Self { spectra }
    }
}

/// Retrieve a required attribute given its name from the provided XML node.
fn required_attr<'a>(node: &'a Node, attr_name: &'static str) -> Result<&'a str, Error> {
    node.attribute(attr_name)
        .ok_or_else(|| Error::MissingAttribute {
            attr_name,
            tag: node.tag_name().name().to_string(),
        })
}

/// Retrieve required text from the given node.
fn required_text<'a>(node: &'a Node) -> Result<&'a str, Error> {
    node.text().ok_or_else(|| Error::MissingText {
        tag: node.tag_name().name().to_string(),
    })
}

fn validate_and_parse(mzxml: Node, only_ms2: bool) -> Result<MzXml, Error> {
    if !mzxml.has_tag_name("mzXML") {
        return Err(Error::NoMzXmlTag);
    }

    let version = {
        let xmlns = mzxml.namespaces()[0].uri();
        let (_, version) = xmlns.rsplit_once('_').ok_or(Error::BadVersion)?;
        version
            .parse::<f32>()
            .map_err(|_| Error::BadVersion)?
            .trunc() as u32
    };

    if version != 2 && version != 3 {
        return Err(Error::BadVersion);
    }

    let msrun = mzxml.first_element_child().unwrap();
    if !msrun.has_tag_name("msRun") {
        return Err(Error::NoMsRunTag);
    }

    // use the required scan count attr to pre-allocate the spectrum vector
    let scan_cnt = required_attr(&msrun, "scanCount")?.parse::<usize>()?;
    let mut spectra = Vec::with_capacity(scan_cnt);

    // parse each scan
    for scan in msrun
        .descendants()
        .filter(|child| child.has_tag_name("scan"))
    {
        // skip all the non MS/MS spectra unless `only_ms2` is false
        let ms_level = required_attr(&scan, "msLevel")?.parse::<usize>()?;
        if only_ms2 && ms_level != 2 {
            continue;
        }

        // the scan number is required
        let scan_num = required_attr(&scan, "num")?.parse::<usize>()?;

        // the precursor tag isn't strictly required, but we will require it here for msLevel 2
        // if `only_ms2` is true.
        let precursor = scan.children().find(|x| x.has_tag_name("precursorMz"));
        if precursor.is_none() && only_ms2 {
            return Err(Error::MissingPrecursorMz(scan_num));
        }
        let precursor_intensity = precursor.and_then(|precursor| {
            precursor
                .attribute("precursorIntensity")
                .map(|intensity| intensity.parse::<f64>())
                .transpose()
                .unwrap()
        });
        let charges = precursor
            .and_then(|precursor| {
                precursor
                    .attribute("precursorCharge")
                    .map(|charge| charge.parse::<i8>())
                    .transpose()
                    .unwrap()
                    .map(|charge| ChargeVec::from(vec![charge]))
            })
            .unwrap_or_default();
        let precursor_scan = precursor.and_then(|precursor| {
            precursor
                .attribute("precursorScanNum")
                .map(|scan_num| scan_num.parse::<usize>())
                .transpose()
                .unwrap()
        });
        let precursor = precursor
            .map(|precursor| required_text(&precursor).unwrap().parse::<f64>().unwrap())
            .unwrap_or_default();

        // get the retention time if it exists
        let mut retention_time = None;
        if let Some(rt) = scan.attribute("retentionTime") {
            if !rt.starts_with("PT") || !rt.ends_with('S') {
                return Err(Error::InvalidRetentionTime(scan_num));
            }
            // skip the PT and S
            let rt = &rt[2..rt.len() - 1];
            retention_time = Some(rt.parse::<f64>()?);
        }

        // now we can parse the peaks
        // first get the tag
        let peaks_node = scan
            .children()
            .find(|x| x.has_tag_name("peaks"))
            .ok_or(Error::MissingPeaks(scan_num))?;

        // the contentType attr is required but is only allowed to be m/z-int
        // however, there's lots of files that don't include this, so we're skipping it

        // use the peak count to pre-allocate the peak vector
        let peak_cnt = required_attr(&scan, "peaksCount")?.parse::<usize>()?;
        let mut peaks = Vec::with_capacity(peak_cnt);

        // check the compression
        let compression = if version == 3 {
            let compression = required_attr(&peaks_node, "compressionType").unwrap_or("none");
            if compression != "none" && compression != "zlib" {
                return Err(Error::InvalidCompressionType {
                    scan_num,
                    value: compression.to_string(),
                });
            }
            compression == "zlib"
        } else {
            false
        };

        // check the precision, error handling for this will be done by fn call
        let precision = required_attr(&peaks_node, "precision")?.parse::<usize>()?;

        // get the actual base64-encoded data, skipping empty spectra
        let data = peaks_node.text();
        if let Some(data) = data {
            let data = decode_data_array(data.as_bytes(), precision, compression, Endianness::Big)
                .map_err(|source| Error::DecodePeaks { scan_num, source })?;
            for mz_int in data.chunks(2) {
                peaks.push(Peak::new(mz_int[0], mz_int[1]));
            }
        }

        // just in case
        peaks.shrink_to_fit();

        spectra.push(Spectrum {
            scan_num,
            precursor,
            ms_level: Some(ms_level),
            precursor_intensity,
            precursor_scan,
            charges,
            retention_time,
            peaks,
        });
    }

    // we over-estimated the # of spectra
    // because there may have been some MS1 spectra we skipped
    spectra.shrink_to_fit();

    Ok(MzXml { spectra })
}

impl FromStr for MzXml {
    type Err = Error;

    fn from_str(contents: &str) -> Result<Self, Self::Err> {
        let document = Document::parse(contents)?;

        // get to required start point
        let mzxml = document.root_element();
        validate_and_parse(mzxml, true)
    }
}

impl MsFormat for MzXml {}

/// A fully parsed mzXML file appropriate for multistage processing.
///
/// The parsed version will contain all spectra regardless of msLevel.
///
/// This should be constructed using its [`FromStr`] or [`MsFormat`] implementations.
pub struct MzXmlMultistage {
    /// The spectra that make up this mzXML file
    pub spectra: Vec<Spectrum>,
}

impl FromStr for MzXmlMultistage {
    type Err = Error;

    fn from_str(contents: &str) -> Result<Self, Self::Err> {
        let document = Document::parse(contents)?;

        // get to required start point
        let mzxml = document.root_element();
        validate_and_parse(mzxml, false).map(|mzxml| mzxml.into())
    }
}

impl From<MzXmlMultistage> for MzXml {
    fn from(mzxml_multistage: MzXmlMultistage) -> Self {
        Self {
            spectra: mzxml_multistage.spectra,
        }
    }
}

impl From<MzXml> for MzXmlMultistage {
    fn from(mzxml: MzXml) -> Self {
        Self {
            spectra: mzxml.spectra,
        }
    }
}

impl From<MzXmlMultistage> for crate::SpectrumCollection {
    fn from(mzxml_multistage: MzXmlMultistage) -> Self {
        MzXml::from(mzxml_multistage).into()
    }
}

impl MsFormat for MzXmlMultistage {}

#[cfg(test)]
mod tests;
