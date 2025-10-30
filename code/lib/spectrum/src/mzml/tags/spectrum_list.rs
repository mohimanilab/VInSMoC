use super::Spectrum;
use serde::Deserialize;

#[derive(Debug, Deserialize, PartialEq)]
pub struct SpectrumList {
    #[serde(rename = "spectrum", default)]
    pub spectra: Vec<Spectrum>,
}
