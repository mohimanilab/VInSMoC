use crate::mzml::Error;
use serde::Deserialize;

use super::cvparam::{CvParam, CvParams};

const RETENTION_TIME_SPECIFIER: &str = "scan start time";

#[derive(Debug, Deserialize, PartialEq)]
pub struct ScanList {
    // From the spec, we can technically have a vector of scans, but for our purposes, we only expect
    // a single one for any given spectrum. Thus, we only take one during deserialization.
    scan: Scan,
}

#[derive(Debug, Deserialize, PartialEq)]
struct Scan {
    #[serde(rename = "cvParam", default)]
    cv_params: Vec<CvParam>,
}

impl ScanList {
    pub fn parse(&self) -> Result<Option<f64>, Error> {
        let param = self.scan.param_with_name(RETENTION_TIME_SPECIFIER);
        if param.is_none() {
            return Ok(None);
        }

        let param = param.unwrap();
        let rt = param
            .value::<f64>()?
            .ok_or_else(|| Error::NoCvParamValue(RETENTION_TIME_SPECIFIER.to_string()))?;

        let rt = match param.units() {
            Some("second") => Some(rt),
            Some("minute") => Some(rt * 60.0),
            _ => return Err(Error::InvalidRetentionTimeUnits),
        };

        Ok(rt)
    }
}

impl CvParams for Scan {
    fn params(&self) -> &Vec<CvParam> {
        &self.cv_params
    }
}
