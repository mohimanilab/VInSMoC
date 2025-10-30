use super::{CvParam, CvParams, SpectrumList};
use serde::Deserialize;

#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Run {
    #[serde(rename = "cvParam", default)]
    cv_params: Vec<CvParam>,
    pub spectrum_list: SpectrumList,
}

impl CvParams for Run {
    fn params(&self) -> &Vec<CvParam> {
        &self.cv_params
    }
}
