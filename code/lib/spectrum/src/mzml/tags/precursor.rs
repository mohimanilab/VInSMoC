use super::cvparam::{CvParam, CvParams};
use crate::{mzml::Error, ChargeVec};
use serde::Deserialize;

const SELECTED_ION_MZ_SPECIFIER: &str = "selected ion m/z";
const CHARGE_STATE_SPECIFIER: &str = "charge state";
const PRECURSOR_INTENSITY_ACCESSION: &str = "MS:1000042";

#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct PrecursorList {
    // Though we only will use the first precursor, deserialization will fail in the case that
    // multiple are included in the file
    precursor: Vec<Precursor>,
}

impl PrecursorList {
    pub fn parse(
        &self,
        multistage: bool,
    ) -> Result<(f64, ChargeVec, Option<String>, Option<f64>), Error> {
        let precursor_id = &self.precursor[0].spectrum_ref;
        let ion = &self.precursor[0].selected_ion_list.selected_ion;

        let param = ion.param_with_name(SELECTED_ION_MZ_SPECIFIER);
        if param.is_none() && !multistage {
            return Err(Error::MissingPrecursor);
        }

        let pepmass = param.and_then(|param| param.value::<f64>().unwrap());
        if pepmass.is_none() && !multistage {
            return Err(Error::NoCvParamValue(SELECTED_ION_MZ_SPECIFIER.to_string()));
        }

        let units = param.and_then(|param| param.units());
        if (units.is_none() || units != Some("m/z")) && !multistage {
            return Err(Error::InvalidPrecursorUnits);
        }

        // Find the precursor intensity
        let precursor_intensity = ion
            .param_with_acc(PRECURSOR_INTENSITY_ACCESSION)
            .map(|x| x.value::<f64>())
            .transpose()
            .unwrap_or_default()
            .flatten();

        // Find the charge state
        let charge = ion.param_with_name(CHARGE_STATE_SPECIFIER);
        if charge.is_none() {
            return Ok((
                pepmass.unwrap_or_default(),
                ChargeVec::new(),
                precursor_id.clone(),
                precursor_intensity,
            ));
        }

        let charge = charge.unwrap().value::<i8>()?;

        let charges = charge
            .map(|x| ChargeVec::from(vec![x]))
            .unwrap_or_else(ChargeVec::new);

        Ok((
            pepmass.unwrap_or_default(),
            charges,
            precursor_id.clone(),
            precursor_intensity,
        ))
    }
}

#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
struct Precursor {
    spectrum_ref: Option<String>,
    // From the spec, a selected ion list is technically optional, but we force it to be present during
    // deserialization to make sure we get a mass to charge ratio
    selected_ion_list: SelectedIonList,
}

#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
struct SelectedIonList {
    // From the spec, we can technically have a vector of selected ions, but for our purposes, we only expect
    // a single one for any given spectrum. Thus, we only take one during deserialization.
    selected_ion: SelectedIon,
}

#[derive(Debug, Deserialize, PartialEq)]
struct SelectedIon {
    #[serde(rename = "cvParam", default)]
    cv_params: Vec<CvParam>,
}

impl CvParams for SelectedIon {
    fn params(&self) -> &Vec<CvParam> {
        &self.cv_params
    }
}
