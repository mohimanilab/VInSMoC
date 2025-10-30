use serde::Deserialize;

use super::{cvparam::CvParam, MzmlDataArrayType};
use crate::{
    mzml::Error,
    xml::{decode_data_array, Endianness},
    Peak,
};

const FLOAT_32_SPECIFIER: &str = "32-bit float";
const FLOAT_64_SPECIFIER: &str = "64-bit float";
const MZ_ARRAY_SPECIFIER: &str = "m/z array";
const INTENSITY_ARRAY_SPECIFIER: &str = "intensity array";
const ZLIB_COMPRESSION_SPECIFIER: &str = "zlib compression";
const NO_COMPRESSION_SPECIFIER: &str = "no compression";

#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct BinDataList {
    #[serde(rename = "binaryDataArray", default)]
    binary_data_arrays: Vec<BinDataArray>,
}

impl BinDataList {
    pub fn parse(&self) -> Result<Vec<Peak>, Error> {
        let mut mz_array = None;
        let mut intensity_array = None;

        // go through all the binary data arrays
        // for the m/z and intensity arrays store them inside this spectrum
        for bin_data_array in self.binary_data_arrays.iter() {
            let (array_type, array) = bin_data_array.decode()?;
            match array_type {
                MzmlDataArrayType::Other => (),
                MzmlDataArrayType::MzArray => mz_array = Some(array),
                MzmlDataArrayType::IntensityArray => intensity_array = Some(array),
            }

            // just leave when we find what we want
            if mz_array.is_some() && intensity_array.is_some() {
                break;
            }
        }

        // Verify that all needed arrays were found
        let mz_array = mz_array.ok_or(Error::NoMzArray)?;
        let intensity_array = intensity_array.ok_or(Error::NoIntensityArray)?;

        let peaks = mz_array
            .into_iter()
            .zip(intensity_array.into_iter())
            .map(Peak::from)
            .collect();

        Ok(peaks)
    }
}

#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
struct BinDataArray {
    binary: String,
    #[serde(rename = "cvParam", default)]
    cv_params: Vec<CvParam>,
}

impl BinDataArray {
    fn decode(&self) -> Result<(MzmlDataArrayType, Vec<f64>), Error> {
        // Flags to tell us whether we found certain attributes
        let mut found_float_type = false;
        let mut found_array_type = false;
        let mut found_compress = false;

        // The actual values of the above attributes
        let mut float_type = 64usize;
        let mut use_zlib = false;
        let mut array_type = MzmlDataArrayType::Other;

        for param in self.cv_params.iter() {
            let name = param.name.as_str();
            // first check if we've already found this
            match name {
                FLOAT_32_SPECIFIER | FLOAT_64_SPECIFIER => {
                    if found_float_type {
                        return Err(Error::MultipleFloat);
                    }
                    found_float_type = true;
                }
                ZLIB_COMPRESSION_SPECIFIER | NO_COMPRESSION_SPECIFIER => {
                    if found_compress {
                        return Err(Error::MultipleCompression);
                    }
                    found_compress = true;
                }
                MZ_ARRAY_SPECIFIER | INTENSITY_ARRAY_SPECIFIER => {
                    if found_array_type {
                        return Err(Error::MultipleArrayDatatype);
                    }
                    found_array_type = true;
                }
                _ => (),
            }

            // now we know we don't have an error, so handle what we found
            match name {
                FLOAT_32_SPECIFIER => float_type = 32usize,
                FLOAT_64_SPECIFIER => float_type = 64usize,
                ZLIB_COMPRESSION_SPECIFIER => use_zlib = true,
                NO_COMPRESSION_SPECIFIER => use_zlib = false,
                MZ_ARRAY_SPECIFIER => array_type = MzmlDataArrayType::MzArray,
                INTENSITY_ARRAY_SPECIFIER => array_type = MzmlDataArrayType::IntensityArray,
                _ => (),
            }
        }

        // Check to make sure that everything required was found
        if !found_float_type {
            return Err(Error::NoFloat);
        }

        if !found_compress {
            return Err(Error::NoCompression);
        }

        if array_type == MzmlDataArrayType::Other {
            return Ok((MzmlDataArrayType::Other, Vec::new()));
        }

        Ok((
            array_type,
            decode_data_array(
                self.binary.as_bytes(),
                float_type,
                use_zlib,
                Endianness::Little,
            )?,
        ))
    }
}
