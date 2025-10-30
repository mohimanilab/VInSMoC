mod bindata;
mod cvparam;
mod mzml;
mod precursor;
mod run;
mod scan;
mod spectrum;
mod spectrum_list;

pub use self::spectrum::Spectrum;
use cvparam::{CvParam, CvParams};
pub use mzml::{Mzml, MzmlMultistage};
use spectrum_list::SpectrumList;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MzmlDataArrayType {
    MzArray,
    IntensityArray,
    Other,
}
