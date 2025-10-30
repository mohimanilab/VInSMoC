use std::{
    io,
    num::{ParseFloatError, ParseIntError},
};

use thiserror::Error;

use crate::xml::DecodeDataError;

/// An error that can occur when parsing an mzXML file.
///
/// Our validation is known to be incomplete. See the [module level docs](crate::mzxml) for more
/// information.
#[derive(Error, Debug)]
pub enum Error {
    /// The input XML file wasn't valid XML
    #[error("couldn't parse general XML structure")]
    Xml(#[from] roxmltree::Error),
    /// The input didn't contain an mzXML tag
    #[error("couldn't find required <mzXML> tag")]
    NoMzXmlTag,
    /// The input had an invalid `xmlns` version string
    #[error("xmlns version string on mzXML tag was invalid")]
    BadVersion,
    /// The input didn't contain an msRun tag
    #[error("couldn't find required <msRun> tag")]
    NoMsRunTag,
    /// The spectra were not centroided
    #[error("spectra are not centroided")]
    Centroiding,
    /// An error doing IO
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Error parsing integer attribute value
    #[error("unable to parse integer-valued attribute")]
    ParseInt(#[from] ParseIntError),
    /// Error parsing integer attribute value
    #[error("unable to parse float-valued attribute")]
    ParseFloat(#[from] ParseFloatError),
    /// Invalid retention time
    #[error("scan {0} had an invalid retention time specifier")]
    InvalidRetentionTime(usize),
    /// A required attribute was missing
    #[error("required attribute {attr_name} on tag {tag} was missing")]
    MissingAttribute {
        /// Name of the attribute we needed but didn't find
        attr_name: &'static str,
        /// Name of the tag we wanted that attribute for
        tag: String,
    },
    /// A tag that requires text didn't have that text
    #[error("tag {tag} didn't have its required text")]
    MissingText {
        /// Name of the tag
        tag: String,
    },
    /// An MS/MS scan didn't have a precursorMz tag
    #[error("scan number {0} had msLevel 2 but no <precursorMz> child")]
    MissingPrecursorMz(usize),
    /// An MS/MS scan didn't have a peaks tag
    #[error("scan number {0} had no <peaks> child")]
    MissingPeaks(usize),
    /// A peaks tag had an invalid compressionType attribute
    #[error("scan number {scan_num} had invalid compressionType {value}")]
    InvalidCompressionType {
        /// The failing scan number
        scan_num: usize,
        /// The compressionType value we actually received
        value: String,
    },
    /// Scan had an invalid contentType attribute
    #[error("scan number {0} <peak> contentType wasn't \"m/z-int\"")]
    InvalidContentType(usize),
    /// Unable to decode the actual base64 peak data
    #[error("scan number {scan_num} had un-decodable <peak> data")]
    DecodePeaks {
        /// The failing scan number
        scan_num: usize,
        #[source]
        /// The cause of this error
        source: DecodeDataError,
    },
    /// Unable to get peaks from the decoded base64 data
    #[error("scan number {scan_num} had {len} floats, this needs to be an even number")]
    UnevenDecodedFloatCount {
        /// The failing scan number
        scan_num: usize,
        /// The number of floats found in the decoded data
        len: usize,
    },
}
