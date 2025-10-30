//! A crate for handling tandem mass spectral data

#![warn(missing_docs)]

/// A single spectral peak
mod peak;
pub use peak::Peak;

/// Traits for general spectra
mod traits;
pub use traits::{Ms, MsFormat, Peaks};

/// A general purpose tandem mass spectrum, not tied to any specific format.
mod spectrum;
pub use crate::spectrum::{Spectrum, SpectrumMultistageMeta};

mod tree;
pub use crate::tree::MultistageTree;

mod spectrum_graph;

pub use mgf::parsers::parse_mgf;

/// A general purpose collection of tandem mass spectra, not tied to any specific format.
mod collection;
pub use collection::SpectrumCollection;

mod errors;
pub use errors::Error;

mod parsed;
pub use parsed::{ParsedCollection, ParsedCollections, ParsedSpectrum, SpectrumPreprocessing};

pub mod splash;
pub use splash::Splash;

/// MGF spectral format parsing
///
/// This module contains types to represent MGF-encoded spectra data and functions to parse MGF
/// data.
///
/// # MGF files
/// The MGF format is an incredibly simple and easy-to-parse format.
/// The specification we followed can be found [here](https://www.matrixscience.com/help/data_file_help.html).
/// Since the specificaiton isn't really a specification but rather an informal description we also
/// referenced the [libmgf grammar](https://github.com/kirchnerlab/libmgf/blob/master/src/Parser.ypp).
///
/// The MGF format consists of several ion blocks containing spectra. Each ion block contains
/// key-value pairs from a controlled vocabulary to describe information about a single spectrum.
/// At the very top level there is a global key-value section that gives information about the
/// whole MGF file.
/// Below is a very simple example interspersed with explanation.
///
/// ```txt
/// CHARGE=1+
/// ```
/// We start with the global parameter section. We parse these values but largely never use them.
/// ```txt
/// BEGIN IONS
/// PEPMASS=200.0
/// CHARGE=1+
/// RTINSECONDS=20
/// ```
/// Next we have the start of a spectrum, marked by `BEGIN IONS`, and the local spectrum
/// parameters. Confusingly, `PEPMASS` typically refers to the precursor ion m/z, which is both not
/// necessarily a peptide and definitely not a mass.
/// ```txt
/// 10.0 21
/// 21.0 22
/// END IONS
/// ```
/// Finally we have the actual peaks, which are presented as space delimited m/z and intensity
/// pairs. The end of the spectrum is marked by `END IONS`.
/// A single MGF file will typically have many spectra and therefore many ion blocks.
///
/// # Implementation Quirks
/// This is most battle-tested parser, so it's fairly complete. However, since people have the
/// unfortunate habit of using the key-value sections of MGF files for random metadata we have
/// chosen to err on the side of tolerance and just discard any invalid key-value pairs -- as long as
/// they follow a reasonable identifier format (i.e. can't start with a number, can't contain
/// spaces) and contain an `=` -- without throwing an error.
///
/// Note since `MSLEVEL` isn't actually a standardized MGF key we cannot rely on it to select
/// tandem spectra from the input. Therefore, the [`Mgf`](mgf::Mgf) type provide by this module
/// requires every spectrum to contain a `PEPMASS` key.
pub mod mgf;

/// mzML spectral format parsing
///
/// This module contains type sto represent mzML-encoded spectral data and functions to parse mzML
/// data.
///
/// # mzML files
/// The mzML format is one of the two commonly used XML-based formats for mass spectra. The
/// specification we followed can be found [here](http://www.peptideatlas.org/tmp/mzML1.1.0.html).
/// This parser is specifically for version 1.1.0, but according to the authors' claims this should
/// be forward compatible with 1.2.X versions as well.
///
/// The mzML format consists of a single top-level mzML tag. This can contain some meta-information
/// but we only care about the `<run>` tag, which contains all of the spectra. Below is a truncated
/// example. Note the actual data is stored as two binary arrays, one for m/z and one for
/// intensity. These are optionally compressed and use different bit-precision. The precursors are
/// not theoretically required by we're only concerned with tandem mass spectra so we skip anything
/// that doesn't have an MS level of 2.
///
/// ```xml
/// <?xml version="1.0" encoding="utf-8"?>
/// <mzML
///     xmlns="http://psi.hupo.org/ms/mzml"
///     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
///     xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"
///     id="C:\Users\ane_q\Desktop\Salinispora_BRAxPTM\LC-MSMS\LCMSMS_BRA-132A_POS_B_1-2_01_9822"
///     version="1.1.0">
///   <!--
///     skipping:
///         - cvList
///         - fileDescription
///         - softwareList
///         - instrumentConfigurationList
///         - dataProcessingList
///     -->
///   <run>
///     <spectrumList count="1">
///      <spectrum index="0" id="scan=1" defaultArrayLength="68">
///        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>
///        <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>
///        <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
///        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
///        <scanList count="1">
///          <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
///          <scan instrumentConfigurationRef="_x0031_">
///            <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="3.18099" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>
///          </scan>
///        </scanList>
///        <precursorList count="1">
///          <precursor spectrumRef="scan=1">
///            <isolationWindow>
///              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.0"/>
///              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="1.0"/>
///            </isolationWindow>
///            <selectedIonList count="1">
///              <selectedIon>
///                <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="3"/>
///                <!-- this is the precursor ion -->
///                <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="352.294799804688" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
///              </selectedIon>
///            </selectedIonList>
///            <activation>
///              <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" value=""/>
///            </activation>
///          </precursor>
///        </precursorList>
///        <binaryDataArrayList count="2">
///          <binaryDataArray encodedLength="364">
///            <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>
///            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
///            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
///            <binary>EyxLQwFyS0MMnk9DvDxZQ8BWWkP+t11DTDVnQ21QZ0PXiGdD1btqQ6IybEMKK3VD/x12Q7wFeEPSt3lD96KAQ/yagUOXkoND1o+IQ+qriEPfHo1D/aKPQ82amEP0lplDpqaaQ7+mo0OmTKRDdqOkQ1C9pEMm4aRDhPykQ/UppUMmSqVDlIClQ5SnpUOkzKVD7jamQxKCpkOKAadDSBinQ0qYp0O4u6dDEp6rQ/Xur0OyGbBDiDKwQwpgsEPujrFD0g6yQ//6skPQurND+Nq0Qw53tUNOardD5xLDQ+TmxkPDks1D5ybOQwEb3EMYH95D7QLfQxPn40PuBuVDA4/yQ+HX+ENaLf5DdEsARIBHDkQ=</binary>
///          </binaryDataArray>
///          <binaryDataArray encodedLength="364">
///            <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>
///            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
///            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>
///            <binary>cJYiRsFq/0RhINhEEhY5RqB0hUS7RwhFTtlHRvfJlEQs33VDj0lPRG2zgkWWzVtGR5jZRDe0cUQsMBNFQfcZRSFO3kTiMchDLgdURVoBM0XaGn9ESdGVRegvGURPC5hFv82gRZ55kEVYR3dFbknoSFeACEdFwJhGxJQFRrFJqUV6UfJElkO7RecJHkW3YmpEndntRMPgM0XOtkFEliNIREMvxUVUGEJFI6/gRKJ8kETEV95Fk6O9RYlvskMG8eVEyFLLQ/Rpq0T/aI9EZGVTRcdy40SniuVDW64HRSwGo0S4W8RDQVqtRGThEUUqImJFhZudRFHFI0VvU25F5eGzRWKDXUTFyjJFUWooRdQ+nEU=</binary>
///          </binaryDataArray>
///        </binaryDataArrayList>
///      </spectrum>
///     </spectrumList>
///   </run>
/// </mzML>
/// ```
///
/// # Retained Information
/// We ignore a lot of the information in the mzML file. In practice we only use the peak lists,
/// precursor mass-to-charge ratio, retention time, and charge state. This means that we mostly
/// focus on the `<spectrum>` tags and ignore most of the other Controlled Vocabulary terms.
pub mod mzml;

/// mzXML spectral format parsing
///
/// This module contains types to represent mzXML-encoded spectral data and functions to parse
/// mzXML data.
///
/// # mzXML files
/// The mzXML format is one of the two commonly used XML-based formats for mass spectra.
/// The specification we followed can be found [here](http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/).
/// This parser is specifically for mzXML version 3.2, but has been tested on version 2.0 as well
/// and appears to work totally fine.
///
/// The mzXML format consists of a single top-level mzXML tag containing various descriptor tags.
/// The most important of these is the msRun tag, which contains all of our spectra.
/// Below is a truncated example.
///
/// ```xml
/// <?xml version="1.0" encoding="ISO-8859-1"?>
/// <!-- We don't parse this initial soup -->
/// <mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2"
///        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
///        xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2 http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/mzXML_idx_3.2.xsd">
///   <!-- All of our data is inside here -->
///   <!-- We can see we have 3265 scans in this file -->
///   <msRun scanCount="3265" startTime="PT2.197S" endTime="PT539.858S">
///     <!-- We have a single scan here -->
///     <!-- The relevant attribute for us the msLevel, which here says this is a MS2 spectrum -->
///     <scan num="2575"
///           scanType="Full"
///           centroided="1"
///           msLevel="2"
///           peaksCount="9"
///           polarity="+"
///           retentionTime="PT2.567S"
///           collisionEnergy="15.9"
///           lowMz="55.934370436424"
///           highMz="118.08568574113"
///           basePeakMz="118.083830426663"
///           basePeakIntensity="661"
///           totIonCurrent="7997"
///           msInstrumentID="1">
///       <!-- We can get the precursor ion here -->
///       <!-- This isn't technically required by the spec, but we require it -->
///       <precursorMz precursorScanNum="2405" precursorIntensity="7997.0" activationMethod="HCD">118.086672610328</precursorMz>
///       <!-- The peak data is a single base64-encoded string -->
///       <!-- In this case we see that it's uncompressed and contains 32-bit floats -->
///       <!-- The encoded data is pairs of floats, alternating b/w m/z and intensity -->
///       <peaks compressionType="none"
///              compressedLen="0"
///              precision="32"
///              byteOrder="network"
///              contentType="m/z-int">
///         Ql+8zEMpU1lCaEJoQ+Ri+kJsSt5Dcyd2Qo/b50LcwrJCkd/oQ64sv0KZ3mVCnUgAQufr6UIm8GRC696JQzZy+0LsK99EKjK3
///       </peaks>
///     </scan>
///   </msRun>
/// </mzXML>
/// ```
///
/// Since this crate is only concerned with *tandem* mass spectra the [`MzXml`](mzxml::MzXml) struct
/// provided by this module will only parse spectra from the input with an MS level of 2. All other
/// spectra are skipped.
///
/// # Implementation Quirks
/// Our mzXML parsing is known to be incomplete. There are many cases where an invalid mzXML file
/// wil pass our parser without an error. Here's a short list of things we don't check.
/// 1. the index
/// 2. parentFile tags
/// 3. msInstrument tags
/// 4. dataProcessing tags, except to make sure that one of them says the file is centroided
/// 5. separation tags
/// 6. spotting tags
/// 7. sha1 tags
/// 8. scanOrigin tags
/// 9. maldi tags
/// 10. any attributes we won't use (polarity, scanType, startMz, etc.)
///
/// We also might fail on a technically valid mzXML file:
/// - If an msLevel 2 scan doesn't have a precursorMz tag we count that as a failure
///
/// # Centroiding
/// We don't have a centroiding implementation, so we recommend that the input data is already
/// centroided. We used to throw an error on un-centroided data, but we now no longer do so.
/// However, it's still recommended that any input data is centroided. One option for this is using
/// [proteowizard](https://github.com/ProteoWizard/pwiz).
pub mod mzxml;

/// Shared parsing for mzML and mzXML files
mod xml;

pub mod draw;

#[cfg(test)]
mod tests;

/// A vector of charges.
///
/// See [`SmallVec`](smallvec::SmallVec) for construction options.
pub type ChargeVec = smallvec::SmallVec<[i8; 2]>;
