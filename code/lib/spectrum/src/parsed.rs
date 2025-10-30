use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::PathBuf,
    str::FromStr,
};

use paths::FilePath;
use serde::ser::SerializeStruct;

use crate::{Peaks, SpectrumCollection};

/// Options for pre-processing mass spectra after they were first parsed.
///
/// This goes through five main steps when pre-processing spectra:
/// 1. Discard a spectrum if it doesn't have at least [`SpectrumPreprocessing::min_peak_count`]
///    peaks.
/// 2. Merge peaks with [`Peaks::merge_peaks`] unless [`SpectrumPreprocessing::no_merge`] is set.
/// 3. Filter peaks with [`Peaks::filter_peaks`] unless [`SpectrumPreprocessing::no_filter`] is
///    set.
/// 4. Applies [`Peaks::filter_top_n_spectrum_peaks`] for each scan, if `max_num_peaks` is set.
/// 5. Sort peaks by increasing m/z.
///
/// There is an additional fifth step, disabled by default, which is peak normalization - this normalizes
/// the peak intensity vector so that its magnitude is 1.0.
///
/// These steps are encapsulated in [`SpectrumPreprocessing::preprocess`].
///
/// # [`clap`]
/// This implements [`clap::Parser`] so it can be used to parse command line arguments.
/// ```
/// use clap::Parser;
/// use spectrum::SpectrumPreprocessing;
///
/// #[derive(Debug, Clone, Default, Parser)]
/// #[clap(name = "bin")]
/// struct Opts {
///     #[clap(flatten)]
///     pre: SpectrumPreprocessing,
/// }
///
/// // long arguments
/// let args = "bin --peak-merge-thresh 0.05 --min-peak-count 1 --peak-filter-window 50.0 --peaks-per-window 5".split_whitespace();
/// let params = Opts::try_parse_from(args).unwrap();
/// assert_eq!(params.pre, SpectrumPreprocessing {
///     peak_merge_thresh: 0.05,
///     min_peak_count: 1,
///     peak_filter_window: 50.0,
///     peaks_per_window: 5,
///     max_num_peaks: None,
///     no_merge: false,
///     no_filter: false,
///     normalize: false,
///     masst_consistent: false,
///     multistage: false,
/// });
///
/// // this has no short arguments
///
/// // --no-merge conflicts with merge parameters
/// let args = "bin --peak-merge-thresh 100.0 --no-merge".split_whitespace();
/// assert!(Opts::try_parse_from(args).is_err());
///
/// // --no-filter conflicts with filter parameters
/// let args = "bin --peak-filter-window 50.0 --peaks-per-window 5 --no-filter".split_whitespace();
/// assert!(Opts::try_parse_from(args).is_err());
///
/// // default CLI args match `Default` implementation
/// let args = "bin".split_whitespace();
/// let params = Opts::try_parse_from(args).unwrap();
/// assert_eq!(params.pre, SpectrumPreprocessing::default());
/// assert_eq!(SpectrumPreprocessing::default(), SpectrumPreprocessing {
///     peak_merge_thresh: 0.05,
///     min_peak_count: 1,
///     peak_filter_window: 50.0,
///     peaks_per_window: 5,
///     max_num_peaks: None,
///     no_merge: false,
///     no_filter: false,
///     normalize: false,
///     masst_consistent: false,
///     multistage: false,
/// });
/// ```
#[derive(Debug, Clone, PartialEq, serde::Serialize, clap::Parser)]
#[clap(help_heading = "Spectra")]
pub struct SpectrumPreprocessing {
    /// Threshold used to merge peaks (Da).
    ///
    /// See [`Peaks::merge_peaks`].
    #[clap(
        long,
        default_value = "0.05",
        conflicts_with = "no-merge",
        long_help = "Threshold used to merge peaks (Da)\n\nPeaks with m/zs within this value of each other will be merged into a single peak with intensity aggregated additively. The minimum number of peaks that can remain after the merge pre-processing step is defined by --min-peak-count. Merging pre-preprocessing can be skipped with --no-merge."
    )]
    pub peak_merge_thresh: f64,
    /// Minimum number of spectral peaks.
    ///
    /// Any spectra with fewer peaks than this will be filtered out. See also [`Peaks::merge_peaks`].
    #[clap(
        long,
        default_value = "1",
        long_help = "Minimum number of spectral peaks\n\nAny spectra with fewer peaks than this will be filtered out. This also serves as the minimum number of peaks that can be returned after the peak merging operation."
    )]
    pub min_peak_count: usize,
    /// Window size for peak filtering (Da).
    ///
    /// See [`Peaks::filter_peaks`].
    #[clap(
        long,
        default_value = "50.0",
        conflicts_with = "no-filter",
        long_help = "Window size for peak filtering (Da)\n\nThis value is used in the filtering pre-processing step of the spectrum. This removes low intensity peaks by keeping only the --peaks-per-window most intense peaks in non-overlapping m/z windows of size --peak-filter-window. Filtering pre-preprocessing can be skipped with --no-filter."
    )]
    pub peak_filter_window: f64,
    /// Number of peaks to keep per window
    ///
    /// See [`Peaks::filter_peaks`].
    #[clap(
        long,
        default_value = "5",
        conflicts_with = "no-filter",
        long_help = "Number of peaks to keep per window\n\nThis value is used in the filtering pre-processing step of the spectrum. This removes low intensity peaks by keeping only the --peaks-per-window most intense peaks in non-overlapping m/z windows of size --peak-filter-window. Filtering pre-preprocessing can be skipped with --no-filter."
    )]
    pub peaks_per_window: usize,
    /// Number of top-intensity peaks per spectrum to keep after peak filtering.
    ///
    /// See [`Peaks::filter_top_n_spectrum_peaks`].
    #[clap(long)]
    pub max_num_peaks: Option<usize>,
    /// Disable peak merging
    #[clap(long)]
    pub no_merge: bool,
    /// Disable peak filtering
    #[clap(long)]
    pub no_filter: bool,
    /// Enable peak intensity normalization
    #[clap(long)]
    pub normalize: bool,
    /// Make preprocessing consistent with MASST, disregarding other spectral processing parameters.
    #[clap(long)]
    pub masst_consistent: bool,
    /// Enable multistage parsing of the spectrum.
    #[clap(skip)]
    pub multistage: bool,
}

impl Default for SpectrumPreprocessing {
    fn default() -> Self {
        Self {
            peak_merge_thresh: 0.05,
            min_peak_count: 1,
            peak_filter_window: 50.0,
            peaks_per_window: 5,
            max_num_peaks: None,
            no_merge: false,
            no_filter: false,
            normalize: false,
            masst_consistent: false,
            multistage: false,
        }
    }
}

impl SpectrumPreprocessing {
    /// Pre-process the provided spectrum.
    ///
    /// Steps in pre-processing are described in [the struct level docs](SpectrumPreprocessing). If
    /// the provided `spectrum` would have been removed during pre-processing this will return
    /// [`None`]. Otherwise, the provided `spectrum` will have the described pre-processing steps
    /// applied to it and then returned as [`Some`]. The pre-processing will happen on the
    /// parameter in place, the signature is to make it easy to filter out spectra using this
    /// method.
    pub fn preprocess(&self, mut spectrum: crate::Spectrum) -> Option<crate::Spectrum> {
        if spectrum.len() < self.min_peak_count {
            return None;
        }

        if self.masst_consistent {
            spectrum.filter_peaks_masst(50.0, 5); // tyler hard coded this?
            spectrum.denormalize_factor = spectrum.normalize_peaks();
            spectrum
                .peaks_mut()
                .sort_unstable_by(|x, y| x.mz.partial_cmp(&y.mz).unwrap());
            return Some(spectrum);
        }

        if !self.no_merge {
            spectrum.merge_peaks(self.peak_merge_thresh, self.min_peak_count);
        }

        if !self.no_filter {
            spectrum.filter_peaks(self.peak_filter_window, self.peaks_per_window);
        }

        if self.normalize {
            spectrum.denormalize_factor = spectrum.normalize_peaks();
        }

        if let Some(num_peaks) = self.max_num_peaks {
            spectrum.filter_top_n_spectrum_peaks(num_peaks);
        }

        spectrum
            .peaks_mut()
            .sort_unstable_by(|x, y| x.mz.partial_cmp(&y.mz).unwrap());

        Some(spectrum)
    }
}

/// A single spectrum collection that has had pre-processing functions applied.
///
/// This wraps a [`SpectrumPreprocessing`] and a path to a spectrum file on disk. The final
/// preprocessed spectrum collection can be retrieved with [`ParsedCollection::collection`].
///
/// See [`SpectrumPreprocessing`] for pre-processing steps. In addition to these steps this sorts
/// all spectra in the collection by ascending scan.
#[derive(Debug, Clone, clap::Parser)]
pub struct ParsedCollection {
    /// Path to a file containing the spectrum collection of interest.
    ///
    /// These are parsed using [`SpectrumCollection::from_path`].
    #[clap(
        short = 's',
        long = "spectra-file",
        long_help = "Path to file containing the spectrum collection of interest."
    )]
    pub path: FilePath,
    /// The pre-processing parameters.
    #[clap(flatten)]
    pub pre: SpectrumPreprocessing,
}

impl ParsedCollection {
    /// Parse a spectrum collection from [`ParsedCollection::path`] and apply pre-processing
    /// according to parameters in [`ParsedCollection::pre`].
    ///
    /// # Errors
    /// If there was an error in parsing the spectrum this returns a [`crate::Error`]. See
    /// [`SpectrumCollection::from_path`].
    pub fn collection(&self) -> Result<SpectrumCollection, crate::Error> {
        let mut collection = if self.pre.multistage {
            SpectrumCollection::from_path_multistage(&self.path)?
        } else {
            SpectrumCollection::from_path(&self.path)?
        };
        collection.preprocess(&self.pre);
        Ok(collection)
    }
}

/// A single scan from a collection that has had pre-processing functions applied.
///
/// This wraps a [`ParsedCollection`] and a scan that exists in that collection. The final
/// preprocessed spectrum can be retrieved with the [`TryFrom`] implementation on
/// [`crate::Spectrum`]. See [`SpectrumPreprocessing`] for pre-processing steps and
/// [`ParsedCollection`] for information about the parsing process.
#[derive(Debug, Clone, clap::Parser)]
pub struct ParsedSpectrum {
    /// Collection that the target scan exists in.
    #[clap(flatten)]
    pub collection: ParsedCollection,
    /// Scan of target spectrum.
    #[clap(long, help_heading = "Spectra")]
    pub scan: usize,
}

impl TryFrom<&ParsedSpectrum> for crate::Spectrum {
    type Error = crate::Error;

    fn try_from(value: &ParsedSpectrum) -> Result<Self, Self::Error> {
        let collection = value.collection.collection()?;
        let scan = collection.into_scan(value.scan);
        let scan = scan.ok_or(crate::Error::MissingScan { scan: value.scan })?;
        Ok(scan)
    }
}

/// An iterator over spectra that have had pre-processing functions applied.
///
/// Internally this merges similar peaks and filters out low intensity peaks using
/// [`Peaks::merge_peaks`] and [`Peaks::filter_peaks`] respectively. It yields both the original
/// input file path of a spectrum file and the preprocessed [`SpectrumCollection`].
///
/// For each collection it
/// 1. Removes scans with fewer than `min_peak_count` peaks.
/// 2. Applies [`Peaks::merge_peaks`] for each scan, unless `no_merge` is set
/// 3. Applies [`Peaks::filter_peaks`] for each scan, unless `no_filter` is set
/// 4. Applies [`Peaks::filter_top_n_spectrum_peaks`] for each scan, if `max_num_peaks` is set.
/// 5. Sorts all peaks in each scan by ascending m/z
/// 6. Sorts collection by ascending scan numbers
///
/// The [`Default`] implementation contains no paths and typically will not be useful.
#[derive(Debug, Default, Clone, clap::Parser)]
#[clap(help_heading = "Spectra")]
pub struct ParsedCollections {
    /// Paths to files containing the spectra of interest.
    ///
    /// These are parsed using [`SpectrumCollection::from_path`]. Text files are treated as newline separated lists of spectra.
    #[clap(
        short = 's',
        required = true,
        multiple = true,
        long = "spectra-files",
        long_help = "Paths to files containing the spectra of interest\n\nTo use more than one spectrum add multiple, space-separated files\nNote that .txt files will be treated as newline separated lists of spectra"
    )]
    pub paths: Vec<FilePath>,
    /// The pre-processing parameters
    #[clap(flatten)]
    pub pre: SpectrumPreprocessing,
    /// The current index of this iterator.
    #[clap(skip = 0usize)]
    idx: usize,
    /// Tracker for when .txt files have been expanded.
    #[clap(skip = false)]
    expanded: bool,
}

impl ParsedCollections {
    /// Parse the list of file paths expanding any .txt files to their constituent spectra files.
    ///
    /// This returns a new list of paths that replaces any .txt files with their line contents.
    /// This operation only needs to be done once per [`ParsedCollections`] instance, which is
    /// tracked via the [`ParsedCollections::expanded`] field.
    ///
    /// # Errors
    /// This will resolve into an error if any of the following conditions occur:
    /// - A path in `paths` has no extension
    /// - A path in `paths` has an extension that isn't valid unicode
    /// - A path in `paths` with a "txt" extension cannot be opened with [`File::open`]
    fn expand_file_paths(paths: Vec<FilePath>) -> Result<Vec<FilePath>, crate::Error> {
        paths
            .into_iter()
            .flat_map(|path| {
                // pull out extension
                let ext = path.extension().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "spectrum file had no extension",
                    )
                });
                if ext.is_err() {
                    return vec![ext.map(|_| path.clone()).map_err(|e| e.into())];
                }
                let ext = ext.unwrap();

                // convert extension to str
                let ext = ext.to_str().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "file extension wasn't valid unicode",
                    )
                });
                if ext.is_err() {
                    return vec![ext.map(|_| path.clone()).map_err(|e| e.into())];
                }
                let ext = ext.unwrap();

                if ext == "txt" {
                    // open txt file
                    let file = File::open(&path);
                    if file.is_err() {
                        return vec![file.map(|_| path.clone()).map_err(|e| e.into())];
                    }
                    let file = file.unwrap();

                    // expand lines into their own paths
                    let buf = BufReader::new(file);
                    buf.lines()
                        .map(|line| {
                            line.and_then(|line| FilePath::from_str(&line))
                                .map_err(|e| e.into())
                        })
                        .collect()
                } else {
                    // otherwise just the original path
                    vec![Ok(path)]
                }
            })
            .collect()
    }

    /// Construct a new [`ParsedCollections`] for the given `paths` using the provided
    /// pre-processing parameters.
    pub fn new(paths: Vec<FilePath>, pre: SpectrumPreprocessing) -> Self {
        Self {
            paths,
            pre,
            idx: 0,
            expanded: false,
        }
    }

    /// Resets the iterator.
    pub fn reset(&mut self) {
        self.idx = 0;
    }
}

impl Iterator for ParsedCollections {
    type Item = Result<(PathBuf, SpectrumCollection), crate::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.expanded {
            match ParsedCollections::expand_file_paths(self.paths.clone()) {
                Ok(paths) => {
                    self.paths = paths;
                    self.expanded = true;
                }
                e => {
                    return Some(
                        e.map(|_| (PathBuf::from(""), SpectrumCollection { spectra: vec![] })),
                    )
                }
            };
        }

        if self.idx >= self.paths.len() {
            return None;
        }

        let mut collection = if self.pre.multistage {
            SpectrumCollection::from_path_multistage(&self.paths[self.idx])
        } else {
            SpectrumCollection::from_path(&self.paths[self.idx])
        };

        if let Ok(inner) = collection.as_mut() {
            inner.preprocess(&self.pre);
        }
        let output =
            Some(collection.map(|collection| (self.paths[self.idx].clone().into(), collection)));
        self.idx += 1;
        output
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.paths.len() - self.idx;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for ParsedCollections {}

impl FromIterator<FilePath> for ParsedCollections {
    fn from_iter<T: IntoIterator<Item = FilePath>>(iter: T) -> Self {
        Self {
            paths: iter.into_iter().collect(),
            pre: SpectrumPreprocessing::default(),
            idx: 0,
            expanded: false,
        }
    }
}

impl serde::Serialize for ParsedCollections {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut s = serializer.serialize_struct("ParsedCollections", 4)?;
        s.serialize_field(
            "paths",
            &self
                .paths
                .iter()
                .map(|path| path.to_str().unwrap())
                .collect::<Vec<&str>>(),
        )?;
        s.serialize_field("pre", &self.pre)?;
        s.serialize_field("idx", &self.idx)?;
        s.serialize_field("expanded", &self.expanded)?;
        s.end()
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use clap::Parser;

    use super::*;

    #[test]
    fn opts() {
        // default params
        let args = "fakebin -s src/lib.rs".split_whitespace();
        let opts = ParsedCollections::try_parse_from(args).unwrap();
        assert_eq!(opts.pre.peak_merge_thresh, 0.05);
        assert_eq!(opts.pre.min_peak_count, 1);
        assert_eq!(opts.pre.peak_filter_window, 50.0);
        assert_eq!(opts.pre.peaks_per_window, 5);
        assert_eq!(&opts.paths[0], &Path::new("src/lib.rs"));
        assert_eq!(opts.pre.no_merge, false);
        assert_eq!(opts.pre.no_filter, false);
        assert_eq!(opts.pre.normalize, false);

        // disable merging
        let args = "fakebin -s src/lib.rs --no-merge --no-filter --normalize".split_whitespace();
        let opts = ParsedCollections::try_parse_from(args).unwrap();
        assert_eq!(opts.pre.no_merge, true);
        assert_eq!(opts.pre.no_filter, true);
        assert_eq!(opts.pre.normalize, true);

        // missing file
        let args = "fakebin".split_whitespace();
        let opts = ParsedCollections::try_parse_from(args);
        assert!(opts.is_err());

        // option but no files
        let args = "fakebin -s".split_whitespace();
        let opts = ParsedCollections::try_parse_from(args);
        assert!(opts.is_err());

        // multiple files
        // option but no files
        let args = "fakebin -s Cargo.toml src/lib.rs".split_whitespace();
        let opts = ParsedCollections::try_parse_from(args);
        assert!(opts.is_ok(), "{:?}", opts);
        let opts = opts.unwrap();
        assert_eq!(&opts.paths[0], &Path::new("Cargo.toml"));
        assert_eq!(&opts.paths[1], &Path::new("src/lib.rs"));
    }
}
