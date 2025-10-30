use std::{
    collections::{BTreeSet, VecDeque},
    fs::File,
    io::{BufReader, Read},
    ops::{Deref, DerefMut},
    path::Path,
    str::FromStr,
};

use balloon::{
    fuzzy::{FuzzyMap, Tol},
    of64,
};
use petgraph::graph::DiGraph;

use itertools::Itertools;
use molecule::Adduct;
use ordered_float::OrderedFloat;

use crate::{errors::Error, spectrum_graph::SpectrumGraph, ChargeVec, Peak};

/// Anything that can be a peak vector.
///
/// This requires that your type can be converted into a vector of peaks by all of these methods
/// (each corresponds to a trait implementation):
/// 1. [An immutable reference](std::ops::Deref)
/// 2. [A mutable reference](std::ops::DerefMut)
/// 3. [A direct cast](std::convert::From)
///
/// If all of these are satisfied this trait provides the ability to
/// 1. [Access peaks immutably](Peaks::peaks)
/// 2. [Access peaks mutably](Peaks::peaks_mut)
/// 3. [Filtering peaks using m/z windows](Peaks::filter_peaks)
/// 4. [Merging close peaks](Peaks::merge_peaks)
/// 5. [Normalizing peaks](Peaks::normalize_peaks)
///
/// # Implementation
/// This specifically requires a [`Vec`] of [`Peak`] because the [`Peaks::filter_peaks`] and
/// [`Peaks::merge_peaks`] methods need to be able to mutate some collection.
/// Most types that implement will actually have an underlying member peak vector.
///
/// # Example
/// This trait is blanket implemented on all types that satisfy its super traits.
/// You should not implement it yourself.
///
/// ```
/// use std::ops::{Deref, DerefMut};
/// use spectrum::{Peak, Peaks};
///
/// #[derive(Clone)]
/// struct MySpectrum {
///     peaks: Vec<Peak>
/// }
///
/// impl Deref for MySpectrum {
///     type Target = Vec<Peak>;
///
///     fn deref(&self) -> &Self::Target {
///         &self.peaks
///     }
/// }
///
/// impl DerefMut for MySpectrum {
///     fn deref_mut(&mut self) -> &mut Self::Target {
///         &mut self.peaks
///     }
/// }
///
/// // you should prefer to implement From instead of Into
/// impl From<MySpectrum> for Vec<Peak> {
///     fn from(spectrum: MySpectrum) -> Self {
///         spectrum.peaks
///     }
/// }
///
/// // note we don't actually implement Peaks, we get it automatically
///
/// let my_spectrum = MySpectrum { peaks: vec![Peak::new(1.0, 2.0), Peak::new(3.0, 4.0)] };
///
/// // we can access the peak vector using the trait
/// assert_eq!(my_spectrum.peaks(), &my_spectrum.peaks);
///
/// // we can filter to keep only the most intense peak in every 10 m/z window
/// let mut filtered_spectrum = my_spectrum.clone();
/// filtered_spectrum.filter_peaks(10.0, 1);
/// assert_eq!(filtered_spectrum.peaks(), &vec![Peak::new(3.0, 4.0)]);
///
/// // we can merge together peaks within 5.0 m/z, making sure we have at least 1 final peak
/// let mut merged_spectrum = my_spectrum.clone();
/// merged_spectrum.merge_peaks(5.0, 1);
/// assert_eq!(merged_spectrum.peaks(), &vec![Peak::new(3.0, 6.0)]);
/// ```
pub trait Peaks: Deref<Target = Vec<Peak>> + DerefMut + Into<Vec<Peak>> {
    /// Get an immutable reference to the underlying peak vector.
    fn peaks(&self) -> &Vec<Peak> {
        self.deref()
    }

    /// Get a mutable reference to the underlying peak vector.
    fn peaks_mut(&mut self) -> &mut Vec<Peak> {
        self.deref_mut()
    }

    // TODO: update docs
    /// Create a datastructure for fast peak m/z existence queries.
    fn peak_set(&self, adducts: &[Adduct]) -> BTreeSet<OrderedFloat<f64>> {
        self.peaks()
            .iter()
            .flat_map(|peak| {
                adducts
                    .iter()
                    .map(move |adduct| adduct.to_mass(peak.mz).into())
            })
            .collect()
    }

    /// Filters out peaks in this spectrum.
    ///
    /// This creates non-overlapping m/z windows of `window_size`.
    /// Within each window on the `num_peaks` most intense peaks are kept.
    /// The rest of the peaks are removed.
    /// Note the windows are in units of m/z, the number of peaks in each window is highly
    /// dependent on the values of [`Peak::mz`].
    ///
    /// This function will also sort the peaks by increasing intensity.
    fn filter_peaks(&mut self, window_size: f64, num_peaks: usize) {
        // if only 1 peak just leave
        if self.len() <= 1 {
            return;
        }

        // sort peaks by m/z, unwrap is safe if no NaNs
        self.peaks_mut()
            .sort_by(|peak1, peak2| peak1.mz.partial_cmp(&peak2.mz).unwrap());

        // we know we have at least 2 peaks and that the peaks are sorted
        let min_mz = self[0].mz;
        let num_windows = ((self.last().unwrap().mz - min_mz) / window_size) as usize + 1;

        let mut new_peaks: Vec<Peak> = Vec::with_capacity(num_windows * num_peaks);
        for window in 1..=num_windows {
            if self.is_empty() {
                break;
            }
            let max_mz = min_mz + ((window as f64) * window_size);

            // find first peak with m/z out of accepted range
            let end_idx = self
                .peaks()
                .iter()
                .position(|x| x.mz > max_mz)
                .unwrap_or_else(|| self.len());

            // no peaks in this window, just continue on
            if end_idx == 0 {
                continue;
            }

            // sort current window by intensity, unwrap is safe if no NaNs
            // add this to new peak vector
            new_peaks.extend(
                self.peaks_mut()
                    .drain(0..end_idx)
                    .sorted_by(|peak1, peak2| {
                        peak2.intensity.partial_cmp(&peak1.intensity).unwrap()
                    })
                    .take(num_peaks),
            );
        }

        new_peaks.shrink_to_fit();
        *self.peaks_mut() = new_peaks;
    }

    /// Merges close-together peaks together in this spectrum.
    ///
    /// Combines peaks within `merge_thresh` m/z of each other into a single peak. When two peaks are
    /// merged the one with lower intensity is merged into the peak with the larger intensity. If
    /// after merging there are fewer than `min_peaks` peaks then the merge operation becomes a
    /// no-op.
    ///
    /// This function will also sort the peaks by increasing intensity.
    fn merge_peaks(&mut self, merge_thresh: f64, min_peaks: usize) {
        self.peaks_mut()
            .sort_by(|peak1, peak2| peak1.mz.partial_cmp(&peak2.mz).unwrap());
        let mut merged_peaks: Vec<Peak> = self.clone();

        // Keep track of the windows of valid peaks
        let mut queue: VecDeque<(usize, Peak)> = VecDeque::with_capacity(self.len());

        for (i, peak) in self.iter().enumerate() {
            // Pop off peaks that are too far behind
            while !queue.is_empty() && queue.front().unwrap().1.mz + merge_thresh < peak.mz {
                queue.pop_front();
            }

            // Push the current peak if need be
            if queue.is_empty() || (queue.back().unwrap().1.mz - peak.mz).abs() > f64::EPSILON {
                queue.push_back((i, *peak));
            }

            // Push on peaks that have higher m_z
            let mut highest_ind: usize = queue.back().unwrap().0 + 1;
            while highest_ind < self.len() && peak.mz + merge_thresh >= self[highest_ind].mz {
                queue.push_back((highest_ind, self[highest_ind]));
                highest_ind += 1;
            }

            let max_ind: usize = queue
                .iter()
                .max_by(|x, y| x.1.intensity.partial_cmp(&y.1.intensity).unwrap())
                .unwrap()
                .0;
            if max_ind != i {
                merged_peaks[max_ind].intensity += peak.intensity;
                merged_peaks[i].mz *= -1.0; // used later for filtering out peaks that were not largest in their range
            }
        }

        let new_peaks = merged_peaks
            .into_iter()
            .filter(|x| x.mz > 0.0)
            .collect::<Vec<Peak>>();

        // Only update if number of peaks remaining is sufficiently large
        if new_peaks.len() >= min_peaks {
            *self.peaks_mut() = new_peaks;
        }
    }

    /// Take the square root of peak intensity vector and then normalize it. This returns the
    /// denominator which is required to de-normalize peak intensities via `[denormalize_peaks]`.
    ///
    /// Normalization enables cosine similiarity calculation. If all peak intensities are zero,
    /// the intensity vector cannot be normalized. In this case we don't panic but instead return the same zeroed intensities.
    fn normalize_peaks(&mut self) -> Option<f64> {
        let mut magnitude = 0.0;
        for peak in self.peaks_mut() {
            magnitude += peak.intensity;
        }
        let denominator = f64::powf(magnitude, 0.5);
        if denominator > 0.0 {
            for peak in self.peaks_mut() {
                peak.intensity = f64::powf(peak.intensity, 0.5) / denominator;
            }
        }
        // it's okay if 0.0 is returned because we'd just multiply each intensity 0 by 0 if we denormalize
        Some(denominator)
    }

    /// Given the denominator returned from `[normalize_peaks]`, undo the normalization applied to each peak to
    /// restore peak intensities to their original intensities.
    fn denormalize_peaks(&mut self, denormalize_factor: f64) {
        for peak in self.peaks_mut() {
            peak.intensity = f64::powf(peak.intensity * denormalize_factor, 2.0);
        }
    }

    /// Filters out peaks in this spectrum according to MASST+ requirement.
    ///
    /// This creates non-overlapping m/z windows of `window_size` based on the implementation of Tyler's MASST+
    /// Within each window on the `num_peaks` most intense peaks are kept.
    /// The rest of the peaks are removed.
    ///
    /// This function will also sort the peaks by increasing intensity.
    fn filter_peaks_masst(&mut self, window_size: f64, num_peaks: usize) {
        if self.is_empty() {
            return;
        }

        // sort peaks by m/z, unwrap is safe if no NaNs
        self.peaks_mut()
            .sort_by(|peak1, peak2| peak1.mz.partial_cmp(&peak2.mz).unwrap());
        let num_windows = ((self.last().unwrap().mz) / window_size) as usize + 1;
        let mut new_peaks = Vec::<Peak>::with_capacity(num_windows * num_peaks);
        let mut window_counter = Vec::<usize>::new();
        window_counter.resize(num_windows, 0);
        for peak in self.peaks() {
            let window_location = (((peak.mz as usize) as f64) / window_size) as usize;
            window_counter[window_location] += 1;
        }
        for &end_idx in &window_counter {
            if end_idx == 0 {
                continue;
            }
            new_peaks.extend(
                self.peaks_mut()
                    .drain(0..end_idx)
                    .sorted_by(|peak1, peak2| {
                        peak2.intensity.partial_cmp(&peak1.intensity).unwrap()
                    })
                    .take(num_peaks),
            );
        }
        new_peaks.shrink_to_fit();
        *self.peaks_mut() = new_peaks;
        // TODO: why even do this?
        self.peaks_mut()
            .sort_by(|peak1, peak2| peak1.mz.partial_cmp(&peak2.mz).unwrap());
    }

    /// Keeps only the `num_peaks` highest-peak-intensity peaks in the spectrum, discarding the
    /// rest.
    fn filter_top_n_spectrum_peaks(&mut self, num_peaks: usize) {
        // sort such that highest intensity peaks are first
        self.peaks_mut()
            .sort_by(|peak1, peak2| peak2.intensity.partial_cmp(&peak1.intensity).unwrap());
        self.peaks_mut().truncate(num_peaks);
        // TODO: why even do this?
        self.peaks_mut()
            .sort_by(|peak1, peak2| peak1.mz.partial_cmp(&peak2.mz).unwrap());
    }

    /// Converts the peaks into a spectrum graph given a list of monomer masses.
    ///
    /// A spectrum graph is a directed acyclic graph where the nodes correspond to peaks, and there
    /// exists a directed edge between two node (peaks) if their mass difference exists in the
    /// provided `monomer_masses`. The runtime of this construction is `O(n^2)`, where `n` is the
    /// number of peaks in the spectrum.
    fn spectrum_graph<'a>(
        &'a self,
        monomer_masses: &[of64],
        monomer_names: &[String],
        precursor_mass: f64,
        tol: f64,
    ) -> Result<SpectrumGraph<'a, Self>, balloon::NanError> {
        let mut graph = DiGraph::new();
        let tol = Tol::Abs(of64::try_from(tol).unwrap());

        let mass_names = monomer_masses
            .iter()
            .cloned()
            .zip(monomer_names.iter().cloned());

        let peak_to_name = FuzzyMap::from_iter(mass_names, tol);

        for peak in self.peaks() {
            graph.add_node(*peak);
        }

        // add water as a starting point, and precursor mass as an ending point
        graph.add_node(Peak::new(1.00794, 1.0));
        graph.add_node(Peak::new(1.00794 + 18.0105646834, 1.0));
        graph.add_node(Peak::new(precursor_mass, 1.0));
        graph.add_node(Peak::new(precursor_mass - 1.00794, 1.0));

        for (i, j) in graph.node_indices().tuple_combinations() {
            let i_mz = of64::try_from(graph[i].mz)?;
            let j_mz = of64::try_from(graph[j].mz)?;
            let (source, target) = if i_mz < j_mz { (i, j) } else { (j, i) };

            // no absolute value, since we want our graph to be a DAG
            let mass_diff = of64::try_from(graph[target].mz - graph[source].mz)?;

            if let Some(name) = peak_to_name.get(mass_diff) {
                graph.add_edge(source, target, name.to_string());
            }
        }

        Ok(SpectrumGraph { graph, peaks: self })
    }
}

impl<T: Deref<Target = Vec<Peak>> + DerefMut + Into<Vec<Peak>>> Peaks for T {}

/// Anything that is considered a tandem mass spectrum.
///
/// This requires that your type implements [`Peaks`], since every mass spectrum should contain a
/// peak list.
///
/// This essentially abstracts getters for the necessary pieces of a *tandem* mass spectrum. We
/// require [a precursor m/z](Ms::pepmass), [a vector of possible charges](Ms::charges), and
/// optionally [a retention time](Ms::retention_time). Using an empty charge vector is acceptable.
pub trait Ms: Peaks {
    /// The precursor m/z of this spectrum.
    ///
    /// This is named confusingly due to the MGF convention of using this for the precursor m/z.
    fn pepmass(&self) -> f64;
    /// The possible charges of this spectrum.
    ///
    /// This value can often be unreliable, so it should be ignored if possible.
    fn charges(&self) -> &ChargeVec;
    /// The retention time for this spectrum.
    ///
    /// This is often just [`None`], but it can be useful information.
    fn retention_time(&self) -> Option<f64>;
    /// The ms level of this spectrum.
    ///
    /// If the original spectrum file didn't contain an ms level then this will return [`None`].
    fn ms_level(&self) -> Option<usize>;
    /// The precursor scan for this spectrum
    ///
    /// If the original spectrum file has an ms level of 1 or otherwise does not have a
    /// precursor scan, this will return [`None`].
    fn precursor_scan(&self) -> Option<usize>;
    /// The precursor intensity for this spectrum.
    ///
    /// If the original spectrum file has an ms level of 1 or otherwise does not have a
    /// precursor intensity, this will return [`None`].
    fn precursor_intensity(&self) -> Option<f64>;
    /// The scan number for this spectrum.
    ///
    /// If the original spectrum file didn't contain a scan then this will return [`None`].
    fn get_scan(&self) -> Option<usize>;

    /// Compute the log ranks of a spectrum.
    ///
    /// The log rank of a peak is defined as the following expression:
    /// ```txt
    /// min(floor(log_2(rank)) + 1, max_rank)
    /// ```
    /// where rank is the ranking of the peak's intensity sorted in descending
    /// order from 1 to n.
    ///
    /// Index `i` of the returned vector is the log rank of peak `i` in the spectrum. Note that the
    /// ranking is no longer valid if the order of the peaks changes.
    ///
    /// Citations:
    /// [\1\]: [Cao, L., Guler, M., Tagirdzhanov, A. et al. MolDiscovery: learning mass spectrometry fragmentation of small molecules. Nat Commun 12, 3718 (2021). https://doi.org/10.1038/s41467-021-23986-0](https://doi.org/10.1038/s41467-021-23986-0)
    fn log_ranks(&self, max_rank: u8) -> Vec<u8> {
        // keep track of the original indices and then sort by intensity
        let mut peak_indices = (0..self.peaks().len()).collect::<Vec<usize>>();
        peak_indices.sort_by(|a, b| {
            self.peaks()[*b]
                .intensity
                .partial_cmp(&self.peaks()[*a].intensity)
                .unwrap()
        });

        let mut log_ranks = vec![0; peak_indices.len()];
        for (rank, original_index) in peak_indices.iter().enumerate() {
            let rank = (rank + 1) as f64;
            log_ranks[*original_index] = std::cmp::min((rank.log2().floor() + 1.0) as u8, max_rank);
        }
        log_ranks
    }
}

/// A marker trait for mass spectra formats.
///
/// Provides functionality to read the format from a file and a path given that it's able to be
/// read from a string already. The implementation of [`FromStr`] for the implementing type must
/// have an error type that can be coerced to a [`Error`].
///
/// Currently we have no streaming parsers, so any file format parser will read the whole file to
/// string and hold it in memory. This ends up being very important for performance in many cases,
/// with the trade-off being possibly elevated memory usage.
pub trait MsFormat: FromStr
where
    Error: From<<Self as FromStr>::Err>,
{
    /// Reads the spectrum collection from a path on disk.
    fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        Ok(std::fs::read_to_string(&path)?.parse()?)
    }

    /// Reads the spectrum collection from anything that implements [`Read`].
    ///
    /// Note this wraps the reader in a buffered reader itself, so caller should not buffer the
    /// reader they pass in.
    fn from_reader<R: Read>(reader: R) -> Result<Self, Error> {
        let mut reader = BufReader::new(reader);
        let mut s = String::new();
        reader.read_to_string(&mut s)?;
        Ok(s.parse()?)
    }

    /// Reads the spectrum collection from a file.
    ///
    /// This is just a convenience wrapper around [`from_reader`](Self::from_reader).
    fn from_file(file: File) -> Result<Self, Error> {
        Self::from_reader(file)
    }

    /// Reads the spectrum collection from a string, unifying error types.
    fn from_str(data: &str) -> Result<Self, Error> {
        Ok(data.parse()?)
    }
}
