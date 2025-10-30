use std::{
    borrow::{Borrow, Cow},
    path::Path,
};

use csv::WriterBuilder;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use spectrum::{SpectrumCollection, SpectrumPreprocessing};

/// Shared outputs for all molecule-spectrum matches to a particular molecule from a particular
/// spectrum file.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct MsmStub<'a> {
    /// Fully canonicalized path to the input spectrum file.
    pub spectrum_file: &'a str,
    /// Index of the matched molecule in the chemical database.
    pub mol_idx: usize,
    /// Name of the matched molecule from the chemical database.
    pub mol_name: &'a str,
    /// SMILES of matched molecule.
    pub smiles: &'a str,
    /// 14-character InChIKey prefix for the matched molecule.
    pub inchikey14: &'a str,
    /// Mass of the matched molecule
    pub mol_mass: f64,
}

/// Keeps track of the score between a particular candidate molecule and query spectrum.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Msm<'a> {
    /// Fully canonicalized path to the input spectrum file.
    pub spectrum_file: Cow<'a, str>,
    /// Scan of the matched spectrum.
    pub scan: usize,
    /// Retention time of matched spectrum, if present in original data.
    pub retention_time: Option<f64>,
    /// The mass to charge ratio of the precursor ion
    pub precursor_mz: f64,
    /// Index of the matched molecule in the chemical database.
    pub mol_idx: usize,
    /// Name of the matched molecule from the chemical database.
    pub mol_name: Cow<'a, str>,
    /// SMILES of matched molecule.
    ///
    /// This is not included in serialization or deserialization for this struct. If you want to
    /// include output SMILES use [`MsmSmiles`].
    #[serde(skip)]
    pub smiles: Option<Cow<'a, str>>,
    /// 14-character InChIKey prefix for the matched molecule.
    pub inchikey14: Cow<'a, str>,
    /// Mass of the output molecule
    pub mol_mass: f64,
    /// Mass difference between the precursor ion and the molecule.
    pub mass_error: f64,
    /// Precursor ion adduct used for this match.
    pub adduct: String,
    /// Assumed charge of precursor ion for this match.
    pub charge: i8,
    /// Score of this match
    pub score: f64,
}

/// Keeps track of the score between a particular candidate molecule and query
/// spectrum, contains SMILES of the molecule.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct MsmSmiles<'a> {
    /// Fully canonicalized path to the input spectrum file.
    pub spectrum_file: Cow<'a, str>,
    /// Scan of the matched spectrum.
    pub scan: usize,
    /// Retention time of matched spectrum, if present in original data.
    pub retention_time: Option<f64>,
    /// The mass to charge ratio of the precursor ion
    pub precursor_mz: f64,
    /// Index of the matched molecule in the chemical database.
    pub mol_idx: usize,
    /// Name of the matched molecule from the chemical database.
    pub mol_name: Cow<'a, str>,
    /// SMILES of matched molecule.
    pub smiles: Cow<'a, str>,
    /// 14-character InChIKey prefix for the matched molecule.
    pub inchikey14: Cow<'a, str>,
    /// Mass of the output molecule
    pub mol_mass: f64,
    /// Mass difference between the precursor ion and the molecule.
    pub mass_error: f64,
    /// Precursor ion adduct used for this match.
    pub adduct: String,
    /// Assumed charge of precursor ion for this match.
    pub charge: i8,
    /// Score of this match
    pub score: f64,
}

/// A [`Msm`] augmented with a p-value.
#[derive(Debug, Clone, Serialize)]
pub struct MsmPValue<'a> {
    /// Fully canonicalized path to the input spectrum file.
    pub spectrum_file: &'a str,
    /// Scan of the matched spectrum.
    pub scan: usize,
    /// Retention time of matched spectrum, if present in original data.
    pub retention_time: Option<f64>,
    /// The mass to charge ratio of the precursor ion
    pub precursor_mz: f64,
    /// Index of the matched molecule in the chemical database.
    pub mol_idx: usize,
    /// Name of the matched molecule from the chemical database.
    pub mol_name: &'a str,
    /// 14-character InChIKey prefix for the matched molecule.
    pub inchikey14: &'a str,
    /// Mass of the output molecule
    pub mol_mass: f64,
    /// Mass difference between the precursor ion and the molecule.
    pub mass_error: f64,
    /// Precursor ion adduct used for this match.
    pub adduct: &'a str,
    /// Assumed charge of precursor ion for this match.
    pub charge: i8,
    /// Score of this match
    pub score: f64,
    /// The p-value of this match.
    ///
    /// If the p-value operation failed this will contain a [`f64::NAN`].
    pub pvalue: f64,
}

impl<'a> Msm<'a> {
    /// Returns the tabular header for this struct using the provided delimiter.
    ///
    /// ```
    /// # use dereplicate::Msm;
    /// assert_eq!(
    ///     Msm::headers(b','),
    ///     "spectrum_file,scan,retention_time,precursor_mz,mol_idx,mol_name,inchikey14,mol_mass,mass_error,adduct,charge,score"
    /// );
    /// ```
    ///
    /// # Panics
    /// If somehow the writer fails to write to an in-memory [`Vec`] buffer, writes invalid UTF-8,
    /// or contains no newlines then this will panic. These are expected not to be possible and
    /// probably indicate a bug in the [`csv`] crate.
    pub fn headers(delimiter: u8) -> String {
        let mut headers = Vec::<u8>::new();
        {
            let mut writer = WriterBuilder::new()
                .delimiter(delimiter)
                .has_headers(true)
                .from_writer(&mut headers);
            writer
                .serialize(Self::default())
                .expect("failed to write example struct");
        }
        let mut headers = String::from_utf8(headers).expect("wrote invalid utf8");
        headers.truncate(headers.find('\n').expect("output had no header"));
        headers
    }

    /// Do a deep copy of this match, returning a fully owned version of it.
    pub fn to_static(&self) -> Msm<'static> {
        Msm {
            spectrum_file: self.spectrum_file.to_string().into(),
            scan: self.scan,
            retention_time: self.retention_time,
            precursor_mz: self.precursor_mz,
            mol_idx: self.mol_idx,
            mol_name: self.mol_name.to_string().into(),
            smiles: self.smiles.as_ref().map(|smiles| smiles.to_string().into()),
            inchikey14: self.inchikey14.to_string().into(),
            mol_mass: self.mol_mass,
            mass_error: self.mass_error,
            adduct: self.adduct.clone(),
            charge: self.charge,
            score: self.score,
        }
    }
}

impl<'a> MsmSmiles<'a> {
    /// Returns the tabular header for this struct using the provided delimiter.
    ///
    /// ```
    /// # use dereplicate::MsmSmiles;
    /// assert_eq!(
    ///     MsmSmiles::headers(b','),
    ///     "spectrum_file,scan,retention_time,precursor_mz,mol_idx,mol_name,smiles,inchikey14,mol_mass,mass_error,adduct,charge,score"
    /// );
    /// ```
    ///
    /// # Panics
    /// If somehow the writer fails to write to an in-memory [`Vec`] buffer, writes invalid UTF-8,
    /// or contains no newlines then this will panic. These are expected not to be possible and
    /// probably indicate a bug in the [`csv`] crate.
    pub fn headers(delimiter: u8) -> String {
        let mut headers = Vec::<u8>::new();
        {
            let mut writer = WriterBuilder::new()
                .delimiter(delimiter)
                .has_headers(true)
                .from_writer(&mut headers);
            writer
                .serialize(Self::default())
                .expect("failed to write example struct");
        }
        let mut headers = String::from_utf8(headers).expect("wrote invalid utf8");
        headers.truncate(headers.find('\n').expect("output had no header"));
        headers
    }

    /// Construct a new [`MsmSmiles`] from a [`Msm`] and a `smiles`.
    pub fn new<'b: 'a>(msm: &'b Msm<'a>, smiles: Cow<'a, str>) -> Self {
        Self {
            spectrum_file: Cow::Borrowed(msm.spectrum_file.borrow()),
            scan: msm.scan,
            retention_time: msm.retention_time,
            precursor_mz: msm.precursor_mz,
            mol_idx: msm.mol_idx,
            mol_name: Cow::Borrowed(msm.mol_name.borrow()),
            inchikey14: Cow::Borrowed(msm.inchikey14.borrow()),
            mol_mass: msm.mol_mass,
            mass_error: msm.mass_error,
            adduct: msm.adduct.clone(),
            charge: msm.charge,
            score: msm.score,
            smiles,
        }
    }

    /// Do a deep copy of this match, returning a fully owned version of it.
    pub fn to_static(&self) -> MsmSmiles<'static> {
        MsmSmiles {
            spectrum_file: self.spectrum_file.to_string().into(),
            scan: self.scan,
            retention_time: self.retention_time,
            precursor_mz: self.precursor_mz,
            mol_idx: self.mol_idx,
            mol_name: self.mol_name.to_string().into(),
            smiles: self.smiles.to_string().into(),
            inchikey14: self.inchikey14.to_string().into(),
            mol_mass: self.mol_mass,
            mass_error: self.mass_error,
            adduct: self.adduct.clone(),
            charge: self.charge,
            score: self.score,
        }
    }
}

impl<'a> From<MsmSmiles<'a>> for Msm<'a> {
    fn from(hit: MsmSmiles<'a>) -> Self {
        Self {
            spectrum_file: hit.spectrum_file,
            scan: hit.scan,
            retention_time: hit.retention_time,
            precursor_mz: hit.precursor_mz,
            mol_idx: hit.mol_idx,
            mol_name: hit.mol_name,
            smiles: Some(hit.smiles),
            inchikey14: hit.inchikey14,
            mol_mass: hit.mol_mass,
            mass_error: hit.mass_error,
            adduct: hit.adduct,
            charge: hit.charge,
            score: hit.score,
        }
    }
}

impl<'a> TryFrom<Msm<'a>> for MsmSmiles<'a> {
    type Error = crate::Error;

    fn try_from(hit: Msm<'a>) -> Result<Self, Self::Error> {
        if hit.smiles.is_none() {
            return Err(crate::Error::MissingSmiles);
        }

        Ok(Self {
            spectrum_file: hit.spectrum_file,
            scan: hit.scan,
            retention_time: hit.retention_time,
            precursor_mz: hit.precursor_mz,
            mol_idx: hit.mol_idx,
            mol_name: hit.mol_name,
            smiles: hit.smiles.unwrap(),
            inchikey14: hit.inchikey14,
            mol_mass: hit.mol_mass,
            mass_error: hit.mass_error,
            adduct: hit.adduct,
            charge: hit.charge,
            score: hit.score,
        })
    }
}

impl<'a> MsmPValue<'a> {
    /// Returns the tabular header for this struct using the provided delimiter.
    ///
    /// This is identical to [`Msm::headers`] except that it includes the extra `pvalue` field at
    /// the end.
    ///
    /// ```
    /// # use dereplicate::{Msm, MsmPValue};
    /// assert_eq!(format!("{},pvalue", Msm::headers(b',')), MsmPValue::headers(b','));
    /// ```
    ///
    /// # Panics
    /// If somehow the writer fails to write to an in-memory [`Vec`] buffer, writes invalid UTF-8,
    /// or contains no newlines then this will panic. These are expected not to be possible and
    /// probably indicate a bug in the [`csv`] crate.
    pub fn headers(delimiter: u8) -> String {
        let mut headers = Msm::headers(delimiter);
        headers.push(char::from(delimiter));
        headers.push_str("pvalue");
        headers
    }

    /// Construct a new [`MsmPValue`] from a [`Msm`] and a `pvalue`.
    pub fn new<'b: 'a>(msm: &'b Msm<'a>, pvalue: Option<f64>) -> Self {
        Self {
            spectrum_file: msm.spectrum_file.borrow(),
            scan: msm.scan,
            retention_time: msm.retention_time,
            precursor_mz: msm.precursor_mz,
            mol_idx: msm.mol_idx,
            mol_name: msm.mol_name.borrow(),
            inchikey14: msm.inchikey14.borrow(),
            mol_mass: msm.mol_mass,
            mass_error: msm.mass_error,
            adduct: &msm.adduct,
            charge: msm.charge,
            score: msm.score,
            pvalue: pvalue.unwrap_or(f64::NAN),
        }
    }
}

impl<'a> From<MsmPValue<'a>> for Msm<'a> {
    fn from(hit: MsmPValue<'a>) -> Self {
        Self {
            spectrum_file: Cow::Borrowed(hit.spectrum_file),
            scan: hit.scan,
            retention_time: hit.retention_time,
            precursor_mz: hit.precursor_mz,
            mol_idx: hit.mol_idx,
            mol_name: Cow::Borrowed(hit.mol_name),
            smiles: None,
            inchikey14: Cow::Borrowed(hit.inchikey14),
            mol_mass: hit.mol_mass,
            mass_error: hit.mass_error,
            adduct: hit.adduct.to_string(),
            charge: hit.charge,
            score: hit.score,
        }
    }
}

/// The result of running a dereplication method.
///
/// This is a wrapper around a vector of owned [`Msm`]s with some utilities for common analysis
/// added as methods.
pub struct DereplicationRun {
    matches: Vec<Msm<'static>>,
}

impl DereplicationRun {
    /// Build a [`DereplicationRun`] from a path to a TSV output from a dereplication method.
    ///
    /// All hits that are below the optional `min_score` will be filtered out. If `min_score` is
    /// [`None`] then no hits will be filtered.
    ///
    /// # Errors
    /// If this path doesn't exist this will return an error. See also errors from
    /// [`csv::Reader::from_path`].
    pub fn from_path(p: impl AsRef<Path>, min_score: Option<f64>) -> Result<Self, csv::Error> {
        let reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(p)?;
        let min_score = min_score.unwrap_or(f64::NEG_INFINITY);
        let matches = reader
            .into_deserialize()
            .filter(|m: &Result<Msm<'static>, csv::Error>| {
                // filter out low scores
                // keep errors so they bubble up on collect
                if let Ok(ref m) = m {
                    m.score >= min_score
                } else {
                    true
                }
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self { matches })
    }

    /// Construct a [`DereplicationRun`] consisting of a single `hit`.
    pub fn singleton(hit: Msm<'static>) -> Self {
        Self { matches: vec![hit] }
    }

    /// Group hits from this run by the spectrum scan they come from, keeping only the `top_k` hits
    /// per spectrum.
    ///
    /// If `top_k` is not provided then all hits for each spectrum will be kept. The return value
    /// is keyed by a [`Msm::spectrum_file`] and [`Msm::scan`] tuple. Its values are a list of
    /// references to the `top_k` [`Msm`]s that came from that spectrum.
    pub fn group_by_spectrum<'a>(
        &'a self,
        top_k: Option<usize>,
    ) -> FxHashMap<(&'a str, usize), Vec<&'a Msm<'static>>> {
        let top_k = top_k.unwrap_or(usize::MAX);
        let mut grouped = FxHashMap::default();

        for row in &self.matches {
            grouped
                .entry((row.spectrum_file.borrow(), row.scan))
                .or_insert_with(Vec::new)
                .push(row);
        }

        for value in grouped.values_mut() {
            value.sort_unstable_by(|m1, m2| m1.score.partial_cmp(&m2.score).unwrap().reverse());
            value.truncate(top_k);
        }

        grouped
    }

    /// Get a cache of parsed versions of all spectral collections in this run.
    ///
    /// The return value is keyed by the string version of a spectrum file's path on disk, as taken
    /// from the internal [`Msm`]s. The values are the parsed and pre-processed
    /// [`SpectrumCollection`]s.
    ///
    /// # Errors
    /// If [`SpectrumCollection::from_path`] fails this will propagate the error.
    pub fn parsed_spectra(
        &self,
        pre: &SpectrumPreprocessing,
    ) -> Result<FxHashMap<String, SpectrumCollection>, crate::Error> {
        let mut output = FxHashMap::default();
        for path in self.matches.iter().map(|m| &m.spectrum_file) {
            let path: &str = path.borrow();
            let entry = output.get_mut(path);

            if entry.is_none() {
                let mut collection = SpectrumCollection::from_path(path)?;
                collection.preprocess(pre);
                output.insert(path.to_string(), collection);
            }
        }
        Ok(output)
    }

    /// Get a reference to the dereplication run's matches.
    pub fn matches(&self) -> &[Msm<'static>] {
        self.matches.as_ref()
    }

    /// Get a mutable reference to the dereplication run's matches.
    pub fn matches_mut(&mut self) -> &mut [Msm<'static>] {
        self.matches.as_mut()
    }
}

impl FromIterator<Msm<'static>> for DereplicationRun {
    fn from_iter<T: IntoIterator<Item = Msm<'static>>>(iter: T) -> Self {
        Self {
            matches: iter.into_iter().collect(),
        }
    }
}

impl IntoIterator for DereplicationRun {
    type Item = Msm<'static>;

    type IntoIter = <Vec<Msm<'static>> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.matches.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::{DereplicationRun, Msm};

    #[test]
    fn headers() {
        assert_eq!(
            Msm::headers(b'\t'),
            String::from("spectrum_file\tscan\tretention_time\tprecursor_mz\tmol_idx\tmol_name\tinchikey14\tmol_mass\tmass_error\tadduct\tcharge\tscore")
        );
        assert_eq!(
            Msm::headers(b','),
            String::from("spectrum_file,scan,retention_time,precursor_mz,mol_idx,mol_name,inchikey14,mol_mass,mass_error,adduct,charge,score")
        );
    }

    #[test]
    fn group_by_spectrum() {
        let run = DereplicationRun::from_path("test_files/matches.tsv", None).unwrap();

        let grouped = run.group_by_spectrum(None);
        assert_eq!(grouped.len(), 1);
        // should be sorted by descending score
        assert_eq!(
            grouped
                .values()
                .next()
                .unwrap()
                .into_iter()
                .map(|msm| msm.score)
                .collect::<Vec<_>>(),
            vec![
                22.0, 21.0, 12.0, 10.0, 10.0, 8.0, 6.0, 5.0, 5.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, 1.0
            ],
        );

        let grouped = run.group_by_spectrum(Some(4));
        assert_eq!(grouped.len(), 1);
        // should be sorted by descending score, only with 4 left
        assert_eq!(
            grouped
                .values()
                .next()
                .unwrap()
                .into_iter()
                .map(|msm| msm.score)
                .collect::<Vec<_>>(),
            vec![22.0, 21.0, 12.0, 10.0],
        );

        // now check if min score is working
        let run = DereplicationRun::from_path("test_files/matches.tsv", Some(5.0)).unwrap();

        let grouped = run.group_by_spectrum(None);
        assert_eq!(grouped.len(), 1);
        assert_eq!(
            grouped
                .values()
                .next()
                .unwrap()
                .into_iter()
                .map(|msm| msm.score)
                .collect::<Vec<_>>(),
            vec![22.0, 21.0, 12.0, 10.0, 10.0, 8.0, 6.0, 5.0, 5.0],
        );

        let grouped = run.group_by_spectrum(Some(6));
        assert_eq!(grouped.len(), 1);
        assert_eq!(
            grouped
                .values()
                .next()
                .unwrap()
                .into_iter()
                .map(|msm| msm.score)
                .collect::<Vec<_>>(),
            vec![22.0, 21.0, 12.0, 10.0, 10.0, 8.0],
        );
    }
}
