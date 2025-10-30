use std::{
    cmp::Ordering,
    ops::{Deref, DerefMut},
};

use serde::{Deserialize, Serialize};

use crate::{ChargeVec, Ms, Peak};

/// A general purpose tandem mass spectrum, not tied to any specific format.
///
/// See also [`crate::Peaks`] and [`crate::Ms`].
// TODO: more docs
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Spectrum {
    pub(crate) peaks: Vec<Peak>,
    pub(crate) pepmass: f64,
    pub(crate) charges: ChargeVec,
    pub(crate) retention_time: Option<f64>,

    pub(crate) multistage_meta: SpectrumMultistageMeta,
    /// The scan of this spectrum
    ///
    /// This field can be used directly instead of [`Ms::get_scan`] because this type is guaranteed
    /// to always have a valid scan.
    pub scan: usize,
    /// The denominator of the normalization applied to this spectrum.
    ///
    /// If no normalization was applied, this is [`None`]
    pub denormalize_factor: Option<f64>,
}

impl Spectrum {
    /// Construct a new spectrum.
    pub fn new(
        mut peaks: Vec<Peak>,
        pepmass: f64,
        charges: ChargeVec,
        retention_time: Option<f64>,
        multistage_meta: SpectrumMultistageMeta,
        scan: usize,
    ) -> Self {
        peaks.sort_unstable_by(|x, y| x.mz.partial_cmp(&y.mz).unwrap());
        Self {
            peaks,
            pepmass,
            charges,
            retention_time,
            multistage_meta,
            scan,
            denormalize_factor: None,
        }
    }

    /// Write this spectrum to `f` in MGF format.
    ///
    /// # Errors
    /// Propagates errors from calls to [`write!`] on `f`.
    pub fn write_mgf<W: std::io::Write>(&self, mut f: W) -> std::io::Result<()> {
        writeln!(f, "BEGIN IONS")?;

        writeln!(f, "PEPMASS={}", self.pepmass)?;
        writeln!(f, "SCANS={}", self.scan)?;
        if let Some(rt) = self.retention_time {
            writeln!(f, "RTINSECONDS={}", rt)?;
        }
        for (idx, charge) in self.charges.iter().enumerate() {
            if idx == 0 {
                write!(f, "CHARGE=")?;
            }
            match charge.cmp(&0) {
                Ordering::Less => write!(f, "{}-", charge.abs())?,
                Ordering::Greater => write!(f, "{}+", charge)?,
                Ordering::Equal => write!(f, "0")?,
            };

            if idx != self.charges.len() - 1 {
                write!(f, ",")?;
            } else {
                writeln!(f)?;
            }
        }

        for peak in &self.peaks {
            writeln!(f, "{} {}", peak.mz, peak.intensity)?;
        }
        writeln!(f, "END IONS")
    }

    /// Calculate the sum of peak intensities for this spectrum.
    pub fn summed_intensity(&self) -> f64 {
        self.peaks.iter().map(|peak| peak.intensity).sum::<f64>()
    }
}

impl Deref for Spectrum {
    type Target = Vec<Peak>;

    fn deref(&self) -> &Self::Target {
        &self.peaks
    }
}

impl DerefMut for Spectrum {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.peaks
    }
}

impl From<Spectrum> for Vec<Peak> {
    fn from(spectrum: Spectrum) -> Self {
        spectrum.peaks
    }
}

impl Ms for Spectrum {
    fn pepmass(&self) -> f64 {
        self.pepmass
    }

    fn charges(&self) -> &ChargeVec {
        &self.charges
    }

    fn retention_time(&self) -> Option<f64> {
        self.retention_time
    }

    fn ms_level(&self) -> Option<usize> {
        self.multistage_meta.ms_level
    }

    fn precursor_scan(&self) -> Option<usize> {
        self.multistage_meta.precursor_scan
    }

    fn precursor_intensity(&self) -> Option<f64> {
        self.multistage_meta.precursor_intensity
    }

    fn get_scan(&self) -> Option<usize> {
        Some(self.scan)
    }
}

/// Metadata relevant to multistage spectrum parsing.
#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpectrumMultistageMeta {
    pub(crate) ms_level: Option<usize>,
    pub(crate) precursor_intensity: Option<f64>,
    pub(crate) precursor_scan: Option<usize>,
}

impl SpectrumMultistageMeta {
    /// Construct a new [`SpectrumMultistageMeta`].
    pub fn new(
        ms_level: Option<usize>,
        precursor_intensity: Option<f64>,
        precursor_scan: Option<usize>,
    ) -> Self {
        Self {
            ms_level,
            precursor_intensity,
            precursor_scan,
        }
    }
}
