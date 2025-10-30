/// A single spectral peak.
///
/// # Construction
/// This can be constructed a few different ways:
/// ```
/// use spectrum::Peak;
/// let peak1 = Peak { mz: 1.0, intensity: 2.0 };
/// let peak2 = Peak::new(1.0, 2.0);
/// let peak3 = Peak::from((1.0, 2.0));
/// assert_eq!(peak1, peak2);
/// assert_eq!(peak2, peak3);
/// ```
///
/// # Matching to Fragments
/// In mass spectrometry a single peak corresponds to some fragment of the input molecular
/// structure. However, instead of directly capturing the mass of a fragment the spectrometer
/// captures the mass-to-charge ratio (m/z) and the relative abundance (usually measured as
/// intensity).
///
/// This means that it's very easy to accidentally be incorrect when matching peaks to fragments.
/// You need to convert the m/z of the peak into an actual mass (by multiplying by peak charge) and
/// then make sure that the fragment is in its ionized form.
///
/// For example, consider [methane](https://pubchem.ncbi.nlm.nih.gov/compound/297). It has a mass
/// of about 16 Dalton. We might be tempted to match it directly to a peak with m/z 16. This might
/// seem ok if the charge of that peak is 1, but there's another catch. Methane is neutral and
/// requires some sort of mass shift to become charged. In the simplest case it can get a charge of
/// `+1` by either losing an electron or by gaining a proton, shifting its mass by the mass of an
/// electron or the mass of a proton respectively. Therefore, to match methane to a peak we need
/// to
/// 1. Convert our peak from m/z to Daltons by multiplying by the absolute value of a charge
/// 2. Convert our methane fragment from a neutral fragment to an ion by shifting its mass
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    serde::Serialize,
    serde::Deserialize,
    rkyv::Serialize,
    rkyv::Deserialize,
    rkyv::Archive,
)]
pub struct Peak {
    /// The mass-to-charge ratio of a fragment
    pub mz: f64,
    /// The relative abundance of a fragment
    pub intensity: f64,
}

impl Peak {
    /// Creates a new peak.
    pub fn new(mz: f64, intensity: f64) -> Self {
        Self { mz, intensity }
    }

    /// Converts this peak's m/z to a mass using the given charge
    pub fn mass(&self, charge: i8) -> f64 {
        self.mz * charge.abs() as f64
    }
}

impl From<(f64, f64)> for Peak {
    fn from(tuple: (f64, f64)) -> Self {
        Self::new(tuple.0, tuple.1)
    }
}
