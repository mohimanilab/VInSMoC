use std::str::FromStr;

use molecule::Adduct;

/// An error tolerance for comparing theoretical and experimental peak masses.
///
/// Error tolerances can either be given in Daltons, which is an absolute error tolerance, or in
/// PPM, which is a relative error tolerance.
///
/// ```
/// # use chemical_db::ErrorTol;
/// let tol = ErrorTol::from_dalton(1.0);
/// assert!(tol.are_equal(2.0, 2.5));
/// assert!(tol.are_equal(2.0, 3.0));
/// assert!(!tol.are_equal(2.0, 3.1));
///
/// // relative tolerances depend on the order of inputs
/// let tol = ErrorTol::from_ppm(10000.0);
/// assert!(tol.are_equal(1000.0, 990.0));
/// assert!(tol.are_equal(990.0, 1000.0));
/// assert!(tol.are_equal(10000.0, 9900.0));
/// assert!(tol.are_equal(9900.0, 100000.0));
/// ```
///
/// # Constructing from strings
/// An error tolerance can be built from a string that consists of the floating point error
/// tolerance optionally suffixed by either "da" for absolute or "ppm" for relative. If no suffix
/// is provided then error tolerances default to absolute.
///
/// ```
/// # use chemical_db::ErrorTol;
/// assert_eq!("1.0ppm".parse::<ErrorTol>().unwrap(), ErrorTol::Ppm(1e-6));
/// assert_eq!("1.0da".parse::<ErrorTol>().unwrap(), ErrorTol::Da(1.0));
/// assert_eq!("1.0".parse::<ErrorTol>().unwrap(), ErrorTol::Da(1.0));
/// assert!("1.0other".parse::<ErrorTol>().is_err())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum ErrorTol {
    /// An absolute error tolerance given and stored in Daltons.
    Da(f64),
    /// A relative error tolerance given as a PPM and stored as `PPM / 1e6`.
    Ppm(f64),
}

impl ErrorTol {
    /// Checks if two floating point values are equal according to this error tolerance.
    pub fn are_equal(&self, expected: f64, observed: f64) -> bool {
        match self {
            ErrorTol::Da(tol) => (expected - observed).abs() <= *tol,
            ErrorTol::Ppm(tol) => (expected - observed).abs() / expected <= tol * expected,
        }
    }

    /// Create an absolute tolerance that is exactly enough to mark `expected` and `observed` as
    /// equal.
    pub fn absolute(expected: f64, observed: f64) -> Self {
        let tol = (expected - observed).abs();
        Self::Da(tol)
    }

    /// Create a relative tolerance that is exactly enough to mark `expected` and `observed` equal.
    pub fn relative(expected: f64, observed: f64) -> Self {
        let tol = (expected - observed).abs() / expected;
        Self::Ppm(tol)
    }

    /// Build a relative error tolerance from the provided `ppm` value.
    pub fn from_ppm(ppm: f64) -> Self {
        Self::Ppm(ppm / 1e6)
    }

    /// Build an exact error tolerance from the provided `dalton` tolerance.
    pub fn from_dalton(dalton: f64) -> Self {
        Self::Da(dalton)
    }

    /// Compute the lower and upper bounds of equal floating point values for this error tolerance.
    ///
    /// The returned tuple is made up of an inclusive lower and inclusive upper bound.
    pub fn bounds(&self, expected: f64) -> (f64, f64) {
        match self {
            ErrorTol::Da(tol) => (expected - tol, expected + tol),
            ErrorTol::Ppm(tol) => {
                let tol = tol * expected;
                (expected - tol, expected + tol)
            }
        }
    }

    /// Retrieve the inner absolute or relative floating point value.
    pub fn internal(&self) -> f64 {
        match *self {
            ErrorTol::Da(x) | ErrorTol::Ppm(x) => x,
        }
    }
}

impl Default for ErrorTol {
    fn default() -> Self {
        Self::from_dalton(0.01)
    }
}

impl PartialEq<f64> for ErrorTol {
    fn eq(&self, other: &f64) -> bool {
        match self {
            ErrorTol::Da(tol) | ErrorTol::Ppm(tol) => tol == other,
        }
    }
}

impl FromStr for ErrorTol {
    type Err = std::num::ParseFloatError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_ascii_lowercase();
        let mut s = s.as_str();

        // check end of string for units
        let mut abs = true;

        // .unwrap()s safe since guarded with .ends_with if checks
        if s.ends_with("ppm") {
            abs = false;
            s = s.strip_suffix("ppm").unwrap().trim();
        } else if s.ends_with("da") {
            s = s.strip_suffix("da").unwrap().trim();
        }

        // parse the float part
        // errors from random non ppm/da extension will show up here
        let tol = s.parse::<f64>()?;

        Ok(if abs {
            Self::from_dalton(tol)
        } else {
            Self::from_ppm(tol)
        })
    }
}

/// Parameters for selecting candidates from a chemical database.
///
/// # [`clap`]
/// This implements [`clap::Parser`] so it can be used to parse command line arguments.
/// ```
/// use clap::Parser;
/// use chemical_db::CandidateFilter;
///
/// // it can be used to add candidate filtering parameters to a CLI
/// #[derive(Debug, Clone, Default, Parser)]
/// #[clap(name = "bin")]
/// struct Opts {
///     #[clap(flatten)]
///     filter: CandidateFilter,
/// }
///
/// // long arguments
/// let args = "bin --precursor-ion-thresh 0.1 --precursor-adducts [M+H]+ --respect-charges".split_whitespace();
/// let params = Opts::try_parse_from(args).unwrap();
/// assert_eq!(params.filter.error_tolerance, 0.1);
/// assert_eq!(params.filter.adducts, vec!["[M+H]+".parse().unwrap()]);
/// assert!(params.filter.respect_charges);
///
/// // short arguments
/// let args = "bin -P 0.1 -A [M+H]+".split_whitespace();
/// let params = Opts::try_parse_from(args).unwrap();
/// assert_eq!(params.filter.error_tolerance, 0.1);
/// assert_eq!(params.filter.adducts, vec!["[M+H]+".parse().unwrap()]);
/// assert!(!params.filter.respect_charges);
///
/// // default CLI args match `Default` implementation
/// let args = "bin".split_whitespace();
/// let params = Opts::try_parse_from(args).unwrap();
/// assert_eq!(params.filter.error_tolerance, CandidateFilter::default().error_tolerance);
/// assert_eq!(&params.filter.adducts, &CandidateFilter::default().adducts);
/// assert_eq!(params.filter.respect_charges, CandidateFilter::default().respect_charges);
/// assert_eq!(params.filter.error_tolerance, 0.01);
/// assert_eq!(
///     params.filter.adducts,
///     vec![
///         "[M+H]+".parse().unwrap(),
///         "[M+2H]++".parse().unwrap(),
///         "[M+3H]+++".parse().unwrap(),
///     ]
/// );
/// assert!(!params.filter.respect_charges);
/// ```
#[derive(Debug, Clone, clap::Parser)]
#[clap(help_heading = "Database")]
pub struct CandidateFilter {
    /// Precursor m/z error tolerance.
    ///
    /// The error tolerance used to match spectral precursor m/zs to candidate molecule masses.
    /// Molecules in the chemical database with masses further away than this threshold are not
    /// scored against that spectrum. Threshold is applied after conversion of precursor m/zs to
    /// masses via the precursor adducts. Error tolerances by default are given as absolute
    /// tolerances in Dalton. The number can be suffixed with "ppm" to convert this tolerance into
    /// a relative tolerance computed using parts-per-million.
    #[clap(
        name = "PRECURSOR-ION-THRESH",
        short = 'P',
        long = "precursor-ion-thresh",
        default_value = "0.01"
    )]
    pub error_tolerance: ErrorTol,
    /// Adducts used for precursor ions
    #[clap(
        name = "PRECURSOR-ADDUCTS",
        short = 'A',
        long = "precursor-adducts",
        default_values = &["[M+H]+", "[M+2H]++", "[M+3H]+++"],
        multiple_values = true,
        takes_value = true,
    )]
    pub adducts: Vec<Adduct>,
    /// Whether or not to respect spectrum charges.
    ///
    /// If this is set then only precursor adducts that match the stated charge of the
    /// spectrum will be considered for candidates. If a spectrum contains no charge information
    /// (an empty charge vector) then that spectrum will receive no candidates. By default this is
    /// disabled and the charge field for the input spectra will be ignored in favor of manual
    /// charge ranges set via the charges of the allowed precursor adducts.
    #[clap(long)]
    pub respect_charges: bool,
}

impl Default for CandidateFilter {
    fn default() -> Self {
        Self {
            error_tolerance: ErrorTol::default(),
            // .unwrap() is safe here because we know for sure these are valid adducts
            adducts: vec![
                "[M+H]+".parse().unwrap(),
                "[M+2H]++".parse().unwrap(),
                "[M+3H]+++".parse().unwrap(),
            ],
            respect_charges: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults() {
        assert_eq!(
            CandidateFilter::default().error_tolerance,
            ErrorTol::Da(0.01)
        );
        assert!(!CandidateFilter::default().respect_charges);
        assert_eq!(
            CandidateFilter::default().adducts,
            vec![
                "[M+H]+".parse().unwrap(),
                "[M+2H]++".parse().unwrap(),
                "[M+3H]+++".parse().unwrap(),
            ]
        );
    }

    #[test]
    fn negative_adducts() {
        use clap::Parser;

        #[derive(Debug, Clone, Default, Parser)]
        struct Opts {
            #[clap(flatten)]
            filter: CandidateFilter,
        }

        let args = "bin -A [M-H]-".split_whitespace();
        let params = Opts::try_parse_from(args).unwrap();
        assert_eq!(params.filter.adducts, vec!["[M-H]-".parse().unwrap()]);

        let args = "bin -A [M-H]- [M-2H]--".split_whitespace();
        let params = Opts::try_parse_from(args).unwrap();
        assert_eq!(
            params.filter.adducts,
            vec!["[M-H]-".parse().unwrap(), "[M-2H]--".parse().unwrap(),]
        );

        let args = "bin -A [M-H]- [M-2H]-- [M+H]+".split_whitespace();
        let params = Opts::try_parse_from(args).unwrap();
        assert_eq!(
            params.filter.adducts,
            vec![
                "[M-H]-".parse().unwrap(),
                "[M-2H]--".parse().unwrap(),
                "[M+H]+".parse().unwrap(),
            ]
        );
    }
}
