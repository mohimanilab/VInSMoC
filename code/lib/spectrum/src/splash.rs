//! Spectral hashes using the SPLASH method.
//!
//! This module contains the [`Splash`] type, which is a spectral hash designed for mass spectra.
//! [`Splash`] contains four blocks separated by a "-".
//!
//! The first block is a version identifier, always made up by the string "splash" followed by a
//! number indicating the type of data encoded (`1` for mass spectra) and the splash version
//! (currently `0`). We only support mass spectra and version `0`, see [`BLOCK1`].
//!
//! The second block is four characters that are either `0-9` or `a-z`. This is the result of
//! 1. Picking only the ten most intense peaks that are at least 10% of the intensity of the most
//!    intense peak
//! 2. Converting into a length 10 wrapping histogram, summing intensities over bins of 5 Da
//! 3. Converting each bin into a base 3 number
//! 4. Converting the whole base 3 length 10 array into a base 36 number (but only using 30 bits)
//!
//! The third block is ten characters that are either `0-9` or `a-z`. This is the result of
//! 1. Creating a length 10 wrapping histogram, summing intensities over bins of 100 Da
//! 2. Converting each bin into a single base 10 digit
//! 3. Outputting these ten digits
//!
//! The fourth and final block is the first 20 characters of a [`sha2::Sha256`] hash of the
//! spectrum. This is what can be used to compare for uniqueness. The spectral peaks have their m/z
//! and intensities scaled and are joined into a string with a pre-defined format.
//!
//! For more information see [the original paper](https://doi.org/10.1038/nbt.3689).

use std::{fmt, str::FromStr};

use itertools::Itertools;
use sha2::Digest;

use crate::{Peak, Peaks};

/// The first '-' separated block of a SPLASH.
///
/// This is hard coded to `splash10` since we only support that version and only work on mass
/// spectra.
///
/// ```
/// assert_eq!(spectrum::splash::BLOCK1, "splash10");
/// ```
pub const BLOCK1: &str = "splash10";

/// A correction value to avoid zeros.
///
/// This was present in the [original implementations of SPLASH](https://github.com/berlinguyinca/spectra-hash/blob/master/cpp/src/splash.cpp).
/// I dislike this approach but am using it to match their versions.
const EPS_CORRECTION: f64 = 1.0e-7;

/// The scaling factor each peak's m/z ratio is multiplied by before being hashed.
const MZ_PRECISION_FACTOR: f64 = 1e6;

/// The scaling factor each peak's intensity is multiplied by before being hashed.
const INTENSITY_PRECISION_FACTOR: f64 = 1e0;

/// A lookup table for base 36 encoded integers.
const BASE_36_TABLE: [u8; 36] = [
    b'0', b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9', b'a', b'b', b'c', b'd', b'e', b'f',
    b'g', b'h', b'i', b'j', b'k', b'l', b'm', b'n', b'o', b'p', b'q', b'r', b's', b't', b'u', b'v',
    b'w', b'x', b'y', b'z',
];

/// A spectral hash value.
///
/// For details on its contents see [the module level docs](crate::splash). This can be constructed
/// using its [`From`] implementation on an implementer of [`Peaks`]. Individual blocks can be
/// accessed with getters and the full version can be accessed with its [`fmt::Display`]
/// implementation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Splash {
    block2: [u8; 4],
    block3: [u8; 10],
    block4: [u8; 20],
}

/// Build block 2 of a SPLASH from a list of peaks.
fn block2(peaks: &[Peak]) -> [u8; 4] {
    const BIN_SIZE: f64 = 5.0;
    const BASE: f64 = 3.0;

    // we work on a filtered peak set for this block
    // up to the 10 most intense peaks, but only peaks that have intensity at least 10% of the
    // most intense peak
    let max_intensity = peaks
        .iter()
        .map(|peak| peak.intensity)
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(0.0);

    let min_kept_intensity = max_intensity * 0.10;
    let mut filter_peaks = peaks
        .iter()
        .copied()
        .filter(|peak| peak.intensity + EPS_CORRECTION >= min_kept_intensity)
        .collect::<Vec<_>>();
    filter_peaks.sort_unstable_by(|x, y| x.intensity.partial_cmp(&y.intensity).unwrap().reverse());
    filter_peaks.truncate(10);

    // we want only the top 10 peaks that are greater than 10% of the most intense peak

    let hist = histogram(&filter_peaks, BIN_SIZE, BASE);

    // we have 10 base 3 numbers (30 bits)
    // we want a single base 36 number

    // first unwrap safe because hist can only contain 0-9,a-z
    // second unwrap safe because the input base to Splash::histogram was 3
    let mut num = usize::from_str_radix(std::str::from_utf8(&hist).unwrap(), BASE as u32).unwrap();
    let mut block = [b'0'; 4];

    // going in reverse automatically gives us padding in the most significant positions
    let mut idx = 3usize;
    while num > 0 {
        block[idx] = BASE_36_TABLE[num % 36];
        num /= 36;
        idx -= 1;
    }

    block
}

/// Build block 3 of a SPLASH from a list of peaks
fn block3(peaks: &[Peak]) -> [u8; 10] {
    const BIN_SIZE: f64 = 100.0;
    const BASE: f64 = 10.0;
    histogram(peaks, BIN_SIZE, BASE)
}

/// Build block 4 of a SPLASH from a list of peaks
fn block4(peaks: &[Peak]) -> [u8; 20] {
    // the spec requires all intensities to be normalized
    // just like the reference implementation we will panic if the maximum intensity peak is 0.0
    let max_intensity = peaks
        .iter()
        .map(|peak| peak.intensity)
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(0.0);

    let mut encoded = String::new();

    // encode with
    // "mz:intensity mz:intensity"
    // normalizing intensities along the way
    for (idx, peak) in peaks
        .iter()
        .sorted_by(|x, y| {
            x.mz.partial_cmp(&y.mz)
                .unwrap()
                .then(x.intensity.partial_cmp(&y.intensity).unwrap().reverse())
        })
        .enumerate()
    {
        let mz = ((peak.mz + EPS_CORRECTION) * MZ_PRECISION_FACTOR) as usize;
        // normalize intensity since this is required by spec
        // intensity becomes % of max intensity peak
        let intensity = peak.intensity / max_intensity * 100.0 + EPS_CORRECTION;
        let intensity = (intensity * INTENSITY_PRECISION_FACTOR) as usize;

        encoded.push_str(&format!("{}:{}", mz, intensity));
        if idx != peaks.len() - 1 {
            encoded.push(' ');
        }
    }

    // hash the encoded string
    // the base16ct::lower::encode_string matches python's hex_digest
    // .unwrap() can't panic since we take the 20 characters
    let mut hasher = sha2::Sha256::new();
    hasher.update(&encoded);
    base16ct::lower::encode_string(&hasher.finalize())[..20]
        .as_bytes()
        .try_into()
        .unwrap()
}

/// Build a histogram of intensities over binned m/zs.
///
/// The final output is always length 10 and overflowed bins wrap to fill from the beginning. The
/// provided `bin_size` is the length of each bin in Daltons. The `base` defines what single digit
/// base each bin should be converted to before finally returning.
fn histogram(peaks: &[Peak], bin_size: f64, base: f64) -> [u8; 10] {
    let mut int_hist = [0.0; 10];
    for &Peak { mz, intensity } in peaks {
        let idx = (mz / bin_size) as usize % 10;
        int_hist[idx] += intensity;
    }

    // .unwrap() safe because int_hist has len 10
    let mut hist = [0u8; 10];
    let max_bin = int_hist
        .iter()
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap();

    for (idx, bin) in int_hist.into_iter().enumerate() {
        hist[idx] = (EPS_CORRECTION + (base - 1.0) * bin / max_bin) as u8;
    }

    for bin in &mut hist {
        *bin = BASE_36_TABLE[*bin as usize];
    }
    hist
}

impl Splash {
    /// Get the SPLASH's first block.
    pub fn block1(&self) -> &'static str {
        BLOCK1
    }

    /// Get the SPLASH's second block.
    pub fn block2(&self) -> &[u8; 4] {
        &self.block2
    }

    /// Get the SPLASH's third block.
    pub fn block3(&self) -> &[u8; 10] {
        &self.block3
    }

    /// Get the SPLASH's fourth block.
    pub fn block4(&self) -> &[u8; 20] {
        &self.block4
    }
}

impl<S: Peaks> From<&S> for Splash {
    fn from(peaks: &S) -> Self {
        Self {
            block2: block2(peaks),
            block3: block3(peaks),
            block4: block4(peaks),
        }
    }
}

impl fmt::Display for Splash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // .unwraps safe here because they're either entries in BASE_36_TABLE or an output from SHA
        // hashing (plus hex encoding), all of which are valid utf8
        let block2 = std::str::from_utf8(&self.block2).unwrap();
        let block3 = std::str::from_utf8(&self.block3).unwrap();
        let block4 = std::str::from_utf8(&self.block4).unwrap();
        write!(f, "{}-{}-{}-{}", BLOCK1, block2, block3, block4)
    }
}

/// An error that happens when parsing a SPLASH from a string.
#[derive(Debug, thiserror::Error)]
pub enum InvalidSplash {
    /// A splash requires four "-" separated blocks, this is raised when there aren't four and it
    /// stores the number of blocks actually found.
    #[error("found {0} - separated blocks, expected 4")]
    BlockLength(usize),
    /// The first block is not the same as our [`BLOCK1`].
    #[error("first block not equal to \"splash10\"")]
    Block1,
    /// A block had the incorrect length.
    ///
    /// This is only enforced for blocks 2-4.
    #[error("block {num} expected {expected} characters but found {found}")]
    BlockLen {
        /// The 1-based block number
        num: u8,
        /// Expected block length
        expected: usize,
        /// Found block length
        found: usize,
    },
    /// A block had invalid characters.
    ///
    /// Blocks 2 and 3 can only have characters that show up in a base 36 number (`0-9` or `a-z`).
    /// If any other characters are found this error is raised. This is not enforced for block 4.
    #[error("block {num} found unexpected character {c}")]
    BlockChars {
        /// The 1-based block number
        num: u8,
        /// The invalid character we found
        c: char,
    },
}

fn valid_byte(byte: u8) -> bool {
    (b'0'..=b'9').contains(&byte) || (b'a'..=b'z').contains(&byte)
}

impl FromStr for Splash {
    type Err = InvalidSplash;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let blocks = s.split('-').map(str::as_bytes).collect::<Vec<_>>();
        if blocks.len() != 4 {
            return Err(InvalidSplash::BlockLength(blocks.len()));
        }

        if blocks[0] != BLOCK1.as_bytes() {
            return Err(InvalidSplash::Block1);
        }

        // parse block 2
        if blocks[1].len() != 4 {
            return Err(InvalidSplash::BlockLen {
                num: 2,
                expected: 4,
                found: blocks[1].len(),
            });
        }
        if let Some(c) = blocks[1].iter().find(|&x| !valid_byte(*x)) {
            return Err(InvalidSplash::BlockChars {
                num: 2,
                c: char::from(*c),
            });
        }
        // .unwrap() safe due to first check
        let block2 = blocks[1].try_into().unwrap();

        // parse block 3
        if blocks[2].len() != 10 {
            return Err(InvalidSplash::BlockLen {
                num: 3,
                expected: 10,
                found: blocks[2].len(),
            });
        }
        if let Some(c) = blocks[2].iter().find(|&x| !valid_byte(*x)) {
            return Err(InvalidSplash::BlockChars {
                num: 3,
                c: char::from(*c),
            });
        }
        // .unwrap() safe due to first check
        let block3 = blocks[2].try_into().unwrap();

        // parse block 4
        if blocks[3].len() != 20 {
            return Err(InvalidSplash::BlockLen {
                num: 4,
                expected: 20,
                found: blocks[3].len(),
            });
        }
        // .unwrap() safe due to first check
        let block4 = blocks[3].try_into().unwrap();

        Ok(Self {
            block2,
            block3,
            block4,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::Spectrum;
    use crate::SpectrumMultistageMeta;

    use super::*;

    macro_rules! pk {
        ($mz:expr, $int:expr) => {
            Peak::new($mz, $int)
        };
    }

    macro_rules! check_splash {
        ($peaks:expr, $splash:literal) => {
            let spectrum = Spectrum::new(
                $peaks,
                0.0,
                vec![].into(),
                None,
                SpectrumMultistageMeta::default(),
                0,
            );
            assert_eq!(Splash::from(&spectrum).to_string(), $splash,);
        };
    }

    #[test]
    fn compare_to_python() {
        // Example python script for spectrum with two peaks:
        // - Peak { mz: 1.0, intensity: 2.0 }
        // - Peak { mz: 3.0, intensity: 4.0 }
        // ```python
        // from splash import Spectrum, SpectrumType, Splash
        // s = Spectrum([(1.0, 2.0), (3.0, 4.0)], SpectrumType.MS)
        // print(Splash().splash(s))
        // ```

        check_splash!(
            vec![pk!(1.0, 2.0), pk!(3.0, 4.0)],
            "splash10-0udi-9000000000-41d31fbe69c73bb50fba"
        );

        check_splash!(
            vec![pk!(0.0, 2.0), pk!(0.0, 4.0)],
            "splash10-0udi-9000000000-70dd0d40b8ac9f896c90"
        );

        // can't compare this since python gives div zero error
        // TODO: open ticket in og repo
        //     vec![pk!(2.0, 0.0), pk!(3.0, 0.0)],

        check_splash!(
            vec![pk!(2.0, 1.0), pk!(3.0, 0.0)],
            "splash10-0udi-9000000000-6dc5faa64f4abd721a2c"
        );

        check_splash!(
            vec![pk!(3.0, 2000.0), pk!(2.0, 1.0), pk!(3.0, 0.0)],
            "splash10-0udi-9000000000-40098b9a3ef770b487fb"
        );

        check_splash!(
            vec![pk!(3.0, 2000.0), pk!(2.0, 1.0), pk!(3.0, 0.0)],
            "splash10-0udi-9000000000-40098b9a3ef770b487fb"
        );

        check_splash!(
            vec![pk!(22.3, 2000.0), pk!(2.444444, 1.0), pk!(3.2222222, 0.0)],
            "splash10-00di-9000000000-15e3444b483cb05aa0f0"
        );

        check_splash!(
            vec![pk!(3.0, 1.0), pk!(2.0, 1.0), pk!(3.0, 1.0)],
            "splash10-0udi-9000000000-51b9139d60f067b1b0e4"
        );
    }
}
