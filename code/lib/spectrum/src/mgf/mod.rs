/// Local and global MGF keys
///
/// Keys are used in in the `KEY=VALUE` entries in MGF files.
mod keys;
pub use keys::{GlobalMgfKey, LocalMgfKey};

/// Parsers for MGF files
pub mod parsers;
use std::{
    convert::TryFrom,
    fmt,
    ops::{Deref, DerefMut, Range},
    str::FromStr,
};

use crate::{spectrum::SpectrumMultistageMeta, Ms, MsFormat};

pub use self::parsers::parse_mgf;

use super::{ChargeVec, Peak};

/// A single spectrum parsed from an MGF file.
///
/// Can only be constructed as part of an MGF file.
#[derive(Debug, PartialEq, Clone)]
pub struct Spectrum {
    peaks: Vec<Peak>,
    pepmass: f64,
    charges: ChargeVec,
    attrs: Vec<LocalMgfKey>,
}

impl fmt::Display for Spectrum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "BEGIN IONS")?;
        for attr in self.attrs.iter() {
            writeln!(f, "{}", attr)?;
        }
        for peak in self.peaks.iter() {
            writeln!(f, "{} {}", peak.mz, peak.intensity)?;
        }
        writeln!(f, "END IONS")
    }
}

// helper macro to generate getters for MgfSpectrum::attrs
macro_rules! attr_getter {
    ($name:ident, str, $variant:path) => {
        attr_getter!($name, str, $variant, stringify!($variant));
    };
    ($name:ident, str, $variant:path, $svariant:expr) => {
        #[doc = "Get the [`"]
        #[doc = $svariant]
        #[doc = "`] attribute."]
        pub fn $name(&self) -> Option<&str> {
            self.attrs.iter().find_map(|key| match key {
                $variant(x) => Some(x.as_str()),
                _ => None,
            })
        }
    };
    ($name:ident, $inner:ty, $variant:path) => {
        attr_getter!($name, $inner, $variant, stringify!($variant));
    };
    ($name:ident, $inner:ty, $variant:path, $svariant:expr) => {
        #[doc = "Get the [`"]
        #[doc = $svariant]
        #[doc = "`] attribute."]
        pub fn $name(&self) -> Option<&$inner> {
            self.attrs.iter().find_map(|key| match key {
                $variant(x) => Some(x),
                _ => None,
            })
        }
    };
}

impl Spectrum {
    attr_getter!(get_title, str, LocalMgfKey::Title);
    attr_getter!(get_comp, str, LocalMgfKey::Comp);
    attr_getter!(get_instrument, str, LocalMgfKey::Instrument);
    attr_getter!(get_it_mods, str, LocalMgfKey::ItMods);
    attr_getter!(get_rt_in_seconds, Range<f64>, LocalMgfKey::RtInSeconds);
    attr_getter!(get_scans, Range<u32>, LocalMgfKey::Scans);
    attr_getter!(get_tolu, str, LocalMgfKey::Tolu);
    attr_getter!(get_seq, str, LocalMgfKey::Seq);
    attr_getter!(get_tag, str, LocalMgfKey::Tag);
    attr_getter!(get_etag, str, LocalMgfKey::Etag);
    attr_getter!(get_tol, f64, LocalMgfKey::Tol);
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
        self.get_rt_in_seconds().map(|range| range.start as f64)
    }

    fn ms_level(&self) -> Option<usize> {
        None
    }

    fn precursor_scan(&self) -> Option<usize> {
        None
    }

    fn precursor_intensity(&self) -> Option<f64> {
        None
    }

    fn get_scan(&self) -> Option<usize> {
        self.get_scans().map(|scans| scans.start as usize)
    }
}

impl From<Spectrum> for crate::Spectrum {
    fn from(mgf: Spectrum) -> Self {
        let retention_time = mgf.retention_time();
        let scan = mgf.get_scan().unwrap_or_default();
        Self {
            pepmass: mgf.pepmass,
            charges: mgf.charges,
            retention_time,
            multistage_meta: SpectrumMultistageMeta::default(),
            peaks: mgf.peaks,
            scan,
            denormalize_factor: None,
        }
    }
}

impl From<crate::Spectrum> for Spectrum {
    fn from(spectrum: crate::Spectrum) -> Self {
        let mut attrs = vec![
            LocalMgfKey::Charge(spectrum.charges.clone()),
            LocalMgfKey::PepMass(spectrum.pepmass),
        ];
        if let Some(rt) = spectrum.retention_time {
            attrs.push(LocalMgfKey::RtInSeconds(rt..rt + f64::EPSILON));
        }
        Self {
            peaks: spectrum.peaks,
            pepmass: spectrum.pepmass,
            charges: spectrum.charges,
            attrs,
        }
    }
}

// This allows creating of a MgfSpectrum struct without pep_mass, charge, and ion_mode as long as
// they are included in the list of all local variables. This is used in the parser to construct
// the spectrum after parsing a block.
impl TryFrom<(Vec<LocalMgfKey>, Vec<Peak>)> for Spectrum {
    type Error = &'static str;

    fn try_from((attrs, peaks): (Vec<LocalMgfKey>, Vec<Peak>)) -> Result<Self, Self::Error> {
        let pepmass: f64 = attrs
            .iter()
            .find_map(|key| match key {
                LocalMgfKey::PepMass(x) => Some(*x),
                _ => None,
            })
            .ok_or("No value in MgfSpectrum for pepmass")?;

        let charges = attrs
            .iter()
            .find_map(|key| match key {
                LocalMgfKey::Charge(x) => Some(x.clone()),
                _ => None,
            })
            .unwrap_or_else(ChargeVec::new);

        Ok(Spectrum {
            peaks,
            pepmass,
            charges,
            attrs,
        })
    }
}

/// A fully parsed MGF file.
#[derive(Debug, PartialEq)]
pub struct Mgf {
    /// The spectra from this MGF file
    pub spectra: Vec<Spectrum>,
    /// The global-level attributes for this MGF file
    pub attrs: Vec<GlobalMgfKey>,
}

impl fmt::Display for Mgf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for attr in self.attrs.iter() {
            writeln!(f, "{}", attr)?;
        }
        for spec in self.spectra.iter() {
            write!(f, "{}", spec)?;
        }
        Ok(())
    }
}

impl Mgf {
    attr_getter!(get_cle, str, GlobalMgfKey::Cle);
    attr_getter!(get_com, str, GlobalMgfKey::Com);
    attr_getter!(get_db, str, GlobalMgfKey::Db);
    attr_getter!(get_format, str, GlobalMgfKey::Format);
    attr_getter!(get_instrument, str, GlobalMgfKey::Instrument);
    attr_getter!(get_it_mods, str, GlobalMgfKey::ItMods);
    attr_getter!(get_itolu, str, GlobalMgfKey::Itolu);
    attr_getter!(get_mass, str, GlobalMgfKey::Mass);
    attr_getter!(get_mods, str, GlobalMgfKey::Mods);
    attr_getter!(get_quantitation, str, GlobalMgfKey::Quantitation);
    attr_getter!(get_report, str, GlobalMgfKey::Report);
    attr_getter!(get_rep_type, str, GlobalMgfKey::RepType);
    attr_getter!(get_search, str, GlobalMgfKey::Search);
    attr_getter!(get_taxonomy, str, GlobalMgfKey::Taxonomy);
    attr_getter!(get_tolu, str, GlobalMgfKey::Tolu);
    attr_getter!(get_user, str, GlobalMgfKey::User);
    attr_getter!(get_user_email, str, GlobalMgfKey::UserEmail);
    attr_getter!(get_username, str, GlobalMgfKey::UserName);
    attr_getter!(get_decoy, bool, GlobalMgfKey::Decoy);
    attr_getter!(get_error_tolerant, bool, GlobalMgfKey::ErrorTolerant);
    attr_getter!(get_pfa, i32, GlobalMgfKey::Pfa);
    attr_getter!(get_itol, f64, GlobalMgfKey::Itol);
    attr_getter!(get_pep_isotope_error, f64, GlobalMgfKey::PepIsotopeError);
    attr_getter!(get_precursor, f64, GlobalMgfKey::Precursor);
    attr_getter!(get_seg, f64, GlobalMgfKey::Seg);
    attr_getter!(get_tol, f64, GlobalMgfKey::Tol);
    attr_getter!(get_frames, str, GlobalMgfKey::Frames);
    attr_getter!(get_charge, ChargeVec, GlobalMgfKey::Charge);

    /// Iterate over all charges.
    ///
    /// If no global charge key exists this will be an empty iterator.
    pub fn charges(&self) -> impl Iterator<Item = i8> {
        let charge_key = self.attrs.iter().find_map(|key| match key {
            GlobalMgfKey::Charge(x) => Some(x),
            _ => None,
        });

        match charge_key {
            Some(x) => x.clone().into_iter(),
            None => ChargeVec::new().into_iter(),
        }
    }
}

// TODO: document this behavior somewhere
impl From<Mgf> for crate::SpectrumCollection {
    fn from(mgf: Mgf) -> Self {
        let spectra = mgf
            .spectra
            .into_iter()
            .enumerate()
            .map(|(idx, spectrum)| {
                let scan = spectrum.get_scan().unwrap_or(idx);
                let mut spectrum = crate::Spectrum::from(spectrum);
                spectrum.scan = scan;
                spectrum
            })
            .collect();

        Self { spectra }
    }
}

impl FromStr for Mgf {
    type Err = crate::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_mgf(s).map_err(|_| crate::Error::Mgf)
    }
}

impl MsFormat for Mgf {}

#[cfg(test)]
mod tests;
