use std::{cmp::Ordering, fmt, ops::Range};

use crate::ChargeVec;

/// Fields for the global parameters given in MGF file headers.
#[allow(missing_docs)]
#[derive(Debug, PartialEq, Clone)]
pub enum GlobalMgfKey {
    Cle(String),
    Com(String),
    Db(String),
    Format(String),
    Instrument(String),
    ItMods(String),
    Itolu(String),
    Mass(String),
    Mods(String),
    Quantitation(String),
    Report(String),
    RepType(String),
    Search(String),
    Taxonomy(String),
    Tolu(String),
    User(String),
    UserEmail(String),
    UserName(String),
    Decoy(bool),
    ErrorTolerant(bool),
    Pfa(i32),
    Itol(f64),
    PepIsotopeError(f64),
    Precursor(f64),
    Seg(f64),
    Tol(f64),
    Charge(ChargeVec),
    Frames(String),
}

/// Fields for local parameters inside an individual MGF file spectrum.
#[allow(missing_docs)]
#[derive(Debug, PartialEq, Clone)]
pub enum LocalMgfKey {
    Title(String),
    Comp(String),
    Instrument(String),
    ItMods(String),
    RtInSeconds(Range<f64>),
    Scans(Range<u32>),
    Tolu(String),
    Seq(String),
    Tag(String),
    Etag(String),
    Tol(f64),
    Charge(ChargeVec),
    PepMass(f64),
}

impl fmt::Display for LocalMgfKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LocalMgfKey::Title(s) => write!(f, "TITLE={}", s),
            LocalMgfKey::Comp(s) => write!(f, "COMP={}", s),
            LocalMgfKey::Instrument(s) => write!(f, "INSTRUMENT={}", s),
            LocalMgfKey::ItMods(s) => write!(f, "IT_MODS={}", s),
            LocalMgfKey::RtInSeconds(x) => {
                if (x.start - x.end).abs() < f64::EPSILON {
                    write!(f, "RTINSECONDS={}", x.start)
                } else {
                    write!(f, "RTINSECONDS={}-{}", x.start, x.end)
                }
            }
            LocalMgfKey::Scans(x) => {
                if x.start >= x.end - 1 {
                    write!(f, "SCANS={}", x.start)
                } else {
                    write!(f, "SCANS={}-{}", x.start, x.end - 1)
                }
            }
            LocalMgfKey::Tolu(s) => write!(f, "TOLU={}", s),
            LocalMgfKey::Seq(s) => write!(f, "SEQ={}", s),
            LocalMgfKey::Tag(s) => write!(f, "TAG={}", s),
            LocalMgfKey::Etag(s) => write!(f, "ETAG={}", s),
            LocalMgfKey::Tol(x) => write!(f, "TOL={}", x),
            LocalMgfKey::Charge(v) => {
                for (idx, chg) in v.into_iter().enumerate() {
                    if idx == 0 {
                        write!(f, "CHARGE=")?;
                    }
                    match chg.cmp(&0) {
                        Ordering::Less => write!(f, "{}-", chg.abs())?,
                        Ordering::Greater => write!(f, "{}+", chg)?,
                        Ordering::Equal => write!(f, "0")?,
                    };

                    if idx < v.len() - 1 {
                        write!(f, ",")?;
                    }
                }
                Ok(())
            }
            LocalMgfKey::PepMass(m) => write!(f, "PEPMASS={}", m),
        }
    }
}

impl fmt::Display for GlobalMgfKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GlobalMgfKey::Cle(s) => write!(f, "CLE={}", s),
            GlobalMgfKey::Com(s) => write!(f, "COM={}", s),
            GlobalMgfKey::Db(s) => write!(f, "DB={}", s),
            GlobalMgfKey::Format(s) => write!(f, "FORMAT={}", s),
            GlobalMgfKey::Instrument(s) => write!(f, "INSTRUMENT={}", s),
            GlobalMgfKey::ItMods(s) => write!(f, "IT_MODS={}", s),
            GlobalMgfKey::Itolu(s) => write!(f, "ITOLU={}", s),
            GlobalMgfKey::Mass(s) => write!(f, "MASS={}", s),
            GlobalMgfKey::Mods(s) => write!(f, "MODS={}", s),
            GlobalMgfKey::Quantitation(s) => write!(f, "QUANTITATION={}", s),
            GlobalMgfKey::Report(s) => write!(f, "REPORT={}", s),
            GlobalMgfKey::RepType(s) => write!(f, "RETYPE={}", s),
            GlobalMgfKey::Search(s) => write!(f, "SEARCH={}", s),
            GlobalMgfKey::Taxonomy(s) => write!(f, "TAXONOMY={}", s),
            GlobalMgfKey::Tolu(s) => write!(f, "TOLU={}", s),
            GlobalMgfKey::User(s) => write!(f, "USER={}", s),
            GlobalMgfKey::UserEmail(s) => write!(f, "USEREMAIL={}", s),
            GlobalMgfKey::UserName(s) => write!(f, "USERNAME={}", s),
            GlobalMgfKey::Decoy(b) => write!(f, "DECOY={}", *b as u8),
            GlobalMgfKey::ErrorTolerant(b) => write!(f, "ERRORTOLERANT={}", *b as u8),
            GlobalMgfKey::Pfa(v) => write!(f, "PFA={}", *v as i32),
            GlobalMgfKey::Itol(v) => write!(f, "ITOL={}", *v as f64),
            GlobalMgfKey::PepIsotopeError(v) => write!(f, "PEP_ISOTOPE_ERROR={}", *v as f64),
            GlobalMgfKey::Precursor(v) => write!(f, "PRECURSOR={}", *v as f64),
            GlobalMgfKey::Seg(v) => write!(f, "SEG={}", *v as f64),
            GlobalMgfKey::Tol(v) => write!(f, "TOL={}", *v as f64),
            GlobalMgfKey::Charge(v) => {
                write!(f, "CHARGE=")?;
                for (idx, chg) in v.into_iter().enumerate() {
                    match chg.cmp(&0) {
                        Ordering::Less => write!(f, "{}-", chg.abs())?,
                        Ordering::Greater => write!(f, "{}+", chg)?,
                        Ordering::Equal => write!(f, "0")?,
                    };

                    if idx != 0 {
                        write!(f, ",")?;
                    }
                }
                Ok(())
            }
            GlobalMgfKey::Frames(s) => write!(f, "FRAMES={}", s),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn global_key_display() {
        let cle_key = GlobalMgfKey::Cle("cle".to_string());
        let com_key = GlobalMgfKey::Com("com".to_string());
        let db_key = GlobalMgfKey::Db("db".to_string());
        let format_key = GlobalMgfKey::Format("format".to_string());
        let instrument_key = GlobalMgfKey::Instrument("instrument".to_string());
        let itmode_key = GlobalMgfKey::ItMods("itmods".to_string());
        let itolu_key = GlobalMgfKey::Itolu("itolu".to_string());
        let mass_key = GlobalMgfKey::Mass("mass".to_string());
        let mod_key = GlobalMgfKey::Mods("mods".to_string());
        let quant_key = GlobalMgfKey::Quantitation("quantitation".to_string());
        let report_key = GlobalMgfKey::Report("report".to_string());
        let rep_key = GlobalMgfKey::RepType("reptype".to_string());
        let search_key = GlobalMgfKey::Search("search".to_string());
        let tax_key = GlobalMgfKey::Taxonomy("taxonomy".to_string());
        let tolu_key = GlobalMgfKey::Tolu("tolu".to_string());
        let user_key = GlobalMgfKey::User("user".to_string());
        let email_key = GlobalMgfKey::UserEmail("useremail".to_string());
        let name_key = GlobalMgfKey::UserName("username".to_string());
        let decoy_key = GlobalMgfKey::Decoy(false);
        let err_key = GlobalMgfKey::ErrorTolerant(false);
        let pfa_key = GlobalMgfKey::Pfa(1);
        let itol_key = GlobalMgfKey::Itol(1.0);
        let peperr_key = GlobalMgfKey::PepIsotopeError(1.0);
        let precursor_key = GlobalMgfKey::Precursor(1.0);
        let seg_key = GlobalMgfKey::Seg(1.0);
        let tol_key = GlobalMgfKey::Tol(1.0);
        let charge_key = GlobalMgfKey::Charge(vec![1].into());
        let frame_key = GlobalMgfKey::Frames("frames".to_string());

        assert_eq!(com_key.to_string(), "COM=com");
        assert_eq!(cle_key.to_string(), "CLE=cle");
        assert_eq!(db_key.to_string(), "DB=db");
        assert_eq!(format_key.to_string(), "FORMAT=format");
        assert_eq!(instrument_key.to_string(), "INSTRUMENT=instrument");
        assert_eq!(itmode_key.to_string(), "IT_MODS=itmods");
        assert_eq!(itolu_key.to_string(), "ITOLU=itolu");
        assert_eq!(mass_key.to_string(), "MASS=mass");
        assert_eq!(mod_key.to_string(), "MODS=mods");
        assert_eq!(quant_key.to_string(), "QUANTITATION=quantitation");
        assert_eq!(report_key.to_string(), "REPORT=report");
        assert_eq!(rep_key.to_string(), "RETYPE=reptype");
        assert_eq!(search_key.to_string(), "SEARCH=search");
        assert_eq!(tax_key.to_string(), "TAXONOMY=taxonomy");
        assert_eq!(tolu_key.to_string(), "TOLU=tolu");
        assert_eq!(user_key.to_string(), "USER=user");
        assert_eq!(email_key.to_string(), "USEREMAIL=useremail");
        assert_eq!(name_key.to_string(), "USERNAME=username");
        assert_eq!(decoy_key.to_string(), "DECOY=0");
        assert_eq!(err_key.to_string(), "ERRORTOLERANT=0");
        assert_eq!(pfa_key.to_string(), "PFA=1");
        assert_eq!(itol_key.to_string(), "ITOL=1");
        assert_eq!(peperr_key.to_string(), "PEP_ISOTOPE_ERROR=1");
        assert_eq!(precursor_key.to_string(), "PRECURSOR=1");
        assert_eq!(seg_key.to_string(), "SEG=1");
        assert_eq!(tol_key.to_string(), "TOL=1");
        assert_eq!(charge_key.to_string(), "CHARGE=1+");
        assert_eq!(frame_key.to_string(), "FRAMES=frames");
    }

    #[test]
    fn local_key_display() {
        let title_key = LocalMgfKey::Title(String::from("title"));
        let comp_key = LocalMgfKey::Comp(String::from("comp"));
        let instrument_key = LocalMgfKey::Instrument(String::from("instrument"));
        let itmod_key = LocalMgfKey::ItMods(String::from("itmods"));
        let rt_seconds_key = LocalMgfKey::RtInSeconds(23.0..44.2);
        let scans_key = LocalMgfKey::Scans(1..2);
        let tolu_key = LocalMgfKey::Tolu(String::from("tolu"));
        let seq_key = LocalMgfKey::Seq(String::from("seq"));
        let tag_key = LocalMgfKey::Tag(String::from("tag"));
        let etag_key = LocalMgfKey::Etag(String::from("etag"));
        let tol_key = LocalMgfKey::Tol(2.0);
        let charge_key = LocalMgfKey::Charge(vec![0, 1].into());
        let pep_key = LocalMgfKey::PepMass(0.0);

        assert_eq!(title_key.to_string(), "TITLE=title");
        assert_eq!(comp_key.to_string(), "COMP=comp");
        assert_eq!(instrument_key.to_string(), "INSTRUMENT=instrument");
        assert_eq!(itmod_key.to_string(), "IT_MODS=itmods");
        assert_eq!(rt_seconds_key.to_string(), "RTINSECONDS=23-44.2");
        assert_eq!(scans_key.to_string(), "SCANS=1");
        assert_eq!(tolu_key.to_string(), "TOLU=tolu");
        assert_eq!(seq_key.to_string(), "SEQ=seq");
        assert_eq!(tag_key.to_string(), "TAG=tag");
        assert_eq!(etag_key.to_string(), "ETAG=etag");
        assert_eq!(tol_key.to_string(), "TOL=2");
        assert_eq!(charge_key.to_string(), "CHARGE=0,1+");
        assert_eq!(pep_key.to_string(), "PEPMASS=0");
    }

    #[test]
    fn display_empty_charge() {
        let charge = LocalMgfKey::Charge(vec![].into());
        assert_eq!(charge.to_string(), "");
    }
}
