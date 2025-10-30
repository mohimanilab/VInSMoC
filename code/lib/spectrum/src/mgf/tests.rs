use approx::{AbsDiffEq, UlpsEq};

use super::*;

use std::fs;

impl UlpsEq for Spectrum {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: <f64 as AbsDiffEq>::Epsilon, max_ulps: u32) -> bool {
        if !f64::ulps_eq(&self.pepmass, &other.pepmass, epsilon, max_ulps) {
            return false;
        }

        for (ind, peak) in self.peaks.iter().enumerate() {
            if !(f64::ulps_eq(&peak.mz, &other.peaks[ind].mz, epsilon, max_ulps)
                && f64::ulps_eq(
                    &peak.intensity,
                    &other.peaks[ind].intensity,
                    epsilon,
                    max_ulps,
                ))
            {
                return false;
            }
        }

        true
    }
}

impl AbsDiffEq for Spectrum {
    type Epsilon = <f64 as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> <f64 as AbsDiffEq>::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: <f64 as AbsDiffEq>::Epsilon) -> bool {
        if !f64::abs_diff_eq(&self.pepmass, &other.pepmass, epsilon) {
            return false;
        }

        for (ind, peak) in self.peaks.iter().enumerate() {
            if !(f64::abs_diff_eq(&peak.mz, &other.peaks[ind].mz, epsilon)
                && f64::abs_diff_eq(&peak.intensity, &other.peaks[ind].intensity, epsilon))
            {
                return false;
            }
        }

        true
    }
}

#[test]
fn local_getters() {
    let spec = Spectrum {
        peaks: Vec::new(),
        pepmass: 0.0,
        charges: vec![0].into(),
        attrs: vec![
            LocalMgfKey::Title(String::from("title")),
            LocalMgfKey::Comp(String::from("comp")),
            LocalMgfKey::Instrument(String::from("instrument")),
            LocalMgfKey::ItMods(String::from("itmods")),
            LocalMgfKey::RtInSeconds(23.0..44.2),
            LocalMgfKey::Scans(1..2),
            LocalMgfKey::Tolu(String::from("tolu")),
            LocalMgfKey::Seq(String::from("seq")),
            LocalMgfKey::Tag(String::from("tag")),
            LocalMgfKey::Etag(String::from("etag")),
            LocalMgfKey::Tol(2.0),
            LocalMgfKey::Charge(vec![0].into()),
            LocalMgfKey::PepMass(0.0),
        ],
    };

    assert_eq!(spec.get_title().unwrap(), "title");
    assert_eq!(spec.get_comp().unwrap(), "comp");
    assert_eq!(spec.get_instrument().unwrap(), "instrument");
    assert_eq!(spec.get_it_mods().unwrap(), "itmods");
    assert_eq!(*spec.get_rt_in_seconds().unwrap(), 23.0..44.2);
    assert_eq!(*spec.get_scans().unwrap(), 1u32..2);
    assert_eq!(spec.get_tolu().unwrap(), "tolu");
    assert_eq!(spec.get_seq().unwrap(), "seq");
    assert_eq!(spec.get_tag().unwrap(), "tag");
    assert_eq!(spec.get_etag().unwrap(), "etag");
    assert_eq!(*spec.get_tol().unwrap(), 2.0);
}

#[test]
fn global_getters() {
    let mgf = Mgf {
        spectra: Vec::new(),
        attrs: vec![
            GlobalMgfKey::Cle("cle".to_string()),
            GlobalMgfKey::Com("com".to_string()),
            GlobalMgfKey::Db("db".to_string()),
            GlobalMgfKey::Format("format".to_string()),
            GlobalMgfKey::Instrument("instrument".to_string()),
            GlobalMgfKey::ItMods("itmods".to_string()),
            GlobalMgfKey::Itolu("itolu".to_string()),
            GlobalMgfKey::Mass("mass".to_string()),
            GlobalMgfKey::Mods("mods".to_string()),
            GlobalMgfKey::Quantitation("quantitation".to_string()),
            GlobalMgfKey::Report("report".to_string()),
            GlobalMgfKey::RepType("reptype".to_string()),
            GlobalMgfKey::Search("search".to_string()),
            GlobalMgfKey::Taxonomy("taxonomy".to_string()),
            GlobalMgfKey::Tolu("tolu".to_string()),
            GlobalMgfKey::User("user".to_string()),
            GlobalMgfKey::UserEmail("useremail".to_string()),
            GlobalMgfKey::UserName("username".to_string()),
            GlobalMgfKey::Decoy(false),
            GlobalMgfKey::ErrorTolerant(false),
            GlobalMgfKey::Pfa(1),
            GlobalMgfKey::Itol(1.0),
            GlobalMgfKey::PepIsotopeError(1.0),
            GlobalMgfKey::Precursor(1.0),
            GlobalMgfKey::Seg(1.0),
            GlobalMgfKey::Tol(1.0),
            GlobalMgfKey::Charge(vec![1].into()),
            GlobalMgfKey::Frames("frames".to_string()),
        ],
    };

    assert_eq!(mgf.get_cle().unwrap(), "cle");
    assert_eq!(mgf.get_com().unwrap(), "com");
    assert_eq!(mgf.get_db().unwrap(), "db");
    assert_eq!(mgf.get_format().unwrap(), "format");
    assert_eq!(mgf.get_instrument().unwrap(), "instrument");
    assert_eq!(mgf.get_it_mods().unwrap(), "itmods");
    assert_eq!(mgf.get_itolu().unwrap(), "itolu");
    assert_eq!(mgf.get_mass().unwrap(), "mass");
    assert_eq!(mgf.get_mods().unwrap(), "mods");
    assert_eq!(mgf.get_quantitation().unwrap(), "quantitation");
    assert_eq!(mgf.get_report().unwrap(), "report");
    assert_eq!(mgf.get_rep_type().unwrap(), "reptype");
    assert_eq!(mgf.get_search().unwrap(), "search");
    assert_eq!(mgf.get_taxonomy().unwrap(), "taxonomy");
    assert_eq!(mgf.get_tolu().unwrap(), "tolu");
    assert_eq!(mgf.get_user().unwrap(), "user");
    assert_eq!(mgf.get_user_email().unwrap(), "useremail");
    assert_eq!(mgf.get_username().unwrap(), "username");
    assert_eq!(*mgf.get_decoy().unwrap(), false);
    assert_eq!(*mgf.get_error_tolerant().unwrap(), false);
    assert_eq!(*mgf.get_pfa().unwrap(), 1i32);
    assert_eq!(*mgf.get_itol().unwrap(), 1.0f64);
    assert_eq!(*mgf.get_pep_isotope_error().unwrap(), 1.0f64);
    assert_eq!(*mgf.get_precursor().unwrap(), 1.0f64);
    assert_eq!(*mgf.get_seg().unwrap(), 1.0f64);
    assert_eq!(*mgf.get_tol().unwrap(), 1.0f64);
    assert_eq!(mgf.charges().collect::<Vec<_>>(), vec![1i8]);
    assert_eq!(mgf.get_frames().unwrap(), "frames");
}

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
fn mgf_output_test1() {
    let contents1 =
        fs::read_to_string("test_files/numberformats.mgf").expect("Unable to read file");
    let result = parse_mgf(&contents1).unwrap();
    let parsed_str = result.to_string();
    let reformat = parse_mgf(&parsed_str);
    println!("{}", parsed_str);
    assert_eq!(reformat, Ok(result));
}

#[test]
fn mgf_output_test2() {
    let contents2 = fs::read_to_string("test_files/tmtx.mgf").expect("Unable to read file");
    let result = parse_mgf(&contents2).unwrap();
    let parsed_str = result.to_string();
    let reformat = parse_mgf(&parsed_str);
    println!("{}", parsed_str);
    assert_eq!(reformat, Ok(result));
}

#[test]
fn mgf_output_test3() {
    let contents3 = fs::read_to_string("test_files/gnps.mgf").expect("Unable to read file");
    let result = parse_mgf(&contents3).unwrap();
    let parsed_str = result.to_string();
    let reformat = parse_mgf(&parsed_str);
    println!("{}", parsed_str);
    assert_eq!(reformat, Ok(result));
}

#[test]
fn mgf_output_test4() {
    let contents4 = fs::read_to_string("test_files/block.mgf").expect("Unable to read file");
    let result = parse_mgf(&contents4).unwrap();
    let parsed_str = result.to_string();
    let reformat = parse_mgf(&parsed_str);
    println!("{}", parsed_str);
    assert_eq!(reformat, Ok(result));
}

#[test]
fn casmi_challenge_002() {
    let contents = fs::read_to_string("test_files/Challenge-002.mgf").expect("unable to read file");
    assert!(parse_mgf(&contents).is_ok());
}
