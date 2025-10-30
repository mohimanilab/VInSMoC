use smallvec::SmallVec;

use crate::mgf::keys::LocalMgfKey;

use super::*;
use std::fs;

#[test]
fn comment_works() {
    assert_eq!(comment("# asdfkajsdlf\n"), Ok(("", " asdfkajsdlf")));
    assert!(comment("# asdfkajsdlf").is_err());
    assert_eq!(comment("     # asdfkajsdlf\n"), Ok(("", " asdfkajsdlf")));
}

#[test]
fn ion_works() {
    assert_eq!(
        ion("110.7 109.3\r\n"),
        Ok((
            "",
            (Peak {
                mz: 110.7,
                intensity: 109.3
            })
        ))
    );
    assert_eq!(
        ion("285.1643 1.324E+06\r\n"),
        Ok((
            "",
            (Peak {
                mz: 285.1643,
                intensity: 1324000.0
            })
        ))
    );
    assert!(ion("1.324E+06\r\n").is_err());
    assert!(ion("285.1643 1.324E+06").is_err());
}

#[test]
fn ions_works() {
    let s1 =
        "232.1765 190.9556\r\n233.1055 212.1309\r\n261.3891 218.4117\r\n285.1643 1.324E+06\r\n";
    assert_eq!(
        ions(s1),
        Ok((
            "",
            (vec![
                Peak {
                    mz: 232.1765,
                    intensity: 190.9556
                },
                Peak {
                    mz: 233.1055,
                    intensity: 212.1309
                },
                Peak {
                    mz: 261.3891,
                    intensity: 218.4117
                },
                Peak {
                    mz: 285.1643,
                    intensity: 1324000.0
                }
            ])
        ))
    );

    let s2 =
        "232.1765\t    190.9556\r\n233.1055    212.1309\r\n261.3891\t218.4117\r\n285.1643 1.324E+06\r\n";
    assert_eq!(
        ions(s2),
        Ok((
            "",
            (vec![
                Peak {
                    mz: 232.1765,
                    intensity: 190.9556
                },
                Peak {
                    mz: 233.1055,
                    intensity: 212.1309
                },
                Peak {
                    mz: 261.3891,
                    intensity: 218.4117
                },
                Peak {
                    mz: 285.1643,
                    intensity: 1324000.0
                }
            ])
        ))
    );

    let s3 = "1.0   2.0    \n3.0   4.0\n";
    assert_eq!(
        ions(s3),
        Ok((
            "",
            vec![
                Peak {
                    mz: 1.0,
                    intensity: 2.0
                },
                Peak {
                    mz: 3.0,
                    intensity: 4.0
                },
            ]
        ))
    );
}

#[test]
fn begin_ions_works() {
    assert_eq!(begin_ions("BEGIN IONS\n"), Ok(("", "BEGIN IONS")));
    assert!(begin_ions("BEGIN IONS").is_err());
}

#[test]
fn block_works() {
    let contents1 = fs::read_to_string("test_files/block.mgf").expect("Unable to read file");

    let peaks: Vec<Peak> = vec![
        Peak {
            mz: 232.1765,
            intensity: 190.9556,
        },
        Peak {
            mz: 233.1055,
            intensity: 212.1309,
        },
        Peak {
            mz: 261.3891,
            intensity: 218.4117,
        },
        Peak {
            mz: 275.1903,
            intensity: 482.2604,
        },
        Peak {
            mz: 285.1643,
            intensity: 1.324E+06,
        },
        Peak {
            mz: 288.1907,
            intensity: 0.0215,
        },
        Peak {
            mz: 320.1273,
            intensity: 552.8523,
        },
        Peak {
            mz: 335.1647,
            intensity: 274.4893,
        },
        Peak {
            mz: 342.2788,
            intensity: 218.6740,
        },
    ];
    let pep_mass = 508.7576;
    let answer = Spectrum {
        peaks,
        pepmass: pep_mass,
        charges: vec![4].into(),
        attrs: vec![
            LocalMgfKey::Title(String::from("MGFp format test")),
            LocalMgfKey::Charge(vec![4].into()),
            LocalMgfKey::PepMass(pep_mass),
        ],
    };

    assert_eq!(answer.get_title(), Some("MGFp format test"));
    assert_eq!(answer.get_tol(), None);
    assert_eq!(answer.get_seq(), None);
    assert_eq!(block(&contents1), Ok(("", answer)));

    let contents2 =
        fs::read_to_string("test_files/missing_pepmass.mgf").expect("Unable to read file");
    assert!(block(&contents2).is_err());
}

#[test]
fn parse_mgf_works() {
    // for single spectra
    let contents1 =
        fs::read_to_string("test_files/numberformats.mgf").expect("Unable to read file");
    let peaks: Vec<Peak> = vec![
        Peak {
            mz: 232.1765,
            intensity: 190.9556,
        },
        Peak {
            mz: 233.1055,
            intensity: 212.1309,
        },
        Peak {
            mz: 261.3891,
            intensity: 218.4117,
        },
        Peak {
            mz: 275.1903,
            intensity: 482.2604,
        },
        Peak {
            mz: 285.1643,
            intensity: 1.324E+06,
        },
        Peak {
            mz: 288.1907,
            intensity: 0.0215,
        },
        Peak {
            mz: 320.1273,
            intensity: 552.8523,
        },
        Peak {
            mz: 335.1647,
            intensity: 274.4893,
        },
        Peak {
            mz: 342.2788,
            intensity: 218.6740,
        },
    ];
    let pep_mass = 508.7576;
    let answer1 = Spectrum {
        peaks,
        pepmass: pep_mass,
        charges: vec![4].into(),
        attrs: vec![
            LocalMgfKey::Title(String::from("MGFp format test")),
            LocalMgfKey::Charge(vec![4].into()),
            LocalMgfKey::PepMass(pep_mass),
            LocalMgfKey::Scans(3..5),
        ],
    };

    let spectra = vec![answer1];
    let spectra_clone = spectra.clone();
    let attrs = vec![];
    assert_eq!(parse_mgf(&contents1), Ok(Mgf { spectra, attrs }));

    // repeat the above test on an identical file with an unsupported attribute
    let spectra = spectra_clone;
    let attrs = vec![];
    let contents1 =
        fs::read_to_string("test_files/unsupported_attr.mgf").expect("Unable to read file");

    assert_eq!(parse_mgf(&contents1), Ok(Mgf { spectra, attrs }));

    // for multiple spectra
    let contents2 = fs::read_to_string("test_files/tmtx.mgf").expect("Unable to read file");

    let peaks1: Vec<Peak> = vec![
        Peak {
            mz: 110.0,
            intensity: 10.4,
        },
        Peak {
            mz: 126.01,
            intensity: 100.00,
        },
        Peak {
            mz: 127.01,
            intensity: 100.00,
        },
        Peak {
            mz: 128.01,
            intensity: 100.00,
        },
        Peak {
            mz: 129.01,
            intensity: 100.00,
        },
        Peak {
            mz: 130.01,
            intensity: 100.00,
        },
        Peak {
            mz: 131.01,
            intensity: 100.00,
        },
        Peak {
            mz: 232.1765,
            intensity: 190.9556,
        },
        Peak {
            mz: 233.1055,
            intensity: 212.1309,
        },
        Peak {
            mz: 261.3891,
            intensity: 218.4117,
        },
        Peak {
            mz: 275.1903,
            intensity: 482.2604,
        },
        Peak {
            mz: 285.1643,
            intensity: 1.324E+06,
        },
        Peak {
            mz: 288.1907,
            intensity: 0.0215,
        },
        Peak {
            mz: 320.1273,
            intensity: 552.8523,
        },
        Peak {
            mz: 335.1647,
            intensity: 274.4893,
        },
        Peak {
            mz: 342.2788,
            intensity: 218.6740,
        },
    ];
    let pep_mass1 = 508.7576;
    let charge1 = 4;
    let spectrum1 = Spectrum {
        peaks: peaks1,
        pepmass: pep_mass1,
        charges: vec![charge1].into(),
        attrs: vec![
            LocalMgfKey::Title(String::from("tmt extraction test 01")),
            LocalMgfKey::Charge(vec![charge1].into()),
            LocalMgfKey::PepMass(pep_mass1),
        ],
    };

    let peaks2: Vec<Peak> = vec![
        Peak {
            mz: 232.1765,
            intensity: 190.9556,
        },
        Peak {
            mz: 233.1055,
            intensity: 212.1309,
        },
        Peak {
            mz: 261.3891,
            intensity: 218.4117,
        },
        Peak {
            mz: 275.1903,
            intensity: 482.2604,
        },
        Peak {
            mz: 155.01,
            intensity: 100.00,
        },
        Peak {
            mz: 156.01,
            intensity: 100.00,
        },
        Peak {
            mz: 157.01,
            intensity: 100.00,
        },
        Peak {
            mz: 158.01,
            intensity: 100.00,
        },
        Peak {
            mz: 160.01,
            intensity: 100.00,
        },
        Peak {
            mz: 285.1643,
            intensity: 1.324E+06,
        },
        Peak {
            mz: 288.1907,
            intensity: 0.0215,
        },
        Peak {
            mz: 320.1273,
            intensity: 552.8523,
        },
        Peak {
            mz: 335.1647,
            intensity: 274.4893,
        },
        Peak {
            mz: 342.2788,
            intensity: 218.6740,
        },
        Peak {
            mz: 126.01,
            intensity: 100.00,
        },
        Peak {
            mz: 127.01,
            intensity: 100.00,
        },
        Peak {
            mz: 129.01,
            intensity: 100.00,
        },
        Peak {
            mz: 130.01,
            intensity: 100.00,
        },
        Peak {
            mz: 131.01,
            intensity: 100.00,
        },
    ];
    let spectrum2 = Spectrum {
        peaks: peaks2,
        pepmass: pep_mass1,
        charges: vec![charge1].into(),
        attrs: vec![
            LocalMgfKey::Title(String::from("tmt extraction test 02")),
            LocalMgfKey::Charge(vec![charge1].into()),
            LocalMgfKey::PepMass(pep_mass1),
        ],
    };

    let peaks3: Vec<Peak> = vec![
        Peak {
            mz: 126.01,
            intensity: 100.00,
        },
        Peak {
            mz: 127.01,
            intensity: 100.00,
        },
        Peak {
            mz: 129.01,
            intensity: 100.00,
        },
        Peak {
            mz: 130.01,
            intensity: 100.00,
        },
        Peak {
            mz: 131.01,
            intensity: 100.00,
        },
        Peak {
            mz: 232.1765,
            intensity: 190.9556,
        },
        Peak {
            mz: 233.1055,
            intensity: 212.1309,
        },
        Peak {
            mz: 261.3891,
            intensity: 218.4117,
        },
        Peak {
            mz: 275.1903,
            intensity: 482.2604,
        },
        Peak {
            mz: 285.1643,
            intensity: 1.324E+06,
        },
        Peak {
            mz: 288.1907,
            intensity: 0.0215,
        },
        Peak {
            mz: 320.1273,
            intensity: 552.8523,
        },
        Peak {
            mz: 335.1647,
            intensity: 274.4893,
        },
        Peak {
            mz: 342.2788,
            intensity: 218.6740,
        },
    ];
    let spectrum3 = Spectrum {
        peaks: peaks3,
        pepmass: pep_mass1,
        charges: vec![charge1].into(),
        attrs: vec![
            LocalMgfKey::Title(String::from("tmt extraction test 03")),
            LocalMgfKey::Charge(vec![charge1].into()),
            LocalMgfKey::PepMass(pep_mass1),
        ],
    };

    let spectra = vec![spectrum1, spectrum2, spectrum3];
    let attrs = vec![
        GlobalMgfKey::Mass(String::from("Monoisotopic")),
        GlobalMgfKey::UserName(String::from("Lou Scene")),
        GlobalMgfKey::Charge(vec![2].into()),
    ];

    let result: Mgf = parse_mgf(&contents2).expect("Did not successfully parse file");
    let result_spectra = result.spectra;
    let result_attrs = result.attrs;

    assert_eq!(spectra.len(), result_spectra.len());
    assert_eq!(spectra[0], result_spectra[0]);
    assert_eq!(spectra[1], result_spectra[1]);
    assert_eq!(spectra[2], result_spectra[2]);

    assert_eq!(attrs, result_attrs);

    assert_eq!(parse_mgf(&contents2), Ok(Mgf { spectra, attrs }));

    let parsed_file = parse_mgf(&contents2).unwrap();
    assert_eq!(parsed_file.get_mass(), Some("Monoisotopic"));
    assert_eq!(parsed_file.get_username(), Some("Lou Scene"));
    assert_eq!(parsed_file.get_tolu(), None);

    let contents3 =
        fs::read_to_string("test_files/invalid_pepmass.mgf").expect("Unable to read file");
    assert!(parse_mgf(&contents3).is_err());

    assert!(parse_mgf(
        &fs::read_to_string("test_files/B3326_R9.mzXML.multiple_charged.mgf").unwrap()
    )
    .is_ok());
}

// 2-pronged bug:
//  1. There were unsupported global keys
//  2. The charge was a vector, not a scalar
#[test]
fn unsupported_global_attr_charge_vec() {
    let content = fs::read_to_string("test_files/cause_error_127.mgf").unwrap();
    let mgf = parse_mgf(&content);
    assert!(mgf.is_ok());
    let mgf = mgf.unwrap();
    let answer = Mgf {
        spectra: vec![Spectrum {
            peaks: vec![Peak {
                mz: 1.0,
                intensity: 1.0,
            }],
            charges: vec![].into(),
            pepmass: 1.0,
            attrs: vec![LocalMgfKey::PepMass(1.0)],
        }],
        attrs: vec![
            GlobalMgfKey::Com("OpenMS_search".to_string()),
            GlobalMgfKey::UserName("OpenMS".to_string()),
            GlobalMgfKey::Format("Mascot generic".to_string()),
            GlobalMgfKey::Tolu("Da".to_string()),
            GlobalMgfKey::Itolu("Da".to_string()),
            // should skip FORMVER
            GlobalMgfKey::Db("MSDB".to_string()),
            GlobalMgfKey::Search("MIS".to_string()),
            GlobalMgfKey::Report("AUTO".to_string()),
            GlobalMgfKey::Cle("Trypsin".to_string()),
            GlobalMgfKey::Mass("monoisotopic".to_string()),
            GlobalMgfKey::Instrument("Default".to_string()),
            GlobalMgfKey::Pfa(1),
            GlobalMgfKey::Tol(3.0),
            GlobalMgfKey::Itol(0.3),
            GlobalMgfKey::Taxonomy("All entries".to_string()),
            GlobalMgfKey::Charge(vec![1, 2, 3].into()),
        ],
    };

    assert_eq!(mgf.spectra.len(), 1);
    assert_eq!(mgf.spectra[0].peaks, answer.spectra[0].peaks);
    assert_eq!(mgf.spectra[0].attrs, answer.spectra[0].attrs);
    assert_eq!(mgf.attrs, answer.attrs);
}

#[test]
fn mgf_parser_failure_199() {
    let contents = fs::read_to_string("test_files/failed_199.mgf").unwrap();
    let mgf = parse_mgf(&contents);
    assert!(mgf.is_ok());
    let mgf = mgf.unwrap();

    assert_eq!(mgf.spectra.len(), 1);
    assert_eq!(mgf.spectra[0].peaks[0], Peak::new(101.0474, 3.6e8));
    assert_eq!(mgf.spectra[0].peaks[8], Peak::new(293.0972, 3.1e10));
    assert!(mgf.attrs.is_empty());
    assert_eq!(
        mgf.spectra[0].attrs,
        vec![
            LocalMgfKey::PepMass(293.0974),
            LocalMgfKey::Scans(1..2),
            LocalMgfKey::RtInSeconds(56.444..56.444),
            LocalMgfKey::Charge(SmallVec::from(vec![1])),
        ]
    );
}

#[test]
fn end_ions_tag() {
    let input = "END IONS\n";
    let (input, output) = end_ions(input).unwrap();
    assert_eq!(input, "");
    assert_eq!(output, "END IONS");

    let input = "END IONS     \n";
    let (input, output) = end_ions(input).unwrap();
    assert_eq!(input, "");
    assert_eq!(output, "END IONS");
}

#[test]
fn begin_ions_tag() {
    let input = "BEGIN IONS\n";
    let (input, output) = begin_ions(input).unwrap();
    assert_eq!(input, "");
    assert_eq!(output, "BEGIN IONS");

    let input = "BEGIN IONS     \n";
    let (input, output) = begin_ions(input).unwrap();
    assert_eq!(input, "");
    assert_eq!(output, "BEGIN IONS");
}

#[test]
fn trailing_whitespace() {
    let contents = r#"BEGIN IONS      
PEPMASS=293.0974  
1.0 2.0    
3.0 4.0     
END IONS
"#;
    eprint!("{}", contents);
    let mgf = parse_mgf(contents);
    assert!(mgf.is_ok());
    let mgf = mgf.unwrap();

    assert_eq!(mgf.spectra.len(), 1);
    assert_eq!(mgf.spectra[0].peaks[0], Peak::new(1.0, 2.0));
    assert_eq!(mgf.spectra[0].peaks[1], Peak::new(3.0, 4.0));
}

#[test]
fn white_space_after_globals() {
    let contents = std::fs::read_to_string("test_files/PLT2_G11.mzML.mgf").unwrap();
    eprint!("{}", contents);
    let mgf = parse_mgf(&contents);
    assert!(mgf.is_ok());
    let mgf = mgf.unwrap();

    assert_eq!(mgf.spectra.len(), 1);
}
