use super::*;
// TODO: use approx to test floating point equality too
extern crate approx;
#[test]
fn hydrogen_works() {
    assert_eq!(PERIODIC_TABLE[0].len(), 2);
    assert_eq!(PERIODIC_TABLE[0][0].atomic_num(), 1);
    assert_eq!(PERIODIC_TABLE[0][0].neutrons(), 0);
    assert_eq!(PERIODIC_TABLE[0][1].atomic_num(), 1);
    assert_eq!(PERIODIC_TABLE[0][1].neutrons(), 1);
    approx::assert_abs_diff_eq!(STANDARD_ATOMIC_WEIGHTS[0], 1.007975, epsilon = 0.00001);
}

macro_rules! check_isotopes_sorted {
    ($a:expr, $b:expr) => {
        for i in 0..($b - 1) {
            assert!($a[i].priority() > $a[i + 1].priority());
        }
    };
}

#[test]
fn periodic_table_sorted() {
    check_isotopes_sorted!(atomic_masses::H, 2);
    check_isotopes_sorted!(atomic_masses::H, 2);
    check_isotopes_sorted!(atomic_masses::HE, 2);
    check_isotopes_sorted!(atomic_masses::LI, 2);
    check_isotopes_sorted!(atomic_masses::BE, 1);
    check_isotopes_sorted!(atomic_masses::B, 2);
    check_isotopes_sorted!(atomic_masses::C, 2);
    check_isotopes_sorted!(atomic_masses::N, 2);
    check_isotopes_sorted!(atomic_masses::O, 3);
    check_isotopes_sorted!(atomic_masses::F, 1);
    check_isotopes_sorted!(atomic_masses::NE, 3);
    check_isotopes_sorted!(atomic_masses::NA, 1);
    check_isotopes_sorted!(atomic_masses::MG, 3);
    check_isotopes_sorted!(atomic_masses::AL, 1);
    check_isotopes_sorted!(atomic_masses::SI, 3);
    check_isotopes_sorted!(atomic_masses::P, 1);
    check_isotopes_sorted!(atomic_masses::S, 4);
    check_isotopes_sorted!(atomic_masses::CL, 2);
    check_isotopes_sorted!(atomic_masses::AR, 3);
    check_isotopes_sorted!(atomic_masses::K, 3);
    check_isotopes_sorted!(atomic_masses::CA, 6);
    check_isotopes_sorted!(atomic_masses::SC, 1);
    check_isotopes_sorted!(atomic_masses::TI, 5);
    check_isotopes_sorted!(atomic_masses::V, 2);
    check_isotopes_sorted!(atomic_masses::CR, 4);
    check_isotopes_sorted!(atomic_masses::MN, 1);
    check_isotopes_sorted!(atomic_masses::FE, 4);
    check_isotopes_sorted!(atomic_masses::CO, 1);
    check_isotopes_sorted!(atomic_masses::NI, 5);
    check_isotopes_sorted!(atomic_masses::CU, 2);
    check_isotopes_sorted!(atomic_masses::ZN, 5);
    check_isotopes_sorted!(atomic_masses::GA, 2);
    check_isotopes_sorted!(atomic_masses::GE, 5);
    check_isotopes_sorted!(atomic_masses::AS, 1);
    check_isotopes_sorted!(atomic_masses::SE, 6);
    check_isotopes_sorted!(atomic_masses::BR, 2);
    check_isotopes_sorted!(atomic_masses::KR, 6);
    check_isotopes_sorted!(atomic_masses::RB, 2);
    check_isotopes_sorted!(atomic_masses::SR, 4);
    check_isotopes_sorted!(atomic_masses::Y, 1);
    check_isotopes_sorted!(atomic_masses::ZR, 5);
    check_isotopes_sorted!(atomic_masses::NB, 1);
    check_isotopes_sorted!(atomic_masses::MO, 7);
    check_isotopes_sorted!(atomic_masses::RU, 7);
    check_isotopes_sorted!(atomic_masses::RH, 1);
    check_isotopes_sorted!(atomic_masses::PD, 6);
    check_isotopes_sorted!(atomic_masses::AG, 2);
    check_isotopes_sorted!(atomic_masses::CD, 8);
    check_isotopes_sorted!(atomic_masses::IN, 2);
    check_isotopes_sorted!(atomic_masses::SN, 10);
    check_isotopes_sorted!(atomic_masses::SB, 2);
    check_isotopes_sorted!(atomic_masses::TE, 8);
    check_isotopes_sorted!(atomic_masses::I, 1);
    check_isotopes_sorted!(atomic_masses::XE, 9);
    check_isotopes_sorted!(atomic_masses::CS, 1);
    check_isotopes_sorted!(atomic_masses::BA, 7);
    check_isotopes_sorted!(atomic_masses::LA, 2);
    check_isotopes_sorted!(atomic_masses::CE, 4);
    check_isotopes_sorted!(atomic_masses::PR, 1);
    check_isotopes_sorted!(atomic_masses::ND, 7);
    check_isotopes_sorted!(atomic_masses::SM, 7);
    check_isotopes_sorted!(atomic_masses::EU, 2);
    check_isotopes_sorted!(atomic_masses::GD, 7);
    check_isotopes_sorted!(atomic_masses::TB, 1);
    check_isotopes_sorted!(atomic_masses::DY, 7);
    check_isotopes_sorted!(atomic_masses::HO, 1);
    check_isotopes_sorted!(atomic_masses::ER, 6);
    check_isotopes_sorted!(atomic_masses::TM, 1);
    check_isotopes_sorted!(atomic_masses::YB, 7);
    check_isotopes_sorted!(atomic_masses::LU, 2);
    check_isotopes_sorted!(atomic_masses::HF, 6);
    check_isotopes_sorted!(atomic_masses::TA, 2);
    check_isotopes_sorted!(atomic_masses::W, 5);
    check_isotopes_sorted!(atomic_masses::RE, 2);
    check_isotopes_sorted!(atomic_masses::OS, 7);
    check_isotopes_sorted!(atomic_masses::IR, 2);
    check_isotopes_sorted!(atomic_masses::PT, 6);
    check_isotopes_sorted!(atomic_masses::AU, 1);
    check_isotopes_sorted!(atomic_masses::HG, 7);
    check_isotopes_sorted!(atomic_masses::TL, 2);
    check_isotopes_sorted!(atomic_masses::PB, 4);
    check_isotopes_sorted!(atomic_masses::BI, 1);
    check_isotopes_sorted!(atomic_masses::TH, 2);
    check_isotopes_sorted!(atomic_masses::PA, 1);
    check_isotopes_sorted!(atomic_masses::U, 3);
    check_isotopes_sorted!(atomic_masses::TC, 2);
    check_isotopes_sorted!(atomic_masses::PM, 1);
    check_isotopes_sorted!(atomic_masses::PO, 1);
    check_isotopes_sorted!(atomic_masses::AT, 1);
    check_isotopes_sorted!(atomic_masses::RN, 1);
    check_isotopes_sorted!(atomic_masses::FR, 1);
    check_isotopes_sorted!(atomic_masses::RA, 1);
    check_isotopes_sorted!(atomic_masses::AC, 1);
    check_isotopes_sorted!(atomic_masses::NP, 1);
    check_isotopes_sorted!(atomic_masses::PU, 1);
    check_isotopes_sorted!(atomic_masses::AM, 1);
    check_isotopes_sorted!(atomic_masses::CM, 1);
    check_isotopes_sorted!(atomic_masses::BK, 1);
    check_isotopes_sorted!(atomic_masses::CF, 1);
    check_isotopes_sorted!(atomic_masses::ES, 1);
    check_isotopes_sorted!(atomic_masses::FM, 1);
    check_isotopes_sorted!(atomic_masses::MD, 1);
    check_isotopes_sorted!(atomic_masses::NO, 1);
    check_isotopes_sorted!(atomic_masses::LR, 1);
    check_isotopes_sorted!(atomic_masses::RF, 1);
    check_isotopes_sorted!(atomic_masses::DB, 2);
    check_isotopes_sorted!(atomic_masses::SG, 2);
    check_isotopes_sorted!(atomic_masses::BH, 3);
    check_isotopes_sorted!(atomic_masses::HS, 3);
    check_isotopes_sorted!(atomic_masses::MT, 3);
    check_isotopes_sorted!(atomic_masses::DS, 2);
    check_isotopes_sorted!(atomic_masses::RG, 2);
    check_isotopes_sorted!(atomic_masses::CN, 1);
    check_isotopes_sorted!(atomic_masses::NH, 1);
    check_isotopes_sorted!(atomic_masses::FL, 1);
    check_isotopes_sorted!(atomic_masses::MC, 3);
    check_isotopes_sorted!(atomic_masses::LV, 1);
    check_isotopes_sorted!(atomic_masses::TS, 1);
    check_isotopes_sorted!(atomic_masses::OG, 1);
}
