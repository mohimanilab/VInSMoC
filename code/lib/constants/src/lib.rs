pub const MASS_WATER: f64 = 18.01534;
pub const MASS_WATER_EXACT: f64 = 18.0105646834;
pub const MASS_PROTON: f64 = 1.00727647;

pub mod amino_acids;
/// Mass of a single electron in Dalton
/// This was taken from [the Wikipedia page](https://en.wikipedia.org/wiki/Electron) on Sat Mar 20 20:33:37 EDT 2021
/// TODO: find a more reliable source for this value
pub const MASS_ELECTRON: f64 = 0.00054858;

pub mod atomic_masses;
pub mod standard_atomic_weights;

pub type Da = f64;
pub struct StableAtomicMass {
    pub atomic_number: u32,
    pub num_neutrons: u32,
    pub isotope_freq: f64,
    pub mass: Da,
}

pub struct RadioactiveAtomicMass {
    pub atomic_number: u32,
    pub num_neutrons: u32,
    pub half_life: f64, // in seconds
    pub mass: Da,
}

pub trait AtomicMass {
    fn atomic_num(&self) -> u32;
    fn neutrons(&self) -> u32;
    fn priority(&self) -> f64;
    fn mass(&self) -> Da;
}

impl AtomicMass for StableAtomicMass {
    fn atomic_num(&self) -> u32 {
        self.atomic_number
    }

    fn neutrons(&self) -> u32 {
        self.num_neutrons
    }

    fn priority(&self) -> f64 {
        self.isotope_freq
    }

    fn mass(&self) -> Da {
        self.mass
    }
}

impl AtomicMass for RadioactiveAtomicMass {
    fn atomic_num(&self) -> u32 {
        self.atomic_number
    }
    fn neutrons(&self) -> u32 {
        self.num_neutrons
    }
    fn priority(&self) -> f64 {
        self.half_life
    }
    fn mass(&self) -> Da {
        self.mass
    }
}

pub const PERIODIC_TABLE_LENGTH: usize = 118;
pub const PERIODIC_TABLE: [&[&dyn AtomicMass]; PERIODIC_TABLE_LENGTH] = [
    &atomic_masses::H,
    &atomic_masses::HE,
    &atomic_masses::LI,
    &atomic_masses::BE,
    &atomic_masses::B,
    &atomic_masses::C,
    &atomic_masses::N,
    &atomic_masses::O,
    &atomic_masses::F,
    &atomic_masses::NE,
    &atomic_masses::NA,
    &atomic_masses::MG,
    &atomic_masses::AL,
    &atomic_masses::SI,
    &atomic_masses::P,
    &atomic_masses::S,
    &atomic_masses::CL,
    &atomic_masses::AR,
    &atomic_masses::K,
    &atomic_masses::CA,
    &atomic_masses::SC,
    &atomic_masses::TI,
    &atomic_masses::V,
    &atomic_masses::CR,
    &atomic_masses::MN,
    &atomic_masses::FE,
    &atomic_masses::CO,
    &atomic_masses::NI,
    &atomic_masses::CU,
    &atomic_masses::ZN,
    &atomic_masses::GA,
    &atomic_masses::GE,
    &atomic_masses::AS,
    &atomic_masses::SE,
    &atomic_masses::BR,
    &atomic_masses::KR,
    &atomic_masses::RB,
    &atomic_masses::SR,
    &atomic_masses::Y,
    &atomic_masses::ZR,
    &atomic_masses::NB,
    &atomic_masses::MO,
    &atomic_masses::TC,
    &atomic_masses::RU,
    &atomic_masses::RH,
    &atomic_masses::PD,
    &atomic_masses::AG,
    &atomic_masses::CD,
    &atomic_masses::IN,
    &atomic_masses::SN,
    &atomic_masses::SB,
    &atomic_masses::TE,
    &atomic_masses::I,
    &atomic_masses::XE,
    &atomic_masses::CS,
    &atomic_masses::BA,
    &atomic_masses::LA,
    &atomic_masses::CE,
    &atomic_masses::PR,
    &atomic_masses::ND,
    &atomic_masses::PM,
    &atomic_masses::SM,
    &atomic_masses::EU,
    &atomic_masses::GD,
    &atomic_masses::TB,
    &atomic_masses::DY,
    &atomic_masses::HO,
    &atomic_masses::ER,
    &atomic_masses::TM,
    &atomic_masses::YB,
    &atomic_masses::LU,
    &atomic_masses::HF,
    &atomic_masses::TA,
    &atomic_masses::W,
    &atomic_masses::RE,
    &atomic_masses::OS,
    &atomic_masses::IR,
    &atomic_masses::PT,
    &atomic_masses::AU,
    &atomic_masses::HG,
    &atomic_masses::TL,
    &atomic_masses::PB,
    &atomic_masses::BI,
    &atomic_masses::PO,
    &atomic_masses::AT,
    &atomic_masses::RN,
    &atomic_masses::FR,
    &atomic_masses::RA,
    &atomic_masses::AC,
    &atomic_masses::TH,
    &atomic_masses::PA,
    &atomic_masses::U,
    &atomic_masses::NP,
    &atomic_masses::PU,
    &atomic_masses::AM,
    &atomic_masses::CM,
    &atomic_masses::BK,
    &atomic_masses::CF,
    &atomic_masses::ES,
    &atomic_masses::FM,
    &atomic_masses::MD,
    &atomic_masses::NO,
    &atomic_masses::LR,
    &atomic_masses::RF,
    &atomic_masses::DB,
    &atomic_masses::SG,
    &atomic_masses::BH,
    &atomic_masses::HS,
    &atomic_masses::MT,
    &atomic_masses::DS,
    &atomic_masses::RG,
    &atomic_masses::CN,
    &atomic_masses::NH,
    &atomic_masses::FL,
    &atomic_masses::MC,
    &atomic_masses::LV,
    &atomic_masses::TS,
    &atomic_masses::OG,
];

pub const STANDARD_ATOMIC_WEIGHTS: [Da; PERIODIC_TABLE_LENGTH] = [
    standard_atomic_weights::H,
    standard_atomic_weights::HE,
    standard_atomic_weights::LI,
    standard_atomic_weights::BE,
    standard_atomic_weights::B,
    standard_atomic_weights::C,
    standard_atomic_weights::N,
    standard_atomic_weights::O,
    standard_atomic_weights::F,
    standard_atomic_weights::NE,
    standard_atomic_weights::NA,
    standard_atomic_weights::MG,
    standard_atomic_weights::AL,
    standard_atomic_weights::SI,
    standard_atomic_weights::P,
    standard_atomic_weights::S,
    standard_atomic_weights::CL,
    standard_atomic_weights::AR,
    standard_atomic_weights::K,
    standard_atomic_weights::CA,
    standard_atomic_weights::SC,
    standard_atomic_weights::TI,
    standard_atomic_weights::V,
    standard_atomic_weights::CR,
    standard_atomic_weights::MN,
    standard_atomic_weights::FE,
    standard_atomic_weights::CO,
    standard_atomic_weights::NI,
    standard_atomic_weights::CU,
    standard_atomic_weights::ZN,
    standard_atomic_weights::GA,
    standard_atomic_weights::GE,
    standard_atomic_weights::AS,
    standard_atomic_weights::SE,
    standard_atomic_weights::BR,
    standard_atomic_weights::KR,
    standard_atomic_weights::RB,
    standard_atomic_weights::SR,
    standard_atomic_weights::Y,
    standard_atomic_weights::ZR,
    standard_atomic_weights::NB,
    standard_atomic_weights::MO,
    standard_atomic_weights::TC,
    standard_atomic_weights::RU,
    standard_atomic_weights::RH,
    standard_atomic_weights::PD,
    standard_atomic_weights::AG,
    standard_atomic_weights::CD,
    standard_atomic_weights::IN,
    standard_atomic_weights::SN,
    standard_atomic_weights::SB,
    standard_atomic_weights::TE,
    standard_atomic_weights::I,
    standard_atomic_weights::XE,
    standard_atomic_weights::CS,
    standard_atomic_weights::BA,
    standard_atomic_weights::LA,
    standard_atomic_weights::CE,
    standard_atomic_weights::PR,
    standard_atomic_weights::ND,
    standard_atomic_weights::PM,
    standard_atomic_weights::SM,
    standard_atomic_weights::EU,
    standard_atomic_weights::GD,
    standard_atomic_weights::TB,
    standard_atomic_weights::DY,
    standard_atomic_weights::HO,
    standard_atomic_weights::ER,
    standard_atomic_weights::TM,
    standard_atomic_weights::YB,
    standard_atomic_weights::LU,
    standard_atomic_weights::HF,
    standard_atomic_weights::TA,
    standard_atomic_weights::W,
    standard_atomic_weights::RE,
    standard_atomic_weights::OS,
    standard_atomic_weights::IR,
    standard_atomic_weights::PT,
    standard_atomic_weights::AU,
    standard_atomic_weights::HG,
    standard_atomic_weights::TL,
    standard_atomic_weights::PB,
    standard_atomic_weights::BI,
    standard_atomic_weights::PO,
    standard_atomic_weights::AT,
    standard_atomic_weights::RN,
    standard_atomic_weights::FR,
    standard_atomic_weights::RA,
    standard_atomic_weights::AC,
    standard_atomic_weights::TH,
    standard_atomic_weights::PA,
    standard_atomic_weights::U,
    standard_atomic_weights::NP,
    standard_atomic_weights::PU,
    standard_atomic_weights::AM,
    standard_atomic_weights::CM,
    standard_atomic_weights::BK,
    standard_atomic_weights::CF,
    standard_atomic_weights::ES,
    standard_atomic_weights::FM,
    standard_atomic_weights::MD,
    standard_atomic_weights::NO,
    standard_atomic_weights::LR,
    standard_atomic_weights::RF,
    standard_atomic_weights::DB,
    standard_atomic_weights::SG,
    standard_atomic_weights::BH,
    standard_atomic_weights::HS,
    standard_atomic_weights::MT,
    standard_atomic_weights::DS,
    standard_atomic_weights::RG,
    standard_atomic_weights::CN,
    standard_atomic_weights::NH,
    standard_atomic_weights::FL,
    standard_atomic_weights::MC,
    standard_atomic_weights::LV,
    standard_atomic_weights::TS,
    standard_atomic_weights::OG,
];

#[cfg(test)]
mod tests;
