use super::*;

pub const H: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 1,
        num_neutrons: 0,
        isotope_freq: 0.999855,
        mass: 1.0078250322,
    },
    &StableAtomicMass {
        atomic_number: 1,
        num_neutrons: 1,
        isotope_freq: 0.000145,
        mass: 2.0141017781,
    },
];

pub const HE: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 2,
        num_neutrons: 2,
        isotope_freq: 0.999998,
        mass: 4.0026032545,
    },
    &StableAtomicMass {
        atomic_number: 2,
        num_neutrons: 1,
        isotope_freq: 0.000002,
        mass: 3.016029322,
    },
];

pub const LI: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 3,
        num_neutrons: 4,
        isotope_freq: 0.9515,
        mass: 7.01600344,
    },
    &StableAtomicMass {
        atomic_number: 3,
        num_neutrons: 3,
        isotope_freq: 0.0485,
        mass: 6.01512289,
    },
];

pub const BE: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 4,
    num_neutrons: 5,
    isotope_freq: 1.0,
    mass: 9.0121831,
}];

pub const B: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 5,
        num_neutrons: 6,
        isotope_freq: 0.8035,
        mass: 11.00930517,
    },
    &StableAtomicMass {
        atomic_number: 5,
        num_neutrons: 5,
        isotope_freq: 0.1965,
        mass: 10.0129369,
    },
];

pub const C: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 6,
        num_neutrons: 6,
        isotope_freq: 0.9894,
        mass: 12.0,
    },
    &StableAtomicMass {
        atomic_number: 6,
        num_neutrons: 7,
        isotope_freq: 0.0106,
        mass: 13.003354835,
    },
];

pub const N: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 7,
        num_neutrons: 7,
        isotope_freq: 0.996205,
        mass: 14.003074004,
    },
    &StableAtomicMass {
        atomic_number: 7,
        num_neutrons: 8,
        isotope_freq: 0.003795,
        mass: 15.000108899,
    },
];

pub const O: [&dyn AtomicMass; 3] = [
    &StableAtomicMass {
        atomic_number: 8,
        num_neutrons: 8,
        isotope_freq: 0.99757,
        mass: 15.994914619,
    },
    &StableAtomicMass {
        atomic_number: 8,
        num_neutrons: 10,
        isotope_freq: 0.002045,
        mass: 17.999159613,
    },
    &StableAtomicMass {
        atomic_number: 8,
        num_neutrons: 9,
        isotope_freq: 0.0003835,
        mass: 16.999131757,
    },
];

pub const F: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 9,
    num_neutrons: 10,
    isotope_freq: 1.0,
    mass: 18.998403163,
}];

pub const NE: [&dyn AtomicMass; 3] = [
    &StableAtomicMass {
        atomic_number: 10,
        num_neutrons: 10,
        isotope_freq: 0.9048,
        mass: 19.99244018,
    },
    &StableAtomicMass {
        atomic_number: 10,
        num_neutrons: 12,
        isotope_freq: 0.0925,
        mass: 21.9913851,
    },
    &StableAtomicMass {
        atomic_number: 10,
        num_neutrons: 11,
        isotope_freq: 0.0027,
        mass: 20.9938467,
    },
];

pub const NA: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 11,
    num_neutrons: 12,
    isotope_freq: 1.0,
    mass: 22.98976928,
}];

pub const MG: [&dyn AtomicMass; 3] = [
    &StableAtomicMass {
        atomic_number: 12,
        num_neutrons: 12,
        isotope_freq: 0.78965,
        mass: 23.9850417,
    },
    &StableAtomicMass {
        atomic_number: 12,
        num_neutrons: 14,
        isotope_freq: 0.11025,
        mass: 25.982593,
    },
    &StableAtomicMass {
        atomic_number: 12,
        num_neutrons: 13,
        isotope_freq: 0.10011,
        mass: 24.985837,
    },
];

pub const AL: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 13,
    num_neutrons: 14,
    isotope_freq: 1.0,
    mass: 26.9815384,
}];

pub const SI: [&dyn AtomicMass; 3] = [
    &StableAtomicMass {
        atomic_number: 14,
        num_neutrons: 14,
        isotope_freq: 0.922545,
        mass: 27.976926535,
    },
    &StableAtomicMass {
        atomic_number: 14,
        num_neutrons: 15,
        isotope_freq: 0.04672,
        mass: 28.976494665,
    },
    &StableAtomicMass {
        atomic_number: 14,
        num_neutrons: 16,
        isotope_freq: 0.030735,
        mass: 29.9737701,
    },
];

pub const P: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 15,
    num_neutrons: 16,
    isotope_freq: 1.0,
    mass: 30.973761998,
}];

pub const S: [&dyn AtomicMass; 4] = [
    &StableAtomicMass {
        atomic_number: 16,
        num_neutrons: 16,
        isotope_freq: 0.9485,
        mass: 31.972071174,
    },
    &StableAtomicMass {
        atomic_number: 16,
        num_neutrons: 18,
        isotope_freq: 0.04365,
        mass: 33.967867,
    },
    &StableAtomicMass {
        atomic_number: 16,
        num_neutrons: 17,
        isotope_freq: 0.00763,
        mass: 32.97145891,
    },
    &StableAtomicMass {
        atomic_number: 16,
        num_neutrons: 20,
        isotope_freq: 0.000158,
        mass: 35.967081,
    },
];

pub const CL: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 17,
        num_neutrons: 18,
        isotope_freq: 0.758,
        mass: 34.9688527,
    },
    &StableAtomicMass {
        atomic_number: 17,
        num_neutrons: 20,
        isotope_freq: 0.242,
        mass: 36.9659026,
    },
];

pub const AR: [&dyn AtomicMass; 3] = [
    &StableAtomicMass {
        atomic_number: 18,
        num_neutrons: 22,
        isotope_freq: 0.968,
        mass: 39.96238312,
    },
    &StableAtomicMass {
        atomic_number: 18,
        num_neutrons: 20,
        isotope_freq: 0.0215,
        mass: 37.962732,
    },
    &StableAtomicMass {
        atomic_number: 18,
        num_neutrons: 18,
        isotope_freq: 0.01035,
        mass: 35.9675451,
    },
];

pub const K: [&dyn AtomicMass; 3] = [
    &StableAtomicMass {
        atomic_number: 19,
        num_neutrons: 20,
        isotope_freq: 0.932581,
        mass: 38.96370649,
    },
    &StableAtomicMass {
        atomic_number: 19,
        num_neutrons: 22,
        isotope_freq: 0.067302,
        mass: 40.96182526,
    },
    &StableAtomicMass {
        atomic_number: 19,
        num_neutrons: 21,
        isotope_freq: 0.000117,
        mass: 39.9639982,
    },
];

pub const CA: [&dyn AtomicMass; 6] = [
    &StableAtomicMass {
        atomic_number: 20,
        num_neutrons: 20,
        isotope_freq: 0.96941,
        mass: 39.9625909,
    },
    &StableAtomicMass {
        atomic_number: 20,
        num_neutrons: 24,
        isotope_freq: 0.02086,
        mass: 43.955481,
    },
    &StableAtomicMass {
        atomic_number: 20,
        num_neutrons: 22,
        isotope_freq: 0.00647,
        mass: 41.958618,
    },
    &StableAtomicMass {
        atomic_number: 20,
        num_neutrons: 28,
        isotope_freq: 0.00187,
        mass: 47.9525229,
    },
    &StableAtomicMass {
        atomic_number: 20,
        num_neutrons: 23,
        isotope_freq: 0.00135,
        mass: 42.958766,
    },
    &StableAtomicMass {
        atomic_number: 20,
        num_neutrons: 26,
        isotope_freq: 0.00004,
        mass: 45.95369,
    },
];

pub const SC: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 21,
    num_neutrons: 24,
    isotope_freq: 1.0,
    mass: 44.955908,
}];

pub const TI: [&dyn AtomicMass; 5] = [
    &StableAtomicMass {
        atomic_number: 22,
        num_neutrons: 26,
        isotope_freq: 0.7372,
        mass: 47.9479409,
    },
    &StableAtomicMass {
        atomic_number: 22,
        num_neutrons: 24,
        isotope_freq: 0.0825,
        mass: 45.952627,
    },
    &StableAtomicMass {
        atomic_number: 22,
        num_neutrons: 25,
        isotope_freq: 0.0744,
        mass: 46.9517577,
    },
    &StableAtomicMass {
        atomic_number: 22,
        num_neutrons: 27,
        isotope_freq: 0.0541,
        mass: 48.9478646,
    },
    &StableAtomicMass {
        atomic_number: 22,
        num_neutrons: 28,
        isotope_freq: 0.0518,
        mass: 49.9447858,
    },
];

pub const V: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 23,
        num_neutrons: 28,
        isotope_freq: 0.9975,
        mass: 50.943957,
    },
    &StableAtomicMass {
        atomic_number: 23,
        num_neutrons: 27,
        isotope_freq: 0.0025,
        mass: 49.947156,
    },
];

pub const CR: [&dyn AtomicMass; 4] = [
    &StableAtomicMass {
        atomic_number: 24,
        num_neutrons: 28,
        isotope_freq: 0.83789,
        mass: 51.940505,
    },
    &StableAtomicMass {
        atomic_number: 24,
        num_neutrons: 29,
        isotope_freq: 0.09501,
        mass: 52.940647,
    },
    &StableAtomicMass {
        atomic_number: 24,
        num_neutrons: 26,
        isotope_freq: 0.04345,
        mass: 49.946041,
    },
    &StableAtomicMass {
        atomic_number: 24,
        num_neutrons: 30,
        isotope_freq: 0.02365,
        mass: 53.938878,
    },
];

pub const MN: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 25,
    num_neutrons: 30,
    isotope_freq: 1.0,
    mass: 54.938043,
}];

pub const FE: [&dyn AtomicMass; 4] = [
    &StableAtomicMass {
        atomic_number: 26,
        num_neutrons: 30,
        isotope_freq: 0.91754,
        mass: 55.934936,
    },
    &StableAtomicMass {
        atomic_number: 26,
        num_neutrons: 28,
        isotope_freq: 0.05845,
        mass: 53.939608,
    },
    &StableAtomicMass {
        atomic_number: 26,
        num_neutrons: 31,
        isotope_freq: 0.02119,
        mass: 56.935392,
    },
    &StableAtomicMass {
        atomic_number: 26,
        num_neutrons: 32,
        isotope_freq: 0.00282,
        mass: 57.933274,
    },
];

pub const CO: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 27,
    num_neutrons: 32,
    isotope_freq: 1.0,
    mass: 58.933194,
}];

pub const NI: [&dyn AtomicMass; 5] = [
    &StableAtomicMass {
        atomic_number: 28,
        num_neutrons: 30,
        isotope_freq: 0.680769,
        mass: 57.935342,
    },
    &StableAtomicMass {
        atomic_number: 28,
        num_neutrons: 32,
        isotope_freq: 0.262231,
        mass: 59.930785,
    },
    &StableAtomicMass {
        atomic_number: 28,
        num_neutrons: 34,
        isotope_freq: 0.036345,
        mass: 61.928345,
    },
    &StableAtomicMass {
        atomic_number: 28,
        num_neutrons: 33,
        isotope_freq: 0.011399,
        mass: 60.931055,
    },
    &StableAtomicMass {
        atomic_number: 28,
        num_neutrons: 36,
        isotope_freq: 0.009256,
        mass: 63.927966,
    },
];

pub const CU: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 29,
        num_neutrons: 34,
        isotope_freq: 0.6915,
        mass: 62.929597,
    },
    &StableAtomicMass {
        atomic_number: 29,
        num_neutrons: 36,
        isotope_freq: 0.3085,
        mass: 64.92779,
    },
];

pub const ZN: [&dyn AtomicMass; 5] = [
    &StableAtomicMass {
        atomic_number: 30,
        num_neutrons: 34,
        isotope_freq: 0.4917,
        mass: 63.929142,
    },
    &StableAtomicMass {
        atomic_number: 30,
        num_neutrons: 36,
        isotope_freq: 0.2773,
        mass: 65.926034,
    },
    &StableAtomicMass {
        atomic_number: 30,
        num_neutrons: 38,
        isotope_freq: 0.1845,
        mass: 67.924844,
    },
    &StableAtomicMass {
        atomic_number: 30,
        num_neutrons: 37,
        isotope_freq: 0.0404,
        mass: 66.927127,
    },
    &StableAtomicMass {
        atomic_number: 30,
        num_neutrons: 40,
        isotope_freq: 0.0061,
        mass: 69.92532,
    },
];

pub const GA: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 31,
        num_neutrons: 38,
        isotope_freq: 0.60108,
        mass: 68.925573,
    },
    &StableAtomicMass {
        atomic_number: 31,
        num_neutrons: 40,
        isotope_freq: 0.39892,
        mass: 70.924702,
    },
];

pub const GE: [&dyn AtomicMass; 5] = [
    &StableAtomicMass {
        atomic_number: 32,
        num_neutrons: 42,
        isotope_freq: 0.3652,
        mass: 73.92117776,
    },
    &StableAtomicMass {
        atomic_number: 32,
        num_neutrons: 40,
        isotope_freq: 0.2745,
        mass: 71.9220758,
    },
    &StableAtomicMass {
        atomic_number: 32,
        num_neutrons: 38,
        isotope_freq: 0.2052,
        mass: 69.924249,
    },
    &StableAtomicMass {
        atomic_number: 32,
        num_neutrons: 41,
        isotope_freq: 0.0776,
        mass: 72.923459,
    },
    &StableAtomicMass {
        atomic_number: 32,
        num_neutrons: 44,
        isotope_freq: 0.0775,
        mass: 75.9214027,
    },
];

pub const AS: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 33,
    num_neutrons: 42,
    isotope_freq: 1.0,
    mass: 74.921595,
}];

pub const SE: [&dyn AtomicMass; 6] = [
    &StableAtomicMass {
        atomic_number: 34,
        num_neutrons: 46,
        isotope_freq: 0.498,
        mass: 79.916522,
    },
    &StableAtomicMass {
        atomic_number: 34,
        num_neutrons: 44,
        isotope_freq: 0.2369,
        mass: 77.917309,
    },
    &StableAtomicMass {
        atomic_number: 34,
        num_neutrons: 42,
        isotope_freq: 0.0923,
        mass: 75.9192137,
    },
    &StableAtomicMass {
        atomic_number: 34,
        num_neutrons: 48,
        isotope_freq: 0.0882,
        mass: 81.916699,
    },
    &StableAtomicMass {
        atomic_number: 34,
        num_neutrons: 43,
        isotope_freq: 0.076,
        mass: 76.9199141,
    },
    &StableAtomicMass {
        atomic_number: 34,
        num_neutrons: 40,
        isotope_freq: 0.0086,
        mass: 73.9224759,
    },
];

pub const BR: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 35,
        num_neutrons: 44,
        isotope_freq: 0.5065,
        mass: 78.918338,
    },
    &StableAtomicMass {
        atomic_number: 35,
        num_neutrons: 46,
        isotope_freq: 0.4935,
        mass: 80.916288,
    },
];

pub const KR: [&dyn AtomicMass; 6] = [
    &StableAtomicMass {
        atomic_number: 36,
        num_neutrons: 48,
        isotope_freq: 0.56987,
        mass: 83.91149773,
    },
    &StableAtomicMass {
        atomic_number: 36,
        num_neutrons: 50,
        isotope_freq: 0.17279,
        mass: 85.91061063,
    },
    &StableAtomicMass {
        atomic_number: 36,
        num_neutrons: 46,
        isotope_freq: 0.11593,
        mass: 81.91348115,
    },
    &StableAtomicMass {
        atomic_number: 36,
        num_neutrons: 47,
        isotope_freq: 0.115,
        mass: 82.91412652,
    },
    &StableAtomicMass {
        atomic_number: 36,
        num_neutrons: 44,
        isotope_freq: 0.02286,
        mass: 79.916378,
    },
    &StableAtomicMass {
        atomic_number: 36,
        num_neutrons: 42,
        isotope_freq: 0.00355,
        mass: 77.920366,
    },
];

pub const RB: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 37,
        num_neutrons: 48,
        isotope_freq: 0.7217,
        mass: 84.91178974,
    },
    &StableAtomicMass {
        atomic_number: 37,
        num_neutrons: 50,
        isotope_freq: 0.2783,
        mass: 86.90918053,
    },
];

pub const SR: [&dyn AtomicMass; 4] = [
    &StableAtomicMass {
        atomic_number: 38,
        num_neutrons: 50,
        isotope_freq: 0.8258,
        mass: 87.90561226,
    },
    &StableAtomicMass {
        atomic_number: 38,
        num_neutrons: 48,
        isotope_freq: 0.0986,
        mass: 85.90926073,
    },
    &StableAtomicMass {
        atomic_number: 38,
        num_neutrons: 49,
        isotope_freq: 0.07,
        mass: 86.9088775,
    },
    &StableAtomicMass {
        atomic_number: 38,
        num_neutrons: 46,
        isotope_freq: 0.0056,
        mass: 83.913419,
    },
];

pub const Y: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 39,
    num_neutrons: 50,
    isotope_freq: 1.0,
    mass: 88.90584,
}];

pub const ZR: [&dyn AtomicMass; 5] = [
    &StableAtomicMass {
        atomic_number: 40,
        num_neutrons: 50,
        isotope_freq: 0.5145,
        mass: 89.9046988,
    },
    &StableAtomicMass {
        atomic_number: 40,
        num_neutrons: 54,
        isotope_freq: 0.1738,
        mass: 93.906313,
    },
    &StableAtomicMass {
        atomic_number: 40,
        num_neutrons: 52,
        isotope_freq: 0.1715,
        mass: 91.9050353,
    },
    &StableAtomicMass {
        atomic_number: 40,
        num_neutrons: 51,
        isotope_freq: 0.1122,
        mass: 90.9056402,
    },
    &StableAtomicMass {
        atomic_number: 40,
        num_neutrons: 56,
        isotope_freq: 0.028,
        mass: 95.9082776,
    },
];

pub const NB: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 41,
    num_neutrons: 52,
    isotope_freq: 1.0,
    mass: 92.90637,
}];

pub const MO: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 42,
        num_neutrons: 56,
        isotope_freq: 0.24292,
        mass: 97.905404,
    },
    &StableAtomicMass {
        atomic_number: 42,
        num_neutrons: 54,
        isotope_freq: 0.16673,
        mass: 95.9046748,
    },
    &StableAtomicMass {
        atomic_number: 42,
        num_neutrons: 53,
        isotope_freq: 0.15873,
        mass: 94.9058374,
    },
    &StableAtomicMass {
        atomic_number: 42,
        num_neutrons: 50,
        isotope_freq: 0.14649,
        mass: 91.906807,
    },
    &StableAtomicMass {
        atomic_number: 42,
        num_neutrons: 58,
        isotope_freq: 0.09744,
        mass: 99.907468,
    },
    &StableAtomicMass {
        atomic_number: 42,
        num_neutrons: 55,
        isotope_freq: 0.09582,
        mass: 96.906017,
    },
    &StableAtomicMass {
        atomic_number: 42,
        num_neutrons: 52,
        isotope_freq: 0.09187,
        mass: 93.905084,
    },
];

pub const RU: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 44,
        num_neutrons: 58,
        isotope_freq: 0.3155,
        mass: 101.90434,
    },
    &StableAtomicMass {
        atomic_number: 44,
        num_neutrons: 60,
        isotope_freq: 0.1862,
        mass: 103.90543,
    },
    &StableAtomicMass {
        atomic_number: 44,
        num_neutrons: 57,
        isotope_freq: 0.1706,
        mass: 100.905573,
    },
    &StableAtomicMass {
        atomic_number: 44,
        num_neutrons: 55,
        isotope_freq: 0.1276,
        mass: 98.90593,
    },
    &StableAtomicMass {
        atomic_number: 44,
        num_neutrons: 56,
        isotope_freq: 0.126,
        mass: 99.904211,
    },
    &StableAtomicMass {
        atomic_number: 44,
        num_neutrons: 52,
        isotope_freq: 0.0554,
        mass: 95.907589,
    },
    &StableAtomicMass {
        atomic_number: 44,
        num_neutrons: 54,
        isotope_freq: 0.0187,
        mass: 97.90529,
    },
];

pub const RH: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 45,
    num_neutrons: 58,
    isotope_freq: 1.0,
    mass: 102.90549,
}];

pub const PD: [&dyn AtomicMass; 6] = [
    &StableAtomicMass {
        atomic_number: 46,
        num_neutrons: 60,
        isotope_freq: 0.2733,
        mass: 105.90348,
    },
    &StableAtomicMass {
        atomic_number: 46,
        num_neutrons: 62,
        isotope_freq: 0.2646,
        mass: 107.903892,
    },
    &StableAtomicMass {
        atomic_number: 46,
        num_neutrons: 59,
        isotope_freq: 0.2233,
        mass: 104.905079,
    },
    &StableAtomicMass {
        atomic_number: 46,
        num_neutrons: 64,
        isotope_freq: 0.1172,
        mass: 109.905173,
    },
    &StableAtomicMass {
        atomic_number: 46,
        num_neutrons: 58,
        isotope_freq: 0.1114,
        mass: 103.90403,
    },
    &StableAtomicMass {
        atomic_number: 46,
        num_neutrons: 56,
        isotope_freq: 0.0102,
        mass: 101.905632,
    },
];

pub const AG: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 47,
        num_neutrons: 60,
        isotope_freq: 0.51839,
        mass: 106.90509,
    },
    &StableAtomicMass {
        atomic_number: 47,
        num_neutrons: 62,
        isotope_freq: 0.48161,
        mass: 108.904756,
    },
];

pub const CD: [&dyn AtomicMass; 8] = [
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 66,
        isotope_freq: 0.28754,
        mass: 113.903365,
    },
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 64,
        isotope_freq: 0.24109,
        mass: 111.902764,
    },
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 63,
        isotope_freq: 0.12795,
        mass: 110.904184,
    },
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 62,
        isotope_freq: 0.1247,
        mass: 109.903008,
    },
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 65,
        isotope_freq: 0.12227,
        mass: 112.904408,
    },
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 68,
        isotope_freq: 0.07512,
        mass: 115.904763,
    },
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 58,
        isotope_freq: 0.01245,
        mass: 105.90646,
    },
    &StableAtomicMass {
        atomic_number: 48,
        num_neutrons: 60,
        isotope_freq: 0.00888,
        mass: 107.904184,
    },
];

pub const IN: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 49,
        num_neutrons: 66,
        isotope_freq: 0.95719,
        mass: 114.90387877,
    },
    &StableAtomicMass {
        atomic_number: 49,
        num_neutrons: 64,
        isotope_freq: 0.04281,
        mass: 112.90406,
    },
];

pub const SN: [&dyn AtomicMass; 10] = [
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 70,
        isotope_freq: 0.3258,
        mass: 119.902202,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 68,
        isotope_freq: 0.2422,
        mass: 117.901607,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 66,
        isotope_freq: 0.1454,
        mass: 115.9017428,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 69,
        isotope_freq: 0.0859,
        mass: 118.903311,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 67,
        isotope_freq: 0.0768,
        mass: 116.902954,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 74,
        isotope_freq: 0.0579,
        mass: 123.905277,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 72,
        isotope_freq: 0.0463,
        mass: 121.90344,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 62,
        isotope_freq: 0.0097,
        mass: 111.904825,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 64,
        isotope_freq: 0.0066,
        mass: 113.9027801,
    },
    &StableAtomicMass {
        atomic_number: 50,
        num_neutrons: 65,
        isotope_freq: 0.0034,
        mass: 114.9033447,
    },
];

pub const SB: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 51,
        num_neutrons: 70,
        isotope_freq: 0.5721,
        mass: 120.90381,
    },
    &StableAtomicMass {
        atomic_number: 51,
        num_neutrons: 72,
        isotope_freq: 0.4279,
        mass: 122.90421,
    },
];

pub const TE: [&dyn AtomicMass; 8] = [
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 78,
        isotope_freq: 0.3408,
        mass: 129.90622275,
    },
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 76,
        isotope_freq: 0.3174,
        mass: 127.904461,
    },
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 74,
        isotope_freq: 0.1884,
        mass: 125.90331,
    },
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 73,
        isotope_freq: 0.0707,
        mass: 124.90443,
    },
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 72,
        isotope_freq: 0.0474,
        mass: 123.90282,
    },
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 70,
        isotope_freq: 0.0255,
        mass: 121.90304,
    },
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 71,
        isotope_freq: 0.0089,
        mass: 122.90427,
    },
    &StableAtomicMass {
        atomic_number: 52,
        num_neutrons: 68,
        isotope_freq: 0.0009,
        mass: 119.90406,
    },
];

pub const I: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 53,
    num_neutrons: 74,
    isotope_freq: 1.0,
    mass: 126.90447,
}];

pub const XE: [&dyn AtomicMass; 9] = [
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 78,
        isotope_freq: 0.26909,
        mass: 131.90415509,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 75,
        isotope_freq: 0.26401,
        mass: 128.90478086,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 77,
        isotope_freq: 0.21232,
        mass: 130.90508414,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 80,
        isotope_freq: 0.10436,
        mass: 133.90539303,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 82,
        isotope_freq: 0.08857,
        mass: 135.90721448,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 76,
        isotope_freq: 0.04071,
        mass: 129.90350935,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 74,
        isotope_freq: 0.0191,
        mass: 127.903531,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 70,
        isotope_freq: 0.00095,
        mass: 123.90589,
    },
    &StableAtomicMass {
        atomic_number: 54,
        num_neutrons: 72,
        isotope_freq: 0.00089,
        mass: 125.9043,
    },
];

pub const CS: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 55,
    num_neutrons: 78,
    isotope_freq: 1.0,
    mass: 132.90545196,
}];

pub const BA: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 56,
        num_neutrons: 82,
        isotope_freq: 0.717,
        mass: 137.905247,
    },
    &StableAtomicMass {
        atomic_number: 56,
        num_neutrons: 81,
        isotope_freq: 0.1123,
        mass: 136.905827,
    },
    &StableAtomicMass {
        atomic_number: 56,
        num_neutrons: 80,
        isotope_freq: 0.0785,
        mass: 135.904576,
    },
    &StableAtomicMass {
        atomic_number: 56,
        num_neutrons: 79,
        isotope_freq: 0.0659,
        mass: 134.905689,
    },
    &StableAtomicMass {
        atomic_number: 56,
        num_neutrons: 78,
        isotope_freq: 0.0242,
        mass: 133.904508,
    },
    &StableAtomicMass {
        atomic_number: 56,
        num_neutrons: 74,
        isotope_freq: 0.0011,
        mass: 129.90632,
    },
    &StableAtomicMass {
        atomic_number: 56,
        num_neutrons: 76,
        isotope_freq: 0.001,
        mass: 131.905061,
    },
];

pub const LA: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 57,
        num_neutrons: 82,
        isotope_freq: 0.9991119,
        mass: 138.90636,
    },
    &StableAtomicMass {
        atomic_number: 57,
        num_neutrons: 81,
        isotope_freq: 0.0008881,
        mass: 137.90712,
    },
];

pub const CE: [&dyn AtomicMass; 4] = [
    &StableAtomicMass {
        atomic_number: 58,
        num_neutrons: 82,
        isotope_freq: 0.8845,
        mass: 139.90545,
    },
    &StableAtomicMass {
        atomic_number: 58,
        num_neutrons: 84,
        isotope_freq: 0.11114,
        mass: 141.90925,
    },
    &StableAtomicMass {
        atomic_number: 58,
        num_neutrons: 80,
        isotope_freq: 0.00251,
        mass: 137.90599,
    },
    &StableAtomicMass {
        atomic_number: 58,
        num_neutrons: 78,
        isotope_freq: 0.00185,
        mass: 135.907129,
    },
];

pub const PR: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 59,
    num_neutrons: 82,
    isotope_freq: 1.0,
    mass: 140.90766,
}];

pub const ND: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 60,
        num_neutrons: 82,
        isotope_freq: 0.27152,
        mass: 141.90773,
    },
    &StableAtomicMass {
        atomic_number: 60,
        num_neutrons: 84,
        isotope_freq: 0.23798,
        mass: 143.91009,
    },
    &StableAtomicMass {
        atomic_number: 60,
        num_neutrons: 86,
        isotope_freq: 0.17189,
        mass: 145.91312,
    },
    &StableAtomicMass {
        atomic_number: 60,
        num_neutrons: 83,
        isotope_freq: 0.12174,
        mass: 142.90982,
    },
    &StableAtomicMass {
        atomic_number: 60,
        num_neutrons: 85,
        isotope_freq: 0.08293,
        mass: 144.91258,
    },
    &StableAtomicMass {
        atomic_number: 60,
        num_neutrons: 88,
        isotope_freq: 0.05756,
        mass: 147.9169,
    },
    &StableAtomicMass {
        atomic_number: 60,
        num_neutrons: 90,
        isotope_freq: 0.05638,
        mass: 149.920902,
    },
];

pub const SM: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 62,
        num_neutrons: 90,
        isotope_freq: 0.2674,
        mass: 151.919739,
    },
    &StableAtomicMass {
        atomic_number: 62,
        num_neutrons: 92,
        isotope_freq: 0.2274,
        mass: 153.92222,
    },
    &StableAtomicMass {
        atomic_number: 62,
        num_neutrons: 85,
        isotope_freq: 0.15,
        mass: 146.9149,
    },
    &StableAtomicMass {
        atomic_number: 62,
        num_neutrons: 87,
        isotope_freq: 0.1382,
        mass: 148.917191,
    },
    &StableAtomicMass {
        atomic_number: 62,
        num_neutrons: 86,
        isotope_freq: 0.1125,
        mass: 147.91483,
    },
    &StableAtomicMass {
        atomic_number: 62,
        num_neutrons: 88,
        isotope_freq: 0.0737,
        mass: 149.917282,
    },
    &StableAtomicMass {
        atomic_number: 62,
        num_neutrons: 82,
        isotope_freq: 0.0308,
        mass: 143.91201,
    },
];

pub const EU: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 63,
        num_neutrons: 90,
        isotope_freq: 0.5219,
        mass: 152.921237,
    },
    &StableAtomicMass {
        atomic_number: 63,
        num_neutrons: 88,
        isotope_freq: 0.4781,
        mass: 150.919857,
    },
];

pub const GD: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 64,
        num_neutrons: 94,
        isotope_freq: 0.2484,
        mass: 157.924112,
    },
    &StableAtomicMass {
        atomic_number: 64,
        num_neutrons: 96,
        isotope_freq: 0.2186,
        mass: 159.927062,
    },
    &StableAtomicMass {
        atomic_number: 64,
        num_neutrons: 92,
        isotope_freq: 0.2047,
        mass: 155.922131,
    },
    &StableAtomicMass {
        atomic_number: 64,
        num_neutrons: 93,
        isotope_freq: 0.1565,
        mass: 156.923968,
    },
    &StableAtomicMass {
        atomic_number: 64,
        num_neutrons: 91,
        isotope_freq: 0.148,
        mass: 154.92263,
    },
    &StableAtomicMass {
        atomic_number: 64,
        num_neutrons: 90,
        isotope_freq: 0.0218,
        mass: 153.920873,
    },
    &StableAtomicMass {
        atomic_number: 64,
        num_neutrons: 88,
        isotope_freq: 0.002,
        mass: 151.919799,
    },
];

pub const TB: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 65,
    num_neutrons: 94,
    isotope_freq: 1.0,
    mass: 158.925354,
}];

pub const DY: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 66,
        num_neutrons: 98,
        isotope_freq: 0.2826,
        mass: 163.929181,
    },
    &StableAtomicMass {
        atomic_number: 66,
        num_neutrons: 96,
        isotope_freq: 0.25475,
        mass: 161.926804,
    },
    &StableAtomicMass {
        atomic_number: 66,
        num_neutrons: 97,
        isotope_freq: 0.24896,
        mass: 162.928737,
    },
    &StableAtomicMass {
        atomic_number: 66,
        num_neutrons: 95,
        isotope_freq: 0.18889,
        mass: 160.926939,
    },
    &StableAtomicMass {
        atomic_number: 66,
        num_neutrons: 94,
        isotope_freq: 0.02329,
        mass: 159.925203,
    },
    &StableAtomicMass {
        atomic_number: 66,
        num_neutrons: 92,
        isotope_freq: 0.00095,
        mass: 157.92441,
    },
    &StableAtomicMass {
        atomic_number: 66,
        num_neutrons: 90,
        isotope_freq: 0.00056,
        mass: 155.924284,
    },
];

pub const HO: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 67,
    num_neutrons: 98,
    isotope_freq: 1.0,
    mass: 164.930328,
}];

pub const ER: [&dyn AtomicMass; 6] = [
    &StableAtomicMass {
        atomic_number: 68,
        num_neutrons: 98,
        isotope_freq: 0.33503,
        mass: 165.930299,
    },
    &StableAtomicMass {
        atomic_number: 68,
        num_neutrons: 100,
        isotope_freq: 0.26978,
        mass: 167.932376,
    },
    &StableAtomicMass {
        atomic_number: 68,
        num_neutrons: 99,
        isotope_freq: 0.22869,
        mass: 166.932054,
    },
    &StableAtomicMass {
        atomic_number: 68,
        num_neutrons: 102,
        isotope_freq: 0.1491,
        mass: 169.93547,
    },
    &StableAtomicMass {
        atomic_number: 68,
        num_neutrons: 96,
        isotope_freq: 0.01601,
        mass: 163.929207,
    },
    &StableAtomicMass {
        atomic_number: 68,
        num_neutrons: 94,
        isotope_freq: 0.00139,
        mass: 161.928787,
    },
];

pub const TM: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 69,
    num_neutrons: 100,
    isotope_freq: 1.0,
    mass: 168.934218,
}];

pub const YB: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 70,
        num_neutrons: 104,
        isotope_freq: 0.31896,
        mass: 173.93886755,
    },
    &StableAtomicMass {
        atomic_number: 70,
        num_neutrons: 102,
        isotope_freq: 0.21754,
        mass: 171.93638666,
    },
    &StableAtomicMass {
        atomic_number: 70,
        num_neutrons: 103,
        isotope_freq: 0.16098,
        mass: 172.93821622,
    },
    &StableAtomicMass {
        atomic_number: 70,
        num_neutrons: 101,
        isotope_freq: 0.14216,
        mass: 170.93633152,
    },
    &StableAtomicMass {
        atomic_number: 70,
        num_neutrons: 106,
        isotope_freq: 0.12887,
        mass: 175.9425747,
    },
    &StableAtomicMass {
        atomic_number: 70,
        num_neutrons: 100,
        isotope_freq: 0.03023,
        mass: 169.93476725,
    },
    &StableAtomicMass {
        atomic_number: 70,
        num_neutrons: 98,
        isotope_freq: 0.00126,
        mass: 167.933889,
    },
];

pub const LU: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 71,
        num_neutrons: 104,
        isotope_freq: 0.97401,
        mass: 174.940777,
    },
    &StableAtomicMass {
        atomic_number: 71,
        num_neutrons: 105,
        isotope_freq: 0.02599,
        mass: 175.942692,
    },
];

pub const HF: [&dyn AtomicMass; 6] = [
    &StableAtomicMass {
        atomic_number: 72,
        num_neutrons: 108,
        isotope_freq: 0.3512,
        mass: 179.94656,
    },
    &StableAtomicMass {
        atomic_number: 72,
        num_neutrons: 106,
        isotope_freq: 0.2728,
        mass: 177.94371,
    },
    &StableAtomicMass {
        atomic_number: 72,
        num_neutrons: 105,
        isotope_freq: 0.1858,
        mass: 176.94323,
    },
    &StableAtomicMass {
        atomic_number: 72,
        num_neutrons: 107,
        isotope_freq: 0.1363,
        mass: 178.94583,
    },
    &StableAtomicMass {
        atomic_number: 72,
        num_neutrons: 104,
        isotope_freq: 0.0524,
        mass: 175.94141,
    },
    &StableAtomicMass {
        atomic_number: 72,
        num_neutrons: 102,
        isotope_freq: 0.00161,
        mass: 173.94005,
    },
];

pub const TA: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 73,
        num_neutrons: 108,
        isotope_freq: 0.9998824,
        mass: 180.948,
    },
    &StableAtomicMass {
        atomic_number: 73,
        num_neutrons: 107,
        isotope_freq: 0.0001176,
        mass: 179.94747,
    },
];

pub const W: [&dyn AtomicMass; 5] = [
    &StableAtomicMass {
        atomic_number: 74,
        num_neutrons: 110,
        isotope_freq: 0.3064,
        mass: 183.950933,
    },
    &StableAtomicMass {
        atomic_number: 74,
        num_neutrons: 112,
        isotope_freq: 0.2843,
        mass: 185.954365,
    },
    &StableAtomicMass {
        atomic_number: 74,
        num_neutrons: 108,
        isotope_freq: 0.265,
        mass: 181.948206,
    },
    &StableAtomicMass {
        atomic_number: 74,
        num_neutrons: 109,
        isotope_freq: 0.1431,
        mass: 182.950224,
    },
    &StableAtomicMass {
        atomic_number: 74,
        num_neutrons: 106,
        isotope_freq: 0.0012,
        mass: 179.94671,
    },
];

pub const RE: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 75,
        num_neutrons: 112,
        isotope_freq: 0.626,
        mass: 186.955752,
    },
    &StableAtomicMass {
        atomic_number: 75,
        num_neutrons: 110,
        isotope_freq: 0.374,
        mass: 184.952958,
    },
];

pub const OS: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 76,
        num_neutrons: 116,
        isotope_freq: 0.4078,
        mass: 191.96148,
    },
    &StableAtomicMass {
        atomic_number: 76,
        num_neutrons: 114,
        isotope_freq: 0.2626,
        mass: 189.958446,
    },
    &StableAtomicMass {
        atomic_number: 76,
        num_neutrons: 113,
        isotope_freq: 0.1615,
        mass: 188.958146,
    },
    &StableAtomicMass {
        atomic_number: 76,
        num_neutrons: 112,
        isotope_freq: 0.1324,
        mass: 187.955837,
    },
    &StableAtomicMass {
        atomic_number: 76,
        num_neutrons: 111,
        isotope_freq: 0.0196,
        mass: 186.95575,
    },
    &StableAtomicMass {
        atomic_number: 76,
        num_neutrons: 110,
        isotope_freq: 0.0159,
        mass: 185.953838,
    },
    &StableAtomicMass {
        atomic_number: 76,
        num_neutrons: 108,
        isotope_freq: 0.0002,
        mass: 183.952493,
    },
];

pub const IR: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 77,
        num_neutrons: 116,
        isotope_freq: 0.6277,
        mass: 192.962924,
    },
    &StableAtomicMass {
        atomic_number: 77,
        num_neutrons: 114,
        isotope_freq: 0.3723,
        mass: 190.960591,
    },
];

pub const PT: [&dyn AtomicMass; 6] = [
    &StableAtomicMass {
        atomic_number: 78,
        num_neutrons: 117,
        isotope_freq: 0.33775,
        mass: 194.964794,
    },
    &StableAtomicMass {
        atomic_number: 78,
        num_neutrons: 116,
        isotope_freq: 0.32864,
        mass: 193.962683,
    },
    &StableAtomicMass {
        atomic_number: 78,
        num_neutrons: 118,
        isotope_freq: 0.25211,
        mass: 195.964955,
    },
    &StableAtomicMass {
        atomic_number: 78,
        num_neutrons: 120,
        isotope_freq: 0.07356,
        mass: 197.9679,
    },
    &StableAtomicMass {
        atomic_number: 78,
        num_neutrons: 114,
        isotope_freq: 0.00782,
        mass: 191.96104,
    },
    &StableAtomicMass {
        atomic_number: 78,
        num_neutrons: 112,
        isotope_freq: 0.00012,
        mass: 189.95995,
    },
];

pub const AU: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 79,
    num_neutrons: 118,
    isotope_freq: 1.0,
    mass: 196.96657,
}];

pub const HG: [&dyn AtomicMass; 7] = [
    &StableAtomicMass {
        atomic_number: 80,
        num_neutrons: 122,
        isotope_freq: 0.2974,
        mass: 201.970644,
    },
    &StableAtomicMass {
        atomic_number: 80,
        num_neutrons: 120,
        isotope_freq: 0.2314,
        mass: 199.968327,
    },
    &StableAtomicMass {
        atomic_number: 80,
        num_neutrons: 119,
        isotope_freq: 0.1694,
        mass: 198.968281,
    },
    &StableAtomicMass {
        atomic_number: 80,
        num_neutrons: 121,
        isotope_freq: 0.1317,
        mass: 200.970303,
    },
    &StableAtomicMass {
        atomic_number: 80,
        num_neutrons: 118,
        isotope_freq: 0.1004,
        mass: 197.966769,
    },
    &StableAtomicMass {
        atomic_number: 80,
        num_neutrons: 124,
        isotope_freq: 0.0682,
        mass: 203.973494,
    },
    &StableAtomicMass {
        atomic_number: 80,
        num_neutrons: 116,
        isotope_freq: 0.0015,
        mass: 195.96583,
    },
];

pub const TL: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 81,
        num_neutrons: 124,
        isotope_freq: 0.70485,
        mass: 204.974427,
    },
    &StableAtomicMass {
        atomic_number: 81,
        num_neutrons: 122,
        isotope_freq: 0.29515,
        mass: 202.972344,
    },
];

pub const PB: [&dyn AtomicMass; 4] = [
    &StableAtomicMass {
        atomic_number: 82,
        num_neutrons: 126,
        isotope_freq: 0.524,
        mass: 207.976652,
    },
    &StableAtomicMass {
        atomic_number: 82,
        num_neutrons: 124,
        isotope_freq: 0.241,
        mass: 205.974465,
    },
    &StableAtomicMass {
        atomic_number: 82,
        num_neutrons: 125,
        isotope_freq: 0.221,
        mass: 206.975897,
    },
    &StableAtomicMass {
        atomic_number: 82,
        num_neutrons: 122,
        isotope_freq: 0.014,
        mass: 203.973043,
    },
];

pub const BI: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 83,
    num_neutrons: 126,
    isotope_freq: 1.0,
    mass: 208.9804,
}];

pub const TH: [&dyn AtomicMass; 2] = [
    &StableAtomicMass {
        atomic_number: 90,
        num_neutrons: 142,
        isotope_freq: 0.9998,
        mass: 232.03805,
    },
    &StableAtomicMass {
        atomic_number: 90,
        num_neutrons: 140,
        isotope_freq: 0.0002,
        mass: 230.033132,
    },
];

pub const PA: [&dyn AtomicMass; 1] = [&StableAtomicMass {
    atomic_number: 91,
    num_neutrons: 140,
    isotope_freq: 1.0,
    mass: 231.03588,
}];

pub const U: [&dyn AtomicMass; 3] = [
    &StableAtomicMass {
        atomic_number: 92,
        num_neutrons: 146,
        isotope_freq: 0.992742,
        mass: 238.05079,
    },
    &StableAtomicMass {
        atomic_number: 92,
        num_neutrons: 143,
        isotope_freq: 0.007204,
        mass: 235.043928,
    },
    &StableAtomicMass {
        atomic_number: 92,
        num_neutrons: 142,
        isotope_freq: 0.000054,
        mass: 234.04095,
    },
];

pub const TC: [&dyn AtomicMass; 2] = [
    &RadioactiveAtomicMass {
        atomic_number: 43,
        num_neutrons: 54,
        half_life: 132853858560000.0,
        mass: 97.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 43,
        num_neutrons: 55,
        half_life: 132538291200000.0,
        mass: 98.0,
    },
];

pub const PM: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 61,
    num_neutrons: 84,
    half_life: 558554227.1999999,
    mass: 145.0,
}];

pub const PO: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 84,
    num_neutrons: 125,
    half_life: 3913035264.0,
    mass: 209.0,
}];

pub const AT: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 85,
    num_neutrons: 125,
    half_life: 29160.0,
    mass: 210.0,
}];

pub const RN: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 86,
    num_neutrons: 136,
    half_life: 330177.6,
    mass: 222.0,
}];

pub const FR: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 87,
    num_neutrons: 136,
    half_life: 1320.0,
    mass: 223.0,
}];

pub const RA: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 88,
    num_neutrons: 138,
    half_life: 50490777600.0,
    mass: 226.0,
}];

pub const AC: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 89,
    num_neutrons: 138,
    half_life: 687053256.1919999,
    mass: 227.0,
}];

pub const NP: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 93,
    num_neutrons: 144,
    half_life: 67657641984000.0,
    mass: 237.0,
}];

pub const PU: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 94,
    num_neutrons: 150,
    half_life: 2524538880000000.0,
    mass: 244.0,
}];

pub const AM: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 95,
    num_neutrons: 148,
    half_life: 232383803904.0,
    mass: 243.0,
}];

pub const CM: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 96,
    num_neutrons: 151,
    half_life: 492285081600000.0,
    mass: 247.0,
}];

pub const BK: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 97,
    num_neutrons: 150,
    half_life: 43548295680.0,
    mass: 247.0,
}];

pub const CF: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 98,
    num_neutrons: 153,
    half_life: 28401062400.0,
    mass: 251.0,
}];

pub const ES: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 99,
    num_neutrons: 153,
    half_life: 40754880.0,
    mass: 252.0,
}];

pub const FM: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 100,
    num_neutrons: 157,
    half_life: 8683200.0,
    mass: 257.0,
}];

pub const MD: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 101,
    num_neutrons: 157,
    half_life: 4449600.0,
    mass: 258.0,
}];

pub const NO: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 102,
    num_neutrons: 157,
    half_life: 3480.0,
    mass: 259.0,
}];

pub const LR: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 103,
    num_neutrons: 159,
    half_life: 14400.0,
    mass: 262.0,
}];

pub const RF: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 104,
    num_neutrons: 163,
    half_life: 9000.0,
    mass: 267.0,
}];

pub const DB: [&dyn AtomicMass; 2] = [
    &RadioactiveAtomicMass {
        atomic_number: 105,
        num_neutrons: 163,
        half_life: 104400.0,
        mass: 268.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 105,
        num_neutrons: 165,
        half_life: 7200.0,
        mass: 270.0,
    },
];

pub const SG: [&dyn AtomicMass; 2] = [
    &RadioactiveAtomicMass {
        atomic_number: 106,
        num_neutrons: 163,
        half_life: 300.0,
        mass: 269.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 106,
        num_neutrons: 165,
        half_life: 186.0,
        mass: 271.0,
    },
];

pub const BH: [&dyn AtomicMass; 3] = [
    &RadioactiveAtomicMass {
        atomic_number: 107,
        num_neutrons: 164,
        half_life: 600.0,
        mass: 271.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 107,
        num_neutrons: 163,
        half_life: 228.0,
        mass: 270.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 107,
        num_neutrons: 167,
        half_life: 60.0,
        mass: 274.0,
    },
];

pub const HS: [&dyn AtomicMass; 3] = [
    &RadioactiveAtomicMass {
        atomic_number: 108,
        num_neutrons: 169,
        half_life: 110.0,
        mass: 277.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 108,
        num_neutrons: 161,
        half_life: 16.0,
        mass: 269.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 108,
        num_neutrons: 162,
        half_life: 9.0,
        mass: 270.0,
    },
];

pub const MT: [&dyn AtomicMass; 3] = [
    &RadioactiveAtomicMass {
        atomic_number: 109,
        num_neutrons: 167,
        half_life: 10.0,
        mass: 276.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 109,
        num_neutrons: 168,
        half_life: 9.0,
        mass: 277.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 109,
        num_neutrons: 169,
        half_life: 7.0,
        mass: 278.0,
    },
];

pub const DS: [&dyn AtomicMass; 2] = [
    &RadioactiveAtomicMass {
        atomic_number: 110,
        num_neutrons: 171,
        half_life: 14.0,
        mass: 281.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 110,
        num_neutrons: 170,
        half_life: 11.0,
        mass: 280.0,
    },
];

pub const RG: [&dyn AtomicMass; 2] = [
    &RadioactiveAtomicMass {
        atomic_number: 111,
        num_neutrons: 171,
        half_life: 96.0,
        mass: 282.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 111,
        num_neutrons: 170,
        half_life: 24.0,
        mass: 281.0,
    },
];

pub const CN: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 112,
    num_neutrons: 173,
    half_life: 32.0,
    mass: 285.0,
}];

pub const NH: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 113,
    num_neutrons: 173,
    half_life: 7.0,
    mass: 286.0,
}];

pub const FL: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 114,
    num_neutrons: 175,
    half_life: 2.4,
    mass: 289.0,
}];

pub const MC: [&dyn AtomicMass; 3] = [
    &RadioactiveAtomicMass {
        atomic_number: 115,
        num_neutrons: 175,
        half_life: 0.41,
        mass: 290.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 115,
        num_neutrons: 174,
        half_life: 0.31,
        mass: 289.0,
    },
    &RadioactiveAtomicMass {
        atomic_number: 115,
        num_neutrons: 173,
        half_life: 0.17,
        mass: 288.0,
    },
];

pub const LV: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 116,
    num_neutrons: 177,
    half_life: 0.08,
    mass: 293.0,
}];

pub const TS: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 117,
    num_neutrons: 177,
    half_life: 0.07,
    mass: 294.0,
}];

pub const OG: [&dyn AtomicMass; 1] = [&RadioactiveAtomicMass {
    atomic_number: 118,
    num_neutrons: 176,
    half_life: 0.00115,
    mass: 294.0,
}];
