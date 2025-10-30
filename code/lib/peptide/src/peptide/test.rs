use approx::assert_abs_diff_eq;
use constants::MASS_WATER;
use molecule::{Mol, MolBuilder};
use petgraph::{algo::is_isomorphic_matching, Graph};
use rustc_hash::FxHashSet;

use crate::{
    amino_acid::{AminoAcid, StdAminoAcid},
    peptide::Peptide,
};

macro_rules! full_mol {
    ($peptide:expr, $smiles:literal) => {
        let real = MolBuilder::from_smiles($smiles).unwrap().build();
        dbg!(&$peptide);
        dbg!($peptide.to_mol().to_smiles());
        assert_eq!(
            real.to_inchi_key().unwrap(),
            $peptide.to_mol().to_inchi_key().unwrap()
        );
    };
}

// this is checked against RDKit's Chem.MolFromFASTA
#[test]
fn to_mol() {
    let empty = Peptide {
        sequence: Vec::new(),
    };
    let empty_mol = empty.to_mol();

    let def_mol = Mol::default();

    assert!(is_isomorphic_matching::<
        &Graph<_, _, _, _>,
        &Graph<_, _, _, _>,
        _,
        _,
    >(
        &empty_mol.clone().graph.into(),
        &def_mol.clone().graph.into(),
        |x, y| x == y,
        |x, y| x == y
    ));

    assert_eq!(empty_mol.hydrogens().clone(), def_mol.hydrogens().clone());

    let pep = vec![StdAminoAcid::R].into_iter().collect::<Peptide>();
    full_mol!(pep, "N=C(N)NCCCC(N)C(=O)O");

    let pep = vec![StdAminoAcid::K, StdAminoAcid::P, StdAminoAcid::G]
        .into_iter()
        .collect::<Peptide>();
    full_mol!(pep, "NCCCCC(N)C(=O)N1CCCC1C(=O)NCC(=O)O");

    let pep = vec![
        StdAminoAcid::P,
        StdAminoAcid::P,
        StdAminoAcid::P,
        StdAminoAcid::P,
    ]
    .into_iter()
    .collect::<Peptide>();
    full_mol!(pep, "O=C(O)C1CCCN1C(=O)C1CCCN1C(=O)C1CCCN1C(=O)C1CCCN1");

    let pep = vec![
        StdAminoAcid::P,
        StdAminoAcid::G,
        StdAminoAcid::P,
        StdAminoAcid::G,
        StdAminoAcid::G,
    ]
    .into_iter()
    .collect::<Peptide>();
    full_mol!(pep, "O=C(O)CNC(=O)CNC(=O)C1CCCN1C(=O)CNC(=O)C1CCCN1");

    let pep = vec![
        StdAminoAcid::G,
        StdAminoAcid::H,
        StdAminoAcid::P,
        StdAminoAcid::H,
        StdAminoAcid::H,
        StdAminoAcid::G,
    ]
    .into_iter()
    .collect::<Peptide>();
    full_mol!(pep, "NCC(=O)NC(CC1=C[NH]C=N1)C(=O)N1CCCC1C(=O)NC(CC1=C[NH]C=N1)C(=O)NC(CC1=C[NH]C=N1)C(=O)NCC(=O)O");

    let pep = vec![
        StdAminoAcid::G,
        StdAminoAcid::G,
        StdAminoAcid::G,
        StdAminoAcid::G,
        StdAminoAcid::G,
    ]
    .into_iter()
    .collect::<Peptide>();
    full_mol!(pep, "NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O");

    let pep = vec![
        StdAminoAcid::G,
        StdAminoAcid::I,
        StdAminoAcid::L,
        StdAminoAcid::K,
        StdAminoAcid::M,
        StdAminoAcid::H,
    ]
    .into_iter()
    .collect::<Peptide>();
    full_mol!(
        pep,
        "CCC(C)C(NC(=O)CN)C(=O)NC(CC(C)C)C(=O)NC(CCCCN)C(=O)NC(CCSC)C(=O)NC(CC1=C[NH]C=N1)C(=O)O"
    );

    let pep = vec![StdAminoAcid::P, StdAminoAcid::E, StdAminoAcid::G]
        .into_iter()
        .collect::<Peptide>();
    full_mol!(pep, "O=C(O)CCC(NC(=O)C1CCCN1)C(=O)NCC(=O)O");

    let pep = vec![
        StdAminoAcid::M,
        StdAminoAcid::A,
        StdAminoAcid::R,
        StdAminoAcid::K,
        StdAminoAcid::S,
    ]
    .into_iter()
    .collect::<Peptide>();
    full_mol!(
        pep,
        "CSCCC(N)C(=O)NC(C)C(=O)NC(CCCNC(=N)N)C(=O)NC(CCCCN)C(=O)NC(CO)C(=O)O"
    );

    let pep = vec![
        StdAminoAcid::L,
        StdAminoAcid::E,
        StdAminoAcid::N,
        StdAminoAcid::I,
        StdAminoAcid::N,
    ]
    .into_iter()
    .collect::<Peptide>();
    full_mol!(
        pep,
        "CCC(C)C(NC(=O)C(CC(N)=O)NC(=O)C(CCC(=O)O)NC(=O)C(N)CC(C)C)C(=O)NC(CC(N)=O)C(=O)O"
    );
}

#[test]
fn get_offset_works() {
    let empty = Peptide {
        sequence: Vec::new(),
    };

    assert!(empty.get_offset(0).is_none());

    let single = Peptide {
        sequence: vec![StdAminoAcid::R.into()],
    };
    assert!(single.get_offset(1).is_none());

    let kpg: Peptide = Peptide {
        sequence: vec![
            StdAminoAcid::K.into(),
            StdAminoAcid::P.into(),
            StdAminoAcid::G.into(),
        ],
    };
    assert_eq!(kpg.get_offset(1).unwrap(), 9);
    assert_eq!(kpg.get_offset(2).unwrap(), 16);
    assert!(kpg.get_offset(3).is_none());
}

#[test]
fn check_peptide_std_mass() {
    let p = Peptide {
        sequence: vec![
            StdAminoAcid::A.into(),
            StdAminoAcid::S.into(),
            StdAminoAcid::T.into(),
            StdAminoAcid::T.into(),
            StdAminoAcid::T.into(),
            StdAminoAcid::N.into(),
            StdAminoAcid::Y.into(),
            StdAminoAcid::T.into(),
        ],
    };
    // If The difference is less than mass of one hydrogen atom, it will pas the test,
    // and difference is assumed as measurement error.
    assert_abs_diff_eq!(
        p.std_mass(),
        p.sequence
            .iter()
            .map(|aa| aa.std_mass() - MASS_WATER)
            .sum::<f64>()
            + MASS_WATER
    );
    assert_abs_diff_eq!(p.exact_mass(), 857.376677, epsilon = 0.1);
}

#[test]
fn sub_peptides() {
    let p = Peptide {
        sequence: vec![
            StdAminoAcid::A.into(),
            StdAminoAcid::S.into(),
            StdAminoAcid::T.into(),
            StdAminoAcid::T.into(),
            StdAminoAcid::T.into(),
            StdAminoAcid::N.into(),
            StdAminoAcid::Y.into(),
            StdAminoAcid::T.into(),
        ],
    };

    // no real length limit
    let ranges = p
        .sub_peptide_ranges(0..p.len() + 1)
        .collect::<FxHashSet<_>>();
    let mut correct_ranges = FxHashSet::default();
    for start in 0..p.len() {
        for end in start..p.len() + 1 {
            correct_ranges.insert(start..end);
        }
    }
    assert_eq!(
        ranges,
        correct_ranges,
        "difference {:?}",
        correct_ranges
            .symmetric_difference(&ranges)
            .collect::<Vec<_>>()
    );

    // exactly len 5
    let ranges = p.sub_peptide_ranges(5..6).collect::<FxHashSet<_>>();
    let mut correct_ranges = FxHashSet::default();
    correct_ranges.insert(0..5);
    correct_ranges.insert(1..6);
    correct_ranges.insert(2..7);
    correct_ranges.insert(3..8);
    assert_eq!(ranges, correct_ranges);

    // lens 3 or 4
    let ranges = p.sub_peptide_ranges(3..5).collect::<FxHashSet<_>>();
    let mut correct_ranges = FxHashSet::default();
    correct_ranges.insert(0..3);
    correct_ranges.insert(0..4);
    correct_ranges.insert(1..4);
    correct_ranges.insert(1..5);
    correct_ranges.insert(2..5);
    correct_ranges.insert(2..6);
    correct_ranges.insert(3..6);
    correct_ranges.insert(3..7);
    correct_ranges.insert(4..7);
    correct_ranges.insert(4..8);
    correct_ranges.insert(5..8);
    assert_eq!(
        ranges,
        correct_ranges,
        "difference {:?}",
        correct_ranges
            .symmetric_difference(&ranges)
            .collect::<Vec<_>>()
    );
}

#[test]
fn to_cyclic() {
    // SMILES here determined by https://www.novoprolabs.com/tools/convert-peptide-to-smiles-string
    use StdAminoAcid::*;

    let reference = MolBuilder::from_smiles("N21[C@@]([H])(CCC1)C(=O)N1[C@@]([H])(CCC1)C(=O)N1[C@@]([H])(CCC1)C(=O)N1[C@@]([H])(CCC1)C2(=O)").unwrap().implicit_hs().build().to_inchi_key().unwrap();
    let ours = vec![P, P, P, P]
        .into_iter()
        .collect::<Peptide>()
        .to_cyclic()
        .to_inchi_key()
        .unwrap();
    assert_eq!(reference, ours, "failed on cyclic PPPP");

    let reference = MolBuilder::from_smiles("N1[C@@]([H])(CO)C(=O)N[C@@]([H])([C@]([H])(CC)C)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])([C@]([H])(CC)C)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CC(C)C)C1(=O)").unwrap().implicit_hs().build().to_inchi_key().unwrap();
    let ours = vec![S, I, C, I, A, L]
        .into_iter()
        .collect::<Peptide>()
        .to_cyclic()
        .to_inchi_key()
        .unwrap();
    assert_eq!(reference, ours, "failed on cyclic SICIAL");

    let reference = MolBuilder::from_smiles("N2[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CC(C)C)C(=O)N1[C@@]([H])(CCC1)C2(=O)").unwrap().implicit_hs().build().to_inchi_key().unwrap();
    let ours = vec![L, L, A, L, L, P]
        .into_iter()
        .collect::<Peptide>()
        .to_cyclic()
        .to_inchi_key()
        .unwrap();
    assert_eq!(reference, ours, "failed on cyclic LLALLP");

    let reference = MolBuilder::from_smiles("N2[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)N[C@@]([H])(CC(C)C)C(=O)N1[C@@]([H])(CCC1)C2(=O)").unwrap().implicit_hs().build().to_inchi_key().unwrap();
    let ours = vec![H, E, L, P]
        .into_iter()
        .collect::<Peptide>()
        .to_cyclic()
        .to_inchi_key()
        .unwrap();
    assert_eq!(reference, ours, "failed on cyclic HELP");
}

#[test]
fn cyclic_lactone_ranges() {
    use StdAminoAcid::*;

    let no_negatives = vec![A, A, A, A].into_iter().collect::<Peptide>();
    assert!(
        no_negatives.cyclic_lactone_ranges(0..5).next().is_none(),
        "failed to ignore region containing no negatively charged amino acids"
    );

    let exactly_third_negative = vec![A, D, A, A].into_iter().collect::<Peptide>();
    assert!(
        exactly_third_negative.cyclic_lactone_ranges(0..1).next() == Some(0..1),
        "failed to recover core upstream 1/3 negatively charged region"
    );

    let lt_third_negative = vec![A, D, A, A, A].into_iter().collect::<Peptide>();
    assert!(
        lt_third_negative
            .cyclic_lactone_ranges(0..1)
            .next()
            .is_none(),
        "failed to reject less than 1/3 negatively charged region"
    );

    let reject_one_region = vec![A, D, A, A, A, A, D, A]
        .into_iter()
        .collect::<Peptide>();
    assert!(
        reject_one_region.cyclic_lactone_ranges(0..1).count() == 1,
        "failed to reject a larger region with less than 1/3 negatively charged residues"
    );

    let real_aip = vec![
        M, E, N, I, F, N, L, F, I, K, F, F, T, T, I, L, E, F, I, G, T, V, A, G, D, S, V, C, A, S,
        Y, F, D, E, P, E, V, P, E, E, L, T, K, L, Y, E,
    ]
    .into_iter()
    .collect::<Peptide>();
    let expected_core = vec![D, S, V, C, A, S, Y, F]
        .into_iter()
        .collect::<Peptide>();
    assert!(
        real_aip
            .cyclic_lactone_ranges(7..8)
            .any(|range| real_aip.sub_peptide(range) == expected_core),
        "failed to recover real_aip core"
    );
}
