use super::*;

#[test]
fn get_atomic_num_works() {
    assert_eq!(Atom::H.atomic_num(), 1);
    assert_eq!(Atom::La.atomic_num(), 57);
    assert_eq!(Atom::Dy.atomic_num(), 66);
    assert_eq!(Atom::Br.atomic_num(), 35);
}

#[test]
fn get_standard_atomic_weight_works() {
    assert_eq!(Atom::H.std_mass(), 1.007975);
    assert_eq!(Atom::La.std_mass(), 138.90547);
    assert_eq!(Atom::Na.std_mass(), 22.98976928);
    assert_eq!(Atom::Au.std_mass(), 196.96657);
}

#[test]
fn ord() {
    assert!(Atom::H < Atom::Bk);
    assert!(Atom::C > Atom::H);
    assert!(Atom::Fm <= Atom::Og);
    assert!(Atom::F >= Atom::O);
}
