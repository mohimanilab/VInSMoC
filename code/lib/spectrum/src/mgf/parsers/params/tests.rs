use super::*;

#[test]
fn pep_mass_works() {
    assert_eq!(
        peptide_mass_local("PEPMASS=483.2930\r\n193.004 290.3\r\n"),
        Ok(("193.004 290.3\r\n", LocalMgfKey::PepMass(483.2930)))
    );
    assert_eq!(
        peptide_mass_local("PEPMASS=483.2930       \n193.004 290.3\r\n"),
        Ok(("193.004 290.3\r\n", LocalMgfKey::PepMass(483.2930)))
    );
    assert_eq!(
        peptide_mass_local("PEPMASS=483.2930\n193.004 290.3\n"),
        Ok(("193.004 290.3\n", LocalMgfKey::PepMass(483.2930)))
    );
    assert_eq!(
        peptide_mass_local("PEPMASS=483.2930     \t \r\n193.004 290.3\n"),
        Ok(("193.004 290.3\n", LocalMgfKey::PepMass(483.2930)))
    );
    assert!(peptide_mass_local("PEPMASS 483.2930\r\n").is_err());
}

#[test]
fn title_works() {
    assert_eq!(
        title_local("TITLE=MGFp format test\nCHARGE=4+"),
        Ok((
            "CHARGE=4+",
            LocalMgfKey::Title("MGFp format test".to_string())
        ))
    );
    assert_eq!(
        title_local("TITLE=MGFp format test      \nCHARGE=4+"),
        Ok((
            "CHARGE=4+",
            LocalMgfKey::Title("MGFp format test      ".to_string())
        ))
    );
}

#[test]
fn charge_works() {
    assert_eq!(
        charge_local("CHARGE=4+\n"),
        Ok(("", LocalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_local("CHARGE=4+    \n"),
        Ok(("", LocalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_local("CHARGE=4\n"),
        Ok(("", LocalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_local("CHARGE=4-\n"),
        Ok(("", LocalMgfKey::Charge(vec![-4].into())))
    );

    assert_eq!(
        charge_local("CHARGE=+4\n"),
        Ok(("", LocalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_local("CHARGE=4\n"),
        Ok(("", LocalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_local("CHARGE=-4\n"),
        Ok(("", LocalMgfKey::Charge(vec![-4].into())))
    );

    assert_eq!(
        charge_local("CHARGE=1,2,3\n"),
        Ok(("", LocalMgfKey::Charge(vec![1, 2, 3].into())))
    );
    assert_eq!(
        charge_local("CHARGE=1,2,3      \n"),
        Ok(("", LocalMgfKey::Charge(vec![1, 2, 3].into())))
    );

    assert_eq!(
        charge_global("CHARGE=4+\n"),
        Ok(("", GlobalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_global("CHARGE=4\n"),
        Ok(("", GlobalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_global("CHARGE=4-\n"),
        Ok(("", GlobalMgfKey::Charge(vec![-4].into())))
    );

    assert_eq!(
        charge_global("CHARGE=+4\n"),
        Ok(("", GlobalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_global("CHARGE=4\n"),
        Ok(("", GlobalMgfKey::Charge(vec![4].into())))
    );
    assert_eq!(
        charge_global("CHARGE=-4\n"),
        Ok(("", GlobalMgfKey::Charge(vec![-4].into())))
    );

    assert_eq!(
        charge_global("CHARGE=1,2,3\n"),
        Ok(("", GlobalMgfKey::Charge(vec![1, 2, 3].into())))
    );
}
