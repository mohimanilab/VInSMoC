use criterion::{black_box, criterion_group, criterion_main, Criterion};

use molecule::parsers::parse_smiles;
use molecule::Mol;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn parse_100_smiles_kekulized(c: &mut Criterion) {
    c.bench_function(
        "parse_smiles 100 kekulized",
        |b| b.iter_with_large_drop(|| {
            parse_smiles(black_box("COC(=O)CC(O)CCC(O)CCCC(C)C")).unwrap();
            parse_smiles(black_box("COC1C(C)OC(OC(C)=O)CC1(C)[N+](=O)[O-]")).unwrap();
            parse_smiles(black_box("CC(O)=NC1=CC(=O)C2(CO)OC2C1O")).unwrap();
            parse_smiles(black_box("CC1=CC2=C(C=C1O)C1=C(C3=CC4=C(C=C(O)C5=CC(C)=C(O)C=C45)C(C)(C)C3=O)C(=O)C(C)(C)C1=CC2=O")).unwrap();
            parse_smiles(black_box("N#CC(N)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("COC1=CC(C=O)=CC=C1OC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC1=CC(CC2=CC=C3C=C(OC)C(O)=NC3=C2OS(=O)(=O)O)=CC=C1O")).unwrap();
            parse_smiles(black_box("C=C1CC(C)CC2CC=CC(CC=CC(=O)OC(C=CC3CC(C)=CCO3)C(O)CC3OC3C(O)C1)O2")).unwrap();
            parse_smiles(black_box("CCCCCCCCCCCCCCCC(=O)OCC1OC(OCC2=C(C)CC(C(C)C3CCC4C5CC6OC67C(O)C=CC(=O)C7(C)C5CCC34C)OC2=O)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("O=C(C=CCCCCCCC=CC1=CC=C2OCOC2=C1)N1CCCCC1")).unwrap();
            parse_smiles(black_box("C=CC1C=C(C)C2C3C(OC4=CC=C(C=C4)CC4(OCC)CC(C(=O)C13)C(O)=N4)C1C(C)CC(C)CC21C")).unwrap();
            parse_smiles(black_box("CCC(C)C1N=C(O)C(CC(C)C)N=C(O)C2CCCN2C(=O)C2CCCN2C1=O")).unwrap();
            parse_smiles(black_box("COC1=CC(CC2CN(C(C)=O)C(CC3=CC=CC=C3)CN2C(C)=O)=CC2=C1OCO2")).unwrap();
            parse_smiles(black_box("COC1C(O)=NC2=CC=C(C=CC3(C)CCC(C(C)(C)O)O3)C(O)=C2C1(O)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("CC1=CC2OC3CC4OC(=O)C=CC=CC(=O)OCCC(C)C(O)C(=O)OCC2(CC1O)C4(C)C31CO1")).unwrap();
            parse_smiles(black_box("CCC1OC2OC3COC(C4=CC=CC=C4)OC3C(=O)C2O1")).unwrap();
            parse_smiles(black_box("COC1=CC2=C(C=C1C1=COC3=C4C=CC(C)(C)OC4=C(OC)C(OC)=C3C1=O)OCO2")).unwrap();
            parse_smiles(black_box("CC(C)CCCCCCCCCCCC1CC(O)=NC(CCC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(O)=NC(CC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(C(C)C)C(=O)O1")).unwrap();
            parse_smiles(black_box("CC(=NNC1=CC=C([N+](=O)[O-])C=C1[N+](=O)[O-])C1=CC=C2OCOC2=C1")).unwrap();
            parse_smiles(black_box("C1=C2CCCCCCCCCCCN3CCC=C(CCCCCCCCCCCN(CC1)C2)C3")).unwrap();
            parse_smiles(black_box("COC1=C(C)C(OC(=O)C2=C(C)C=C(O)C=C2O)=C(C)C(C)=C1C(=O)OC1=C(C)C(C)=C(C(=O)O)C(O)=C1C")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C=C2OC(=O)C=C2C2=CC(Br)=C(O)C(Br)=C2)C=C1")).unwrap();
            parse_smiles(black_box("C=C1C(=O)OC2C=C(C)C(OC(C)=O)CC=C(COC(C)=O)CC(OC(=O)C(C)=CCO)C12")).unwrap();
            parse_smiles(black_box("CC1OC(N2C3=CC=CC=C3C3=CC=C4C5=CC(Cl)=CC=C5[NH]C4=C32)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC1OCC(N=C(C)O)C(O)C1O")).unwrap();
            parse_smiles(black_box("C=CC1(C)CC(OC(=O)C(=C)CO)C(C(=C)C(=O)OC)C(O)C1C(=C)CO")).unwrap();
            parse_smiles(black_box("NCC(CCC(N)C(=O)O)OC1OC(CO)C(O)C(O)C1OC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("CCCCCCCCC=CCCCCCCCC(=O)OCC(COC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)OC(=O)CCCCCCCC=CCCCCCCCC")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C2OC3=CC(OC)=CC(O)=C3C(=O)C2OC2OC(C)C(O)C(O)C2O)C=C1")).unwrap();
            parse_smiles(black_box("CN=C(O)N=C(O)OC")).unwrap();
            parse_smiles(black_box("CC1OC(OC2CCC3(C)C4CCC5(C)C(C6=CC(=O)OC6)CCC5(O)C4CCC3(O)C2)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("O=C(O)CCCCCCCCC=CC(O)=NC1CC2(O)C(=O)C(C1O)C(O)C1OC12")).unwrap();
            parse_smiles(black_box("CONCCCCCCCCCCCCC1=CC=CN=C1")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C(=O)O)C(C(=O)N2CCC3=CC4=C(C=C3C2=O)OCO4)=C1OC")).unwrap();
            parse_smiles(black_box("CNC(C)C1C(O)CC2(C)C3CCC4C(=CC3=CCC12C)CCC(N(C)C(C)=O)C4(C)C")).unwrap();
            parse_smiles(black_box("C1=CC=C(CC2CNC(CC3=CC=CC=C3)CN2)C=C1")).unwrap();
            parse_smiles(black_box("COC1=CC(=O)C2(C3=CC(=O)C4=C(O)C=C(O)C=C4O3)C(CC=C(C)C2CC2(O)C(C)(O)CCC3C(C)(C)CC(O)CC32C)C1=O")).unwrap();
            parse_smiles(black_box("O=C1OC(C(O)CO)C(=NNC2=CC=C([N+](=O)[O-])C=C2[N+](=O)[O-])C1=NNC1=CC=C([N+](=O)[O-])C=C1[N+](=O)[O-]")).unwrap();
            parse_smiles(black_box("CCCCCC=CCC=CCCCC(O)CC(=O)OC1CN(C)C(C(OC2OC(CN)C(O)C2OS(=O)(=O)O)C2OC(N3C=CC(O)=NC3=O)C(O)C2O)C(=O)N(C)C1C(=O)O")).unwrap();
            parse_smiles(black_box("CC1CCC(OC(=O)C2=CC=CC=C2)C2(C)C(OC(=O)C3=CC=CC=C3)C(O)C3C(O)C12OC3(C)C")).unwrap();
            parse_smiles(black_box("CCC1CN2C=CC3=C4[C](C=CC=C4OC)NC3=C2CC1C(=COC)C(=O)OC")).unwrap();
            parse_smiles(black_box("C=CC1C2=CCOC(=O)C2=COC1OC1OC(COC(C)=O)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("C=CC(C)(O)CCC=C(C)C(=O)OCC1OC(OC2OC=C(C(=O)O)C3CCC(C)(O)C23)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC1=C(C=CC(=O)O)C(C(=O)O)=NC2=C1C=CO2")).unwrap();
            parse_smiles(black_box("CC=C(C)C1C(C=CC=C(C)C(O)C(C)CO)C2CC(Cl)C(C)(O)CC2C2OC21C")).unwrap();
            parse_smiles(black_box("COC1=CC(=O)C(=O)C2=C1C1=C(O2)C2=C(OC1)C(OC(C)=O)=C(OC)C=C2")).unwrap();
            parse_smiles(black_box("CCCCOC1C2=C(OC(C)C1C)C1=C(OC(C)(C)C=C1)C1=C2OC(=O)C=C1C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C2=CC(=O)C3=C(O)C=C(OC4OC(COC(=O)CC(=O)OCC5OC(OC6=CC(O)=C7C(=O)C=C(C8=CC=C(OC)C=C8)OC7=C6)C(O)C(O)C5O)C(O)C(O)C4O)C=C3O2)C=C1")).unwrap();
            parse_smiles(black_box("COC(=O)CCS")).unwrap();
            parse_smiles(black_box("CC1(CO)CCC2(C)CCC3(C)C4=C(CCC3(C)C2C1)C1(C)CCC(OC(=O)C2=CC=CC=C2)C(C)(C)C1CC4=O")).unwrap();
            parse_smiles(black_box("COC1=C(O)C2=C(CCN2C(=N)O)C2=C1[NH]C(C(=O)O)=C2")).unwrap();
            parse_smiles(black_box("COC1OC(C)(C2CC(C)=C(C)C(=O)O2)C2CCC3C4CC=C5C(O)C=CC(=O)C5(C)C4CCC132")).unwrap();
            parse_smiles(black_box("CCOC(=O)C1=CC(Br)=C[NH]1")).unwrap();
            parse_smiles(black_box("CCC(O)=NC1C(CO)OC(OC(C(O)CN)C(N)C(O)CO)C(O)C1O")).unwrap();
            parse_smiles(black_box("CCCC1N=C(O)C(N=C(O)C(C)N=C(O)CN=C(O)C2=CC(O)=CC=C2O)CN=C(O)C(CS(=O)(=O)O)N=C(O)C=CC2=CSC(=N2)C(CC2=CC=C(O)C=C2)NC(=O)C(=O)C(C(C)CC)N=C1O")).unwrap();
            parse_smiles(black_box("COCC(OC)C(O)C(OC)C(OC)C(=N)O")).unwrap();
            parse_smiles(black_box("CC1OC(OC2C3OCC(O3)C(O)C2O)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("C=C1CC(OC)(C(O)C(O)=NCC2CC(O)C(C)(C)C(CC(=O)CCCC=CC=CC=CC(O)=NC(CCCNC(=N)N)C(=O)O)O2)OC(C)C1C")).unwrap();
            parse_smiles(black_box("CN(CCC1=CN(C2OC(CO)C(O)C(O)C2O)C=N1)C(=O)C=CC1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("C=C1C(OC(=O)C=CC2=CC=CC=C2)CC(OC(C)=O)C2(COC(C)=O)C(OC(C)=O)C(OC(C)=O)C3(O)C4(C)OCC3(C)C(CC4=O)C(OC(C)=O)C12")).unwrap();
            parse_smiles(black_box("CCC(C(=O)C(C)C(O)C(C)C1OC(O)(CC(=O)O)C(CC)C(O)C1C)C1OC(C)(C2(O)OC(CC)(C(O)CC)CC2C)CC1C")).unwrap();
            parse_smiles(black_box("COC1=CC(C=CC(=O)OCC2OC(OC3=CC(C4=C(OC5OC(CO)C(O)C(O)C5OC(=O)C=CC5=CC=C(O)C(OC)=C5)C(=O)C5=C(O)C=C(O)C=C5O4)=CC=C3OC3OC(CO)C(O)C(O)C3O)C(O)C(O)C2O)=CC=C1O")).unwrap();
            parse_smiles(black_box("CCCCCCCCCCCCOC(=O)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("N=C(O)C1=C(C2=C[NH]C3=CC(O)=CC=C23)N=CC2=C1C1=CC(O)=CC=C1[NH]2")).unwrap();
            parse_smiles(black_box("OCC1OC(OC2C3COC2C(O)C(O)O3)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCO")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C(O)C(C)O)C=C1")).unwrap();
            parse_smiles(black_box("COC(=O)C1=CC2=CC(O)=C(O)C=C2C(C2=CC(=O)OC(C(=O)O)=C2)C1C(=O)O")).unwrap();
            parse_smiles(black_box("CCOC1=CC(C=O)=CC=C1OCC1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("O=C(C=CC1=CC=C(O)C(O)=C1)OCC1OC(OC2=CC=CC=C2C(=O)O)C(O)C(OC(=O)C=CC2=CC=C(O)C(O)=C2)C1O")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C2=C(O)C(O)C(O)(C3=CC=C(O)C(Cl)=C3)C2=O)C=C1Cl")).unwrap();
            parse_smiles(black_box("CC1CC2=CC=CC=C2C(=O)CC2=C(Cl)C(O)=CC(O)=C2C(=O)O1")).unwrap();
            parse_smiles(black_box("CCCCCCCCC=CCCCCCCCCCCCCC(O)C(O)=NC(CO)C(O)C(O)CCCC=CCCCCCCCCC")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C(=O)CC#N)C=C1")).unwrap();
            parse_smiles(black_box("COC1=CC(C)C2CC3OC(OC4OC(CO)C(O)C(O)C4O)CC4C(C)(O)C(OC(C)=O)C(O)C(C2(C)C1=O)C34C")).unwrap();
            parse_smiles(black_box("COC(=O)CCCCCCCC1C=CC(=O)C1CC=CC(C)=O")).unwrap();
            parse_smiles(black_box("CC=C(C)C(=O)OC1C(OC(=O)C(C)=CC)C2(CO)C(O)CC3(C)C(=CCC4C5(C)CCC(O)C(C)(C)C5CCC43C)C2CC1(C)C")).unwrap();
            parse_smiles(black_box("CC(CCC1(O)OC2CC3C4CCC5CC(OC6OC(CO)C(OC7OC(CO)C(O)C(OC8OCC(O)C(O)C8O)C7O)C(O)C6O)C(O)CC5(C)C4CCC3(C)C2C1C)COC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("CC(C)CCCCCC=CCCCCCCCCOC(=O)C(C)C")).unwrap();
            parse_smiles(black_box("CC1(C)CC(O)C2(C)CCC3(C)C(=CCC4C5(C)CCC(OC6OC(C(=O)O)C(O)C(O)C6OC6OCC(O)C(O)C6OC6OC(CO)C(O)C(O)C6O)C(C)(CO)C5CCC43C)C2C1")).unwrap();
            parse_smiles(black_box("CCCCCCCCCCCCCCCCCC(=O)OC12C(O)C(OC)C3(O)CC(C1C3OC(=O)C1=CC=CC=C1)C13C(OC)CCC4(COC)CN(C)C1C2C(OC)C43")).unwrap();
            parse_smiles(black_box("CC=C(C)C(=O)OC1CC(C)(C)CC2C3=CCC4C5(C)CCC(OC6OC(C(=O)O)C(O)C(OC7OC(CO)C(O)C(O)C7O)C6O)C(C)(C)C5CCC4(C)C3(C)CC(O)C12CO")).unwrap();
            parse_smiles(black_box("CCC(C)CC(C)CCCCCCCCCCC(O)C(C)N")).unwrap();
            parse_smiles(black_box("CC(=O)OC1CCC2(C)C3CCC(=CC(=O)N(C)CCO)C(C)C3C(O)CC2C1(C)C")).unwrap();
            parse_smiles(black_box("OCC1C(C2=CC=C(O)C(O)=C2)OC(C2=CC=C(O)C(O)=C2)C1CO")).unwrap();
            parse_smiles(black_box("CCC(C)C=C(C)C=CC(O)C(C)(O)C(O)=NCC(O)=NC(C(O)=NC(C(O)=NC(C(O)=NC1C(O)=NC(COC)C(O)=NCC(O)=NC(C)C(O)=NC(C(C)O)C(O)=NC(C(OC)C2=CC=C(O)C=C2)C(=O)N2CCCCC2C(=O)OC1C(C)C)C(C)C(C)C(=N)O)C(C)N)C(C)O")).unwrap();
            parse_smiles(black_box("CC1=CC(O)N(CCCCN=C(O)C=CC2=CC=CC=C2)C1=O")).unwrap();
            parse_smiles(black_box("NC(CCC1=CN2C=CCCC2C2CCC(C(=O)O)N=C12)C(=O)O")).unwrap();
            parse_smiles(black_box("O=C1C=C2C=CCCC23C2=CC4=C(C=C2CCN13)OCO4")).unwrap();
            parse_smiles(black_box("C=CC(C)(CCC=C(C)CCC=C(C)CCC=C(C)COC1OC(COC2OC(C)C(O)C(O)C2O)C(OC2OC(C)C(O)C(O)C2OC(=O)C(C)=CCCC(C)=CCCC(C)=CCCC(C)(C=C)OC2OC(CO)C(O)C(O)C2OC2OC(CO)C(O)C(O)C2O)C(O)C1O)OC1OC(CO)C(O)C(O)C1OC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC(=O)C1=CC(C(=O)O)=CO1")).unwrap();
            parse_smiles(black_box("COC1CC(OC2CCC3(C)C(=CCC4(O)C3CC(OC(=O)C=C(C)C(C)C)C3(C)C(C(C)=O)CCC43O)C2)OC(C)C1OC1CC(OC)C(OC2CC(OC)C(OC3OC(CO)C(O)C(O)C3O)C(C)O2)C(C)O1")).unwrap();
            parse_smiles(black_box("OC(=NCCC1=CC=CC=C1)C(O)C(O)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("C=C1CCC(O)C(C)(CCC=C(C)C)C1CCC1(C)CCCC2=C1CCC(C)C2(C)CC1=CC(OS(=O)(=O)O)=CC=C1O")).unwrap();
            parse_smiles(black_box("CS(=O)(=O)C=CC#CC=C1C=CC2(CCCO2)O1")).unwrap();
            parse_smiles(black_box("CC1=CC(=O)C2=C(O)C3=C(C=C2O1)OC(C)(C)C(OC(=O)C(C)(O)C(C)O)C3")).unwrap();
            parse_smiles(black_box("O=C(C=CC1=CC=C(O)C=C1)OC1OC2COC(=O)C3=CC(O)=C(O)C(O)=C3C3=C(C=C(O)C(O)=C3O)C(=O)OC2C(OC(=O)C2=CC(O)=C(O)C(O)=C2)C1O")).unwrap();
            parse_smiles(black_box("C=CC(C)(C)C1=CC(O)=C(O)C2=C1OC1=C(O)C3=C(C=C1C2=O)C(O)C(C(C)(C)O)O3")).unwrap();
            parse_smiles(black_box("COC1=CC(C)=C2C(=O)C=C(CC(C)=O)OC2=C1C1OC(CO)C(O)C(O)C1OC(=O)C=CC1=CC=C(O)C(O)=C1")).unwrap();
            parse_smiles(black_box("CC(=O)C1=C(C)C=C2C(=C1O)C(=O)C1=C(OC3OC(CO)C(O)C(O)C3O)C=C(O)C=C1C2C1C2=CC(O)=CC(OC3OC(CO)C(O)C(O)C3O)=C2C(=O)C2=C1C=C(C)C(C(C)=O)=C2O")).unwrap();
            parse_smiles(black_box("CC=C(C)C(=O)OC1C=C(C)CC(OC(=O)C2=CC=C(OC)C=C2)C2C(C(C)(C)O)CCC12C")).unwrap();
        }));
}

// #[bench]
// fn parse_100_smiles(b: &mut test::Bencher) {
fn parse_100_smiles(c: &mut Criterion) {
    c.bench_function(
        "parse_smiles 100",
        |b| b.iter_with_large_drop(|| {
            parse_smiles(black_box("COC(=O)CC(O)CCC(O)CCCC(C)C")).unwrap();
            parse_smiles(black_box("COC1C(C)OC(OC(C)=O)CC1(C)[N+](=O)[O-]")).unwrap();
            parse_smiles(black_box("CC(O)=NC1=CC(=O)C2(CO)OC2C1O")).unwrap();
            parse_smiles(black_box("CC1=CC2=C(C=C1O)C1=C(C3=CC4=C(C=C(O)C5=CC(C)=C(O)C=C45)C(C)(C)C3=O)C(=O)C(C)(C)C1=CC2=O")).unwrap();
            parse_smiles(black_box("N#CC(N)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("COC1=CC(C=O)=CC=C1OC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC1=CC(CC2=CC=C3C=C(OC)C(O)=NC3=C2OS(=O)(=O)O)=CC=C1O")).unwrap();
            parse_smiles(black_box("C=C1CC(C)CC2CC=CC(CC=CC(=O)OC(C=CC3CC(C)=CCO3)C(O)CC3OC3C(O)C1)O2")).unwrap();
            parse_smiles(black_box("CCCCCCCCCCCCCCCC(=O)OCC1OC(OCC2=C(C)CC(C(C)C3CCC4C5CC6OC67C(O)C=CC(=O)C7(C)C5CCC34C)OC2=O)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("O=C(C=CCCCCCCC=CC1=CC=C2OCOC2=C1)N1CCCCC1")).unwrap();
            parse_smiles(black_box("C=CC1C=C(C)C2C3C(OC4=CC=C(C=C4)CC4(OCC)CC(C(=O)C13)C(O)=N4)C1C(C)CC(C)CC21C")).unwrap();
            parse_smiles(black_box("CCC(C)C1N=C(O)C(CC(C)C)N=C(O)C2CCCN2C(=O)C2CCCN2C1=O")).unwrap();
            parse_smiles(black_box("COC1=CC(CC2CN(C(C)=O)C(CC3=CC=CC=C3)CN2C(C)=O)=CC2=C1OCO2")).unwrap();
            parse_smiles(black_box("COC1C(O)=NC2=CC=C(C=CC3(C)CCC(C(C)(C)O)O3)C(O)=C2C1(O)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("CC1=CC2OC3CC4OC(=O)C=CC=CC(=O)OCCC(C)C(O)C(=O)OCC2(CC1O)C4(C)C31CO1")).unwrap();
            parse_smiles(black_box("CCC1OC2OC3COC(C4=CC=CC=C4)OC3C(=O)C2O1")).unwrap();
            parse_smiles(black_box("COC1=CC2=C(C=C1C1=COC3=C4C=CC(C)(C)OC4=C(OC)C(OC)=C3C1=O)OCO2")).unwrap();
            parse_smiles(black_box("CC(C)CCCCCCCCCCCC1CC(O)=NC(CCC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(O)=NC(CC(C)C)C(O)=NC(CC(=O)O)C(O)=NC(CC(C)C)C(O)=NC(C(C)C)C(=O)O1")).unwrap();
            parse_smiles(black_box("CC(=NNC1=CC=C([N+](=O)[O-])C=C1[N+](=O)[O-])C1=CC=C2OCOC2=C1")).unwrap();
            parse_smiles(black_box("C1=C2CCCCCCCCCCCN3CCC=C(CCCCCCCCCCCN(CC1)C2)C3")).unwrap();
            parse_smiles(black_box("COC1=C(C)C(OC(=O)C2=C(C)C=C(O)C=C2O)=C(C)C(C)=C1C(=O)OC1=C(C)C(C)=C(C(=O)O)C(O)=C1C")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C=C2OC(=O)C=C2C2=CC(Br)=C(O)C(Br)=C2)C=C1")).unwrap();
            parse_smiles(black_box("C=C1C(=O)OC2C=C(C)C(OC(C)=O)CC=C(COC(C)=O)CC(OC(=O)C(C)=CCO)C12")).unwrap();
            parse_smiles(black_box("CC1OC(N2C3=CC=CC=C3C3=CC=C4C5=CC(Cl)=CC=C5[NH]C4=C32)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC1OCC(N=C(C)O)C(O)C1O")).unwrap();
            parse_smiles(black_box("C=CC1(C)CC(OC(=O)C(=C)CO)C(C(=C)C(=O)OC)C(O)C1C(=C)CO")).unwrap();
            parse_smiles(black_box("NCC(CCC(N)C(=O)O)OC1OC(CO)C(O)C(O)C1OC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("CCCCCCCCC=CCCCCCCCC(=O)OCC(COC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)OC(=O)CCCCCCCC=CCCCCCCCC")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C2OC3=CC(OC)=CC(O)=C3C(=O)C2OC2OC(C)C(O)C(O)C2O)C=C1")).unwrap();
            parse_smiles(black_box("CN=C(O)N=C(O)OC")).unwrap();
            parse_smiles(black_box("CC1OC(OC2CCC3(C)C4CCC5(C)C(C6=CC(=O)OC6)CCC5(O)C4CCC3(O)C2)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("O=C(O)CCCCCCCCC=CC(O)=NC1CC2(O)C(=O)C(C1O)C(O)C1OC12")).unwrap();
            parse_smiles(black_box("CONCCCCCCCCCCCCC1=CC=CN=C1")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C(=O)O)C(C(=O)N2CCC3=CC4=C(C=C3C2=O)OCO4)=C1OC")).unwrap();
            parse_smiles(black_box("CNC(C)C1C(O)CC2(C)C3CCC4C(=CC3=CCC12C)CCC(N(C)C(C)=O)C4(C)C")).unwrap();
            parse_smiles(black_box("C1=CC=C(CC2CNC(CC3=CC=CC=C3)CN2)C=C1")).unwrap();
            parse_smiles(black_box("COC1=CC(=O)C2(C3=CC(=O)C4=C(O)C=C(O)C=C4O3)C(CC=C(C)C2CC2(O)C(C)(O)CCC3C(C)(C)CC(O)CC32C)C1=O")).unwrap();
            parse_smiles(black_box("O=C1OC(C(O)CO)C(=NNC2=CC=C([N+](=O)[O-])C=C2[N+](=O)[O-])C1=NNC1=CC=C([N+](=O)[O-])C=C1[N+](=O)[O-]")).unwrap();
            parse_smiles(black_box("CCCCCC=CCC=CCCCC(O)CC(=O)OC1CN(C)C(C(OC2OC(CN)C(O)C2OS(=O)(=O)O)C2OC(N3C=CC(O)=NC3=O)C(O)C2O)C(=O)N(C)C1C(=O)O")).unwrap();
            parse_smiles(black_box("CC1CCC(OC(=O)C2=CC=CC=C2)C2(C)C(OC(=O)C3=CC=CC=C3)C(O)C3C(O)C12OC3(C)C")).unwrap();
            parse_smiles(black_box("CCC1CN2C=CC3=C4[C](C=CC=C4OC)NC3=C2CC1C(=COC)C(=O)OC")).unwrap();
            parse_smiles(black_box("C=CC1C2=CCOC(=O)C2=COC1OC1OC(COC(C)=O)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("C=CC(C)(O)CCC=C(C)C(=O)OCC1OC(OC2OC=C(C(=O)O)C3CCC(C)(O)C23)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC1=C(C=CC(=O)O)C(C(=O)O)=NC2=C1C=CO2")).unwrap();
            parse_smiles(black_box("CC=C(C)C1C(C=CC=C(C)C(O)C(C)CO)C2CC(Cl)C(C)(O)CC2C2OC21C")).unwrap();
            parse_smiles(black_box("COC1=CC(=O)C(=O)C2=C1C1=C(O2)C2=C(OC1)C(OC(C)=O)=C(OC)C=C2")).unwrap();
            parse_smiles(black_box("CCCCOC1C2=C(OC(C)C1C)C1=C(OC(C)(C)C=C1)C1=C2OC(=O)C=C1C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C2=CC(=O)C3=C(O)C=C(OC4OC(COC(=O)CC(=O)OCC5OC(OC6=CC(O)=C7C(=O)C=C(C8=CC=C(OC)C=C8)OC7=C6)C(O)C(O)C5O)C(O)C(O)C4O)C=C3O2)C=C1")).unwrap();
            parse_smiles(black_box("COC(=O)CCS")).unwrap();
            parse_smiles(black_box("CC1(CO)CCC2(C)CCC3(C)C4=C(CCC3(C)C2C1)C1(C)CCC(OC(=O)C2=CC=CC=C2)C(C)(C)C1CC4=O")).unwrap();
            parse_smiles(black_box("COC1=C(O)C2=C(CCN2C(=N)O)C2=C1[NH]C(C(=O)O)=C2")).unwrap();
            parse_smiles(black_box("COC1OC(C)(C2CC(C)=C(C)C(=O)O2)C2CCC3C4CC=C5C(O)C=CC(=O)C5(C)C4CCC132")).unwrap();
            parse_smiles(black_box("CCOC(=O)C1=CC(Br)=C[NH]1")).unwrap();
            parse_smiles(black_box("CCC(O)=NC1C(CO)OC(OC(C(O)CN)C(N)C(O)CO)C(O)C1O")).unwrap();
            parse_smiles(black_box("CCCC1N=C(O)C(N=C(O)C(C)N=C(O)CN=C(O)C2=CC(O)=CC=C2O)CN=C(O)C(CS(=O)(=O)O)N=C(O)C=CC2=CSC(=N2)C(CC2=CC=C(O)C=C2)NC(=O)C(=O)C(C(C)CC)N=C1O")).unwrap();
            parse_smiles(black_box("COCC(OC)C(O)C(OC)C(OC)C(=N)O")).unwrap();
            parse_smiles(black_box("CC1OC(OC2C3OCC(O3)C(O)C2O)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("C=C1CC(OC)(C(O)C(O)=NCC2CC(O)C(C)(C)C(CC(=O)CCCC=CC=CC=CC(O)=NC(CCCNC(=N)N)C(=O)O)O2)OC(C)C1C")).unwrap();
            parse_smiles(black_box("CN(CCC1=CN(C2OC(CO)C(O)C(O)C2O)C=N1)C(=O)C=CC1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("C=C1C(OC(=O)C=CC2=CC=CC=C2)CC(OC(C)=O)C2(COC(C)=O)C(OC(C)=O)C(OC(C)=O)C3(O)C4(C)OCC3(C)C(CC4=O)C(OC(C)=O)C12")).unwrap();
            parse_smiles(black_box("CCC(C(=O)C(C)C(O)C(C)C1OC(O)(CC(=O)O)C(CC)C(O)C1C)C1OC(C)(C2(O)OC(CC)(C(O)CC)CC2C)CC1C")).unwrap();
            parse_smiles(black_box("COC1=CC(C=CC(=O)OCC2OC(OC3=CC(C4=C(OC5OC(CO)C(O)C(O)C5OC(=O)C=CC5=CC=C(O)C(OC)=C5)C(=O)C5=C(O)C=C(O)C=C5O4)=CC=C3OC3OC(CO)C(O)C(O)C3O)C(O)C(O)C2O)=CC=C1O")).unwrap();
            parse_smiles(black_box("CCCCCCCCCCCCOC(=O)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("N=C(O)C1=C(C2=C[NH]C3=CC(O)=CC=C23)N=CC2=C1C1=CC(O)=CC=C1[NH]2")).unwrap();
            parse_smiles(black_box("OCC1OC(OC2C3COC2C(O)C(O)O3)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCO")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C(O)C(C)O)C=C1")).unwrap();
            parse_smiles(black_box("COC(=O)C1=CC2=CC(O)=C(O)C=C2C(C2=CC(=O)OC(C(=O)O)=C2)C1C(=O)O")).unwrap();
            parse_smiles(black_box("CCOC1=CC(C=O)=CC=C1OCC1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("O=C(C=CC1=CC=C(O)C(O)=C1)OCC1OC(OC2=CC=CC=C2C(=O)O)C(O)C(OC(=O)C=CC2=CC=C(O)C(O)=C2)C1O")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C2=C(O)C(O)C(O)(C3=CC=C(O)C(Cl)=C3)C2=O)C=C1Cl")).unwrap();
            parse_smiles(black_box("CC1CC2=CC=CC=C2C(=O)CC2=C(Cl)C(O)=CC(O)=C2C(=O)O1")).unwrap();
            parse_smiles(black_box("CCCCCCCCC=CCCCCCCCCCCCCC(O)C(O)=NC(CO)C(O)C(O)CCCC=CCCCCCCCCC")).unwrap();
            parse_smiles(black_box("COC1=CC=C(C(=O)CC#N)C=C1")).unwrap();
            parse_smiles(black_box("COC1=CC(C)C2CC3OC(OC4OC(CO)C(O)C(O)C4O)CC4C(C)(O)C(OC(C)=O)C(O)C(C2(C)C1=O)C34C")).unwrap();
            parse_smiles(black_box("COC(=O)CCCCCCCC1C=CC(=O)C1CC=CC(C)=O")).unwrap();
            parse_smiles(black_box("CC=C(C)C(=O)OC1C(OC(=O)C(C)=CC)C2(CO)C(O)CC3(C)C(=CCC4C5(C)CCC(O)C(C)(C)C5CCC43C)C2CC1(C)C")).unwrap();
            parse_smiles(black_box("CC(CCC1(O)OC2CC3C4CCC5CC(OC6OC(CO)C(OC7OC(CO)C(O)C(OC8OCC(O)C(O)C8O)C7O)C(O)C6O)C(O)CC5(C)C4CCC3(C)C2C1C)COC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("CC(C)CCCCCC=CCCCCCCCCOC(=O)C(C)C")).unwrap();
            parse_smiles(black_box("CC1(C)CC(O)C2(C)CCC3(C)C(=CCC4C5(C)CCC(OC6OC(C(=O)O)C(O)C(O)C6OC6OCC(O)C(O)C6OC6OC(CO)C(O)C(O)C6O)C(C)(CO)C5CCC43C)C2C1")).unwrap();
            parse_smiles(black_box("CCCCCCCCCCCCCCCCCC(=O)OC12C(O)C(OC)C3(O)CC(C1C3OC(=O)C1=CC=CC=C1)C13C(OC)CCC4(COC)CN(C)C1C2C(OC)C43")).unwrap();
            parse_smiles(black_box("CC=C(C)C(=O)OC1CC(C)(C)CC2C3=CCC4C5(C)CCC(OC6OC(C(=O)O)C(O)C(OC7OC(CO)C(O)C(O)C7O)C6O)C(C)(C)C5CCC4(C)C3(C)CC(O)C12CO")).unwrap();
            parse_smiles(black_box("CCC(C)CC(C)CCCCCCCCCCC(O)C(C)N")).unwrap();
            parse_smiles(black_box("CC(=O)OC1CCC2(C)C3CCC(=CC(=O)N(C)CCO)C(C)C3C(O)CC2C1(C)C")).unwrap();
            parse_smiles(black_box("OCC1C(C2=CC=C(O)C(O)=C2)OC(C2=CC=C(O)C(O)=C2)C1CO")).unwrap();
            parse_smiles(black_box("CCC(C)C=C(C)C=CC(O)C(C)(O)C(O)=NCC(O)=NC(C(O)=NC(C(O)=NC(C(O)=NC1C(O)=NC(COC)C(O)=NCC(O)=NC(C)C(O)=NC(C(C)O)C(O)=NC(C(OC)C2=CC=C(O)C=C2)C(=O)N2CCCCC2C(=O)OC1C(C)C)C(C)C(C)C(=N)O)C(C)N)C(C)O")).unwrap();
            parse_smiles(black_box("CC1=CC(O)N(CCCCN=C(O)C=CC2=CC=CC=C2)C1=O")).unwrap();
            parse_smiles(black_box("NC(CCC1=CN2C=CCCC2C2CCC(C(=O)O)N=C12)C(=O)O")).unwrap();
            parse_smiles(black_box("O=C1C=C2C=CCCC23C2=CC4=C(C=C2CCN13)OCO4")).unwrap();
            parse_smiles(black_box("C=CC(C)(CCC=C(C)CCC=C(C)CCC=C(C)COC1OC(COC2OC(C)C(O)C(O)C2O)C(OC2OC(C)C(O)C(O)C2OC(=O)C(C)=CCCC(C)=CCCC(C)=CCCC(C)(C=C)OC2OC(CO)C(O)C(O)C2OC2OC(CO)C(O)C(O)C2O)C(O)C1O)OC1OC(CO)C(O)C(O)C1OC1OC(CO)C(O)C(O)C1O")).unwrap();
            parse_smiles(black_box("COC(=O)C1=CC(C(=O)O)=CO1")).unwrap();
            parse_smiles(black_box("COC1CC(OC2CCC3(C)C(=CCC4(O)C3CC(OC(=O)C=C(C)C(C)C)C3(C)C(C(C)=O)CCC43O)C2)OC(C)C1OC1CC(OC)C(OC2CC(OC)C(OC3OC(CO)C(O)C(O)C3O)C(C)O2)C(C)O1")).unwrap();
            parse_smiles(black_box("OC(=NCCC1=CC=CC=C1)C(O)C(O)C1=CC=CC=C1")).unwrap();
            parse_smiles(black_box("C=C1CCC(O)C(C)(CCC=C(C)C)C1CCC1(C)CCCC2=C1CCC(C)C2(C)CC1=CC(OS(=O)(=O)O)=CC=C1O")).unwrap();
            parse_smiles(black_box("CS(=O)(=O)C=CC#CC=C1C=CC2(CCCO2)O1")).unwrap();
            parse_smiles(black_box("CC1=CC(=O)C2=C(O)C3=C(C=C2O1)OC(C)(C)C(OC(=O)C(C)(O)C(C)O)C3")).unwrap();
            parse_smiles(black_box("O=C(C=CC1=CC=C(O)C=C1)OC1OC2COC(=O)C3=CC(O)=C(O)C(O)=C3C3=C(C=C(O)C(O)=C3O)C(=O)OC2C(OC(=O)C2=CC(O)=C(O)C(O)=C2)C1O")).unwrap();
            parse_smiles(black_box("C=CC(C)(C)C1=CC(O)=C(O)C2=C1OC1=C(O)C3=C(C=C1C2=O)C(O)C(C(C)(C)O)O3")).unwrap();
            parse_smiles(black_box("COC1=CC(C)=C2C(=O)C=C(CC(C)=O)OC2=C1C1OC(CO)C(O)C(O)C1OC(=O)C=CC1=CC=C(O)C(O)=C1")).unwrap();
            parse_smiles(black_box("CC(=O)C1=C(C)C=C2C(=C1O)C(=O)C1=C(OC3OC(CO)C(O)C(O)C3O)C=C(O)C=C1C2C1C2=CC(O)=CC(OC3OC(CO)C(O)C(O)C3O)=C2C(=O)C2=C1C=C(C)C(C(C)=O)=C2O")).unwrap();
            parse_smiles(black_box("CC=C(C)C(=O)OC1C=C(C)CC(OC(=O)C2=CC=C(OC)C=C2)C2C(C(C)(C)O)CCC12C")).unwrap();
        }));
}

fn parse_all_dnp_smiles(c: &mut Criterion) {
    c.bench_function("parse_smiles all DNP", |b| {
        b.iter(|| {
            let file = File::open("dnp.smi").unwrap();
            let reader = BufReader::new(file);
            let _mols: Vec<Mol> = reader
                .lines()
                .map(|x| parse_smiles(black_box(x.unwrap().trim())).unwrap())
                .collect();
        })
    });
}

criterion_group!(
    benches,
    parse_100_smiles_kekulized,
    parse_100_smiles,
    parse_all_dnp_smiles
);
criterion_main!(benches);
