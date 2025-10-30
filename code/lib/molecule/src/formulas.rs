use super::Atom;
use crate::{parsers::parse_formula, Mol};
use hash::Hash;
use itertools::Itertools;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{error::Error, fmt, hash, str::FromStr};

#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize)]
pub struct ChemFormula {
    pub cnts: FxHashMap<Atom, u32>,
}

impl ChemFormula {
    pub fn std_mass(&self) -> f64 {
        self.cnts
            .iter()
            .map(|(atom, count)| atom.std_mass() * (*count as f64))
            .sum()
    }

    pub fn exact_mass(&self) -> f64 {
        self.cnts
            .iter()
            .map(|(atom, count)| atom.exact_mass() * (*count as f64))
            .sum()
    }
}

#[derive(Debug)]
pub struct ParseFormulaError {
    inner_cause: String,
}

impl fmt::Display for ParseFormulaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "error parsing molecular formula: {}", &self.inner_cause)
    }
}

impl Error for ParseFormulaError {}

impl<'a> FromStr for ChemFormula {
    type Err = ParseFormulaError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match parse_formula(s) {
            Ok(x) => Ok(x),
            Err(x) => Err(ParseFormulaError {
                inner_cause: nom::error::convert_error(s, x),
            }),
        }
    }
}

/// When converting ChemFormulas to strings we make sure the atoms are sorted in order of
/// increasing atomic number. Additionally, any 1-count atoms have the 1 omitted for cleanliness of
/// output. The strings for equivalent formulas should always be the same due to this sorting.
impl fmt::Display for ChemFormula {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (atom, count) in self.cnts.iter().sorted_by(|(x, _), (y, _)| x.cmp(y)) {
            if *count == 1 {
                write!(f, "{}", atom)?;
            } else {
                write!(f, "{}{}", atom, count)?;
            }
        }
        Ok(())
    }
}

#[allow(clippy::derive_hash_xor_eq)]
impl Hash for ChemFormula {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        self.to_string().hash(state);
    }
}

impl From<&Mol> for ChemFormula {
    fn from(mol: &Mol) -> Self {
        let mut cnts = FxHashMap::<Atom, u32>::default();
        for atom in mol.graph.node_indices().map(|x| mol.graph[x]) {
            let curr_cnt = cnts.entry(atom).or_insert_with(|| 0);
            *curr_cnt += 1;
        }

        let hcount = cnts.entry(Atom::H).or_insert_with(|| 0);
        *hcount += mol.hydrogens.values().map(|x| *x as u32).sum::<u32>();

        Self { cnts }
    }
}

impl From<Atom> for ChemFormula {
    fn from(atom: Atom) -> Self {
        let mut cnts = FxHashMap::default();
        cnts.insert(atom, 1);
        Self { cnts }
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::hash_map::DefaultHasher, hash::Hasher};

    use super::*;
    use crate::parsers::parse_smiles;
    use approx::assert_relative_eq;

    macro_rules! form {
        ($($key:ident , $value:literal);+) => {
            {
                let mut map = FxHashMap::<Atom, u32>::default();
                $(
                    map.insert(Atom::$key, $value);
                 )*

                ChemFormula { cnts: map }
            }
        }
    }

    #[test]
    fn display() {
        assert_eq!("C".parse::<ChemFormula>().unwrap().to_string(), "C");
        assert_eq!("COH".parse::<ChemFormula>().unwrap().to_string(), "HCO");
        assert_eq!(
            "C5H8O4".parse::<ChemFormula>().unwrap().to_string(),
            "H8C5O4"
        );
        assert_eq!(
            "C5O4H8".parse::<ChemFormula>().unwrap().to_string(),
            "H8C5O4"
        );
        assert_eq!(
            "O4H8C5".parse::<ChemFormula>().unwrap().to_string(),
            "H8C5O4"
        );
        assert_eq!(
            "C10H6N2OH".parse::<ChemFormula>().unwrap().to_string(),
            "H7C10N2O"
        );
    }

    // taken from rust docs
    fn calculate_hash(form: &ChemFormula) -> u64 {
        let mut s = DefaultHasher::new();
        form.hash(&mut s);
        s.finish()
    }

    #[test]
    fn hash_eq_agreement() {
        // C14H28O5S
        let mols = vec![
            "CCCCCCCCCCCCOC(=O)CS(=O)(=O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C14H28O5S"
        let mols = vec![
            "CCCCCCCCCCCCOC(=O)CS(=O)(=O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
            "CCCCCCCCSC1C(C(C(C(O1)CO)O)O)O",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C20H20N3O"
        let mols = vec![
            "CC[N+](=C1C=CC2=NC3=C(C=C(C4=CC=CC=C43)N)OC2=C1)CC",
            "CC1=[N+](N(C(=O)C1N=CC=CC2=CC=CC=C2)C3=CC=CC=C3)C",
            "CC1=C(C=C(C=C1)C2=[NH+]NC(=C2)C=CC(=O)C3=CC=C(C=C3)N)C",
            "CC1=CC(=CC=C1)C[N-]C2CCN(C2=O)CC3=CC(=CC=C3)C#N",
            "CC1=CC=C(C=C1)C[N-]C2CCN(C2=O)CC3=CC(=CC=C3)C#N",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C15H15ClN2O2"
        let mols = vec![
            "CN(C)C(=O)NC1=CC=C(C=C1)OC2=CC=C(C=C2)Cl",
            "CN1CCN(CC1)C2=C(C(=O)C3=CC=CC=C3C2=O)Cl",
            "CCOC1=CC=CC=C1NC(=O)NC2=CC=C(C=C2)Cl",
            "CC(=O)NC(C1CC1)C2=CC(=C3C=CC=NC3=C2O)Cl",
            "C1CCC(C1)NC(=O)C2=NOC(=C2)C3=CC(=CC=C3)Cl",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C13H10N2O"
        let mols = vec![
            "C1=CC=C2C(=C1)N=C(O2)C3=CC=C(C=C3)N",
            "C1=CC=C2C(=C1)NC(=N2)C3=CC(=CC=C3)O",
            "C1=CC=C2C(=C1)C3=C(C=CC(=C3)N)NC2=O",
            "COC1=CC2=NC3=CC=CC=C3N=C2C=C1",
            "C1=CC=C2C(=C1)C(=O)NC3=CC=CC=C3N2",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C4H11NO"
        let mols = vec!["CC(C)(CO)N", "CN(C)CCO", "CCN(CC)O", "CCNCCO", "COCCCN"]
            .into_iter()
            .map(|x| (&parse_smiles(x).unwrap()).into())
            .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C21H27F4NO"
        let mols = vec![
            "C1CC(CC2(C1)CCN(CC2)CCCC(=O)C3=CC=C(C=C3)F)C(F)(F)F",
            "CCC(CC(=C(C1=C(C=C(C=C1F)OCCF)F)NC)C(=CC)C=CC)F",
            "CC1=NC2=C(C=C1)C(=C(C=C2)OF)C(F)(F)F.CC(C)(C)C1CCCCC1",
            "CC1=CC(=CN(C1)CC(CC(C)(C)C2=C(C=CC(=C2)F)C)(C(F)(F)F)O)C",
            "CCCCCC1CCC(CC1)CCC(OC2=CC(=C(C(=C2)F)C#N)F)(F)F",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C4H6O3"
        let mols = vec![
            "CC(=O)OC(=O)C",
            "CC1COC(=O)O1",
            "CCC(=O)C(=O)O",
            "CC(=O)C(=O)OC",
            "C(CC(=O)O)C=O",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C18H17ClN2O2"
        let mols = vec![
            "CC1=NN=C(C2=CC(=C(C=C2C1)OC)OC)C3=CC(=CC=C3)Cl",
            "CC1CN2CC(=O)NC3=C(C2(O1)C4=CC=CC=C4)C=C(C=C3)Cl",
            "CC1CN2CC(=O)NC3=C(C2(O1)C4=CC=CC=C4)C=C(C=C3)Cl",
            "C1CNCC1C(=O)C2=CC=C(C=C2)NC(=O)C3=CC=C(C=C3)Cl",
            "C1C2=C(C=CC(=C2NC1=O)Cl)OCCN3CC4=CC=CC=C4C3",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C15H26N2"
        let mols = vec![
            "C1CCN2CC3CC(C2C1)CN4C3CCCC4",
            "C1CCN2CC3CC(C2C1)CN4C3CCCC4",
            "C1CCN2CC3CC(C2C1)CN4C3CCCC4",
            "C1CCN2CC3CC(C2C1)CN4C3CCCC4",
            "C1CCN2CC3CC(C2C1)CN4C3CCCC4",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C5H12O2"
        let mols = vec![
            "C(CCO)CCO",
            "CC(C)(CO)CO",
            "CCCOCCO",
            "CC(C)OCCO",
            "CCOCOCC",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C13H10BrN"
        let mols = vec![
            "C1=CC=C(C=C1)C=NC2=CC=C(C=C2)Br",
            "C1C2=C(C=CC(=C2)N)C3=C1C=C(C=C3)Br",
            "C1=CC=C(C=C1)N=CC2=CC=C(C=C2)Br",
            "C1=CC=C(C=C1)C=NC2=CC(=CC=C2)Br",
            "CN1C2=C(C=C(C=C2)Br)C3=CC=CC=C31",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C10H15N4O9P"
        let mols = vec![
            "C1=NC(=C(N1C2C(C(C(O2)COP(=O)(O)O)O)O)NC=O)C(=O)N",
            "C1=NC(=C(N1C2C(C(C(O2)COP(=O)(O)O)O)O)NC=O)C(=O)N",
            "C1=NC(=C(N1C2C(C(C(O2)COP(=O)(O)O)O)O)NC=O)C(=O)N",
            "C1=NC(=C(N1C2C(C(C(O2)COP(=O)(O)O)O)O)NC=O)C(=O)N",
            "C1=C(NC(=N1)C2C(C(C(O2)C(NC=O)OP(=O)(O)O)O)O)C(=O)N",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C20H24N4"
        let mols = vec![
            "CN(C)C1=CC=C(C=C1)CCNCC2=CC3=C(C=C2)C=CC(=N3)N",
            "CC1=C2C=C(C=CC2=NN1)NC3CCCN(C3)CC4=CC=CC=C4",
            "C1CC(=CCC1N2CCN(CC2)C3=CC=CC=N3)C4=CC=CC=N4",
            "CC1=C(C=C(C=C1)N2CCN(CC2)CC3=NC4=CC=CC=C4N3)C",
            "CC1=CC2=C(C=C1C)N(C(=N2)N3CCNCC3)CC4=CC=CC=C4",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C17H24Cl3N"
        let mols = vec![
            "CCC1CCN2CCCCC2C1C3=CC(=C(C=C3)Cl)Cl.Cl",
            "CCC1CCN2CCCCC2C1C3=CC(=C(C=C3)Cl)Cl.Cl",
            "C[NH+](C)CC1C2CCC(C1C3=CC(=C(C=C3)Cl)Cl)CC2.[Cl-]",
            "CN(C)CC1C2CCC(C1C3=CC(=C(C=C3)Cl)Cl)CC2.Cl",
            "C[NH+](C)CC1C2CCC(C1C3=CC(=C(C=C3)Cl)Cl)CC2.[Cl-]",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C7H7N5O"
        let mols = vec![
            "CC1=CN=C2C(=N1)C(=O)NC(=N2)N",
            "CC(=O)NC1=NC=NC2=C1NC=N2",
            "C1=COC(=C1)C2=NC(=NC(=N2)N)N",
            "CC1=CN=C2C(=O)NC(=NC2=N1)N",
            "C1=CC(=CC=C1C(=O)NN)N=[N+]=[N-]",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C13H7FO"
        let mols = vec![
            "C1=CC=C2C(=C1)C3=C(C2=O)C=C(C=C3)F",
            "C1=CC=C2C(=C1)C3=C(C2=O)C=CC=C3F",
            "C1=CC=C2C(=C1)C3=C(C2=O)C=CC(=C3)F",
            "C1=CC=C2C(=C1)C3=C(C2=O)C(=CC=C3)F",
            "C1=CC=C2C(=C1)C=C3C2=CC=C(C3=O)F",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C6H13NS2"
        let mols = vec![
            "CC1NC(SC(S1)C)C",
            "CCN(CC)C(=S)SC",
            "CC1NC(SC(S1)C)C",
            "CN(C)C1CSCSC1",
            "CN(C)C1SCCCS1",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C5H7N"
        let mols = vec![
            "CCC=CC#N",
            "CC=C(C)C#N",
            "CC(C=C)C#N",
            "CCC=CC#N",
            "CC1=CNC=C1",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C17H24ClN3O"
        let mols = vec![
            "CN(C)CCN(CC1=CC=C(C=C1)OC)C2=CC=CC=N2.Cl",
            "CCN(CC)CCCNCC1=CC(=C2C=CC=NC2=C1O)Cl",
            "CC1CN(CCN1)CC(=O)N2CC(C3=C2C=C(C=C3)Cl)(C)C",
            "C1CCC(CC1)NC(=O)N2CCN(CC2)C3=CC=C(C=C3)Cl",
            "C1CCN(CC1)C2CCN(CC2)C(=O)C3=CC(=C(C=C3)N)Cl",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C4HCl2N3O2"
        let mols = vec![
            "C1=NC(=C(C(=N1)Cl)[N+](=O)[O-])Cl",
            "C1=C(C(=NC(=N1)Cl)Cl)[N+](=O)[O-]",
            "C1=C(N=C(N=C1Cl)[N+](=O)[O-])Cl",
            "C1=C(N=C(N=C1Cl)Cl)[N+](=O)[O-]",
            "C1(=NNN=C1C(=O)Cl)C(=O)Cl",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C21H27FN2O2"
        let mols = vec![
            "COC1=CC=CC=C1N2CCN(CC2)CCCC(C3=CC=C(C=C3)F)O",
            "CC1=C(C2=CC=CC=C2N1CCF)C(=O)OC3CC4CCCC(C3)N4C",
            "COC1=CC=CC=C1N2CCN(CC2)CCCCOC3=CC=C(C=C3)F",
            "C1C2CC3(CC1CC(C2)(C3)F)CNCC4=CC=C(C=C4)C=CC(=O)NO",
            "CC(C)(C)C1=CC=C(C=C1)CNC(=O)CC2=CC(=C(C=C2)OCCN)F",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C6H7NO2"
        let mols = vec![
            "CCOC(=O)C(=C)C#N",
            "CCN1C(=O)C=CC1=O",
            "C=CC(=O)OCCC#N",
            "C1=CC(=C[N+](=C1)[O-])CO",
            "CC(=O)OC(C=C)C#N",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C23H29NO3"
        let mols = vec![
            "CCCOC(C1=CC=CC=C1)(C2=CC=CC=C2)C(=O)OC3CCN(CC3)C",
            "CCC(C1=CC=CC=C1)C(=O)OCCN2CCOC(C2C)C3=CC=CC=C3",
            "CCOC(=O)C1(CCN(CC1)CCOCC2=CC=CC=C2)C3=CC=CC=C3",
            "CCOC(=O)C1(CCN(CC1)CCC(C2=CC=CC=C2)O)C3=CC=CC=C3",
            "CC1(CCCC2(C1CCC3(C2(CC4=CNC5=C4C3=CC(=O)C5=O)O)C)C)C",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C11H17NO2"
        let mols = vec![
            "CC(CC1=C(C=CC(=C1)OC)OC)N",
            "CC1=CC=C(C=C1)N(CCO)CCO",
            "CNCCC1=CC(=C(C=C1)OC)OC",
            "CC(CC1=CC(=C(C=C1)OC)OC)N",
            "CC1=CC(=CC=C1)N(CCO)CCO",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C18H18ClNS"
        let mols = vec![
            "CN(C)CCC=C1C2=CC=CC=C2SC3=C1C=C(C=C3)Cl",
            "CN(C)CCC=C1C2=CC=CC=C2SC3=C1C=C(C=C3)Cl",
            "CN(C)CCC=C1C2=CC=CC=C2SC3=C1C=C(C=C3)Cl",
            "CC1=NC(=CC=C1)C#CC(C(C)C)SC2=CC(=CC=C2)Cl",
            "CCC1(C2=C(C=CC(=C2)C3=CC(=CC=C3)Cl)NC1=S)CC",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C9H15N5O4S"
        let mols = vec![
            "C1CCN(CC1)C2=NC(=[N+](C(=C2)N)OS(=O)(=O)[O-])N",
            "C1CCN(CC1)C2=NC(=N)N(C(=C2)N)OS(=O)(=O)O",
            "CCS(=O)(=O)N1C(=NC(=N1)C(=O)N2CCOCC2)N",
            "CC(=O)C1CN(CCN1[N+]2=CC(=N)O[N-]2)S(=O)(=O)C",
            "CCS(=O)(=O)N1CCN(C(C1)C=O)[N+]2=CC(=N)O[N-]2",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C11H20ClN5"
        let mols = vec![
            "CCN(CC)C1=NC(=NC(=N1)Cl)N(CC)CC",
            "CC(C)(C)NC1=NC(=NC(=N1)Cl)NC(C)(C)C",
            "CCC(C)NC1=NC(=NC(=N1)Cl)NC(C)CC",
            "CCCCNC1=NC(=NC(=N1)Cl)NCCCC",
            "CCN(CC)C1=NC(=NC(=N1)NC(C)(C)C)Cl",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C4H8N2O2"
        let mols = vec![
            "C1COCCN1N=O",
            "CC(=NO)C(=NO)C",
            "C1CN(C(=O)C1N)O",
            "CC(=O)NNC(=O)C",
            "C(CC(=O)N)C(=O)N",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C8H7NO"
        let mols = vec![
            "C1C2=CC=CC=C2NC1=O",
            "C1=CC(=CC=C1CC#N)O",
            "CC1=NC2=CC=CC=C2O1",
            "C1=CC(=CC2=C1C=CN2)O",
            "C1=CC2=C(C=CN2)C(=C1)O",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
        // "C35H58O5"
        let mols = vec![
            "CC(C)C1CCC2(C1C3CCC4C5(CCC(C(C5CCC4(C3(CC2)C)C)(C)C)OC(=O)CCCC(=O)O)C)CO",
            "CCCCC=CCCCCCCCC(=O)OC(CO)COC(=O)CCCCC=CCC=CCC=CCC=CCC",
            "CCCCC=CCCCCCCCC(=O)OCC(CO)OC(=O)CCCCC=CCC=CCC=CCC=CCC",
            "CC(=C)C1CCC2(C1C3CCC4C5(CCC(C(C5CCC4(C3(CC2)C)C)(C)C)O)C)C(=O)OCCOCCOC",
            "CCCCCCCCC=CCCCC(=O)OCC(COC(=O)CCCCC=CCC=CCC=CCC=CCC)O",
        ]
        .into_iter()
        .map(|x| (&parse_smiles(x).unwrap()).into())
        .collect::<Vec<ChemFormula>>();
        for (i, j) in mols
            .iter()
            .tuple_combinations::<(&ChemFormula, &ChemFormula)>()
        {
            assert_eq!(i, j);
            assert_eq!(i.to_string(), j.to_string());
            assert_eq!(calculate_hash(i), calculate_hash(j));
        }
    }

    #[test]
    fn chemformula_fromstr_works() {
        assert!("abcdefghi".parse::<ChemFormula>().is_err());
        assert_eq!("C".parse::<ChemFormula>().unwrap(), form! {C , 1});
        assert_eq!("COH".parse::<ChemFormula>().unwrap(), form! {C,1; O,1; H,1});
        assert_eq!(
            "C5H8O4".parse::<ChemFormula>().unwrap(),
            form! {C,5; O,4; H,8}
        );
        assert_eq!(
            "C10H6N2OH".parse::<ChemFormula>().unwrap(),
            form! {C,10; O,1; H,7; N,2}
        );
    }

    #[test]
    fn chemformula_exact_mass_works() {
        // compares to PubChem Monoisotopic Mass
        // https://pubchem.ncbi.nlm.nih.gov/compounds/15804
        assert_relative_eq!(
            "C14H28O5S".parse::<ChemFormula>().unwrap().exact_mass(),
            308.165745,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/16939
        assert_relative_eq!(
            "C20H20N3O".parse::<ChemFormula>().unwrap().exact_mass(),
            318.160637,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/16115
        assert_relative_eq!(
            "C15H15ClN2O2".parse::<ChemFormula>().unwrap().exact_mass(),
            290.082205,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/31631
        assert_relative_eq!(
            "C13H10N2O".parse::<ChemFormula>().unwrap().exact_mass(),
            210.079313,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/1672
        assert_relative_eq!(
            "C4H11NO".parse::<ChemFormula>().unwrap().exact_mass(),
            89.084064,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/30761
        assert_relative_eq!(
            "C21H27F4NO".parse::<ChemFormula>().unwrap().exact_mass(),
            385.202877,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/296
        assert_relative_eq!(
            "C4H6O3".parse::<ChemFormula>().unwrap().exact_mass(),
            102.031694,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/4617
        assert_relative_eq!(
            "C18H17ClN2O2".parse::<ChemFormula>().unwrap().exact_mass(),
            328.097855,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/3966
        assert_relative_eq!(
            "C15H26N2".parse::<ChemFormula>().unwrap().exact_mass(),
            234.209599,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/18053
        assert_relative_eq!(
            "C5H12O2".parse::<ChemFormula>().unwrap().exact_mass(),
            104.08373,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/23126
        assert_relative_eq!(
            "C13H10BrN".parse::<ChemFormula>().unwrap().exact_mass(),
            258.99966,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/1011
        assert_relative_eq!(
            "C10H15N4O9P".parse::<ChemFormula>().unwrap().exact_mass(),
            366.057665,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/18894
        assert_relative_eq!(
            "C20H24N4".parse::<ChemFormula>().unwrap().exact_mass(),
            320.200097,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/32539
        assert_relative_eq!(
            "C17H24Cl3N".parse::<ChemFormula>().unwrap().exact_mass(),
            347.097433,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/20813
        assert_relative_eq!(
            "C7H7N5O".parse::<ChemFormula>().unwrap().exact_mass(),
            177.06506,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/15193
        assert_relative_eq!(
            "C13H7FO".parse::<ChemFormula>().unwrap().exact_mass(),
            198.048093,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/30079
        assert_relative_eq!(
            "C6H13NS2".parse::<ChemFormula>().unwrap().exact_mass(),
            163.048942,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/7304
        assert_relative_eq!(
            "C5H7N".parse::<ChemFormula>().unwrap().exact_mass(),
            81.057849,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/18029
        assert_relative_eq!(
            "C17H24ClN3O".parse::<ChemFormula>().unwrap().exact_mass(),
            321.16079,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/20312
        assert_relative_eq!(
            "C4HCl2N3O2".parse::<ChemFormula>().unwrap().exact_mass(),
            192.944582,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/13333
        assert_relative_eq!(
            "C21H27FN2O2".parse::<ChemFormula>().unwrap().exact_mass(),
            358.205656,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/25919
        assert_relative_eq!(
            "C6H7NO2".parse::<ChemFormula>().unwrap().exact_mass(),
            125.047678,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/4942
        assert_relative_eq!(
            "C23H29NO3".parse::<ChemFormula>().unwrap().exact_mass(),
            367.214744,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/7553
        assert_relative_eq!(
            "C11H17NO2".parse::<ChemFormula>().unwrap().exact_mass(),
            195.125929,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/2729
        assert_relative_eq!(
            "C18H18ClNS".parse::<ChemFormula>().unwrap().exact_mass(),
            315.084849,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/4202
        assert_relative_eq!(
            "C9H15N5O4S".parse::<ChemFormula>().unwrap().exact_mass(),
            289.084475,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/11380
        assert_relative_eq!(
            "C11H20ClN5".parse::<ChemFormula>().unwrap().exact_mass(),
            257.140723,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/1232
        assert_relative_eq!(
            "C4H8N2O2".parse::<ChemFormula>().unwrap().exact_mass(),
            116.058578,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/10199
        assert_relative_eq!(
            "C8H7NO".parse::<ChemFormula>().unwrap().exact_mass(),
            133.052764,
            max_relative = 0.000001
        );
        // https://pubchem.ncbi.nlm.nih.gov/compounds/5676
        assert_relative_eq!(
            "C35H58O5".parse::<ChemFormula>().unwrap().exact_mass(),
            558.428425,
            max_relative = 0.000001
        );
    }
}
