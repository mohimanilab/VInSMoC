#[derive(Debug, Eq, PartialEq, Clone)]
pub enum MolecularNode {
    NoBracket(OrganicAtom),
    Bracketed(Atom, u8, i8),
    RingClosure(u8),
}

impl From<OrganicAtom> for MolecularNode {
    fn from(atom: OrganicAtom) -> Self {
        Self::NoBracket(atom)
    }
}

impl From<(Atom, u8, i8)> for MolecularNode {
    fn from(bracketed_atom: (Atom, u8, i8)) -> Self {
        Self::Bracketed(bracketed_atom.0, bracketed_atom.1, bracketed_atom.2)
    }
}

impl From<u8> for MolecularNode {
    fn from(ring_closure: u8) -> Self {
        Self::RingClosure(ring_closure)
    }
}

// /// Checks if a character is a valid uppercase organic atom.
// /// In the SMILES specification only the organic atomic subset can be represented without
// /// surrounding the atom in square brackets.
// /// This intentionally does does not work for the lowercase second-characters of atomic symbols.
// fn is_organic_char_upper(c: char) -> bool {
//     "BCNOPSFI".contains(c)
// }

// /// Checks if character could appear in the atomic symbol of an organic atom.
// /// In the SMILES specification only the organic atomic subset can be represented without
// /// surrounding the atom in square brackets.
// fn is_organic_char(c: char) -> bool {
//     "BCNOPSFIlr".contains(c)
// }


fn molecule_node(input: &str) -> IResult<&str, MolecularNode> {
    alt((into(organic_atom), into(bracketed_atom), into(ring_closure)))(input)
}

fn molecular_path(input: &str) -> IResult<&str, Vec<(MolecularNode, BondType)>> {
    many1(tuple((molecule_node, bond)))(input)
}

type MolecularPath = Vec<(MolecularNode, BondType)>;
fn path_or_branch(input: &str) -> IResult<&str, MolecularPath> {
    eprintln!("input path_or_branch: {}", input);
    fold_many0(
        alt((molecular_path, molecular_branch)),
        Vec::new(),
        |mut v, mut other| v.into_iter().chain(other.into_iter()).collect(),
    )(input)
}

fn molecular_branch(input: &str) -> IResult<&str, MolecularPath> {
    eprintln!("input mol branch {}", input);
    delimited(
        nom::character::complete::char('('),
        map(opt(path_or_branch), |x| x.unwrap_or(Vec::new())),
        nom::character::complete::char(')'),
    )(input)
}

