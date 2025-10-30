use std::fmt;

use super::{charge_sign, forced_charge, remaining_hydrogens};
use crate::atoms::{Atom, InorganicAtom, OrganicAtom};
use crate::{BondType, Mol, MolGraph, NodeIndex};
use nom::error::{context, convert_error};
use parsers::take_n_try_parse;
use regex::Regex;

use rustc_hash::{FxHashMap, FxHashSet};

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{char, digit0},
    combinator::{all_consuming, into, map, map_res, opt, peek, value},
    error::VerboseError,
    multi::{fold_many0, many0, many1, separated_list1},
    sequence::{delimited, preceded, tuple},
    IResult,
};

/// Information used to produce a fully-specified Smiles Syntax Atom Node.
///
/// The individual components are:
/// - The [`Atom`].
/// - Whether the atom is aromatic.
/// - The number of hydrogens bound to the atom.
/// - The charge on the atom.
type AtomSpecifiers = ((Atom, bool), u8, i8);

#[derive(Debug, Eq, PartialEq, Clone)]
enum SmilesSyntaxNode {
    // the actual atom, the number of hydrogens attached, its charge, and whether or not it's
    // aromatic
    Atom(Atom, Option<u8>, i8, bool),
    RingClosure(u8),
    Bond(BondType),
    BranchBegin,
    BranchEnd,
}

impl From<(OrganicAtom, bool)> for SmilesSyntaxNode {
    fn from(atom: (OrganicAtom, bool)) -> Self {
        let (org_atom, is_aromatic) = atom;
        SmilesSyntaxNode::Atom(org_atom.into(), None, 0, is_aromatic)
    }
}

impl From<AtomSpecifiers> for SmilesSyntaxNode {
    fn from(bracketed_atom: AtomSpecifiers) -> Self {
        let ((atom, is_aromatic), hydrogens, charge) = bracketed_atom;
        SmilesSyntaxNode::Atom(atom, Some(hydrogens), charge, is_aromatic)
    }
}

impl From<BondType> for SmilesSyntaxNode {
    fn from(bond_type: BondType) -> Self {
        SmilesSyntaxNode::Bond(bond_type)
    }
}

impl From<u8> for SmilesSyntaxNode {
    fn from(ring_closure: u8) -> Self {
        SmilesSyntaxNode::RingClosure(ring_closure)
    }
}

/// Parses an aliphatic organic atom (one that does not need to be surrounded by square brackets
/// in the SMILES spec) from a string.
///
/// Matches any of: B, C, N, O, P, S, F, Cl, Br, or I.
fn aliphatic_organic_atom(input: &str) -> IResult<&str, (OrganicAtom, bool), VerboseError<&str>> {
    // Note that the two-letter tags must be tried first to avoid parsing issues
    let atom_tags = (
        tag("Br"),
        tag("Cl"),
        tag("C"),
        tag("N"),
        tag("O"),
        tag("P"),
        tag("S"),
        tag("F"),
        tag("I"),
        tag("B"),
    );

    let (input, atom) = map_res(alt(atom_tags), |s: &str| s.parse::<OrganicAtom>())(input)?;

    Ok((input, (atom, false)))
}

/// Parses an aromatic organic atom from a string assuming the lower-case convention covered in:
/// https://www.chemeurope.com/en/encyclopedia/Simplified_molecular_input_line_entry_specification.html.
///
/// Matches any of: c, n, o, s.
fn aromatic_organic_atom(input: &str) -> IResult<&str, (OrganicAtom, bool), VerboseError<&str>> {
    let atom_tags = (tag("c"), tag("n"), tag("o"), tag("s"));

    let (input, atom) = map_res(alt(atom_tags), |s: &str| {
        s.to_uppercase().parse::<OrganicAtom>()
    })(input)?;

    Ok((input, (atom, true)))
}

/// Parses any organic atom by attempting the aliphatic and aromatic parsers in serial.
fn organic_atom(input: &str) -> IResult<&str, (OrganicAtom, bool), VerboseError<&str>> {
    alt((aliphatic_organic_atom, aromatic_organic_atom))(input)
}

/// Parses any inorganic atom (those that cannot be represented without [] in the
/// SMILES spec) from a string.
fn inorganic_atom(input: &str) -> IResult<&str, InorganicAtom, VerboseError<&str>> {
    alt((
        // We attempt to take two characters first to deal with two-character atoms
        // that have one-character atom substrings
        take_n_try_parse::<InorganicAtom>(2),
        take_n_try_parse::<InorganicAtom>(1),
    ))(input)
}

/// Parses a hydrogen from the current string.
fn hydrogen(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    tag("H")(input)
}

/// Parses number of hydrogens from a string.
/// Intended for use in parsing hydrogen counts from [] surrounded atoms in a SMILES string.
/// Hydrogen counts are not expected to exceed 255, and therefore are parsed into a u8.
/// If provided an empty string it is assumed that the input has no hydrogens, and a count of 0 is
/// returned.
fn hydrogen_count(input: &str) -> IResult<&str, u8, VerboseError<&str>> {
    // Find 0 or more digits preceded by a hydrogen.
    // If no digits are found, the hydrogen count is assumed to be 1.
    let (input, count) = opt(map_res(preceded(hydrogen, digit0), |s: &str| match s {
        "" => Ok(1),
        x => x.parse::<u8>(),
    }))(input)?;

    // Original parser was optional.
    // If it could not match the preceding hydrogen ('H') hydrogen count is assumed to be 0.
    Ok((input, count.unwrap_or(0u8)))
}

/// Parses a signed charge from a string.
/// Intended for use in parsing charges from [] surrounded atoms in a SMILES string.
/// Charges are never expected to exceed +/- 128 and are therefore parsed into an i8.
/// Charges strings for SMILES can be specified in a few ways.
/// First, if the charge string is an empty string the charge is assumed to be 0.
/// Next, if the charge is not empty it must begin with its sign character, either + or -.
/// After the sign character there can be a number, denoting the absolute value of the
/// charge. Alternatively, there can be multiple of the same sign character, in which case the
/// absolute charge is the number of times that sign character appears. For example, ++ and +2 both
/// denote the absolute charge 2 with a positive sign.
fn charge(input: &str) -> IResult<&str, i8, VerboseError<&str>> {
    // first find out the charge sign character
    let (input, charge_char) = opt(peek(charge_sign))(input)?;

    // if there is no such charcter, our chage must be empty and therefore we should return a 0
    // charge
    if charge_char.is_none() {
        return Ok((input, 0));
    }

    forced_charge(input)
}

// BEGIN CHIRAL GROUP
// These are to match chirality, but results are currently ignored

fn chiral_class(input: &str) -> IResult<&str, (), VerboseError<&str>> {
    let (input, _) = alt((tag("AL"), tag("SP"), tag("TB"), tag("OH")))(input)?;
    Ok((input, ()))
}

fn chiral_order(input: &str) -> IResult<&str, (), VerboseError<&str>> {
    let (input, _) = digit0(input)?;
    Ok((input, ()))
}

fn chirality(input: &str) -> IResult<&str, (), VerboseError<&str>> {
    let (input, _) = tuple((
        many1(nom::character::complete::char('@')),
        opt(tuple((chiral_class, chiral_order))),
    ))(input)?;
    Ok((input, ()))
}

// END CHIRAL GROUP

/// Parses out a smiles atom written in its most explicit, bracketed form (e.g. `[CH4]` instead of the implicit `C`).
fn bracketed_atom(input: &str) -> IResult<&str, AtomSpecifiers, VerboseError<&str>> {
    let (input, (curr_atom, _, curr_hs, curr_charge)) = delimited(
        nom::character::complete::char('['),
        tuple((
            alt((
                map(inorganic_atom, |atom| (atom.into(), false)),
                map(organic_atom, |(atom, is_aromatic)| {
                    (atom.into(), is_aromatic)
                }),
                value((Atom::H, false), hydrogen),
            )),
            opt(chirality), // matches chiral center specifiers but just ignores them for now
            hydrogen_count,
            charge,
        )),
        nom::character::complete::char(']'),
    )(input)?;
    Ok((input, (curr_atom, curr_hs, curr_charge)))
}

/// Uses `into` functions to automatically convert atom parser output to a [`SmilesSyntaxNode`].
fn atom(input: &str) -> IResult<&str, SmilesSyntaxNode, VerboseError<&str>> {
    alt((into(organic_atom), into(bracketed_atom)))(input)
}

/// Takes one character and tries to parse it into a bond type.
fn bond_type(input: &str) -> IResult<&str, BondType, VerboseError<&str>> {
    take_n_try_parse::<BondType>(1)(input)
}

/// Uses `into` functions to automatically convert bond parser output to a [`SmilesSyntaxNode`].
fn bond(input: &str) -> IResult<&str, SmilesSyntaxNode, VerboseError<&str>> {
    into(bond_type)(input)
}

/// Attempt to parse a ring closure id (in smiles, rings are specified by marking two atoms
/// joined by a back edge in the DFS traversal with the same digit, e.g. `C1CCCCC1` for cyclohexane).
fn ring_closure_digit(input: &str) -> IResult<&str, u8, VerboseError<&str>> {
    alt((
        take_n_try_parse::<u8>(1),
        preceded(tag("%"), take_n_try_parse::<u8>(2)),
    ))(input)
}

/// Uses `into` functions to automatically convert ring closure parser output to a [`SmilesSyntaxNode`].
fn ring_closure(input: &str) -> IResult<&str, SmilesSyntaxNode, VerboseError<&str>> {
    into(ring_closure_digit)(input)
}

fn ring_bond(input: &str) -> IResult<&str, Vec<SmilesSyntaxNode>, VerboseError<&str>> {
    context(
        "ring bond",
        map(tuple((opt(bond), ring_closure)), |(b, r)| match b {
            Some(x) => vec![x, r],
            None => vec![r],
        }),
    )(input)
}

/// Parse all branches coming out of the current position in the chain by repeatedly calling [`branch`].
fn branched_atom(input: &str) -> IResult<&str, Vec<SmilesSyntaxNode>, VerboseError<&str>> {
    let (input, mut output) = map(atom, |a| vec![a])(input)?;
    let (input, mut ringbonds) =
        map(many0(ring_bond), |v| v.into_iter().flatten().collect())(input)?;

    // Once we have the atom and its ring, keep trying to parse outgoing branches until impossible and append each one to an accumulator
    let (input, mut branches) = fold_many0(branch, Vec::new(), |mut acc, mut item| {
        acc.append(&mut item);
        acc
    })(input)?;

    output.append(&mut ringbonds);
    output.append(&mut branches);
    Ok((input, output))
}

fn branch(input: &str) -> IResult<&str, Vec<SmilesSyntaxNode>, VerboseError<&str>> {
    let (input, mut bond_vec) = preceded(
        tag("("),
        map(opt(bond), |one_bond| {
            one_bond.map(|b| vec![b]).unwrap_or_default()
        }),
    )(input)?;
    let (input, mut inner_chain) = atom_chain(input)?;

    // add markers for branch start and end, combine whole thing into one flat vector
    let mut full_branch = Vec::with_capacity(2 + inner_chain.len() + bond_vec.len());
    full_branch.push(SmilesSyntaxNode::BranchBegin);
    full_branch.append(&mut bond_vec);
    full_branch.append(&mut inner_chain);
    full_branch.push(SmilesSyntaxNode::BranchEnd);

    value(full_branch, tag(")"))(input)
}

fn bond_vec(input: &str) -> IResult<&str, Vec<SmilesSyntaxNode>, VerboseError<&str>> {
    map(bond, |b| vec![b])(input)
}

/// Combine the previous parsers and apply them repeatedly to parse a whole formula
fn atom_chain(input: &str) -> IResult<&str, Vec<SmilesSyntaxNode>, VerboseError<&str>> {
    map(
        many1(alt((
            branched_atom,
            map(
                tuple((alt((bond_vec, ring_bond)), branched_atom)),
                |(b, a)| b.into_iter().chain(a.into_iter()).collect(),
            ),
        ))),
        |v| v.into_iter().flatten().collect(),
    )(input)
}

/// Possible errors when constructing a molecule from a SMILES string.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub enum SmilesParseError {
    /// Input SMILES was syntactically invalid
    Parsing(String),
    /// Input SMILES open a ring closure but never closed it
    UnmatchedRingClosure,
    /// Input SMILES began with a ring closure or bond
    FirstNodeNotAtom,
    /// Input SMILES contained an organic atom with too many bonds
    HydrogenInference,
    /// Input SMILES has inconsistent bonds connecting ring closures
    InconsistentRingBonds(u8),
}

impl fmt::Display for SmilesParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Error parsing SMILES: ")?;
        match self {
            Self::Parsing(cause) => write!(f, "unable to parse input\n{}", cause),
            Self::UnmatchedRingClosure => write!(f, "unmatched ring closures"),
            Self::FirstNodeNotAtom => write!(f, "first element is not an atom"),
            Self::HydrogenInference => write!(f, "organic atom has too many bonds"),
            Self::InconsistentRingBonds(ring) => {
                write!(f, "ring bond {} has inconsistent bonds", ring)
            }
        }
    }
}

impl std::error::Error for SmilesParseError {}

pub fn parse_smiles(smiles: &str) -> Result<Mol, SmilesParseError> {
    let smiles: String = smiles.chars().filter(|x| *x != '/' && *x != '\\').collect();
    let parse_trees = match all_consuming(separated_list1(char('.'), atom_chain))(&smiles) {
        Ok((_, y)) => Ok(y),
        Err(nom::Err::Failure(e) | nom::Err::Error(e)) => {
            Err(SmilesParseError::Parsing(convert_error(smiles.as_str(), e)))
        }
        _ => unreachable!(),
    }?;

    // make sure ring closures are matching
    let mut ring_closure_cnts: FxHashMap<&u8, usize> = FxHashMap::default();
    for parse_node in parse_trees.iter().flatten() {
        let ring_closure_cnt = match parse_node {
            SmilesSyntaxNode::RingClosure(x) => ring_closure_cnts.entry(x).or_insert(0),
            _ => continue,
        };
        *ring_closure_cnt += 1;
    }
    if !ring_closure_cnts.values().all(|x| x % 2 == 0) {
        return Err(SmilesParseError::UnmatchedRingClosure);
    }

    let num_atoms = parse_trees
        .iter()
        .flatten()
        .filter(|x| matches!(x, SmilesSyntaxNode::Atom(_, _, _, _)))
        .count();
    let mut charges: Vec<i8> = Vec::with_capacity(num_atoms);

    let mut implicit_h_tracker: FxHashSet<NodeIndex> = FxHashSet::default();
    let mut arom_track: FxHashSet<NodeIndex> = FxHashSet::default();
    let mut graph = MolGraph::with_capacity(num_atoms, num_atoms);

    for parse_tree in parse_trees {
        let mut last_atom_idx;
        let mut last_node_was_atom;
        let mut last_bond_type = None;
        let mut ring_tracker: FxHashMap<u8, (NodeIndex, Option<BondType>)> = FxHashMap::default();
        let mut branch_stack: Vec<NodeIndex> = Vec::new();

        let mut it = parse_tree.into_iter();

        // special case for first item in the SMILES, it must be a node
        let first_node = it.next().unwrap();
        match first_node {
            SmilesSyntaxNode::Atom(atom, h_count_opt, charge, aromatic) => {
                last_atom_idx = graph.add_node(atom);
                charges.push(charge);
                match h_count_opt {
                    Some(x) => {
                        charges.resize(charges.len() + x as usize, 0);
                        for _ in 0..x {
                            let h_idx = graph.add_node(Atom::H);
                            graph.add_edge(last_atom_idx, h_idx, BondType::Single);
                        }
                    }
                    None => {
                        implicit_h_tracker.insert(last_atom_idx);
                    }
                };
                if aromatic {
                    arom_track.insert(last_atom_idx);
                }

                last_node_was_atom = true;
            }
            _ => return Err(SmilesParseError::FirstNodeNotAtom),
        }

        for parse_node in it {
            match parse_node {
                // this node represents an atom in the target molecule
                SmilesSyntaxNode::Atom(atom, h_count_opt, charge, aromatic) => {
                    let curr_idx = graph.add_node(atom);
                    charges.push(charge);
                    match h_count_opt {
                        Some(x) => {
                            charges.resize(charges.len() + x as usize, 0);
                            for _ in 0..x {
                                let h_idx = graph.add_node(Atom::H);
                                graph.add_edge(curr_idx, h_idx, BondType::Single);
                            }
                        }
                        None => {
                            implicit_h_tracker.insert(curr_idx);
                        }
                    };
                    if aromatic {
                        arom_track.insert(curr_idx);
                    }

                    // there was an implicit bond between previous atom and this one, add it in
                    if last_node_was_atom {
                        // if both atoms are aromatic then the implicit bond is aromatic, otherwise
                        // it's a single bond
                        let bond_type = match arom_track.contains(&last_atom_idx) && aromatic {
                            true => BondType::Aromatic,
                            false => BondType::Single,
                        };
                        graph.add_edge(last_atom_idx, curr_idx, bond_type);
                    } else if graph.node_count() > 1 {
                        // there was an explicit bond before this
                        graph.add_edge(last_atom_idx, curr_idx, last_bond_type.unwrap());
                    }

                    last_node_was_atom = true;
                    last_atom_idx = curr_idx;
                    last_bond_type = None;
                }
                SmilesSyntaxNode::Bond(bond_type) => {
                    last_node_was_atom = false;
                    last_bond_type = Some(bond_type);
                }
                SmilesSyntaxNode::RingClosure(ring_number) => {
                    match ring_tracker.remove(&ring_number) {
                        // we found the first instance of this ring already
                        Some((neighbor, bond_type)) => {
                            let bond_type = match bond_type {
                                // it had an explicit bond type
                                Some(bt) => {
                                    if let Some(last_bt) = last_bond_type {
                                        if last_bt != bt {
                                            return Err(SmilesParseError::InconsistentRingBonds(
                                                ring_number,
                                            ));
                                        }
                                    }
                                    bt
                                }
                                // it had an implicit bond type
                                None => {
                                    if let Some(bt) = last_bond_type {
                                        bt
                                    } else {
                                        match arom_track.contains(&last_atom_idx)
                                            && arom_track.contains(&neighbor)
                                        {
                                            true => BondType::Aromatic,
                                            false => BondType::Single,
                                        }
                                    }
                                }
                            };
                            graph.add_edge(neighbor, last_atom_idx, bond_type);
                        }
                        // this is the first time we found this ring
                        None => {
                            ring_tracker.insert(ring_number, (last_atom_idx, last_bond_type));
                        }
                    };
                    last_node_was_atom = true;
                    last_bond_type = None;
                }
                SmilesSyntaxNode::BranchBegin => {
                    branch_stack.push(last_atom_idx);
                    last_node_was_atom = true;
                }
                SmilesSyntaxNode::BranchEnd => {
                    last_atom_idx = branch_stack.pop().unwrap();
                    last_node_was_atom = true;
                }
            }
        }
    }

    let charges: FxHashMap<NodeIndex, i8> = (0..graph.node_count())
        .map(NodeIndex::new)
        .zip(charges.into_iter())
        .filter(|(_, charge)| *charge != 0)
        .collect();

    let hydrogens: FxHashMap<NodeIndex, u8> = FxHashMap::default();
    let mut mol = Mol {
        graph,
        hydrogens,
        charges,
        aromatic_atoms: arom_track.clone(),
    };
    for atom_idx in implicit_h_tracker {
        // FIXME: use OrganicAtom inference functions
        let mut h_count =
            match remaining_hydrogens(&mol.graph[atom_idx], mol.total_bonds(atom_idx), None) {
                Ok(x) => x,
                Err(_) => {
                    return Err(SmilesParseError::HydrogenInference);
                }
            };

        if arom_track.contains(&atom_idx) && h_count > 0 {
            h_count -= 1;
        }
        mol.hydrogens.insert(atom_idx, h_count);
    }

    Ok(mol)
}

pub fn parse_indexed_smiles(
    smiles: &str,
) -> (Result<Mol, SmilesParseError>, FxHashMap<usize, usize>) {
    let mut idx_map = FxHashMap::default();
    let smiles_regex = Regex::new(r"\[([^\[:]+):(\d+)]").unwrap();

    for (idx, mapping) in smiles_regex.find_iter(smiles).enumerate() {
        let caps = smiles_regex.captures(mapping.as_str()).unwrap();
        idx_map.insert(idx, caps[2].parse::<usize>().unwrap());
    }
    let result_smiles = smiles_regex.replace_all(smiles, "$1");
    let mol = parse_smiles(&result_smiles);
    (mol, idx_map)
}

#[cfg(test)]
mod tests;
