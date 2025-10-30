use std::{error::Error, fmt};

use crate::atoms::Atom;
use crate::{BondType, Mol, MolGraph, NodeIndex};

use nom::{
    bytes::complete::tag,
    character::complete::{alpha1, digit1, line_ending, not_line_ending, space1},
    combinator::{all_consuming, complete, map, map_opt, map_res, opt, value},
    error::{context, convert_error, VerboseError},
    multi::{count, many0, many1},
    number::complete::float,
    sequence::{delimited, preceded, separated_pair, terminated, tuple},
    Err as NomErr, IResult,
};
use petgraph::graph::node_index;
use rustc_hash::{FxHashMap, FxHashSet};

// Used to make transfer of data from line to node cleaner
struct AtomLine {
    pub atom: Atom,
    pub charge: i8,
    pub hydrogen_count: Option<u8>,
}

fn remove_line(input: &str) -> IResult<&str, &str, VerboseError<&str>> {
    terminated(not_line_ending, line_ending)(input)
}

fn parse_header_block(input: &str) -> IResult<&str, Vec<&str>, VerboseError<&str>> {
    count(remove_line, 4)(input)
}

fn parse_counts_line(input: &str) -> IResult<&str, (usize, usize), VerboseError<&str>> {
    context(
        "counts line",
        map_res(
            preceded(
                tag("M  V30 COUNTS "),
                separated_pair(digit1, space1, digit1),
            ),
            |(atom_str, bond_str): (&str, &str)| {
                atom_str
                    .parse::<usize>()
                    .and_then(|atom_cnt| Ok((atom_cnt, bond_str.parse::<usize>()?)))
            },
        ),
    )(input)
}

// Returns the (Atom type, charge, hydrogen count (if provided))
fn parse_atom_line(input: &str) -> IResult<&str, AtomLine, VerboseError<&str>> {
    // Parses the atom type
    let (input, atom) = preceded(tuple((tag("M  V30 "), digit1, space1)), alpha1)(input)?;

    // Ignores fields relating to atom coordinates and atom-atom mapping
    let (input, _) = count(tuple((space1, float)), 4)(input)?;

    // Parses optional fields for HCOUNT= and CHG= fields
    let (input, option_fields) =
        many0(preceded(space1, separated_pair(alpha1, tag("="), digit1)))(input)?;
    let mut charge = 0;
    let mut hydrogen_count = None;
    for (field, num_str) in option_fields {
        match field {
            "CHG" => charge = num_str.parse::<i8>().unwrap(),
            "HCOUNT" => hydrogen_count = Some(num_str.parse::<u8>().unwrap()),
            _ => (),
        };
    }

    let (input, _) = line_ending(input)?;

    Ok((
        input,
        AtomLine {
            atom: atom.parse::<Atom>().unwrap(),
            charge,
            hydrogen_count,
        },
    ))
}

// Returns the indices of the atom being connected as well as the BondType
fn parse_bond_line(input: &str) -> IResult<&str, (usize, usize, BondType), VerboseError<&str>> {
    // Parses for the string indicating edge type
    let (input, edge_str) =
        delimited(tuple((tag("M  V30 "), digit1, space1)), digit1, space1)(input)?;

    // Converts edge string to BondType
    let edge = match edge_str {
        "1" => BondType::Single,
        "2" => BondType::Double,
        "3" => BondType::Triple,
        "4" => BondType::Aromatic,
        _ => unreachable!(), // It is possible that a bond is ambiguous (edge code 5-8) but we assume that they won't be used.
    };

    // Parses for atom indices being connexted by bond
    let (input, (atom1, atom2)) = map_res(
        separated_pair(digit1, space1, digit1),
        |(atom1_str, atom2_str): (&str, &str)| {
            atom1_str
                .parse::<usize>()
                .and_then(|atom1| Ok((atom1, atom2_str.parse::<usize>()?)))
        },
    )(input)?;

    let (rem, _) = remove_line(input)?;

    Ok((rem, (atom1, atom2, edge)))
}

fn parse_sgroup_block(input: &str) -> IResult<&str, (), VerboseError<&str>> {
    value(
        (),
        tuple((
            tag("M  V30 BEGIN SGROUP"),
            line_ending,
            many1(remove_line),
            tag("M  V30 END SGROUP"),
        )),
    )(input)
}

fn parse_three_dim_block(input: &str) -> IResult<&str, (), VerboseError<&str>> {
    value(
        (),
        tuple((
            tag("M  V30 BEGIN OBJ3D"),
            line_ending,
            many1(remove_line),
            tag("M  V30 END OBJ3D"),
        )),
    )(input)
}

fn parse_connection_table_block(input: &str) -> IResult<&str, Mol, VerboseError<&str>> {
    let (input, _) = tuple((tag("M  V30 BEGIN CTAB"), line_ending))(input)?;
    let (input, (num_atoms, num_bonds)) = parse_counts_line(input)?;

    // Parse Atom block starting here
    let (input, _) = context(
        "start atom block",
        tuple((remove_line, tag("M  V30 BEGIN ATOM"), line_ending)),
    )(input)?;

    let (input, atom_info): (&str, Vec<AtomLine>) =
        context("atom lines", count(parse_atom_line, num_atoms))(input)?;

    // Add vertices to nodes of graph as well as their charges and hydrogen count
    let mut graph = MolGraph::with_capacity(num_atoms, num_bonds);
    let mut charges = FxHashMap::<NodeIndex, i8>::default();
    let mut hydrogens = FxHashMap::<NodeIndex, u8>::default();
    let aromatic_atoms = FxHashSet::<NodeIndex>::default();
    for single_atom in atom_info {
        let node_idx = graph.add_node(single_atom.atom);
        if single_atom.charge != 0 {
            charges.insert(node_idx, single_atom.charge);
        }
        if let Some(hcount) = single_atom.hydrogen_count {
            if hcount == 0 {
                continue;
            }
            hydrogens.insert(node_idx, hcount);
        }
    }

    let (input, _) = tuple((tag("M  V30 END ATOM"), line_ending))(input)?;

    // Parse Bond block starting here
    let input = if num_bonds > 0 {
        let (input, _) = context(
            "start bond block",
            tuple((tag("M  V30 BEGIN BOND"), line_ending)),
        )(input)?;

        let mut running_input = input;
        for _ in 0..num_bonds {
            let (input, (idx1, idx2, edge_type)) = context(
                "bond indices",
                map_opt(
                    map(
                        context("bond line", parse_bond_line),
                        |(idx1, idx2, edge_type)| {
                            (node_index(idx1 - 1), node_index(idx2 - 1), edge_type)
                        },
                    ),
                    |(idx1, idx2, et)| {
                        graph
                            .node_weight(idx1)
                            .and(graph.node_weight(idx2))
                            .map(|_| (idx1, idx2, et))
                    },
                ),
            )(running_input)?;

            // add edges to mol graph
            graph.add_edge(idx1, idx2, edge_type);
            running_input = input;
        }

        let (input, _) = context(
            "end bond block",
            tuple((tag("M  V30 END BOND"), line_ending)),
        )(running_input)?;
        input
    } else {
        input
    };

    let mol = Mol {
        graph,
        hydrogens,
        charges,
        aromatic_atoms,
    };

    // Rest of Mol file does not really matter (until we need to add 3D structure)
    let (input, _) = context("SGROUP block", opt(parse_sgroup_block))(input)?;
    let (input, _) = context("OBJ3D block", opt(parse_three_dim_block))(input)?;

    let (input, _) = tuple((tag("M  V30 END CTAB"), line_ending))(input)?;

    Ok((input, mol))
}

#[derive(Debug)]
pub struct MolFileParseError {
    trace: String,
}

impl fmt::Display for MolFileParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.trace)
    }
}

impl Error for MolFileParseError {}

pub fn parse_mol_file(input: &str) -> Result<Mol, MolFileParseError> {
    let result = all_consuming(complete(tuple((
        parse_header_block,
        parse_connection_table_block,
        tag("M  END"),
        opt(line_ending),
    ))))(input);

    match result {
        Ok((_, (_, mol, _, _))) => Ok(mol),
        Err(x) => match x {
            NomErr::Error(x) | NomErr::Failure(x) => Err(MolFileParseError {
                trace: convert_error(input, x),
            }),
            _ => unreachable!(),
        },
    }
}

#[cfg(test)]
mod tests;
