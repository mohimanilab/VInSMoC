//! Utility methods for working with molecules.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
};

use anyhow::Result;
use chemical_db::{ChemicalDb, DatabaseRow};
use clap::Parser;
use molecule::{Atom, ChemFormula, Mol, MolBuilder, Point};
use paths::FilePath;
use petgraph::{
    dot::{Config, Dot},
    visit::{EdgeRef, IntoEdgeReferences, IntoNodeReferences, NodeRef},
};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use serde_json::json;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Utility methods for working with molecules.
///
/// Information is extracted via the available subcommands. No matter the subcommand user must
/// provide the SMILES of the molecule we want to examine. A single smiles can be provided with the
/// --smiles file. A file containing one SMILES on each line can be provided with --file. If
/// neither of these are given then one SMILES will be read at a time from STDIN.
///
/// Some subcommands produce multi-line outputs (atoms, bonds, dot). These are not intended for
/// multi-molecule inputs and therefore will only apply to the first SMILES (however it was input).
#[derive(Debug, Parser)]
#[clap(name = "mol", version = env!("CARGO_PKG_VERSION"))]
struct Opts {
    /// SMILES of a single input molecule.
    #[clap(long, short = 's')]
    smiles: Option<String>,
    /// Path to a file containing a single SMILES on each line.
    #[clap(long, short = 'f', conflicts_with = "smiles")]
    file: Option<FilePath>,
    #[clap(long, short = 'i')]
    implicit_hs: bool,
    #[clap(subcommand)]
    sub: Subcommand,
}

#[derive(Debug, Parser)]
enum Subcommand {
    /// Print out the number of non-hydrogen atoms in the molecule.
    HeavyAtomCount,
    /// Print out the number of atoms represented as molecular graph nodes.
    ExplicitAtomCount,
    /// Print out the number of atoms in the molecule.
    AtomCount,
    /// Print out the chemical formula of the molecule.
    Formula,
    /// Print out the exact mass of the molecule.
    ExactMass {
        #[clap(long, default_value = "3")]
        precision: u8,
    },
    /// Print out the InChIKey of the molecule.
    InchiKey,
    /// Print out the ECFP of the molecule.
    ///
    /// This is a hash based on the ECFP of the provided --radius and modified to be similar to the
    /// InChIKey in that it is a 14-character ASCII string.
    Ecfp {
        /// The radius to use when computed the ECFP hash.
        #[clap(long, default_value = "8")]
        radius: u8,
    },
    /// Print out the number of bonds in the molecule.
    BondCount,
    /// Print out each atom in the molecule.
    ///
    /// The order of atoms corresponds to the order of atoms in the internal representation of the
    /// molecule. Only atoms that receive their own node in the molecular graph will be printed.
    /// When both --hydrogens and --charges are enabled this will print out a three column CSV with
    /// the first column being the atomic symbol, the second column being the number of implicit
    /// hydrogens connected to that atom, and the third column being the charge of that atom.
    Atoms {
        /// Show implicit hydrogen counts
        #[clap(long)]
        hydrogens: bool,
        /// Show atomic charges
        #[clap(long)]
        charges: bool,
    },
    /// Print out the multiplicity of each bond in the molecule.
    ///
    /// The order of bonds corresponds to the order of edges in the internal representation of the
    /// molecular graph. Each bond is weighted by a bond type and each bond type will be on its
    /// own line. A bond type is represented by the SMILES bond type specification character. If
    /// the --atoms option is set each bond will be flanked by the atomic symbols of the atoms it
    /// connects. The atom will the smaller atomic number will always appear on the left side of
    /// the bond.
    ///
    /// See also <https://daylight.com/dayhtml/doc/theory/theory.smiles.html> Section 3.2.2
    Bonds {
        /// Show atomic symbols on either of the bond
        #[clap(long)]
        atoms: bool,
    },
    /// Print out the graphviz DOT representation of the molecular graph.
    ///
    /// Nodes will be labeled with their atoms. Edges will be labeled with their bond type. If
    /// --hydrogens is enabled then nodes will also be labeled with hydrogen counts. If --charges
    /// is enabled then nodes will further also be labeled with their charge.
    Dot {
        /// Show implicit hydrogen counts
        #[clap(long)]
        hydrogens: bool,
        /// Show atomic charges
        #[clap(long)]
        charges: bool,
    },
    /// Compute 2D coordinates of the molecular graph for drawing.
    ///
    /// Output is JSON with two top level keys for atoms and bonds. The atoms are a JSON array of
    /// the atomic symbol, the atom's x coordinate, and the atom's y coordinate. The edges are a
    /// JSON array of the indices of the two atoms in the bond and the multiplicity of the bond as
    /// an integer.
    Coordinates {
        /// Pretty print the JSON output
        #[clap(long)]
        pretty: bool,
    },
    /// Construct a chemical database CSV for input to a dereplication method.
    ///
    /// This will print a CSV with four columns: name,exact_mass,smiles,inchikey14. This is the
    /// format accepted by our dereplication methods. There will be one row for each input
    /// molecule. The name column will be the same as the inchikey14 column unless overridden with
    /// the --names option.
    Database {
        /// Use InChIKey prefixes for the molecular hashes.
        ///
        /// The default behavior is to use ECFP hashes for molecular hashes as they're faster to
        /// compute.
        #[clap(long)]
        inchikey: bool,
        /// A list of whitespace-separated names for the molecules in this database.
        ///
        /// If omitted the name column with be set to be the same as the inchikey14 column. If
        /// given this list must match the length of the input molecules. Only the minimum of the
        /// compound count and name count rows will be reported in the final output.
        #[clap(long)]
        names: Option<Vec<String>>,
        /// Sort molecules by ascending exact mass before printing.
        ///
        /// The default behavior simply outputs molecules as they appear in the input SMILES list.
        #[clap(long)]
        sort: bool,
    },
    /// Get all molecules / BGCs containing the input smiles.
    ///
    /// This will print a CSV with two columns: BGC,smiles, such that all the outputted smiles
    /// contain the input smiles as a substructure. The input `database` should be a comma
    /// separated CSV of name,smiles for each BGC (where name is some identifier for the BGC).
    ///
    /// # Panics
    /// This function panics if `database` is not a valid CSV file.
    Substruct {
        /// A database file, which should be a csv of name,smiles for each BGC.
        #[clap(long, short = 'd')]
        database: paths::FilePath,
    },
    /// Print out the input molecule as a V3000 MOL file.
    ///
    /// This only works with single-molecule inputs.
    MolFile,
}

static mut STDIN: Option<std::io::Stdin> = None;

fn setup_input(
    smiles: Option<&String>,
    file: Option<&FilePath>,
) -> Result<Box<dyn Iterator<Item = Result<Mol>>>> {
    if let Some(smiles) = smiles {
        let mol = MolBuilder::from_smiles(smiles);
        return Ok(Box::new(std::iter::once_with(|| Ok(mol?.build()))));
    }
    if let Some(file) = file {
        let file = File::open(file)?;
        let reader = BufReader::new(file);
        return Ok(Box::new(
            reader
                .lines()
                .map(|smiles| Ok(MolBuilder::from_smiles(&smiles?)?.build())),
        ));
    }

    let reader = unsafe {
        STDIN = Some(std::io::stdin());
        STDIN.as_ref()
    }
    .unwrap();
    let reader: std::io::StdinLock<'static> = reader.lock();
    Ok(Box::new(reader.lines().map(|smiles| {
        Ok(MolBuilder::from_smiles(&smiles?)?.build())
    })))
}

fn print_all<T, F>(it: Box<dyn Iterator<Item = Result<Mol>>>, f: F) -> Result<()>
where
    F: Fn(&Mol) -> T,
    T: std::fmt::Display,
{
    let stdout = std::io::stdout();
    let mut stdout = BufWriter::new(stdout.lock());

    for mol in it {
        let mol = mol?;
        writeln!(stdout, "{}", f(&mol))?;
    }

    Ok(())
}

fn main() -> Result<()> {
    let opts = Opts::parse();

    let mut mol_iter = setup_input(opts.smiles.as_ref(), opts.file.as_ref())?;

    match opts.sub {
        Subcommand::HeavyAtomCount => {
            print_all(mol_iter, |mol| {
                mol.graph
                    .node_weights()
                    .filter(|&&atom| atom != Atom::H)
                    .count()
            })?;
        }
        Subcommand::ExplicitAtomCount => {
            print_all(mol_iter, |mol| mol.graph.node_count())?;
        }
        Subcommand::AtomCount => {
            print_all(mol_iter, |mol| {
                mol.graph.node_count() + mol.hydrogens().values().copied().sum::<u8>() as usize
            })?;
        }
        Subcommand::Formula => {
            #[allow(clippy::redundant_closure)]
            print_all(mol_iter, |mol| ChemFormula::from(mol))?;
        }
        Subcommand::ExactMass { precision } => {
            print_all(mol_iter, |mol| {
                format!("{1:.0$}", precision as usize, mol.exact_mass())
            })?;
        }
        Subcommand::InchiKey => {
            print_all(mol_iter, |mol| {
                mol.to_inchi_key().expect("failure to generate inchikey")
            })?;
        }
        Subcommand::Ecfp { radius } => {
            print_all(mol_iter, |mol| {
                let hash = mol.ecfp_hash(radius);
                String::from_utf8(hash.to_vec()).expect("ECFP hash wasn't valid utf8")
            })?;
        }
        Subcommand::BondCount => {
            print_all(mol_iter, |mol| {
                mol.graph.edge_count() + mol.hydrogens().values().copied().sum::<u8>() as usize
            })?;
        }
        Subcommand::Atoms { hydrogens, charges } => {
            if let Some(mol) = mol_iter.next() {
                let mol = mol?;
                for i in mol.graph.node_indices() {
                    print!("{}", mol.graph[i]);
                    if hydrogens {
                        print!(",{}", mol.h_count(&i));
                    }
                    if charges {
                        print!(",{}", mol.charge(&i));
                    }
                    println!();
                }
            }
            if mol_iter.next().is_some() {
                eprintln!("WARNING: the atoms subcommand only applies to the first molecule, but you provided more than one input molecule");
            }
        }
        Subcommand::Bonds { atoms } => {
            if let Some(mol) = mol_iter.next() {
                let mol = mol?;
                for b in mol.graph.edge_references() {
                    let mult = b.weight();

                    if atoms {
                        let mut atom1 = mol.graph[b.source()];
                        let mut atom2 = mol.graph[b.target()];
                        if atom2.atomic_num() < atom1.atomic_num() {
                            std::mem::swap(&mut atom1, &mut atom2);
                        }
                        println!("{}{}{}", atom1, mult, atom2);
                    } else {
                        println!("{}", mult);
                    }
                }
            }
            if mol_iter.next().is_some() {
                eprintln!("WARNING: the bonds subcommand only applies to the first molecule, but you provided more than one input molecule");
            }
        }
        Subcommand::Dot { hydrogens, charges } => {
            if let Some(mol) = mol_iter.next() {
                let mol = mol?;
                println!(
                    "{}",
                    Dot::with_attr_getters(
                        &mol.graph,
                        &[Config::EdgeNoLabel, Config::NodeNoLabel],
                        &|_, edge| format!("label = \"{}\"", edge.weight()),
                        &|_, (idx, atom)| {
                            let mut label = atom.to_string();
                            if hydrogens {
                                label += &format!(" h:{}", mol.h_count(&idx));
                            }
                            if charges {
                                label += &format!(" c:{}", mol.charge(&idx));
                            }
                            format!("label = \"{}\"", label)
                        }
                    )
                );
            }
            if mol_iter.next().is_some() {
                eprintln!("WARNING: the dot subcommand only applies to the first molecule, but you provided more than one input molecule");
            }
        }
        Subcommand::Coordinates { pretty } => {
            if let Some(mol) = mol_iter.next() {
                let mol = mol?;
                let coords = mol.coordinates();
                let atoms = mol
                    .graph
                    .node_references()
                    .zip(&coords)
                    .map(|(node, Point { x, y })| {
                        json!({
                            "symbol": node.weight().to_string(),
                            "x": x,
                            "y": y
                        })
                    })
                    .collect::<Vec<_>>();
                let bonds = mol
                    .graph
                    .edge_references()
                    .map(|edge| {
                        json!({
                            "atom1": edge.source().index(),
                            "atom2": edge.target().index(),
                            "mult": *edge.weight() as usize
                        })
                    })
                    .collect::<Vec<_>>();

                let mol_json = json!({
                    "atoms": atoms,
                    "bonds": bonds,
                });

                if pretty {
                    println!("{}", serde_json::to_string_pretty(&mol_json)?);
                } else {
                    println!("{}", serde_json::to_string(&mol_json)?);
                }
            }
            if mol_iter.next().is_some() {
                eprintln!("WARNING: the coordinates subcommand only applies to the first molecule, but you provided more than one input molecule");
            }
        }
        Subcommand::Database {
            inchikey,
            names,
            sort,
        } => {
            let mols = mol_iter.collect::<Result<Vec<Mol>>>()?;
            let mut db = if let Some(names) = names {
                mols.into_par_iter()
                    .zip(names)
                    .map(|(mol, name)| -> Result<DatabaseRow> {
                        if inchikey {
                            let inchikey = mol.to_inchi_key().map_err(|code| {
                                anyhow::format_err!("inchikey conversion error code {}", code)
                            })?;
                            Ok(DatabaseRow::new(
                                name,
                                mol.to_smiles(),
                                mol.exact_mass(),
                                inchikey.as_bytes().try_into()?,
                            ))
                        } else {
                            Ok(DatabaseRow::from_named_mol(name, &mol))
                        }
                    })
                    .collect::<Result<ChemicalDb>>()?
            } else {
                mols.into_par_iter()
                    .map(|mol| {
                        if inchikey {
                            let inchikey = mol.to_inchi_key().map_err(|code| {
                                anyhow::format_err!("inchikey conversion error code {}", code)
                            })?;
                            Ok(DatabaseRow::new(
                                inchikey.clone(),
                                mol.to_smiles(),
                                mol.exact_mass(),
                                inchikey.as_bytes().try_into()?,
                            ))
                        } else {
                            Ok(DatabaseRow::from_mol(&mol))
                        }
                    })
                    .collect::<Result<ChemicalDb>>()?
            };

            if sort {
                db.par_sort_rows();
            }

            // dump output
            let stdout = std::io::stdout();
            let stdout = stdout.lock();
            let mut stdout = BufWriter::new(stdout);
            db.write_csv(&mut stdout)?;
            stdout.flush()?;
        }
        Subcommand::Substruct { database } => {
            if let Some(mol) = mol_iter.next() {
                let mut pattern_mol = mol?;
                pattern_mol.explicit_hs();
                let f = File::open(database)?;
                let f = BufReader::new(f);
                for line in f.lines() {
                    let line = line.unwrap();
                    // tokens = [name, smiles]
                    let tokens = line.split(',').collect::<Vec<_>>();
                    let mut db_mol = MolBuilder::from_smiles(tokens[1]).unwrap().build();
                    db_mol.explicit_hs();
                    let isomorphisms = db_mol.ri_subgraph_isomorphism(&pattern_mol);
                    if !isomorphisms.is_empty() {
                        println!("{},{}", tokens[0], tokens[1]);
                    }
                }
            }

            if mol_iter.next().is_some() {
                eprintln!("WARNING: the substruct subcommand only searches for the first molecule, but you provided more than one input molecule");
            }
        }
        Subcommand::MolFile => {
            if let Some(mol) = mol_iter.next() {
                let mol = mol?;
                mol.to_mol_file(&mut std::io::stdout())?;
            }
            if mol_iter.next().is_some() {
                eprintln!("WARNING: the atoms subcommand only applies to the first molecule, but you provided more than one input molecule");
            }
        }
    }

    Ok(())
}
