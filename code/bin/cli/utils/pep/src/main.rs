//! Utility methods for working with peptides.

use std::convert::TryFrom;

use clap::Parser;
use peptide::{Peptide, StdAminoAcid};

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Utility methods for working with peptides
///
/// Information is extracted via the avaiable subcommands. No matter the subcommand the user must
/// provide a single peptide using the --peptide flag. This should be a string whose characters are
/// the single letter amino acids codes for the 20 standard amino acids.
#[derive(Debug, Parser)]
#[clap(name = "pep", version = env!("CARGO_PKG_VERSION"))]
struct Opts {
    /// The amino acid sequence of peptide using the single-letter amino acid code.
    #[clap(long, short = 'p', help_heading = "Input")]
    peptide: String,
    #[clap(subcommand)]
    sub: Subcommand,
}

#[derive(Debug, Parser)]
enum Subcommand {
    /// Print out the smiles of this peptide
    Smiles,
    /// Print out the smiles of the head-to-tail cyclized version of this peptide.
    CyclicSmiles,
    /// Print out the V3000 MOL file for this peptide
    MolFile,
    /// Shuffle amino acids in this peptide and print out the shuffled version
    Shuffle,
}

fn main() -> anyhow::Result<()> {
    let opts = Opts::parse();

    let mut pep = opts
        .peptide
        .chars()
        .map(StdAminoAcid::try_from)
        .collect::<Result<Peptide, _>>()?;

    match opts.sub {
        Subcommand::Smiles => println!("{}", pep.to_mol().to_smiles()),
        Subcommand::CyclicSmiles => println!("{}", pep.to_cyclic().to_smiles()),
        Subcommand::MolFile => pep.to_mol().to_mol_file(&mut std::io::stdout())?,
        Subcommand::Shuffle => {
            let mut rng = rand::thread_rng();
            pep.shuffle(&mut rng);
            println!(
                "{}",
                pep.into_iter()
                    .map(|aa| char::from(aa.as_std().unwrap()))
                    .collect::<String>()
            );
        }
    }

    Ok(())
}
