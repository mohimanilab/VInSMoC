use anyhow::Context;
use clap::Parser;
use modification::SubgraphMod;
use molecule::MolBuilder;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Utility for applying modifications to a molecule.
///
/// For our seq2x methods we store molecule modifications on disk in a directory containing two
/// files: a mod.mol file describing a molecular substructure and a mod.txt that modifies that
/// structure with graph operations (add/delete a node/edge). This utility applies the modification
/// at the specified directory to a molecule encoded as SMILES.
#[derive(Debug, Parser)]
#[clap(name = "modify", version = env!("CARGO_PKG_VERSION"))]
struct Opts {
    /// Path to a directory containing the target modification.
    ///
    /// This directory should contain a mod.mol file containing the target pattern and a mod.txt
    /// file containing the modification steps.
    #[clap(long, short = 'm')]
    modification: paths::DirPath,
    /// SMILES of molecule to be modified.
    #[clap(long, short = 's')]
    smiles: String,
    /// Print out the original, unmodified SMILES as well.
    #[clap(long)]
    include_original: bool,
    /// Print out the atom indices of the modification site.
    ///
    /// This converts the output into a TSV where the first column is the modified SMILES and the
    /// second column is a semicolon-separated list of atom indices describing the modified site.
    #[clap(long)]
    include_site: bool,
}

fn main() -> anyhow::Result<()> {
    let opts = Opts::parse();

    let mol_file = opts.modification.join("mod.mol");
    let mod_file = opts.modification.join("mod.txt");
    let subgraph_mod = SubgraphMod::from_mol_and_mod_file(mol_file, mod_file, "modification")?;

    let mol = MolBuilder::from_smiles(&opts.smiles)
        .with_context(|| "invalid input SMILES")?
        .implicit_hs()
        .build();

    if opts.include_original {
        if opts.include_site {
            println!("{}\t", mol.to_smiles());
        } else {
            println!("{}", mol.to_smiles());
        }
    }

    for node_map in subgraph_mod.conversions(&mol) {
        let node_map = node_map
            .into_iter()
            .map(|(x, y)| (x.index() as i32, y.index() as i32))
            .collect();
        let mod_seq = subgraph_mod
            .mod_seq
            .clone()
            .convert(&node_map)
            .with_context(|| "failed to convert indices from subgraph isomorphism")?;

        let mut mol = mol.clone();
        mod_seq
            .modify(&mut mol)
            .with_context(|| "failed to apply modification")?;

        if opts.include_site {
            let site = mod_seq
                .modified_site(&mol)?
                .into_iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>();
            let site = site.join(";");
            println!("{}\t{}", mol.to_smiles(), site);
        } else {
            println!("{}", mol.to_smiles());
        }
    }

    Ok(())
}
