use clap::Parser;
use molecule::MolBuilder;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Utility methods for working with chemical reactions.
#[derive(Debug, Parser)]
#[clap(name = "rxn", version = env!("CARGO_PKG_VERSION"))]
struct Opts {
    #[clap(subcommand)]
    sub: Subcommand,
}

#[derive(Debug, Parser)]
enum Subcommand {
    /// Print out the submolecule pattern in MOL V3000 format.
    ///
    /// Note that the MOL format uses 1-based indexing but the modifications use 0-based indexing.
    /// User provides a reactant and product using the indexed SMILES format.
    CenterPattern {
        /// The reactant indexed SMILES
        #[clap(long)]
        reactant: String,
        /// The product indexed SMILES
        #[clap(long)]
        product: String,
    },
    /// Print out the modifications to convert the reactant to the product.
    ///
    /// Modifications are sequential graph-based changes to a molecule.
    ///
    /// Note that node indices here are 0-based and the MOL pattern has 1-based indices.
    CenterMods {
        /// The reactant indexed SMILES
        #[clap(long)]
        reactant: String,
        /// The product indexed SMILES
        #[clap(long)]
        product: String,
    },
}

fn main() -> anyhow::Result<()> {
    let opts = Opts::parse();

    match opts.sub {
        Subcommand::CenterPattern { reactant, product } => {
            let reactant = MolBuilder::from_indexed_smiles(&reactant)?;
            let product = MolBuilder::from_indexed_smiles(&product)?;
            let subgraph_mod = modification::reaction(&reactant, &product, "".into());
            let (mol, _) = subgraph_mod.pattern.reset_indices();
            let mut stdout = std::io::stdout();
            mol.to_mol_file(&mut stdout)?;
        }
        Subcommand::CenterMods { reactant, product } => {
            let reactant = MolBuilder::from_indexed_smiles(&reactant)?;
            let product = MolBuilder::from_indexed_smiles(&product)?;
            let subgraph_mod = modification::reaction(&reactant, &product, "".into());
            let (_, idx_map) = subgraph_mod.pattern.reset_indices();
            let idx_map = idx_map
                .into_iter()
                .map(|(k, v)| (v.index() as i32, k.index() as i32))
                .collect();
            let mods = subgraph_mod.mod_seq.convert(&idx_map)?;
            print!("{}", mods);
        }
    }

    Ok(())
}
