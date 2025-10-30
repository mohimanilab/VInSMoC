use std::collections::HashMap;

use clap::Parser;
use molecule::Atom;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Utility methods for working with atoms.
#[derive(Debug, Parser)]
#[clap(name = "atom", version = env!("CARGO_PKG_VERSION"))]
struct Opts {
    #[clap(subcommand)]
    sub: Subcommand,
}

#[derive(Debug, Parser)]
enum Subcommand {
    /// Print out the exact mass of this atom
    ExactMass {
        #[clap(long, short = 'a', required = true)]
        atom: Atom,
    },
    /// Print out the standard mass of this atom
    StdMass {
        #[clap(long, short = 'a', required = true)]
        atom: Atom,
    },
    /// Print out the periodic table of masses as JSON.
    PeriodicTable,
}

fn main() -> anyhow::Result<()> {
    let opts = Opts::parse();

    match opts.sub {
        Subcommand::ExactMass { atom } => {
            println!("{}", atom.exact_mass());
        }
        Subcommand::StdMass { atom } => {
            println!("{}", atom.std_mass());
        }
        Subcommand::PeriodicTable => {
            let mut masses = HashMap::new();
            for atomic_num in 1..=118 {
                let atom = Atom::try_from(atomic_num)?;
                masses.insert(atom, atom.exact_mass());
            }
            serde_json::to_writer_pretty(std::io::stdout(), &masses)?;
            // otherwise we wont have a newline
            println!();
        }
    }

    Ok(())
}
