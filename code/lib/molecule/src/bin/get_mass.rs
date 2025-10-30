use molecule::{parsers::parse_smiles, ChemFormula};

use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, StdinLock, Write};

enum InputSource<'a> {
    Stdin(StdinLock<'a>),
    File(BufReader<File>),
}

impl Read for InputSource<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Stdin(x) => x.read(buf),
            Self::File(x) => x.read(buf),
        }
    }
}

impl BufRead for InputSource<'_> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::Stdin(x) => x.fill_buf(),
            Self::File(x) => x.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::Stdin(x) => x.consume(amt),
            Self::File(x) => x.consume(amt),
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // setup input
    let args: Vec<String> = std::env::args().skip(1).collect();
    let stdin = io::stdin();
    let stdin_lock = stdin.lock();
    let readers: Vec<InputSource> = match args.len() {
        0 => vec![InputSource::Stdin(stdin_lock)],
        _ => args
            .iter()
            .map(|file| {
                InputSource::File(BufReader::new(
                    File::open(file).unwrap_or_else(|_| panic!("Could not open file {}", file)),
                ))
            })
            .collect(),
    };

    // setup output
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    for reader in readers {
        for line in reader.lines() {
            dbg!(&line);
            let line = line.expect("Failed to read line from input");
            match parse_smiles(&line) {
                Ok(mol) => {
                    let exact_mass = mol.exact_mass();
                    let charge = mol.formal_charge();
                    let formula = ChemFormula::from(&mol);
                    writeln!(handle, "{} {} {}", exact_mass, charge, formula)?;
                }
                Err(_) => {
                    eprintln!("Failed to parse smiles {}, writing empty line", line);
                    writeln!(handle)?;
                }
            };
        }
    }

    Ok(())
}
