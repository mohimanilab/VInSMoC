use std::env;

// Converts a peptide string into a .MOL file.

fn main() {
    let sequence = env::args()
        // skip the name of the binary (which is always arg 0) and get the next argument, this returns an option
        .nth(1)
        // that we require to be Some, not None, so we can unwrap here
        .unwrap();
    // throw away the first return value, which we won't use and keep just the second return value in the output tuple
    let (_, peptide) = peptide::parsers::parse_one_letter_peptide(&sequence).unwrap();
    let mol = peptide.to_mol();
    mol.to_mol_file(&mut std::io::stdout()).unwrap();
}
