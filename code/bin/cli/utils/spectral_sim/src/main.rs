use clap::Parser;
use paths::FilePath;
use spectrum::{Peaks, SpectrumCollection, SpectrumPreprocessing};

/// Compute a distance matrix between spectra in the two spectrum files.
///
/// Distances are defined as the dot products between spectra.
#[derive(Debug, Parser)]
struct Opts {
    /// Path to the first spectrum file.
    #[clap(long)]
    spectra1: FilePath,
    #[clap(long)]
    /// Path to the second spectrum file.
    spectra2: FilePath,
    #[clap(flatten)]
    pre: SpectrumPreprocessing,
}

fn dot_product(s1: &spectrum::Spectrum, s2: &spectrum::Spectrum) -> f64 {
    let mut s1 = s1.peaks().iter();
    let mut s2 = s2.peaks().iter();

    let mut peak1 = s1.next();
    let mut peak2 = s2.next();
    let mut product = 0.0;

    while peak1.is_some() && peak2.is_some() {
        let inner1 = peak1.unwrap();
        let inner2 = peak2.unwrap();
        if (inner1.mz - inner2.mz).abs() <= 0.01 {
            product += inner1.intensity * inner2.intensity;
        }

        if inner1.mz < inner2.mz {
            peak1 = s1.next();
        } else {
            peak2 = s2.next();
        }
    }

    product
}

fn main() -> anyhow::Result<()> {
    let mut opts = Opts::parse();

    opts.pre.normalize = true;

    let mut coll1 = SpectrumCollection::from_path(&opts.spectra1)?;
    coll1.preprocess(&opts.pre);
    let mut coll2 = SpectrumCollection::from_path(&opts.spectra2)?;
    coll2.preprocess(&opts.pre);

    for s1 in &coll1.spectra {
        let sims = coll2
            .spectra
            .iter()
            .map(|s2| dot_product(s1, s2).to_string())
            .collect::<Vec<_>>();
        println!("{}", sims.join(","))
    }

    Ok(())
}
