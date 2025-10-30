use std::borrow::Borrow;

use anyhow::Context;
use clap::Parser;
use dereplicate::{DereplicationRun, Msm};
use paths::FilePath;
use rustc_hash::FxHashMap;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Debug, serde::Deserialize)]
struct Solution {
    spectrum_file: String,
    scan: usize,
    inchikey14: String,
}

#[derive(Debug, Parser)]
enum Subcommand {
    /// Computes to the top K accuracy for some set of values for K for this search, printing out a
    /// CSV with the values for K and the corresponding accuracy.
    ///
    /// The top K accuracy is the percent of spectrum for which the correct molecule is in the top
    /// K molecules. In the case of ties in scoring the ranks of all molecules with the same score
    /// are set to the average of their original ranks.
    TopK {
        /// The values of K for which we want to compute top K accuracy.
        ///
        /// Multiple whitespace-delimited values can be provided.
        #[clap(long, short = 'k', multiple_values = true)]
        ks: Vec<u8>,
    },
    /// Return the scans for which the correct molecule is within the top k highest scores hits.
    ///
    /// Correct scan numbers will be reported with one on each line.
    CorrectScans {
        /// The number of top-scoring hits that are considered when tracking correct scans.
        #[clap(long, short = 'k', default_value = "1")]
        k: u8,
    },
    /// Computes the top N accuracy for some set of values for N for this search, printing out a
    /// CSV with the values for N and the corresponding accuracy.
    TopN {
        #[clap(long, short = 'n', multiple_values = true)]
        ns: Vec<usize>,
    },
}

#[derive(Debug, Parser)]
struct Opts {
    #[clap(subcommand)]
    sub: Subcommand,
    /// Path to TSV file output from a dereplication method.
    #[clap(long)]
    hits: FilePath,
    /// Path to CSV file containing spectrum_file, scan, and inchikey14 columns in any order.
    ///
    /// This file associates spectra to the InChIKey of their compounds. All subcommands require
    /// some reference solution.
    #[clap(long)]
    solutions: FilePath,
}

/// An identifier for a spectrum, consisting of the spectrum file name and target scan.
type SpecId = (String, usize);

/// Parse the input data.
///
/// The `hits` are parsed into a [`DereplicationRun`] with no minimum score. The `solutions` are
/// first parsed from the on-disk CSV into a [`Solution`], then these are converted to a hash map
/// from spectrum identifiers (file path and scan tuples) to InChIKeys.
///
/// The returned value is a tuple of the parsed run and the solutions hash map.
///
/// # Errors
/// IO and CSV deserialization errors are propagated.
fn parse_data(
    hits: FilePath,
    solutions: FilePath,
) -> anyhow::Result<(DereplicationRun, FxHashMap<SpecId, String>)> {
    // parse hits TSV
    let hits = DereplicationRun::from_path(hits, None)?;

    // parse reference solution CSV
    let mut solutions =
        csv::Reader::from_path(solutions).with_context(|| "failed to open solutions CSV")?;
    let solutions = solutions
        .deserialize::<Solution>()
        .collect::<Result<Vec<_>, _>>()
        .with_context(|| "failed to parse solutions CSV")?;
    let solutions = solutions
        .into_iter()
        .map(|sol| ((sol.spectrum_file, sol.scan), sol.inchikey14))
        .collect::<FxHashMap<_, _>>();

    Ok((hits, solutions))
}

/// Compute the ranks of the provided `hits`.
///
/// Ranks are computed according to position when `hits` is sorted by descending [`Msm::score`]. In
/// the case of ties in score all hits with that score are assigned the average rank. For example,
/// if hits ranked `3` and `4` both have a score of `10.0` then both will receive a rank of `3.5`.
///
/// This function sorts the provided `hits` in place by descending score. The returned value is a
/// vector of the final ranks, with ties broken as described above. This vector is co-indexed with
/// the sorted version of `hits`, meaning that the rank of `hits[0]` is `ranks(hits)[0]`.
///
/// The generic parameter `T` allows this function to work on both slices of [`Msm`]s and slices of
/// references to [`Msm`]s.
fn ranks<T: Borrow<Msm<'static>>>(hits: &mut [T]) -> Vec<f64> {
    hits.sort_by(|m1, m2| {
        m1.borrow()
            .score
            .partial_cmp(&m2.borrow().score)
            .unwrap()
            .reverse()
    });
    let mut ranks = (1..=hits.len()).map(|x| x as f64).collect::<Vec<_>>();

    let mut start = 0;
    let mut end = 0;
    let mut prev_score = hits[0].borrow().score;
    for hit in hits.iter().skip(1) {
        if hit.borrow().score == prev_score {
            end += 1;
        } else {
            let avg = ranks[start..=end].iter().sum::<f64>();
            let avg = avg / ranks[start..=end].len() as f64;
            ranks[start..=end].fill(avg);
            start = end + 1;
            end = start;
        }
        prev_score = hit.borrow().score;
    }
    let avg = ranks[start..=end].iter().sum::<f64>();
    let avg = avg / ranks[start..=end].len() as f64;
    ranks[start..=end].fill(avg);

    ranks
}

/// Computes a mapping from spectrum identifiers to the rank of the correct match for that
/// spectrum.
///
/// If no correct hit exists for a spectrum its value in the returned map will be [`f64::MAX`].
fn best_rank_per_spectrum<'a>(
    hits: &'a DereplicationRun,
    solutions: &FxHashMap<SpecId, String>,
) -> FxHashMap<(&'a str, usize), f64> {
    let grouped = hits.group_by_spectrum(None);

    grouped
        .into_iter()
        .map(|(k, mut v)| {
            let ranks = ranks(&mut v);
            let rank = match v.into_iter().position(|hit| {
                hit.inchikey14 == solutions[&(hit.spectrum_file.to_string(), hit.scan)]
            }) {
                Some(idx) => ranks[idx],
                None => f64::MAX,
            };
            (k, rank)
        })
        .collect::<FxHashMap<_, _>>()
}

fn main() -> anyhow::Result<()> {
    let opts = Opts::parse();
    let (mut hits, solutions) = parse_data(opts.hits, opts.solutions)?;

    match opts.sub {
        Subcommand::TopK { ks } => {
            let grouped = best_rank_per_spectrum(&hits, &solutions);
            println!("k,accuracy");
            for k in ks {
                let correct = grouped.values().filter(|&&x| x <= k as f64).count() as f64;
                let correct = correct / solutions.len() as f64;
                println!("{},{}", k, correct);
            }
        }
        Subcommand::CorrectScans { k } => {
            let grouped = best_rank_per_spectrum(&hits, &solutions);
            for (_, scan) in grouped
                .into_iter()
                .filter(|(_, v)| *v <= k as f64)
                .map(|(key, _)| key)
            {
                println!("{}", scan);
            }
        }
        Subcommand::TopN { ns } => {
            let ranks = ranks(hits.matches_mut());

            println!("n,accuracy");
            for n in ns {
                let mut cutoff_idx = ranks.len();

                // find the upper bound hit index (first index that exceeds target rank)
                // falls back to the full hit list just in case n exceeds number of hits
                for (idx, &rank) in ranks.iter().enumerate() {
                    if rank > n as f64 {
                        cutoff_idx = idx;
                        break;
                    }
                }

                // check how many of those hits are correct
                let mut accuracy = hits.matches()[..cutoff_idx]
                    .iter()
                    .filter(|&hit| {
                        hit.inchikey14 == solutions[&(hit.spectrum_file.to_string(), hit.scan)]
                    })
                    .count() as f64;
                accuracy /= cutoff_idx as f64;
                println!("{},{}", n, accuracy);
            }
        }
    };

    Ok(())
}
