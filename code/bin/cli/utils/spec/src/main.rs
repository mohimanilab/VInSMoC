//! Utility methods for working with spectra.

use anyhow::{Context, Result};
use balloon::{
    fuzzy::{FuzzyMap, Tol},
    of64,
};
use clap::Parser;
use constants::MASS_WATER_EXACT;
use itertools::Itertools;
use paths::NewDirPath;
use peptide::AminoAcid;
use rustc_hash::FxHashSet;
use spectrum::{Ms, ParsedCollection, Peaks, Spectrum, Splash};
use std::io::BufRead;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Utility methods for working with spectra.
///
/// All commands require a path to a spectrum file in MGF, mzML, or mzXML format. This is done at
/// the very top level before any subcommands are specified. Any spectrum pre-processing steps can
/// be specified there too, although these are most commonly left to defaults.
///
/// For sub-commands that support it, a subset of the scans in the input spectrum file can be
/// specified. This subset can take the form of spectrum indices (0-indexed and only including MS2
/// scans) or scan numbers. If both are specified then the sub-command will be run on the scans
/// specified by index first and then the scans specified by scan number. In both cases the outputs
/// will be in the same order as the inputs. For example, --indices 0 1 and --scans 0 1 will output
/// on scan with index 0, then scan with index 1, then scan with scan number 0, and finally scan
/// with scan number 1.
#[derive(Debug, clap::Parser)]
#[clap(name = "spec")]
struct Opts {
    #[clap(subcommand)]
    sub: SubCommand,
    #[clap(flatten)]
    spec: ParsedCollection,
}

#[derive(Debug, clap::Parser)]
struct SpectrumSubset {
    /// Scan numbers of scans of interest.
    ///
    /// Scan numbers are taken directly from the spectrum file for mzML and mzXML formatted inputs.
    /// For MGF inputs scan numbers are assigned, starting with 0, to all MS2 scans. The scan
    /// subset specified by --scans will be operated on after those specified by --indices. If
    /// neither --indices nor --scans are given then all MS2 scans in the input file will be
    /// operated on.
    #[clap(long, multiple_values = true, takes_value = true)]
    scans: Vec<usize>,
    /// Indices of scans of interest.
    ///
    /// Scan indices are 0-based indices starting from the first MS2 scan and counting upwards
    /// (only considered MS2 scans). In the case of MGF files they are the same as scan numbers.
    /// The scan subset specified by --indices will be operated on before those specified by
    /// --scans. If neither --indices nor --scans are given then all MS2 scans in the input file
    /// will be operated on.
    #[clap(long, multiple_values = true, takes_value = true)]
    indices: Vec<usize>,
}

#[derive(Debug, clap::Parser)]
enum SubCommand {
    /// Print out spectrum in MGF format.
    Print {
        #[clap(flatten)]
        subset: SpectrumSubset,
        /// Print in JSON format instead of MGF.
        #[clap(long)]
        json: bool,
    },
    /// Draw spectra as SVGs.
    ///
    /// Drawn spectra will be written to SVG files named after their scan numbers in the provided
    /// --dir.
    Draw {
        #[clap(flatten)]
        subset: SpectrumSubset,
        /// Directory in which the drawn spectrum SVGs will be created.
        #[clap(long)]
        dir: NewDirPath,
        /// The background color for the spectrum drawings.
        #[clap(long, short = 'b', default_value = "RGBA(0,0,0,0.0)")]
        background_color: draw::Color,
        /// The number of most instens peaks to draw.
        ///
        /// If not provided then all peaks will be drawn. Otherwise only the most intense -n peaks
        /// will be drawn in all the final images.
        #[clap(long, short = 'n')]
        num_peaks: Option<usize>,
        /// Height of SVG images.
        #[clap(long, default_value = "500")]
        height: u32,
        /// Width of SVG images.
        #[clap(long, default_value = "1000")]
        width: u32,
    },
    /// Convert spectrum index to scan number.
    ///
    /// A spectrum index is a 0-based index, excluding MS1 spectra. A scan is something specified
    /// in the collection format being used. Note, for MGF-formatted spectra these are actually
    /// equivalent.
    IndexToScan {
        /// Index of spectrum.
        #[clap(long)]
        index: usize,
    },
    /// Convert spectrum scan number to index.
    ///
    /// A spectrum index is a 0-based index, excluding MS1 spectra. A scan is something specified
    /// in the collection format being used. Note, for MGF-formatted spectra these are actually
    /// equivalent.
    ScanToIndex {
        /// Scan of spectrum.
        #[clap(long)]
        scan: usize,
    },
    /// Count the number of peaks in the target spectra.
    PeakCount {
        #[clap(flatten)]
        subset: SpectrumSubset,
    },
    /// Get the precursor m/z of the target spectra.
    Precursor {
        #[clap(flatten)]
        subset: SpectrumSubset,
    },
    /// Get the retention time of the target spectra.
    RetentionTime {
        #[clap(flatten)]
        subset: SpectrumSubset,
    },
    /// Get the charges of the target spectra.
    Charge {
        #[clap(flatten)]
        subset: SpectrumSubset,
    },
    /// Compute the SPLASH spectral hash of the target spectra.
    Splash {
        #[clap(flatten)]
        subset: SpectrumSubset,
    },
    /// Print all peaks in the target spectra with m/z and intensity separated by the specified
    /// --separator.
    Peaks {
        #[clap(flatten)]
        subset: SpectrumSubset,
        /// The separator between peak m/z and intensity values.
        #[clap(long, required = true, default_value = ",")]
        separator: char,
    },
    /// Print the m/zs of all peaks in the target spectra.
    Mzs {
        #[clap(flatten)]
        subset: SpectrumSubset,
    },
    /// Print the intensities of all peaks in the target spectra.
    Intensities {
        #[clap(flatten)]
        subset: SpectrumSubset,
    },
    /// Print the sequence tags found given the peaks in the target spectra.
    Tags {
        #[clap(flatten)]
        subset: SpectrumSubset,
        /// The length of tags.
        ///
        /// If 0, the maximal tag length will be inferred.
        #[clap(long, default_value = "3")]
        tag_length: usize,
        /// The floating point precision.
        #[clap(long, default_value = "1e-3")]
        precision: f64,
        /// An optional configuration file.
        ///
        /// The file must be specified as a TSV, with monomer name in the first column, and monomer
        /// mass in the second column.
        #[clap(long)]
        monomer_config: Option<std::path::PathBuf>,
        /// Sets the mass loss, or the reduction in mass when a monomer is disconnected.
        ///
        /// Should be either H2O or H2.
        #[clap(long, default_value = "H2O", possible_values = &["H2O", "H2"])]
        mass_loss: String,
    },
}

impl SubCommand {
    pub fn scans(&self) -> Option<&[usize]> {
        use SubCommand::*;
        match self {
            Print { subset, .. }
            | Draw { subset, .. }
            | PeakCount { subset }
            | Precursor { subset }
            | RetentionTime { subset }
            | Charge { subset }
            | Splash { subset }
            | Peaks { subset, .. }
            | Mzs { subset }
            | Tags { subset, .. }
            | Intensities { subset } => (!subset.scans.is_empty()).then(|| subset.scans.as_ref()),
            SubCommand::IndexToScan { .. } => None,
            SubCommand::ScanToIndex { .. } => None,
        }
    }

    pub fn indices(&self) -> Option<&[usize]> {
        use SubCommand::*;
        match self {
            Print { subset, .. }
            | Draw { subset, .. }
            | PeakCount { subset }
            | Precursor { subset }
            | RetentionTime { subset }
            | Charge { subset }
            | Splash { subset }
            | Peaks { subset, .. }
            | Mzs { subset }
            | Tags { subset, .. }
            | Intensities { subset } => {
                (!subset.indices.is_empty()).then(|| subset.indices.as_ref())
            }
            SubCommand::IndexToScan { .. } => None,
            SubCommand::ScanToIndex { .. } => None,
        }
    }
}

/// Execute `op` on all relevant scans on the spectrum in `opts`.
///
/// If neither `indices` or `scans` are specified this runs on all of the scans. Otherwise it runs
/// on the scans specified, running `op` first on the scans specified by index and then on scans
/// specified by scan number.
fn run_on_collection<F>(opts: &Opts, op: F) -> Result<()>
where
    F: Fn(&Spectrum) -> Result<()>,
{
    let collection = opts.spec.collection()?;
    if opts.sub.scans().is_none() && opts.sub.indices().is_none() {
        for scan in &collection.spectra {
            op(scan)?;
        }
    } else {
        if let Some(indices) = opts.sub.indices() {
            for &index in indices {
                let scan = &collection.spectra.get(index).with_context(|| {
                    format!(
                        "index {} doesn't exist in {}",
                        index,
                        opts.spec.path.display()
                    )
                })?;
                op(scan)?;
            }
        }

        if let Some(scans) = opts.sub.scans() {
            for &scan in scans {
                let scan = collection.get_scan(scan).with_context(|| {
                    format!(
                        "scan {} doesn't exist in {}",
                        scan,
                        opts.spec.path.display()
                    )
                })?;
                op(scan)?;
            }
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let opts = Opts::parse();

    match opts.sub {
        SubCommand::Print { json, .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                if json {
                    serde_json::to_writer_pretty(std::io::stdout(), scan)?;
                } else {
                    scan.write_mgf(std::io::stdout())?;
                }
                Ok(())
            };

            run_on_collection(&opts, op)?;
        }
        SubCommand::Draw {
            ref dir,
            background_color,
            num_peaks,
            width,
            height,
            ..
        } => {
            dir.create()?;
            let op = |scan: &Spectrum| -> Result<()> {
                let scan_num = scan.scan;
                let file = dir.join(format!("scan{}.svg", scan_num));
                let mut drawing =
                    spectrum::draw::Drawing::new(scan).background_color(&background_color);
                if let Some(num_peaks) = num_peaks {
                    drawing = drawing.limit_peaks(num_peaks);
                }
                drawing.draw_svg(&file, (width, height))?;
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::IndexToScan { index } => {
            println!(
                "{}",
                opts.spec
                    .collection()?
                    .spectra
                    .get(index)
                    .with_context(|| format!(
                        "index {} doesn't exist in {}",
                        index,
                        opts.spec.path.display()
                    ))?
                    .scan
            );
        }
        SubCommand::ScanToIndex { scan } => {
            println!(
                "{}",
                opts.spec
                    .collection()?
                    .spectra
                    .iter()
                    .position(|spec| spec.scan == scan)
                    .with_context(|| format!(
                        "scan {} doesn't exist in {}",
                        scan,
                        opts.spec.path.display()
                    ))?
            );
        }
        SubCommand::PeakCount { .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                println!("{}", scan.peaks().len());
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::Precursor { .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                println!("{}", scan.pepmass());
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::RetentionTime { .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                if let Some(rt) = scan.retention_time() {
                    print!("{}", rt);
                }
                println!();
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::Charge { .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                println!(
                    "{}",
                    scan.charges()
                        .iter()
                        .map(|charge| charge.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                );
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::Splash { .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                println!("{}", Splash::from(scan));
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::Peaks { separator, .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                for peak in scan.peaks() {
                    println!("{}{}{}", peak.mz, separator, peak.intensity);
                }
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::Mzs { .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                for peak in scan.peaks() {
                    println!("{}", peak.mz);
                }
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::Intensities { .. } => {
            let op = |scan: &Spectrum| -> Result<()> {
                for peak in scan.peaks() {
                    println!("{}", peak.intensity);
                }
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
        SubCommand::Tags {
            tag_length,
            precision,
            ref monomer_config,
            ref mass_loss,
            ..
        } => {
            let op = |scan: &Spectrum| -> Result<()> {
                // Need to get the residual masses for all standard and nonstandard amino acids
                let all_aas = peptide::StdOrNonStdAa::get_all_aas();

                let mass_loss = match mass_loss.as_str() {
                    "H2O" => MASS_WATER_EXACT,
                    "H2" => 2.015650,
                    _ => unreachable!(),
                };

                let (monomer_masses, monomer_names) = if let Some(monomer_config) =
                    monomer_config.as_ref()
                {
                    let mut monomer_masses = Vec::new();
                    let mut monomer_names = Vec::new();
                    let reader = std::io::BufReader::new(std::fs::File::open(monomer_config)?);
                    for line in reader.lines() {
                        let line = line?;
                        let mut split = line.split('\t');
                        let name = split.next().expect("missing monomer name");
                        let mass = split.next().expect("missing monomer mass");
                        let mass = mass
                            .parse::<f64>()
                            .expect("failed to parse monomer mass as f64");

                        if let Ok(mass) = of64::try_from(mass - mass_loss) {
                            monomer_masses.push(mass);
                            monomer_names.push(name.to_string());
                        }
                    }
                    (monomer_masses, monomer_names)
                } else {
                    let monomer_masses = all_aas
                        .iter()
                        .filter_map(|aa| of64::try_from(aa.exact_mass() - MASS_WATER_EXACT).ok())
                        .collect::<Vec<_>>();

                    let monomer_names = all_aas.iter().map(|aa| aa.to_string()).collect::<Vec<_>>();

                    (monomer_masses, monomer_names)
                };

                let sg = scan.spectrum_graph(
                    &monomer_masses,
                    &monomer_names,
                    scan.pepmass(),
                    precision,
                )?;

                let paths = if tag_length != 0 {
                    sg.sequence_tags(tag_length)
                } else {
                    // infer the maximal tag length
                    let mut tags = sg.sequence_tags(1);

                    for i in 1..sg.size() {
                        let new_tags = sg.sequence_tags(i + 1);
                        if new_tags.is_empty() {
                            break;
                        }
                        tags = new_tags;
                    }

                    tags
                };

                let tol = Tol::Abs(of64::try_from(precision).unwrap());
                let mut monomer_mass_to_aas: FuzzyMap<Vec<&str>> = FuzzyMap::new(tol);
                for (mass, aa_name) in monomer_masses.iter().zip(monomer_names.iter()) {
                    if let Some(aas_with_mass) = monomer_mass_to_aas.get_mut(*mass) {
                        aas_with_mass.push(aa_name);
                    } else {
                        monomer_mass_to_aas.insert(*mass, vec![aa_name]);
                    }
                }

                let mut all_tags = FxHashSet::default();

                // print the outputs for each sequence tag
                for path in paths {
                    // amino acid tags, and mass peaks for them. this outputs a single tag
                    // with arbitrary choices for which aas are used for a mass shift. we will
                    // print out all possible amino acids at a given position with a given mass
                    let (_, masses) = path;

                    let mass_diffs = masses.windows(2).map(|w| w[1] - w[0]).collect::<Vec<_>>();

                    // print out all amino acids at each position with the same mass,
                    // using `monomer_mass_to_aas`, which gives a list of all amino acids with a
                    // given mass

                    let empty_vec = vec![];

                    // Generate combinations using iterators
                    let combinations = mass_diffs
                        .iter()
                        .map(|mass_diff| {
                            monomer_mass_to_aas
                                .get(of64::try_from(*mass_diff).ok().unwrap())
                                .unwrap_or(&empty_vec)
                        })
                        .multi_cartesian_product();

                    // Print all combinations
                    for aas in combinations {
                        let aa_tag = aas.iter().join(",");
                        if !all_tags.insert(aa_tag.clone()) {
                            continue;
                        }

                        print!("{}\t", scan.scan);
                        print!("{}\t", scan.pepmass());
                        print!("{}\t", scan.retention_time().unwrap_or(-1.0));
                        print!("{}\t", scan.precursor_intensity().unwrap_or(-1.0));
                        print!("{}\t", aa_tag);
                        println!("{}", masses.iter().join(","));
                    }
                }
                Ok(())
            };
            run_on_collection(&opts, op)?;
        }
    }
    Ok(())
}
