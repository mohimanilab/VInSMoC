#![warn(missing_docs)]

//! A crate to handle MC (Monte Carlo) and MCMC (Markov Chain Monte Carlo) sampling methods and compute pvalues given scoring functions.

use balloon::hist::SaturatingHist;
use rand::Rng;

mod msdpr;

const MC_MAX_SAMPLE_NUM: usize = 20;

/// A generic trait to define how to modify an implementer for an MC step.
///
/// The generic type `State` defines the state needed to revert a modify step, in case it is rejected
/// in MCMC sampling. The trait provides a default implementation for i.i.d sampling. Since `State`
/// is an associated type implementers of `Mutable` must define type wrappers if they want
/// different states.
pub trait Mutable: Sized {
    /// Any memory needed to revert a single modification step on this type.
    ///
    /// When [`Mutable::modify`] is called it returns a state to track all information needed to
    /// revert that modification with [`Mutable::revert`].
    type State;

    /// A function to modify [`Self`], returns `State`.
    ///
    /// This is used as an MCMC step. It should represent a slight modification to [`Self`], and return
    /// an instance of `State` which is needed to revert the step.
    fn modify<R: Rng + ?Sized>(&mut self, rng: &mut R) -> Self::State;

    /// A function to revert self to undo a MCMC step.
    ///
    /// This should always be called with the `State` that is returned by [`Self::modify`].
    fn revert(&mut self, state: Self::State);

    /// A default implementation to draw an i.i.d sample from the distribution over the
    /// implementor.
    ///
    /// This is used to draw an i.i.d sample for MC simulation. The base implementation is simply
    /// calling [`Self::modify`] a random number of times, where each call to [`Self::modify`]
    /// results in a slight modification to [`Self`]. See [`Mutable::iterated_sample`] for
    /// information about the return value.
    fn sample<R: Rng + ?Sized>(&mut self, rng: &mut R) -> Vec<Self::State> {
        let rand_num_iter = rng.gen_range(1..MC_MAX_SAMPLE_NUM);
        self.iterated_sample(rand_num_iter, rng)
    }

    /// A function to draw a pseudo i.i.d sample from the distribution over the implementer,
    /// returning a list of states visited.
    ///
    /// This is used to draw a pseudo i.i.d sample for MC simulation. The base implementation takes
    /// a [`usize`] `n`, and calls [`Self::modify`] `n` times, where each call to [`Self::modify`]
    /// results in a slight modification to [`Self`]. All states visited are returned and this
    /// sample can be undone by calling [`Mutable::revert`] on the output in reverse order.
    fn iterated_sample<R: Rng + ?Sized>(&mut self, n: usize, rng: &mut R) -> Vec<Self::State> {
        let mut states = Vec::with_capacity(n);
        for _ in 0..n {
            states.push(self.modify(rng));
        }
        states
    }

    /// Perform a single step of the Metropolis-Hastings algorithm.
    ///
    /// This will [`Mutable::modify`] the element once, compare the new score to the old score, and
    /// accept the new state with probability `new_score / old_score`. Optionally `log_weights` can
    /// be provided to weight scores when computing the ratio, which is relevant for [`MsDpr`]
    /// significance computation methods.
    fn mh_step<S: Scorer<Self>, R: Rng + ?Sized>(
        &mut self,
        scorer: &S,
        old_score: f64,
        rng: &mut R,
        log_weights: Option<&SaturatingHist<f64>>,
    ) -> f64 {
        let state = self.modify(rng);
        let candidate_score = scorer.score(self);
        let ratio = match log_weights {
            Some(log_weights) => {
                let num = *log_weights.get(candidate_score);
                let denom = *log_weights.get(old_score);
                (num - denom).exp().min(1.0)
            }
            None => (candidate_score / old_score).min(1.0),
        };

        // accept new candidate with prob of ratio
        // otherwise revert and reject
        if rng.gen_bool(ratio) {
            candidate_score
        } else {
            self.revert(state);
            old_score
        }
    }
}

/// Something that can score some data.
///
/// A `Scorer` should contain all the state outside of the data point necessary for computing the
/// score. For example, when scoring molecules against mass spectra and computing significance
/// given a single spectrum the scorer should contain that spectrum.
///
/// The `data` used when scoring is behind a mutable reference to allow for more complex scoring
/// methods that mutate the type they are scoring. This may not always be necessary, but requiring
/// a mutable reference gives more options for implementers of this trait.
///
/// This trait is already implemented on closures that take in a `&mut T`, so a new type does not
/// always need to be defined in order to score data points.
pub trait Scorer<T> {
    /// Compute the score of `data`.
    ///
    /// This score must return a floating point value, but otherwise does not have any strong
    /// restrictions to allow for maximum flexibility when computing p-values.
    fn score(&self, data: &mut T) -> f64;
}

impl<F, T> Scorer<T> for F
where
    F: Fn(&mut T) -> f64,
{
    fn score(&self, data: &mut T) -> f64 {
        self(data)
    }
}

/// An extension to a [`Scorer`] for when scores have a theoretical lower and upper bound that can
/// be computed from a data point.
///
/// For example, if our score is a probability we can easily define [`BoundedScorer::min_score`] to
/// be `0.0` and [`BoundedScorer::max_score`] to be `1.0`. Bounded scores are not typically
/// required to compute significance, except in the case of [`MsDpr`].
pub trait BoundedScorer<T>: Scorer<T> {
    /// The minimum possible score that this scorer could compute, given something that has been
    /// [`Mutable::modify`]ed from `data`.
    fn min_score(&self, data: &mut T) -> f64;

    /// The maximum possible score that this scorer could compute, given something that has been
    /// [`Mutable::modify`]ed from `data`.
    fn max_score(&self, data: &mut T) -> f64;

    /// An increment between two scores to use as a bucket size.
    fn score_step(&self, data: &mut T) -> f64;
}

/// A method for computing p-values.
///
/// The p-values will be based on mutating a type `T` and scoring with scorer `S`.
pub trait Method<T, S>
where
    S: Scorer<T>,
    T: Mutable,
{
    /// Compute the significance of `data` when scored with `scorer`.
    ///
    /// If there was some failure in the computation of the p-value (such as failure to converge),
    /// this should return [`None`].
    fn pvalue<R: Rng + ?Sized>(self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64>;
}

/// Monte Carlo method for estimation of significance.
///
/// The Monte Carlo method generates random samples from a distribution to approximate that
/// distribution. In this case we use [`Mutable::sample`] [`MonteCarlo::iter`] times to approximate
/// the distribution of scores. The [`Method::pvalue`] is then the number of samples that got a
/// score at least as high as the starting score.
///
/// P-value computation using this method is infallible.
pub struct MonteCarlo {
    /// Number of iterations to perform
    pub iter: u32,
}

impl<T, S> Method<T, S> for MonteCarlo
where
    S: Scorer<T>,
    T: Mutable + Clone,
{
    fn pvalue<R: Rng + ?Sized>(self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64> {
        // Keep track of how many samples score higher than the start point
        let mut num_samples_higher_score = 0;
        let start_score = scorer.score(data);

        for _ in 0..self.iter {
            let states = data.sample(rng);
            let score = scorer.score(data);
            if score >= start_score {
                num_samples_higher_score += 1;
            }
            for state in states.into_iter().rev() {
                data.revert(state);
            }
        }

        // +1, since by definition score(&mut start) >= start_score
        Some(f64::from(num_samples_higher_score + 1) / f64::from(self.iter + 1))
    }
}

/// A modified version of Monte Carlo significance estimation that looks at samples around a
/// target.
///
/// Monte Carlo significance estimation relies on the ability to draw IID random samples from the
/// space of inputs, functionality that we expose with [`Mutable::sample`]. Now, in some cases it's
/// expensive (or maybe even not possible) to generate true IID samples. In those cases
/// [`Mutable::sample`] falls back to calling [`Mutable::modify`] in a loop. The local Monte Carlo
/// significance computation method provides an alternative formulation that only looks at some
/// number of single mutations of the original input. This is identical to [`MonteCarlo`] except
/// that the iterations consist of single calls to [`Mutable::modify`] instead of calls to
/// [`Mutable::sample`].
///
/// P-value computation using this method is infallible.
///
/// # See Also
/// - [`MonteCarlo`]
pub struct LocalMonteCarlo {
    /// Number of iterations to perform.
    pub iter: u32,
}

impl<T, S> Method<T, S> for LocalMonteCarlo
where
    S: Scorer<T>,
    T: Mutable,
{
    fn pvalue<R: Rng + ?Sized>(self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64> {
        let mut num_samples_higher_score = 0;
        let start_score = scorer.score(data);

        for _ in 0..self.iter {
            let state = data.modify(rng);
            let score = scorer.score(data);
            if score >= start_score {
                num_samples_higher_score += 1;
            }
            data.revert(state);
        }

        // +1, since by definition score(&mut start) >= start_score
        Some(f64::from(num_samples_higher_score + 1) / f64::from(self.iter + 1))
    }
}

/// The Metropolis-Hastings algorithm for computing p-values.
///
/// The Metropolis-Hastings algorithm begins with a starting state and at each step
/// [`Mutable::modify`]s it once. If the new state has a better score it's instantly accepted, if
/// not there is still a probability that the worse new state will be accepted. For more
/// information on this see [`Mutable::mh_step`].
///
/// This method also optionally allows for `log_weights` of scores to be provided. These weights
/// will be used in each step of the this algorithm for computing acceptance probabilities. If not
/// provided the scores will simply be compared instead.
///
/// This implementation of Metropolis-Hastings has no special stopping or convergence conditions,
/// instead stopping after [`MetropolisHastings::iter`] steps have been performed.
///
/// P-value computation using this method is infallible.
pub struct MetropolisHastings {
    /// The number of steps to perform before stopping.
    pub iter: u32,
    /// The optional logarithm weights for scores used when computing acceptance probabilities.
    pub log_weights: Option<SaturatingHist<f64>>,
}

impl<T, S> Method<T, S> for MetropolisHastings
where
    S: Scorer<T>,
    T: Mutable,
{
    fn pvalue<R: Rng + ?Sized>(self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64> {
        let mut num_samples_higher_score = 0u32;
        let mut total_prob = 0.0;
        let start_score = scorer.score(data);
        let mut curr_score = start_score;

        let normalizer = self
            .log_weights
            .as_ref()
            .map(|weights| weights.values().map(|x| (-x).exp()).sum::<f64>());

        for _ in 0u32..self.iter {
            curr_score = data.mh_step(scorer, curr_score, rng, self.log_weights.as_ref());
            if curr_score >= start_score {
                if let Some(weights) = self.log_weights.as_ref() {
                    total_prob += (-weights.get(curr_score)).exp() / normalizer.unwrap();
                } else {
                    num_samples_higher_score += 1;
                }
            }
        }

        Some(match &self.log_weights {
            Some(_) => total_prob / f64::from(self.iter),
            // +1, since by definition score(&mut start) >= start_score
            None => f64::from(num_samples_higher_score + 1) / f64::from(self.iter + 1),
        })
    }
}

/// The MS-DPR method of p-value computation.
///
/// This is a weighted Metropolis-Hastings using Wang-Landau sampling to compute the weights. This
/// stops after performing [`MsDpr::iter`] steps of Metropolis-Hastings after computing the score
/// weights. If you would like a more sophisticated (but more expensive) stopping criterion use
/// [`MsDprTavc`].
///
/// P-value computation with this method is fallible. If Wang-Landau Sampling fails to converge to
/// a flat histogram at any point this will return the placeholder p-value of [`None`].
pub struct MsDpr {
    /// Number of Metropolis-Hastings steps to perform between convergence checks.
    pub iter: u32,
    /// Upper bound on number of Metropolis-Hastings steps done in each iteration of Wang-Landau
    /// Sampling.
    ///
    /// If the histogram is not flat after this many iterations we immediately move on to the next
    /// constant factor.
    pub max_wang_landau_iter: usize,
}

impl<T, S> Method<T, S> for MsDpr
where
    S: BoundedScorer<T>,
    T: Mutable + Clone,
{
    fn pvalue<R: Rng + ?Sized>(self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64> {
        let weights = msdpr::WangLandauHist::new(
            scorer.min_score(data),
            scorer.max_score(data),
            scorer.score_step(data),
        )
        .into_log_weights(data, scorer, rng, self.max_wang_landau_iter)?;

        MetropolisHastings {
            iter: self.iter,
            log_weights: Some(weights),
        }
        .pvalue(data, scorer, rng)
    }
}

/// The MS-DPR method of p-value computation with time average variance constant stopping
/// condition.
///
/// This is a weighted Metropolis-Hastings with a more sophisticated stopping criterion. For
/// weights it uses Wang-Landau sampling to compute log score weights to redistribute probability
/// when performing [`Mutable::mh_step`]s. For its stopping criterion it computes a confidence
/// internal using a time-average variance constant and stops when that interval is smaller than
/// some percentage of an estimate of the sample variance.
///
/// P-value computation with this method is fallible. If Wang-Landau Sampling fails to converge to
/// a flat histogram at any point this will return the placeholder p-value of [`None`].
pub struct MsDprTavc {
    /// Number of Metropolis-Hastings steps to perform between convergence checks.
    pub iter: u32,
    /// Upper bound on number of Metropolis-Hastings steps done in each iteration of Wang-Landau
    /// Sampling.
    ///
    /// If the histogram is not flat after this many iterations we immediately move on to the next
    /// constant factor.
    pub max_wang_landau_iter: usize,
}

impl<T, S> Method<T, S> for MsDprTavc
where
    S: BoundedScorer<T>,
    T: Mutable + Clone,
{
    fn pvalue<R: Rng + ?Sized>(self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64> {
        let weights = msdpr::WangLandauHist::new(
            scorer.min_score(data),
            scorer.max_score(data),
            scorer.score_step(data),
        )
        .into_log_weights(data, scorer, rng, self.max_wang_landau_iter)?;

        let og_score = scorer.score(data);
        let mut traj = msdpr::TrajectoryInfo::new(&weights, og_score);

        // burn in for 5 times the number of iterations
        let mut curr_score = og_score;
        for _ in 0..5 * self.iter as usize {
            curr_score = data.mh_step(scorer, curr_score, rng, Some(&weights));
            traj.add(curr_score);
        }

        loop {
            for _ in 0..self.iter {
                curr_score = data.mh_step(scorer, curr_score, rng, Some(&weights));
                traj.add(curr_score);
            }

            if traj.stop_iteration(0.02, 0.0) {
                break;
            }
        }
        Some(traj.mean())
    }
}

/// A [`Method`] that picks between [`MonteCarlo`], [`LocalMonteCarlo`], [`MetropolisHastings`],
/// [`MsDpr`], and [`MsDprTavc`].
#[derive(Debug, Clone, Copy, clap::Parser)]
#[clap(help_heading = "P-Value")]
pub struct MethodSelection {
    /// Number of steps to perform.
    ///
    /// When performing Metropolis-Hastings for p-values this is simply the number of steps of
    /// mutation applied to the original structure. If performing MS-DPR with constant iterations
    /// this is the number of Metropolis-Hastings steps after Wang-Landau sampling to compute score
    /// weights. If performing MS-DPR with time average variance constant stopping conditions this
    /// is the number of Metropolis-Hastings steps between convergence checks. If performing Monte
    /// Carlo simulations this is just the number of samples generated.
    #[clap(long, short = 'i', default_value = "100000")]
    iter: u32,
    /// Use Monte Carlo.
    ///
    /// This overrides the default Metropolis-Hastings p-value implementation and switches it to
    /// Monte Carlo with a constant number of iterations instead. It is an error to enable Monte
    /// Carlo and MS-DPR or Monte Carlo and Local Monte Carlo at the same time.
    #[clap(long, conflicts_with_all = &["msdpr", "local-monte-carlo"])]
    monte_carlo: bool,
    /// Use Local Monte Carlo.
    ///
    /// This overrides the default Metropolis-Hastings p-value implementation and switches it to
    /// Local Monte Carlo with a constant number of iterations instead. Local Monte Carlo is
    /// identical to Monte Carlo except that instead of IID samples it uses single-step mutations
    /// of the original input. It is an error to enable Local Monte Carlo and MS-DPR or Local Monte
    /// Carlo and Monte Carlo at the same time.
    #[clap(long, conflicts_with_all = &["msdpr", "monte-carlo"])]
    local_monte_carlo: bool,
    /// Use the MS-DPR method.
    ///
    /// This overrides the default Metropolis-Hastings p-value implementation and switches it to
    /// MS-DPR with a constant number of iterations. This is essentially weighted
    /// Metropolis-Hastings where the weights are computed using Wang-Landau sampling. This is
    /// typically slower than simple Metropolis-Hastings but can get p-values for rarer events.
    #[clap(long)]
    msdpr: bool,
    /// Enable time average variance constant stopping conditions for MS-DPR.
    ///
    /// This has no effect unless perform MS-DPR for p-values. This overrides the stopping
    /// condition and switches it to a check of a confidence interval based on a time average
    /// variance constant estimation via overlapping estimation with a trapezoidal selection rule.
    /// This is more expensive than constant iteration MS-DPR but is more theoretically sound. In
    /// practice this does not have a strong effect and is not recommended due to its slower
    /// runtime.
    #[clap(long)]
    tavc: bool,
    /// The maximum number of Metropolis-Hastings steps to use for a single constant factor during
    /// Wang-Landau sampling.
    ///
    /// Wang-Landau sampling goes through exponentially decreasing constant factors and runs
    /// Metropolis-Hastings simulations for these constant factors, re-weighting until the observed
    /// score histogram is relatively flat. Now, given a sufficiently effective mutation function
    /// this should always converge. We have found that there are some cases where convergence
    /// takes a very long time. To address this we introduce an upper bound to the number of
    /// iterations to cause early halting. P-values that are a result of this early halting are not
    /// typically very informative, but this allows us to move on when computing many p-values at
    /// once instead of waiting for a very long time for a single p-value.
    ///
    /// This has no effect unless --msdpr is set.
    #[clap(long, default_value = "100000")]
    max_wang_landau_iter: usize,
    /// The cutoff to do p-value consistency checking by re-running.
    ///
    /// By random chance it's possible that you can get lower p-values that you should. These
    /// spurious significant hits can be highly undesirable when accuracy is paramount. Providing a
    /// floating point value to this flag turns on consistency checking mode. In this mode any
    /// p-value that is below the specified floating point value is re-run --rerun-iter times
    /// and the median is reported. In the worst case this can cause large slowdowns to it's
    /// disabled by default.
    #[clap(long)]
    rerun_cutoff: Option<f64>,
    /// The number of times to check consistency of a low p-value by re-running.
    ///
    /// This flag is meaningless unless --rerun-cutoff is provided. When a p-value is below
    /// the --rerun-cutoff this flag is the number of times the p-value should be
    /// re-computed. The final reported value in this setting is the median of all runs of the same
    /// p-value. See --rerun-cutoff for more information on consistency checking.
    #[clap(long, default_value = "10", requires = "rerun-cutoff", value_parser = clap::value_parser!(u32).range(2..))]
    rerun_iter: u32,
}

impl Default for MethodSelection {
    fn default() -> Self {
        Self {
            iter: 100000,
            monte_carlo: false,
            local_monte_carlo: false,
            msdpr: false,
            tavc: false,
            max_wang_landau_iter: 100000,
            rerun_cutoff: None,
            rerun_iter: 10,
        }
    }
}

/// Get a string literal with the description of [`MethodSelection`] flags and their usage.
#[macro_export]
macro_rules! method_selection_about {
    () => {
        r#"P-value simulation by default is done with Metropolis-Hastings for --iter steps. The --monte-carlo flag switches this over to Monte Carlo simulation for --iter samples. The --local-monte-carlo flag uses a variant of Monte Carlo that only samples one mutation from the original input molecule, again for --iter samples. The --msdpr flag switches over to the MS-DPR algorithm, which first computes score weights via Wang-Landau sampling and then runs --iter steps of weighted Metropolis-Hastings. Adding the --tavc flag in addition to --msdpr switches the stopping condition from a constant number of Metropolis-Hastings steps post Wang-Landau sampling to a stopping condition based on the Time Average Variance Constant. This stopping condition is more theoretically sound and interesting, but in practice does not seem to increase consistency or accuracy while significantly increasing runtimes. If --msdpr is set you can control the upper bound on the number of Metropolis-Hastings steps for a single constant factor in Wang-Landau sampling with the --max-wang-landau-iter flag. Higher values will result in more consistent p-values, but may cause runtime to increase.

If you suspect that p-values are highly inconsistent in your runs the --rerun-iter and --rerun-cutoff flags can help, at the cost of runtime overhead. If --rerun-cutoff is set then all p-values that are below this cutoff are run again for --rerun-iter, with the final p-value being the median of all these reruns.
"#;
    }
}

impl MethodSelection {
    /// Use the [`MetropolisHastings`] method with `iter` steps.
    ///
    /// This does no consistency checks.
    pub fn metropolis_hastings(iter: u32) -> Self {
        Self {
            iter,
            ..Self::default()
        }
    }

    /// Use the [`MonteCarlo`] method with `iter` samples.
    ///
    /// This does no consistency checks.
    pub fn monte_carlo(iter: u32) -> Self {
        Self {
            iter,
            monte_carlo: true,
            ..Self::default()
        }
    }

    /// Use the [`LocalMonteCarlo`] method with `iter` samples.
    ///
    /// This does no consistency checks.
    pub fn local_monte_carlo(iter: u32) -> Self {
        Self {
            iter,
            local_monte_carlo: true,
            ..Self::default()
        }
    }

    /// Use the MS-DPR method with `iter` steps.
    ///
    /// When `tavc` is enabled this is [`MsDprTavc`] with `iter` Metropolis-Hastings steps
    /// between convergence checks. Otherwise, this is [`MsDpr`] with `iter` Metrpolis-Hatings
    /// steps after Wang-Landau sampling. The upper bound on number of Metroplis-Hastings steps per
    /// constant factor during Wang-Landau sampling can be controlled with `max_wang_landau_iter`.
    ///
    /// This does no consistency checks.
    pub fn ms_dpr(iter: u32, tavc: bool, max_wang_landau_iter: usize) -> Self {
        Self {
            iter,
            msdpr: true,
            tavc,
            max_wang_landau_iter,
            ..Self::default()
        }
    }

    /// Use p-value consistency checking with this method.
    ///
    /// Consistency checking does `iter` re-runs of all p-values below `cutoff`. The final p-value
    /// becomes the median of these re-runs. If consistency was already enabled this will overwrite
    /// the consistency parameters. See also [`MethodSelection::inconsistent`].
    ///
    /// If `iter` isn't at least `2` this will do nothing.
    ///
    /// Consistency checking can be a nice backup measure when p-values are highly sensitive to
    /// initialization, but does add a large cost. First, in the worst case runtime with be `iter`
    /// times slower. This case is unlikely, but even if no p-values pass the consistency check
    /// `cutoff` every time a p-value is computed there will be a extra clone of the data.
    pub fn consistent(&mut self, cutoff: f64, iter: u32) -> &mut Self {
        if iter < 2 {
            return self;
        }
        self.rerun_cutoff = Some(cutoff);
        self.rerun_iter = iter;
        self
    }

    /// Disable p-value consistency checking with this method.
    ///
    /// Resets this method to the default consistency checking, which is none. See
    /// [`MethodSelection::consistent`] for a description of consistency checking.
    pub fn inconsistent(&mut self) -> &mut Self {
        self.rerun_cutoff = None;
        self.rerun_iter = Self::default().rerun_iter;
        self
    }

    /// If consistency checks are enabled return the cutoff and number of extra iterations.
    pub fn consistency(&self) -> Option<(f64, u32)> {
        self.rerun_cutoff.map(|cutoff| (cutoff, self.rerun_iter))
    }

    fn one_pvalue<T, S, R>(&self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64>
    where
        T: Mutable + Clone,
        S: BoundedScorer<T>,
        R: Rng + ?Sized,
    {
        if self.monte_carlo {
            return MonteCarlo { iter: self.iter }.pvalue(data, scorer, rng);
        }

        if self.local_monte_carlo {
            return LocalMonteCarlo { iter: self.iter }.pvalue(data, scorer, rng);
        }

        if self.msdpr {
            if self.tavc {
                MsDprTavc {
                    iter: self.iter,
                    max_wang_landau_iter: self.max_wang_landau_iter,
                }
                .pvalue(data, scorer, rng)
            } else {
                MsDpr {
                    iter: self.iter,
                    max_wang_landau_iter: self.max_wang_landau_iter,
                }
                .pvalue(data, scorer, rng)
            }
        } else {
            MetropolisHastings {
                iter: self.iter,
                log_weights: None,
            }
            .pvalue(data, scorer, rng)
        }
    }
}

impl<T: Mutable + Clone, S: BoundedScorer<T>> Method<T, S> for MethodSelection {
    /// This p-value implementation selects one of [`MetropolisHastings`], [`MsDpr`], and
    /// [`MsDprTavc`] based on the parameters of [`MethodSelection`].
    ///
    /// The different constructors influence which method is used for the p-value computation here.
    /// - [`MethodSelection::metropolis_hastings`] uses [`MetropolisHastings`]
    /// - [`MethodSelection::monte_carlo`] uses [`MonteCarlo`]
    /// - [`MethodSelection::local_monte_carlo`] uses [`LocalMonteCarlo`]
    /// - [`MethodSelection::ms_dpr`] uses [`MsDprTavc`] if `tavc` was set, otherwise it uses
    /// [`MsDpr`]
    ///
    /// If rerun-based consistency checks were enabled with [`MethodSelection::consistent`] then
    /// checks will be performed on the specified cutoff for the specified number of iterations.
    fn pvalue<R: Rng + ?Sized>(self, data: &mut T, scorer: &S, rng: &mut R) -> Option<f64> {
        // consistency checks require that we have a fresh backup `data`
        // however, this clone is wasted computation if we aren't doing consistency checks.
        // we can't even just compute the pvalue and check the cutoff, because we need the clone
        // first to avoid messing up `data`
        // so, we'll store in here if there's a chance we'll need consistency checks
        let consistency_backup = self.rerun_cutoff.is_some().then(|| data.clone());

        let pvalue = self.one_pvalue(data, scorer, rng)?;

        // we're either not doing consistency or we're too insignificant
        if self.rerun_cutoff.is_none() || self.rerun_cutoff.unwrap() < pvalue {
            return Some(pvalue);
        }

        // recompute as many pvalues as asked for
        let iter = self.rerun_iter as usize;
        let mut pvalues = Vec::with_capacity(iter);
        for _ in 0..iter {
            *data = consistency_backup.as_ref().unwrap().clone();
            pvalues.push(self.one_pvalue(data, scorer, rng)?);
        }
        pvalues.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap());

        // get median
        let mid = iter / 2;
        Some(if iter % 2 == 0 {
            (pvalues[mid - 1] + pvalues[mid]) / 2.0
        } else {
            pvalues[mid]
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use rand::rngs::mock::StepRng;

    // Test a mutable pvalue function on Method

    #[derive(Clone)]
    struct FloatWrapper {
        val: f32,
        incr: f32,
    }

    impl Mutable for FloatWrapper {
        type State = Option<f32>;

        fn modify<R: Rng + ?Sized>(&mut self, _rng: &mut R) -> Option<f32> {
            self.incr += 1.0;
            self.val = self.incr;
            Some(self.val)
        }

        fn revert(&mut self, state: Option<f32>) {
            match state {
                Some(prev_val) => {
                    self.val = prev_val;
                }
                None => {}
            }
        }
    }

    #[test]
    fn try_mutable_mc() {
        // rng.gen_range(1..x) always returns 1
        let mut rng = StepRng::new(1, 0);
        let mut start = FloatWrapper {
            val: 800.0,
            incr: 0.0,
        };
        let mc = MonteCarlo { iter: 1000 };
        let pval = mc.pvalue(
            &mut start,
            &|x: &mut FloatWrapper| x.val as f64 / 1000.0,
            &mut rng,
        );
        // Draw values in range from 0 to 1000, so expect 200 values over 800
        assert_abs_diff_eq!(pval.unwrap(), 0.2, epsilon = 1e-2);
    }

    #[test]
    fn try_mutable_metropolis_hastings() {
        let mut rng = StepRng::new(0, 1000);
        let mut start = FloatWrapper {
            val: 800.0,
            incr: 0.0,
        };
        let mc = MetropolisHastings {
            iter: 1000,
            log_weights: None,
        };
        let pval = mc.pvalue(
            &mut start,
            &(|x: &mut FloatWrapper| x.val as f64 / 1000.0),
            &mut rng,
        );
        // Draw values in range from 0 to 1000, so expect 200 values over 800
        assert_abs_diff_eq!(pval.unwrap(), 0.2, epsilon = 1e-2);
    }

    #[test]
    fn try_capture_scope() {
        let mut rng = StepRng::new(0, 1000);
        let mut start = FloatWrapper {
            val: 800.0,
            incr: 0.0,
        };
        let mc = MetropolisHastings {
            iter: 10000,
            log_weights: None,
        };
        // Verify that scope can be captured by closures
        let scale = 1.0;
        let pval = mc.pvalue(
            &mut start,
            &|x: &mut FloatWrapper| x.val as f64 / 10000.0 * scale,
            &mut rng,
        );
        // Draw values in range of 0 to 10000, so expect 9200 values over 800
        assert_abs_diff_eq!(pval.unwrap(), 0.92, epsilon = 1e-3);
    }
}
