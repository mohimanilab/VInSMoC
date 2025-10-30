use balloon::hist::SaturatingHist;
use rand::Rng;

use crate::{Mutable, Scorer};

/// The histogram of score counts used when computing log score weights using the Wang-Landau
/// algorithm.
///
/// This is mainly used for its [`WangLandauHist::into_log_weights`] method to consume this type
/// and receive log score weights for use in the [`super::MsDpr`] method for p-value computation.
/// That implementation primarily follows the implementation in the original C++ codebase, which in
/// turn is based on [the updated MS-DPR paper](https://doi.org/10.4230/LIPIcs.WABI.2017.14) with
/// minor discrepancies in constants used. In cases of disagreement we follow the values from the
/// codebase.
///
/// # The Algorithm
/// The purpose of Wang-Landau sampling is to reduce the number of samples necessary for a MC
/// algorithm. It accomplishes this by weighting energy (or score) values by their frequency and
/// computing a Metropolis-Hastings acceptance ratio with those weights. Ideally we find weights
/// that perfectly correct for the probability of each score (`weight(score) = 1 / Pr(score)`). To
/// do so we try to find a sufficiently flat score histogram. The final score weights are then
/// scaled by some factor `C` (the original paper uses `f`). This factor is decreased exponentially
/// each iteration. We finish using a given factor when the histogram is sufficiently flat or we
/// hit our upper bound on iterations.
///
/// The associated constants for this type store the relevant values for factors, flatness
/// conditions, and iteration upper bounds. See also [`WangLandauHist::into_log_weights`] and
/// [`WangLandauHist::is_flat`].
///
/// # Relevant Papers
/// - [Updated MS-DPR](https://doi.org/10.4230/LIPIcs.WABI.2017.14)
/// - [Wang-Landau Sampling](https://doi.org/10.1119/1.1707017)
pub struct WangLandauHist {
    inner: SaturatingHist<u32>,
}

impl WangLandauHist {
    /// The starting factor for Wang-Landau sampling.
    ///
    /// This value is taken from [the updated MS-DPR](https://doi.org/10.4230/LIPIcs.WABI.2017.14)
    /// paper, where it is set to `0.6.exp()`. However, we store the logarithm instead of the
    /// actual exponentiated value.
    pub const C_START: f64 = 0.6;

    /// The final factor for Wang-Landau sampling.
    ///
    /// This value is taken from [the updated MS-DPR](https://doi.org/10.4230/LIPIcs.WABI.2017.14)
    /// paper, where it is set to `0.0000367.exp()`. However, we store the logarithm instead of the
    /// actual exponentiated value.
    pub const C_END: f64 = 0.0000367;

    /// The lower bound multiplicative factor for histogram flatness.
    ///
    /// This value was described in the paper to be `70%`. However, the old codebase uses `80%`
    /// instead so we follow the codebase.
    ///
    /// See [`WangLandauHist::is_flat`].
    pub const FLAT_LOWER: f64 = 0.7;

    /// The upper bound multiplicative factor for histogram flatness.
    ///
    /// This value was described in the paper to be `130%`. However, the old codebase uses `120%`
    /// instead so we follow the codebase.
    ///
    /// See [`WangLandauHist::is_flat`].
    pub const FLAT_UPPER: f64 = 1.3;

    /// Construct an empty histogram with buckets `step` apart, supporting scores between
    /// `min_score` and `max_score`.
    ///
    /// If any scores we compute leave the bounds of `min_score` and `max_score` then are truncated
    /// to `min_score` or `max_score`. It is therefore not necessarily catastrophic to incorrectly
    /// estimate these bounds, but if estimates are extremely off it can impact the reliability of
    /// the weights.
    ///
    /// See [`SaturatingHist::new`] also.
    pub fn new(min_score: f64, max_score: f64, step: f64) -> Self {
        Self {
            inner: SaturatingHist::new(min_score, max_score, step),
        }
    }

    /// Set all counts in the histogram to zero.
    pub fn zero(&mut self) {
        for cnt in self.inner.values_mut() {
            *cnt = 0;
        }
    }

    /// Increment the count of the bucket corresponding to `score` by 1.
    pub fn increment(&mut self, score: f64) {
        *self.inner.get_mut(score) += 1;
    }

    /// Get the number of buckets in this histogram.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Check if the score count histogram is sufficiently flat.
    ///
    /// We don't expect that we can ever get a perfectly flat distribution of scores, even though
    /// that is our goal, so we give some wiggle room for counting a histogram as flat. We first
    /// compute the expected count in the case where our histogram is perfectly flat. Our histogram
    /// is then sufficiently flat if all buckets have counts that between
    /// [`WangLandauHist::FLAT_LOWER`] percent and [`WangLandauHist::FLAT_UPPER`] percent of the
    /// expected perfectly flat value.
    ///
    /// One potential pitfall is that the histogram could easily appear sufficiently flat when not
    /// too many values have been added to it. For example, at the start the histogram is actually
    /// perfectly flat at all `0`s. To remedy this a histogram cannot be considered sufficiently
    /// flat if there are any buckets with counts lower than `20`. This value was taken from the
    /// old codebase and is not documented in [the paper](https://doi.org/10.4230/LIPIcs.WABI.2017.14).
    pub fn is_flat(&self) -> bool {
        let total_buckets = self.len() as f64;
        let expected_count = self.inner.values().sum::<u32>() as f64 / total_buckets;
        let count_range = (Self::FLAT_LOWER * expected_count)..=(Self::FLAT_UPPER * expected_count);

        for count in self.inner.values().copied().map(f64::from) {
            if count < 20.0 || !count_range.contains(&count) {
                return false;
            }
        }

        true
    }

    /// Compute the log score weights using Wang-Landau sampling.
    ///
    /// This implementation is based on the pseudocode in `Algorithm 3` of [this
    /// paper](https://doi.org/10.4230/LIPIcs.WABI.2017.14). Starting with the factor
    /// [`WangLandauHist::C_START`] we decrease the factor by a factor of 2 every iteration
    /// (equivalent to taking the root, since we're working in log). In each factor-iteration we
    /// perform weighted [`Mutable::mh_step`]s until either we hit `max_iter` or our internal count
    /// histogram [`WangLandauHist::is_flat`].
    ///
    /// The returned values can be used as the optional weights in [`Mutable::mh_step`].
    pub fn into_log_weights<T: Mutable + Clone, S: Scorer<T>, R: Rng + ?Sized>(
        mut self,
        data: &T,
        scorer: &S,
        rng: &mut R,
        max_iter: usize,
    ) -> Option<SaturatingHist<f64>> {
        let mut c = Self::C_START;
        let mut weights =
            SaturatingHist::new(self.inner.min(), self.inner.max(), self.inner.step());

        while c > Self::C_END {
            self.zero();

            // this guarantees we start with a clean slate for every WL factor
            // we could also undo all the modifys, but that is prolly more expensive than just
            // cloning
            let mut data = data.clone();
            let mut score = scorer.score(&mut data);

            let mut flat = false;

            for _ in 0..max_iter {
                score = data.mh_step(scorer, score, rng, Some(&weights));
                *weights.get_mut(score) -= c;
                self.increment(score);

                if self.is_flat() {
                    flat = true;
                    break;
                }
            }

            if !flat {
                return None;
            }

            // guarantee our scores are <= 0
            // this is from the original codebase, the purpose seems to be to guarantee that
            // negating and taking the exponent gives a valid probability.
            let max_weight = *weights
                .values()
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap();
            for weight in weights.values_mut() {
                *weight -= max_weight;
            }

            c /= 2.0;
        }

        Some(weights)
    }
}

/// An estimator for the time-average variance constant using a trapezoidal selection rule.
///
/// This estimator is from [1].
///
/// By the Markov Chain central limit theorem the difference between the sample mean and the true
/// mean (modulo a factor based on number of samples) will approach a normal distribution with `0`
/// mean and some variance. This variance is called the time-average variance constant (TAVC) [2].
///
/// We use the TAVC to determine a confidence interval for our stopping condition in MS-DPR [3].
/// See [`TrajectoryInfo`] for more details on how this plays into stopping conditions.
///
/// # The Algorithm
///
/// ## Characteristics
/// This is a batch means estimator for TAVC, meaning that observations are broken up in to
/// batches. Specifically this method uses overlapping batches that are shaped like trapezoids
/// (hence the name trapezoidal selection rule for batches). Additionally this is a recursive
/// estimator, meaning that we *do not* need to store the whole trajectory to estimate TAVC. This
/// is very beneficial for long running processes that can have millions of observations.
///
/// ## Estimation Process
/// The final estimator is the sum of squared sums of observations in each batch divided by the
/// sums of all batch sizes. This, however, only works in the case of zero mean observation sets.
/// We use correction for non-zero means shown on page 612 in [1].
///
/// The recursive process stores two sums for the current batch: one total and one offset. At each
/// step when we move to a new step we remove the offset from the batch sum to get the new batch
/// sum. Details are shown in Section 3.1 of [1].
///
/// The last remaining aspect is to select when to start a new batch. This is done with two
/// sequences. First, `a` marks which observations split batches into triangles.
/// ```txt
///        a2
///        |---
///        ****
///        ***
///        **
///        *
///    ****
///    ***
///    **
///    *
///    |---
///    a1
/// ```
/// Second, `b` marks how many samples to keep as a rectangle to the left of the triangle.
/// ```txt
/// b=3 a
///     |---
/// *** ****
/// *** ***
/// *** **
/// *** *
/// ```
/// We use the suggested definition of `a_1 = 1` and `a_k = ⌊k^(1.5)⌋`. Additionally, we set
/// `b_k = a_k - a_{k-1}`.
///
/// # Using the Estimator
/// The proper initial state will be set with the [`Default`] implementation for this type. Then,
/// as observations appear, they can be accounted for with [`TavcEstimator::update`]. The current
/// sample mean can be retrieved with [`TavcEstimator::mean`] and the estimate of the TAVC can be
/// retrieved with [`TavcEstimator::variance`].
///
///
/// # Relevant Papers
/// 1. [The method we use](https://doi.org/10.1007/s11222-015-9548-7)
/// 2. [Its predecessor](https://doi.org/10.1214/08-AAP587)
/// 3. [Updated MS-DPR](https://doi.org/10.4230/LIPIcs.WABI.2017.14)
#[derive(Debug)]
struct TavcEstimator {
    /// The number of samples we've observed so far.
    ///
    /// Each time we call [`TavcEstimator::update`] we increment this count by one. This is used in
    /// conjunction with [`TavcEstimator::running_sum`] to compute the mean. This variable is `n`
    /// in the original paper and codebase.
    num_obs: usize,
    /// The total value of all observations seen so far.
    ///
    /// Each time we call [`TavcEstimator::update`] the observation parameter is added to our
    /// running sum for use in computation of the mean. Note, this isn't how it's descried in the
    /// paper but we're trading fewer operations for slightly more memory used. This variable is
    /// `X_sum` in the original codebase.
    running_sum: f64,
    /// The sample mean of the observations.
    ///
    /// This is computed with [`TavcEstimator::running_sum`] and [`TavcEstimator::num_obs`] and
    /// cached her to save some operations at the cost of memory. This variable is X bar in the
    /// original text and `mu` in the old codebase.
    mean: f64,
    /// The index of the current batch we are on.
    ///
    /// This is used to compute the boundaries of the `a` and `b` series, so it needs to be stored.
    /// In the paper and the original codebase this variable is `k`.
    batch_idx: usize,
    /// The sum of observations in the current batch.
    ///
    /// This is batch [`TavcEstimator::batch_idx`]. In the original paper this variable is `W`. It
    /// is computed by adding in new observations until we spill over into a new batch, in which
    /// case we correct for spurious overlap with the offset `W_del`.
    batch_sum: f64,
    /// An offset to prevent spurious overlap when we move over to a new batch.
    ///
    /// This is relative to the batch [`TavcEstimator::batch_idx`]. In the original paper this
    /// variable is called `W_del`. It is reset when we are in the rectangular region of the
    /// trapezoid.
    batch_sum_offset: f64,
    /// The numerator of the variance estimation when zero mean.
    ///
    /// This is described in the struct-level docs. It is the sum of the squares of all the batch
    /// sums seen so far. While this is named as the numerator of the estimation, it isn't quite
    /// complete because it doesn't account for the non-zero mean correction. In the original paper
    /// this variable is `V`.
    est_numer: f64,
    /// The denominator of the variance estimation when zero mean.
    ///
    /// This is described in the struct-level docs. It is the sum of all the batch sizes seen so
    /// far. While this is named as the denominator of the estimation, it isn't quite complete
    /// because it doesn't account for the non-zero mean correction. In the original paper this
    /// variable is `v`.
    est_denom: f64,
    /// The size of the current batch.
    ///
    /// This relates to the batch [`TavcEstimator::batch_idx`]. The original paper calls this value
    /// `B`.
    current_batch_size: usize,
    /// The sample that will be the start of the next batch.
    ///
    /// This is computed using the trapezoidal selection rule in [`TavcEstimator::a_tsr`],
    /// described in the struct-level docs. The original paper calls this value `a_new`.
    next_batch_start: usize,
    /// The width of the rectangular part of the trapezoid for the next batch.
    ///
    /// This is computed using the trapezoidal selection rule in [`TavcEstimator::b_tsr`]. The
    /// original paper calls this value `b_new`.
    next_additional_samples: usize,
    /// The order one coefficient for the non-zero mean correction.
    ///
    /// This is the sum of all batches weighted by the sizes of the batches. It is called `H` in
    /// the original paper.
    batch_size_weighted_sum: f64,
    /// The order two coefficient for the non-zero mean correction.
    ///
    /// This is the sum of the squares of all the batch sizes. It is called `h` in the original
    /// paper.
    squared_batch_size_sum: f64,
    /// An estimate of the posterior variance without normalization.
    ///
    /// This is computed by just summing the difference between each new observation and the
    /// running sample mean. This is used in the stopping condition for [`super::MsDpr`].
    lambda_sum: f64,
}

impl Default for TavcEstimator {
    fn default() -> Self {
        Self {
            num_obs: 0,
            running_sum: 0.0,
            mean: 0.0,
            batch_idx: 1,
            batch_sum: 0.0,
            batch_sum_offset: 0.0,
            est_numer: 0.0,
            est_denom: 0.0,
            current_batch_size: 0,
            next_batch_start: 1,
            next_additional_samples: 0,
            batch_size_weighted_sum: 0.0,
            squared_batch_size_sum: 0.0,
            lambda_sum: 0.0,
        }
    }
}

impl TavcEstimator {
    pub fn update(&mut self, new_obs: f64) {
        // simple addition for running sums, also updating mean
        self.num_obs += 1;
        self.running_sum += new_obs;
        self.batch_sum += new_obs;
        self.mean = self.running_sum / self.num_obs as f64;
        self.lambda_sum += (new_obs - self.mean) * (new_obs - self.mean);

        if self.num_obs == self.next_batch_start {
            // we're at a batch boundary now
            // increment batch counter, remove offset from sum to avoid overcounting
            // we are also now just the rectangle part of the trapezoid plus one more sample
            // finally, we need to compute the new next boundaries
            self.batch_idx += 1;
            self.batch_sum -= self.batch_sum_offset;
            self.current_batch_size = self.next_additional_samples + 1;
            self.next_batch_start = Self::a_tsr(self.batch_idx);
            self.next_additional_samples = Self::b_tsr(self.batch_idx);
        } else {
            // same batch, track the batch size
            self.current_batch_size += 1;
        }

        // update the batch sum correction offset
        // next_real_start is the first sample in the rectangular part of the trapezoid.
        // let next_real_start = self.next_batch_start - self.next_additional_samples;
        let next_real_start = self
            .next_batch_start
            .checked_sub(self.next_additional_samples)
            .unwrap();
        if self.num_obs == next_real_start - 1 {
            self.batch_sum_offset = self.batch_sum;
        } else if self.num_obs == next_real_start {
            self.batch_sum_offset = self.batch_sum - new_obs;
        }

        // finally compute necessary parts for the final variance estimation
        self.est_numer += self.batch_sum * self.batch_sum;
        self.batch_size_weighted_sum += self.current_batch_size as f64 * self.batch_sum;
        self.est_denom += self.current_batch_size as f64;
        self.squared_batch_size_sum += (self.current_batch_size * self.current_batch_size) as f64;
    }

    /// Get the first sample number for the triangular part of the trapezoidal batch.
    fn a_tsr(batch_num: usize) -> usize {
        if batch_num == 1 {
            return 1;
        }

        let out = (batch_num as f64).powf(1.5).floor() as usize;
        out.max(batch_num)
    }

    /// Get the width of the rectangular part of the trapezoidal batch.
    fn b_tsr(batch_num: usize) -> usize {
        Self::a_tsr(batch_num).saturating_sub(Self::a_tsr(batch_num - 1) + 1)
    }

    /// Retrieve the current sample mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Retrieve the current estimate for the TAVC.
    pub fn variance(&self) -> f64 {
        let mean = self.mean();
        (self.est_numer - 2.0 * self.batch_size_weighted_sum * mean
            + self.squared_batch_size_sum * mean * mean)
            / self.est_denom
    }

    /// An estimate of the posterior variance.
    pub fn lambda(&self) -> f64 {
        self.lambda_sum / self.num_obs as f64
    }
}

/// A tracker for the trajectory of a [`super::MsDpr`] simulation.
///
/// This tracks observations and uses a [`TavcEstimator`] to decide when to terminate simulation.
/// This is done by comparing the confidence interval that can be computed from the TAVC to the
/// sample variance scaled by some user-specified ε.
///
/// # Relevant Papers
/// - [Updated MS-DPR](https://doi.org/10.4230/LIPIcs.WABI.2017.14)
/// - [Stopping Criterion](https://www.jstor.org/stable/24311039)
#[derive(Debug)]
pub struct TrajectoryInfo<'a> {
    /// Log score weights computed using Wang-Landau sampling.
    weights: &'a SaturatingHist<f64>,
    /// The minimum score to count a sample as having a better score than our target.
    threshold: f64,
    /// A normalizing constant for score weights.
    ///
    /// Score weights are an estimate of the inverse log probability of seeing that score. To get
    /// to a final probability we need to normalize these, since use in Metropolis-Hastings does
    /// not require normalization. This is the sum of the negative exponentiation of each weight.
    normalizer: f64,
    /// The estimator for TAVC.
    tavc: TavcEstimator,
}

impl<'a> TrajectoryInfo<'a> {
    /// Construct a new trajectory information tracker with the provided log score weights and
    /// score threshold.
    pub fn new(weights: &'a SaturatingHist<f64>, threshold: f64) -> Self {
        let normalizer = weights
            .values()
            .map(|log_prob| (-log_prob).exp())
            .sum::<f64>();
        Self {
            weights,
            threshold,
            normalizer,
            tavc: TavcEstimator::default(),
        }
    }

    /// Add a new observation to this trajectory tracker.
    pub fn add(&mut self, value: f64) {
        if value >= self.threshold {
            let weight = self.weights.get(value);
            self.tavc.update((-weight).exp() / self.normalizer);
        } else {
            self.tavc.update(0.0);
        }
    }

    /// Get the current sample mean.
    pub fn mean(&self) -> f64 {
        self.tavc.mean()
    }

    /// Determine if the sampling should be terminated.
    pub fn stop_iteration(&self, target_eps: f64, threshold: f64) -> bool {
        let se = (self.tavc.variance() / self.tavc.num_obs as f64).sqrt();
        // double the limit of the confidence interval for 95% of std normal
        let w = 2.0 * 1.96 * se;

        if se == 0.0 || self.tavc.lambda() == 0.0 {
            return false;
        }

        let curr_eps = w / self.tavc.lambda().sqrt();

        if threshold == 0.0 {
            curr_eps <= target_eps
        } else {
            (threshold < self.tavc.mean() - w / 2.0) || (threshold > self.tavc.mean() + w / 2.0)
        }
    }
}
