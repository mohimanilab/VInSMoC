//! This crate defines the data API and contains implementations of machine learning techniques and
//! methods built on this data API.

#![warn(missing_docs, clippy::all)]

use ordered_float::OrderedFloat;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub mod data;
use data::{DataPoint, DataView, Labeled};
pub mod clustering;
pub mod counting;
pub mod trees;

/// Type alias for ordered float comparisons
pub(crate) type Of64 = OrderedFloat<f64>;

/// Trait for calculating distance, used in models and applications such as clustering.
pub trait Distance: Sized {
    /// Calculate the distance between `self` and `other`.
    ///
    /// This should be a symmetric distance measure, i.e. distance from `self` to `other` should be
    /// the same as the distance between `other` and `self`.
    fn distance(&self, other: &Self) -> Of64;
}

/// Trait for calculating a centroid which represents an aggregate across a set of points.
pub trait Centroid: Sized {
    /// Calculate a centroid instance of [`Self`] by aggregating the input points.
    fn centroid<'a, I>(points: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
        Self: 'a;
}

/// Trait for a classification model.
pub trait Model<L: data::CategoricalFeature>: Sync + Send + Sized {
    /// Each model should have an associated set of hyperparameters.
    ///
    /// The hyperparameters are required to implement [`Default`].
    type HyperParams: Default;

    /// Trains on a [`DataView`] and returns a fitted [`Model`] according to the
    /// specified hyperparameters.
    fn train<D>(dataview: DataView<D, L>, hyperparameters: &Self::HyperParams) -> Self
    where
        D: DataPoint<usize> + Labeled<L>;

    /// Trains on a [`DataView`] and returns a fitted [`Model`] using default
    /// hyperparameters.
    fn train_default<D>(dataview: DataView<D, L>) -> Self
    where
        D: DataPoint<usize> + Labeled<L>,
    {
        Self::train(dataview, &Self::HyperParams::default())
    }

    /// Function to classify a [`DataPoint`].
    ///
    /// The return value the label for the dataset.
    fn classify<D: DataPoint<usize>>(&self, point: &D) -> L;

    /// Default implementation to classify a list of [`DataPoint`]'s, given a slice of points.
    fn classify_all<D: DataPoint<usize>>(&self, points: &[D]) -> Vec<L> {
        points.iter().map(|point| self.classify(point)).collect()
    }

    /// Default implementation to classify a list of [`DataPoint`]'s, given anything which
    /// implements [`IntoIterator`].
    fn classify_iter<'a, D: 'a + DataPoint<usize>, I>(&'a self, points: I) -> Vec<L>
    where
        I: IntoIterator<Item = &'a D>,
    {
        points
            .into_iter()
            .map(|point| self.classify(point))
            .collect()
    }

    /// Default implementation to classify a list of [`DataPoint`]'s in parallel.
    ///
    /// The recommended usage of this is for cases when we have a huge list of datapoints.
    fn classify_all_parallel<D: DataPoint<usize> + Sync + Send>(&self, points: &[D]) -> Vec<L> {
        points
            .into_par_iter()
            .map(|point| self.classify(point))
            .collect()
    }
}
