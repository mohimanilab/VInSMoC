//! Data handling for machine learning applications.

use std::marker::PhantomData;

use rand::{prelude::SliceRandom, Rng};
use rustc_hash::{FxHashMap, FxHashSet};
use uuid::Uuid;

use crate::trees::SplitQuality;

mod categorical;
pub use categorical::CategoricalFeature;
mod frozen;
pub use frozen::{FeatureSet, FeatureUnit};
mod numeric;
pub use numeric::{NumericFeature, NumericValue, Precision};
mod onehot;
pub use onehot::OneHot;

/// Trait to denote that a datapoint is labeled, i.e. the implementor has an associated label for
/// the purpose of classification tasks.
pub trait Labeled<L: categorical::CategoricalFeature> {
    /// Retrieve the label of the datapoint.
    fn label(&self) -> L;
}

/// A trait to mark a type as a data point consisting of a feature set and a unique identifier.
pub trait DataPoint<T: frozen::FeatureUnit> {
    /// Retrieve a reference to the feature set of this data point.
    fn features(&self) -> &frozen::FeatureSet<T>;

    /// Retrieve the unique identifier for this data point.
    ///
    /// This function must be consistent, i.e. calling [`DataPoint::id`] multiple times on the same
    /// data point must return the same unique ID.
    fn id(&self) -> Uuid;
}

/// A collection of labeled data points with categorical feature sets.
///
/// This is generic over the data point itself, `D`, and the labels for those data points, `L`.
#[derive(Clone)]
pub struct CategoricalDataset<D, L> {
    /// A map between the unique id for a datapoint and the datapoint itself
    pub points: FxHashMap<Uuid, D>,
    label: PhantomData<L>,
}

impl<L: categorical::CategoricalFeature, D: DataPoint<usize>> FromIterator<D>
    for CategoricalDataset<D, L>
{
    /// Create a new [`CategoricalDataset`] from an iterator of datapoints with type `D`.
    fn from_iter<T: IntoIterator<Item = D>>(datapoints: T) -> Self {
        let points = datapoints
            .into_iter()
            .map(|x| (x.id(), x))
            .collect::<FxHashMap<_, _>>();

        Self {
            points,
            label: PhantomData::<L>,
        }
    }
}

impl<L: categorical::CategoricalFeature, D: DataPoint<usize>> CategoricalDataset<D, L> {
    /// Create a new [`CategoricalDataset`] from a list of datapoints with type `D`.
    pub fn new(datapoints: Vec<D>) -> Self {
        let points = datapoints
            .into_iter()
            .map(|x| (x.id(), x))
            .collect::<FxHashMap<_, _>>();
        Self {
            points,
            label: PhantomData::<L>,
        }
    }

    /// Create a full [`DataView`] for this [`CategoricalDataset`].
    pub fn full(&self) -> DataView<D, L> {
        DataView::full(self)
    }
}

impl<D: DataPoint<usize>, L: categorical::CategoricalFeature> CategoricalDataset<D, L> {
    /// Get the number of data points in this dataset.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Check if this dataset has no data points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Get the number of features a particular data point in this data set has.
    ///
    /// # Note
    /// This assumes that all data points in the data set have the same number of features. It is
    /// possible to manufacture a [`DataPoint`] implementation that does not follow this
    /// restriction. In those cases this method's output can't be trusted.
    ///
    /// # Panics
    /// If this dataset has no points this will panic.
    pub fn num_features(&self) -> usize {
        self.points
            .values()
            .next()
            .expect("dataset had no points")
            .features()
            .len()
    }

    /// Get the number of possible distinct labels for points in this dataset.
    ///
    /// This is based off of the [`CategoricalFeature`] implementation of the label. If the dataset
    /// does not use all possible labels this will still return the number of *possible* labels for
    /// a dataset over points of this type.
    pub fn num_labels(&self) -> usize {
        L::NUM_CATEGORIES
    }

    /// Construct views into the dataset grouped by distinct values taken on by the feature at
    /// `feature_idx`.
    ///
    /// The returned map is keyed by observed values of the feature at `feature_idx`. If a possible
    /// value of that feature is not observed it will not have a key in the returned map, meaning
    /// that all values in the returned map have at least one point in them.
    pub fn group_by_feature(&self, feature_idx: usize) -> FxHashMap<usize, DataView<D, L>> {
        let mut groups = FxHashMap::default();
        groups.reserve(self.num_features());

        for (id, point) in &self.points {
            let view = groups
                .entry(point.features()[feature_idx])
                .or_insert_with(|| DataView::empty(self));
            view.add_point(*id);
        }
        groups
    }

    /// Construct views into the dataset grouped by distinct label values.
    ///
    /// The returned map is keyed by observed values of the labels. If a possible label value is
    /// never observed in this dataset it will not appear as a key in the returned map, meaning
    /// that all values in the returned map have at least one point in them.
    pub fn group_by_label(&self) -> FxHashMap<usize, DataView<D, L>>
    where
        D: Labeled<L>,
    {
        let mut groups = FxHashMap::default();
        groups.reserve(self.num_labels());
        for (id, point) in &self.points {
            let view = groups
                .entry(point.label().index())
                .or_insert_with(|| DataView::empty(self));
            view.add_point(*id);
        }
        groups
    }

    /// Split this dataset into training and testing views.
    ///
    /// Roughly `train_percent` of the points will be in the training set, with the remainder in
    /// the testing set. The first view returned is the training set and the second view is the
    /// testing set.
    pub fn train_test_split<'a, R: Rng + ?Sized>(
        &'a self,
        train_percent: f64,
        rng: &mut R,
    ) -> (DataView<'a, D, L>, DataView<'a, D, L>) {
        let mut train = DataView::empty(self);
        let mut test = DataView::empty(self);

        for id in self.points.keys().copied() {
            if rng.gen_bool(train_percent) {
                train.add_point(id);
            } else {
                test.add_point(id);
            }
        }
        (train, test)
    }

    /// Sample with replacement from the dataset.
    pub fn bootstrap<'a, R: Rng + ?Sized>(
        &'a self,
        num_samples: usize,
        rng: &mut R,
    ) -> DataView<'a, D, L> {
        let mut bootstrap_sample = DataView::empty(self);
        let ids = self.points.keys().copied().collect::<Vec<_>>();

        for _ in 0..num_samples {
            let uuid = ids.choose(rng).expect("dataset had no points");
            bootstrap_sample.add_point(*uuid);
        }

        bootstrap_sample
    }
}

/// Mask over a [`CategoricalDataset`]
pub struct DataView<'db, D, L> {
    /// Reference to the dataset
    dataset: &'db CategoricalDataset<D, L>,
    /// Vec of keys contained in this view.
    keys: FxHashSet<Uuid>,
}

impl<'db, D: DataPoint<usize>, L: categorical::CategoricalFeature> DataView<'db, D, L> {
    /// Construct an empty view into the provided `dataset`.
    pub fn empty(dataset: &'db CategoricalDataset<D, L>) -> Self {
        Self {
            dataset,
            keys: FxHashSet::default(),
        }
    }

    /// Construct a view into the provided `dataset` containing all of its points.
    pub fn full(dataset: &'db CategoricalDataset<D, L>) -> Self {
        Self {
            dataset,
            keys: dataset.points.keys().copied().collect(),
        }
    }

    /// Construct a view into the provided `dataset` containing points with IDs in `keys`.
    pub fn new<I: IntoIterator<Item = Uuid>>(
        dataset: &'db CategoricalDataset<D, L>,
        keys: I,
    ) -> Self {
        Self {
            dataset,
            keys: keys.into_iter().collect(),
        }
    }

    /// Get a list of the data points in this view.
    pub fn datapoints(&self) -> Vec<&D> {
        self.dataset
            .points
            .iter()
            .filter(|(uuid, _point)| self.keys.contains(uuid))
            .map(|(_x, y)| y)
            .collect()
    }

    /// Get a reference to the underlying dataset
    pub fn dataset(&self) -> &'db CategoricalDataset<D, L> {
        self.dataset
    }

    /// Get the number of features for each datapoint
    pub fn num_features(&self) -> usize {
        self.dataset.num_features()
    }

    /// Get the number of data points in this view.
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Check if the data view is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Add a point with the given `id` to this view.
    pub fn add_point(&mut self, id: Uuid) {
        self.keys.insert(id);
    }

    /// Returns the best value split for a given feature.
    ///
    /// The best value split is the split into all datapoints with that feature
    /// value and all datapoints without that feature value with best criterion.
    pub(crate) fn partition_by_feature(
        &self,
        feature_idx: usize,
        criterion: SplitQuality,
    ) -> Option<(usize, f64, Self, Self)>
    where
        D: Labeled<L>,
    {
        // If there is only one point, or all of the same point, there is no more splitting to be
        // done
        if self.len() == 1 || self.entropy() == 0.0 {
            return None;
        }

        let mut best_split = None;
        let mut best_score = f64::MAX;

        // For each feature value that appears
        for (value, points_with_value) in self.group_by_feature(feature_idx) {
            // Split into points with given feature value, and points without that feature value
            let points_without_value = Self {
                dataset: self.dataset,
                keys: &self.keys - &points_with_value.keys,
            };
            // Calculate weighted score based on specified criteria
            let score = match criterion {
                SplitQuality::Entropy => {
                    points_with_value.len() as f64 * points_with_value.entropy()
                        + points_without_value.len() as f64 * points_without_value.entropy()
                }
                SplitQuality::Impurity => {
                    points_with_value.len() as f64 * points_with_value.impurity()
                        + points_without_value.len() as f64 * points_without_value.impurity()
                }
            };
            if score < best_score {
                best_score = score;
                best_split = Some((value, score, points_with_value, points_without_value));
            }
        }

        best_split
    }

    /// Group points in this view by their values at feature `feature_idx`.
    ///
    /// See [`CategoricalDataset::group_by_feature`].
    pub fn group_by_feature(&self, feature_idx: usize) -> FxHashMap<usize, Self> {
        let mut groups = FxHashMap::default();
        groups.reserve(self.dataset.num_features());

        for (id, point) in &self.dataset.points {
            if !self.keys.contains(id) {
                continue;
            }
            let view = groups
                .entry(point.features()[feature_idx])
                .or_insert_with(|| DataView::empty(self.dataset));
            view.add_point(*id);
        }
        groups.shrink_to_fit();
        groups
    }

    /// Group points in this view by their labels.
    ///
    /// See [`CategoricalDataset::group_by_label`].
    pub fn group_by_label(&self) -> FxHashMap<usize, Self>
    where
        D: Labeled<L>,
    {
        let mut groups = FxHashMap::default();
        groups.reserve(self.dataset.num_labels());
        for (id, point) in &self.dataset.points {
            if !self.keys.contains(id) {
                continue;
            }
            let view = groups
                .entry(point.label().index())
                .or_insert_with(|| DataView::empty(self.dataset));
            view.add_point(*id);
        }
        groups.shrink_to_fit();
        groups
    }

    /// Compute the entropy over labels for this data view.
    ///
    /// See [definition of entropy](https://en.wikipedia.org/wiki/Entropy_(information_theory)).
    /// This implementation using base 2 logarithms.
    pub fn entropy(&self) -> f64
    where
        D: Labeled<L>,
    {
        let total = self.len() as f64;
        self.group_by_label()
            .values()
            .map(|x| x.len())
            .map(|len| {
                let p = len as f64 / total;
                -p * p.log2()
            })
            .sum()
    }

    /// Compute the Gini impurity over labels for this data view.
    ///
    /// See [definition of
    /// impurity](https://en.wikipedia.org/wiki/Decision_tree_learning#Gini_impurity).
    pub fn impurity(&self) -> f64
    where
        D: Labeled<L>,
    {
        let total = self.len() as f64;
        let one_minus_impurity = self
            .group_by_label()
            .values()
            .map(|x| x.len())
            .map(|len| {
                let p = len as f64 / total;
                p * p
            })
            .sum::<f64>();
        1.0f64 - one_minus_impurity
    }
}
