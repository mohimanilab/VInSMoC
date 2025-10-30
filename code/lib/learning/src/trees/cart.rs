use super::SplitQuality;
use fixedbitset::FixedBitSet;
use petgraph::graph::NodeIndex;
use rustc_hash::FxHashMap;

use crate::{
    data::{CategoricalFeature, DataPoint, DataView, Labeled},
    Of64,
};

/// Training Node Context is a mask over the training data points and the available features.
pub struct TrainingNodeContext<'db, D, L> {
    pub available_features: FixedBitSet,
    pub depth: usize,
    pub data: DataView<'db, D, L>,
}

pub(crate) type BestSplit<'db, D, L> = (usize, f64, DataView<'db, D, L>, DataView<'db, D, L>);

impl<'db, D: DataPoint<usize> + Labeled<L>, L: CategoricalFeature> TrainingNodeContext<'db, D, L> {
    pub fn root(dataview: DataView<'db, D, L>) -> Self {
        let num_features = dataview.dataset().num_features();
        let mut available_features = FixedBitSet::with_capacity(num_features);
        available_features.toggle_range(..);
        Self {
            available_features,
            depth: 0,
            data: dataview,
        }
    }

    /// Returns (Option(feature_idx, split_value, score, left, right))
    pub fn select_feature<'a>(
        &'a mut self,
        split_criterion: SplitQuality,
    ) -> (usize, Option<BestSplit<'db, D, L>>) {
        // all features should be available, on each round
        // if bagging, then sample a subset of the features
        self.available_features
            .ones()
            .map(|feature_idx| {
                (
                    feature_idx,
                    self.data.partition_by_feature(feature_idx, split_criterion),
                )
            })
            .min_by_key(|(_feature, grouped)| match grouped.as_ref() {
                Some((_, score, _, _)) => Of64::from(*score),
                None => Of64::from(f64::MAX),
            })
            .unwrap()
    }
}

/// Map between decision tree node and the corresponding training context.
pub struct TrainingContext<'db, D, L> {
    pub data: FxHashMap<NodeIndex, TrainingNodeContext<'db, D, L>>,
}
