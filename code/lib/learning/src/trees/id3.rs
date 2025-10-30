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
    pub data: DataView<'db, D, L>,
}

impl<'db, D: DataPoint<usize> + Labeled<L>, L: CategoricalFeature> TrainingNodeContext<'db, D, L> {
    pub fn root(dataview: DataView<'db, D, L>) -> Self {
        let num_features = dataview.dataset().num_features();
        let mut available_features = FixedBitSet::with_capacity(num_features);
        available_features.toggle_range(..);
        Self {
            available_features,
            data: dataview,
        }
    }

    pub fn get_prob_dist(&self) -> FxHashMap<usize, Of64> {
        // data: &DataView<'db, D, L>) -> FxHashMap<usize, Of64> {
        let tot = self.data.len() as f64;
        self.data
            .group_by_label()
            .iter()
            .map(|(&k, v)| (k, Of64::from(v.len() as f64 / tot)))
            .collect()
    }

    pub fn select_feature<'a>(
        &'a self,
        split_type: SplitQuality,
    ) -> (usize, FxHashMap<usize, DataView<'db, D, L>>) {
        self.available_features
            .ones()
            .map(|feature| (feature, self.data.group_by_feature(feature)))
            .max_by_key(|(_, grouped)| {
                // Each grouped is a potential data split by a feature
                match split_type {
                    SplitQuality::Entropy => Of64::from(
                        grouped
                            .values()
                            .map(|x| x.len() as f64 * x.entropy())
                            .sum::<f64>(),
                    ),
                    SplitQuality::Impurity => Of64::from(
                        grouped
                            .values()
                            .map(|x| x.len() as f64 * x.impurity())
                            .sum::<f64>(),
                    ),
                }
                // The value that is returned is the feature which leads to the best data split,
                // and then a hashmap of feature value to the dataview it corresponds to.
            })
            .unwrap()
    }
}

/// Map between node and training context.
///
/// Each node acts as a split on our data. Therefore we maintain a map of the node index, and the
/// training node context it corresponds to. This way we can quickly get the probability
/// distribution of a node.
pub struct TrainingContext<'db, D, L> {
    pub data: FxHashMap<NodeIndex, TrainingNodeContext<'db, D, L>>,
}
