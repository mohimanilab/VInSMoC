//! Implementations of various Decision Tree and Random Forest machine learning methods.

use std::{hash::Hash, marker::PhantomData};

use crate::{Labeled, Model};
use fixedbitset::FixedBitSet;
use petgraph::{
    graph::{DiGraph, NodeIndex},
    visit::EdgeRef,
};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

use crate::{
    data::{CategoricalFeature, DataPoint, DataView},
    Of64,
};

mod cart;
mod id3;

/// Algorithm used to train the decision tree
#[non_exhaustive]
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum TrainingAlgorithm {
    /// The Iterative Dichotomiser 3 (ID3) algorithm for training decision trees.
    ///
    /// This algorithm at each iteration greedily selects the feature index whose entropy is
    /// minimized. This feature index is now forbidden from splitting again in that rooted subtree.
    /// It continues until either there is only one label left or no feature remain.
    ///
    /// # Relevant Papers
    /// See the [original paper](https://doi.org/10.1007/BF00116251).
    Id3,
    /// The Classification And Regression Trees (CART) algorithm for training decision trees.
    ///
    /// 1. Find each feature’s best split. For each feature with K different values there exist
    ///    K-1 possible splits. Find the split which optimizes the splitting criterion. The
    ///    resulting set of splits contains the best splits (one for each feature).
    ///
    /// 2. Find the node's best split across all features. Among the best splits from Step i, find
    ///    the best one.
    ///
    /// 3. Split the node using best node split from Step ii and repeat from Step i until stopping
    ///    criterion is satisfied
    ///
    /// # Relevant Papers
    /// See the [original book](https://doi.org/10.1201/9781315139470).
    Cart,
}

/// An enum to describe the different ways to judge how good a split is.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum SplitQuality {
    /// Uses entropy of the labels, see
    /// [here](https://en.wikipedia.org/wiki/Entropy_(information_theory)) for more details.
    Entropy,
    /// Uses Gini impurity of the labels, see
    /// [here](https://en.wikipedia.org/wiki/Decision_tree_learning#Gini_impurity) for more
    /// details.
    Impurity,
}

/// Structure for the RandomForest hyperparameters.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct RandomForestHyperparameters {
    /// Number of [`DecisionTree`] classifiers to use.
    pub num_classifiers: usize,
    /// Hyperparameters which apply to all [`DecisionTree`]'s.
    pub tree_hyperparameters: DecisionTreeHyperparameters,
}

impl Default for RandomForestHyperparameters {
    /// Default values for RandomForest hyperparameters.
    ///
    /// The default value of `num_classifiers` is `100usize`. For the default value of the tree
    /// hyperparameters, see [`DecisionTreeHyperparameters::default`].
    fn default() -> Self {
        Self {
            num_classifiers: 100,
            tree_hyperparameters: DecisionTreeHyperparameters::default(),
        }
    }
}

/// Implementation of RandomForest, which is a collection of [`DecisionTree`]'s.
///
/// Each [`DecisionTree`] is trained on a bootstrapped sample of the training data. Each
/// bootstrapped sample is the same size as the original training data.
#[derive(Serialize, Deserialize, Debug)]
pub struct RandomForest<L> {
    classifiers: Vec<DecisionTree<L>>,
}

impl<L: CategoricalFeature + Hash + Eq + Copy + Sync + Send> RandomForest<L> {
    /// Parallel implementation of training.
    ///
    /// This is a trivial parallelized implementation of random forest training. Since each tree is
    /// trained on a bootstrapped dataset, training can be done completely independently for each
    /// classifier.
    pub fn train_parallel<D>(
        dataview: DataView<D, L>,
        hyperparameters: &RandomForestHyperparameters,
    ) -> Self
    where
        D: DataPoint<usize> + Sync + Send + Labeled<L>,
    {
        let mut classifiers = vec![];
        (0..hyperparameters.num_classifiers)
            .into_par_iter()
            .map(|_x| {
                let dataset = dataview.dataset();
                let mut rng = rand::thread_rng();
                let clf = DecisionTree::train(
                    dataset.bootstrap(dataset.len(), &mut rng),
                    &hyperparameters.tree_hyperparameters,
                );
                clf
            })
            .collect_into_vec(&mut classifiers);

        RandomForest { classifiers }
    }
}

impl<L: CategoricalFeature + Hash + Eq> Model<L> for RandomForest<L> {
    type HyperParams = RandomForestHyperparameters;

    /// Perform training, given the training data and a set of hyperparameters.
    fn train<D>(dataview: DataView<D, L>, hyperparameters: &Self::HyperParams) -> Self
    where
        D: DataPoint<usize> + Labeled<L>,
    {
        let dataset = dataview.dataset();
        let classifiers = (0..hyperparameters.num_classifiers)
            .into_iter()
            .map(|_i| {
                let mut rng = rand::thread_rng();
                let clf = DecisionTree::train(
                    dataset.bootstrap(dataset.len(), &mut rng),
                    &hyperparameters.tree_hyperparameters,
                );
                clf
            })
            .collect::<Vec<_>>();

        Self { classifiers }
    }

    /// Perform classification, given an input [`DataPoint`].
    fn classify<D: DataPoint<usize>>(&self, point: &D) -> L {
        let mut votes = FxHashMap::default();
        for clf in self.classifiers.iter() {
            let prediction = clf.classify(point);
            let entry = votes.entry(prediction).or_insert(0usize);
            *entry += 1;
        }

        votes
            .iter()
            .max_by(|a, b| a.1.cmp(b.1))
            .map(|(k, _v)| *k)
            .unwrap()
    }
}

impl<L: CategoricalFeature + Hash + Eq + Copy> RandomForest<L> {
    /// Number of classifiers used for RandomForest
    pub fn num_classifiers(&self) -> usize {
        self.classifiers.len()
    }

    // Helper function to return the label distribution of the inner or leaf node that `point` gets
    // mapped to.
    fn distribution<D: DataPoint<usize>>(&self, point: &D) -> FxHashMap<usize, Of64> {
        let mut prediction_freqs = FxHashMap::default();
        for clf in self.classifiers.iter() {
            // Unwrap is infallible, since cart always goes to a leaf node, and id3 will keep the
            // hashmap populated on each inner node
            for (k, v) in clf.distribution(point).unwrap().iter() {
                let entry = prediction_freqs
                    .entry(*k)
                    .or_insert_with(|| Of64::from(0.0f64));
                *entry += *v
            }
        }

        prediction_freqs
    }

    /// Return the top k predictions for a given datapoint.
    ///
    /// This function should be used if ranked predictions are required.
    pub fn predict_top_k<D: DataPoint<usize>>(&self, point: &D, k: usize) -> Vec<L> {
        let prediction_freqs = self.distribution(point);
        let mut keys = prediction_freqs.keys().collect::<Vec<_>>();
        keys.sort_by(|a, b| {
            prediction_freqs[a]
                .partial_cmp(&prediction_freqs[b])
                .unwrap()
        });
        let mut ordered_labels = keys
            .iter()
            .map(|k| L::from_index(**k).unwrap())
            .collect::<Vec<_>>();
        ordered_labels.truncate(k);
        ordered_labels
    }
}

/// Decision tree structure, which consists of a graph of decision nodes and associated
/// hyperparameters.
#[derive(Debug, Deserialize, Serialize)]
pub struct DecisionTree<L> {
    tree: DiGraph<Node, usize>,
    algorithm: TrainingAlgorithm,
    label: PhantomData<L>,
}

/// Structure to represent hyperparameters when building a decision tree.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct DecisionTreeHyperparameters {
    /// The maximal depth for a tree, beyond which further splits are not considered.
    pub max_depth: usize,
    /// The [`SplitQuality`] rule by which the best split is determined.
    pub criterion: SplitQuality,
    /// The [`TrainingAlgorithm`] used to train the decision tree.
    pub algorithm: TrainingAlgorithm,
}

impl Default for DecisionTreeHyperparameters {
    /// Return default hyperparameters for the decision tree.
    ///
    /// Default values are as follows: `max_depth` is `20usize`, `criterion` is
    /// [`SplitQuality::Entropy`], and `algorithm is [`TrainingAlgorithm::Cart`].
    fn default() -> Self {
        Self {
            max_depth: 20usize,
            criterion: SplitQuality::Entropy,
            algorithm: TrainingAlgorithm::Cart,
        }
    }
}

enum Binary {
    Left,
    Right,
}

#[derive(Debug, Deserialize, Serialize)]
enum Node {
    /// Variant to denote an inner node in the decision tree, and the various attributes it needs
    /// to store.
    Inner {
        feature_idx: usize,
        distribution: Option<FxHashMap<usize, Of64>>,
        feature_split_val: Option<usize>,
    },
    /// Variant to denote a leaf node in the decision tree, and the various attributes it needs
    /// to store.
    Leaf {
        label: usize,
        distribution: FxHashMap<usize, Of64>,
    },
}

impl<L: CategoricalFeature> Model<L> for DecisionTree<L> {
    type HyperParams = DecisionTreeHyperparameters;

    /// Perform training, given a the training data and a set of hyperparameters.
    fn train<D>(dataview: DataView<D, L>, hyperparameters: &Self::HyperParams) -> Self
    where
        D: DataPoint<usize> + Labeled<L>,
    {
        let mut tree = Self::new(hyperparameters);
        match hyperparameters.algorithm {
            TrainingAlgorithm::Id3 => {
                tree.train_id3(dataview, hyperparameters);
            }
            TrainingAlgorithm::Cart => {
                tree.train_cart(dataview, hyperparameters);
            }
        }
        tree
    }

    /// Perform classification on an input [`DataPoint`].
    fn classify<D: DataPoint<usize>>(&self, point: &D) -> L {
        let features = point.features();
        let mut curr_node = NodeIndex::new(0);
        match self.algorithm {
            TrainingAlgorithm::Id3 => {
                while let Node::Inner {
                    feature_idx,
                    distribution: _,
                    feature_split_val: _,
                } = self.tree[curr_node]
                {
                    let target_value = features[feature_idx];
                    curr_node = match self
                        .tree
                        .edges(curr_node)
                        .find(|e| *e.weight() == target_value)
                    {
                        Some(x) => x.target(),
                        None => break,
                    }
                }
            }
            TrainingAlgorithm::Cart => {
                while let Node::Inner {
                    feature_idx,
                    distribution: _,
                    feature_split_val,
                } = self.tree[curr_node]
                {
                    // If node is inner, then feature split val should always be Some
                    let feature_split_val = feature_split_val.unwrap();
                    let target_value = features[feature_idx];
                    // If target value is equal to split value for the given feature, go to left
                    // subchild. Otherwise, if it is not equal, go to right subchild.
                    let mut which_child = Binary::Right as usize;
                    if target_value == feature_split_val {
                        which_child = Binary::Left as usize;
                    }

                    curr_node = match self
                        .tree
                        .edges(curr_node)
                        .find(|e| *e.weight() == which_child)
                    {
                        Some(x) => x.target(),
                        None => break,
                    }
                }
            }
        }

        self.tree[curr_node].classify().unwrap()
    }
}

impl<L: CategoricalFeature> DecisionTree<L> {
    fn new(hyperparameters: &DecisionTreeHyperparameters) -> Self {
        Self {
            tree: DiGraph::default(),
            algorithm: hyperparameters.algorithm,
            label: PhantomData,
        }
    }

    fn train_cart<D>(
        &mut self,
        dataview: DataView<D, L>,
        hyperparameters: &DecisionTreeHyperparameters,
    ) where
        D: DataPoint<usize> + Labeled<L>,
    {
        let mut available_features = FixedBitSet::with_capacity(dataview.num_features());
        available_features.toggle_range(..);

        let mut context = cart::TrainingContext {
            data: FxHashMap::default(),
        };
        let root_idx = self.tree.add_node(Node::Inner {
            feature_idx: 0,
            distribution: None,
            feature_split_val: None,
        });

        context
            .data
            .insert(root_idx, cart::TrainingNodeContext::root(dataview));

        let mut to_visit = vec![root_idx];
        while let Some(parent_node) = to_visit.pop() {
            let mut current_context = context.data.remove(&parent_node).unwrap();

            let (feature_idx, best_split) =
                current_context.select_feature(hyperparameters.criterion);
            match best_split {
                // for each value, we look at the datapoints which have that value
                Some((split_val, _, left_split, right_split)) => {
                    // Update the current node to be the best it can be
                    // If split would make depth too much, make this node a leaf and do not recurse
                    if current_context.depth > hyperparameters.max_depth {
                        // current node should be a leaf
                        self.tree[parent_node] = Node::new_leaf(&current_context.data);
                        continue;
                    }

                    self.tree[parent_node] = Node::Inner {
                        feature_idx,
                        distribution: None,
                        feature_split_val: Some(split_val),
                    };

                    // add left split as left child, right split as right child
                    // left child is dummy, right child is dummy
                    // will be updated in a subsequent iteration
                    let left_child = self.tree.add_node(Node::Inner {
                        feature_idx: 0,
                        distribution: None,
                        feature_split_val: None,
                    });

                    context.data.insert(
                        left_child,
                        cart::TrainingNodeContext {
                            available_features: available_features.clone(),
                            data: left_split,
                            depth: current_context.depth + 1,
                        },
                    );

                    to_visit.push(left_child);
                    self.tree
                        .add_edge(parent_node, left_child, Binary::Left as usize);

                    let right_child = self.tree.add_node(Node::Inner {
                        feature_idx: 0,
                        distribution: None,
                        feature_split_val: None,
                    });

                    context.data.insert(
                        right_child,
                        cart::TrainingNodeContext {
                            available_features: available_features.clone(),
                            data: right_split,
                            depth: current_context.depth + 1,
                        },
                    );

                    to_visit.push(right_child);
                    self.tree
                        .add_edge(parent_node, right_child, Binary::Right as usize);
                }
                None => {
                    // There is no way to split this node so as to increase entropy; therefore it
                    // is a leaf node. Don't push anything to stack
                    self.tree[parent_node] = Node::new_leaf(&current_context.data);
                }
            }
        }
    }

    fn train_id3<D>(
        &mut self,
        dataview: DataView<D, L>,
        hyperparameters: &DecisionTreeHyperparameters,
    ) where
        D: DataPoint<usize> + Labeled<L>,
    {
        let dataset_ref = dataview.dataset();
        let mut available_features = FixedBitSet::with_capacity(dataset_ref.num_features());
        available_features.toggle_range(..);

        let mut context = id3::TrainingContext {
            data: FxHashMap::default(),
        };
        let root_idx = self.tree.add_node(Node::Inner {
            feature_idx: 0,
            distribution: Some(FxHashMap::default()),
            feature_split_val: None,
        });
        context
            .data
            .insert(root_idx, id3::TrainingNodeContext::root(dataview));

        let mut to_visit = vec![root_idx];
        while let Some(parent_node) = to_visit.pop() {
            // Remove the current parent node, and select the best
            let current_context = context.data.remove(&parent_node).unwrap();

            // We choose the best split based on the current context
            let (best_feature, grouped) = current_context.select_feature(hyperparameters.criterion);

            let prob_dist = current_context.get_prob_dist();

            // inner node should keep track of its best feature, and its probability distribution
            self.tree[parent_node] = Node::Inner {
                feature_idx: best_feature,
                distribution: Some(prob_dist),
                feature_split_val: None,
            };

            let mut child_features = current_context.available_features.clone();
            child_features.set(best_feature, false);
            let num_avail_features = child_features.count_ones(..);
            let curr_depth = child_features.len() - num_avail_features;
            let no_more_features = num_avail_features == 0;

            // for each value, we look at the datapoints which have that value
            for (value, data) in grouped {
                let leaf = Node::new_leaf(&data);

                // termination conditions:
                // - only one label in data subset
                // - no more feature to explore
                // - max depth has been reached (add to training context)
                if curr_depth >= hyperparameters.max_depth
                    || leaf.as_leaf().unwrap().len() == 1
                    || no_more_features
                {
                    let node = self.tree.add_node(leaf);
                    self.tree.add_edge(parent_node, node, value);
                    continue;
                }

                let node = self.tree.add_node(Node::Inner {
                    feature_idx: 0,
                    distribution: Some(FxHashMap::default()),
                    feature_split_val: None,
                });
                context.data.insert(
                    node,
                    id3::TrainingNodeContext {
                        data,
                        available_features: child_features.clone(),
                    },
                );
                to_visit.push(node);
                self.tree.add_edge(parent_node, node, value);
            }
        }
    }

    fn distribution<D: DataPoint<usize>>(&self, point: &D) -> Option<&FxHashMap<usize, Of64>> {
        let features = point.features();

        let mut curr_node = NodeIndex::new(0);

        while let Node::Inner {
            feature_idx,
            distribution: _,
            feature_split_val: _,
        } = self.tree[curr_node]
        {
            let target_value = features[feature_idx];
            curr_node = match self
                .tree
                .edges(curr_node)
                .find(|e| *e.weight() == target_value)
            {
                Some(x) => x.target(),
                // curr_node points to an inner node
                None => break,
            }
        }

        // if the algorithm is cart, then it will always point to a leaf node
        self.tree[curr_node].distribution()
    }
}

impl Node {
    /// Returns the underlying data label distribution for the `Node`.
    pub fn distribution(&self) -> Option<&FxHashMap<usize, Of64>> {
        match self {
            Node::Leaf {
                label: _,
                distribution,
            } => Some(distribution),
            Node::Inner {
                feature_idx: _,
                distribution,
                feature_split_val: _,
            } => distribution.as_ref(),
        }
    }

    /// Performs classification, and returns the index of the predicted class.
    pub fn classify_index(&self) -> usize {
        // Unwrap is safe because self will be Node::Inner only if TrainingAlgorithm::Id3 is used
        let distribution = match self {
            Node::Leaf {
                label: _,
                distribution,
            } => distribution,
            Node::Inner {
                feature_idx: _,
                distribution,
                feature_split_val: _,
            } => distribution.as_ref().unwrap(),
        };

        distribution
            .iter()
            .max_by_key(|(_, v)| **v)
            .map(|(k, _)| *k)
            .expect("leaf node had no datapoints")
    }

    /// Performs classification, and returns the label of the predicted class.
    pub fn classify<L: CategoricalFeature>(&self) -> Option<L> {
        L::from_index(self.classify_index())
    }

    /// Creates a new leaf node.
    pub fn new_leaf<D, L: CategoricalFeature>(data: &DataView<D, L>) -> Self
    where
        D: DataPoint<usize> + Labeled<L>,
    {
        let tot = data.len() as f64;
        let probs = data
            .group_by_label()
            .iter()
            .map(|(&k, v)| (k, Of64::from(v.len() as f64 / tot)))
            .collect::<FxHashMap<usize, Of64>>();
        // There should be at least one point, so unwrap should be infallible
        let label = probs.iter().max_by(|a, b| a.1.cmp(b.1)).unwrap();
        let label = *label.0;
        Self::Leaf {
            label,
            distribution: probs,
        }
    }

    /// Get reference to the underlying distribution, if it has type [`Self::Leaf`].
    fn as_leaf(&self) -> Option<&FxHashMap<usize, Of64>> {
        if let Self::Leaf {
            label: _,
            distribution,
        } = self
        {
            Some(distribution)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(test)]
    use crate::data::CategoricalDataset;
    use crate::data::FeatureSet;

    #[derive(Clone, Copy)]
    struct Coordinate {
        c: usize,
    }

    impl CategoricalFeature for Coordinate {
        const NUM_CATEGORIES: usize = 10;

        fn index(&self) -> usize {
            self.c
        }

        fn from_index(idx: usize) -> Option<Self> {
            Some(Coordinate { c: idx })
        }
    }

    #[derive(Clone)]
    struct TwoDimensionalDiscreteDataPoint {
        features: FeatureSet<usize>,
        label: Label,
        id: uuid::Uuid,
    }

    #[derive(Copy, Clone, Debug, Eq, PartialEq)]
    struct Label(usize);

    impl CategoricalFeature for Label {
        const NUM_CATEGORIES: usize = 2;

        fn index(&self) -> usize {
            self.0
        }

        fn from_index(index: usize) -> Option<Self> {
            Some(Label(index))
        }
    }

    impl TwoDimensionalDiscreteDataPoint {
        #[cfg(test)]
        pub fn new(x: usize, y: usize, label: usize) -> Self {
            let c1 = Coordinate { c: x };
            let c2 = Coordinate { c: y };
            let coords = [c1, c2];

            TwoDimensionalDiscreteDataPoint {
                features: FeatureSet::from_iter(coords),
                label: Label(label),
                id: uuid::Uuid::new_v4(),
            }
        }
    }

    impl Labeled<Label> for TwoDimensionalDiscreteDataPoint {
        fn label(&self) -> Label {
            self.label
        }
    }

    impl DataPoint<usize> for TwoDimensionalDiscreteDataPoint {
        fn features(&self) -> &FeatureSet<usize> {
            &self.features
        }

        fn id(&self) -> uuid::Uuid {
            self.id
        }
    }

    #[test]
    fn test_best_split() {
        // As a heuristic, pick features to split on which don't have overlaps in feature val
        let p1 = TwoDimensionalDiscreteDataPoint::new(0, 9, 1);
        let p2 = TwoDimensionalDiscreteDataPoint::new(1, 8, 1);
        let p3 = TwoDimensionalDiscreteDataPoint::new(2, 7, 1);
        let p4 = TwoDimensionalDiscreteDataPoint::new(3, 6, 0);
        let p5 = TwoDimensionalDiscreteDataPoint::new(4, 5, 0);
        let p6 = TwoDimensionalDiscreteDataPoint::new(5, 4, 0);
        let p7 = TwoDimensionalDiscreteDataPoint::new(6, 3, 0);
        let p8 = TwoDimensionalDiscreteDataPoint::new(7, 2, 0);
        let p9 = TwoDimensionalDiscreteDataPoint::new(8, 1, 1);
        let p10 = TwoDimensionalDiscreteDataPoint::new(9, 0, 1);
        let datapoints = vec![p1, p2, p3, p4, p5, p6, p7, p8, p9, p10];
        let dataset = CategoricalDataset::new(datapoints.clone());

        let dataview = DataView::full(&dataset);

        let hyperparams = DecisionTreeHyperparameters {
            max_depth: 10,
            criterion: SplitQuality::Entropy,
            algorithm: TrainingAlgorithm::Cart,
        };
        let dtree = DecisionTree::<Label>::train(dataview, &hyperparams);
        for p in datapoints {
            assert_eq!(dtree.classify(&p), p.label());
        }
    }
}
