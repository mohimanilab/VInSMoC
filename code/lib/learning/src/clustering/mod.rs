//! A module for clustering algorithms.
//!
//! So far we support hierarchical clustering and KMeans clustering.

use crate::{Centroid, Distance, Of64};
use itertools::Itertools;
use petgraph::graph::{DiGraph, NodeIndex};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

/// An enum to describe the different types of clustering algorithms that we can use for
/// hierarchical clustering models.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum Linkage {
    /// Max linkage.
    ///
    /// Here the linkage distance between two clusters is calculated as the maximum difference
    /// between any constituent points within those clusters.
    Max,
    /// Min linkage.
    ///
    /// Here the linkage distance between two clusters is calculated as the minimum difference
    /// between any constituent points within those clusters.
    Min,
    /// Average linkage.
    ///
    /// Here the linkage distance between two clusters is calculated as the average distance across
    /// all pairs of points, such that one point is in the first cluster and the other point is in
    /// the second cluster.
    Average,
}

impl Default for Linkage {
    fn default() -> Self {
        Self::Max
    }
}

impl Linkage {
    /// Calculate linkage distance between two slices of points which define a distance function.
    ///
    /// # Panics
    /// This function panics if either `set1` or `set2` are empty.
    pub fn linkage_distance<D: Distance>(&self, set1: &[&D], set2: &[&D]) -> Of64 {
        match *self {
            // the largest distance between i and j for i in set1 and j in set2
            Self::Max => set1
                .iter()
                .map(|i| set2.iter().map(|j| i.distance(j)).max().unwrap())
                .max()
                .unwrap(),
            // the smallest distance between i and j for i in set1 and j in set2
            Self::Min => set1
                .iter()
                .map(|i| set2.iter().map(|j| i.distance(j)).min().unwrap())
                .min()
                .unwrap(),
            Self::Average => {
                let mut cum_sum = Of64::from(0f64);
                set1.iter().for_each(|i| {
                    set2.iter().for_each(|j| {
                        cum_sum += i.distance(j);
                    });
                });

                cum_sum / Of64::from(f64::from(set1.len() as u32 * set2.len() as u32))
            }
        }
    }
}

// type alias for nodes in the hierarchical clustering graph
type HcNodeIndex = NodeIndex<u16>;

/// A helper type to denote a pair of nodes.
///
/// This type maintains the invariant that the first node index will always be smaller than the
/// second node index.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
struct NodePair(HcNodeIndex, HcNodeIndex);

impl NodePair {
    fn new(node1: HcNodeIndex, node2: HcNodeIndex) -> Self {
        if node1.index() < node2.index() {
            Self(node1, node2)
        } else {
            Self(node2, node1)
        }
    }
}

/// A datatype to perform hierarchical clustering.
///
/// This structure stores the best clustering for all possible values of k, where k denotes the
/// number of clusters.
pub struct HierarchicalClustering<'db, D: Distance> {
    linkage_type: Linkage,
    clusters: DiGraph<Vec<&'db D>, (), u16>,
    leaves: FxHashSet<HcNodeIndex>,
    /// A map from k to the best set of k clusters that covers all the data points
    history: FxHashMap<usize, Vec<HcNodeIndex>>,
    distance_cache: FxHashMap<NodePair, Of64>,
    sink: Option<HcNodeIndex>,
}

impl<'db, D: Distance + std::fmt::Debug> HierarchicalClustering<'db, D> {
    /// Initialize the clustering model given a slice of `datapoints` and `linkage_type`.
    pub fn init(linkage_type: Linkage, datapoints: &'db [D]) -> Self {
        let n = datapoints.len();
        let mut clusters = DiGraph::with_capacity(datapoints.len(), datapoints.len());
        let mut history = FxHashMap::default();
        let mut leaves = FxHashSet::default();
        let mut distance_cache = FxHashMap::default();

        history.reserve(2 * n - 1);
        history.insert(n, Vec::with_capacity(datapoints.len()));
        leaves.reserve(n);
        // calculated as (2n - 1)^2 - (2n - 1), since we have 2n - 1 nodes, but distance from a
        // node to itself is always zero
        distance_cache.reserve(4 * n * n);

        datapoints.iter().for_each(|point| {
            let point_slice = vec![point];
            let node_idx = clusters.add_node(point_slice);

            history.entry(n).or_insert_with(Vec::new).push(node_idx);
            leaves.insert(node_idx);
        });

        Self {
            linkage_type,
            clusters,
            leaves,
            history,
            distance_cache,
            sink: None,
        }
    }

    /// Get the best set of `k` clusters such that each datapoint belongs to exactly one cluster
    /// and all datapoints belong to a cluster.
    pub fn get_clusters(&self, k: usize) -> Option<impl Iterator<Item = &Vec<&D>>> {
        self.history.get(&k).map(|node_inds| {
            node_inds
                .iter()
                .map(|node_ind| self.clusters.node_weight(*node_ind).unwrap())
        })
    }

    /// Learn the hierarchical clusters.
    ///
    /// This learns the best clustering for all choices of `k` clusters.
    pub fn learn_clusters(&mut self) {
        // while there are still leaves to compute
        while let Some(NodePair(ind1, ind2)) = self.get_closest_cluster_pair() {
            // both of these are removed from leaves
            self.leaves.remove(&ind1);
            self.leaves.remove(&ind2);

            // create a slice over both elements
            let slice1 = &self.clusters[ind1];
            let slice2 = &self.clusters[ind2];
            let mut slice12 = slice1.clone();
            slice12.extend(slice2.iter().copied());

            let node_idx = self.clusters.add_node(slice12);
            self.clusters.add_edge(ind1, node_idx, ());
            self.clusters.add_edge(ind2, node_idx, ());
            self.leaves.insert(node_idx);

            self.history.insert(
                self.leaves.len(),
                self.leaves.iter().copied().collect::<Vec<_>>(),
            );

            if self.leaves.len() == 1 {
                self.sink = Some(node_idx);
            }
        }
    }

    /// Get the closest pair of nodes in this clustering scheme.
    fn get_closest_cluster_pair(&mut self) -> Option<NodePair> {
        if self.leaves.len() <= 1 {
            return None;
        }

        // Otherwise there are at least two leaves, so we shall unwrap accordingly
        let mut ind1 = None;
        let mut ind2 = None;
        let mut best_distance = None;

        // calculate index i and j which are least distance
        // iterate over the node indices in the leaves and skip ones that are the same and return
        // the pair which should be merged according to the linkage type
        for (i, j) in self.leaves.iter().tuple_combinations() {
            let slice_i = &self.clusters[*i];
            let slice_j = &self.clusters[*j];
            // no two node indices should have same slice of point
            let node_pair = NodePair::new(*i, *j);
            self.distance_cache
                .entry(node_pair)
                .or_insert_with(|| self.linkage_type.linkage_distance(slice_i, slice_j));

            let distance = self.distance_cache[&node_pair];

            match best_distance.as_ref() {
                None => {
                    ind1 = Some(i);
                    ind2 = Some(j);
                    best_distance = Some(distance);
                }
                Some(prev_distance) => {
                    use std::cmp::Ordering;
                    match (&distance).cmp(prev_distance) {
                        Ordering::Less => {
                            ind1 = Some(i);
                            ind2 = Some(j);
                            best_distance = Some(distance);
                        }
                        Ordering::Greater => {}
                        Ordering::Equal => {
                            // if the distance is the same, then connect the smaller clusters
                            // unwrap safe since they must be Some if prev_distance is Some
                            let curr_size = self.clusters[*ind1.unwrap()].len()
                                + self.clusters[*ind2.unwrap()].len();
                            let new_size = self.clusters[*i].len() + self.clusters[*j].len();

                            if new_size < curr_size {
                                ind1 = Some(i);
                                ind2 = Some(j);
                                best_distance = Some(distance);
                            }
                        }
                    };
                }
            }
        }

        // at this point there must be at least two nodes, so both of these must be `Some`
        Some(NodePair::new(*ind1.unwrap(), *ind2.unwrap()))
    }
}

/// An enum to describe the different types of KMeans initialization we can use.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum KMeansInitialization {
    /// A simple initialization where we just assign points to clusters in the order they are
    /// provided in.
    ///
    /// More specifically, we assign point n to cluster n % k.
    Simple,
}

impl Default for KMeansInitialization {
    fn default() -> Self {
        KMeansInitialization::Simple
    }
}

/// A structure to perform K-Means clustering.
pub struct KMeans<'a, D: Distance> {
    /// The constant k which denotes the number of clusters
    k: usize,
    /// The dataset over which we perform the clustering
    dataset: &'a [D],
    /// The splits of the dataset
    splits: Vec<Vec<&'a D>>,
    /// The current centroids for `splits`
    centroids: Vec<D>,
}

impl<'a, D: Distance + Centroid + std::fmt::Debug> KMeans<'a, D> {
    /// Initialize a new KMeans clustering object.
    ///
    /// # Panics
    /// This function panics if `k` is zero.
    pub fn init(k: usize, dataset: &'a [D], init: KMeansInitialization) -> Self {
        assert!(k > 0, "cannot have 0 clusters");

        let split_capacity = dataset.len() / k;
        let mut splits: Vec<Vec<&D>> = vec![Vec::with_capacity(split_capacity); k];
        let mut centroids = Vec::with_capacity(k);

        match init {
            KMeansInitialization::Simple => {
                dataset.iter().enumerate().for_each(|(idx, point)| {
                    splits[idx % k].push(point);
                });

                for split in &splits {
                    centroids.push(D::centroid(split.iter().copied()))
                }
            }
        }

        Self {
            k,
            dataset,
            splits,
            centroids,
        }
    }

    /// Return the value of k used to initialize this kmeans object
    pub fn k(&self) -> usize {
        self.k
    }

    /// FIXME: need better stopping criteria
    pub fn learn_clusters(&mut self, n: usize) {
        for _i in 0..n {
            self.update_centroids();
            self.update_splits();
        }
    }

    /// Recalculate the centroids for clustering k
    fn update_centroids(&mut self) {
        // calculate the centroids based on the current split assignments
        self.centroids
            .iter_mut()
            .zip(self.splits.iter())
            .for_each(|(centroid, split)| {
                // centroid should be calculated over slice of references of D
                *centroid = D::centroid(split.iter().copied());
            });
    }

    /// Calculate the closest centroid for a given datapoint.
    ///
    /// The return value of this function is guaranteed to be less than or equal to k.
    fn get_closest_centroid(&self, datapoint: &D) -> usize {
        // unwrap at the end is safe since we panic if k = 0 during init
        self.centroids
            .iter()
            .position_min_by_key(|centroid| centroid.distance(datapoint))
            .unwrap()
    }

    /// Recalculate the dataview splits
    fn update_splits(&mut self) {
        // first clear the current splits
        self.splits.iter_mut().for_each(|dataview| {
            dataview.clear();
        });

        // for each datapoint in the dataset, add it to the split corresponding to the centroid
        // that it is closest to
        self.dataset.iter().for_each(|datapoint| {
            let split_idx = self.get_closest_centroid(datapoint);
            self.splits[split_idx].push(datapoint);
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Pow;

    #[derive(Clone, Debug)]
    struct Point {
        x: f64,
        y: f64,
    }

    impl Point {
        fn new(x: f64, y: f64) -> Self {
            Self { x, y }
        }
    }

    impl PartialEq for Point {
        fn eq(&self, other: &Self) -> bool {
            self.x == other.x && self.y == other.y
        }
    }

    impl Distance for Point {
        fn distance(&self, other: &Self) -> Of64 {
            let x_dist = (self.x - other.x).pow(2u16);
            let y_dist = (self.y - other.y).pow(2u16);
            Of64::from((x_dist + y_dist).sqrt())
        }
    }

    impl Centroid for Point {
        fn centroid<'a, I>(points: I) -> Self
        where
            I: IntoIterator<Item = &'a Self>,
            Self: 'a,
        {
            let mut cum_x = 0f64;
            let mut cum_y = 0f64;

            let mut num_points = 0u32;
            for point in points {
                cum_x += point.x;
                cum_y += point.y;
                num_points += 1;
            }

            Self {
                x: cum_x / num_points as f64,
                y: cum_y / num_points as f64,
            }
        }
    }

    #[test]
    fn test_k_means() {
        let ds: Vec<Point> = vec![
            Point::new(0., 0.),
            Point::new(10., 10.),
            Point::new(1., 0.),
            Point::new(11., 10.),
            Point::new(2., 3.),
        ];

        let mut kmeans = KMeans::init(2, &ds, KMeansInitialization::default());
        kmeans.learn_clusters(10);
        let split_one = if kmeans.splits[0].len() == 3 {
            &kmeans.splits[0]
        } else {
            &kmeans.splits[1]
        };

        assert!(split_one.contains(&&ds[0]));
        assert!(split_one.contains(&&ds[2]));
        assert!(split_one.contains(&&ds[4]));

        let split_two = if kmeans.splits[1].len() != 3 {
            &kmeans.splits[1]
        } else {
            &kmeans.splits[0]
        };

        assert!(split_two.contains(&&ds[1]));
        assert!(split_two.contains(&&ds[3]));
    }

    #[test]
    fn test_hierarchical() {
        let ds: Vec<Point> = vec![
            Point::new(2., 0.),
            Point::new(0., 0.),
            Point::new(16., 0.),
            Point::new(14., 0.),
        ];

        let mut hc = HierarchicalClustering::init(Linkage::Average, &ds);
        hc.learn_clusters();

        // one cluster should have the first two points, and the second should have the second two
        // points
        let clusters = hc.get_clusters(2).unwrap().collect::<Vec<_>>();
        let order1 = clusters[0].contains(&&ds[0])
            && clusters[0].contains(&&ds[1])
            && clusters[1].contains(&&ds[2])
            && clusters[1].contains(&&ds[3]);
        let order2 = clusters[1].contains(&&ds[0])
            && clusters[1].contains(&&ds[1])
            && clusters[0].contains(&&ds[2])
            && clusters[0].contains(&&ds[3]);
        assert!(order1 || order2);
    }
}
