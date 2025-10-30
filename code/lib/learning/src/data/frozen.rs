//! Sets of features for datapoints.

use super::{categorical::CategoricalFeature, numeric::NumericFeature, onehot::OneHot};

/// A marker trait for types that can be a single feature for [`FeatureSet`].
///
/// This trait cannot be implemented on your own types.
pub trait FeatureUnit: sealed::Sealed {}
mod sealed {
    pub trait Sealed {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
    impl Sealed for usize {}
    impl Sealed for bool {}
}

impl FeatureUnit for f32 {}
impl FeatureUnit for f64 {}
impl FeatureUnit for usize {}
impl FeatureUnit for bool {}

/// A [`FeatureSet`] is a flattened collection of features that are either all numeric or all
/// categorical.
///
/// A [`FeatureSet<usize>`] is a collection of categorical features. Both [`FeatureSet<f32>`] and
/// [`FeatureSet<f64>`] are collections of numeric features.
///
/// # Construction
/// A feature set can be built from an iterator over categorical features, an iterator over numeric
/// features, a single categorical feature, or a single numeric feature. It can also be built via
/// its [`Default`] implementation and the various `push_` methods.
#[derive(Clone, Default)]
pub struct FeatureSet<T: FeatureUnit> {
    features: Vec<T>,
}

impl<T: FeatureUnit> std::ops::Deref for FeatureSet<T> {
    type Target = Vec<T>;

    fn deref(&self) -> &Self::Target {
        &self.features
    }
}

impl<T: FeatureUnit> FeatureSet<T> {
    /// Get a reference to the frozen's features.
    pub fn features(&self) -> &Vec<T> {
        &self.features
    }

    /// Get the number of features in this feature set.
    pub fn len(&self) -> usize {
        self.features.len()
    }

    /// Check if this feature set is empty or not.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl FeatureSet<f32> {
    /// Add a numeric feature to this feature set.
    pub fn push_numeric<T: NumericFeature>(&mut self, value: T) {
        self.features.extend(value.value().into_singles());
    }

    /// Add a categorical feature to this feature set via 1-hot encoding.
    pub fn push_onehot<T: CategoricalFeature>(&mut self, value: T) {
        let oh = OneHot::<_, f32>::new(value);
        self.features.extend(oh.value().into_singles());
    }

    /// Add a categorical feature to this feature set by convert its index into a numeric value.
    ///
    /// This should generally be avoided unless your feature is ordinal. For categorical features
    /// that are not ordinal the numeric representation of that feature will be entirely
    /// meaningless. If your feature is categorical but not ordinal prefer
    /// [`FeatureSet::push_onehot`].
    pub fn push_categorical<T: CategoricalFeature>(&mut self, value: T) {
        self.features.push(value.index() as f32);
    }
}

impl FeatureSet<f64> {
    /// Add a numeric feature to this feature set.
    pub fn push_numeric<T: NumericFeature>(&mut self, value: T) {
        self.features.extend(value.value().into_doubles());
    }

    /// Add a categorical feature to this feature set via 1-hot encoding.
    pub fn push_onehot<T: CategoricalFeature>(&mut self, value: T) {
        let oh = OneHot::<_, f64>::new(value);
        self.features.extend(oh.value().into_doubles());
    }

    /// Add a categorical feature to this feature set by convert its index into a numeric value.
    ///
    /// This should generally be avoided unless your feature is ordinal. For categorical features
    /// that are not ordinal the numeric representation of that feature will be entirely
    /// meaningless. If your feature is categorical but not ordinal prefer
    /// [`FeatureSet::push_onehot`].
    pub fn push_categorical<T: CategoricalFeature>(&mut self, value: T) {
        self.features.push(value.index() as f64);
    }
}

impl FeatureSet<usize> {
    /// Add a categorical feature to this feature set.
    pub fn push<T: CategoricalFeature>(&mut self, value: T) {
        self.features.push(value.index());
    }
}

impl<T: CategoricalFeature> From<T> for FeatureSet<usize> {
    fn from(cat: T) -> Self {
        Self {
            features: vec![cat.index()],
        }
    }
}

impl<T: CategoricalFeature> FromIterator<T> for FeatureSet<usize> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self {
            features: iter.into_iter().map(|x| x.index()).collect(),
        }
    }
}

impl<T: NumericFeature> From<T> for FeatureSet<f32> {
    fn from(num: T) -> Self {
        Self {
            features: num.value().into_singles().collect(),
        }
    }
}

impl<T: NumericFeature> FromIterator<T> for FeatureSet<f32> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self {
            features: iter
                .into_iter()
                .flat_map(|x| x.value().into_singles())
                .collect(),
        }
    }
}

impl<T: NumericFeature> From<T> for FeatureSet<f64> {
    fn from(num: T) -> Self {
        Self {
            features: num.value().into_doubles().collect(),
        }
    }
}

impl<T: NumericFeature> FromIterator<T> for FeatureSet<f64> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self {
            features: iter
                .into_iter()
                .flat_map(|x| x.value().into_doubles())
                .collect(),
        }
    }
}
