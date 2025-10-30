//! Categorical features.

/// A trait marking a type as a categorical feature.
///
/// Categorical features can take one of a predefined set of values. We also introduce the
/// additional constraint that these predefined values can be represented as indices from 0 to
/// [`CategoricalFeature::NUM_CATEGORIES`]. This is intended for implementation on enums.
pub trait CategoricalFeature: Sync + Send + Copy + Sized {
    /// The total possible values of this categorical feature.
    const NUM_CATEGORIES: usize;

    /// Convert this feature into its canonical feature index.
    ///
    /// This index is used for equality comparisons, so if your implementation does not map to
    /// unique indices for each logically distinct category it will be incorrect.
    fn index(&self) -> usize;

    /// Construct this feature from its index.
    ///
    /// If the `index` is invalid this function should return [`None`]. This function should be the
    /// direct inverse of [`CategoricalFeature::index`], otherwise the implementation is incorrect.
    fn from_index(index: usize) -> Option<Self>;

    /// Construct this feature from its index, panicking on invalid indices.
    ///
    /// The default implementation simply unwraps a [`CategoricalFeature::from_index`] call.
    ///
    /// # Panics
    /// If the `index` is not valid this panics.
    fn from_index_unchecked(index: usize) -> Self {
        Self::from_index(index).expect("invalid categorical feature index")
    }
}

impl CategoricalFeature for bool {
    const NUM_CATEGORIES: usize = 2;

    fn index(&self) -> usize {
        *self as usize
    }

    fn from_index(index: usize) -> Option<Self> {
        match index {
            0 => Some(false),
            1 => Some(true),
            _ => None,
        }
    }
}

impl<T: CategoricalFeature> CategoricalFeature for Option<T> {
    const NUM_CATEGORIES: usize = T::NUM_CATEGORIES + 1;

    fn index(&self) -> usize {
        match self {
            Some(x) => x.index(),
            None => Self::NUM_CATEGORIES - 1,
        }
    }

    fn from_index(index: usize) -> Option<Self> {
        match index {
            x if x == Self::NUM_CATEGORIES - 1 => Some(None),
            x => T::from_index(x).map(Some),
        }
    }
}
