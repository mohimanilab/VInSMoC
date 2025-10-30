//! One-hot encoding for categorical features.

use std::marker::PhantomData;

use super::{
    categorical::CategoricalFeature,
    numeric::{NumericFeature, NumericValue},
};

/// A one-hot encoded feature.
///
/// A one-hot encoding serves as a bridge between a categorical feature and a numeric feature.
/// Categorical features should not use their index identifier in general as a numeric feature, as
/// these index identifiers are largely meaningless except in the special case where the feature is
/// ordinal. A more general method is create a vector of length
/// [`CategoricalFeature::NUM_CATEGORIES`] that is zero valued everywhere except for the index
/// identifying this feature value, where it's one.
///
/// The precision of the numeric one-hot encoding can be selected with the `F` generic parameter.
///
/// # Example
/// ```
/// use learning::data::{
///     CategoricalFeature,
///     NumericFeature,
///     NumericValue::{F32Vec, F64Vec},
///     OneHot
/// };
///
/// #[derive(Clone, Copy)]
/// #[repr(u8)]
/// enum Feat {
///     A = 0,
///     B,
///     C,
///     D
/// }
///
/// impl CategoricalFeature for Feat {
///     const NUM_CATEGORIES: usize = 4;
///
///     fn index(&self) -> usize {
///         *self as usize
///     }
///
///     fn from_index(index: usize) -> Option<Self> {
///         match index {
///             0 => Some(Feat::A),
///             1 => Some(Feat::B),
///             2 => Some(Feat::C),
///             3 => Some(Feat::D),
///             _ => None
///         }
///     }
/// }
///
/// assert_eq!(OneHot::<_, f32>::new(Feat::A).value(), F32Vec(vec![1.0f32, 0.0, 0.0, 0.0]));
/// assert_eq!(OneHot::<_, f32>::new(Feat::B).value(), F32Vec(vec![0.0f32, 1.0, 0.0, 0.0]));
/// assert_eq!(OneHot::<_, f32>::new(Feat::C).value(), F32Vec(vec![0.0f32, 0.0, 1.0, 0.0]));
/// assert_eq!(OneHot::<_, f32>::new(Feat::D).value(), F32Vec(vec![0.0f32, 0.0, 0.0, 1.0]));
/// assert_eq!(OneHot::<_, f64>::new(Feat::A).value(), F64Vec(vec![1.0f64, 0.0, 0.0, 0.0]));
/// assert_eq!(OneHot::<_, f64>::new(Feat::B).value(), F64Vec(vec![0.0f64, 1.0, 0.0, 0.0]));
/// assert_eq!(OneHot::<_, f64>::new(Feat::C).value(), F64Vec(vec![0.0f64, 0.0, 1.0, 0.0]));
/// assert_eq!(OneHot::<_, f64>::new(Feat::D).value(), F64Vec(vec![0.0f64, 0.0, 0.0, 1.0]));
/// ```
pub struct OneHot<T, F = f32> {
    inner: T,
    marker: PhantomData<F>,
}

impl<T: CategoricalFeature> OneHot<T, f32> {
    /// Construct a new 1-hot encoded feature from a categorical feature.
    pub fn new(inner: T) -> Self {
        Self {
            inner,
            marker: PhantomData,
        }
    }
}

impl<T: CategoricalFeature> OneHot<T, f64> {
    /// Construct a new 1-hot encoded feature from a categorical feature.
    pub fn new(inner: T) -> Self {
        Self {
            inner,
            marker: PhantomData,
        }
    }
}

impl<T: CategoricalFeature> From<T> for OneHot<T, f32> {
    fn from(inner: T) -> Self {
        Self::new(inner)
    }
}

impl<T: CategoricalFeature> From<T> for OneHot<T, f64> {
    fn from(inner: T) -> Self {
        Self::new(inner)
    }
}

impl<T: CategoricalFeature> NumericFeature for OneHot<T, f32> {
    fn value(&self) -> NumericValue {
        let mut one_hot = vec![0.0f32; T::NUM_CATEGORIES];
        one_hot[self.inner.index()] = 1.0;
        one_hot.into()
    }
}

impl<T: CategoricalFeature> NumericFeature for OneHot<T, f64> {
    fn value(&self) -> NumericValue {
        let mut one_hot = vec![0.0f64; T::NUM_CATEGORIES];
        one_hot[self.inner.index()] = 1.0;
        one_hot.into()
    }
}
