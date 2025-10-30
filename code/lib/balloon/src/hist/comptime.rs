//! Compile-time versions of [`crate::hist::Hist`] and [`crate::hist::SaturatingHist`].
//!
//! These versions of the general-purpose histograms are backed by an array instead of a vector,
//! meaning that their constructors are constant functions.

use std::iter::FusedIterator;

pub use saturating::SaturatingHist;

/// A collection of bounded floating point keys associated with some value, backed by an array.
///
/// This histogram has buckets for all possible floating point values within the minimum and
/// maximum bounds provided during its construction. In exchange for the potentially high memory
/// usage because of the eager bucket initialization lookup is some minor arithmetic and an array
/// indexing operation, making it performant. This particular histogram rejects values that are
/// outside of its bounds.
///
/// If you would like a more flexible histogram with regards to bounds see [`SaturatingHist`]. If
/// you want dynamically allocated storage see [`crate::hist::Hist`].
#[derive(Debug, Clone)]
pub struct Hist<T, const N: usize> {
    values: [T; N],
    min: f64,
    max: f64,
    step: f64,
}

impl<T, const N: usize> TryFrom<super::Hist<T>> for Hist<T, N> {
    type Error = super::Hist<T>;

    fn try_from(hist: super::Hist<T>) -> Result<Self, Self::Error> {
        hist.values
            .try_into()
            .map(|values| Self {
                values,
                min: hist.min,
                max: hist.max,
                step: hist.step,
            })
            .map_err(|values| super::Hist {
                values,
                min: hist.min,
                max: hist.max,
                step: hist.step,
            })
    }
}

impl<T: Default + Copy, const N: usize> Hist<T, N> {
    /// Construct a histogram with provided bounds and bucket size.
    ///
    /// The buckets will be filled with the [`Default`] value for `T`.
    pub fn new(min: f64, max: f64, step: f64) -> Self {
        let values = [T::default(); N];
        Self {
            values,
            min,
            max,
            step,
        }
    }
}

impl<T, const N: usize> Hist<T, N> {
    /// Construct a histogram with the provided `values` as bucket values.
    pub const fn from_values(values: [T; N], min: f64, max: f64, step: f64) -> Self {
        Self {
            values,
            min,
            max,
            step,
        }
    }

    /// Compute the index of `key` in the underlying buffer.
    ///
    /// If the value of `key` is outside of the bounds of this histogram this method returns
    /// [`None`]. Otherwise this returns the index of the histogram bucket for `key`.
    pub fn index(&self, key: f64) -> Option<usize> {
        if key < self.min || key > self.max {
            return None;
        }
        Some(((key - self.min) / self.step).round() as usize)
    }

    /// Get the value inside the bucket at `index`.
    ///
    /// If the index is out of bounds this returns [`None`].
    pub fn bucket(&self, index: usize) -> Option<&T> {
        self.values.get(index)
    }

    /// Get the number of buckets in this histogram.
    pub const fn len(&self) -> usize {
        self.values.len()
    }

    /// Check if this histogram has no buckets.
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get a reference to the value in the bucket for `key`.
    ///
    /// If the `key` is out of bounds this returns [`None`].
    pub fn get(&self, key: f64) -> Option<&T> {
        let index = self.index(key)?;
        Some(&self.values[index])
    }

    /// Get a mutable reference to the value in the bucket for `key`.
    ///
    /// If the `key` is out of bounds this returns [`None`].
    pub fn get_mut(&mut self, key: f64) -> Option<&mut T> {
        let index = self.index(key)?;
        Some(&mut self.values[index])
    }

    /// Get the histogram's bucket values.
    pub const fn raw_values(&self) -> &[T] {
        &self.values
    }

    /// Get the histogram's minimum value.
    pub const fn min(&self) -> f64 {
        self.min
    }

    /// Get the histogram's maximum value.
    pub const fn max(&self) -> f64 {
        self.max
    }

    /// Get the step between buckets in the histogram.
    pub const fn step(&self) -> f64 {
        self.step
    }

    /// Get an iterator over references to values for each bucket in this histogram.
    pub fn values(&self) -> Values<'_, T> {
        Values {
            inner: self.values.iter(),
        }
    }

    /// Get an iterator over mutable references to values for each bucket in this histogram.
    pub fn values_mut(&mut self) -> ValuesMut<'_, T> {
        ValuesMut {
            inner: self.values.iter_mut(),
        }
    }

    /// Convert into an iterator over owned values for each bucket in this histogram.
    pub fn into_values(self) -> IntoValues<T, N> {
        IntoValues {
            inner: self.values.into_iter(),
        }
    }
}

/// An iterator over references to values in a [`Hist`].
///
/// Constructed using [`Hist::values`].
pub struct Values<'a, T> {
    inner: std::slice::Iter<'a, T>,
}
inheriting_iter! { Values<'a, T>, &'a T, 'a, T }
inheriting_iter! { @all Values<'a, T>, 'a, T }

/// An iterator over mutable references to values in a [`Hist`].
///
/// Constructed using [`Hist::values_mut`].
pub struct ValuesMut<'a, T> {
    inner: std::slice::IterMut<'a, T>,
}
inheriting_iter! { ValuesMut<'a, T>, &'a mut T, 'a, T }
inheriting_iter! { @all ValuesMut<'a, T>, 'a, T }

/// An iterator over owned values in a [`Hist`].
///
/// Constructed using [`Hist::into_values`].
pub struct IntoValues<T, const N: usize> {
    inner: std::array::IntoIter<T, N>,
}

impl<T, const N: usize> Iterator for IntoValues<T, N> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<T, const N: usize> FusedIterator for IntoValues<T, N> {}

impl<T, const N: usize> ExactSizeIterator for IntoValues<T, N> {
    fn len(&self) -> usize {
        self.inner.len()
    }
}

impl<T, const N: usize> DoubleEndedIterator for IntoValues<T, N> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner.next_back()
    }
}

/// Bounded floating point collections that snap out of bounds queries to the closest bound.
pub mod saturating {

    use std::iter::FusedIterator;

    use super::Hist;

    /// A [`Hist`] that has infallible queries and insertions by replacing out of bounds queries with
    /// the closest bound.
    #[derive(Debug, Clone)]
    pub struct SaturatingHist<T, const N: usize> {
        hist: Hist<T, N>,
    }

    impl<T, const N: usize> TryFrom<crate::hist::SaturatingHist<T>> for SaturatingHist<T, N> {
        type Error = crate::hist::SaturatingHist<T>;

        fn try_from(hist: crate::hist::SaturatingHist<T>) -> Result<Self, Self::Error> {
            Ok(Self {
                hist: super::Hist::try_from(hist.hist)
                    .map_err(|hist| crate::hist::SaturatingHist { hist })?,
            })
        }
    }

    impl<T: Default + Copy, const N: usize> SaturatingHist<T, N> {
        /// Construct a histogram with provided bounds and bucket size.
        ///
        /// The buckets will be filled with the [`Default`] value for `T`.
        pub fn new(min: f64, max: f64, step: f64) -> Self {
            Self {
                hist: Hist::new(min, max, step),
            }
        }
    }

    impl<T, const N: usize> SaturatingHist<T, N> {
        /// Construct a histogram with the provided `values` as bucket values.
        pub const fn from_values(values: [T; N], min: f64, max: f64, step: f64) -> Self {
            Self {
                hist: Hist::from_values(values, min, max, step),
            }
        }

        /// Compute the index of `key` in the underlying buffer.
        ///
        /// If the value of `key` is smaller than the lower bound of this histogram this returns `0`.
        /// If the value of `key` is larger than the upper bound of this histogram this returns the
        /// index of the last bucket. Otherwise this returns the index of the histogram bucket for
        /// `key`.
        pub fn index(&self, key: f64) -> usize {
            if key < self.hist.min {
                0
            } else if key > self.hist.max {
                self.hist.values.len() - 1
            } else {
                ((key - self.hist.min) / self.hist.step).round() as usize
            }
        }

        /// Get the value of bucket at `index`.
        ///
        /// # Panics
        /// If `index` is out of bounds this panics.
        pub const fn bucket(&self, index: usize) -> &T {
            &self.hist.values[index]
        }

        /// Get the number of buckets in this histogram.
        pub const fn len(&self) -> usize {
            self.hist.len()
        }

        /// Check if this histogram has no buckets.
        pub const fn is_empty(&self) -> bool {
            self.len() == 0
        }

        /// Insert `value` into bucket of `key`.
        ///
        /// This returns the old value in that bucket.
        pub fn insert(&mut self, key: f64, mut value: T) -> T {
            let index = self.index(key);
            std::mem::swap(&mut self.hist.values[index], &mut value);
            value
        }

        /// Get a reference to the value of the bucket of `key`.
        pub fn get(&self, key: f64) -> &T {
            &self.hist.values[self.index(key)]
        }

        /// Get a mutable reference to the value of the bucket of `key`.
        pub fn get_mut(&mut self, key: f64) -> &mut T {
            let index = self.index(key);
            &mut self.hist.values[index]
        }

        /// Get an iterator over references to values for each bucket in this histogram.
        pub fn values(&self) -> Values<'_, T> {
            Values {
                inner: self.hist.values(),
            }
        }

        /// Get an iterator over mutable references to values for each bucket in this histogram.
        pub fn values_mut(&mut self) -> ValuesMut<'_, T> {
            ValuesMut {
                inner: self.hist.values_mut(),
            }
        }

        /// Convert into an iterator over owned values for each bucket in this histogram.
        pub fn into_values(self) -> IntoValues<T, N> {
            IntoValues {
                inner: self.hist.into_values(),
            }
        }

        /// Get the histogram's minimum value.
        pub const fn min(&self) -> f64 {
            self.hist.min()
        }

        /// Get the histogram's maximum value.
        pub const fn max(&self) -> f64 {
            self.hist.max()
        }

        /// Get the step between buckets in the histogram.
        pub const fn step(&self) -> f64 {
            self.hist.step()
        }
    }

    /// An iterator over references to values in a [`SaturatingHist`].
    ///
    /// Constructed using [`SaturatingHist::values`].
    pub struct Values<'a, T> {
        inner: super::Values<'a, T>,
    }
    inheriting_iter! { Values<'a, T>, &'a T, 'a, T }
    inheriting_iter! { @all Values<'a, T>, 'a, T }

    /// An iterator over mutable references to values in a [`SaturatingHist`].
    ///
    /// Constructed using [`SaturatingHist::values_mut`].
    pub struct ValuesMut<'a, T> {
        inner: super::ValuesMut<'a, T>,
    }
    inheriting_iter! { ValuesMut<'a, T>, &'a mut T, 'a, T }
    inheriting_iter! { @all ValuesMut<'a, T>, 'a, T }

    /// An iterator over owned values in a [`SaturatingHist`].
    ///
    /// Constructed using [`SaturatingHist::into_values`].
    pub struct IntoValues<T, const N: usize> {
        inner: super::IntoValues<T, N>,
    }

    impl<T, const N: usize> Iterator for IntoValues<T, N> {
        type Item = T;

        fn next(&mut self) -> Option<Self::Item> {
            self.inner.next()
        }

        fn size_hint(&self) -> (usize, Option<usize>) {
            self.inner.size_hint()
        }
    }

    impl<T, const N: usize> FusedIterator for IntoValues<T, N> {}

    impl<T, const N: usize> ExactSizeIterator for IntoValues<T, N> {
        fn len(&self) -> usize {
            self.inner.len()
        }
    }

    impl<T, const N: usize> DoubleEndedIterator for IntoValues<T, N> {
        fn next_back(&mut self) -> Option<Self::Item> {
            self.inner.next_back()
        }
    }
}
