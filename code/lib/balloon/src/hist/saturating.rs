//! Bounded floating point collections that snap out of bounds queries to the closest bound.

use rkyv::{out_field, ser::ScratchSpace, Archive, Archived, Fallible};

use super::Hist;

/// A [`Hist`] that has infallible queries and insertions by replacing out of bounds queries with
/// the closest bound.
#[derive(Debug, Clone)]
pub struct SaturatingHist<T> {
    pub(crate) hist: Hist<T>,
}

impl<T: Default + Clone> SaturatingHist<T> {
    /// Construct a histogram with provided bounds and bucket size.
    ///
    /// The buckets will be filled with the [`Default`] value for `T`.
    pub fn new(min: f64, max: f64, step: f64) -> Self {
        Self {
            hist: Hist::new(min, max, step),
        }
    }
}

impl<T: Default> SaturatingHist<T> {
    /// Set all bins in this histogram to the default value of `T`.
    pub fn clear(&mut self) {
        self.hist.clear();
    }
}

impl<T> SaturatingHist<T> {
    /// Compute the index of `key` in the underlying buffer.
    ///
    /// If the value of `key` is smaller than the lower bound of this histogram this returns `0`.
    /// If the value of `key` is larger than the upper bound of this histogram this returns the
    /// index of the last bucket. Otherwise this returns the index of the histogram bucket for
    /// `key`.
    ///
    /// # Panics
    /// If the provided key is [`f64::NAN`] this panics.
    pub fn index(&self, key: f64) -> usize {
        assert!(!key.is_nan(), "attempted to index with NaN key");
        if key < self.hist.min {
            0
        } else if key > self.hist.max {
            self.hist.values.len() - 1
        } else {
            ((key - self.hist.min) / self.hist.step).round() as usize
        }
    }

    /// Compute the center of the bucket at the given index and return that as the key to the bucket
    pub fn key(&self, index: usize) -> Option<f64> {
        self.hist.key(index)
    }

    /// Get the number of buckets in this histogram.
    pub fn len(&self) -> usize {
        self.hist.len()
    }

    /// Check if this histogram has no buckets.
    pub fn is_empty(&self) -> bool {
        self.hist.is_empty()
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

    /// Get an iterator over bucket keys for each bucket in this histogram.
    ///
    /// Bucket keys are computed using [`Hist::key`].
    pub fn keys(&self) -> Keys<'_, T> {
        Keys {
            inner: self.hist.keys(),
        }
    }

    /// Get an iterator over mutable references to values for each bucket in this histogram.
    pub fn values_mut(&mut self) -> ValuesMut<'_, T> {
        ValuesMut {
            inner: self.hist.values_mut(),
        }
    }

    /// Convert into an iterator over owned bucket keys for each bucket in this histogram.
    pub fn into_keys(self) -> IntoKeys<T> {
        IntoKeys {
            inner: self.hist.into_keys(),
        }
    }

    /// Convert into an iterator over owned values for each bucket in this histogram.
    pub fn into_values(self) -> IntoValues<T> {
        IntoValues {
            inner: self.hist.into_values(),
        }
    }

    /// Get an iterator over bucket keys and values as tuples.
    pub fn iter(&self) -> Iter<'_, T> {
        Iter {
            inner: self.hist.iter(),
        }
    }

    /// Get the histogram's minimum value.
    pub fn min(&self) -> f64 {
        self.hist.min()
    }

    /// Get the histogram's maximum value.
    pub fn max(&self) -> f64 {
        self.hist.max()
    }

    /// Get the step between buckets in the histogram.
    pub fn step(&self) -> f64 {
        self.hist.step()
    }
}

/// An iterator over references to values in a [`SaturatingHist`].
///
/// Constructed using [`SaturatingHist::values`].
#[derive(Debug, Clone)]
pub struct Values<'a, T> {
    inner: super::Values<'a, T>,
}
inheriting_iter! { Values<'a, T>, &'a T, 'a, T }
inheriting_iter! { @all Values<'a, T>, 'a, T }

/// An iterator over the floating point bucket keys in a [`SaturatingHist`].
///
/// Constructed using [`SaturatingHist::keys`].
#[derive(Debug, Clone)]
pub struct Keys<'a, T> {
    inner: super::Keys<'a, T>,
}
inheriting_iter! { Keys<'a, T>, f64, 'a, T }
inheriting_iter! { @all Keys<'a, T>, 'a, T }

/// An iterator over mutable references to values in a [`SaturatingHist`].
///
/// Constructed using [`SaturatingHist::values_mut`].
#[derive(Debug)]
pub struct ValuesMut<'a, T> {
    inner: super::ValuesMut<'a, T>,
}
inheriting_iter! { ValuesMut<'a, T>, &'a mut T, 'a, T }
inheriting_iter! { @all ValuesMut<'a, T>, 'a, T }

/// An iterator over owned floating point bucket keys in a [`SaturatingHist`].
///
/// Constructed using [`SaturatingHist::into_keys`].
#[derive(Debug, Clone)]
pub struct IntoKeys<T> {
    inner: super::IntoKeys<T>,
}
inheriting_iter! { IntoKeys<T>, f64, T }
inheriting_iter! { @all IntoKeys<T>, T }

/// An iterator over owned values in a [`SaturatingHist`].
///
/// Constructed using [`SaturatingHist::into_values`].
#[derive(Debug, Clone)]
pub struct IntoValues<T> {
    inner: super::IntoValues<T>,
}
inheriting_iter! { IntoValues<T>, T, T }
inheriting_iter! { @all IntoValues<T>, T }

/// An iterator over bucket keys and values in a [`SaturatingHist`].
///
/// Constructed using [`SaturatingHist::iter`].
#[derive(Debug, Clone)]
pub struct Iter<'a, T> {
    inner: super::Iter<'a, T>,
}
inheriting_iter! { Iter<'a, T>, (f64, &'a T), 'a, T }
inheriting_iter! { @all Iter<'a, T>, 'a, T }

/// An iterator over owned bucket keys and owned values in a [`SaturatingHist`].
///
/// Constructed using [`IntoIterator`].
#[derive(Debug, Clone)]
pub struct IntoIter<T> {
    inner: super::IntoIter<T>,
}
inheriting_iter! { IntoIter<T>, (f64, T), T }
inheriting_iter! { @all IntoIter<T>, T }

impl<T> IntoIterator for SaturatingHist<T> {
    type Item = (f64, T);

    type IntoIter = IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIter {
            inner: self.hist.into_iter(),
        }
    }
}

/// An [`rkyv`] compatible version of [`SaturatingHist`]
pub struct ArchivedSaturatingHist<T: Archive> {
    pub(crate) hist: Archived<Hist<T>>,
}

impl<T: Archive> Archive for SaturatingHist<T> {
    type Archived = ArchivedSaturatingHist<T>;

    type Resolver = <Hist<T> as Archive>::Resolver;

    unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
        let (fp, fo) = out_field!(out.hist);
        self.hist.resolve(pos + fp, resolver, fo);
    }
}

impl<S, T> rkyv::Serialize<S> for SaturatingHist<T>
where
    S: rkyv::ser::Serializer + ScratchSpace + ?Sized,
    T: Archive + rkyv::Serialize<S>,
{
    fn serialize(
        &self,
        serializer: &mut S,
    ) -> Result<Self::Resolver, <S as rkyv::Fallible>::Error> {
        self.hist.serialize(serializer)
    }
}

impl<T, D> rkyv::Deserialize<SaturatingHist<T>, D> for ArchivedSaturatingHist<T>
where
    T: Archive,
    T::Archived: rkyv::Deserialize<T, D>,
    D: Fallible + ?Sized,
{
    fn deserialize(
        &self,
        deserializer: &mut D,
    ) -> Result<SaturatingHist<T>, <D as Fallible>::Error> {
        Ok(SaturatingHist {
            hist: self.hist.deserialize(deserializer)?,
        })
    }
}

impl<T: Archive> ArchivedSaturatingHist<T> {
    /// Compute the index of `key` in the underlying buffer.
    ///
    /// If the value of `key` is smaller than the lower bound of this histogram this returns `0`.
    /// If the value of `key` is larger than the upper bound of this histogram this returns the
    /// index of the last bucket. Otherwise this returns the index of the histogram bucket for
    /// `key`.
    ///
    /// # Panics
    /// If the provided key is [`f64::NAN`] this panics.
    pub fn index(&self, key: f64) -> usize {
        assert!(!key.is_nan(), "attempted to index with NaN key");
        if key < self.hist.min {
            0
        } else if key > self.hist.max {
            self.hist.values.len() - 1
        } else {
            ((key - self.hist.min) / self.hist.step).round() as usize
        }
    }
}
