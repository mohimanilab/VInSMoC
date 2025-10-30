//! Histograms for bounded collections of floating point values.

use std::{
    iter::{Enumerate, Zip},
    ops::Range,
};

use rkyv::{
    out_field,
    ser::ScratchSpace,
    vec::{ArchivedVec, VecResolver},
    Archive, Archived, Fallible,
};

pub mod comptime;
pub mod histmap;
pub mod saturating;

pub use saturating::SaturatingHist;

/// A collection of bounded floating point keys associated with some value.
///
/// This histogram has buckets for all possible floating point values within the minimum and
/// maximum bounds provided during its construction. In exchange for the potentially high memory
/// usage because of the eager bucket initialization lookup is some minor arithmetic and an array
/// indexing operation, making it performant. This particular histogram rejects values that are
/// outside of its bounds.
///
/// If you would like a more flexible histogram with regards to bounds see [`SaturatingHist`].
#[derive(Debug, Clone)]
pub struct Hist<T> {
    values: Vec<T>,
    min: f64,
    max: f64,
    step: f64,
}

impl<T: Default + Clone> Hist<T> {
    /// Construct a histogram with provided bounds and bucket size.
    ///
    /// The buckets will be filled with the [`Default`] value for `T`.
    pub fn new(min: f64, max: f64, step: f64) -> Self {
        let mut out = Self {
            values: vec![],
            min,
            max,
            step,
        };
        let max_index = out.index(max).unwrap();
        let values = vec![T::default(); max_index + 1];
        out.values = values;
        out
    }
}

impl<T: Default> Hist<T> {
    /// Set all bins in this histogram to the default value of `T`.
    pub fn clear(&mut self) {
        for v in &mut self.values {
            *v = T::default();
        }
    }
}

impl<T> Hist<T> {
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

    /// Compute the center of the bucket at the given index and return that as the key to the bucket
    pub fn key(&self, index: usize) -> Option<f64> {
        if index >= self.len() {
            return None;
        }
        Some(self.min + ((index as f64) * self.step))
    }

    /// Get the value inside the bucket at `index`.
    ///
    /// If the index is out of bounds this returns [`None`].
    pub fn bucket(&self, index: usize) -> Option<&T> {
        self.values.get(index)
    }

    /// Insert `value` into the bucket for `key`.
    ///
    /// This returns [`None`] if the `key` is out of bounds. Otherwise, it returns the value that
    /// was previously in that bucket.
    pub fn insert(&mut self, key: f64, mut value: T) -> Option<T> {
        let index = self.index(key)?;
        std::mem::swap(&mut self.values[index], &mut value);
        Some(value)
    }

    /// Get the number of buckets in this histogram.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Check if this histogram has no buckets.
    pub fn is_empty(&self) -> bool {
        self.values.len() == 0
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

    /// Get an iterator over references to values for each bucket in this histogram.
    pub fn values(&self) -> Values<'_, T> {
        Values {
            inner: self.values.iter(),
        }
    }

    /// Get an iterator over bucket keys for each bucket in this histogram.
    ///
    /// Bucket keys are computed using [`Hist::key`].
    pub fn keys(&self) -> Keys<'_, T> {
        let inner = 0..self.values.len();
        Keys { hist: self, inner }
    }

    /// Get an iterator over mutable references to values for each bucket in this histogram.
    pub fn values_mut(&mut self) -> ValuesMut<'_, T> {
        ValuesMut {
            inner: self.values.iter_mut(),
        }
    }

    /// Convert into an iterator over owned bucket keys for each bucket in this histogram.
    pub fn into_keys(self) -> IntoKeys<T> {
        let inner = 0..self.values.len();
        IntoKeys { hist: self, inner }
    }

    /// Convert into an iterator over owned values for each bucket in this histogram.
    pub fn into_values(self) -> IntoValues<T> {
        IntoValues {
            inner: self.values.into_iter(),
        }
    }

    /// Get an iterator over bucket keys and values as tuples.
    pub fn iter(&self) -> Iter<'_, T> {
        Iter {
            inner: self.keys().zip(self.values()),
        }
    }

    /// Get the histogram's minimum value.
    pub fn min(&self) -> f64 {
        self.min
    }

    /// Get the histogram's maximum value.
    pub fn max(&self) -> f64 {
        self.max
    }

    /// Get the step between buckets in the histogram.
    pub fn step(&self) -> f64 {
        self.step
    }
}

/// An iterator over references to values in a [`Hist`].
///
/// Constructed using [`Hist::values`].
#[derive(Debug, Clone)]
pub struct Values<'a, T> {
    inner: std::slice::Iter<'a, T>,
}
inheriting_iter! { Values<'a, T>, &'a T, 'a, T }
inheriting_iter! { @all Values<'a, T>, 'a, T }

/// An iterator over all floating point bucket keys in a [`Hist`].
///
/// Constructed using [`Hist::keys`].
#[derive(Debug, Clone)]
pub struct Keys<'a, T> {
    hist: &'a Hist<T>,
    inner: Range<usize>,
}

impl<'a, T> Iterator for Keys<'a, T> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        let index = self.inner.next()?;
        self.hist.key(index)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

inheriting_iter! { @exact Keys<'a, T>, 'a, T }
inheriting_iter! { @fused Keys<'a, T>, 'a, T }

impl<'a, T> DoubleEndedIterator for Keys<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let index = self.inner.next_back()?;
        self.hist.key(index)
    }
}

/// An iterator over mutable references to values in a [`Hist`].
///
/// Constructed using [`Hist::values_mut`].
#[derive(Debug)]
pub struct ValuesMut<'a, T> {
    inner: std::slice::IterMut<'a, T>,
}
inheriting_iter! { ValuesMut<'a, T>, &'a mut T, 'a, T }
inheriting_iter! { @all ValuesMut<'a, T>, 'a, T }

/// An iterator over owned floating point bucket keys in a [`Hist`].
///
/// Constructed using [`Hist::into_keys`].
#[derive(Debug, Clone)]
pub struct IntoKeys<T> {
    hist: Hist<T>,
    inner: Range<usize>,
}

impl<T> Iterator for IntoKeys<T> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        let index = self.inner.next()?;
        self.hist.key(index)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}
inheriting_iter! { @exact IntoKeys<T>, T }
inheriting_iter! { @fused IntoKeys<T>, T }

impl<'a, T> DoubleEndedIterator for IntoKeys<T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let index = self.inner.next_back()?;
        self.hist.key(index)
    }
}

/// An iterator over owned values in a [`Hist`].
///
/// Constructed using [`Hist::into_values`].
#[derive(Debug, Clone)]
pub struct IntoValues<T> {
    inner: std::vec::IntoIter<T>,
}
inheriting_iter! { IntoValues<T>, T, T }
inheriting_iter! { @all IntoValues<T>, T }

/// An iterator over bucket keys and values in a [`Hist`].
///
/// Constructed using [`Hist::iter`].
#[derive(Debug, Clone)]
pub struct Iter<'a, T> {
    inner: Zip<Keys<'a, T>, Values<'a, T>>,
}
inheriting_iter! { Iter<'a, T>, (f64, &'a T), 'a, T }
inheriting_iter! { @all Iter<'a, T>, 'a, T }

/// An iterator over owned bucket keys and owned values in a [`Hist`].
///
/// Constructed using [`IntoIterator`].
#[derive(Debug, Clone)]
pub struct IntoIter<T> {
    inner: Enumerate<std::vec::IntoIter<T>>,
    min: f64,
    step: f64,
}

impl<T> IntoIter<T> {
    /// Compute the center of the bucket at the given index and return that as the key to the bucket
    ///
    /// Assumes a valid index will be given.
    fn key(&self, index: usize) -> f64 {
        self.min + ((index as f64) * self.step)
    }
}

impl<T> Iterator for IntoIter<T> {
    type Item = (f64, T);

    fn next(&mut self) -> Option<Self::Item> {
        let (index, value) = self.inner.next()?;
        let key = self.key(index);
        Some((key, value))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

inheriting_iter! { @exact IntoIter<T>, T }
inheriting_iter! { @fused IntoIter<T>, T }

impl<T> DoubleEndedIterator for IntoIter<T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let (index, value) = self.inner.next_back()?;
        let key = self.key(index);
        Some((key, value))
    }
}

impl<T> IntoIterator for Hist<T> {
    type Item = (f64, T);

    type IntoIter = IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        let inner = self.values.into_iter().enumerate();
        IntoIter {
            inner,
            min: self.min,
            step: self.step,
        }
    }
}

/// An [`rkyv`] compatible version of [`Hist`].
pub struct ArchivedHist<T: Archive> {
    values: Archived<Vec<T>>,
    min: Archived<f64>,
    max: Archived<f64>,
    step: Archived<f64>,
}

impl<T: Archive> Archive for Hist<T> {
    type Archived = ArchivedHist<T>;

    type Resolver = VecResolver;

    unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
        let (fp, fo) = out_field!(out.values);
        ArchivedVec::resolve_from_len(self.values.len(), pos + fp, resolver, fo);

        let (fp, fo) = out_field!(out.min);
        self.min.resolve(pos + fp, (), fo);

        let (fp, fo) = out_field!(out.max);
        self.max.resolve(pos + fp, (), fo);

        let (fp, fo) = out_field!(out.step);
        self.step.resolve(pos + fp, (), fo);
    }
}

impl<S, T> rkyv::Serialize<S> for Hist<T>
where
    S: rkyv::ser::Serializer + ScratchSpace + ?Sized,
    T: Archive + rkyv::Serialize<S>,
{
    fn serialize(
        &self,
        serializer: &mut S,
    ) -> Result<Self::Resolver, <S as rkyv::Fallible>::Error> {
        self.values.serialize(serializer)
    }
}

impl<T, D> rkyv::Deserialize<Hist<T>, D> for ArchivedHist<T>
where
    T: Archive,
    T::Archived: rkyv::Deserialize<T, D>,
    D: Fallible + ?Sized,
{
    fn deserialize(&self, deserializer: &mut D) -> Result<Hist<T>, <D as Fallible>::Error> {
        Ok(Hist {
            values: self.values.deserialize(deserializer)?,
            min: self.min.deserialize(deserializer)?,
            max: self.max.deserialize(deserializer)?,
            step: self.step.deserialize(deserializer)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contruction_bucket_count() {
        // 0.0 bucket and 1.0 bucket
        let hist = Hist::<()>::new(0.0, 1.0, 1.0);
        assert_eq!(hist.len(), 2);

        // 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 buckets
        let hist = Hist::<()>::new(1.0, 15.0, 1.0);
        assert_eq!(hist.len(), 15);

        let hist = Hist::<()>::new(0.0, 15.0, 1.0);
        assert_eq!(hist.len(), 16);

        let hist = Hist::<()>::new(-1.0, 20.0, 1.0);
        assert_eq!(hist.len(), 22);

        let hist = Hist::<()>::new(0.0, 10.0, 2.0);
        assert_eq!(hist.len(), 6);

        let hist = Hist::<()>::new(0.0, 9.0, 2.0);
        assert_eq!(hist.len(), 6);

        let hist = Hist::<()>::new(0.0, 8.0, 2.0);
        assert_eq!(hist.len(), 5);
    }

    #[test]
    fn get_boundaries() {
        let mut hist = Hist::<usize>::new(1.0, 15.0, 1.0);
        *hist.get_mut(1.0).unwrap() = 100;
        *hist.get_mut(15.0).unwrap() = 100;

        assert_eq!(hist.get(1.0).unwrap(), &100);
        assert_eq!(hist.values[0], 100);
        assert_eq!(hist.get(15.0).unwrap(), &100);
        assert_eq!(hist.values[14], 100);
        for i in 2..15 {
            assert_eq!(hist.get(i as f64).unwrap(), &0);
            assert_eq!(hist.values[i - 1], 0);
        }
    }

    #[test]
    fn get_out_of_bounds() {
        let hist = Hist::<()>::new(1.0, 15.0, 1.0);
        assert!(hist.get(0.99).is_none());
        assert!(hist.get(1.0).is_some());
        assert!(hist.get(-1.0).is_none());
        assert!(hist.get(15.001).is_none());
        assert!(hist.get(15.0).is_some());

        let hist = Hist::<()>::new(-1.0, 8.0, 1.5);
        assert!(hist.get(-1.1).is_none());
        assert!(hist.get(-1.0).is_some());
        assert!(hist.get(8.0).is_some());
        assert!(hist.get(8.1).is_none());
        assert!(hist.get(4.0).is_some());
    }

    #[test]
    fn iterators() {
        let mut hist = Hist::<usize>::new(1.0, 15.0, 1.0);
        hist.insert(2.0, 1);
        hist.insert(5.0, 4);
        hist.insert(3.0, 2);
        hist.insert(4.0, 3);

        // 4 elements, 15 buckets, 11 zeros
        assert_eq!(hist.values().len(), 15);
        assert_eq!(
            hist.values().collect::<Vec<&usize>>(),
            vec![&0, &1, &2, &3, &4, &0, &0, &0, &0, &0, &0, &0, &0, &0, &0]
        );
        assert_eq!(
            hist.clone().into_values().collect::<Vec<usize>>(),
            vec![0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );
        assert_eq!(hist.keys().len(), 15);
        assert_eq!(
            hist.keys().collect::<Vec<f64>>(),
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
        );
        assert_eq!(
            hist.clone().into_keys().collect::<Vec<f64>>(),
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
        );
        assert_eq!(hist.iter().len(), 15);
        assert_eq!(
            hist.iter().collect::<Vec<(f64, &usize)>>(),
            vec![
                (1.0, &0),
                (2.0, &1),
                (3.0, &2),
                (4.0, &3),
                (5.0, &4),
                (6.0, &0),
                (7.0, &0),
                (8.0, &0),
                (9.0, &0),
                (10.0, &0),
                (11.0, &0),
                (12.0, &0),
                (13.0, &0),
                (14.0, &0),
                (15.0, &0),
            ]
        );
        assert_eq!(
            hist.clone().into_iter().collect::<Vec<(f64, usize)>>(),
            vec![
                (1.0, 0),
                (2.0, 1),
                (3.0, 2),
                (4.0, 3),
                (5.0, 4),
                (6.0, 0),
                (7.0, 0),
                (8.0, 0),
                (9.0, 0),
                (10.0, 0),
                (11.0, 0),
                (12.0, 0),
                (13.0, 0),
                (14.0, 0),
                (15.0, 0),
            ]
        );
    }
}
