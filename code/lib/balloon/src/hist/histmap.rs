//! A wrapper for [`SaturatingHist`] with a map-like interface similar to
//! [`FuzzyMap`](crate::fuzzy::map::FuzzyMap).

use std::ops::Range;

use rkyv::{
    option::ArchivedOption,
    out_field,
    ser::{ScratchSpace, Serializer},
    Archive, Archived, Deserialize, Fallible, Serialize,
};

use super::saturating::SaturatingHist;

/// The struct for a [`HistMap`].
///
/// A [`HistMap`] wraps a [`SaturatingHist`] and holds a tolerance. The tolerance determines the bucket intervals
/// of the underlying [`SaturatingHist`], allowing fuzzy matching of [`f64`] keys.
///
/// # Tolerance
/// The bucket intervals will be the same as the tolerance. For example, a tolerance of 0.01 on a [`HistMap`] with
/// lower bound of 0.0 will result in buckets with centers at 0.00, 0.01, 0.02, ...
///
/// Unlike [`crate::fuzzy::map::FuzzyMap`], the buckets include bottom bound but not the top bound. For example, a
/// bucket at 1.01 in a HistMap with tolerance 0.01 will contain 1.00 but not 1.02.
///
/// # Order-insensitivity
/// Unlike [`crate::fuzzy::map::FuzzyMap`], a [`HistMap`] is not order-sensitive and
/// therefore allows self-consistency on fuzzy key matching regardless of the order of insertion.
///
/// For example, Suppose a [`crate::fuzzy::map::FuzzyMap`] is created with tolerance 1.0. If (1.2, 1) and (1.3, 2)
/// are inserted in that order, the final bucket will have range [0.2, 2.2], whereas if (1.3, 2) is inserted before
/// (1.2, 1), the final bucket will have range [0.3, 2.3]. Therefore, [`crate::fuzzy::map::FuzzyMap`] is
/// order-sensitive.
///
/// On the other hand, [`HistMap`] is not order-sensitive. Suppose a [`HistMap`] is created with tolerance 1.0. Then
/// all the bucket centers will be pre-defined. In the case where the lower bound of the [`HistMap`] is 0.0, the bucket
/// centers will be 0.0, 1.0, 2.0, ... Therefore, (1.2, 1) and (1.3, 2) will both be in the bucket with range
/// [0.5, 1.5) regardless of insertion order.
///
/// # Polymorphism
/// HistMap allows for a generic type `T` to be stored in each bucket. The underlying [`SaturatingHist`] stores a
/// [`f64`] as the key and a tuple of [`f64`] and `Option<T>` as the value. The [`f64`] in the tuple is the center
/// of the bucket and the `Option<T>` is the value stored in the bucket. The `Option<T>` is `None` if no value has been
/// inserted into the bucket, which corresponds to an inexistent key in a map. The `Option<T>` is `Some(T)` if a value
/// has been inserted into the bucket, which corresponds to an existent key-value pair in a map.
#[derive(Debug, Clone)]
pub struct HistMap<T> {
    pub(crate) hist: SaturatingHist<(f64, Option<T>)>,
}

impl<T: Clone> HistMap<T> {
    /// Creates a new [`HistMap`] like [`HistMap::new`] but with a specified range.
    ///
    /// The bucket intervals will be the same as the tolerance. The `min` and `max` arguments specify the minimum and
    /// maximum bucket centers.
    pub fn new(min: f64, max: f64, step: f64) -> Self {
        let mut hist = SaturatingHist::new(min, max, step);
        for i in 0..hist.len() {
            // .unwrap() safe since we're guaranteeing valid indices with the loop
            let center = hist.key(i).unwrap();
            hist.insert(center, (center, None));
        }
        Self { hist }
    }
}

impl<T> HistMap<T> {
    /// Get the bucket index for a given key.
    ///
    /// This will return the same result as getting the index of the underlying histogram via
    /// [`SaturatingHist::index`].
    pub fn index(&self, key: f64) -> usize {
        self.hist.index(key)
    }

    /// Get the center of the bucket corresponding to the given key.
    ///
    /// Note that keys that are out of range will be clamped to the nearest bucket center.
    pub fn bucket_center(&self, key: f64) -> f64 {
        self.hist.get(key).0
    }

    /// Get the value in the bucket associated with this key, if any.
    ///
    /// Note that since values are placed into buckets, the key may not be exactly the same as the key used to insert
    /// the value. Therefore, only values within the corresponding bucket will be returned, rather than values inserted
    /// with any key within the tolerance.
    ///
    /// # Example
    /// ```
    /// # use balloon::{fuzzy::Tol, hist::histmap::HistMap, of64};
    /// let mut hmap = HistMap::<usize>::new(0.0, 3000.0, 1.0);
    ///
    /// hmap.insert(1.0, 1);
    /// assert_eq!(hmap.get(1.0), Some(&1));
    /// assert_eq!(hmap.get(0.5), Some(&1));
    /// assert_eq!(hmap.get(0.9), Some(&1));
    /// assert_eq!(hmap.get(1.499), Some(&1));
    /// ```
    pub fn get(&self, key: f64) -> Option<&T> {
        self.hist.get(key).1.as_ref()
    }

    /// Get the bucket center, which serves as a key, and value in the bucket associated with this key, if any.
    ///
    /// Note that the bucket that is being checked has same caveat as [`HistMap::get`].
    pub fn get_entry(&self, key: f64) -> &(f64, Option<T>) {
        self.hist.get(key)
    }

    /// Get a mutable reference to the value in the bucket associated with this key, if any.
    ///
    /// Note that the bucket that is being checked has same caveat as [`HistMap::get`].
    pub fn get_mut(&mut self, key: f64) -> Option<&mut T> {
        self.hist.get_mut(key).1.as_mut()
    }

    /// Check if the map has something in the bucket associated with the given key.
    ///
    /// Note that the bucket that is being checked has same caveat as [`HistMap::get`].
    pub fn contains_key(&self, key: f64) -> bool {
        self.get(key).is_some()
    }

    /// Get all values associated with with this key.
    ///
    /// For a [`HistMap`], this would be the bucket for the key and its two neighboring buckets, in which entries inserted
    /// with keys within the tolerance of the bucket center could be present.
    pub fn get_all_relevant(&self, key: f64) -> impl Iterator<Item = &T> {
        self.get_all_relevant_keys(key)
            .flat_map(move |k| self.get(*k))
    }

    /// Get an iterator over all bucket centers that could be associated with this key.
    ///
    /// This would include the closest bucket center and its two neighboring buckets, if they exist. For example, if the
    /// tolerance is 1.0 and the [`HistMap`] starts at 0.0, then calling this with 1.2 would return the iterator
    /// [0.0, 1.0, 2.0], while calling this with 0.1 will return [0.0, 1.0].
    pub fn get_all_relevant_keys(&self, key: f64) -> impl Iterator<Item = &f64> + '_ {
        let center_i = self.index(key) as i32;
        [center_i - 1, center_i, center_i + 1]
            .into_iter()
            .filter(|i| i >= &0 && i < &(self.hist.len() as i32))
            .filter(|i| self.hist.hist.values[*i as usize].1.is_some())
            .map(|i| &self.hist.hist.values[i as usize].0)
    }

    /// Get the bucket center and value associated with this key, if any.
    ///
    /// This access the same bucket as [`HistMap::get`], but returns a reference to the bucket center as well.
    pub fn get_key_value(&self, key: f64) -> Option<(&f64, &T)> {
        match self.get_entry(key) {
            (k, Some(v)) => Some((k, v)),
            _ => None,
        }
    }

    /// Get the bucket center and value associated with this key, if any.
    ///
    /// This access the same bucket as [`HistMap::get`], but returns the bucket center and a mutable references to
    /// the value.
    pub fn get_key_value_mut(&mut self, key: f64) -> Option<(f64, &mut T)> {
        let center = self.bucket_center(key);
        self.get_mut(key).map(|v| (center, v))
    }

    /// Get all keys and values that could match with the current key.
    ///
    /// This would return entries the bucket for the given key and its two neighboring buckets, which are guaranteed
    /// to contain all potential entries inserted with keys within the tolerance of the given key (though not all
    /// returned entries returned are guaranteed to be within the tolerance of given key).
    pub fn get_all_relevant_key_value(&self, key: f64) -> impl Iterator<Item = (&f64, &T)> {
        self.get_all_relevant_keys(key)
            .zip(self.get_all_relevant(key))
    }

    /// Insert a value into the map with the given key.
    ///
    /// If a value already exists in the bucket associated with the key, it will be replaced and returned. Otherwise,
    /// `None` will be returned
    pub fn insert(&mut self, key: f64, value: T) -> Option<T> {
        let (_, old_value) = self
            .hist
            .insert(key, (self.bucket_center(key), Some(value)));
        old_value
    }

    // TODO into_keys, into_values omitted

    /// Get an iterator over all key-value pairs in the map.
    pub fn iter(&self) -> Iter<'_, T> {
        Iter {
            inner: self.hist.hist.values.iter(),
        }
    }

    /// Get an iterator over all keys in the map.
    pub fn keys(&self) -> Keys<'_, T> {
        Keys {
            hist: self,
            inner: 0..self.hist_len(),
        }
    }

    /// Get an iterator over all values in the map.
    pub fn values(&self) -> impl Iterator<Item = &T> + '_ {
        self.hist
            .values()
            .filter(|(_, v)| v.is_some())
            .map(|(_, v)| v.as_ref().unwrap())
    }

    /// Get the number of buckets in the histogram.
    pub fn hist_len(&self) -> usize {
        self.hist.len()
    }

    /// Get the number of entries in the map.
    pub fn len(&self) -> usize {
        self.hist.values().filter(|(_, v)| v.is_some()).count()
    }

    /// Check if all there is no value inserted into the map.
    pub fn is_empty(&self) -> bool {
        self.hist
            .hist
            .values
            .iter()
            .all(|(_, bucket)| bucket.is_none())
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

/// An iterator over all floating point bucket keys in a [`HistMap`].
#[derive(Debug, Clone)]
pub struct Keys<'a, T> {
    hist: &'a HistMap<T>,
    inner: Range<usize>,
}

impl<'a, T> Iterator for Keys<'a, T> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        for index in &mut self.inner {
            let (key, value) = &self.hist.hist.hist.values[index];
            if value.is_some() {
                return Some(*key);
            }
        }
        None
    }
}
inheriting_iter! { @fused Keys<'a, T>, 'a, T }

impl<'a, T> DoubleEndedIterator for Keys<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        while let Some(index) = self.inner.next_back() {
            let (key, value) = &self.hist.hist.hist.values[index];
            if value.is_some() {
                return Some(*key);
            }
        }
        None
    }
}

/// An iterator over references to the key-value pairs of a [`HistMap`].
#[derive(Debug, Clone)]
pub struct Iter<'a, T> {
    inner: std::slice::Iter<'a, (f64, Option<T>)>,
}

impl<'a, T> Iterator for Iter<'a, T> {
    type Item = (f64, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        for (key, value) in &mut self.inner {
            if let Some(value) = value.as_ref() {
                return Some((*key, value));
            }
        }
        None
    }
}
inheriting_iter! { @exact Iter<'a, T>, 'a, T }
inheriting_iter! { @fused Iter<'a, T>, 'a, T }

/// An iterator over copies of the key-value pairs of a [`HistMap`].
#[derive(Debug, Clone)]
pub struct IntoIter<T> {
    inner: std::vec::IntoIter<(f64, Option<T>)>,
}

impl<T> Iterator for IntoIter<T> {
    type Item = (f64, T);

    fn next(&mut self) -> Option<Self::Item> {
        for (key, value) in &mut self.inner {
            if let Some(value) = value {
                return Some((key, value));
            }
        }
        None
    }
}
inheriting_iter! { @exact IntoIter<T>, T }
inheriting_iter! { @fused IntoIter<T>, T }

impl<T> IntoIterator for HistMap<T> {
    type Item = (f64, T);

    type IntoIter = IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIter {
            inner: self.hist.hist.values.into_iter(),
        }
    }
}

type HistMapHist<T> = SaturatingHist<(f64, Option<T>)>;

/// An archived version of [`HistMap`].
pub struct ArchivedHistMap<T: Archive> {
    hist: Archived<HistMapHist<T>>,
}

impl<T: Archive> Archive for HistMap<T> {
    type Archived = ArchivedHistMap<T>;
    type Resolver = <HistMapHist<T> as Archive>::Resolver;

    unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
        let (fp, fo) = out_field!(out.hist);
        self.hist.resolve(pos + fp, resolver, fo);
    }
}

impl<S, T> Serialize<S> for HistMap<T>
where
    S: ScratchSpace + Serializer + ?Sized,
    T: Archive + Serialize<S>,
{
    fn serialize(&self, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        self.hist.serialize(serializer)
    }
}

impl<T, D> Deserialize<HistMap<T>, D> for ArchivedHistMap<T>
where
    T: Archive,
    T::Archived: Deserialize<T, D>,
    D: Fallible + ?Sized,
{
    fn deserialize(&self, deserializer: &mut D) -> Result<HistMap<T>, D::Error> {
        let hist = self.hist.deserialize(deserializer)?;
        Ok(HistMap { hist })
    }
}

impl<T: Archive> ArchivedHistMap<T> {
    /// Returns the archived minimum bucket of the histogram.
    pub fn min(&self) -> Archived<f64> {
        self.hist.hist.min
    }

    /// Returns the archived maximum bucket of the histogram.
    pub fn max(&self) -> Archived<f64> {
        self.hist.hist.max
    }

    /// Returns the archived step size of the histogram.
    pub fn step(&self) -> Archived<f64> {
        self.hist.hist.step
    }

    /// get index in array for given key for an archived histmap.
    ///
    /// This is not public because the vectors are not exposed.
    fn index(&self, key: f64) -> usize {
        if key < self.min() {
            0
        } else if key > self.max() {
            self.hist.hist.values.len() - 1
        } else {
            ((key - self.min()) / self.step()).round() as usize
        }
    }

    /// Returns an iterator over all archived key-value pairs in the map.
    pub fn get_all_key_value(
        &self,
        key: Archived<f64>,
    ) -> impl Iterator<Item = (&Archived<f64>, &<T as Archive>::Archived)> {
        let center = self.index(key);
        let lower = center.saturating_sub(1);
        let upper = (center + 1).min(self.hist.hist.values.len() - 1);

        (lower..=upper).filter_map(|index| {
            let (key, value) = &self.hist.hist.values[index];
            match value {
                ArchivedOption::None => None,
                ArchivedOption::Some(value) => Some((key, value)),
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{of64, ArchivedOf64};

    #[test]
    fn contruction_bucket_count() {
        let f = |x: f64| f64::try_from(x).unwrap();

        let hmap = HistMap::<()>::new(0.0, 3000.0, 1.0);
        assert_eq!(hmap.hist_len(), 3001);
        assert_eq!(hmap.hist.len(), 3001);
        assert_eq!(hmap.hist.hist.values[0].0, f(0.0));
        assert_eq!(hmap.hist.hist.values[1].0, f(1.0));
        assert_eq!(hmap.hist.hist.values[2].0, f(2.0));

        let hmap = HistMap::<()>::new(0.0, 3000.0, 5.0);
        assert_eq!(hmap.hist_len(), 601);
        assert_eq!(hmap.hist.len(), 601);

        let hmap = HistMap::<()>::new(0.0, 3000.0, 1001.0);
        assert_eq!(hmap.hist_len(), 4);
        assert_eq!(hmap.hist.len(), 4);
    }

    #[test]
    fn test_index() {
        let hmap = HistMap::<()>::new(0.0, 3000.0, 1.0);
        assert_eq!(hmap.index(0.1), 0);
        assert_eq!(hmap.index(0.9), 1);
        assert_eq!(hmap.index(1.0), 1);
        assert_eq!(hmap.index(1.1), 1);
        assert_eq!(hmap.index(1.5), 2);
        assert_eq!(hmap.index(1.9), 2);
        assert_eq!(hmap.index(120.344), 120);
    }

    #[test]
    fn test_bucket_center() {
        let hmap = HistMap::<()>::new(0.0, 3000.0, 1.0);
        assert_eq!(hmap.bucket_center(0.1), 0.0);
        assert_eq!(hmap.bucket_center(0.9), 1.0);
        assert_eq!(hmap.bucket_center(1.0), 1.0);
        assert_eq!(hmap.bucket_center(1.1), 1.0);
        assert_eq!(hmap.bucket_center(1.9), 2.0);
        assert_eq!(hmap.bucket_center(120.344), 120.0);
    }

    #[test]
    fn test_get() {
        let mut hmap = HistMap::<usize>::new(0.0, 3000.0, 1.0);
        hmap.insert(0.2, 1);
        hmap.insert(1.0, 2);
        hmap.insert(2.87, 3);
        hmap.insert(187.5, 4);
        hmap.insert(3000.499, 5);

        assert_eq!(hmap.get(0.3).unwrap(), &1);
        assert_eq!(hmap.get_key_value(0.3).unwrap(), (&0.0, &1));

        assert_eq!(hmap.get(0.0).unwrap(), &1);
        assert_eq!(hmap.get_key_value(0.3).unwrap(), (&0.0, &1));

        // one tolerance below min should still be accepted.
        assert_eq!(hmap.get(-0.4).unwrap(), &1);
        assert_eq!(hmap.get_key_value(0.3).unwrap(), (&0.0, &1));

        assert_eq!(hmap.get(0.9).unwrap(), &2);
        assert_eq!(hmap.get_key_value(0.9).unwrap(), (&1.0, &2));

        assert_eq!(hmap.get(1.3).unwrap(), &2);
        assert_eq!(hmap.get_key_value(1.3).unwrap(), (&1.0, &2));

        assert_eq!(hmap.get(2.9).unwrap(), &3);
        assert_eq!(hmap.get_key_value(2.9).unwrap(), (&3.0, &3));

        assert_eq!(hmap.get(3.11).unwrap(), &3);
        assert_eq!(hmap.get_key_value(3.11).unwrap(), (&3.0, &3));

        assert_eq!(hmap.get(187.999).unwrap(), &4);
        assert_eq!(hmap.get_key_value(187.999).unwrap(), (&188.0, &4));

        // one tolerance above max should still be accepted.
        assert_eq!(hmap.get(3000.4).unwrap(), &5);
        assert_eq!(hmap.get_key_value(3000.4).unwrap(), (&3000.0, &5));
    }

    #[test]
    fn test_archive() {
        let f = |x: f64| of64::try_from(x).unwrap();
        let af = |x: f64| ArchivedOf64::from(of64::try_from(x).unwrap());

        let mut hmap = HistMap::<of64>::new(0.0, 3000.0, 1.0);
        hmap.insert(0.2, f(1.0));
        hmap.insert(1.0, f(2.0));
        hmap.insert(2.87, f(3.0));
        hmap.insert(187.5, f(4.0));

        let serialized = rkyv::to_bytes::<_, 256>(&hmap).unwrap();
        let archived = unsafe { rkyv::archived_root::<HistMap<of64>>(&serialized[..]) };
        let deserialized: HistMap<of64> = archived.deserialize(&mut rkyv::Infallible).unwrap();

        assert_eq!(deserialized.hist.hist.min, hmap.hist.hist.min);
        assert_eq!(deserialized.hist.hist.max, hmap.hist.hist.max);
        assert_eq!(
            deserialized.hist.hist.values.len(),
            hmap.hist.hist.values.len()
        );

        for (a, b) in deserialized
            .hist
            .hist
            .values
            .iter()
            .zip(hmap.hist.hist.values.iter())
        {
            assert_eq!(a, b);
        }

        for i in 0..100 {
            for ((ak, av), (bk, bv)) in deserialized
                .get_all_relevant_key_value(0.4 * i as f64)
                .zip(archived.get_all_key_value((0.4 * i as f64) as Archived<f64>))
            {
                assert_eq!(&(*ak as Archived<f64>), bk, "failed k on {}", i);
                assert_eq!(&af(av.into_inner()), bv, "failed v on {}", i);
            }
        }

        assert_eq!(deserialized.get(0.3).unwrap(), &f(1.0));
        assert_eq!(deserialized.get_key_value(0.3).unwrap(), (&0.0, &f(1.0)));

        assert_eq!(deserialized.get(0.0).unwrap(), &f(1.0));
        assert_eq!(deserialized.get_key_value(0.3).unwrap(), (&0.0, &f(1.0)));

        assert_eq!(deserialized.get(0.9).unwrap(), &f(2.0));
        assert_eq!(deserialized.get_key_value(0.9).unwrap(), (&1.0, &f(2.0)));

        assert_eq!(deserialized.get(1.3).unwrap(), &f(2.0));
        assert_eq!(deserialized.get_key_value(1.3).unwrap(), (&1.0, &f(2.0)));

        assert_eq!(deserialized.get(2.9).unwrap(), &f(3.0));
        assert_eq!(deserialized.get_key_value(2.9).unwrap(), (&3.0, &f(3.0)));

        assert_eq!(deserialized.get(3.11).unwrap(), &f(3.0));
        assert_eq!(deserialized.get_key_value(3.11).unwrap(), (&3.0, &f(3.0)));

        assert_eq!(deserialized.get(187.999).unwrap(), &f(4.0));
        assert_eq!(
            deserialized.get_key_value(187.999).unwrap(),
            (&188.0, &f(4.0))
        );
    }

    #[test]
    fn test_get_mirror_fuzzymap() {
        let mut hmap = HistMap::<usize>::new(0.0, 3000.0, 2.0);

        hmap.insert(2.0, 10);
        hmap.insert(2.5, 12);
        hmap.insert(3.0, 15);
        hmap.insert(4.1, 17);

        assert_eq!(hmap.get(1.8).unwrap(), &12);
        assert_eq!(hmap.get(2.2).unwrap(), &12);
        assert_eq!(hmap.get(4.9).unwrap(), &17);
        assert_eq!(hmap.get(6.3), None);
        assert_eq!(hmap.get(-0.3), None);
    }

    #[test]
    fn test_get_all() {
        let mut hmap = HistMap::<usize>::new(0.0, 3000.0, 1.0);
        hmap.insert(0.2, 1);
        hmap.insert(1.0, 2);
        hmap.insert(2.1, 3);
        hmap.insert(2.87, 4);
        hmap.insert(187.5, 5);
        hmap.insert(3000.0, 6); // max bucket

        assert_eq!(
            hmap.get_all_relevant(1.1).collect::<Vec<&usize>>(),
            vec![&1, &2, &3]
        );
        assert_eq!(
            hmap.get_all_relevant_keys(1.1).collect::<Vec<&f64>>(),
            vec![&0.0, &1.0, &2.0]
        );
        assert_eq!(
            hmap.get_all_relevant_key_value(1.1)
                .collect::<Vec<(&f64, &usize)>>(),
            vec![(&0.0, &1), (&1.0, &2), (&2.0, &3)]
        );

        assert_eq!(
            hmap.get_all_relevant(0.3).collect::<Vec<&usize>>(),
            vec![&1, &2]
        );
        assert_eq!(
            hmap.get_all_relevant_keys(0.3).collect::<Vec<&f64>>(),
            vec![&0.0, &1.0]
        );
        assert_eq!(
            hmap.get_all_relevant_key_value(0.3)
                .collect::<Vec<(&f64, &usize)>>(),
            vec![(&0.0, &1), (&1.0, &2)]
        );

        assert_eq!(
            hmap.get_all_relevant(187.6).collect::<Vec<&usize>>(),
            vec![&5]
        );
        assert_eq!(
            hmap.get_all_relevant_keys(187.6).collect::<Vec<&f64>>(),
            vec![&188.0]
        );
        assert_eq!(
            hmap.get_all_relevant_key_value(187.6)
                .collect::<Vec<(&f64, &usize)>>(),
            vec![(&188.0, &5)]
        );

        assert_eq!(hmap.get_all_relevant(18.0).collect::<Vec<&usize>>(), {
            let target: Vec<&usize> = vec![];
            target
        });
        assert_eq!(
            hmap.get_all_relevant(2999.9).collect::<Vec<&usize>>(),
            vec![&6]
        );
        assert_eq!(
            hmap.get_all_relevant(3000.0).collect::<Vec<&usize>>(),
            vec![&6]
        );
        assert_eq!(
            hmap.get_all_relevant(3000.5).collect::<Vec<&usize>>(),
            vec![&6]
        );
    }
}
