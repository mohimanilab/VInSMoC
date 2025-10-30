//! Features that can be represented as continuous numeric values.

/// A single numeric value.
///
/// A numeric value can either be a scalar or a vector, in either single or double precision.
#[derive(Debug, Clone, PartialEq)]
pub enum NumericValue {
    /// A scalar, single precision numeric value
    F32(f32),
    /// A scalar, double precision numeric value
    F64(f64),
    /// A vector, single precision numeric value
    F32Vec(Vec<f32>),
    /// A vector, double precision numeric value
    F64Vec(Vec<f64>),
}

impl From<f32> for NumericValue {
    fn from(inner: f32) -> Self {
        Self::F32(inner)
    }
}

impl From<f64> for NumericValue {
    fn from(inner: f64) -> Self {
        Self::F64(inner)
    }
}

impl From<Vec<f32>> for NumericValue {
    fn from(inner: Vec<f32>) -> Self {
        Self::F32Vec(inner)
    }
}

impl From<Vec<f64>> for NumericValue {
    fn from(inner: Vec<f64>) -> Self {
        Self::F64Vec(inner)
    }
}

/// The precision of a numeric value.
#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub enum Precision {
    /// Single precision, corresponds to 32-bit floating points.
    Single,
    /// Double precision, corresponds to 64-bit floating points.
    Double,
}

impl NumericValue {
    /// Combine this numeric value with an `other` value.
    ///
    /// First, this matches precision by promoting single precision values into double precision.
    /// Afterwards this returns a new [`NumericValue`] that is guaranteed to be a vector. The
    /// values in `self` are guaranteed to precede the values in `other`.
    ///
    /// # Examples
    /// ```
    /// use learning::data::NumericValue;
    ///
    /// // fusion always preserves order
    /// let one = NumericValue::from(1.0f32);
    /// let two = NumericValue::from(2.0f32);
    /// assert_eq!(one.clone().fuse(two.clone()), NumericValue::from(vec![1.0f32, 2.0]));
    /// assert_eq!(two.clone().fuse(one.clone()), NumericValue::from(vec![2.0f32, 1.0]));
    ///
    /// // if both aren't single precision the lower precision is promoted
    /// let one = NumericValue::from(1.0f32);
    /// let two = NumericValue::from(2.0f64);
    /// assert_eq!(one.clone().fuse(two.clone()), NumericValue::from(vec![1.0f64, 2.0]));
    /// assert_eq!(two.clone().fuse(one.clone()), NumericValue::from(vec![2.0f64, 1.0]));
    /// ```
    pub fn fuse(self, other: Self) -> Self {
        // promote precision if mismatching
        let (this, other) = match (self.precision(), other.precision()) {
            (Precision::Double, Precision::Single) => (self, other.into_double()),
            (Precision::Single, Precision::Double) => (self.into_double(), other),
            _ => (self, other),
        };

        match (this, other) {
            (Self::F32(f1), Self::F32(f2)) => Self::F32Vec(vec![f1, f2]),
            (Self::F64(f1), Self::F64(f2)) => Self::F64Vec(vec![f1, f2]),
            (Self::F32(f), Self::F32Vec(mut v)) => {
                v.insert(0, f);
                Self::F32Vec(v)
            }
            (Self::F32Vec(mut v), Self::F32(f)) => {
                v.push(f);
                Self::F32Vec(v)
            }
            (Self::F64(f), Self::F64Vec(mut v)) => {
                v.insert(0, f);
                Self::F64Vec(v)
            }
            (Self::F64Vec(mut v), Self::F64(f)) => {
                v.push(f);
                Self::F64Vec(v)
            }
            (Self::F32Vec(mut v1), Self::F32Vec(mut v2)) => {
                v1.append(&mut v2);
                Self::F32Vec(v1)
            }
            (Self::F64Vec(mut v1), Self::F64Vec(mut v2)) => {
                v1.append(&mut v2);
                Self::F64Vec(v1)
            }
            _ => unreachable!(),
        }
    }

    /// Convert this numeric value into double precision.
    ///
    /// This is a no-op if the original value was already double precision.
    pub fn into_double(self) -> Self {
        match self {
            Self::F32(f) => Self::F64(f64::from(f)),
            Self::F32Vec(v) => Self::F64Vec(v.into_iter().map(f64::from).collect()),
            x => x,
        }
    }

    /// Convert this numeric value into single precision.
    ///
    /// This is a no-op if the original value was already single precision. Otherwise the double
    /// precision values will be truncated to 32 bits.
    pub fn into_single(self) -> Self {
        match self {
            Self::F64(f) => Self::F32(f as f32),
            Self::F64Vec(v) => Self::F32Vec(v.into_iter().map(|x| x as f32).collect()),
            x => x,
        }
    }

    /// Convert a scalar numeric value into a one element vector.
    ///
    /// This is a no-op if the original value was already a vector.
    pub fn into_vector(self) -> Self {
        match self {
            Self::F32(f) => Self::F32Vec(vec![f]),
            Self::F64(f) => Self::F64Vec(vec![f]),
            _ => self,
        }
    }

    /// Retrieve the precision of this numeric value.
    pub fn precision(&self) -> Precision {
        match self {
            Self::F32(_) | Self::F32Vec(_) => Precision::Single,
            Self::F64(_) | Self::F64Vec(_) => Precision::Double,
        }
    }

    /// Check if this numeric value is a vector.
    pub fn is_vector(&self) -> bool {
        matches!(self, Self::F32Vec(_) | Self::F64Vec(_))
    }

    /// Convert this numeric value into an iterator of [`f32`]s.
    pub fn into_singles(self) -> impl Iterator<Item = f32> {
        match self.into_single() {
            NumericValue::F32(f) => vec![f].into_iter(),
            NumericValue::F32Vec(v) => v.into_iter(),
            _ => unreachable!(),
        }
    }

    /// Convert this numeric value into an iterator of [`f64`]s.
    pub fn into_doubles(self) -> impl Iterator<Item = f64> {
        match self.into_double() {
            NumericValue::F64(f) => vec![f].into_iter(),
            NumericValue::F64Vec(v) => v.into_iter(),
            _ => unreachable!(),
        }
    }
}

impl FromIterator<f32> for NumericValue {
    fn from_iter<T: IntoIterator<Item = f32>>(iter: T) -> Self {
        Self::F32Vec(iter.into_iter().collect())
    }
}

impl FromIterator<f64> for NumericValue {
    fn from_iter<T: IntoIterator<Item = f64>>(iter: T) -> Self {
        Self::F64Vec(iter.into_iter().collect())
    }
}

/// A feature that takes a continuous numeric value, either scalar or vector.
pub trait NumericFeature {
    /// Get the [`NumericValue`] representing the continuous feature value.
    fn value(&self) -> NumericValue;
}
