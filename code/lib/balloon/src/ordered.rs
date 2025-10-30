use rkyv::{out_field, ser::Serializer, Archive, Archived, Deserialize, Resolver, Serialize};
use std::{
    fmt, hash,
    ops::{Add, Mul, Sub},
};

// TODO: revisit mathematical op traits
// one option is to just define nothing and force people to go through the std f64 and reconstruct
// a of64. Slightly less ergonomic for users, but much easier to enforce safety without a ton of
// runtime checks

/// A floating point type that cannot be NaN.
///
/// This allows us to have a full instead of partial equivalence relationship. It also gives us the
/// ability to compute orderings and hash.
#[repr(transparent)]
#[allow(non_camel_case_types)]
#[derive(Debug, Deserialize, Clone, Copy)]
pub struct of64(f64);

impl Default for of64 {
    fn default() -> of64 {
        of64(0.0)
    }
}

impl std::fmt::Display for of64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<of64> for f64 {
    fn from(f: of64) -> Self {
        f.0
    }
}

impl From<&of64> for f64 {
    fn from(f: &of64) -> Self {
        f.0
    }
}

impl TryFrom<f64> for of64 {
    type Error = NanError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl TryFrom<&f64> for of64 {
    type Error = NanError;

    fn try_from(value: &f64) -> Result<Self, Self::Error> {
        Self::new(*value)
    }
}

impl AsRef<f64> for of64 {
    fn as_ref(&self) -> &f64 {
        &self.0
    }
}

// TODO: * can produces NaNs even if no inputs are NaN
// http://www.hlam.ece.ufl.edu/EEL4712/Labs/Lab6/IEEEStandard754FP.pdf
// ± inf * 0 := NaN
impl Mul for of64 {
    type Output = of64;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul<&of64> for of64 {
    type Output = of64;

    fn mul(self, rhs: &of64) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul for &of64 {
    type Output = of64;

    fn mul(self, rhs: Self) -> Self::Output {
        of64(self.0 * rhs.0)
    }
}

impl Mul<of64> for &of64 {
    type Output = of64;

    fn mul(self, rhs: of64) -> Self::Output {
        of64(self.0 * rhs.0)
    }
}

// TODO: - can produce NaNs even if no inputs are NaN
// see link above
// inf - inf := NaN
impl Sub for of64 {
    type Output = of64;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Add for of64 {
    type Output = of64;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

/// An error for when a [`of64`] is constructed with a [`f64::NAN`].
#[derive(Debug)]
pub struct NanError;

impl fmt::Display for NanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "attempted to create of64 with a NaN")
    }
}

impl std::error::Error for NanError {}

macro_rules! of64_const {
    ($name:ident) => {
        #[doc = concat!("See [`f64::", stringify!($name), "`]")]
        pub const $name: Self = unsafe { Self::new_unchecked(f64::$name) };
    };
}

impl of64 {
    of64_const! { EPSILON }
    of64_const! { MIN }
    of64_const! { MIN_POSITIVE }
    of64_const! { MAX }

    /// Attempt to convert a [`f64`] into a [`of64`].
    ///
    /// # Errors
    /// If `inner` is [`f64::NAN`] this will return an error.
    pub fn new(inner: f64) -> Result<Self, NanError> {
        if inner.is_nan() {
            return Err(NanError);
        }

        Ok(Self(inner))
    }

    /// Convert a [`f64`] into a [`of64`] without checking for NaNs.
    ///
    /// # Safety
    /// - `inner` must not be NaN
    pub const unsafe fn new_unchecked(inner: f64) -> Self {
        Self(inner)
    }

    /// Returned the wrapped `f64`.
    pub fn into_inner(self) -> f64 {
        self.0
    }

    /// Compute the absolute value of this ordered float.
    pub fn abs(&self) -> of64 {
        Self(self.0.abs())
    }
}

impl PartialEq for of64 {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for of64 {}

impl PartialOrd for of64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for of64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // SAFETY: this can only be unsafe if the user uses an unsafe constructor, otherwise all
        // constructors checked that the inner value is not NaN, so the partial comparison should
        // never be None
        unsafe { self.partial_cmp(other).unwrap_unchecked() }
    }
}

impl hash::Hash for of64 {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        self.0.to_bits().hash(state);
    }
}

/// [`ArchivedOf64`] implements [`Archive`] for the ordered floating point type [`of64`]
///
/// This type can be created by serializing [`of64`] and casting the resulting byte buffer.
#[derive(Debug, Clone, Copy)]
pub struct ArchivedOf64(Archived<f64>);

// This is a convenience function to archive an ordered float without going through rkyv
impl From<of64> for ArchivedOf64 {
    fn from(f: of64) -> ArchivedOf64 {
        ArchivedOf64(f.0)
    }
}

impl TryFrom<f64> for ArchivedOf64 {
    type Error = NanError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        of64::new(value).map(|x| x.into())
    }
}

impl ArchivedOf64 {
    /// This function converts the inner floating point to its absolute value.
    pub fn abs(&self) -> ArchivedOf64 {
        Self(self.0.abs())
    }

    /// This function returns the inner float
    pub fn inner(&self) -> f64 {
        self.0
    }
}

impl Archive for of64 {
    type Archived = ArchivedOf64;
    type Resolver = Resolver<f64>;

    unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
        let (fp, fo) = out_field!(out.0);
        self.0.resolve(pos + fp, resolver, fo)
    }
}

impl<S: Serializer + ?Sized> Serialize<S> for of64 {
    fn serialize(&self, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        self.0.serialize(serializer)
    }
}

impl Ord for ArchivedOf64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // SAFETY: this can only be unsafe if the user uses an unsafe constructor, otherwise all
        // constructors checked that the inner value is not NaN, so the partial comparison should
        // never be None
        unsafe { self.partial_cmp(other).unwrap_unchecked() }
    }
}

impl Eq for ArchivedOf64 {}

impl PartialOrd for ArchivedOf64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl PartialEq for ArchivedOf64 {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Mul for ArchivedOf64 {
    type Output = ArchivedOf64;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul<&ArchivedOf64> for ArchivedOf64 {
    type Output = ArchivedOf64;

    fn mul(self, rhs: &ArchivedOf64) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul for &ArchivedOf64 {
    type Output = ArchivedOf64;

    fn mul(self, rhs: Self) -> Self::Output {
        ArchivedOf64(self.0 * rhs.0)
    }
}

impl Mul<ArchivedOf64> for &ArchivedOf64 {
    type Output = ArchivedOf64;

    fn mul(self, rhs: ArchivedOf64) -> Self::Output {
        ArchivedOf64(self.0 * rhs.0)
    }
}

impl Sub for ArchivedOf64 {
    type Output = ArchivedOf64;

    fn sub(self, rhs: ArchivedOf64) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Add for ArchivedOf64 {
    type Output = ArchivedOf64;

    fn add(self, rhs: ArchivedOf64) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}
