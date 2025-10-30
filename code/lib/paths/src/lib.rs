//! Utility types to represent particular kinds of paths on a filesystem.
//!
//! This crate provides thin wrappers around [`PathBuf`](std::path::PathBuf) that ensure upon
//! construction that
//!
//! - [The path corresponds to an existing file](crate::FilePath)
//! - [The path corresponds to an new file](crate::NewFilePath)
//! - [The path corresponds to an existing directory](crate::DirPath)
//! - [The path corresponds to a new directory](crate::NewDirPath)
//!
//! # Trait Implementations
//! All the provided types implement converters from the standard library path types. During
//! construction the invariants for each type are enforced, returning a [`std::io::Error`] on
//! failure.
//!
//! All types implement [`std::str::FromStr`], allowing for `structopt` compatibility.
//!
//! All types implement equality with anything that implements `AsRef<Path>`. Equality is done by
//! comparing paths.
//!
//! ## Note
//! For standard library types equality is only in one direction.
//! ```
//! # use paths::NewFilePath;
//! # use std::{path::Path, convert::TryFrom};
//! // this works
//! NewFilePath::try_from(Path::new("fakefile")).unwrap() == Path::new("fakefile");
//! ```
//! ```compile_fail
//! # use paths::NewFilePath;
//! # use std::{path::Path, convert::TryFrom};
//! // but this doesn't
//! Path::new("fakefile") == NewFilePath::try_from(Path::new("fakefile")).unwrap();
//! ```
//! But for our own types equality is bidirectional.
//! ```
//! # use paths::{NewFilePath, NewDirPath};
//! # use std::{path::Path, str::FromStr};
//! NewFilePath::from_str("fakefile").unwrap() == NewDirPath::from_str("fakefile").unwrap();
//! NewDirPath::from_str("fakefile").unwrap() == NewFilePath::from_str("fakefile").unwrap();
//! ```
//!
//! # Usage
//! See struct docs for usage examples.

#![warn(missing_docs)]

use std::{
    borrow::Borrow,
    convert::TryFrom,
    fs::{DirBuilder, File},
    io::{Error, ErrorKind, Result},
    ops::Deref,
    path::{Path, PathBuf},
    str::FromStr,
};

/// A path that can only be an existing file.
///
/// # Example
/// ```
/// # fn tmp() -> std::io::Result<()> {
/// use tempfile::NamedTempFile;
/// use paths::FilePath;
/// use std::{convert::TryFrom, str::FromStr, path::{Path, PathBuf}};
///
/// // make a file that we know exists
/// let temp_file = NamedTempFile::new()?;
/// let temp_path = temp_file.path();
///
/// // we can create this from a Path
/// let file_path = FilePath::try_from(temp_path)?;
/// assert_eq!(&file_path, &temp_path);
///
/// // but if the file doesn't exist that will give an error
/// let fake_path = "thisreallyshouldneverexistonyourfs";
/// assert!(fake_path.parse::<FilePath>().is_err());
/// assert!(FilePath::try_from(Path::new(fake_path)).is_err());
/// assert!(FilePath::try_from(PathBuf::from(fake_path)).is_err());
/// # Ok(())
/// # }
/// # tmp().unwrap();
/// ```
#[derive(Debug, Clone, Eq, Hash, PartialOrd, Ord, Default)]
pub struct FilePath(PathBuf);

/// A path that can only be a new file.
///
/// # Example
/// ```
/// # fn tmp() -> std::io::Result<()> {
/// use tempfile::NamedTempFile;
/// use paths::NewFilePath;
/// use std::{convert::TryFrom, str::FromStr, path::{Path, PathBuf}};
///
/// // make a file that we know exists
/// let temp_file = NamedTempFile::new()?;
/// let temp_path = temp_file.path();
///
/// // this will fail, the file already exists
/// assert!(NewFilePath::try_from(temp_path).is_err());
///
/// // but if the file doesn't exist then we can create it
/// let fake_path = "thisreallyshouldneverexistonyourfs";
/// assert_eq!(fake_path.parse::<NewFilePath>()?, Path::new(fake_path));
/// let fake_path = Path::new(fake_path);
/// assert_eq!(NewFilePath::try_from(fake_path)?, fake_path);
/// assert_eq!(NewFilePath::try_from(PathBuf::from(fake_path))?, fake_path);
/// # Ok(())
/// # }
/// # tmp().unwrap();
/// ```
#[derive(Debug, Clone, Eq, Hash, PartialOrd, Ord, Default)]
pub struct NewFilePath(PathBuf);

/// A path that can only be an existing directory.
///
/// # Example
/// ```
/// # fn tmp() -> std::io::Result<()> {
/// use tempfile::TempDir;
/// use paths::DirPath;
/// use std::{convert::TryFrom, str::FromStr, path::{Path, PathBuf}};
///
/// // make a directory that we know exists
/// let temp_dir = TempDir::new()?;
/// let temp_path = temp_dir.path();
///
/// // we can create this from a Path
/// let file_path = DirPath::try_from(temp_path)?;
/// assert_eq!(&file_path, &temp_path);
///
/// // but if the dir doesn't exist that will give an error
/// let fake_path = "thisreallyshouldneverexistonyourfs";
/// assert!(fake_path.parse::<DirPath>().is_err());
/// assert!(DirPath::try_from(Path::new(fake_path)).is_err());
/// assert!(DirPath::try_from(PathBuf::from(fake_path)).is_err());
/// # Ok(())
/// # }
/// # tmp().unwrap();
/// ```
#[derive(Debug, Clone, Eq, Hash, PartialOrd, Ord, Default)]
pub struct DirPath(PathBuf);

/// A path that can only be a new directory
///
/// # Example
/// ```
/// # fn tmp() -> std::io::Result<()> {
/// use tempfile::TempDir;
/// use paths::NewDirPath;
/// use std::{convert::TryFrom, str::FromStr, path::{Path, PathBuf}};
///
/// // make a directory that we know exists
/// let temp_dir = TempDir::new()?;
/// let temp_path = temp_dir.path();
///
/// // this will fail, the file already exists
/// assert!(NewDirPath::try_from(temp_path).is_err());
///
/// // but if the dir doesn't exist then we can create it
/// let fake_path = "thisreallyshouldneverexistonyourfs";
/// assert_eq!(fake_path.parse::<NewDirPath>()?, Path::new(fake_path));
/// let fake_path = Path::new(fake_path);
/// assert_eq!(NewDirPath::try_from(fake_path)?, fake_path);
/// assert_eq!(NewDirPath::try_from(PathBuf::from(fake_path))?, fake_path);
/// # Ok(())
/// # }
/// # tmp().unwrap();
/// ```
#[derive(Debug, Clone, Eq, Hash, PartialOrd, Ord, Default)]
pub struct NewDirPath(PathBuf);

macro_rules! shared_traits {
    ($typ:ident) => {
        impl AsRef<Path> for $typ {
            fn as_ref(&self) -> &Path {
                &self.0
            }
        }

        impl AsRef<PathBuf> for $typ {
            fn as_ref(&self) -> &PathBuf {
                &self.0
            }
        }

        impl Deref for $typ {
            type Target = PathBuf;

            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl Borrow<Path> for $typ {
            fn borrow(&self) -> &Path {
                &self.0
            }
        }

        impl TryFrom<&Path> for $typ {
            type Error = Error;

            fn try_from(path: &Path) -> Result<Self> {
                Self::try_from(path.to_path_buf())
            }
        }

        impl FromStr for $typ {
            type Err = Error;

            fn from_str(s: &str) -> Result<Self> {
                let path = PathBuf::from(s);
                Self::try_from(path)
            }
        }

        impl From<$typ> for PathBuf {
            fn from(path: $typ) -> Self {
                path.0
            }
        }

        /// Equality is done by conversion to a [`Path`](std::path::Path).
        impl<T: AsRef<Path>> PartialEq<T> for $typ {
            fn eq(&self, other: &T) -> bool {
                &self.0 == other.as_ref()
            }
        }
    };
}

shared_traits! {FilePath}
shared_traits! {NewFilePath}
shared_traits! {DirPath}
shared_traits! {NewDirPath}

fn check_new_path(path: &Path) -> Result<()> {
    if path.exists() {
        return Err(Error::new(
            ErrorKind::AlreadyExists,
            format!("path {} already exists", path.display()),
        ));
    }

    let parent = path.parent().ok_or_else(|| {
        Error::new(
            ErrorKind::NotFound,
            format!("path {} has no parent directory", path.display()),
        )
    })?;

    // empty paths never exist, so replace it with the current directory
    // the current directory should really always exist
    if parent == Path::new("") {
        return Ok(());
    }

    if !parent.exists() {
        return Err(Error::new(
            ErrorKind::NotFound,
            format!("parent directory {} does not exist", parent.display()),
        ));
    }

    Ok(())
}

impl TryFrom<PathBuf> for FilePath {
    type Error = Error;

    fn try_from(path: PathBuf) -> Result<Self> {
        if !path.is_file() {
            return Err(Error::new(
                ErrorKind::NotFound,
                format!("path {} was not a file", path.display()),
            ));
        }
        Ok(Self(path))
    }
}

impl TryFrom<PathBuf> for NewFilePath {
    type Error = Error;

    fn try_from(path: PathBuf) -> Result<Self> {
        check_new_path(&path)?;
        Ok(Self(path))
    }
}

impl TryFrom<PathBuf> for DirPath {
    type Error = Error;

    fn try_from(path: PathBuf) -> Result<Self> {
        if !path.is_dir() {
            return Err(Error::new(
                ErrorKind::NotFound,
                format!("path {} was not a directory", path.display()),
            ));
        }

        Ok(Self(path))
    }
}

impl TryFrom<PathBuf> for NewDirPath {
    type Error = Error;

    fn try_from(path: PathBuf) -> Result<Self> {
        check_new_path(&path)?;
        Ok(Self(path))
    }
}

impl NewFilePath {
    /// Creates the new file.
    ///
    /// See [`File::create`].
    pub fn create(&self) -> Result<File> {
        File::create(&self.0)
    }
}

impl NewDirPath {
    /// Creates the new directory.
    ///
    /// See [`DirBuilder::create`].
    pub fn create(&self) -> Result<()> {
        DirBuilder::new().create(&self.0)
    }
}
