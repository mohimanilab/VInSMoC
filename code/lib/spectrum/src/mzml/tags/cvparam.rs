use serde::Deserialize;
use std::{fmt, str::FromStr};

pub trait CvParams {
    /// Get a reference to the controlled variable parameters for this type
    fn params(&self) -> &Vec<CvParam>;

    /// Retrieve parameter with given accession.
    ///
    /// Returns [`None`] if no such [`CvParam`] exists.
    fn param_with_acc<S: AsRef<str>>(&self, accession: S) -> Option<&CvParam> {
        self.params()
            .iter()
            .find(|param| param.accession == accession.as_ref())
    }

    /// Retrieve parameter with given name.
    ///
    /// Returns [`None`] if no such [`CvParam`] exists.
    fn param_with_name<S: AsRef<str>>(&self, name: S) -> Option<&CvParam> {
        self.params()
            .iter()
            .find(|param| param.name == name.as_ref())
    }
}

#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct CvParam {
    pub accession: String,
    pub cv_ref: String,
    pub name: String,
    pub unit_accession: Option<String>,
    pub unit_cv_ref: Option<String>,
    pub unit_name: Option<String>,
    pub value: Option<String>,
}

impl CvParam {
    /// Retrieve the value of this parameter as the given type.
    ///
    /// If there was no `value` attribute then this will pass and return [`None`].
    /// If the `value` attribute existing but couldn't be parsed as the given type, then this will
    /// be an [`Err`].
    pub fn value<T>(&self) -> Result<Option<T>, T::Err>
    where
        T: FromStr + fmt::Debug,
        T::Err: std::error::Error,
    {
        if self.value.is_none() {
            return Ok(None);
        }

        let value_str = self.value.as_ref().unwrap().as_str();
        let value = value_str.parse::<T>()?;

        Ok(Some(value))
    }

    /// Get the name of units of this parameter.
    ///
    /// If no units were specified then this returns [`None`].
    pub fn units(&self) -> Option<&str> {
        self.unit_name.as_deref()
    }
}
