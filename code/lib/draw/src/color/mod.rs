//! Colors that can be parsed from strings.

use std::str::FromStr;

use nom::error::convert_error;
pub use plotters::style::Color as PlottersColor;
use plotters::style::{HSLColor, RGBAColor};
use plotters_backend::BackendColor;

mod parser;

// TODO: HSL conversion in plotters seems to have bugs
// examine further

/// A color that can be parsed from a string.
///
/// The color is stored in RGBA format for convenience.
///
/// String representations of colors map to RGB, RGBA, or HSL color representations. The formats
/// for each of these are a case-insensitive format specifier followed by the numeric components
/// for that format.
///
/// # Examples
/// ```
/// use draw::Color;
///
/// assert_eq!(
///     "RGB(10,10,10)".parse::<Color>().unwrap(),
///     Color { r: 10, g: 10, b: 10, alpha: 1.0 }
/// );
/// assert_eq!(
///     "RGBA(10,10,10,2.0)".parse::<Color>().unwrap(),
///     Color { r: 10, g: 10, b: 10, alpha: 2.0 }
/// );
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color {
    /// The red channel.
    pub r: u8,
    /// The green channel.
    pub g: u8,
    /// The blue channel.
    pub b: u8,
    /// The alpha channel.
    pub alpha: f64,
}

impl From<(u8, u8, u8)> for Color {
    fn from((r, g, b): (u8, u8, u8)) -> Self {
        Self {
            r,
            g,
            b,
            alpha: 1.0,
        }
    }
}

impl From<(u8, u8, u8, f64)> for Color {
    fn from((r, g, b, alpha): (u8, u8, u8, f64)) -> Self {
        Self { r, g, b, alpha }
    }
}

impl From<(f64, f64, f64)> for Color {
    fn from((h, s, l): (f64, f64, f64)) -> Self {
        use plotters::style::Color;
        dbg!(h, s, l, HSLColor(h, s, l).to_rgba());
        let RGBAColor(r, g, b, alpha) = HSLColor(h, s, l).to_rgba();
        Self { r, g, b, alpha }
    }
}

impl plotters::style::Color for Color {
    fn to_backend_color(&self) -> BackendColor {
        BackendColor {
            alpha: self.alpha,
            rgb: (self.r, self.g, self.b),
        }
    }
}

/// An error that occurred when parsing a color.
#[derive(Debug, thiserror::Error)]
#[error("failed to parse color specifier:\n{0}")]
pub struct ParseColorError(String);

impl FromStr for Color {
    type Err = ParseColorError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match parser::parse(s) {
            Err(nom::Err::Error(e) | nom::Err::Failure(e)) => {
                Err(ParseColorError(convert_error(s, e)))
            }
            Err(nom::Err::Incomplete(_)) => unreachable!(),
            Ok((_, c)) => Ok(c),
        }
    }
}
