//! Utilities for creating visualizations of mass spectra.

use std::path::Path;

use plotters::{
    backend::{DrawingBackend, SVGBackend},
    chart::ChartBuilder,
    coord::Shift,
    drawing::DrawingAreaErrorKind,
    prelude::{DrawingArea, IntoDrawingArea, LabelAreaPosition},
    style::{colors, AsRelative, Color, RGBAColor},
};
use rustc_hash::FxHashSet;

use crate::Peaks;

/// Settings for drawing a mass spectrum.
///
/// This holds on to a reference to any [`Peaks`] implementer and follows a builder pattern to
/// configure how that spectrum should be drawn. To consume to builder to create a drawing use the
/// [`Drawing::draw_on`] or the utility [`Drawing::draw_svg`].
///
/// # Drawings
/// Spectra are drawn as a plot with [`Peak::mz`]s on the x-axis and [`Peak::intensity`]s on the
/// y-axis. Each peak is drawn as a thin bar on these axes.
///
/// # Configurable Parameters
/// | Parameter | Description | Default |
/// | ---- | ---- | ---- |
/// | [background color] | Background of full drawing | Transparent |
/// | [peak colors] | Color that each peak is drawn | Black for all peaks |
/// | [mz labels] | Show the m/z axis tick labels | True |
/// | [intensity labels] | Show the intensity axis tick labels | False |
/// | [number of peaks] | Control number of most intense peaks drawn | All peaks in spectrum |
///
/// [background color]: Drawing::background_color
/// [peak colors]: Drawing::color_peaks
/// [mz labels]: Drawing::hide_mz_labels
/// [intensity labels]: Drawing::show_intensity_labels
/// [number of peaks]: Drawing::limit_peaks
/// [`Peak::mz`]: crate::Peak::mz
/// [`Peak::intensity`]: crate::Peak::intensity
#[derive(Debug, Clone)]
pub struct Drawing<'a, S> {
    spectrum: &'a S,
    background_color: RGBAColor,
    peak_colors: Vec<RGBAColor>,
    mz_labels: bool,
    intensity_labels: bool,
    num_peaks: usize,
}

impl<'a, S: Peaks> Drawing<'a, S> {
    /// Create a new spectrum drawing.
    ///
    /// The provided `spectrum` will be drawn when [`Drawing::draw_on`] is called.
    pub fn new(spectrum: &'a S) -> Self {
        Self {
            spectrum,
            background_color: colors::TRANSPARENT,
            peak_colors: vec![colors::BLACK.to_rgba(); spectrum.peaks().len()],
            mz_labels: true,
            intensity_labels: false,
            num_peaks: spectrum.peaks().len(),
        }
    }

    /// Reset all drawing settings to the defaults.
    pub fn reset(mut self) -> Self {
        self.background_color = colors::TRANSPARENT;
        self.peak_colors = vec![colors::BLACK.to_rgba(); self.spectrum.peaks().len()];
        self.mz_labels = true;
        self.intensity_labels = false;
        self.num_peaks = self.spectrum.peaks().len();
        self
    }

    /// Set the background color of the drawing.
    ///
    /// This background color will be applied to the full drawing area passed to
    /// [`Drawing::draw_on`]. The default background color is transparent.
    pub fn background_color<C: Color>(mut self, color: &C) -> Self {
        self.background_color = color.to_rgba();
        self
    }

    /// Set the colors of each peak in the spectrum drawing.
    ///
    /// The `colors` iterator parameter should have an element for each peak in the spectrum that
    /// this drawing was constructed with. The first color in `colors` will be the color that the
    /// first peak in the spectrum will be drawn with. The default colors are black for all peaks.
    pub fn color_peaks<'b, I: IntoIterator<Item = &'b C>, C: Color + 'b>(
        mut self,
        colors: I,
    ) -> Self {
        self.peak_colors = colors.into_iter().map(Color::to_rgba).collect();
        self
    }

    /// Get a mutable reference to the color of the peak at `peak_index`.
    ///
    /// This is a utility to make coloring a small number of peaks more convenient.
    ///
    /// # Panics
    /// If `peak_index` is out of bounds this will panic.
    pub fn peak_color_mut(&mut self, peak_index: usize) -> &mut RGBAColor {
        &mut self.peak_colors[peak_index]
    }

    /// Display the tick labels for peak m/zs on the x-axis.
    ///
    /// This is already the default functionality.
    pub fn show_mz_labels(mut self) -> Self {
        self.mz_labels = true;
        self
    }

    /// Do not display the tick labels for peak m/zs on the x-axis.
    ///
    /// This is not the default functionality. By default the m/z tick labels will be shown and
    /// calling this method disables them.
    pub fn hide_mz_labels(mut self) -> Self {
        self.mz_labels = false;
        self
    }

    /// Display the tick labels for peak intensities on the y-axis.
    ///
    /// This is not the default functionality. By default intensity tick labels will not be shown
    /// and call this method enables them.
    pub fn show_intensity_labels(mut self) -> Self {
        self.intensity_labels = true;
        self
    }

    /// Do not display the tick labels for peak intensities on the y-axis.
    ///
    /// This is the default functionality.
    pub fn hide_intensity_labels(mut self) -> Self {
        self.intensity_labels = false;
        self
    }

    /// Limit drawn peaks to the `count` most intense peaks.
    ///
    /// All peaks outside of the `count` most intense will not be drawn. If you set peak colors
    /// using [`Drawing::color_peaks`] those colors will still be respected, however this
    /// will cause only the `count` most intense peaks to be drawn with their colors. The default
    /// functionality is to draw all peaks, regardless of intensity.
    pub fn limit_peaks(mut self, count: usize) -> Self {
        self.num_peaks = count;
        self
    }

    /// Using the provided settings draw this spectrum on the provided `area`.
    ///
    /// # Errors
    /// This propagates errors from the [`plotters`] API.
    pub fn draw_on<DB: DrawingBackend>(
        self,
        area: &DrawingArea<DB, Shift>,
    ) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
        // fill in the background
        area.fill(&self.background_color)?;

        // determine the bounds for m/z and intensity axes
        let mut min_mz = f64::INFINITY;
        let mut max_mz = 0.0;
        let mut max_int = 0.0;

        for peak in self.spectrum.peaks() {
            if peak.mz < min_mz {
                min_mz = peak.mz;
            }
            if peak.mz > max_mz {
                max_mz = peak.mz;
            }
            if peak.intensity > max_int {
                max_int = peak.intensity;
            }
        }

        // pick the indices we use according to `self.num_peaks`
        let mut desc_intensity = (0..self.spectrum.peaks().len()).collect::<Vec<_>>();
        desc_intensity.sort_by(|&x, &y| {
            self.spectrum.peaks()[x]
                .intensity
                .partial_cmp(&self.spectrum.peaks()[y].intensity)
                .unwrap()
                .reverse()
        });
        let target_peak_indices = desc_intensity
            .into_iter()
            .take(self.num_peaks)
            .collect::<FxHashSet<_>>();

        let mut chart_builder = &mut ChartBuilder::on(area);

        // reserve space for labels
        // if they want ticks we need to save more space on that axis' label

        let bottom_percent = match self.mz_labels {
            true => 15.percent(),
            false => 10.percent(),
        };
        let left_percent = match self.intensity_labels {
            true => 15.percent(),
            false => 10.percent(),
        };

        chart_builder = chart_builder
            .margin(5u32)
            .set_label_area_size(LabelAreaPosition::Left, left_percent)
            .set_label_area_size(LabelAreaPosition::Bottom, bottom_percent);

        // keep label sizes consistent, this uses a percentage of the smaller dimension
        let desc_style = ("serif", 8.percent(), &colors::BLACK);

        // build and configure the chart and mesh
        let mut chart = chart_builder.build_cartesian_2d(min_mz..max_mz, 0.0..max_int)?;
        let mut mesh = &mut chart.configure_mesh();
        mesh.disable_mesh()
            .axis_desc_style(desc_style)
            .x_desc("m/z")
            .y_desc("Intensity");

        // if users don't turn on tick labels don't draw them
        if !self.mz_labels {
            mesh = mesh.x_labels(0);
        }
        if !self.intensity_labels {
            mesh = mesh.y_labels(0);
        }
        mesh.draw()?;

        // draw the peaks
        chart.draw_series(
            self.spectrum
                .peaks()
                .iter()
                .zip(&self.peak_colors)
                .enumerate()
                .filter_map(|(idx, (peak, color))| {
                    if !target_peak_indices.contains(&idx) {
                        return None;
                    }
                    Some(plotters::prelude::PathElement::new(
                        [(peak.mz, 0.0), (peak.mz, peak.intensity)],
                        color,
                    ))
                }),
        )?;

        Ok(())
    }

    /// Using the provided settings draw this spectrum into an SVG at `path` with `size`.
    ///
    /// This is a utility that simple wraps the construction of a [`plotters::backend::SVGBackend`]
    /// and a call to [`Drawing::draw_on`].
    pub fn draw_svg<P: AsRef<Path> + ?Sized>(
        self,
        path: &P,
        size: (u32, u32),
    ) -> Result<(), DrawingAreaErrorKind<<SVGBackend<'_> as DrawingBackend>::ErrorType>> {
        let backend = SVGBackend::new(path, size);
        let area = backend.into_drawing_area();
        self.draw_on(&area)
    }
}

#[cfg(test)]
mod tests {
    use crate::{ParsedCollection, ParsedSpectrum, Spectrum};

    use super::*;

    fn spec() -> Spectrum {
        let spectrum = Spectrum::try_from(&ParsedSpectrum {
            collection: ParsedCollection {
                path: Path::new("./test_files/gnps.mgf").try_into().unwrap(),
                pre: Default::default(),
            },
            scan: 1,
        })
        .unwrap();
        spectrum
    }

    #[test]
    fn reset() {
        let spectrum = spec();
        let drawing = Drawing::new(&spectrum)
            .hide_mz_labels()
            .show_intensity_labels()
            .color_peaks([&colors::RED; 56])
            .limit_peaks(30)
            .background_color(&colors::BLUE)
            .reset();

        assert!(drawing.mz_labels);
        assert!(!drawing.intensity_labels);
        assert_eq!(drawing.background_color, colors::TRANSPARENT);
        assert_eq!(drawing.num_peaks, 56);
        assert_eq!(drawing.peak_colors, [colors::BLACK.to_rgba(); 56]);
    }

    #[test]
    fn background_color() {
        let spectrum = spec();
        let mut drawing = Drawing::new(&spectrum);
        assert_eq!(drawing.background_color, colors::TRANSPARENT);
        drawing = drawing.background_color(&colors::RED);
        assert_eq!(drawing.background_color, colors::RED.to_rgba());
        drawing = drawing.background_color(&colors::GREEN);
        assert_eq!(drawing.background_color, colors::GREEN.to_rgba());
    }

    #[test]
    fn peak_colors() {
        let spectrum = spec();
        let mut drawing = Drawing::new(&spectrum);
        assert_eq!(drawing.peak_colors, vec![colors::BLACK.to_rgba(); 56]);
        drawing = drawing.color_peaks([&colors::BLUE; 56]);
        assert_eq!(drawing.peak_colors, vec![colors::BLUE.to_rgba(); 56]);
        *drawing.peak_color_mut(10) = colors::GREEN.to_rgba();
        assert_eq!(drawing.peak_colors[10], colors::GREEN.to_rgba());
        assert!(drawing
            .peak_colors
            .iter()
            .take(10)
            .all(|&x| x == colors::BLUE.to_rgba()));
        assert!(drawing
            .peak_colors
            .iter()
            .skip(11)
            .all(|&x| x == colors::BLUE.to_rgba()));
    }

    #[test]
    fn mz_labels() {
        let spectrum = spec();
        let mut drawing = Drawing::new(&spectrum);
        assert!(drawing.mz_labels);
        drawing = drawing.hide_mz_labels();
        assert!(!drawing.mz_labels);
        drawing = drawing.show_mz_labels();
        assert!(drawing.mz_labels);
    }

    #[test]
    fn intensity_labels() {
        let spectrum = spec();
        let mut drawing = Drawing::new(&spectrum);
        assert!(!drawing.intensity_labels);
        drawing = drawing.hide_intensity_labels();
        assert!(!drawing.intensity_labels);
        drawing = drawing.show_intensity_labels();
        assert!(drawing.intensity_labels);
    }

    #[test]
    fn num_peaks() {
        let spectrum = spec();
        let mut drawing = Drawing::new(&spectrum);
        assert_eq!(drawing.num_peaks, 56);
        drawing = drawing.limit_peaks(10);
        assert_eq!(drawing.num_peaks, 10);
        drawing = drawing.limit_peaks(40);
        assert_eq!(drawing.num_peaks, 40);
    }

    #[test]
    fn draw() {
        let p = tempfile::NamedTempFile::new().unwrap();
        let area = SVGBackend::new(p.path(), (400, 500)).into_drawing_area();
        let spectrum = Spectrum::try_from(&ParsedSpectrum {
            collection: ParsedCollection {
                path: Path::new("./test_files/gnps.mgf").try_into().unwrap(),
                pre: Default::default(),
            },
            scan: 1,
        })
        .unwrap();
        assert!(Drawing::new(&spectrum)
            .limit_peaks(10)
            .show_mz_labels()
            .hide_intensity_labels()
            .draw_on(&area)
            .is_ok());
    }
}
