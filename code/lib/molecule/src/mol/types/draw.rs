use std::ops::{Add, Div, Sub};

use petgraph::visit::{EdgeRef, IntoEdgeReferences};

use super::Mol;

#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct Point {
    pub x: f32,
    pub y: f32,
}

impl From<(f32, f32)> for Point {
    fn from((x, y): (f32, f32)) -> Self {
        Self { x, y }
    }
}

impl From<&Point> for (i32, i32) {
    fn from(point: &Point) -> Self {
        (point.x as i32, point.y as i32)
    }
}

impl Div<f32> for Point {
    type Output = Point;

    fn div(self, rhs: f32) -> Self::Output {
        Point {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Add<&Point> for &Point {
    type Output = Point;

    fn add(self, rhs: &Point) -> Self::Output {
        Point {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub<&Point> for &Point {
    type Output = Point;

    fn sub(self, rhs: &Point) -> Self::Output {
        Point {
            x: self.x - rhs.y,
            y: self.y - rhs.y,
        }
    }
}

impl Point {
    pub fn centroided(points: &[Self]) -> Self {
        debug_assert!(!points.is_empty());
        let min_x = points
            .iter()
            .map(|point| point.x)
            .min_by(|x1, x2| x1.partial_cmp(x2).unwrap())
            .unwrap();
        let min_y = points
            .iter()
            .map(|point| point.y)
            .min_by(|y1, y2| y1.partial_cmp(y2).unwrap())
            .unwrap();

        let max_x = points
            .iter()
            .map(|point| point.x)
            .max_by(|x1, x2| x1.partial_cmp(x2).unwrap())
            .unwrap();
        let max_y = points
            .iter()
            .map(|point| point.y)
            .max_by(|y1, y2| y1.partial_cmp(y2).unwrap())
            .unwrap();

        Self {
            x: (min_x + max_x) / 2.0,
            y: (min_y + max_y) / 2.0,
        }
    }

    /// Translates point so the center is located at the top left corner of the canvas, rather than
    /// the center of the canvas.
    ///
    /// This is useful for when your drawing system has an origin in the top left corner (e.g.
    /// HTML5 canvas). The shifting is done by setting the point at `(min_x, max_y)` to `(0, 0)`
    /// and shifting all coordinates based on that. The `top_left` parameter is this reference
    /// point.
    ///
    /// The origin (`o`) and coordinate system moves as follows.
    /// ```txt
    ///      Original              Translated
    ///    ▲ #############       ▲ #############
    ///    | #           #       0 # o         #
    /// +y | #           #       | #           #
    ///    | #           #       | #           #
    ///    0 #     o     #    +y | #           #
    ///    | #           #       | #           #
    /// -y | #           #       | #           #
    ///    | #           #       | #           #
    ///    ▼ #############       ▼ #############
    ///      ◀-----0-----▶         ◀0----------▶
    ///       -x     +x                 +x
    /// ```
    pub fn translate(&mut self, top_left: &Point) {
        self.y = -self.y + top_left.y;
        self.x -= top_left.x;
        debug_assert!(self.x >= 0.0 && self.y >= 0.0);
    }

    /// Normalize x and y coordinates by provided `max` point.
    ///
    /// Applying this to a collection of points converts them to percentages that can easily be
    /// scaled to any canvas size. Note there may be negative values if [`Point::translate`] is not called
    /// first.
    pub fn normalize(&mut self, max: &Point) {
        self.x /= max.x;
        self.y /= max.y;
    }

    pub fn scale(&mut self, (width, height): (f32, f32), (padding_x, padding_y): (f32, f32)) {
        self.x = (width - padding_x * width) * self.x + padding_x / 2.0 * width;
        self.y = (height - padding_x * height) * self.y + padding_y / 2.0 * height;
    }
}

impl Mol {
    pub fn coordinates(&self) -> Vec<Point> {
        let atoms = self
            .graph
            .node_weights()
            .map(|x| x.atomic_num() as u8)
            .collect::<Vec<u8>>();
        let bonds = self
            .graph
            .edge_references()
            .map(|edge| {
                [
                    edge.source().index() as u16,
                    edge.target().index() as u16,
                    *edge.weight() as u16,
                ]
            })
            .collect::<Vec<[u16; 3]>>();
        let mut coords = unsafe { coordgen::gen_coords_unchecked(&atoms, &bonds) }
            .into_iter()
            .map(Point::from)
            .collect::<Vec<_>>();

        // translate first
        let x = coords
            .iter()
            .map(|point| point.x)
            .min_by(|p1, p2| p1.partial_cmp(p2).unwrap())
            .unwrap();
        let y = coords
            .iter()
            .map(|point| point.y)
            .max_by(|p1, p2| p1.partial_cmp(p2).unwrap())
            .unwrap();
        let top_left = Point { x, y };
        for point in coords.iter_mut() {
            point.translate(&top_left);
        }

        // then normalize
        let x = coords
            .iter()
            .map(|point| point.x)
            .max_by(|p1, p2| p1.partial_cmp(p2).unwrap())
            .unwrap();
        let y = coords
            .iter()
            .map(|point| point.y)
            .max_by(|p1, p2| p1.partial_cmp(p2).unwrap())
            .unwrap();
        let max = Point { x, y };
        for point in coords.iter_mut() {
            point.normalize(&max);
        }

        coords
    }
}
