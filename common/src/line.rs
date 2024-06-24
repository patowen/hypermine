use std::ops::Neg;

use crate::{
    dodeca::{Side, Vertex},
    math::{self, lorentz_normalize, mip, origin, translate, translate_along},
    Plane,
};

/// A hyperbolic line
#[derive(Debug, Copy, Clone)]
pub struct Line<N: na::RealField> {
    vector: na::Vector4<N>,
    point: na::Vector4<N>,
}

impl From<Side> for Line<f64> {
    /// A surface overlapping with a particular dodecahedron side
    fn from(side: Side) -> Self {
        Self {
            vector: *side.normal_f64(),
            point: *&origin(),
        }
    }
}

impl<N: na::RealField + Copy> From<na::Unit<na::Vector3<N>>> for Line<N> {
    /// A line passing through the origin
    fn from(x: na::Unit<na::Vector3<N>>) -> Self {
        Self {
            vector: x.into_inner().push(na::zero()),
            point: *&origin(),
        }
    }
}

impl<N: na::RealField + Copy> Neg for Line<N> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            vector: -self.vector,
            point: self.point,
        }
    }
}

impl<N: na::RealField + Copy> Line<N> {
    /// Hyperbolic vector identifying the line
    pub fn vector(&self) -> &na::Vector4<N> {
        &self.vector
    }

    /// Hyperbolic point identifying the line
    pub fn point(&self) -> &na::Vector4<N> {
        &self.point
    }

    /// Shortest distance between the line and a point
    pub fn distance_to(&self, point: &na::Vector4<N>) -> N {
        (math::mip(point, &self.point).powi(2) - math::mip(point, &self.vector).powi(2))
            .sqrt()
            .acosh()
    }
    /// Like `distance_to`, but using chunk coordinates for a chunk in the same node space
    pub fn from_points(&point1: &na::Vector4<N>, &point2: &na::Vector4<N>) -> Self {
        let point1 = math::lorentz_normalize(&point1);
        let point_diff = point2 - point1;
        let a = math::lorentz_normalize(&(point_diff + point1 * math::mip(&point1, &point_diff)));
        Self {
            vector: a,
            point: point1,
        }
    }
}

impl Line<f64> {
    /// Like `distance_to`, but using chunk coordinates for a chunk in the same node space
    pub fn distance_to_chunk(&self, chunk: Vertex, coord: &na::Vector3<f64>) -> f64 {
        let pos = lorentz_normalize(&(chunk.chunk_to_node_f64() * coord.push(1.0)));
        self.distance_to(&pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{origin, translate_along};
    use approx::*;

    #[test]
    fn distance_sanity() {
        for &axis in &[
            na::Vector3::x_axis(),
            na::Vector3::y_axis(),
            na::Vector3::z_axis(),
        ] {
            for &axis2 in &[
                na::Vector3::x_axis(),
                na::Vector3::y_axis(),
                na::Vector3::z_axis(),
            ] {
                for &distance in &[-1.5f32, 0.0f32, 1.5f32] {
                    if axis != axis2 {
                        let line = Line::from(axis);
                        assert_abs_diff_eq!(
                            line.distance_to(&(translate_along(&(*axis2 * distance)) * origin())),
                            distance.abs()
                        );
                    }
                }
            }
        }
    }
}
