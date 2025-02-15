//! Tools for processing the geometry of a right dodecahedron

use data::*;
use serde::{Deserialize, Serialize};

use crate::math::{MIsometry, MVector};
use crate::voxel_math::ChunkAxisPermutation;

/// Sides of a right dodecahedron
///
/// These sides are arranged based on the following adjacency graph, although it
/// is recommended not to hardcode side names in other code:
/// ```nocode
///          A
/// (D-) I-E-B-C-D (-I)
///  (K-) L-G-F-H-K (-L)
///           J
/// ```
/// The above adjacency graph can be read as a world map, where side A is at the
/// north pole, and side J is at the south pole.
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub enum Side {
    A,
    B,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    J,
    K,
    L,
}

impl Side {
    pub const COUNT: usize = 12;

    pub const VALUES: [Self; Self::COUNT] = [
        Self::A,
        Self::B,
        Self::C,
        Self::D,
        Self::E,
        Self::F,
        Self::G,
        Self::H,
        Self::I,
        Self::J,
        Self::K,
        Self::L,
    ];

    #[inline]
    pub fn from_index(x: usize) -> Self {
        use Side::*;
        const VALUES: [Side; Side::COUNT] = [A, B, C, D, E, F, G, H, I, J, K, L];
        VALUES[x]
    }

    pub fn iter() -> impl ExactSizeIterator<Item = Self> {
        Self::VALUES.iter().copied()
    }

    /// Whether `self` and `other` share an edge
    ///
    /// `false` when `self == other`.
    #[inline]
    pub fn adjacent_to(self, other: Side) -> bool {
        ADJACENT[self as usize][other as usize]
    }

    /// Outward normal vector of this side
    #[inline]
    pub fn normal(self) -> &'static MVector<f32> {
        &SIDE_NORMALS_F32[self as usize]
    }

    /// Outward normal vector of this side
    #[inline]
    pub fn normal_f64(self) -> &'static MVector<f64> {
        &SIDE_NORMALS_F64[self as usize]
    }

    /// Reflection across this side
    #[inline]
    pub fn reflection(self) -> &'static MIsometry<f32> {
        &REFLECTIONS_F32[self as usize]
    }

    /// Reflection across this side
    #[inline]
    pub fn reflection_f64(self) -> &'static MIsometry<f64> {
        &REFLECTIONS_F64[self as usize]
    }

    /// Whether `p` is opposite the dodecahedron across the plane containing `self`
    #[inline]
    pub fn is_facing(self, p: &MVector<f32>) -> bool {
        let r = self.reflection().row(3).clone_owned();
        (r * na::Vector4::from(*p)).x < p.w
    }
}

/// Vertices of a right dodecahedron
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Vertex {
    A,
    B,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    J,
    K,
    L,
    M,
    N,
    O,
    P,
    Q,
    R,
    S,
    T,
}

impl Vertex {
    pub const COUNT: usize = 20;

    pub const VALUES: [Self; Self::COUNT] = [
        Self::A,
        Self::B,
        Self::C,
        Self::D,
        Self::E,
        Self::F,
        Self::G,
        Self::H,
        Self::I,
        Self::J,
        Self::K,
        Self::L,
        Self::M,
        Self::N,
        Self::O,
        Self::P,
        Self::Q,
        Self::R,
        Self::S,
        Self::T,
    ];

    pub fn iter() -> impl ExactSizeIterator<Item = Self> {
        Self::VALUES.iter().copied()
    }

    /// Vertex shared by three sides, if any
    #[inline]
    pub fn from_sides(a: Side, b: Side, c: Side) -> Option<Self> {
        SIDES_TO_VERTEX[a as usize][b as usize][c as usize]
    }

    /// Sides incident to this vertex, in canonical order.
    ///
    /// This canonical order determines the X, Y, and Z axes of the chunk
    /// corresponding to the vertex.
    #[inline]
    pub fn canonical_sides(self) -> [Side; 3] {
        VERTEX_SIDES[self as usize]
    }

    /// Vertices adjacent to this vertex in canonical order.
    ///
    /// The canonical order of adjacent vertices is based on the canonical order
    /// of sides incident to the vertex, as each of the three adjacent vertices
    /// corresponds to one of the three sides. As for which side, when two
    /// vertices are adjacent, they share two out of three sides of the
    /// dodecahedron. The side they do _not_ share is the side they correspond
    /// to.
    ///
    /// Put another way, anything leaving a chunk in the negative-X direction
    /// will end up crossing `canonical_sides()[0]`, while anything leaving a
    /// chunk in the positive-X direction will end up arriving at
    /// `adjacent_vertices()[0]`.
    #[inline]
    pub fn adjacent_vertices(self) -> [Vertex; 3] {
        ADJACENT_VERTICES[self as usize]
    }

    /// Chunk axes permutations for vertices adjacent to this vertex in
    /// canonical order.
    ///
    /// The chunks of two adjacent vertices meet at a plane. When swiching
    /// reference frames from one vertex to another, it is necessary to reflect
    /// about this plane and then apply the permutation returned by this
    /// function.
    #[inline]
    pub fn chunk_axis_permutations(self) -> &'static [ChunkAxisPermutation; 3] {
        &CHUNK_AXIS_PERMUTATIONS[self as usize]
    }

    /// For each vertex of the cube dual to this dodecahedral vertex, provides an iterator of at
    /// most 3 steps to reach the corresponding graph node, and binary coordinates of the vertex in
    /// question with respect to the origin vertex of the cube.
    pub fn dual_vertices(
        self,
    ) -> impl ExactSizeIterator<Item = ([bool; 3], impl ExactSizeIterator<Item = Side>)> {
        let [a, b, c] = self.canonical_sides();
        let verts = [
            ([Side::A; 3], 0, [false, false, false]),
            ([c, Side::A, Side::A], 1, [false, false, true]),
            ([b, Side::A, Side::A], 1, [false, true, false]),
            ([b, c, Side::A], 2, [false, true, true]),
            ([a, Side::A, Side::A], 1, [true, false, false]),
            ([a, c, Side::A], 2, [true, false, true]),
            ([a, b, Side::A], 2, [true, true, false]),
            ([a, b, c], 3, [true, true, true]),
        ];
        (0..8).map(move |i| {
            let (sides, len, coords) = verts[i];
            (coords, (0..len).map(move |i| sides[i]))
        })
    }

    /// Transform from euclidean chunk coordinates to hyperbolic node space
    pub fn chunk_to_node(self) -> na::Matrix4<f32> {
        na::Matrix4::from(*self.dual_to_node())
            * na::Matrix4::new_scaling(1.0 / Self::dual_to_chunk_factor())
    }

    /// Transform from euclidean chunk coordinates to hyperbolic node space
    pub fn chunk_to_node_f64(self) -> na::Matrix4<f64> {
        na::Matrix4::from(*self.dual_to_node_f64())
            * na::Matrix4::new_scaling(1.0 / Self::dual_to_chunk_factor_f64())
    }

    /// Transform from hyperbolic node space to euclidean chunk coordinates
    pub fn node_to_chunk(self) -> na::Matrix4<f32> {
        na::Matrix4::new_scaling(Self::dual_to_chunk_factor())
            * na::Matrix4::from(*self.node_to_dual())
    }

    /// Transform from hyperbolic node space to euclidean chunk coordinates
    pub fn node_to_chunk_f64(self) -> na::Matrix4<f64> {
        na::Matrix4::new_scaling(Self::dual_to_chunk_factor_f64())
            * na::Matrix4::from(*self.node_to_dual_f64())
    }

    /// Transform from cube-centric coordinates to dodeca-centric coordinates
    pub fn dual_to_node(self) -> &'static MIsometry<f32> {
        &DUAL_TO_NODE_F32[self as usize]
    }

    /// Transform from cube-centric coordinates to dodeca-centric coordinates
    pub fn dual_to_node_f64(self) -> &'static MIsometry<f64> {
        &DUAL_TO_NODE_F64[self as usize]
    }

    /// Transform from dodeca-centric coordinates to cube-centric coordinates
    pub fn node_to_dual(self) -> &'static MIsometry<f32> {
        &NODE_TO_DUAL_F32[self as usize]
    }

    /// Transform from dodeca-centric coordinates to cube-centric coordinates
    pub fn node_to_dual_f64(self) -> &'static MIsometry<f64> {
        &NODE_TO_DUAL_F64[self as usize]
    }

    /// Scale factor used in conversion from cube-centric coordinates to euclidean chunk coordinates.
    /// Scaling the x, y, and z components of a vector in cube-centric coordinates by this value
    /// and dividing them by the w coordinate will yield euclidean chunk coordinates.
    pub fn dual_to_chunk_factor() -> f32 {
        *DUAL_TO_CHUNK_FACTOR_F32
    }

    /// Scale factor used in conversion from cube-centric coordinates to euclidean chunk coordinates.
    /// Scaling the x, y, and z components of a vector in cube-centric coordinates by this value
    /// and dividing them by the w coordinate will yield euclidean chunk coordinates.
    pub fn dual_to_chunk_factor_f64() -> f64 {
        *DUAL_TO_CHUNK_FACTOR_F64
    }

    /// Scale factor used in conversion from euclidean chunk coordinates to cube-centric coordinates.
    /// Scaling the x, y, and z components of a vector in homogeneous euclidean chunk coordinates by this value
    /// and lorentz-normalizing the result will yield cube-centric coordinates.
    pub fn chunk_to_dual_factor() -> f32 {
        *CHUNK_TO_DUAL_FACTOR_F32
    }

    /// Scale factor used in conversion from euclidean chunk coordinates to cube-centric coordinates.
    /// Scaling the x, y, and z components of a vector in homogeneous euclidean chunk coordinates by this value
    /// and lorentz-normalizing the result will yield cube-centric coordinates.
    pub fn chunk_to_dual_factor_f64() -> f64 {
        *CHUNK_TO_DUAL_FACTOR_F64
    }

    /// Convenience method for `self.chunk_to_node().determinant() < 0`.
    pub fn parity(self) -> bool {
        CHUNK_TO_NODE_PARITY[self as usize]
    }
}

pub const BOUNDING_SPHERE_RADIUS_F64: f64 = 1.2264568712514068;
pub const BOUNDING_SPHERE_RADIUS: f32 = BOUNDING_SPHERE_RADIUS_F64 as f32;

mod data {
    use std::array;
    use std::sync::LazyLock;

    use crate::dodeca::{Side, Vertex};
    use crate::math::{self, MIsometry, MVector};
    use crate::voxel_math::ChunkAxisPermutation;

    /// Whether two sides share an edge
    pub static ADJACENT: LazyLock<[[bool; Side::COUNT]; Side::COUNT]> = LazyLock::new(|| {
        Side::VALUES.map(|side0| {
            Side::VALUES.map(|side1| {
                let cosh_distance = (*side0.reflection_f64() * *side1.reflection_f64())[(3, 3)];
                // Possile cosh_distances: 1, 4.23606 = 2+sqrt(5), 9.47213 = 5+2*sqrt(5), 12.70820 = 6+3*sqrt(5);
                // < 2.0 indicates identical faces; < 5.0 indicates adjacent faces; > 5.0 indicates non-adjacent faces
                (2.0..5.0).contains(&cosh_distance)
            })
        })
    });

    /// Vector corresponding to the outer normal of each side
    pub static SIDE_NORMALS_F64: LazyLock<[MVector<f64>; Side::COUNT]> = LazyLock::new(|| {
        let phi = libm::sqrt(1.25) + 0.5; // golden ratio
        let template_normal = MVector::new(1.0, phi, 0.0, libm::sqrt(phi)).normalized();
        let signed_template_normals = {
            let n = template_normal;
            [
                MVector::new(n.x, n.y, n.z, n.w),
                MVector::new(-n.x, n.y, -n.z, n.w),
                MVector::new(n.x, -n.y, -n.z, n.w),
                MVector::new(-n.x, -n.y, n.z, n.w),
            ]
        };

        Side::VALUES.map(|side| {
            let signed_template_normal = signed_template_normals[side as usize / 3];
            math::tuv_to_xyz((3 - side as usize % 3) % 3, signed_template_normal)
        })
    });

    /// Transform that moves from a neighbor to a reference node, for each side
    pub static REFLECTIONS_F64: LazyLock<[MIsometry<f64>; Side::COUNT]> =
        LazyLock::new(|| SIDE_NORMALS_F64.map(|r| MIsometry::reflection(&r)));

    /// Sides incident to a vertex, in canonical order
    pub static VERTEX_SIDES: LazyLock<[[Side; 3]; Vertex::COUNT]> = LazyLock::new(|| {
        let mut result: Vec<[Side; 3]> = Vec::new();
        // Kind of a hack, but working this out by hand isn't any fun.
        for a in 0..Side::COUNT {
            for b in (a + 1)..Side::COUNT {
                for c in (b + 1)..Side::COUNT {
                    if !ADJACENT[a][b] || !ADJACENT[b][c] || !ADJACENT[c][a] {
                        continue;
                    }
                    result.push([
                        Side::from_index(a),
                        Side::from_index(b),
                        Side::from_index(c),
                    ]);
                }
            }
        }
        result
            .try_into()
            .expect("All vertices should be initialized.")
    });

    // Which vertices are adjacent to other vertices and opposite the canonical sides
    pub static ADJACENT_VERTICES: LazyLock<[[Vertex; 3]; Vertex::COUNT]> = LazyLock::new(|| {
        Vertex::VALUES.map(|vertex| {
            array::from_fn(|result_index| {
                // Reorder the canonical sides so that the third element is the
                // one that must change to find the correct adjacent vertex.
                let adjacent_sides = math::tuv_to_xyz(2 - result_index, vertex.canonical_sides());
                // Try every possible side to find an adjacent vertex.
                for test_side in Side::iter() {
                    if test_side == adjacent_sides[2] {
                        continue;
                    }
                    if let Some(adjacent_vertex) =
                        Vertex::from_sides(adjacent_sides[0], adjacent_sides[1], test_side)
                    {
                        return adjacent_vertex;
                    }
                }
                panic!("No suitable vertex found");
            })
        })
    });

    // Which transformations have to be done after a reflection to switch reference frames from one vertex
    // to one of its adjacent vertices (ordered similarly to ADJACENT_VERTICES)
    pub static CHUNK_AXIS_PERMUTATIONS: LazyLock<[[ChunkAxisPermutation; 3]; Vertex::COUNT]> =
        LazyLock::new(|| {
            array::from_fn(|vertex| {
                array::from_fn(|result_index| {
                    let mut test_sides = VERTEX_SIDES[vertex];
                    // Keep modifying the result_index'th element of test_sides until its three elements are all
                    // adjacent to a single vertex (determined using `Vertex::from_sides`).
                    for side in Side::iter() {
                        if side == VERTEX_SIDES[vertex][result_index] {
                            continue;
                        }
                        test_sides[result_index] = side;
                        let Some(adjacent_vertex) =
                            Vertex::from_sides(test_sides[0], test_sides[1], test_sides[2])
                        else {
                            continue;
                        };
                        // Compare the natural permutation of sides after a reflection from `vertex` to `adjacent_vertex`
                        // to the canonical permutation of the sides for `adjacent_vertex`.
                        return ChunkAxisPermutation::from_permutation(
                            test_sides,
                            adjacent_vertex.canonical_sides(),
                        );
                    }
                    panic!("No suitable vertex found");
                })
            })
        });

    /// Transform that converts from cube-centric coordinates to dodeca-centric coordinates
    pub static DUAL_TO_NODE_F64: LazyLock<[MIsometry<f64>; Vertex::COUNT]> = LazyLock::new(|| {
        let mip_origin_normal = MVector::origin().mip(&SIDE_NORMALS_F64[0]); // This value is the same for every side
        let mut result = [MIsometry::identity(); Vertex::COUNT];
        for (i, map) in result.iter_mut().enumerate() {
            let [a, b, c] = VERTEX_SIDES[i];
            let vertex_position = (MVector::origin()
                - (*a.normal_f64() + *b.normal_f64() + *c.normal_f64()) * mip_origin_normal)
                .normalized();
            *map = MIsometry::from_columns_unchecked(&[
                -*a.normal_f64(),
                -*b.normal_f64(),
                -*c.normal_f64(),
                vertex_position,
            ]);
        }
        result
    });

    /// Transform that converts from dodeca-centric coordinates to cube-centric coordinates
    pub static NODE_TO_DUAL_F64: LazyLock<[MIsometry<f64>; Vertex::COUNT]> =
        LazyLock::new(|| DUAL_TO_NODE_F64.map(|m| m.inverse()));

    pub static DUAL_TO_CHUNK_FACTOR_F64: LazyLock<f64> =
        LazyLock::new(|| (2.0 + 5.0f64.sqrt()).sqrt());

    pub static CHUNK_TO_DUAL_FACTOR_F64: LazyLock<f64> =
        LazyLock::new(|| 1.0 / *DUAL_TO_CHUNK_FACTOR_F64);

    /// Vertex shared by 3 sides
    pub static SIDES_TO_VERTEX: LazyLock<
        [[[Option<Vertex>; Side::COUNT]; Side::COUNT]; Side::COUNT],
    > = LazyLock::new(|| {
        let mut result = [[[None; Side::COUNT]; Side::COUNT]; Side::COUNT];
        let mut vertex = Vertex::iter();
        // Kind of a hack, but working this out by hand isn't any fun.
        for a in 0..Side::COUNT {
            for b in (a + 1)..Side::COUNT {
                for c in (b + 1)..Side::COUNT {
                    if !Side::from_index(a).adjacent_to(Side::from_index(b))
                        || !Side::from_index(b).adjacent_to(Side::from_index(c))
                        || !Side::from_index(c).adjacent_to(Side::from_index(a))
                    {
                        continue;
                    }
                    let v = Some(vertex.next().unwrap());
                    result[a][b][c] = v;
                    result[a][c][b] = v;
                    result[b][a][c] = v;
                    result[b][c][a] = v;
                    result[c][a][b] = v;
                    result[c][b][a] = v;
                }
            }
        }
        assert_eq!(vertex.next(), None);
        result
    });

    /// Whether the determinant of the cube-to-node transform is negative
    pub static CHUNK_TO_NODE_PARITY: LazyLock<[bool; Vertex::COUNT]> = LazyLock::new(|| {
        let mut result = [false; Vertex::COUNT];

        for v in Vertex::iter() {
            result[v as usize] = v.dual_to_node().parity();
        }
        result
    });

    pub static SIDE_NORMALS_F32: LazyLock<[MVector<f32>; Side::COUNT]> =
        LazyLock::new(|| SIDE_NORMALS_F64.map(|n| n.to_f32()));

    pub static REFLECTIONS_F32: LazyLock<[MIsometry<f32>; Side::COUNT]> =
        LazyLock::new(|| REFLECTIONS_F64.map(|n| n.to_f32()));

    pub static DUAL_TO_NODE_F32: LazyLock<[MIsometry<f32>; Vertex::COUNT]> =
        LazyLock::new(|| DUAL_TO_NODE_F64.map(|n| n.to_f32()));

    pub static NODE_TO_DUAL_F32: LazyLock<[MIsometry<f32>; Vertex::COUNT]> =
        LazyLock::new(|| NODE_TO_DUAL_F64.map(|n| n.to_f32()));

    pub static DUAL_TO_CHUNK_FACTOR_F32: LazyLock<f32> =
        LazyLock::new(|| *DUAL_TO_CHUNK_FACTOR_F64 as f32);

    pub static CHUNK_TO_DUAL_FACTOR_F32: LazyLock<f32> =
        LazyLock::new(|| *CHUNK_TO_DUAL_FACTOR_F64 as f32);
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;

    #[test]
    fn vertex_sides_consistent() {
        use std::collections::HashSet;
        let triples = VERTEX_SIDES.iter().collect::<HashSet<_>>();
        assert_eq!(triples.len(), Vertex::COUNT);
        for &triple in VERTEX_SIDES.iter() {
            let mut sorted = triple;
            sorted.sort_unstable();
            assert_eq!(triple, sorted);
            assert!(ADJACENT[triple[0] as usize][triple[1] as usize]);
            assert!(ADJACENT[triple[1] as usize][triple[2] as usize]);
            assert!(ADJACENT[triple[2] as usize][triple[0] as usize]);
        }
    }

    #[test]
    fn sides_to_vertex() {
        for v in Vertex::iter() {
            let [a, b, c] = v.canonical_sides();
            assert_eq!(v, Vertex::from_sides(a, b, c).unwrap());
            assert_eq!(v, Vertex::from_sides(a, c, b).unwrap());
            assert_eq!(v, Vertex::from_sides(b, a, c).unwrap());
            assert_eq!(v, Vertex::from_sides(b, c, a).unwrap());
            assert_eq!(v, Vertex::from_sides(c, a, b).unwrap());
            assert_eq!(v, Vertex::from_sides(c, b, a).unwrap());
        }
    }

    #[test]
    fn adjacent_chunk_axis_permutations() {
        // Assumptions for this test to be valid. If any assertions in this section fail, the test itself
        // needs to be modified
        assert_eq!(Vertex::A.canonical_sides(), [Side::A, Side::B, Side::C]);
        assert_eq!(Vertex::B.canonical_sides(), [Side::A, Side::B, Side::E]);

        assert_eq!(Vertex::F.canonical_sides(), [Side::B, Side::C, Side::F]);
        assert_eq!(Vertex::J.canonical_sides(), [Side::C, Side::F, Side::H]);

        // Test cases

        // Variables with name vertex_?_canonical_sides_reflected refer to the canonical sides
        // of a particular vertex after a reflection that moves it to another vertex.
        // For instance, vertex_a_canonical_sides_reflected is similar to Vertex::A.canonical_sides(),
        // but one of the sides is changed to match Vertex B, but the order of the other two sides is left alone.
        let vertex_a_canonical_sides_reflected = [Side::A, Side::B, Side::E];
        let vertex_b_canonical_sides_reflected = [Side::A, Side::B, Side::C];
        assert_eq!(
            Vertex::A.chunk_axis_permutations()[2],
            ChunkAxisPermutation::from_permutation(
                vertex_a_canonical_sides_reflected,
                Vertex::B.canonical_sides()
            )
        );
        assert_eq!(
            Vertex::B.chunk_axis_permutations()[2],
            ChunkAxisPermutation::from_permutation(
                vertex_b_canonical_sides_reflected,
                Vertex::A.canonical_sides()
            )
        );

        let vertex_f_canonical_sides_reflected = [Side::H, Side::C, Side::F];
        let vertex_j_canonical_sides_reflected = [Side::C, Side::F, Side::B];
        assert_eq!(
            Vertex::F.chunk_axis_permutations()[0],
            ChunkAxisPermutation::from_permutation(
                vertex_f_canonical_sides_reflected,
                Vertex::J.canonical_sides()
            )
        );
        assert_eq!(
            Vertex::J.chunk_axis_permutations()[2],
            ChunkAxisPermutation::from_permutation(
                vertex_j_canonical_sides_reflected,
                Vertex::F.canonical_sides()
            )
        );
    }

    #[test]
    fn side_is_facing() {
        for side in Side::iter() {
            assert!(!side.is_facing(&MVector::origin()));
            assert!(side.is_facing(&(*side.reflection() * MVector::origin())));
        }
    }

    #[test]
    fn radius() {
        let corner = *Vertex::A.dual_to_node_f64() * MVector::origin();
        assert_abs_diff_eq!(
            BOUNDING_SPHERE_RADIUS_F64,
            corner.distance(&MVector::origin()),
            epsilon = 1e-10
        );
        let phi = (1.0 + 5.0f64.sqrt()) / 2.0; // Golden ratio
        assert_abs_diff_eq!(
            BOUNDING_SPHERE_RADIUS_F64,
            (1.5 * phi).sqrt().asinh(),
            epsilon = 1e-10
        );
    }

    #[test]
    fn chunk_to_node() {
        // Chunk coordinates of (1, 1, 1) should be at the center of a dodecahedron.
        let mut chunk_corner_in_node_coordinates =
            Vertex::A.chunk_to_node_f64() * na::Vector4::new(1.0, 1.0, 1.0, 1.0);
        chunk_corner_in_node_coordinates /= chunk_corner_in_node_coordinates.w;
        assert_abs_diff_eq!(
            chunk_corner_in_node_coordinates,
            na::Vector4::new(0.0, 0.0, 0.0, 1.0),
            epsilon = 1e-10
        );
    }

    #[test]
    fn node_to_chunk() {
        assert_abs_diff_eq!(
            Vertex::A.chunk_to_node_f64().try_inverse().unwrap(),
            Vertex::A.node_to_chunk_f64(),
            epsilon = 1e-10
        );
    }
}
