//! Tools for processing the geometry of a right dodecahedron

use std::ops::{Index, IndexMut};

use data::*;
use enum_map::{Enum, EnumMap};
use serde::{Deserialize, Serialize};

use crate::math::{self, MIsometry, MVector};
use crate::voxel_math::ChunkAxisPermutation;

pub const BOUNDING_SPHERE_RADIUS_F64: f64 = 1.2264568712514068;
pub const BOUNDING_SPHERE_RADIUS: f32 = BOUNDING_SPHERE_RADIUS_F64 as f32;

/// Sides of a right dodecahedron
#[derive(
    Debug, Copy, Clone, Enum, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize,
)]
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

    pub fn iter() -> impl ExactSizeIterator<Item = Self> {
        Self::VALUES.into_iter()
    }

    /// Whether `self` and `other` share an edge
    ///
    /// `false` when `self == other`.
    #[inline]
    pub fn adjacent_to(self, other: Side) -> bool {
        adjacent()[self as usize][other as usize]
    }

    /// Outward normal vector of this side
    #[inline]
    pub fn normal(self) -> &'static MVector<f32> {
        &side_normals_f32()[self as usize]
    }

    /// Outward normal vector of this side
    #[inline]
    pub fn normal_f64(self) -> &'static MVector<f64> {
        &side_normals_f64()[self as usize]
    }

    /// Reflection across this side
    #[inline]
    pub fn reflection(self) -> &'static MIsometry<f32> {
        &reflections_f32()[self as usize]
    }

    /// Reflection across this side
    #[inline]
    pub fn reflection_f64(self) -> &'static MIsometry<f64> {
        &reflections_f64()[self as usize]
    }

    /// Whether `p` is opposite the dodecahedron across the plane containing `self`
    #[inline]
    pub fn is_facing(self, p: &MVector<f32>) -> bool {
        let r = na::convert::<_, na::RowVector4<f32>>(self.reflection().row(3).clone_owned());
        (r * na::Vector4::from(*p)).x < p.w
    }
}

/// Vertices of a right dodecahedron
#[derive(Debug, Copy, Clone, Enum, Eq, PartialEq, Hash, Serialize, Deserialize)]
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
        Self::VALUES.into_iter()
    }

    /// Vertex shared by three sides, if any
    #[inline]
    pub fn from_sides(a: Side, b: Side, c: Side) -> Option<Self> {
        sides_to_vertex()[a as usize][b as usize][c as usize]
    }

    /// Sides incident to this vertex, in canonical order.
    ///
    /// This canonical order determines the X, Y, and Z axes of the chunk
    /// corresponding to the vertex.
    #[inline]
    pub fn canonical_sides(self) -> [Side; 3] {
        vertex_sides()[self as usize]
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
        adjacent_vertices()[self as usize]
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
        &chunk_axis_permutations()[self as usize]
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
        &dual_to_node_f32()[self as usize]
    }

    /// Transform from cube-centric coordinates to dodeca-centric coordinates
    pub fn dual_to_node_f64(self) -> &'static MIsometry<f64> {
        &dual_to_node_f64()[self as usize]
    }

    /// Transform from dodeca-centric coordinates to cube-centric coordinates
    pub fn node_to_dual(self) -> &'static MIsometry<f32> {
        &node_to_dual_f32()[self as usize]
    }

    /// Transform from dodeca-centric coordinates to cube-centric coordinates
    pub fn node_to_dual_f64(self) -> &'static MIsometry<f64> {
        &node_to_dual_f64()[self as usize]
    }

    /// Scale factor used in conversion from cube-centric coordinates to euclidean chunk coordinates.
    /// Scaling the x, y, and z components of a vector in cube-centric coordinates by this value
    /// and dividing them by the w coordinate will yield euclidean chunk coordinates.
    pub fn dual_to_chunk_factor() -> f32 {
        dual_to_chunk_factor_f32()
    }

    /// Scale factor used in conversion from cube-centric coordinates to euclidean chunk coordinates.
    /// Scaling the x, y, and z components of a vector in cube-centric coordinates by this value
    /// and dividing them by the w coordinate will yield euclidean chunk coordinates.
    pub fn dual_to_chunk_factor_f64() -> f64 {
        dual_to_chunk_factor_f64()
    }

    /// Scale factor used in conversion from euclidean chunk coordinates to cube-centric coordinates.
    /// Scaling the x, y, and z components of a vector in homogeneous euclidean chunk coordinates by this value
    /// and lorentz-normalizing the result will yield cube-centric coordinates.
    pub fn chunk_to_dual_factor() -> f32 {
        chunk_to_dual_factor_f32()
    }

    /// Scale factor used in conversion from euclidean chunk coordinates to cube-centric coordinates.
    /// Scaling the x, y, and z components of a vector in homogeneous euclidean chunk coordinates by this value
    /// and lorentz-normalizing the result will yield cube-centric coordinates.
    pub fn chunk_to_dual_factor_f64() -> f64 {
        chunk_to_dual_factor_f64()
    }

    /// Convenience method for `self.chunk_to_node().determinant() < 0`.
    pub fn parity(self) -> bool {
        chunk_to_node_parity()[self as usize]
    }
}

struct HalfEdge(u8);

impl HalfEdge {
    const COUNT: usize = 60;
}

impl From<u8> for HalfEdge {
    fn from(value: u8) -> Self {
        assert!(value < Self::COUNT as u8);
        Self(value)
    }
}

struct HalfEdgeList<T>([T; HalfEdge::COUNT]);

impl<T> Index<HalfEdge> for HalfEdgeList<T> {
    type Output = T;

    fn index(&self, index: HalfEdge) -> &Self::Output {
        &self.0[index.0 as usize]
    }
}

impl<T> IndexMut<HalfEdge> for HalfEdgeList<T> {
    fn index_mut(&mut self, index: HalfEdge) -> &mut Self::Output {
        &mut self.0[index.0 as usize]
    }
}

struct SideData {
    adjacency_map: EnumMap<Side, bool>,

    normal_f64: MVector<f64>,
    reflection_f64: MIsometry<f64>,

    normal_f32: MVector<f32>,
    reflection_f32: MIsometry<f32>,
}

struct VertexData {
    incident_sides: [Side; 3],
    adjacent_vertices: [Vertex; 3],
    chunk_axis_permutations: [ChunkAxisPermutation; 3],

    dual_to_node_f64: MIsometry<f64>,
    node_to_dual_f64: MIsometry<f64>,

    dual_to_node_f32: MIsometry<f32>,
    node_to_dual_f32: MIsometry<f32>,
}

struct HalfEdgeData {
    next: HalfEdge,
    previous: HalfEdge,
    twin: HalfEdge,
    side: Side,
    vertex: Vertex,
}

struct Dodecahedron {
    sides: EnumMap<Side, SideData>,
    vertices: EnumMap<Vertex, VertexData>,
    half_edges: HalfEdgeList<HalfEdgeData>,
}

impl Dodecahedron {
    fn new() -> Self {
        let phi = libm::sqrt(1.25) + 0.5; // golden ratio
        let reference_normal = MVector::new(1.0, phi, 0.0, libm::sqrt(phi)).lorentz_normalize();

        let side_data: EnumMap<Side, _> = EnumMap::from_fn(|side: Side| {
            let rectangle = side.rectangle();
            let sign_options = [1.0, -1.0];
            let short_side_sign = side.short_side_sign();
            let long_side_sign = side.long_side_sign();
            let side_parity_sign = (short_side_sign + long_side_sign) % 2;

            let normal_f64 = math::tuv_to_xyz(
                rectangle,
                MVector::new(
                    reference_normal.x * sign_options[short_side_sign],
                    reference_normal.y * sign_options[long_side_sign],
                    reference_normal.z * sign_options[side_parity_sign],
                    reference_normal.w,
                ),
            );

            let reflection_f64 = normal_f64.reflect();
            /*
            Around (0,  1,  p): (0, -1,  p), ( p, 0,  1), ( 1,  p, 0), (-1,  p, 0), (-p, 0,  1)
            Around (0, -1,  p): (0,  1,  p), (-p, 0,  1), (-1, -p, 0), ( 1, -p, 0), ( p, 0,  1)
            Around (0,  1, -p): (0, -1, -p), (-p, 0, -1), (-1,  p, 0), ( 1,  p, 0), ( p, 0, -1)
            Around (0, -1, -p): (0,  1, -p), ( p, 0, -1), ( 1, -p, 0), (-1, -p, 0), (-p, 0, -1)
            */

            let adjacent_sides = [
                Side::packed_index(rectangle, short_side_sign + 1, long_side_sign),
                Side::packed_index(rectangle + 1, long_side_sign, side_parity_sign),
                Side::packed_index(rectangle + 2, side_parity_sign, short_side_sign),
                Side::packed_index(rectangle + 2, side_parity_sign + 1, short_side_sign),
                Side::packed_index(rectangle + 1, long_side_sign, side_parity_sign + 1),
            ];
            let mut adjacency_map = EnumMap::from_fn(|_| false);
            for adjacent_side in adjacent_sides {
                adjacency_map[adjacent_side] = true;
            }

            SideData {
                adjacency_map,

                normal_f32: normal_f64.to_f32(),
                normal_f64,

                reflection_f32: reflection_f64.to_f32(),
                reflection_f64,
            }
        });

        todo!()
    }
}

impl Side {
    fn packed_index(rectangle: usize, short_side_sign: usize, long_side_sign: usize) -> Self {
        Self::VALUES[(rectangle % 3) + (short_side_sign % 2) * 3 + (long_side_sign % 2) * 6]
    }

    fn rectangle(self) -> usize {
        (self as usize) % 3
    }

    fn short_side_sign(self) -> usize {
        ((self as usize) / 3) % 2
    }

    fn long_side_sign(self) -> usize {
        ((self as usize) / 6) % 2
    }
}

mod data {
    use std::array;
    use std::sync::OnceLock;

    use crate::dodeca::{Side, Vertex};
    use crate::math::{MIsometry, MVector};
    use crate::voxel_math::ChunkAxisPermutation;

    /// Whether two sides share an edge
    pub fn adjacent() -> &'static [[bool; Side::COUNT]; Side::COUNT] {
        static LOCK: OnceLock<[[bool; Side::COUNT]; Side::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| {
            let mut result = [[false; Side::COUNT]; Side::COUNT];
            for (i, side) in result.iter_mut().enumerate() {
                for (j, is_adjacent) in side.iter_mut().enumerate() {
                    let cosh_distance = (reflections_f64()[i] * reflections_f64()[j])[(3, 3)];
                    // Possile cosh_distances: 1, 4.23606 = 2+sqrt(5), 9.47213 = 5+2*sqrt(5), 12.70820 = 6+3*sqrt(5);
                    // < 2.0 indicates identical faces; < 5.0 indicates adjacent faces; > 5.0 indicates non-adjacent faces
                    *is_adjacent = (2.0..5.0).contains(&cosh_distance);
                }
            }
            result
        })
    }

    /// Vector corresponding to the outer normal of each side
    pub fn side_normals_f64() -> &'static [MVector<f64>; Side::COUNT] {
        static LOCK: OnceLock<[MVector<f64>; Side::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| {
            let phi = libm::sqrt(1.25) + 0.5; // golden ratio
            let f = MVector::new(1.0, phi, 0.0, libm::sqrt(phi)).lorentz_normalize();

            let mut result: [MVector<f64>; Side::COUNT] = [MVector::zero(); Side::COUNT];
            let mut i = 0;
            for (x, y, z, w) in [
                (f.x, f.y, f.z, f.w),
                (-f.x, f.y, -f.z, f.w),
                (f.x, -f.y, -f.z, f.w),
                (-f.x, -f.y, f.z, f.w),
            ] {
                for (x, y, z, w) in [(x, y, z, w), (y, z, x, w), (z, x, y, w)] {
                    result[i] = MVector::new(x, y, z, w);
                    i += 1;
                }
            }
            result
        })
    }

    /// Transform that moves from a neighbor to a reference node, for each side
    pub fn reflections_f64() -> &'static [MIsometry<f64>; Side::COUNT] {
        static LOCK: OnceLock<[MIsometry<f64>; Side::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| side_normals_f64().map(|r| r.reflect()))
    }

    /// Sides incident to a vertex, in canonical order
    pub fn vertex_sides() -> &'static [[Side; 3]; Vertex::COUNT] {
        static LOCK: OnceLock<[[Side; 3]; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| {
            let mut result = [[Side::A; 3]; Vertex::COUNT];
            let mut vertex = 0;
            // Kind of a hack, but working this out by hand isn't any fun.
            for a in 0..Side::COUNT {
                for b in (a + 1)..Side::COUNT {
                    for c in (b + 1)..Side::COUNT {
                        if !adjacent()[a][b] || !adjacent()[b][c] || !adjacent()[c][a] {
                            continue;
                        }
                        result[vertex] = [Side::VALUES[a], Side::VALUES[b], Side::VALUES[c]];
                        vertex += 1;
                    }
                }
            }
            assert_eq!(vertex, 20);
            result
        })
    }

    // Which vertices are adjacent to other vertices and opposite the canonical sides
    pub fn adjacent_vertices() -> &'static [[Vertex; 3]; Vertex::COUNT] {
        static LOCK: OnceLock<[[Vertex; 3]; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| {
            let mut result = [[Vertex::A; 3]; Vertex::COUNT];

            for (i, triple) in result.iter_mut().enumerate() {
                for result_index in 0..3 {
                    let mut test_sides = vertex_sides()[i];
                    // Keep modifying the result_index'th element of test_sides until its three elements are all
                    // adjacent to a single vertex. That vertex is the vertex we're looking for.
                    for side in Side::iter() {
                        if side == vertex_sides()[i][result_index] {
                            continue;
                        }
                        test_sides[result_index] = side;
                        if let Some(adjacent_vertex) =
                            Vertex::from_sides(test_sides[0], test_sides[1], test_sides[2])
                        {
                            triple[result_index] = adjacent_vertex;
                        }
                    }
                }
            }
            result
        })
    }

    // Which transformations have to be done after a reflection to switch reference frames from one vertex
    // to one of its adjacent vertices (ordered similarly to ADJACENT_VERTICES)
    pub fn chunk_axis_permutations() -> &'static [[ChunkAxisPermutation; 3]; Vertex::COUNT] {
        static LOCK: OnceLock<[[ChunkAxisPermutation; 3]; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| {
            array::from_fn(|vertex| {
                array::from_fn(|result_index| {
                    let mut test_sides = vertex_sides()[vertex];
                    // Keep modifying the result_index'th element of test_sides until its three elements are all
                    // adjacent to a single vertex (determined using `Vertex::from_sides`).
                    for side in Side::iter() {
                        if side == vertex_sides()[vertex][result_index] {
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
        })
    }

    /// Transform that converts from cube-centric coordinates to dodeca-centric coordinates
    pub fn dual_to_node_f64() -> &'static [MIsometry<f64>; Vertex::COUNT] {
        static LOCK: OnceLock<[MIsometry<f64>; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| {
            let mip_origin_normal = MVector::origin().mip(&side_normals_f64()[0]); // This value is the same for every side
            let mut result = [MIsometry::identity(); Vertex::COUNT];
            for (i, map) in result.iter_mut().enumerate() {
                let [a, b, c] = vertex_sides()[i];
                let vertex_position = (MVector::origin()
                    - (*a.normal_f64() + *b.normal_f64() + *c.normal_f64()) * mip_origin_normal)
                    .lorentz_normalize();
                *map = MIsometry::from_columns_unchecked(&[
                    -*a.normal_f64(),
                    -*b.normal_f64(),
                    -*c.normal_f64(),
                    vertex_position,
                ]);
            }
            result
        })
    }

    /// Transform that converts from dodeca-centric coordinates to cube-centric coordinates
    pub fn node_to_dual_f64() -> &'static [MIsometry<f64>; Vertex::COUNT] {
        static LOCK: OnceLock<[MIsometry<f64>; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| dual_to_node_f64().map(|m| m.mtranspose()))
    }

    pub fn dual_to_chunk_factor_f64() -> f64 {
        static LOCK: OnceLock<f64> = OnceLock::new();
        *LOCK.get_or_init(|| (2.0 + 5.0f64.sqrt()).sqrt())
    }

    pub fn chunk_to_dual_factor_f64() -> f64 {
        static LOCK: OnceLock<f64> = OnceLock::new();
        *LOCK.get_or_init(|| 1.0 / dual_to_chunk_factor_f64())
    }

    /// Vertex shared by 3 sides
    pub fn sides_to_vertex() -> &'static [[[Option<Vertex>; Side::COUNT]; Side::COUNT]; Side::COUNT]
    {
        static LOCK: OnceLock<[[[Option<Vertex>; Side::COUNT]; Side::COUNT]; Side::COUNT]> =
            OnceLock::new();
        LOCK.get_or_init(|| {
            let mut result = [[[None; Side::COUNT]; Side::COUNT]; Side::COUNT];
            let mut vertex = Vertex::iter();
            // Kind of a hack, but working this out by hand isn't any fun.
            for a in 0..Side::COUNT {
                for b in (a + 1)..Side::COUNT {
                    for c in (b + 1)..Side::COUNT {
                        if !Side::VALUES[a].adjacent_to(Side::VALUES[b])
                            || !Side::VALUES[b].adjacent_to(Side::VALUES[c])
                            || !Side::VALUES[c].adjacent_to(Side::VALUES[a])
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
        })
    }

    /// Whether the determinant of the cube-to-node transform is negative
    pub fn chunk_to_node_parity() -> &'static [bool; Vertex::COUNT] {
        static LOCK: OnceLock<[bool; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| {
            let mut result = [false; Vertex::COUNT];

            for v in Vertex::iter() {
                result[v as usize] = v.dual_to_node().parity();
            }
            result
        })
    }

    pub fn side_normals_f32() -> &'static [MVector<f32>; Side::COUNT] {
        static LOCK: OnceLock<[MVector<f32>; Side::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| side_normals_f64().map(|n| n.to_f32()))
    }

    pub fn reflections_f32() -> &'static [MIsometry<f32>; Side::COUNT] {
        static LOCK: OnceLock<[MIsometry<f32>; Side::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| reflections_f64().map(|n| n.to_f32()))
    }

    pub fn dual_to_node_f32() -> &'static [MIsometry<f32>; Vertex::COUNT] {
        static LOCK: OnceLock<[MIsometry<f32>; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| dual_to_node_f64().map(|n| n.to_f32()))
    }

    pub fn node_to_dual_f32() -> &'static [MIsometry<f32>; Vertex::COUNT] {
        static LOCK: OnceLock<[MIsometry<f32>; Vertex::COUNT]> = OnceLock::new();
        LOCK.get_or_init(|| node_to_dual_f64().map(|n| n.to_f32()))
    }

    pub fn dual_to_chunk_factor_f32() -> f32 {
        static LOCK: OnceLock<f32> = OnceLock::new();
        *LOCK.get_or_init(|| dual_to_chunk_factor_f64() as f32)
    }

    pub fn chunk_to_dual_factor_f32() -> f32 {
        static LOCK: OnceLock<f32> = OnceLock::new();
        *LOCK.get_or_init(|| chunk_to_dual_factor_f64() as f32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math;
    use approx::*;

    #[test]
    fn vertex_sides_consistent() {
        use std::collections::HashSet;
        let triples = vertex_sides().iter().collect::<HashSet<_>>();
        assert_eq!(triples.len(), Vertex::COUNT);
        for &triple in vertex_sides() {
            let mut sorted = triple;
            sorted.sort_unstable();
            assert_eq!(triple, sorted);
            assert!(adjacent()[triple[0] as usize][triple[1] as usize]);
            assert!(adjacent()[triple[1] as usize][triple[2] as usize]);
            assert!(adjacent()[triple[2] as usize][triple[0] as usize]);
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
            math::distance(&corner, &MVector::origin()),
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

    struct ManifoldFinder {
        oriented_sides: [[Side; 5]; 12],
        gluings: [Option<(Side, usize)>; 12],
    }

    impl ManifoldFinder {
        fn new() -> Self {
            let mut oriented_sides: [[Side; 5]; 12] = [[Side::A; 5]; 12];
            for side in Side::iter() {
                for adjacent_side in Side::iter() {
                    if !side.adjacent_to(adjacent_side) {
                        continue;
                    }
                    // Use adjacent_side
                    oriented_sides[side as usize][0] = adjacent_side;
                    let mut current_side = adjacent_side;
                    for i in 1..5 {
                        // Find appropriate vertex
                        for vertex in Vertex::iter() {
                            let Some(first_index) =
                                vertex.canonical_sides().iter().position(|s| *s == side)
                            else {
                                continue;
                            };
                            let Some(second_index) = vertex
                                .canonical_sides()
                                .iter()
                                .position(|s| *s == current_side)
                            else {
                                continue;
                            };
                            let diff = (second_index + 3 - first_index) % 3;
                            // Make sure we're going the right way
                            if (diff == 2) != vertex.parity() {
                                continue;
                            }
                            let Some(third_index) = vertex
                                .canonical_sides()
                                .iter()
                                .position(|s| *s != side && *s != current_side)
                            else {
                                panic!();
                            };
                            current_side = vertex.canonical_sides()[third_index];
                            oriented_sides[side as usize][i] = current_side;
                            break;
                        }
                    }
                    break;
                }
            }
            ManifoldFinder {
                oriented_sides,
                gluings: Default::default(),
            }
        }

        fn check_edge(side: Side, side_edge: usize) {
            // TODO: The dodeca data structure isn't very good for this.
        }
    }

    #[test]
    fn find_manifold() {
        let mut manifold_finder = ManifoldFinder::new();
        println!("{:?}", manifold_finder.oriented_sides);
    }
}
