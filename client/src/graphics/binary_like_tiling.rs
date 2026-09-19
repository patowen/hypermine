use ash::vk;
use common::{
    graph::{Graph, NodeId},
    math::{MIsometry, MPoint, MVector, PermuteXYZ, sqr},
    proto::Position,
    worldgen,
};
use fxhash::FxHashMap;
use libm::{coshf, logf, powf, sinhf, sqrtf, tanhf};

use crate::graphics::{
    Mesh, Meshes,
    asset_loader::AssetLoadContext,
    meshes::{MeshGeometryDefinition, Vertex},
};

fn pseudo_chunk_to_mvector_simple(pseudo_chunk: na::Vector3<f32>) -> MVector<f32> {
    let factor = 1.0 / sqrtf(1.0 - sqr(pseudo_chunk.x) - sqr(pseudo_chunk.y));
    //MVector::new(voxel.x, voxel.y, factor * tanhf(voxel.z), 1.0)

    // This is already pre-scaled
    MVector::new(
        pseudo_chunk.x * coshf(pseudo_chunk.z) * factor,
        pseudo_chunk.y * coshf(pseudo_chunk.z) * factor,
        sinhf(pseudo_chunk.z),
        coshf(pseudo_chunk.z) * factor,
    )
}

fn pseudo_chunk_to_isometry(pseudo_chunk: na::Vector3<f32>, boost: f32) -> MIsometry<f32> {
    let w = pseudo_chunk_to_mvector_boosted(pseudo_chunk, boost).to_point_unchecked();
    let x = pseudo_chunk_to_mvector_boosted_partial_x(pseudo_chunk, boost).normalized_direction();
    let z = pseudo_chunk_to_mvector_boosted_partial_z(pseudo_chunk, boost).to_direction_unchecked();
    // TODO: There might be a better way to get `y`, especially since voxel_to_mvector_boosted_partial_y will be
    // in the wrong direction.
    let mut y = pseudo_chunk_to_mvector_boosted_partial_y(pseudo_chunk, boost);
    y -= MVector::from(x) * y.mip(&x);
    let y = y.normalized_direction();
    MIsometry::from_columns_unchecked(&[x, y, z], w)
}

fn chunk_to_isometry(
    klein_coords: na::Vector2<f32>,
    chunk: na::Vector3<f32>,
    boost: f32,
) -> MIsometry<f32> {
    let pseudo_chunk = chunk_to_pseudo_chunk(klein_coords, chunk, boost);
    let pseudo_chunk_partial_x = chunk_to_pseudo_chunk_partial_x(klein_coords, chunk, boost);
    let w = pseudo_chunk_to_mvector_boosted(pseudo_chunk, boost).to_point_unchecked();
    let x = (pseudo_chunk_to_mvector_boosted_partial_x(pseudo_chunk, boost)
        * pseudo_chunk_partial_x.x
        + pseudo_chunk_to_mvector_boosted_partial_y(pseudo_chunk, boost)
            * pseudo_chunk_partial_x.y)
        .normalized_direction();
    let z = pseudo_chunk_to_mvector_boosted_partial_z(pseudo_chunk, boost).to_direction_unchecked();
    let mut y = pseudo_chunk_to_mvector_boosted_partial_y(pseudo_chunk, boost);
    y -= MVector::from(x) * y.mip(&x);
    let y = y.normalized_direction();
    MIsometry::from_columns_unchecked(&[x, y, z], w)
}

fn chunk_to_mvector_simple(
    klein_coords: na::Vector2<f32>,
    chunk: na::Vector3<f32>,
    boost: f32,
) -> MVector<f32> {
    pseudo_chunk_to_isometry(
        na::Vector3::new(
            klein_coords.x * coshf(boost),
            klein_coords.y * coshf(boost),
            0.0,
        ),
        boost,
    )
    .inverse()
        * pseudo_chunk_to_mvector_boosted(
            na::Vector3::new(
                klein_coords.x * coshf(boost) + chunk.x,
                klein_coords.y * coshf(boost) + chunk.y,
                chunk.z,
            ),
            boost,
        )
}

fn chunk_to_pseudo_chunk(
    klein_coords: na::Vector2<f32>,
    chunk: na::Vector3<f32>,
    boost: f32,
) -> na::Vector3<f32> {
    let origin = na::Vector3::new(klein_coords[0], klein_coords[1], 1.0);
    let origin_norm = sqrtf(-(sqr(origin[0]) + sqr(origin[1]) - sqr(origin[2])));
    let normalized_origin = origin / origin_norm;
    let z = normalized_origin;
    let mut y = z.cross(&na::Vector3::x());
    y.z *= -1.0;
    y /= sqrtf(sqr(y.x) + sqr(y.y) - sqr(y.z));
    let mut x = y.cross(&z);
    x.z *= -1.0;
    let mut conversion = na::Matrix3::from_columns(&[x, y, z]).try_inverse().unwrap();
    conversion *= na::Matrix3::new_translation(&klein_coords); // Applying skew

    // This is equivalent (but numerically more stable) to `new_scale(coshf(boost)) * conversion * new_scale(1.0/coshf(boost))`
    conversion[(2, 0)] /= coshf(boost);
    conversion[(2, 1)] /= coshf(boost);

    let horizontal_coords = conversion * chunk.xy().push(1.0);
    na::Vector3::new(
        horizontal_coords.x / horizontal_coords.z,
        horizontal_coords.y / horizontal_coords.z,
        chunk.z,
    )
}

fn chunk_to_pseudo_chunk_partial_x(
    klein_coords: na::Vector2<f32>,
    chunk: na::Vector3<f32>,
    boost: f32,
) -> na::Vector3<f32> {
    let origin = na::Vector3::new(klein_coords[0], klein_coords[1], 1.0);
    let origin_norm = sqrtf(-(sqr(origin[0]) + sqr(origin[1]) - sqr(origin[2])));
    let normalized_origin = origin / origin_norm;
    let z = normalized_origin;
    let mut y = z.cross(&na::Vector3::x());
    y.z *= -1.0;
    y /= sqrtf(sqr(y.x) + sqr(y.y) - sqr(y.z));
    let mut x = y.cross(&z);
    x.z *= -1.0;
    let mut conversion = na::Matrix3::from_columns(&[x, y, z]).try_inverse().unwrap();
    conversion *= na::Matrix3::new_translation(&klein_coords); // Applying skew

    // This is equivalent (but numerically more stable) to `new_scale(coshf(boost)) * conversion * new_scale(1.0/coshf(boost))`
    conversion[(2, 0)] /= coshf(boost);
    conversion[(2, 1)] /= coshf(boost);

    let horizontal_coords = conversion * chunk.xy().push(1.0);
    let horizontal_coords_partial_x = conversion * na::Vector3::x();
    na::Vector3::new(
        (horizontal_coords_partial_x.x * horizontal_coords.z
            - horizontal_coords.x * horizontal_coords_partial_x.z)
            / sqr(horizontal_coords.z),
        (horizontal_coords_partial_x.y * horizontal_coords.z
            - horizontal_coords.y * horizontal_coords_partial_x.z)
            / sqr(horizontal_coords.z),
        0.0,
    )
}

/// Computes
/// `translation_along([0, 0, -boost]) * voxel_to_mvector_simple([x / cosh(boost), y / cosh(boost), boost + z])`
fn pseudo_chunk_to_mvector_boosted(pseudo_chunk: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        pseudo_chunk.x * factor * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        pseudo_chunk.y * factor * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        sinhf(pseudo_chunk.z) * (1.0 - scaled_factor_delta * sqr(tanhf(boost)))
            - scaled_factor_delta * tanhf(boost) * coshf(pseudo_chunk.z),
        coshf(pseudo_chunk.z) * (1.0 + scaled_factor_delta)
            + scaled_factor_delta * tanhf(boost) * sinhf(pseudo_chunk.z),
    )
}

fn pseudo_chunk_to_mvector_boosted_partial_x(
    pseudo_chunk: na::Vector3<f32>,
    boost: f32,
) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let dist_squared_partial_x = pseudo_chunk.x * 2.0;
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor_radicand_partial_x = -dist_squared_partial_x / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let factor_partial_x = -0.5 * factor_radicand_partial_x * factor / factor_radicand;
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    let scaled_factor_delta_partial_x = (dist_squared_partial_x
        * (factor_radicand + sqrtf(factor_radicand))
        - dist_squared
            * (factor_radicand_partial_x
                + 0.5 / sqrtf(factor_radicand) * factor_radicand_partial_x))
        / sqr(factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        (factor + pseudo_chunk.x * factor_partial_x)
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        pseudo_chunk.y
            * factor_partial_x
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        sinhf(pseudo_chunk.z) * (-scaled_factor_delta_partial_x * sqr(tanhf(boost)))
            - scaled_factor_delta_partial_x * tanhf(boost) * coshf(pseudo_chunk.z),
        coshf(pseudo_chunk.z) * (scaled_factor_delta_partial_x)
            + scaled_factor_delta_partial_x * tanhf(boost) * sinhf(pseudo_chunk.z),
    )
}

fn pseudo_chunk_to_mvector_boosted_partial_y(
    pseudo_chunk: na::Vector3<f32>,
    boost: f32,
) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let dist_squared_partial_y = pseudo_chunk.y * 2.0;
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor_radicand_partial_y = -dist_squared_partial_y / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let factor_partial_y = -0.5 * factor_radicand_partial_y * factor / factor_radicand;
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    let scaled_factor_delta_partial_y = (dist_squared_partial_y
        * (factor_radicand + sqrtf(factor_radicand))
        - dist_squared
            * (factor_radicand_partial_y
                + 0.5 / sqrtf(factor_radicand) * factor_radicand_partial_y))
        / sqr(factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        pseudo_chunk.x
            * factor_partial_y
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        (factor + pseudo_chunk.y * factor_partial_y)
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        sinhf(pseudo_chunk.z) * (-scaled_factor_delta_partial_y * sqr(tanhf(boost)))
            - scaled_factor_delta_partial_y * tanhf(boost) * coshf(pseudo_chunk.z),
        coshf(pseudo_chunk.z) * (scaled_factor_delta_partial_y)
            + scaled_factor_delta_partial_y * tanhf(boost) * sinhf(pseudo_chunk.z),
    )
}

fn pseudo_chunk_to_mvector_boosted_partial_z(
    pseudo_chunk: na::Vector3<f32>,
    boost: f32,
) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        pseudo_chunk.x * factor * (sinhf(pseudo_chunk.z) + coshf(pseudo_chunk.z) * tanhf(boost)),
        pseudo_chunk.y * factor * (sinhf(pseudo_chunk.z) + coshf(pseudo_chunk.z) * tanhf(boost)),
        coshf(pseudo_chunk.z) * (1.0 - scaled_factor_delta * sqr(tanhf(boost)))
            - scaled_factor_delta * tanhf(boost) * sinhf(pseudo_chunk.z),
        sinhf(pseudo_chunk.z) * (1.0 + scaled_factor_delta)
            + scaled_factor_delta * tanhf(boost) * coshf(pseudo_chunk.z),
    )
}

fn add_quad(
    chunk: &BltChunk,
    layout: &BltLayout,
    transform: &MIsometry<f32>,
    geometry: &mut MeshGeometryDefinition,
    points: [na::Vector3<i32>; 4],
    texture: usize,
) {
    let vertices: Vec<_> = points
        .into_iter()
        .enumerate()
        .map(|(i, point)| {
            let len = geometry.vertices.len();
            geometry.vertices.push(Vertex {
                position: common::dodeca::Side::A.reflection()
                    * common::dodeca::Vertex::A.dual_to_node()
                    * (transform * chunk.point_from_voxel(layout, point.cast())).tuv_to_xyz(1),
                texcoords: na::Vector3::new((i & 1) as f32, ((i >> 1) & 1) as f32, texture as f32),
                normal: common::math::MDirection::x(),
            });
            len as u32
        })
        .collect();
    geometry.indices.extend(&[
        vertices[0],
        vertices[1],
        vertices[2],
        vertices[1],
        vertices[3],
        vertices[2],
    ]);
}

fn add_voxel(
    chunk: &BltChunk,
    layout: &BltLayout,
    transform: &MIsometry<f32>,
    geometry: &mut MeshGeometryDefinition,
    coords: na::Vector3<i32>,
) {
    for x_axis in 0..3 {
        let t = na::Vector3::x().tuv_to_xyz(x_axis);
        let u = na::Vector3::y().tuv_to_xyz(x_axis);
        let v = na::Vector3::z().tuv_to_xyz(x_axis);
        add_quad(
            chunk,
            layout,
            transform,
            geometry,
            [coords, coords + t, coords + u, coords + t + u],
            x_axis,
        );
        add_quad(
            chunk,
            layout,
            transform,
            geometry,
            [
                coords + v,
                coords + u + v,
                coords + t + v,
                coords + t + u + v,
            ],
            x_axis,
        );
    }
}

pub struct BltChunkSurface {
    pub geometry: MeshGeometryDefinition,
}

impl skid_steer::Source for BltChunkSurface {
    type Output = Mesh;

    async fn load<'a>(self, context: &'a skid_steer::Context<'a>) -> Option<Mesh> {
        let ctx: &AssetLoadContext = context.get().unwrap();

        let colors = ctx.load_cached(crate::graphics::PngArray {
            path: "materials".into(),
            size: common::world::Material::COUNT - 1,
        });

        Mesh::from_definition(ctx, self.geometry, colors).await
    }

    fn free(mut output: Self::Output, context: &skid_steer::Context) {
        let ctx: &AssetLoadContext = context.get().unwrap();
        unsafe { output.destroy(ctx.device()) };
    }
}

#[derive(Clone, Copy)]
struct BltChunkId(u32);

impl std::fmt::Debug for BltChunkId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self.0, f)
    }
}

#[derive(Clone, Copy)]
struct QuadIndex(u8);

impl QuadIndex {
    fn x(self) -> u8 {
        self.0 & 1
    }

    fn y(self) -> u8 {
        (self.0 >> 1) & 1
    }

    pub const VALUES: [Self; 4] = [QuadIndex(0), QuadIndex(1), QuadIndex(2), QuadIndex(3)];

    fn neighbor(self, side_index: SideIndex) -> QuadIndexNeighbor {
        QuadIndexNeighbor {
            index: QuadIndex(self.0 ^ (1 << side_index.coordinate())),
            different_inner_chunk: (self.0 >> side_index.coordinate()) & 1 == side_index.extreme(),
            return_side: side_index.opposite(),
        }
    }
}

impl std::fmt::Debug for QuadIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self.0, f)
    }
}

struct QuadIndexNeighbor {
    index: QuadIndex,
    different_inner_chunk: bool,
    return_side: SideIndex,
}

#[derive(Default)]
struct QuadIndexMap<T>([T; 4]);

impl<T> std::ops::Index<QuadIndex> for QuadIndexMap<T> {
    type Output = T;

    fn index(&self, index: QuadIndex) -> &Self::Output {
        &self.0[index.0 as usize]
    }
}

impl<T> std::ops::IndexMut<QuadIndex> for QuadIndexMap<T> {
    fn index_mut(&mut self, index: QuadIndex) -> &mut Self::Output {
        &mut self.0[index.0 as usize]
    }
}

impl<T: std::fmt::Debug> std::fmt::Debug for QuadIndexMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self.0, f)
    }
}

#[derive(Clone, Copy)]
struct SideIndex(u8);

impl SideIndex {
    fn coordinate(self) -> u8 {
        self.0 & 1
    }

    fn extreme(self) -> u8 {
        (self.0 >> 1) & 1
    }

    fn opposite(self) -> Self {
        SideIndex(self.0 ^ 2)
    }

    pub const VALUES: [Self; 4] = [SideIndex(0), SideIndex(1), SideIndex(2), SideIndex(3)];
}

impl std::fmt::Debug for SideIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self.0, f)
    }
}

#[derive(Default)]
struct SideIndexMap<T>([T; 4]);

impl<T> std::ops::Index<SideIndex> for SideIndexMap<T> {
    type Output = T;

    fn index(&self, index: SideIndex) -> &Self::Output {
        &self.0[index.0 as usize]
    }
}

impl<T> std::ops::IndexMut<SideIndex> for SideIndexMap<T> {
    fn index_mut(&mut self, index: SideIndex) -> &mut Self::Output {
        &mut self.0[index.0 as usize]
    }
}

impl<T: std::fmt::Debug> std::fmt::Debug for SideIndexMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self.0, f)
    }
}

pub struct BltGraph {
    chunks: Vec<BltChunk>,
    root_chunk: BltChunkId,
    layout: BltLayout,
    meshes: FxHashMap<NodeId, Vec<skid_steer::Asset<Mesh>>>,
    shadow_graph: Graph,
    loader: skid_steer::Loader,
}

impl BltGraph {
    pub fn new(loader: skid_steer::Loader) -> Self {
        let mut result = BltGraph {
            chunks: vec![BltChunk::new_central()],
            root_chunk: BltChunkId(0),
            layout: BltLayout::default(),
            meshes: FxHashMap::default(),
            shadow_graph: Graph::new(12),
            loader,
        };
        result.init_chunk_mesh(result.root_chunk);
        result
    }

    pub fn get_meshes(&self, node: NodeId) -> &[skid_steer::Asset<Mesh>] {
        self.meshes.get(&node).map_or_default(|x| x.as_slice())
    }

    pub fn initialize_for_test(&mut self) {
        let mut current = self.root_chunk;
        for _ in 0..3 {
            current = self.ensure_outer(current, QuadIndex(0));
        }
    }

    fn new_chunk(&mut self, chunk: BltChunk) -> BltChunkId {
        let id = BltChunkId(self.chunks.len() as u32);
        self.chunks.push(chunk);
        id
    }

    fn init_chunk_mesh(&mut self, chunk: BltChunkId) {
        let mut geometry = MeshGeometryDefinition {
            vertices: Vec::new(),
            indices: Vec::new(),
        };
        for x in (0..(self.layout.horizontal_size as i32)).step_by(2) {
            for y in (0..(self.layout.horizontal_size as i32)).step_by(2) {
                for z in (0..(self.layout.outer_vertical_size as i32)).step_by(2) {
                    add_voxel(
                        self.chunk(chunk),
                        &self.layout,
                        &self.chunk(chunk).position.local,
                        &mut geometry,
                        na::Vector3::new(x, y, z),
                    );
                }
            }
        }
        self.meshes
            .entry(self.chunk(chunk).position.node)
            .or_default()
            .push(self.loader.load(BltChunkSurface { geometry }));
    }

    fn ensure_outer(&mut self, inner: BltChunkId, index: QuadIndex) -> BltChunkId {
        if let Some(outer) = self.chunk(inner).outer_neighbors[index] {
            return outer;
        }
        let mut position = self.chunk(inner).position;
        // TODO: Need to change position's node, as well as expanding the graph. Also need additional parents for numerical stability
        position.local *= self.chunk(inner).outer_isometry(&self.layout, index);
        let outer = self.new_chunk(self.chunk(inner).new_outer(&self.layout, index, position));
        self.chunk_mut(inner).outer_neighbors[index] = Some(outer);
        self.chunk_mut(outer).inner_neighbor = Some(inner);
        for side_index in SideIndex::VALUES {
            let neighbor = index.neighbor(side_index);
            if neighbor.different_inner_chunk {
                if let Some(side_of_inner) = self.chunk(inner).side_neighbors[side_index]
                    && let Some(side_of_outer) =
                        self.chunk(side_of_inner).outer_neighbors[neighbor.index]
                {
                    self.chunk_mut(outer).side_neighbors[side_index] = Some(side_of_outer);
                    self.chunk_mut(side_of_outer).side_neighbors[neighbor.return_side] =
                        Some(outer);
                }
            } else {
                if let Some(side_of_outer) = self.chunk(inner).outer_neighbors[neighbor.index] {
                    self.chunk_mut(outer).side_neighbors[side_index] = Some(side_of_outer);
                    self.chunk_mut(side_of_outer).side_neighbors[neighbor.return_side] =
                        Some(outer);
                }
            }
        }
        self.init_chunk_mesh(outer);
        outer
    }

    fn ensure_side(&mut self, current: BltChunkId, side_index: SideIndex) -> Option<BltChunkId> {
        if let Some(side) = self.chunk(current).side_neighbors[side_index] {
            return Some(side);
        }
        if self.chunk(current).boost == 0.0 {
            return None;
        }
        let neighbor = self
            .chunk(current)
            .inner_neighbor_index
            .neighbor(side_index);
        let parent = self.chunk(current).inner_neighbor.unwrap();
        let side_parent = if neighbor.different_inner_chunk {
            self.ensure_side(parent, side_index)?
        } else {
            parent
        };
        Some(self.ensure_outer(side_parent, neighbor.index))
    }

    fn chunk(&self, chunk: BltChunkId) -> &BltChunk {
        &self.chunks[chunk.0 as usize]
    }

    fn chunk_mut(&mut self, chunk: BltChunkId) -> &mut BltChunk {
        &mut self.chunks[chunk.0 as usize]
    }
}

#[derive(Debug)]
struct BltLayout {
    horizontal_size: u8,
    central_vertical_size: u8,
    outer_vertical_size: u8,
    central_voxel_width: f32,
    voxel_height: f32,
}

impl Default for BltLayout {
    fn default() -> Self {
        Self {
            horizontal_size: 11,
            central_vertical_size: 11,
            outer_vertical_size: 11,
            central_voxel_width: 0.7 / 11.0,
            voxel_height: logf(2.0) / 11.0,
        }
    }
}

#[derive(Debug)]
struct BltChunk {
    inner_neighbor: Option<BltChunkId>,
    inner_neighbor_index: QuadIndex,
    outer_neighbors: QuadIndexMap<Option<BltChunkId>>,
    side_neighbors: SideIndexMap<Option<BltChunkId>>, // Most significant bit: extreme. Least significant bit: axis
    klein_coords: na::Vector2<f32>,
    voxel_width_factor: f32,
    boost: f32,
    position: Position,
}

impl BltChunk {
    fn new_central() -> Self {
        BltChunk {
            inner_neighbor: None,
            inner_neighbor_index: QuadIndex(0),
            outer_neighbors: QuadIndexMap::default(),
            side_neighbors: SideIndexMap::default(),
            klein_coords: na::Vector2::zeros(),
            voxel_width_factor: 1.0,
            boost: 0.0,
            position: Position::origin(),
        }
    }

    fn point_from_voxel(&self, layout: &BltLayout, voxel_coords: na::Vector3<f32>) -> MPoint<f32> {
        self.point_from_chunk(na::Vector3::new(
            voxel_coords[0] * layout.central_voxel_width * self.voxel_width_factor,
            voxel_coords[1] * layout.central_voxel_width * self.voxel_width_factor,
            voxel_coords[2] * layout.voxel_height,
        ))
    }

    fn point_from_chunk(&self, chunk_coords: na::Vector3<f32>) -> MPoint<f32> {
        pseudo_chunk_to_mvector_boosted(
            chunk_to_pseudo_chunk(self.klein_coords, chunk_coords, self.boost),
            self.boost,
        )
        .to_point_unchecked()
    }

    fn outer_isometry(&self, layout: &BltLayout, index: QuadIndex) -> MIsometry<f32> {
        let scale = layout.central_voxel_width
            * layout.horizontal_size as f32
            * self.voxel_width_factor
            * 0.5;
        let chunk_pos = na::Vector3::new(
            scale * index.x() as f32,
            scale * index.y() as f32,
            layout.voxel_height * layout.outer_vertical_size as f32,
        );
        self.isometry_from_chunk(chunk_pos)
    }

    fn isometry_from_chunk(&self, chunk: na::Vector3<f32>) -> MIsometry<f32> {
        let result = chunk_to_isometry(self.klein_coords, chunk, self.boost);
        println!("sanity check: {:?}", result.inverse() * result);
        println!("result: {:?}", result);
        result
    }

    fn new_outer(&self, layout: &BltLayout, index: QuadIndex, position: Position) -> Self {
        let new_boost = self.boost + layout.voxel_height * layout.outer_vertical_size as f32;
        let scale_factor = coshf(new_boost) / coshf(self.boost); // Make computation numerically table
        let displacement_scale = layout.central_voxel_width
            * layout.horizontal_size as f32
            * self.voxel_width_factor
            * 0.5
            / coshf(self.boost);
        let displacement = na::Vector2::new(
            displacement_scale * index.x() as f32,
            displacement_scale * index.y() as f32,
        );
        println!("{:?}", self.klein_coords + displacement);
        BltChunk {
            inner_neighbor: None,
            inner_neighbor_index: index,
            outer_neighbors: QuadIndexMap::default(),
            side_neighbors: SideIndexMap::default(),
            klein_coords: self.klein_coords + displacement,
            voxel_width_factor: self.voxel_width_factor * scale_factor * 0.5,
            boost: new_boost,
            position,
        }
    }
}

#[cfg(test)]
mod tests {
    use common::math::MIsometry;

    use super::*;

    #[test]
    fn example() {
        let example = na::vector![0.1, 0.2, 2.0];
        let boost = 1.0;
        println!(
            "simple: {:?}",
            MIsometry::translation_along(&(na::Vector3::z() * -boost))
                * pseudo_chunk_to_mvector_simple(na::Vector3::new(
                    example.x / coshf(boost),
                    example.y / coshf(boost),
                    boost + example.z
                ))
        );
        println!(
            "boosted: {:?}",
            pseudo_chunk_to_mvector_boosted(example, boost)
        );

        println!(
            "derivative approx 1: {:?}",
            (pseudo_chunk_to_mvector_boosted(example + na::Vector3::x() * 0.01, boost)
                - pseudo_chunk_to_mvector_boosted(example + na::Vector3::x() * -0.01, boost))
                / 0.02
        );
        println!(
            "derivative approx 2: {:?}",
            (pseudo_chunk_to_mvector_boosted(example + na::Vector3::x() * 0.001, boost)
                - pseudo_chunk_to_mvector_boosted(example + na::Vector3::x() * -0.001, boost))
                / 0.002
        );
        println!(
            "derivative: {:?}",
            pseudo_chunk_to_mvector_boosted_partial_x(example, boost)
        );
        BltChunk::new_central().isometry_from_chunk(example);
    }

    #[test]
    fn test_chunk_to_mvector() {
        let klein_coords = na::Vector2::new(0.2, 0.3);
        let chunk = na::Vector3::new(0.25, 0.35, 0.0);
        let boost = 1.1;
        println!("{:?}", chunk_to_mvector_simple(klein_coords, chunk, boost));
        println!(
            "{:?}",
            pseudo_chunk_to_mvector_boosted(
                chunk_to_pseudo_chunk(klein_coords, chunk, boost),
                boost
            )
        );
    }

    #[test]
    fn test_graph_structure() {
        let mut graph = BltGraph::new(skid_steer::Loader::new());
        let a = graph.ensure_outer(graph.root_chunk, QuadIndex(3));
        let b = graph.ensure_outer(a, QuadIndex(0));
        let c = graph.ensure_side(b, SideIndex(0));
        for i in 0..graph.chunks.len() {
            println!(
                "{}: {{ inner_neighbor: {:?}, inner_neighbor_index: {:?}, outer_neighbors: {:?}, side_neighbors: {:?} }}",
                i,
                graph.chunks[i].inner_neighbor,
                graph.chunks[i].inner_neighbor_index,
                graph.chunks[i].outer_neighbors,
                graph.chunks[i].side_neighbors
            );
        }
    }
}
