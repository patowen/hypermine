use std::collections::VecDeque;

use common::{
    dodeca::Side,
    graph::{Graph, NodeId},
    math::{MDirection, MIsometry, MPoint, MVector, PermuteXYZ, sqr},
    proto::Position,
    traversal,
};
use fxhash::{FxHashMap, FxHashSet};
use libm::{coshf, logf, sinhf, sqrtf, tanhf};

use crate::graphics::{
    Mesh,
    asset_loader::AssetLoadContext,
    meshes::{MeshGeometryDefinition, Vertex},
};

#[cfg(test)]
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

#[cfg(test)]
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

#[cfg(test)]
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

fn pseudo_chunk_to_chunk(
    klein_coords: na::Vector2<f32>,
    pseudo_chunk: na::Vector3<f32>,
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

    let horizontal_coords = conversion.try_inverse().unwrap() * pseudo_chunk.xy().push(1.0);
    na::Vector3::new(
        horizontal_coords.x / horizontal_coords.z,
        horizontal_coords.y / horizontal_coords.z,
        pseudo_chunk.z,
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

fn point_to_pseudo_chunk_boosted(point: MPoint<f32>, boost: f32) -> na::Vector3<f32> {
    /*
    In simple mode:
        sinh(Z) = point.mip([0, 0, 1, 0])
    In full mode:
        sinh(Z+boost) = point.mip([0, 0, cosh(boost), -sinh(boost)])
        cosh(Z)sinh(boost) + sinh(Z)cosh(boost) = point.z * cosh(boost) + point.w * sinh(boost)
        cosh(Z)tanh(boost) + sinh(Z) = point.z + point.w * tanh(boost)
    */
    let tanhf_boost = tanhf(boost);
    let point_mip_up = point.z + point.w * tanhf(boost);
    let z = logf(
        (point_mip_up + sqrtf(sqr(point_mip_up) - sqr(tanhf_boost) + 1.0)) / (tanhf_boost + 1.0),
    );
    let squared_xy_dist_times_factor =
        (sqr(point.x) + sqr(point.y)) / sqr(coshf(z) + sinhf(z) * tanhf(boost));
    let squared_xy_dist = 1.0 / (1.0 / squared_xy_dist_times_factor + 1.0 / sqr(coshf(boost)));
    let restoring_factor = sqrtf(squared_xy_dist / (sqr(point.x) + sqr(point.y)));
    na::Vector3::new(point.x * restoring_factor, point.y * restoring_factor, z)
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
    //let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
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
    //let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
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
                position: transform * chunk.point_from_voxel(layout, point.cast()),
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

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
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

    #[allow(unused)]
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
    chunks: Vec<BltChunkWithPosition>,
    root_chunk: BltChunkId,
    layout: BltLayout,
    meshes: FxHashMap<NodeId, Vec<skid_steer::Asset<Mesh>>>,
    shadow_graph: Graph,
    loader: skid_steer::Loader,
    current_chunk: BltChunkId,
}

impl BltGraph {
    pub fn new(loader: skid_steer::Loader) -> Self {
        let initial_transform = common::dodeca::Side::A.reflection()
            * common::dodeca::Vertex::A.dual_to_node()
            * MIsometry::from_columns_unchecked(
                &[MDirection::y(), MDirection::z(), MDirection::x()],
                MPoint::w(),
            );
        let mut result = BltGraph {
            chunks: Vec::new(),
            root_chunk: BltChunkId(0),
            layout: BltLayout::default(),
            meshes: FxHashMap::default(),
            shadow_graph: Graph::new(12),
            loader,
            current_chunk: BltChunkId(0),
        };
        result.root_chunk = result.new_chunk(
            BltChunk::new_central(),
            Position {
                node: NodeId::ROOT,
                local: initial_transform,
            },
        );
        result.current_chunk = result.root_chunk;
        result.init_chunk_mesh(result.root_chunk);
        result
    }

    pub fn get_meshes(&self, node: NodeId) -> &[skid_steer::Asset<Mesh>] {
        self.meshes.get(&node).map_or_default(|x| x.as_slice())
    }

    pub fn fill_radius(&mut self, external_graph: &Graph, mut position: Position, radius: f32) {
        // Note: For simplicity, we fill up to one chunk past the radius, since the logic can be that
        // if we're still within the radius, we continue to expand.
        let mut best_chunk = self.current_chunk;
        let mut best_chunk_cosh_distance = f32::INFINITY;

        while !self.shadow_graph.contains(position.node) {
            let side = external_graph
                .primary_parent_side(position.node)
                .expect("not root");
            position.node = external_graph
                .neighbor(position.node, side)
                .expect("parent");
            position.local = side.reflection() * position.local;
        }

        traversal::ensure_nearby(&mut self.shadow_graph, &position, radius);
        let valid_shadow_nodes: FxHashMap<NodeId, MIsometry<f32>> =
            traversal::nearby_nodes(&self.shadow_graph, &position, radius)
                .into_iter()
                .collect();

        let mut pending = VecDeque::<BltChunkId>::new();
        let mut visited: FxHashSet<BltChunkId> = FxHashSet::default();

        pending.push_back(self.current_chunk);
        visited.insert(self.current_chunk);

        while let Some(blt_chunk_id) = pending.pop_front() {
            let Some(transform) = valid_shadow_nodes.get(&self.chunk_position(blt_chunk_id).node)
            else {
                continue;
            };
            if transform.m44 < best_chunk_cosh_distance {
                best_chunk_cosh_distance = transform.m44;
                best_chunk = blt_chunk_id;
            }
            for i in QuadIndex::VALUES {
                let neighbor = self.ensure_outer(blt_chunk_id, i);
                if visited.insert(neighbor) {
                    pending.push_back(neighbor);
                }
            }
            for i in SideIndex::VALUES {
                let Some(neighbor) = self.ensure_side(blt_chunk_id, i) else {
                    continue;
                };
                if visited.insert(neighbor) {
                    pending.push_back(neighbor);
                }
            }
        }

        self.current_chunk = best_chunk
    }

    pub fn ensure_position(&mut self, mut position: Position, external_graph: &Graph) {
        /*println!(
            "initial: {:?}: {:?}",
            external_graph.debug_node_path(position.node),
            position.local * MPoint::origin()
        );*/
        while !self.shadow_graph.contains(position.node) {
            let side = external_graph
                .primary_parent_side(position.node)
                .expect("not root");
            position.node = external_graph
                .neighbor(position.node, side)
                .expect("parent");
            position.local = side.reflection() * position.local;
        }
        /*println!(
            "shadow: {:?}: {:?}",
            self.shadow_graph.debug_node_path(position.node),
            position.local * MPoint::origin()
        );*/
        let radius = 4.0;
        for (node, transform) in traversal::nearby_nodes(&self.shadow_graph, &position, radius) {
            if node == self.chunk_position(self.current_chunk).node {
                position.node = node;
                position.local = transform.inverse() * position.local;
            }
        }
        if position.node != self.chunk_position(self.current_chunk).node {
            tracing::warn!("Could not find where position is relative to current_chunk");
            return;
        }
        /*println!(
            "matching: {:?}: {:?}",
            self.shadow_graph.debug_node_path(position.node),
            position.local * MPoint::origin()
        );*/
        let mut point = self.chunk_position(self.current_chunk).local.inverse()
            * position.local
            * MPoint::origin();
        for i in 0..10 {
            let voxel = self
                .chunk(self.current_chunk)
                .voxel_from_point(&self.layout, point);
            //println!("{:?} -> {:?}", point, voxel);
            if voxel.x < -0.5
                && let Some(new_chunk) = self.ensure_side(self.current_chunk, SideIndex(0))
            {
                point = self
                    .chunk(self.current_chunk)
                    .side_isometry(&self.layout, SideIndex(0))
                    .inverse()
                    * point;
                self.current_chunk = new_chunk;
            } else if voxel.y < -0.5
                && let Some(new_chunk) = self.ensure_side(self.current_chunk, SideIndex(1))
            {
                point = self
                    .chunk(self.current_chunk)
                    .side_isometry(&self.layout, SideIndex(1))
                    .inverse()
                    * point;
                self.current_chunk = new_chunk;
            } else if voxel.x > self.layout.horizontal_size as f32 + 0.5
                && let Some(new_chunk) = self.ensure_side(self.current_chunk, SideIndex(2))
            {
                point = self
                    .chunk(self.current_chunk)
                    .side_isometry(&self.layout, SideIndex(2))
                    .inverse()
                    * point;
                self.current_chunk = new_chunk;
            } else if voxel.y > self.layout.horizontal_size as f32 + 0.5
                && let Some(new_chunk) = self.ensure_side(self.current_chunk, SideIndex(3))
            {
                point = self
                    .chunk(self.current_chunk)
                    .side_isometry(&self.layout, SideIndex(3))
                    .inverse()
                    * point;
                self.current_chunk = new_chunk;
            } else if voxel.z > self.layout.outer_vertical_size as f32 + 0.5 {
                let x_beyond = voxel.x > self.layout.horizontal_size as f32 * 0.5;
                let y_beyond = voxel.y > self.layout.horizontal_size as f32 * 0.5;
                let quad_index =
                    QuadIndex((if x_beyond { 1 } else { 0 }) | (if y_beyond { 2 } else { 0 }));
                let new_chunk = self.ensure_outer(self.current_chunk, quad_index);
                point = self
                    .chunk(self.current_chunk)
                    .outer_isometry(&self.layout, quad_index)
                    .inverse()
                    * point;
                self.current_chunk = new_chunk;
            } else if voxel.z < -0.5
                && let Some(new_chunk) = self.chunk(self.current_chunk).inner_neighbor
            {
                point = self
                    .chunk(self.current_chunk)
                    .inner_isometry(&self.layout)
                    .inverse()
                    * point;
                self.current_chunk = new_chunk;
            } else {
                break;
            }
            if i == 9 {
                tracing::warn!("Taking longer than expected to reach position");
            }
        }
    }

    fn new_chunk(&mut self, chunk: BltChunk, input_position: Position) -> BltChunkId {
        let id = BltChunkId(self.chunks.len() as u32);
        self.chunks.push(BltChunkWithPosition::from_chunk(
            chunk,
            &self.layout,
            input_position,
            &mut self.shadow_graph,
        ));
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
                        &self.chunk_position(chunk).local,
                        &mut geometry,
                        na::Vector3::new(x, y, z),
                    );
                }
            }
        }
        self.meshes
            .entry(self.chunk_position(chunk).node)
            .or_default()
            .push(self.loader.load(BltChunkSurface { geometry }));
    }

    fn ensure_outer(&mut self, inner: BltChunkId, index: QuadIndex) -> BltChunkId {
        if let Some(outer) = self.chunk(inner).outer_neighbors[index] {
            return outer;
        }
        let mut position = *self.chunk_position(inner);
        // TODO: Need additional parents for numerical stability
        position.local *= self.chunk(inner).outer_isometry(&self.layout, index);
        position.local = position.local.renormalized();
        let outer = self.new_chunk(self.chunk(inner).new_outer(&self.layout, index), position);
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

    #[allow(unused)]
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
        &self.chunks[chunk.0 as usize].chunk
    }

    fn chunk_position(&self, chunk: BltChunkId) -> &Position {
        &self.chunks[chunk.0 as usize].position
    }

    fn chunk_mut(&mut self, chunk: BltChunkId) -> &mut BltChunk {
        &mut self.chunks[chunk.0 as usize].chunk
    }
}

#[derive(Debug)]
struct BltLayout {
    horizontal_size: u8,
    #[allow(unused)]
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
struct BltChunkWithPosition {
    chunk: BltChunk,
    position: Position,
}

impl BltChunkWithPosition {
    pub fn from_chunk(
        chunk: BltChunk,
        layout: &BltLayout,
        input_position: Position,
        shadow_graph: &mut Graph,
    ) -> Self {
        let mut position = input_position;
        let center_point = chunk.center_point(layout);
        // TODO: Need additional parents for numerical stability
        'outer: loop {
            for side in Side::iter() {
                if side.is_facing(&(position.local * center_point)) {
                    position.local = side.reflection() * position.local;
                    position.node = shadow_graph.ensure_neighbor(position.node, side);
                    continue 'outer;
                }
            }
            break;
        }
        BltChunkWithPosition { chunk, position }
    }
}

#[derive(Debug)]
struct BltChunk {
    inner_neighbor: Option<BltChunkId>,
    #[allow(unused)]
    inner_neighbor_index: QuadIndex,
    outer_neighbors: QuadIndexMap<Option<BltChunkId>>,
    side_neighbors: SideIndexMap<Option<BltChunkId>>, // Most significant bit: extreme. Least significant bit: axis
    klein_coords: na::Vector2<f32>,
    voxel_width_factor: f32,
    boost: f32,
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
        }
    }

    fn point_from_voxel(&self, layout: &BltLayout, voxel_coords: na::Vector3<f32>) -> MPoint<f32> {
        self.point_from_chunk(na::Vector3::new(
            voxel_coords[0] * layout.central_voxel_width * self.voxel_width_factor,
            voxel_coords[1] * layout.central_voxel_width * self.voxel_width_factor,
            voxel_coords[2] * layout.voxel_height,
        ))
    }

    fn voxel_from_point(&self, layout: &BltLayout, point: MPoint<f32>) -> na::Vector3<f32> {
        let chunk = self.chunk_from_point(layout, point);
        na::Vector3::new(
            chunk[0] / (layout.central_voxel_width * self.voxel_width_factor),
            chunk[1] / (layout.central_voxel_width * self.voxel_width_factor),
            chunk[2] / layout.voxel_height,
        )
    }

    fn center_point(&self, layout: &BltLayout) -> MPoint<f32> {
        self.point_from_chunk(na::Vector3::new(
            layout.horizontal_size as f32
                * 0.5
                * layout.central_voxel_width
                * self.voxel_width_factor,
            layout.horizontal_size as f32
                * 0.5
                * layout.central_voxel_width
                * self.voxel_width_factor,
            layout.outer_vertical_size as f32 * 0.5 * layout.voxel_height,
        ))
    }

    fn point_from_chunk(&self, chunk_coords: na::Vector3<f32>) -> MPoint<f32> {
        pseudo_chunk_to_mvector_boosted(
            chunk_to_pseudo_chunk(self.klein_coords, chunk_coords, self.boost),
            self.boost,
        )
        .to_point_unchecked()
    }

    fn chunk_from_point(&self, layout: &BltLayout, point: MPoint<f32>) -> na::Vector3<f32> {
        pseudo_chunk_to_chunk(
            self.klein_coords,
            point_to_pseudo_chunk_boosted(point, self.boost),
            self.boost,
        )
    }

    fn inner_isometry(&self, layout: &BltLayout) -> MIsometry<f32> {
        let scale =
            -layout.central_voxel_width * layout.horizontal_size as f32 * self.voxel_width_factor;
        let index = self.inner_neighbor_index;
        let chunk_pos = na::Vector3::new(
            scale * index.x() as f32,
            scale * index.y() as f32,
            -layout.voxel_height * layout.outer_vertical_size as f32,
        );
        self.isometry_from_chunk(chunk_pos)
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

    fn side_isometry(&self, layout: &BltLayout, index: SideIndex) -> MIsometry<f32> {
        let chunk_coord = (if index.extreme() == 0 { -1.0 } else { 1.0 })
            * layout.central_voxel_width
            * layout.horizontal_size as f32
            * self.voxel_width_factor;
        self.isometry_from_chunk(
            (if index.coordinate() == 0 {
                na::Vector3::x()
            } else {
                na::Vector3::y()
            }) * chunk_coord,
        )
    }

    fn isometry_from_chunk(&self, chunk: na::Vector3<f32>) -> MIsometry<f32> {
        let result = chunk_to_isometry(self.klein_coords, chunk, self.boost);
        // println!("sanity check: {:?}", result.inverse() * result);
        // println!("result: {:?}", result);
        result
    }

    fn new_outer(&self, layout: &BltLayout, index: QuadIndex) -> Self {
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
        // println!("{:?}", self.klein_coords + displacement);
        BltChunk {
            inner_neighbor: None,
            inner_neighbor_index: index,
            outer_neighbors: QuadIndexMap::default(),
            side_neighbors: SideIndexMap::default(),
            klein_coords: self.klein_coords + displacement,
            voxel_width_factor: self.voxel_width_factor * scale_factor * 0.5,
            boost: new_boost,
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
        let boosted_result = pseudo_chunk_to_mvector_boosted(example, boost);
        println!("boosted: {:?}", boosted_result);
        println!(
            "Getting back original coordinates: {:?}",
            point_to_pseudo_chunk_boosted(boosted_result.to_point_unchecked(), boost)
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
        println!(
            "{:?}",
            pseudo_chunk_to_chunk(
                klein_coords,
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
        let _c = graph.ensure_side(b, SideIndex(0));
        for i in 0..(graph.chunks.len() as u32) {
            println!(
                "{}: {{ inner_neighbor: {:?}, inner_neighbor_index: {:?}, outer_neighbors: {:?}, side_neighbors: {:?} }}",
                i,
                graph.chunk(BltChunkId(i)).inner_neighbor,
                graph.chunk(BltChunkId(i)).inner_neighbor_index,
                graph.chunk(BltChunkId(i)).outer_neighbors,
                graph.chunk(BltChunkId(i)).side_neighbors
            );
        }
    }
}
