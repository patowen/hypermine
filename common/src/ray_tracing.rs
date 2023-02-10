use std::collections::{HashSet, VecDeque};

use crate::{
    dodeca::{self, Vertex},
    math,
    node::{Chunk, ChunkId, DualGraph, VoxelData},
    world::Material,
};

pub fn trace_ray(
    graph: &DualGraph,
    dimension: usize,
    chunk_ray_tracer: &impl ChunkRayTracer,
    chunk: ChunkId,
    transform: na::Matrix4<f32>, // TODO: Consider where this transformation gets applied
    ray: Ray,
    tanh_distance: f32,
) -> RayStatus {
    // A collision check is assumed to be a miss until a collision is found.
    // This `status` variable gets updated over time before being returned.
    let mut status = RayStatus {
        result: RayTracingResult::Miss,
        tanh_distance,
    };

    // Start a breadth-first search of the graph's chunks, performing collision checks in each relevant chunk.
    // The `chunk_queue` contains ordered pairs containing the `ChunkId` and the transformation needed to switch
    // from the original node coordinates to the current chunk's node coordinates.
    let mut visited_chunks: HashSet<ChunkId> = HashSet::new();
    let mut chunk_queue: VecDeque<(ChunkId, na::Matrix4<f32>)> = VecDeque::new();
    chunk_queue.push_back((chunk, na::Matrix4::identity()));

    // Precalculate the chunk boundaries for collision purposes. If the collider goes outside these bounds,
    // the corresponding neighboring chunk will also be used for collision checking.
    let klein_lower_boundary = chunk_ray_tracer.max_radius().tanh();
    let klein_upper_boundary =
        ((Vertex::chunk_to_dual_factor() as f32).atanh() - chunk_ray_tracer.max_radius()).tanh();

    // Breadth-first search loop
    while let Some((chunk, node_transform)) = chunk_queue.pop_front() {
        let node = graph.get(chunk.node).as_ref().unwrap();
        let Chunk::Populated {
                voxels: ref voxel_data,
                ..
            } = node.chunks[chunk.vertex] else {
                // Collision checking on unpopulated chunk
                status.result = RayTracingResult::Inconclusive;
                return status;
            };
        let local_ray = chunk.vertex.node_to_dual().cast::<f32>() * node_transform * &ray;

        // Check collision within a single chunk
        chunk_ray_tracer.trace_ray(
            &RtChunkContext {
                dimension,
                dimension_f32: dimension as f32,
                chunk,
                transform: transform
                    * math::mtranspose(&node_transform)
                    * chunk.vertex.dual_to_node().cast(),
                voxel_data,
                ray: &local_ray,
                bounding_box: CubicVoxelRegion::from_ray_segment_and_radius(
                    dimension,
                    &local_ray,
                    status.tanh_distance,
                    chunk_ray_tracer.max_radius(),
                ),
            },
            &mut status,
        );

        // Compute the Klein-Beltrami coordinates of the ray segment's endpoints. To check whether neighboring chunks
        // are needed, we need to check whether the endpoints of the line segments lie outside the boundaries of the square
        // bounded by `klein_lower_boundary` and `klein_upper_boundary`.
        let klein_ray_start = na::Point3::from_homogeneous(local_ray.position).unwrap();
        let klein_ray_end =
            na::Point3::from_homogeneous(local_ray.ray_point(status.tanh_distance)).unwrap();

        // Add neighboring chunks as necessary, using one coordinate at a time.
        for coord in 0..3 {
            // Check for neighboring nodes
            if klein_ray_start[coord] <= klein_lower_boundary
                || klein_ray_end[coord] <= klein_lower_boundary
            {
                let side = chunk.vertex.canonical_sides()[coord];
                let next_node_transform = side.reflection().cast::<f32>() * node_transform;
                // Crude check to ensure that the neighboring chunk's node can be in the path of the ray. For simplicity, this
                // check treats each node as a sphere and assumes the ray is pointed directly towards its center. The check is
                // needed because chunk generation uses this approximation, and this check is not guaranteed to pass near corners.
                let ray_node_distance = (next_node_transform * ray.position)
                    .xyz()
                    .magnitude()
                    .acosh();
                let ray_length = status.tanh_distance.atanh();
                if ray_node_distance - ray_length
                    > dodeca::BOUNDING_SPHERE_RADIUS as f32 + chunk_ray_tracer.max_radius()
                {
                    // Ray cannot intersect node
                    continue;
                }
                // If we have to do collision checking on nodes that don't exist in the graph, we cannot have a conclusive result.
                let Some(neighbor) = graph.neighbor(chunk.node, side) else {
                    // Collision checking on nonexistent node
                    status.result = RayTracingResult::Inconclusive;
                    return status;
                };
                // Assuming everything goes well, add the new chunk to the queue.
                let next_chunk = (neighbor, chunk.vertex).into();
                if visited_chunks.insert(next_chunk) {
                    chunk_queue.push_back((next_chunk, next_node_transform));
                }
            }

            // Check for neighboring chunks within the same node
            if klein_ray_start[coord] >= klein_upper_boundary
                || klein_ray_end[coord] >= klein_upper_boundary
            {
                let vertex = chunk.vertex.adjacent_vertices()[coord];
                let next_chunk = (chunk.node, vertex).into();
                if visited_chunks.insert(next_chunk) {
                    chunk_queue.push_back((next_chunk, node_transform));
                }
            }
        }
    }

    status
}

pub trait ChunkRayTracer {
    fn trace_ray(&self, ctx: &RtChunkContext, status: &mut RayStatus);
    fn max_radius(&self) -> f32;
}

pub struct RtChunkContext<'a> {
    pub dimension: usize,
    pub dimension_f32: f32,
    pub chunk: ChunkId,
    pub transform: na::Matrix4<f32>,
    pub voxel_data: &'a VoxelData,
    pub ray: &'a Ray,
    pub bounding_box: CubicVoxelRegion,
}

impl RtChunkContext<'_> {
    // Also allows access to margins
    pub fn get_voxel(&self, coords: [usize; 3]) -> Material {
        let dimension_with_margin = self.dimension + 2;
        assert!(coords[0] < dimension_with_margin);
        assert!(coords[1] < dimension_with_margin);
        assert!(coords[2] < dimension_with_margin);
        self.voxel_data.get(
            coords[0]
                + coords[1] * dimension_with_margin
                + coords[2] * dimension_with_margin.pow(2),
        )
    }
}

pub struct RayStatus {
    pub tanh_distance: f32,
    pub result: RayTracingResult,
}

#[derive(Debug)]
pub enum RayTracingResult {
    Miss,
    Intersection(RayTracingIntersection),
    Inconclusive,
}

#[derive(Debug)]
pub struct RayTracingIntersection {
    pub chunk: ChunkId,
    pub normal: na::Vector4<f32>,
}

impl RayStatus {
    pub fn update(
        &mut self,
        context: &RtChunkContext<'_>,
        tanh_distance: f32,
        normal: na::Vector4<f32>,
    ) {
        self.tanh_distance = tanh_distance;
        self.result = RayTracingResult::Intersection(RayTracingIntersection {
            chunk: context.chunk,
            normal: context.transform * normal,
        });
    }
}

pub struct Ray {
    pub position: na::Vector4<f32>,
    pub direction: na::Vector4<f32>,
}

impl Ray {
    pub fn new(position: na::Vector4<f32>, direction: na::Vector4<f32>) -> Ray {
        Ray {
            position,
            direction,
        }
    }

    /// Returns a point along this ray tanh_distance units away from the origin
    pub fn ray_point(&self, tanh_distance: f32) -> na::Vector4<f32> {
        self.position + self.direction * tanh_distance
    }
}

impl std::ops::Mul<&Ray> for na::Matrix4<f32> {
    type Output = Ray;

    #[inline]
    fn mul(self, rhs: &Ray) -> Self::Output {
        Ray {
            position: self * rhs.position,
            direction: self * rhs.direction,
        }
    }
}

pub struct CubicVoxelRegion {
    bounds: [[usize; 2]; 3],
}

impl CubicVoxelRegion {
    pub fn from_ray_segment_and_radius(
        dimension: usize,
        ray: &Ray,
        tanh_distance: f32,
        radius: f32,
    ) -> CubicVoxelRegion {
        let float_dimension = dimension as f32;
        let voxel_start = na::Point3::from_homogeneous(ray.position).unwrap()
            * Vertex::dual_to_chunk_factor() as f32
            * float_dimension;
        let voxel_end = na::Point3::from_homogeneous(ray.ray_point(tanh_distance)).unwrap()
            * Vertex::dual_to_chunk_factor() as f32
            * float_dimension;
        let max_voxel_radius = radius * Vertex::dual_to_chunk_factor() as f32 * float_dimension;
        let bounds = [0, 1, 2].map(|coord| {
            if !voxel_start[coord].is_finite() || !voxel_end[coord].is_finite() {
                return [0, dimension + 1];
            }
            let result_min =
                (voxel_start[coord].min(voxel_end[coord]) - max_voxel_radius).max(-0.5) + 1.0;
            let result_max = (voxel_start[coord].max(voxel_end[coord]) + max_voxel_radius)
                .min(dimension as f32 + 0.5)
                + 1.0;

            if result_min > result_max {
                // Empty range
                return [1, 0];
            }

            [result_min.floor() as usize, result_max.floor() as usize]
        });

        CubicVoxelRegion { bounds }
    }

    /// Creates an iterator over voxels, represented as ordered triples
    pub fn grid_point_iterator(
        &self,
        coord0: usize,
        coord1: usize,
        coord2: usize,
    ) -> impl Iterator<Item = (usize, usize, usize)> {
        let bounds = self.bounds;
        (bounds[coord0][0]..bounds[coord0][1]).flat_map(move |i| {
            (bounds[coord1][0]..bounds[coord1][1])
                .flat_map(move |j| (bounds[coord2][0]..bounds[coord2][1]).map(move |k| (i, j, k)))
        })
    }

    /// Creates an iterator over voxel lines, represented as ordered pairs determining the line's two fixed coordinates
    pub fn grid_line_iterator(
        &self,
        coord0: usize,
        coord1: usize,
    ) -> impl Iterator<Item = (usize, usize)> {
        let bounds = self.bounds;
        (bounds[coord0][0]..bounds[coord0][1])
            .flat_map(move |i| (bounds[coord1][0]..bounds[coord1][1]).map(move |j| (i, j)))
    }

    /// Creates an iterator over voxel planes, represented as integers determining the plane's fixed coordinate
    pub fn grid_plane_iterator(&self, coord: usize) -> impl Iterator<Item = usize> {
        self.bounds[coord][0]..self.bounds[coord][1]
    }
}
