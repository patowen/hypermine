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
    start_chunk: ChunkId,
    start_node_transform: na::Matrix4<f32>,
    ray: Ray,
    tanh_distance: f32,
) -> Result<RayTracingResult, RayTracingError> {
    // A collision check is assumed to be a miss until a collision is found.
    // This `result` variable gets updated over time before being returned.
    let mut result = RayTracingResult {
        tanh_distance,
        intersection: None,
    };

    // Start a breadth-first search of the graph's chunks, performing collision checks in each relevant chunk.
    // The `chunk_queue` contains ordered pairs containing the `ChunkId` and the transformation needed to switch
    // from the original node coordinates to the current chunk's node coordinates.
    let mut visited_chunks: HashSet<ChunkId> = HashSet::new();
    let mut chunk_queue: VecDeque<(ChunkId, na::Matrix4<f32>)> = VecDeque::new();
    chunk_queue.push_back((start_chunk, start_node_transform));

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
                return Err(RayTracingError::OutOfBounds);
            };
        let local_ray = chunk.vertex.node_to_dual().cast::<f32>() * node_transform * &ray;

        let bounding_box = CubicVoxelRegion::from_ray_segment_and_radius(
            dimension,
            &local_ray,
            result.tanh_distance,
            chunk_ray_tracer.max_radius(),
        );

        // Check collision within a single chunk
        if let Some(bounding_box) = bounding_box {
            chunk_ray_tracer.trace_ray(
                &RtChunkContext {
                    dimension,
                    dimension_f32: dimension as f32,
                    chunk,
                    transform: math::mtranspose(&node_transform)
                        * chunk.vertex.dual_to_node().cast(),
                    voxel_data,
                    ray: &local_ray,
                    bounding_box,
                },
                &mut result,
            );
        }

        // Compute the Klein-Beltrami coordinates of the ray segment's endpoints. To check whether neighboring chunks
        // are needed, we need to check whether the endpoints of the line segments lie outside the boundaries of the square
        // bounded by `klein_lower_boundary` and `klein_upper_boundary`.
        let klein_ray_start = na::Point3::from_homogeneous(local_ray.position).unwrap();
        let klein_ray_end =
            na::Point3::from_homogeneous(local_ray.ray_point(result.tanh_distance)).unwrap();

        // Add neighboring chunks as necessary, using one coordinate at a time.
        for axis in 0..3 {
            // Check for neighboring nodes
            if klein_ray_start[axis] <= klein_lower_boundary
                || klein_ray_end[axis] <= klein_lower_boundary
            {
                let side = chunk.vertex.canonical_sides()[axis];
                let next_node_transform = side.reflection().cast::<f32>() * node_transform;
                // Crude check to ensure that the neighboring chunk's node can be in the path of the ray. For simplicity, this
                // check treats each node as a sphere and assumes the ray is pointed directly towards its center. The check is
                // needed because chunk generation uses this approximation, and this check is not guaranteed to pass near corners.
                let ray_node_distance = (next_node_transform * ray.position)
                    .xyz()
                    .magnitude()
                    .acosh();
                let ray_length = result.tanh_distance.atanh();
                if ray_node_distance - ray_length
                    > dodeca::BOUNDING_SPHERE_RADIUS as f32 + chunk_ray_tracer.max_radius()
                {
                    // Ray cannot intersect node
                    continue;
                }
                // If we have to do collision checking on nodes that don't exist in the graph, we cannot have a conclusive result.
                let Some(neighbor) = graph.neighbor(chunk.node, side) else {
                    // Collision checking on nonexistent node
                    return Err(RayTracingError::OutOfBounds);
                };
                // Assuming everything goes well, add the new chunk to the queue.
                let next_chunk = (neighbor, chunk.vertex).into();
                if visited_chunks.insert(next_chunk) {
                    chunk_queue.push_back((next_chunk, next_node_transform));
                }
            }

            // Check for neighboring chunks within the same node
            if klein_ray_start[axis] >= klein_upper_boundary
                || klein_ray_end[axis] >= klein_upper_boundary
            {
                let vertex = chunk.vertex.adjacent_vertices()[axis];
                let next_chunk = (chunk.node, vertex).into();
                if visited_chunks.insert(next_chunk) {
                    chunk_queue.push_back((next_chunk, node_transform));
                }
            }
        }
    }

    Ok(result)
}

pub trait ChunkRayTracer {
    fn trace_ray(&self, ctx: &RtChunkContext, result: &mut RayTracingResult);
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

pub struct RayTracingResult {
    pub tanh_distance: f32,
    pub intersection: Option<RayTracingIntersection>,
}

#[derive(Debug)]
pub enum RayTracingError {
    OutOfBounds,
}

#[derive(Debug)]
pub struct RayTracingIntersection {
    pub chunk: ChunkId,
    pub normal: na::Vector4<f32>,
}

impl RayTracingResult {
    pub fn update(
        &mut self,
        context: &RtChunkContext<'_>,
        tanh_distance: f32,
        normal: na::Vector4<f32>,
    ) {
        self.tanh_distance = tanh_distance;
        self.intersection = Some(RayTracingIntersection {
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

/// Represents a discretized region in the voxel grid contained by an axis-aligned bounding box.
pub struct CubicVoxelRegion {
    // The bounds are of the form [[x_min, x_max], [y_min, y_max], [z_min, z_max]], using voxel coordinates with margins.
    // Any voxel that intersects the cube of interest is included in these bounds. By adding or subtracting 1 in the right
    // places, these bounds can be used to find other useful info related to the cube of interset, such as what grid points
    // it contains.
    bounds: [[usize; 2]; 3],
}

impl CubicVoxelRegion {
    /// Returns a bounding box that is guaranteed to cover a given radius around a ray segment. Returns None if the
    /// bounding box lies entirely outside the chunk, including its margins.
    pub fn from_ray_segment_and_radius(
        dimension: usize,
        ray: &Ray,
        tanh_distance: f32,
        radius: f32,
    ) -> Option<CubicVoxelRegion> {
        let dimension_f32 = dimension as f32;
        let grid_start = na::Point3::from_homogeneous(ray.position).unwrap()
            * Vertex::dual_to_chunk_factor() as f32
            * dimension_f32;
        let grid_end = na::Point3::from_homogeneous(ray.ray_point(tanh_distance)).unwrap()
            * Vertex::dual_to_chunk_factor() as f32
            * dimension_f32;
        let max_grid_radius = radius * Vertex::dual_to_chunk_factor() as f32 * dimension_f32;
        let mut bounds = [[0; 2]; 3];
        for axis in 0..3 {
            let grid_min = grid_start[axis].min(grid_end[axis]) - max_grid_radius;
            let grid_max = grid_start[axis].max(grid_end[axis]) + max_grid_radius;
            let voxel_min = (grid_min + 1.0).floor().max(0.0);
            let voxel_max = (grid_max + 1.0).floor().min(dimension_f32 + 1.0);

            // This will happen when voxel_min is greater than dimension+1 or voxel_max is less than 0, which
            // occurs when the cube is out of range.
            if voxel_min > voxel_max {
                return None;
            }

            // We convert to usize here instead of earlier because out-of-range voxel coordinates can be negative.
            bounds[axis] = [voxel_min.floor() as usize, voxel_max.floor() as usize];
        }

        Some(CubicVoxelRegion { bounds })
    }

    /// Creates an iterator over grid points contained in the region, represented as ordered triples
    pub fn grid_point_iterator(
        &self,
        axis0: usize,
        axis1: usize,
        axis2: usize,
    ) -> impl Iterator<Item = (usize, usize, usize)> {
        let bounds = self.bounds;
        (bounds[axis0][0]..bounds[axis0][1]).flat_map(move |i| {
            (bounds[axis1][0]..bounds[axis1][1])
                .flat_map(move |j| (bounds[axis2][0]..bounds[axis2][1]).map(move |k| (i, j, k)))
        })
    }

    /// Creates an iterator over grid lines intersecting the region, represented as ordered pairs determining the line's two fixed coordinates
    pub fn grid_line_iterator(
        &self,
        axis0: usize,
        axis1: usize,
    ) -> impl Iterator<Item = (usize, usize)> {
        let bounds = self.bounds;
        (bounds[axis0][0]..bounds[axis0][1])
            .flat_map(move |i| (bounds[axis1][0]..bounds[axis1][1]).map(move |j| (i, j)))
    }

    /// Creates an iterator over grid planes intersecting the region, represented as integers determining the plane's fixed coordinate
    pub fn grid_plane_iterator(&self, axis: usize) -> impl Iterator<Item = usize> {
        self.bounds[axis][0]..self.bounds[axis][1]
    }
}
