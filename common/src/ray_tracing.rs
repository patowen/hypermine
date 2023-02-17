use std::collections::{HashSet, VecDeque};

use crate::{
    dodeca::{self, Vertex},
    math,
    node::{Chunk, ChunkId, DualGraph, VoxelData},
    world::Material,
};

/// Performs ray tracing against the voxels in the `DualGraph`. This function is suitable for collision checking with
/// any different collider shape. Each collider would have a different implemenation of `ChunkRayTracer`, which describes
/// how the collision checking math works within a single chunk.
///
/// The `start_node_transform` parameter determines which coordinate system the `ray` parameter and any resulting hit
/// normals are given in. Specifically, the `start_node_transform` matrix converts this coordinate system to the coordinate
/// system of `start_chunk`'s node.
///
/// The `tanh_distance` is the hyperbolic tangent of the distance along the ray to check for hits.
pub fn trace_ray(
    graph: &DualGraph,
    dimension: usize,
    chunk_ray_tracer: &impl ChunkRayTracer,
    start_chunk: ChunkId,
    start_node_transform: na::Matrix4<f32>,
    ray: &Ray,
    tanh_distance: f32,
) -> Result<RayEndpoint, RayTracingError> {
    // A collision check is assumed to be a miss until a collision is found.
    // This `endpoint` variable gets updated over time before being returned.
    let mut endpoint = RayEndpoint {
        tanh_distance,
        hit: None,
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
        let local_ray = chunk.vertex.node_to_dual().cast::<f32>() * node_transform * ray;

        let dual_to_grid_factor = Vertex::dual_to_chunk_factor() as f32 * dimension as f32;

        let bounding_box = VoxelAABB::from_ray_segment_and_radius(
            dimension,
            dual_to_grid_factor,
            &local_ray,
            endpoint.tanh_distance,
            chunk_ray_tracer.max_radius(),
        );

        // Check collision within a single chunk
        if let Some(bounding_box) = bounding_box {
            chunk_ray_tracer.trace_ray(
                &ChunkRayTracingContext {
                    dimension,
                    dual_to_grid_factor,
                    chunk,
                    transform: math::mtranspose(&node_transform)
                        * chunk.vertex.dual_to_node().cast(),
                    voxel_data,
                    ray: &local_ray,
                    bounding_box,
                },
                &mut endpoint,
            );
        }

        // Compute the Klein-Beltrami coordinates of the ray segment's endpoints. To check whether neighboring chunks
        // are needed, we need to check whether the endpoints of the line segments lie outside the boundaries of the square
        // bounded by `klein_lower_boundary` and `klein_upper_boundary`.
        let klein_ray_start = na::Point3::from_homogeneous(local_ray.position).unwrap();
        let klein_ray_end =
            na::Point3::from_homogeneous(local_ray.ray_point(endpoint.tanh_distance)).unwrap();

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
                let ray_length = endpoint.tanh_distance.atanh();
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
                let next_chunk = ChunkId::new(neighbor, chunk.vertex);
                if visited_chunks.insert(next_chunk) {
                    chunk_queue.push_back((next_chunk, next_node_transform));
                }
            }

            // Check for neighboring chunks within the same node
            if klein_ray_start[axis] >= klein_upper_boundary
                || klein_ray_end[axis] >= klein_upper_boundary
            {
                let vertex = chunk.vertex.adjacent_vertices()[axis];
                let next_chunk = ChunkId::new(chunk.node, vertex);
                if visited_chunks.insert(next_chunk) {
                    chunk_queue.push_back((next_chunk, node_transform));
                }
            }
        }
    }

    Ok(endpoint)
}

/// Wraps any logic needed for a particular collider to perform ray tracing within a given chunk
pub trait ChunkRayTracer {
    /// Performs ray tracing with a single chunk. If an intersection is found, `endpoint` is updated
    /// to reflect this intersection, overriding any old intersections that may have been found.
    fn trace_ray(&self, ctx: &ChunkRayTracingContext, endpoint: &mut RayEndpoint);

    /// Returns the radius of the sphere of influence of the collider. `trace_ray` might not be called on
    /// chunks outside this sphere of influence.
    fn max_radius(&self) -> f32;
}

/// Contains all the immutable data needed for `ChunkRayTracer` to perform its logic
pub struct ChunkRayTracingContext<'a> {
    pub dimension: usize,
    pub dual_to_grid_factor: f32,
    pub chunk: ChunkId,
    pub transform: na::Matrix4<f32>,
    pub voxel_data: &'a VoxelData,
    pub ray: &'a Ray,
    pub bounding_box: VoxelAABB,
}

impl ChunkRayTracingContext<'_> {
    /// Convenience function to get data from a single voxel given its coordinates.
    /// Each coordinate runs from `0` to `dimension + 1` inclusive, as margins are included.
    /// To get data within the chunk, use coordinates in the range from `1` to `dimension` inclusive.
    pub fn get_voxel(&self, coords: [usize; 3]) -> Material {
        let dimension_with_margin = self.dimension + 2;
        debug_assert!(coords[0] < dimension_with_margin);
        debug_assert!(coords[1] < dimension_with_margin);
        debug_assert!(coords[2] < dimension_with_margin);
        self.voxel_data.get(
            coords[0]
                + coords[1] * dimension_with_margin
                + coords[2] * dimension_with_margin.pow(2),
        )
    }
}

/// Where a cast ray ended, and all information about the hit at the end if such a hit occurred.
pub struct RayEndpoint {
    /// The tanh of the length of the resulting ray segment so far. As new intersections are found, the
    /// ray segment gets shorter each time.
    pub tanh_distance: f32,

    /// Information about the intersection at the end of the ray segment. If this is `None`, there
    /// are no intersections.
    pub hit: Option<RayHit>,
}

#[derive(Debug)]
pub enum RayTracingError {
    OutOfBounds,
}

/// Information about the intersection at the end of a ray segment.
pub struct RayHit {
    /// Which chunk in the graph the hit occurred
    pub chunk: ChunkId,

    /// Represents the normal vector of the hit surface in the original coordinate system
    /// of the ray tracing. To get the actual normal vector, project it so that it is orthogonal
    /// to the ray's endpoint in Lorentz space.
    pub normal: na::Vector4<f32>,
}

impl RayEndpoint {
    /// Convenience function to report a new hit found when ray tracing. The `normal` parameter
    /// should be provided in the chunk's "dual" coordinate system.
    pub fn update(
        &mut self,
        context: &ChunkRayTracingContext<'_>,
        tanh_distance: f32,
        normal: na::Vector4<f32>,
    ) {
        self.tanh_distance = tanh_distance;
        self.hit = Some(RayHit {
            chunk: context.chunk,
            normal: context.transform * normal,
        });
    }
}

/// A ray in hyperbolic space. The fields must be lorentz normalized, with `mip(position, position) == -1`,
/// `mip(direction, direction) == 1`, and `mip(position, direction) == 0`.
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

    /// Returns a point along this ray `atanh(tanh_distance)` units away from the origin. This point
    /// is _not_ lorentz normalized.
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
pub struct VoxelAABB {
    // The bounds are of the form [[x_min, x_max], [y_min, y_max], [z_min, z_max]], using voxel coordinates with margins.
    // Any voxel that intersects the cube of interest is included in these bounds. By adding or subtracting 1 in the right
    // places, these bounds can be used to find other useful info related to the cube of interset, such as what grid points
    // it contains.
    bounds: [[usize; 2]; 3],
}

impl VoxelAABB {
    /// Returns a bounding box that is guaranteed to cover a given radius around a ray segment. Returns None if the
    /// bounding box lies entirely outside the chunk, including its margins.
    pub fn from_ray_segment_and_radius(
        dimension: usize,
        dual_to_grid_factor: f32,
        ray: &Ray,
        tanh_distance: f32,
        radius: f32,
    ) -> Option<VoxelAABB> {
        // Convert the ray to grid coordinates
        let grid_start = na::Point3::from_homogeneous(ray.position).unwrap() * dual_to_grid_factor;
        let grid_end = na::Point3::from_homogeneous(ray.ray_point(tanh_distance)).unwrap()
            * dual_to_grid_factor;
        // Convert the radius to grid coordinates using a crude conservative estimate
        let max_grid_radius = radius * dual_to_grid_factor;
        let mut bounds = [[0; 2]; 3];
        for axis in 0..3 {
            let grid_min = grid_start[axis].min(grid_end[axis]) - max_grid_radius;
            let grid_max = grid_start[axis].max(grid_end[axis]) + max_grid_radius;
            let voxel_min = (grid_min + 1.0).floor().max(0.0);
            let voxel_max = (grid_max + 1.0).floor().min(dimension as f32 + 1.0);

            // This will happen when voxel_min is greater than dimension+1 or voxel_max is less than 0, which
            // occurs when the cube is out of range.
            if voxel_min > voxel_max {
                return None;
            }

            // We convert to usize here instead of earlier because out-of-range voxel coordinates can be negative.
            bounds[axis] = [voxel_min.floor() as usize, voxel_max.floor() as usize];
        }

        Some(VoxelAABB { bounds })
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
