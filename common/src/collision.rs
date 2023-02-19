use std::collections::{HashSet, VecDeque};

use crate::{
    dodeca::{self, Vertex},
    math,
    node::{Chunk, ChunkId, DualGraph, VoxelData},
    sphere_collider::chunk_sphere_cast,
    world::Material,
};

/// Performs sphere casting (swept collision query) against the voxels in the `DualGraph`
///
/// The `start_node_transform` parameter determines which coordinate system the `ray` parameter and any resulting hit
/// normals are given in. Specifically, the `start_node_transform` matrix converts this coordinate system to the coordinate
/// system of `start_chunk`'s node.
///
/// The `tanh_distance` is the hyperbolic tangent of the distance along the ray to check for hits.
pub fn sphere_cast(
    graph: &DualGraph,
    dimension: usize,
    collider_radius: f32,
    start_chunk: ChunkId,
    start_node_transform: na::Matrix4<f32>,
    ray: &Ray,
    tanh_distance: f32,
) -> Result<CastEndpoint, SphereCastError> {
    // A collision check is assumed to be a miss until a collision is found.
    // This `endpoint` variable gets updated over time before being returned.
    let mut endpoint = CastEndpoint {
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
    let klein_lower_boundary = collider_radius.tanh();
    let klein_upper_boundary =
        ((Vertex::chunk_to_dual_factor() as f32).atanh() - collider_radius).tanh();

    // Breadth-first search loop
    while let Some((chunk, node_transform)) = chunk_queue.pop_front() {
        let node = graph.get(chunk.node).as_ref().unwrap();
        let Chunk::Populated {
                voxels: ref voxel_data,
                ..
            } = node.chunks[chunk.vertex] else {
                // Collision checking on unpopulated chunk
                return Err(SphereCastError::OutOfBounds);
            };
        let local_ray = chunk.vertex.node_to_dual().cast::<f32>() * node_transform * ray;

        let dual_to_grid_factor = Vertex::dual_to_chunk_factor() as f32 * dimension as f32;

        let bounding_box = VoxelAABB::from_ray_segment_and_radius(
            dimension,
            dual_to_grid_factor,
            &local_ray,
            endpoint.tanh_distance,
            collider_radius,
        );

        // Check collision within a single chunk
        if let Some(bounding_box) = bounding_box {
            chunk_sphere_cast(
                &ChunkSphereCastContext {
                    dimension,
                    dual_to_grid_factor,
                    chunk,
                    transform: math::mtranspose(&node_transform)
                        * chunk.vertex.dual_to_node().cast(),
                    voxel_data,
                    collider_radius,
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
                    > dodeca::BOUNDING_SPHERE_RADIUS as f32 + collider_radius
                {
                    // Ray cannot intersect node
                    continue;
                }
                // If we have to do collision checking on nodes that don't exist in the graph, we cannot have a conclusive result.
                let Some(neighbor) = graph.neighbor(chunk.node, side) else {
                    // Collision checking on nonexistent node
                    return Err(SphereCastError::OutOfBounds);
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

/// Contains all the immutable data needed for `chunk_sphere_cast` to perform its logic
pub struct ChunkSphereCastContext<'a> {
    pub dimension: usize,
    pub dual_to_grid_factor: f32,
    pub chunk: ChunkId,
    pub transform: na::Matrix4<f32>,
    pub voxel_data: &'a VoxelData,
    pub collider_radius: f32,
    pub ray: &'a Ray,
    pub bounding_box: VoxelAABB,
}

impl ChunkSphereCastContext<'_> {
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
#[derive(Debug)]
pub struct CastEndpoint {
    /// The tanh of the length of the resulting ray segment so far. As new intersections are found, the
    /// ray segment gets shorter each time.
    pub tanh_distance: f32,

    /// Information about the intersection at the end of the ray segment. If this is `None`, there
    /// are no intersections.
    pub hit: Option<CastHit>,
}

#[derive(Debug)]
pub enum SphereCastError {
    OutOfBounds,
}

/// Information about the intersection at the end of a ray segment.
#[derive(Debug)]
pub struct CastHit {
    /// Which chunk in the graph the hit occurred
    pub chunk: ChunkId,

    /// Represents the normal vector of the hit surface in the original coordinate system
    /// of the sphere casting. To get the actual normal vector, project it so that it is orthogonal
    /// to the endpoint in Lorentz space.
    pub normal: na::Vector4<f32>,
}

impl CastEndpoint {
    /// Convenience function to report a new hit found when sphere casting. The `normal` parameter
    /// should be provided in the chunk's "dual" coordinate system.
    pub fn update(
        &mut self,
        context: &ChunkSphereCastContext<'_>,
        tanh_distance: f32,
        normal: na::Vector4<f32>,
    ) {
        self.tanh_distance = tanh_distance;
        self.hit = Some(CastHit {
            chunk: context.chunk,
            normal: context.transform * normal,
        });
    }
}

/// A ray in hyperbolic space. The fields must be lorentz normalized, with `mip(position, position) == -1`,
/// `mip(direction, direction) == 1`, and `mip(position, direction) == 0`.
#[derive(Debug)]
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

    /// Iterator over grid points contained in the region, represented as ordered triples
    pub fn grid_points(
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

    /// Iterator over grid lines intersecting the region, represented as ordered pairs determining the line's two fixed coordinates
    pub fn grid_lines(&self, axis0: usize, axis1: usize) -> impl Iterator<Item = (usize, usize)> {
        let bounds = self.bounds;
        (bounds[axis0][0]..bounds[axis0][1])
            .flat_map(move |i| (bounds[axis1][0]..bounds[axis1][1]).map(move |j| (i, j)))
    }

    /// Iterator over grid planes intersecting the region, represented as integers determining the plane's fixed coordinate
    pub fn grid_planes(&self, axis: usize) -> impl Iterator<Item = usize> {
        self.bounds[axis][0]..self.bounds[axis][1]
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        dodeca::{Side, Vertex},
        graph::NodeId,
        node::populate_fresh_nodes,
        proto::Position,
        traversal::{ensure_nearby, nearby_nodes},
    };

    use super::*;

    #[test]
    fn sphere_cast_example() {
        let dimension: usize = 12;
        let dual_to_grid_factor = Vertex::dual_to_chunk_factor() as f32 * dimension as f32;
        let mut graph = DualGraph::new();
        let graph_radius = 3.0;

        // Set up a graph with void chunks
        ensure_nearby(&mut graph, &Position::origin(), graph_radius);
        populate_fresh_nodes(&mut graph);
        for (node, _) in nearby_nodes(&graph, &Position::origin(), graph_radius) {
            for vertex in dodeca::Vertex::iter() {
                graph[ChunkId::new(node, vertex)] = Chunk::Populated {
                    voxels: VoxelData::Solid(Material::Void),
                    surface: None,
                };
            }
        }

        // Populate an arbitrary nearby chunk
        let chosen_chunk = ChunkId::new(graph.neighbor(NodeId::ROOT, Side::G).unwrap(), Vertex::I);
        let Chunk::Populated { voxels, .. } = graph.get_chunk_mut(chosen_chunk).unwrap() else {
            panic!("All chunks should be populated.");
        };

        // Populate (3, 4, 5) with dirt.
        voxels.data_mut(dimension as u8)[3 + 4 * (dimension + 2) + 5 * (dimension + 2).pow(2)] =
            Material::Dirt;

        let dirt_position = Side::G.reflection().cast()
            * (Vertex::I.chunk_to_node().cast()
                * math::lorentz_normalize(&na::Vector4::new(
                    2.5 / dual_to_grid_factor,
                    3.5 / dual_to_grid_factor,
                    4.5 / dual_to_grid_factor,
                    1.0,
                )));

        let ray_position = math::origin();
        let ray_direction = dirt_position - ray_position;

        let ray = Ray::new(
            ray_position,
            math::lorentz_normalize(
                &(ray_direction + ray_position * math::mip(&ray_position, &ray_direction)),
            ),
        );

        let tanh_distance = (-math::mip(&ray_position, &dirt_position)).acosh().tanh();

        let endpoint = sphere_cast(
            &graph,
            dimension,
            0.02,
            ChunkId::new(NodeId::ROOT, Vertex::A),
            na::Matrix4::identity(),
            &ray,
            tanh_distance,
        )
        .expect("conclusive collision result");

        assert!(endpoint.hit.is_some(), "no collision detected");
        assert_eq!(
            endpoint.hit.as_ref().unwrap().chunk,
            chosen_chunk,
            "collision occurred in wrong chunk"
        );
        assert!(
            math::mip(&endpoint.hit.as_ref().unwrap().normal, &ray.direction) < 0.0,
            "normal is facing the wrong way"
        );
    }

    /// Any voxel AABB should at least cover a capsule-shaped region consisting of all points
    /// `radius` units away from the ray's line segment. This region consists of two spheres
    /// and a cylinder. We only test planes because covered lines and points are a strict subset.
    #[test]
    fn voxel_aabb_coverage() {
        let dimension = 12;
        let dual_to_grid_factor = Vertex::dual_to_chunk_factor() as f32 * dimension as f32;

        // Pick an arbitrary ray by transforming the positive-x-axis ray.
        let ray = na::Rotation3::from_euler_angles(0.1, 0.2, 0.3).to_homogeneous()
            * math::translate_along(&na::Vector3::new(0.2, 0.3, 0.1))
            * &Ray::new(na::Vector4::w(), na::Vector4::x());

        let tanh_distance = 0.2;
        let radius = 0.1;

        let aabb = VoxelAABB::from_ray_segment_and_radius(
            dimension,
            dual_to_grid_factor,
            &ray,
            tanh_distance,
            radius,
        )
        .unwrap();

        // Test planes in all 3 axes. For variable names and further comments, we use a tuv coordinate system,
        // which is a permuted xyz coordinate system.
        for t_axis in 0..3 {
            let covered_planes: HashSet<_> = aabb.grid_planes(t_axis).collect();

            // Check that all uv-aligned planes that should be covered are covered
            let ray_end = math::lorentz_normalize(&ray.ray_point(tanh_distance));
            for t in 0..=dimension {
                if covered_planes.contains(&t) {
                    continue;
                }

                let mut plane_normal = na::Vector4::zeros();
                plane_normal[t_axis] = 1.0;
                plane_normal[3] = t as f32 / dual_to_grid_factor;
                let plane_normal = math::lorentz_normalize(&plane_normal);

                // Get the sinh of the signed distance from the ray segment's endpoints to the plane
                let ray_start_sinh_displacement = math::mip(&ray.position, &plane_normal);
                let ray_end_sinh_displacement = math::mip(&ray_end, &plane_normal);

                // Ensure that both ray endpoints are far enough away on the same side of the plane
                assert!(
                    ray_start_sinh_displacement.min(ray_end_sinh_displacement) > radius.sinh()
                        || ray_start_sinh_displacement.max(ray_end_sinh_displacement)
                            < -radius.sinh(),
                    "Plane not covered: axis={}, t={}",
                    t_axis,
                    t,
                );
            }
        }

        // We do not test that the right lines and points are covered because the current
        // implementation guarantees that if the right planes are covered, so are the right
        // lines and points. Actually testing for this would be too complicated.
    }
}
