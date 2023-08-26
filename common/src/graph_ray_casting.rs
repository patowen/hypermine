use std::collections::VecDeque;

use fxhash::FxHashSet;

use crate::{
    chunk_ray_casting::chunk_ray_cast,
    dodeca::Vertex,
    graph_collision::Ray,
    math,
    node::{Chunk, ChunkId, ChunkLayout, DualGraph},
    proto::Position,
};

/// Performs sphere casting (swept collision query) against the voxels in the `DualGraph`
///
/// The `ray` parameter and any resulting hit normals are given in the local coordinate system of `position`.
///
/// The `tanh_distance` is the hyperbolic tangent of the cast_distance, or the distance along the ray to check for hits.
///
/// This function may return a `SphereCastError` if not enough chunks are generated, even if the ray never reaches an
/// ungenerated chunk. To prevent these errors, make sure that the distance between the ray's start point and the center of
/// the closest node with ungenerated chunks is greater than `cast_distance + collider_radius + dodeca::BOUNDING_SPHERE_RADIUS`
pub fn ray_cast(
    graph: &DualGraph,
    layout: &ChunkLayout,
    position: &Position,
    ray: &Ray,
    tanh_distance: f32,
) -> Result<Option<GraphCastHit>, RayCastError> {
    // A collision check is assumed to be a miss until a collision is found.
    // This `hit` variable gets updated over time before being returned.
    let mut hit: Option<GraphCastHit> = None;

    // Pick the vertex closest to position.local as the vertex of the chunk to use to start collision checking
    let start_vertex = Vertex::iter()
        .map(|v| {
            (
                v,
                (v.node_to_dual().cast::<f32>() * position.local * math::origin()).w,
            )
        })
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap()
        .0;

    // Start a breadth-first search of the graph's chunks, performing collision checks in each relevant chunk.
    // The `chunk_queue` contains ordered pairs containing the `ChunkId` and the transformation needed to switch
    // from the original node coordinates to the current chunk's node coordinates.
    let mut visited_chunks = FxHashSet::<ChunkId>::default();
    let mut chunk_queue: VecDeque<(ChunkId, na::Matrix4<f32>)> = VecDeque::new();
    chunk_queue.push_back((ChunkId::new(position.node, start_vertex), position.local));

    // Precalculate the chunk boundaries for collision purposes. If the collider goes outside these bounds,
    // the corresponding neighboring chunk will also be used for collision checking.
    let klein_lower_boundary = 0.0;
    let klein_upper_boundary = Vertex::chunk_to_dual_factor() as f32;

    // Breadth-first search loop
    while let Some((chunk, node_transform)) = chunk_queue.pop_front() {
        let Chunk::Populated {
            voxels: ref voxel_data,
            ..
        } = graph[chunk]
        else {
            // Collision checking on unpopulated chunk
            return Err(RayCastError::OutOfBounds);
        };
        let local_ray = chunk.vertex.node_to_dual().cast::<f32>() * node_transform * ray;

        // Check collision within a single chunk
        let current_tanh_distance = hit.as_ref().map_or(tanh_distance, |hit| hit.tanh_distance);
        hit = chunk_ray_cast(voxel_data, layout, &local_ray, current_tanh_distance).map_or(
            hit,
            |hit| {
                Some(GraphCastHit {
                    tanh_distance: hit.tanh_distance,
                    chunk,
                    voxel_coords: hit.voxel_coords,
                    face_axis: hit.face_axis,
                    face_direction: hit.face_direction,
                })
            },
        );

        // Compute the Klein-Beltrami coordinates of the ray segment's endpoints. To check whether neighboring chunks
        // are needed, we need to check whether the endpoints of the line segments lie outside the boundaries of the square
        // bounded by `klein_lower_boundary` and `klein_upper_boundary`.
        let klein_ray_start = na::Point3::from_homogeneous(local_ray.position).unwrap();
        let klein_ray_end =
            na::Point3::from_homogeneous(local_ray.ray_point(current_tanh_distance)).unwrap();

        // Add neighboring chunks as necessary based on a conservative AABB check, using one coordinate at a time.
        for axis in 0..3 {
            // Check for neighboring nodes
            if klein_ray_start[axis] <= klein_lower_boundary
                || klein_ray_end[axis] <= klein_lower_boundary
            {
                let side = chunk.vertex.canonical_sides()[axis];
                let next_node_transform = side.reflection().cast::<f32>() * node_transform;
                // If we have to do collision checking on nodes that don't exist in the graph, we cannot have a conclusive result.
                let Some(neighbor) = graph.neighbor(chunk.node, side) else {
                    // Collision checking on nonexistent node
                    return Err(RayCastError::OutOfBounds);
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

    Ok(hit)
}

#[derive(Debug)]
pub enum RayCastError {
    OutOfBounds,
}

/// Information about the intersection at the end of a ray segment.
#[derive(Debug)]
pub struct GraphCastHit {
    /// The tanh of the distance traveled along the ray to result in this hit.
    pub tanh_distance: f32,

    /// Which chunk in the graph the hit occurred in
    pub chunk: ChunkId,

    /// The coordinates of the block that was hit, including margins.
    pub voxel_coords: [usize; 3],

    /// Which of the three axes is orthogonal to the face of the block that was hit.
    pub face_axis: u32,

    // Either +1 or -1, depending on whether the outside of the face that was hit was in the positive or
    // negative direction in `face_axis`.
    pub face_direction: i32,
}
