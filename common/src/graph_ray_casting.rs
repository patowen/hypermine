use crate::{
    chunk_ray_casting::chunk_ray_cast,
    collision_math::Ray,
    graph::Graph,
    node::{Chunk, ChunkId, Coords},
    proto::Position,
    traversal::RayTraverser,
};

/// Performs ray casting against the voxels in the `DualGraph`
///
/// The `ray` parameter and any resulting hit normals are given in the local coordinate system of `position`.
///
/// The `tanh_distance` is the hyperbolic tangent of the cast_distance, or the distance along the ray to check for hits.
///
/// This function may return a `RayCastError` if not enough chunks are generated, even if the ray never reaches an
/// ungenerated chunk. To prevent these errors, make sure that the distance between the ray's start point and the center of
/// the closest node with ungenerated chunks is greater than `cast_distance + dodeca::BOUNDING_SPHERE_RADIUS`
pub fn ray_cast(
    graph: &Graph,
    position: &Position,
    ray: &Ray,
    mut tanh_distance: f32,
) -> Result<Option<GraphCastHit>, RayCastError> {
    // A ray cast is assumed to be a miss until a collision is found.
    // This `hit` variable gets updated over time before being returned.
    let mut hit: Option<GraphCastHit> = None;

    let mut traverser = RayTraverser::new(graph, *position, ray, 0.0);
    while let Some((chunk, transform)) = traverser.next(tanh_distance) {
        let Some(chunk) = chunk else {
            // Ray reached chunk outside of graph
            return Err(RayCastError::OutOfBounds);
        };
        let Chunk::Populated {
            voxels: ref voxel_data,
            ..
        } = graph[chunk]
        else {
            // Ray reached unpopulated chunk
            return Err(RayCastError::OutOfBounds);
        };

        hit = chunk_ray_cast(
            voxel_data,
            graph.layout(),
            &(transform * ray),
            tanh_distance,
        )
        .map_or(hit, |hit| {
            tanh_distance = hit.tanh_distance;
            Some(GraphCastHit {
                tanh_distance: hit.tanh_distance,
                chunk,
                voxel_coords: hit.voxel_coords,
                face_axis: hit.face_axis,
                face_direction: hit.face_direction,
            })
        });
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
    pub voxel_coords: Coords,

    /// Which of the three axes is orthogonal to the face of the block that was hit.
    pub face_axis: u32,

    // Either +1 or -1, depending on whether the outside of the face that was hit was in the positive or
    // negative direction in `face_axis`.
    pub face_direction: i8,
}