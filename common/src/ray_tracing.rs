use std::collections::{HashSet, VecDeque};

use crate::{
    dodeca::Vertex,
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
    let mut status = RayStatus {
        result: RayTracingResult::Miss,
        tanh_distance,
    };

    let mut visited_chunks: HashSet<ChunkId> = HashSet::new();
    let mut chunk_queue: VecDeque<(ChunkId, na::Matrix4<f32>)> = VecDeque::new();
    chunk_queue.push_back((chunk, na::Matrix4::identity()));

    const EPSILON: f32 = 1e-5;
    let klein_boundary0 = (chunk_ray_tracer.max_radius() + EPSILON).tanh();
    let klein_boundary1 = ((Vertex::chunk_to_dual_factor() as f32).atanh()
        - (chunk_ray_tracer.max_radius() + EPSILON))
        .tanh();

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
        chunk_ray_tracer.trace_ray(
            &RtChunkContext {
                dimension,
                chunk,
                transform: transform
                    * math::mtranspose(&node_transform)
                    * chunk.vertex.dual_to_node().cast(),
                voxel_data,
                ray: &local_ray,
            },
            &mut status,
        );

        let klein_ray_start = na::Point3::from_homogeneous(local_ray.position).unwrap();
        let klein_ray_end =
            na::Point3::from_homogeneous(local_ray.point(status.tanh_distance)).unwrap();

        // If pos or pos+dir*max_t lies beyond the chunk boundary, with a buffer to account for radius, repeat
        // ray tracing with the neighboring chunk unless it has already been visited. We start at vertex
        // AB for simplicity even if that's not where pos is, although this should be optimized later.
        for coord_boundary in 0..3 {
            // Check for neighboring nodes
            if klein_ray_start[coord_boundary] <= klein_boundary0
                || klein_ray_end[coord_boundary] <= klein_boundary0
            {
                let side = chunk.vertex.canonical_sides()[coord_boundary];
                let Some(neighbor) = graph.neighbor(chunk.node, side) else {
                    // Collision checking on nonexistent node
                    status.result = RayTracingResult::Inconclusive;
                    return status;
                };
                let next_chunk = (neighbor, chunk.vertex).into();
                if visited_chunks.insert(next_chunk) {
                    chunk_queue
                        .push_back((next_chunk, side.reflection().cast::<f32>() * node_transform));
                }
            }

            // Check for neighboring chunks within the same node
            if klein_ray_start[coord_boundary] >= klein_boundary1
                || klein_ray_end[coord_boundary] >= klein_boundary1
            {
                let vertex = chunk.vertex.adjacent_vertices()[coord_boundary];
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
    pub chunk: ChunkId,
    pub transform: na::Matrix4<f32>,
    pub voxel_data: &'a VoxelData,
    pub ray: &'a Ray,
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
    pub fn point(&self, tanh_distance: f32) -> na::Vector4<f32> {
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
