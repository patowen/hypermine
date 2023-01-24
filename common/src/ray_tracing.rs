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
    ray: na::Matrix4x2<f32>,
    tanh_length: f32,
) {
    let mut status = RayStatus {
        result: RayTracingResult::Miss,
        tanh_length,
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
                return;
            };
        let local_ray = chunk.vertex.node_to_dual().cast::<f32>() * node_transform * ray;
        chunk_ray_tracer.trace_ray(
            &RtChunkContext {
                dimension,
                chunk,
                transform: transform
                    * math::mtranspose(&node_transform)
                    * chunk.vertex.dual_to_node().cast(),
                voxel_data,
                ray: local_ray,
            },
            &mut status,
        );

        // If pos or pos+dir*max_t lies beyond the chunk boundary, with a buffer to account for radius, repeat
        // ray tracing with the neighboring chunk unless it has already been visited. We start at vertex
        // AB for simplicity even if that's not where pos is, although this should be optimized later.
        for coord_boundary in 0..3 {
            let klein_pos0_val = local_ray[(coord_boundary, 0)] / local_ray[(coord_boundary, 3)];
            let klein_pos1_val = (local_ray[(coord_boundary, 0)]
                + local_ray[(coord_boundary, 1)] * status.tanh_length)
                / (local_ray[(3, 0)] + local_ray[(3, 1)] * status.tanh_length);

            // Check for neighboring nodes. TODO: The use of unwrap here will cause a crash if you arrive at an ungenerated chunk.
            if klein_pos0_val <= klein_boundary0 || klein_pos1_val <= klein_boundary0 {
                let side = chunk.vertex.canonical_sides()[coord_boundary];
                let next_chunk = (graph.neighbor(chunk.node, side).unwrap(), chunk.vertex).into();
                if visited_chunks.insert(next_chunk) {
                    chunk_queue
                        .push_back((next_chunk, side.reflection().cast::<f32>() * transform));
                }
            }

            // Check for neighboring chunks within the same node
            if klein_pos0_val >= klein_boundary1 || klein_pos1_val >= klein_boundary1 {
                let vertex = chunk.vertex.adjacent_vertices()[coord_boundary];
                let next_chunk = (chunk.node, vertex).into();
                if visited_chunks.insert(next_chunk) {
                    chunk_queue.push_back((next_chunk, transform));
                }
            }
        }
    }
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
    pub ray: na::Matrix4x2<f32>,
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
    pub tanh_length: f32,
    pub result: RayTracingResult,
}

pub enum RayTracingResult {
    Miss,
    Intersection(RayTracingIntersection),
    Inconclusive,
}

pub struct RayTracingIntersection {
    pub chunk: ChunkId,
    pub normal: na::Vector4<f32>,
}

impl RayStatus {
    pub fn update(
        &mut self,
        context: &RtChunkContext<'_>,
        tanh_length: f32,
        normal: na::Vector4<f32>,
    ) {
        self.tanh_length = tanh_length;
        self.result = RayTracingResult::Intersection(RayTracingIntersection {
            chunk: context.chunk,
            normal: context.transform * normal,
        });
    }
}
