use crate::ray_tracing::{ChunkRayTracer, RayStatus, RtChunkContext};

pub struct SphereCollider {}

impl ChunkRayTracer for SphereCollider {
    fn trace_ray(&self, ctx: &RtChunkContext, status: &mut RayStatus) {
        todo!()
    }

    fn max_radius(&self) -> f32 {
        todo!()
    }
}

impl SphereCollider {
    fn find_vertex_collision(
        &self,
        ctx: &RtChunkContext,
        bbox: [[usize; 2]; 3],
        status: &mut RayStatus,
    ) {
        for i in bbox[0][0]..bbox[0][1] {
            for j in bbox[1][0]..bbox[1][1] {
                for k in bbox[2][0]..bbox[2][1] {
                }
            }
        }
    }
}
