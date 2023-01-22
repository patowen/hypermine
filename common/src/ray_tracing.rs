use crate::{node::{ChunkId, VoxelData}, world::Material};

pub struct RtChunkContext<'a> {
    pub dimension: usize,
    pub chunk: ChunkId,
    pub transform: na::Matrix4<f32>,
    pub voxel_data: &'a VoxelData,
    pub ray: na::Matrix2x4<f32>,
}

impl RtChunkContext<'_> {
    pub fn get_voxel(&self, coords: [usize; 3]) -> Material {
        assert!(coords[0] < self.dimension);
        assert!(coords[1] < self.dimension);
        assert!(coords[2] < self.dimension);
        let dimension_with_margin = self.dimension + 2;
        self.voxel_data.get(
            (coords[0] + 1)
                + (coords[1] + 1) * dimension_with_margin
                + (coords[2] + 1) * dimension_with_margin.pow(2),
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
