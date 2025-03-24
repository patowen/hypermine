use libm::{acosf, cosf, sinf, sqrtf};
use rand::Rng;
use rand_distr::Poisson;
use rand_pcg::Pcg64Mcg;

use crate::{
    dodeca::{Side, Vertex},
    graph::{Graph, NodeId},
    math::MVector,
    node::VoxelData,
    voxel_math::Coords,
    world::Material,
};

#[derive(Clone, Copy)]
pub struct Horosphere {
    pub owner: NodeId,
    pub vector: MVector<f32>, // TODO: Explain (equation of horosphere is `mip(h,p) == -1`)
}

impl Horosphere {
    pub fn propagate(&self, side: Side) -> Option<Horosphere> {
        if self.vector.mip(side.normal()) < 1.0 {
            Some(Horosphere {
                owner: self.owner,
                vector: *side.reflection() * self.vector,
            })
        } else {
            None
        }
    }

    pub fn average_with(&mut self, other: Horosphere, other_weight: f32) {
        assert!(self.owner == other.owner);
        self.vector = self.vector * (1.0 - other_weight) + other.vector * other_weight;
    }

    pub fn renormalize(&mut self) {
        self.vector.w = self.vector.xyz().norm();
    }

    pub fn chunk_data(&self, vertex: Vertex) -> HorosphereChunk {
        HorosphereChunk {
            vector: *vertex.node_to_dual() * self.vector,
        }
    }
}

pub struct HorosphereChunk {
    pub vector: MVector<f32>,
}

impl HorosphereChunk {
    pub fn generate(&self, voxels: &mut VoxelData, chunk_size: u8) {
        for x in 0..chunk_size {
            for y in 0..chunk_size {
                for z in 0..chunk_size {
                    // TODO: This math should be in ChunkLayout
                    let pos = MVector::new(
                        x as f32 + 0.5,
                        y as f32 + 0.5,
                        z as f32 + 0.5,
                        chunk_size as f32 * Vertex::dual_to_chunk_factor(),
                    )
                    .normalized_point();
                    if pos.mip(&self.vector) > -1.0 {
                        voxels.data_mut(chunk_size)[Coords([x, y, z]).to_index(chunk_size)] =
                            Material::RedSandstone;
                    }
                }
            }
        }
    }
}

pub fn get_random_candidate_horosphere(
    rng: &mut Pcg64Mcg,
    graph: &Graph,
    node: NodeId,
) -> Option<MVector<f32>> {
    for _ in 0..rng.sample(Poisson::new(6.0).unwrap()) as u32 {
        let horosphere_pos = random_horosphere(rng);
        if is_horosphere_valid(graph, node, &horosphere_pos) {
            return Some(horosphere_pos);
        }
    }
    None
}

fn is_horosphere_valid(graph: &Graph, node: NodeId, horosphere: &MVector<f32>) -> bool {
    // TODO: This needs an explanation.
    Side::iter().all(|s| s.normal().mip(&horosphere) < 1.0)
        && (graph.descenders(node)).all(|(s, _)| s.normal().mip(horosphere) < -1.0)
}

fn random_horosphere(rng: &mut Pcg64Mcg) -> MVector<f32> {
    // TODO: This needs an explanation.
    let vertex_w = sqrtf(2.0) * (3.0 + sqrtf(5.0)) / 4.0; // w-coordinate of every vertex in dodeca-coordinates
    let max_w = sqrtf(3.0) * (vertex_w + sqrtf(vertex_w * vertex_w - 1.0)); // Maximum possible w-coordinate of valid horosphere

    let w = sqrtf(rng.random::<f32>()) * max_w;
    let phi = acosf(rng.random::<f32>() * 2.0 - 1.0);
    let theta = rng.random::<f32>() * std::f32::consts::TAU;
    MVector::new(
        w * sinf(phi) * cosf(theta),
        w * sinf(phi) * sinf(theta),
        w * cosf(phi),
        w,
    )
}
