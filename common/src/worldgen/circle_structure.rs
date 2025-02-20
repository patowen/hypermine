use libm::{acosf, cosf, sinf, sqrtf};
use rand::Rng;
use rand_distr::Poisson;
use rand_pcg::Pcg64Mcg;

use crate::{
    dodeca::Side,
    graph::{Graph, NodeId},
    math::MVector,
};

pub struct Horosphere {
    pub owner: NodeId,
    pub pos: MVector<f32>, // TODO: Explain
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
