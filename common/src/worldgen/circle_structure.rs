use libm::{acosf, cosf, sinf, sqrtf};
use rand::{Rng, SeedableRng};
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

#[derive(Clone)]
pub struct Horosphere {
    pub owner: NodeId,
    pub vector: MVector<f32>, // TODO: Explain (equation of horosphere is `mip(h,p) == -1`)
}

impl Horosphere {
    pub fn create_from_parents(graph: &Graph, node_id: NodeId) -> Option<Horosphere> {
        let mut horospheres_to_average_iter =
            graph
                .descenders(node_id)
                .filter_map(|(parent_side, parent_id)| {
                    (graph.node_state(parent_id).horosphere.as_ref())
                        .filter(|h| h.should_propagate(parent_side))
                        .map(|h| h.propagate(parent_side))
                });

        let mut horosphere = horospheres_to_average_iter.next()?;
        let mut count = 1;
        for other in horospheres_to_average_iter {
            count += 1;
            horosphere.average_with(other, 1.0 / count as f32);
        }

        horosphere.renormalize();
        Some(horosphere)
    }

    pub fn maybe_create_fresh(graph: &Graph, node_id: NodeId) -> Option<Horosphere> {
        let spice = graph.hash_of(node_id) as u64;
        let mut rng = rand_pcg::Pcg64Mcg::seed_from_u64(spice.wrapping_add(42));
        for _ in 0..rng.sample(Poisson::new(6.0).unwrap()) as u32 {
            let horosphere_pos = Self::random_horosphere_pos(&mut rng);
            if Self::is_horosphere_pos_valid(graph, node_id, &horosphere_pos) {
                return Some(Horosphere {
                    owner: node_id,
                    vector: horosphere_pos,
                });
            }
        }
        None
    }

    pub fn should_propagate(&self, side: Side) -> bool {
        self.vector.mip(side.normal()) > -1.0
    }

    pub fn propagate(&self, side: Side) -> Horosphere {
        Horosphere {
            owner: self.owner,
            vector: side.reflection() * self.vector,
        }
    }

    pub fn average_with(&mut self, other: Horosphere, other_weight: f32) {
        // TODO: This assertion can fail. May need fuzz testing.
        if self.owner != other.owner {
            tracing::error!("average_with failed. {:?} != {:?}", self.owner, other.owner);
            panic!("average_with failed");
        }
        self.vector = self.vector * (1.0 - other_weight) + other.vector * other_weight;
    }

    pub fn renormalize(&mut self) {
        self.vector.w = self.vector.xyz().norm();
    }

    /// Returns whether the structure is freshly created, or whether it's a
    /// reference to a structure created earlier on in the node graph.
    pub fn is_fresh(&self, node_id: NodeId) -> bool {
        self.owner == node_id
    }

    /// If self and other have to compete to exist as an actual structure,
    /// returns whether self wins.
    pub fn has_priority(&self, other: &Horosphere, node_id: NodeId) -> bool {
        // If both structures are fresh, use the w-coordinate as an arbitrary
        // tie-breaker to decide which horosphere should win.
        !self.is_fresh(node_id) || (other.is_fresh(node_id) && self.vector.w < other.vector.w)
    }

    /// Based on other nodes in the graph, determines whether the structure
    /// should generate. If false, it means that another horosphere elsewhere
    /// would interfere, and generation should not proceed.
    pub fn should_generate(&self, graph: &Graph, node_id: NodeId) -> bool {
        if !self.is_fresh(node_id) {
            // The horosphere is propagated and so is already proven to exist.
            return true;
        }

        let length = graph.length(node_id);
        for (parent_side, parent_id) in graph.descenders(node_id) {
            for sibling_side in Side::iter().filter(|s| s.adjacent_to(parent_side)) {
                let sibling_id = graph.neighbor(parent_id, sibling_side).unwrap();
                if graph.length(sibling_id) != length {
                    continue;
                }
                let Some(sibling_horosphere) = graph
                    .partial_node_state(sibling_id)
                    .candidate_horosphere
                    .as_ref()
                else {
                    continue;
                };
                if !self.has_priority(sibling_horosphere, node_id)
                    // Check that these structures can interfere by seeing if they share a node in common.
                    && sibling_horosphere.should_propagate(parent_side)
                    && self.should_propagate(sibling_side)
                {
                    return false;
                }
            }
        }
        true
    }

    pub fn chunk_data(&self, vertex: Vertex) -> HorosphereChunk {
        HorosphereChunk {
            vector: vertex.node_to_dual() * self.vector,
        }
    }

    /// Returns whether the given horosphere position could represent a horosphere generated by the
    /// given node. The requirement is that a horosphere must be bounded by all of the node's descenders
    /// (as otherwise, a parent node would own the horosphere), and the horosphere must not be fully
    /// behind any of the other dodeca sides (as otherwise, a child node would own the horosphere). Note
    /// that the horosphere does not necessarily need to intersect the dodeca to be valid.
    fn is_horosphere_pos_valid(graph: &Graph, node_id: NodeId, horosphere: &MVector<f32>) -> bool {
        // TODO: This needs an explanation.
        Side::iter().all(|s| s.normal().mip(&horosphere) < 1.0)
            && (graph.descenders(node_id)).all(|(s, _)| s.normal().mip(horosphere) < -1.0)
    }

    /// Returns a vector representing a random horosphere close enough to the origin that it would
    /// intersect with the sphere that a dodeca is inscribed in.
    fn random_horosphere_pos(rng: &mut Pcg64Mcg) -> MVector<f32> {
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
