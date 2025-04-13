use libm::{acosf, cosf, sinf, sqrtf};
use rand::{Rng, SeedableRng};
use rand_distr::Poisson;
use rand_pcg::Pcg64Mcg;

use crate::{
    dodeca::{Side, Vertex},
    graph::{Graph, NodeId},
    math::MVector,
    node::VoxelData,
    peer_traverser::PeerTraverser,
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
        // TODO: Consider adding epsilon to ensure floating point precision
        // doesn't cause `average_with` to fail
        self.vector.mip(side.normal()) > -1.0
    }

    pub fn propagate(&self, side: Side) -> Horosphere {
        Horosphere {
            owner: self.owner,
            vector: side.reflection() * self.vector,
        }
    }

    pub fn average_with(&mut self, other: Horosphere, other_weight: f32) {
        if self.owner != other.owner {
            panic!("Tried to average two unrelated horospheres");
        }
        self.vector = self.vector * (1.0 - other_weight) + other.vector * other_weight;
    }

    pub fn renormalize(&mut self) {
        self.vector.w = self.vector.xyz().norm();
    }

    /// Returns whether the horosphere is freshly created, or whether it's a
    /// reference to a horosphere created earlier on in the node graph.
    pub fn is_fresh(&self, node_id: NodeId) -> bool {
        self.owner == node_id
    }

    /// If self and other have to compete to exist as an actual horosphere,
    /// returns whether self wins.
    pub fn has_priority(&self, other: &Horosphere, node_id: NodeId) -> bool {
        // If both horospheres are fresh, use the w-coordinate as an arbitrary
        // tie-breaker to decide which horosphere should win.
        !self.is_fresh(node_id) || (other.is_fresh(node_id) && self.vector.w < other.vector.w)
    }

    /// Based on other nodes in the graph, determines whether the horosphere
    /// should generate. If false, it means that another horosphere elsewhere
    /// would interfere, and generation should not proceed.
    pub fn should_generate(&self, graph: &Graph, node_id: NodeId) -> bool {
        if !self.is_fresh(node_id) {
            // The horosphere is propagated and so is already proven to exist.
            return true;
        }

        let mut peers = PeerTraverser::new(node_id);
        while let Some(peer) = peers.next(graph) {
            let Some(peer_horosphere) = graph
                .partial_node_state(peer.node())
                .candidate_horosphere
                .as_ref()
            else {
                continue;
            };
            if !self.has_priority(peer_horosphere, node_id)
                // Check that these horospheres can interfere by seeing if they share a node in common.
                && peer_horosphere.should_propagate_through_path(peer.path_from_peer())
                && self.should_propagate_through_path(peer.path_from_base())
            {
                return false;
            }
        }
        true
    }

    fn should_propagate_through_path(&self, mut path: impl ExactSizeIterator<Item = Side>) -> bool {
        let mut current_horosphere = self.clone();
        while let Some(side) = path.next() {
            if !current_horosphere.should_propagate(side) {
                return false;
            }
            if path.len() == 0 {
                return true;
            }
            current_horosphere = current_horosphere.propagate(side);
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

#[cfg(test)]
#[allow(unused)]
mod test {
    use std::collections::BTreeMap;

    use super::*;
    use crate::{
        dodeca::Side,
        math::MIsometry,
        proto::Position,
        traversal::{self, ensure_nearby, nearby_nodes},
    };
    use Side::{A, B, C, D, E, F, G, H, I, J, K, L};
    use fxhash::{FxHashMap, FxHashSet};

    #[test]
    #[rustfmt::skip]
    fn average_with_test() {
        /*
        [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, G, E, C]
        !=
        [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, E, F, B]
        for
        [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, G, E, F, C, B]
        */


        /*
        [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, E, G, C]
        !=
        [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, E, F, B]
        for
        [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, E, G, C, F, B]

        Sibling relationship path: CBGF
        */
        let mut graph = Graph::new(12);
        let mut node_id = NodeId::ROOT;
        for side in [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, G, E, C] {
            node_id = graph.ensure_neighbor(node_id, side);
            graph.ensure_node_state(node_id);
            ensure_nearby(&mut graph, &Position { node: node_id, local: MIsometry::identity() }, 3.0);
        }

        let mut node_id = NodeId::ROOT;
        for side in [B, D, E, B, D, E, C, B, I, C, B, I, D, E, C, E, C, B, A, L, E, F, B] {
            node_id = graph.ensure_neighbor(node_id, side);
            graph.ensure_node_state(node_id);
            ensure_nearby(&mut graph, &Position { node: node_id, local: MIsometry::identity() }, 3.0);
        }

        let all_nodes = nearby_nodes(&graph, &Position::origin(), f32::INFINITY);
        for (node_id, _) in all_nodes {
            println!("{:?}", graph.node_path(node_id));
            graph.ensure_node_state(node_id);
        }
    }

    fn is_shortlex(last: Side, rest: &[Side]) -> bool {
        let Some(second_last) = rest.last() else {
            // One-element node strings are always shortlex.
            return true;
        };
        if last == *second_last {
            // Backtracking is not valid (short part of shortlex)
            return false;
        }
        if last.adjacent_to(*second_last) && (last as usize) < (*second_last as usize) {
            // Unnecessarily having a higher side index first is not valid (lex part of shortlex)
            return false;
        }
        true
    }

    // Like is_shortlex, but without the lex part. Slightly more complicated because it's no longer
    // sufficient to check for immediate backtracking
    fn is_shortest_path(last: Side, rest: &[Side]) -> bool {
        for &side in rest.iter().rev() {
            if last == side {
                // Backtracking discovered
                return false;
            }
            if !last.adjacent_to(side) {
                // Chain of adjacencies has ended, so there's no way to shorten the path.
                return true;
            }
        }
        true
    }

    fn is_interfering_path(
        last: Side,
        rest: &[Side],
        depth: usize,
        graph: &Graph,
        node: NodeId,
    ) -> bool {
        if !is_shortest_path(last, rest) {
            return false;
        }
        if rest.len() >= depth {
            for &side in &rest[0..depth] {
                if !last.adjacent_to(side) && last != side {
                    return false;
                }
            }
        }
        if rest.len() < depth {
            let ancestor_node = (rest.iter().chain(Some(last).iter()))
                .try_fold(node, |current_node, side| {
                    graph.neighbor(current_node, *side)
                });
            let Some(ancestor_node) = ancestor_node else {
                return false;
            };
            if graph.length(ancestor_node) + rest.len() as u32 + 1 != graph.length(node) {
                return false;
            }
        }
        true
    }

    // Returns false if we were unable to increment the path because we exhausted all possibilities
    fn increment_path(
        path: &mut [Side],
        mut always_increment: bool,
        validity_check: &impl Fn(Side, &[Side]) -> bool,
    ) -> bool {
        if path.is_empty() {
            // Empty paths always pass validity checks, so looping is only necessary
            // if incrementing was requested.
            return !always_increment;
        }
        let (last, rest) = path.split_last_mut().unwrap();
        while !validity_check(*last, rest) || always_increment {
            always_increment = false;
            *last = Side::VALUES[(*last as usize + 1) % Side::VALUES.len()];
            if *last == Side::A && !increment_path(rest, true, validity_check) {
                return false;
            }
        }
        true
    }

    struct NodeString {
        path: Vec<Side>,
        initialized: bool,
    }

    impl NodeString {
        fn new(len: usize) -> Self {
            NodeString {
                path: vec![Side::A; len],
                initialized: false,
            }
        }

        fn increment(&mut self, validity_check: &impl Fn(Side, &[Side]) -> bool) -> bool {
            if !self.initialized {
                for i in 0..self.path.len() {
                    if !increment_path(&mut self.path[0..=i], false, validity_check) {
                        return false;
                    }
                }
                self.initialized = true;
                true
            } else {
                increment_path(&mut self.path, true, validity_check)
            }
        }
    }

    #[test]
    fn enumerate_siblings() {
        // If BF...A.... and CG.... reach the same place, then the CG string must have an A in it.
        // Transforming one path to the other will still require transferring B and F past C and G
        // So, strings x and y can interfere iff everything in x commutes with everything in y.
        // However, x and y are only worth considering if they don't share a common parent (other than)
        // the root. This is a bit complicated because CG and GC don't share a parent, but BF and FB do (because they commute).

        // Actually, checking for redundancies between BF and CG is tricky, but there's a simpler way to think about it:
        // These two paths are redundant iff FBCG can be shortened.
        // Perhaps, the problem is as simple as follows:
        // Find all shortlex node strings of a particular even length where the first half commutes with the second half.

        let mut num_sibling_nodes: BTreeMap<usize, u32> = BTreeMap::new();

        let mut graph = Graph::new(1);
        ensure_nearby(&mut graph, &Position::origin(), 5.0);
        let base_nodes = nearby_nodes(&graph, &Position::origin(), 5.0);
        for (base_node, _) in base_nodes {
            /*let mut base_node = NodeId::ROOT;
            for side in [Side::A, Side::B, Side::C, Side::D] {
                base_node = graph.ensure_neighbor(base_node, side);
            }*/

            let depth = 1;

            let mut path = NodeString::new(depth * 2);
            let mut sibling_nodes = FxHashSet::default();
            //println!("List of depth-{} interference patterns", depth);
            while path
                .increment(&|last, rest| is_interfering_path(last, rest, depth, &graph, base_node))
            {
                let found_node_id = path.path.iter().fold(base_node, |current_node, side| {
                    graph.ensure_neighbor(current_node, *side)
                });
                if graph.length(found_node_id) != graph.length(base_node) {
                    continue;
                }
                if sibling_nodes.insert(found_node_id) {
                    //println!("{:?}: {:?}", path.path, found_node_id);
                }
            }
            //println!("End of list");
            //println!("{:?}: {}", base_node, sibling_nodes.len());
            *num_sibling_nodes.entry(sibling_nodes.len()).or_default() += 1;
        }

        println!("Sibling node count histogram: {:?}", num_sibling_nodes);
    }

    // BF vs CG
}
