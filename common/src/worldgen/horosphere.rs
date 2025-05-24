use libm::{cosf, sinf, sqrtf};
use rand::{Rng, SeedableRng};
use rand_distr::Poisson;
use rand_pcg::Pcg64Mcg;

use crate::{
    dodeca::{Side, Vertex},
    graph::{Graph, NodeId},
    math::MVector,
    node::VoxelData,
    peer_traverser,
    voxel_math::Coords,
    world::Material,
};

/// Whether an assortment of random horospheres should be added to world generation. This is a temporary
/// option until large structures that fit with the theme of the world are introduced.
/// For code simplicity, this is made into a constant instead of a configuration option.
const HOROSPHERES_ENABLED: bool = true;

/// Represents a node's reference to a particular horosphere. As a general rule, for any give horosphere,
/// every node in the convex hull of nodes containing the horosphere will have a `HorosphereNode`
/// referencing it. The unique node in this convex hull with the smallest depth in the graph is the owner
/// of the horosphere, where it is originally generated.
#[derive(Copy, Clone)]
pub struct HorosphereNode {
    /// The node that originally created the horosphere. All parts of the horosphere will
    /// be in a node with this as an ancestor, and all HorosphereNodes with the same `owner` correspond
    /// to the same horosphere.
    owner: NodeId,

    /// The vector representing the horosphere with respect to the node storing this `HorosphereNode`. A vector
    /// `point` is in this horosphere if `point.mip(&self.pos) == -1`. This vector should always have
    /// the invariant `self.pos.mip(&self.pos) == 0`, behaving much like a "light-like" vector
    /// in Minkowski space. One consequence of this invariant is that this vector's length is always
    /// proportional to its w-coordinate. If the w-coordinate is 1, the horosphere intersects the origin.
    /// If it's less than 1, the horosphere contains the origin, and if it's greater than 1, the origin
    /// is outside the horosphere. The vector points in the direction of the horosphere's ideal point.
    ///
    /// TODO: If a player traverses too far inside a horosphere, this vector will underflow, preventing
    /// the horosphere from generating properly. Fixing this requires using logic similar to `Plane` to
    /// increase the range of magnitudes the vector can take.
    pos: MVector<f32>,
}

impl HorosphereNode {
    /// Returns the `HorosphereNode` for the given node, either by propagating an existing parent
    /// `HorosphereNode` or by randomly generating a new one.
    pub fn new(graph: &Graph, node_id: NodeId) -> Option<HorosphereNode> {
        if !HOROSPHERES_ENABLED {
            return None;
        }
        HorosphereNode::create_from_parents(graph, node_id)
            .or_else(|| HorosphereNode::maybe_create_fresh(graph, node_id))
    }

    /// Propagates `HorosphereNode` information from the given parent nodes to this child node. Returns
    /// `None` if there's no horosphere to propagate, either because none of the parent nodes have a
    /// horosphere associated with them, or because any existing horosphere is outside the range
    /// of this node.
    fn create_from_parents(graph: &Graph, node_id: NodeId) -> Option<HorosphereNode> {
        // Rather than selecting an arbitrary parent HorosphereNode, we average all of them. This
        // is important because otherwise, the propagation of floating point precision errors could
        // create a seam. This ensures that all errors average out, keeping the horosphere smooth.
        let mut horospheres_to_average_iter =
            graph
                .descenders(node_id)
                .filter_map(|(parent_side, parent_id)| {
                    graph
                        .node_state(parent_id)
                        .horosphere
                        .as_ref()
                        .and_then(|h| h.propagate(parent_side))
                });

        let mut horosphere = horospheres_to_average_iter.next()?;
        let mut count = 1;
        for other in horospheres_to_average_iter {
            // Take an average of all HorosphereNodes in this iterator, giving each of them equal weight
            // by keeping track of a moving average with a weight that changes over time to make the
            // numbers work out the same way.
            count += 1;
            horosphere.average_with(other, 1.0 / count as f32);
        }

        horosphere.renormalize();
        Some(horosphere)
    }

    /// Create a `HorosphereNode` corresponding to a freshly created horosphere with the given node as its owner,
    /// if one should be created. This function is called on every node that doesn't already have a horosphere
    /// associated with it, so this function has control over how frequent the horospheres should be.
    fn maybe_create_fresh(graph: &Graph, node_id: NodeId) -> Option<HorosphereNode> {
        const HOROSPHERE_DENSITY: f32 = 6.0;

        let spice = graph.hash_of(node_id) as u64;
        let mut rng = rand_pcg::Pcg64Mcg::seed_from_u64(spice.wrapping_add(42));
        for _ in 0..rng.sample(Poisson::new(HOROSPHERE_DENSITY).unwrap()) as u32 {
            // This logic is designed to create an average of "HOROSPHERE_DENSITY" horosphere candiates
            // in the region determined by `random_horosphere_pos` and then filters the resulting
            // list of candiates to only ones where the current node is the suitable owner for them.
            // Filtering instead of rejection sampling ensures a uniform distribution of horosphere
            // even though different nodes have different-sized regions for valid horospheres.

            // However, we do return early to ensure that after filtering, we only take the first
            // horosphere if there is one, since a node can have at most one horosphere.
            let horosphere_pos = random_horosphere_pos(&mut rng);
            if is_horosphere_pos_valid(graph, node_id, &horosphere_pos) {
                return Some(HorosphereNode {
                    owner: node_id,
                    pos: horosphere_pos,
                });
            }
        }
        None
    }

    /// Returns an estimate of the `HorosphereNode` corresponding to the node adjacent to the current node
    /// at the given side, or `None` if the horosphere is no longer relevant after crossing the given side.
    /// The estimates given by multiple nodes may be used to produce the actual `HorosphereNode`.
    fn propagate(&self, side: Side) -> Option<HorosphereNode> {
        // If the horosphere is entirely behind the plane bounded by the given side, it is no longer relevant.
        // The relationship between a horosphere and a directed plane given by a normal vector can be determined with
        // the Minkowski inner product (mip) between their respective vectors. If it's positive, the ideal point is in front
        // of the plane, and if it's negative, the ideal point is behind it. The horosphere intersects the plane
        // exactly when the mip is between -1 and 1. Therefore, if it's less than -1, the horosphere is entirely
        // behind the plane.
        if self.pos.mip(side.normal()) <= -1.0 {
            return None;
        }

        Some(HorosphereNode {
            owner: self.owner,
            pos: side.reflection() * self.pos,
        })
    }

    /// Takes the weighted average of the coordinates of this horosphere with the coordinates of the other horosphere.
    fn average_with(&mut self, other: HorosphereNode, other_weight: f32) {
        if self.owner != other.owner {
            // If this panic is triggered, it may mean that two horospheres were generated that interfere
            // with each other. The logic in `should_generate` should prevent this, so this would be a sign
            // of a bug in that function's implementation.
            panic!("Tried to average two unrelated HorosphereNodes");
        }
        self.pos = self.pos * (1.0 - other_weight) + other.pos * other_weight;
    }

    /// Ensures that the horosphere invariant holds (`pos.mip(&pos) == 0`), as numerical error can otherwise propagate,
    /// potentially making the surface behave more like a sphere or an equidistant surface.
    fn renormalize(&mut self) {
        self.pos.w = self.pos.xyz().norm();
    }

    /// Returns whether the horosphere is freshly created, instead of a
    /// reference to a horosphere created earlier on in the node graph.
    fn is_fresh(&self, node_id: NodeId) -> bool {
        self.owner == node_id
    }

    /// If `self` and `other` would propagate to the same node, to avoid interference, only one of these
    /// two horospheres can generate. This function determines whether `self` should be the one to generate.
    fn has_priority(&self, other: &HorosphereNode, node_id: NodeId) -> bool {
        // If both horospheres are fresh, use the w-coordinate as an arbitrary
        // tie-breaker to decide which horosphere should win.
        !self.is_fresh(node_id) || (other.is_fresh(node_id) && self.pos.w < other.pos.w)
    }

    /// Based on other nodes in the graph, determines whether the horosphere
    /// should generate. If false, it means that another horosphere elsewhere
    /// would interfere, and generation should not proceed.
    pub fn should_generate(&self, graph: &Graph, node_id: NodeId) -> bool {
        if !self.is_fresh(node_id) {
            // The horosphere is propagated and so is already proven to exist.
            return true;
        }

        for peer in peer_traverser::expect_peer_nodes(graph, node_id) {
            let Some(peer_horosphere) = graph
                .partial_node_state(peer.node())
                .candidate_horosphere
                .as_ref()
            else {
                continue;
            };
            if !self.has_priority(peer_horosphere, node_id)
                // Check that these horospheres can interfere by seeing if they share a node in common.
                && peer_horosphere.should_propagate_through_path(peer.peer_to_shared())
                && self.should_propagate_through_path(peer.base_to_shared())
            {
                return false;
            }
        }
        true
    }

    /// This function returns whether the horosphere will propagate all the way to the node given
    /// by the sequence of sides.
    fn should_propagate_through_path(&self, path: impl ExactSizeIterator<Item = Side>) -> bool {
        let mut current_horosphere = *self;
        for side in path {
            let Some(horosphere) = current_horosphere.propagate(side) else {
                return false;
            };
            current_horosphere = horosphere;
        }
        true
    }
}

/// Returns whether the given horosphere position could represent a horosphere generated by the
/// given node. The requirement is that a horosphere must be bounded by all of the node's descenders
/// (as otherwise, a parent node would own the horosphere), and the horosphere must not be fully
/// behind any of the other dodeca sides (as otherwise, a child node would own the horosphere). Note
/// that the horosphere does not necessarily need to intersect the dodeca to be valid.
fn is_horosphere_pos_valid(graph: &Graph, node_id: NodeId, horosphere_pos: &MVector<f32>) -> bool {
    // See `propagate`'s implementation for an explanation of what the mip between a horosphere position and
    // a plane's normal signifies.
    Side::iter().all(|s| s.normal().mip(&horosphere_pos) < 1.0)
        && (graph.descenders(node_id)).all(|(s, _)| s.normal().mip(horosphere_pos) < -1.0)
}

/// Returns a vector representing a uniformly random horosphere within a certain distance to the origin.
/// This distance is set up to ensure that `is_horosphere_pos_valid` would always return false if it were any futher,
/// ensuring that this function does not artificially restrict which horospheres can be created.
fn random_horosphere_pos(rng: &mut Pcg64Mcg) -> MVector<f32> {
    // Pick a w-coordinate whose probability density function is `p(w) = w`. By trial and error,
    // this seems to produce horospheres with a uniform and isotropic distribution.
    // TODO: Find a rigorous explanation for this. We would want to show that the probability density is unchanged
    // when an isometry is applied.
    let w = sqrtf(rng.random::<f32>()) * MAX_OWNED_HOROSPHERE_W;

    // Uniformly pick spherical coordinates from a unit sphere
    let cos_phi = rng.random::<f32>() * 2.0 - 1.0;
    let sin_phi = sqrtf(1.0 - cos_phi * cos_phi);
    let theta = rng.random::<f32>() * std::f32::consts::TAU;

    // Construct the resulting vector.
    MVector::new(
        w * sin_phi * cosf(theta),
        w * sin_phi * sinf(theta),
        w * cos_phi,
        w,
    )
}

/// The maximum node-centric w-coordinate a horosphere can have such that the node in question
/// is still the owner of the horosphere.
// See `test_max_owned_horosphere_w()` for how this is computed.
const MAX_OWNED_HOROSPHERE_W: f32 = 5.9047837;

/// Represents a chunks's reference to a particular horosphere.
pub struct HorosphereChunk {
    /// The vector representing the horosphere in the perspective of the relevant chunk.
    /// See `HorosphereNode::pos` for details.
    pub pos: MVector<f32>,
}

impl HorosphereChunk {
    /// Creates a `HorosphereChunk` based on a `HorosphereNode`
    pub fn new(horosphere_node: &HorosphereNode, vertex: Vertex) -> Self {
        HorosphereChunk {
            pos: vertex.node_to_dual() * horosphere_node.pos,
        }
    }

    /// Rasterizes the horosphere chunk into the given `VoxelData`
    pub fn generate(&self, voxels: &mut VoxelData, chunk_size: u8) {
        for z in 0..chunk_size {
            for y in 0..chunk_size {
                for x in 0..chunk_size {
                    let pos = MVector::new(
                        x as f32 + 0.5,
                        y as f32 + 0.5,
                        z as f32 + 0.5,
                        chunk_size as f32 * Vertex::dual_to_chunk_factor(),
                    )
                    .normalized_point();
                    if pos.mip(&self.pos) > -1.0 {
                        voxels.data_mut(chunk_size)[Coords([x, y, z]).to_index(chunk_size)] =
                            Material::RedSandstone;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::math::MPoint;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_max_owned_horosphere_w() {
        // This tests that `MAX_OWNED_HOROSPHERE_W` is set to the correct value.

        // The worst case scenario would be a horosphere located directly in the direction of a dodeca's vertex.
        // This is because the horosphere can be outside the dodeca, tangent to each of the planes that extend the
        // dodeca's sides adjancent to that vertex. If that horosphere were brought any closer, it would intersect
        // all three of those planes, making it impossible for any child node to own the dodeca and forcing the node
        // in focus to own it.

        // First, find an arbitrary horosphere in the direction of a vertex.
        let example_vertex = Vertex::A;
        let example_vertex_pos = example_vertex.dual_to_node() * MPoint::origin();
        let mut horosphere_pos = MVector::from(example_vertex_pos);
        horosphere_pos.w = horosphere_pos.xyz().norm();

        // Then, scale the horosphere so that it's mip with each of the sides of the vertex is 1, making it tangent.
        horosphere_pos /= horosphere_pos.mip(example_vertex.canonical_sides()[0].normal());
        for side in example_vertex.canonical_sides() {
            assert_abs_diff_eq!(horosphere_pos.mip(side.normal()), 1.0, epsilon = 1.0e-6);
        }

        // Finally, compare that horosphere's w-coordinate to `MAX_OWNED_HOROSPHERE_W`
        assert_abs_diff_eq!(horosphere_pos.w, MAX_OWNED_HOROSPHERE_W, epsilon = 1.0e-6);
    }
}
