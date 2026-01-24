use libm::{coshf, sqrtf};

use crate::{
    dodeca::{Side, Vertex},
    graph::{Graph, NodeId},
    math::{self, MDirection, MIsometry, MPoint, MVector},
    node::VoxelData,
    voxel_math::Coords,
    world::Material,
};

const CASTLE_SEED: u64 = 934820348209; // TODO

#[derive(Copy, Clone)]
pub struct CastleNode {
    cylinder: StraightWallCylinder,
}

impl CastleNode {
    pub fn new(graph: &Graph, node_id: NodeId) -> Option<CastleNode> {
        CastleNode::create_from_parents(graph, node_id)
            .or_else(|| CastleNode::maybe_create_fresh(graph, node_id))
    }

    fn create_from_parents(graph: &Graph, node_id: NodeId) -> Option<CastleNode> {
        // Rather than selecting an arbitrary parent CastleNode, we average all of them. This
        // is important because otherwise, the propagation of floating point precision errors could
        // create a seam. This ensures that all errors average out, keeping the horosphere smooth.
        let mut cylinders_to_average_iter =
            graph
                .parents(node_id)
                .filter_map(|(parent_side, parent_id)| {
                    graph
                        .node_state(parent_id)
                        .castle
                        .as_ref()
                        .and_then(|c| c.propagate(parent_side))
                });

        let mut castle_node = cylinders_to_average_iter.next()?;
        let mut count = 1;
        for other in cylinders_to_average_iter {
            // Take an average of all HorosphereNodes in this iterator, giving each of them equal weight
            // by keeping track of a moving average with a weight that changes over time to make the
            // numbers work out the same way.
            count += 1;
            castle_node.average_with(other, 1.0 / count as f32);
        }

        castle_node.cylinder.renormalize();
        // TODO: castle_node.tighten_region_bounds();
        Some(castle_node)
    }

    fn maybe_create_fresh(graph: &Graph, node_id: NodeId) -> Option<CastleNode> {
        if node_id != NodeId::ROOT {
            return None;
        }

        Some(CastleNode {
            cylinder: StraightWallCylinder {
                center: Vertex::A.dual_to_node()
                    * MIsometry::translation_along(&na::Vector3::new(0.0, 30.0 + 1.5, 0.0))
                    * MPoint::origin(),
                axis: *Side::A.normal(),
                center_radius: 30.0,
            },
        })
    }

    fn propagate(&self, side: Side) -> Option<CastleNode> {
        // TODO: Don't propagate beyond the already-computed bounds of the `HorosphereNode`.

        Some(CastleNode {
            cylinder: side.reflection() * self.cylinder,
        })
    }

    fn average_with(&mut self, other: CastleNode, other_weight: f32) {
        // TODO: Add implementation (blank implementation works but puts full weight to `self`, making it numerically unstable)
    }
}

pub struct CastleChunk {
    pub cylinder: StraightWallCylinder,
}

impl CastleChunk {
    /// Creates a `HorosphereChunk` based on a `HorosphereNode`
    pub fn new(castle_node: &CastleNode, vertex: Vertex) -> Self {
        CastleChunk {
            cylinder: vertex.node_to_dual() * castle_node.cylinder,
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
                    if self.cylinder.contains_point(&pos) {
                        voxels.data_mut(chunk_size)[Coords([x, y, z]).to_index(chunk_size)] =
                            Material::RedSandstone;
                    }
                }
            }
        }
    }
}

#[derive(Copy, Clone)]
pub struct StraightWallCylinder {
    center: MPoint<f32>,
    axis: MDirection<f32>,
    center_radius: f32,
}

impl StraightWallCylinder {
    pub fn contains_point(&self, point: &MPoint<f32>) -> bool {
        /*
        projected_point =
        (p - a*<p,a>) / sqrt(-<p-a*<p,a>, p-a*<p,a>>)
        (p - a*<p,a>) / sqrt(- [<p,p> + (<a,a> - 2)*<p,a>*<p,a>])
        (p - a*<p,a>) / sqrt(- [-1 - <p,a>*<p,a>])
        (p - a*<p,a>) / sqrt(1 + <p,a>*<p,a>)

        -projected_point mip center =
        (-<p,c> + <a,c>*<p,a>) / sqrt(1 + <p,a>*<p,a>)
        (-<p,c> + 0*<p,a>) / sqrt(1 + <p,a>*<p,a>)
        -<p,c> / sqrt(1 + <p,a>*<p,a>)
         */

        let cosh_horizontal_distance = -point.mip(&self.center) / sqrtf(1.0 + math::sqr(point.mip(&self.axis)));
        cosh_horizontal_distance < coshf(self.center_radius) && cosh_horizontal_distance.acosh().fract() < 0.1
    }

    pub fn renormalize(&mut self) {
        //self.center = self.center.as_ref().normalized_point();
        //self.axis = self.axis.as_ref().normalized_direction();
    }
}

impl std::ops::Mul<StraightWallCylinder> for &MIsometry<f32> {
    type Output = StraightWallCylinder;

    fn mul(self, rhs: StraightWallCylinder) -> Self::Output {
        StraightWallCylinder {
            center: self * rhs.center,
            axis: self * rhs.axis,
            center_radius: rhs.center_radius,
        }
    }
}
