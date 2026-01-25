use libm::{coshf, sqrtf};

use crate::{
    dodeca::{Side, Vertex},
    graph::{Graph, NodeId},
    math::{self, MDirection, MIsometry, MPoint, MVector, sqr},
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
    pub fn new(graph: &Graph, node_id: NodeId, up: &MDirection<f32>) -> Option<CastleNode> {
        CastleNode::create_from_parents(graph, node_id, up)
            .or_else(|| CastleNode::maybe_create_fresh(graph, node_id, up))
    }

    fn create_from_parents(
        graph: &Graph,
        node_id: NodeId,
        up: &MDirection<f32>,
    ) -> Option<CastleNode> {
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
                        .and_then(|c| c.propagate(parent_side, up))
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

    fn maybe_create_fresh(
        graph: &Graph,
        node_id: NodeId,
        up: &MDirection<f32>,
    ) -> Option<CastleNode> {
        if node_id != NodeId::ROOT {
            return None;
        }

        Some(CastleNode {
            cylinder: StraightWallCylinder {
                axis_point: Vertex::A.dual_to_node()
                    * MIsometry::translation_along(&na::Vector3::new(0.0, 0.5, 0.0))
                    * MPoint::origin(),
                axis_direction: Vertex::A.dual_to_node() * MDirection::x(),
                tangent_plane: Vertex::A.dual_to_node() * MDirection::y(),
                axis_to_tangent: Vertex::A.dual_to_node()
                    * MIsometry::translation_along(&na::Vector3::new(0.0, 0.5, 0.0))
                    * -MDirection::y(),
            },
        })
    }

    fn propagate(&self, side: Side, up: &MDirection<f32>) -> Option<CastleNode> {
        // TODO: Don't propagate beyond the already-computed bounds of the `HorosphereNode`.

        let mut new_cylinder = side.reflection() * self.cylinder;
        new_cylinder.set_axis(*up);

        Some(CastleNode {
            cylinder: new_cylinder,
        })
    }

    fn average_with(&mut self, other: CastleNode, other_weight: f32) {
        // TODO: Add implementation (blank implementation works but puts full weight to `self`, making it numerically unstable)
    }

    pub fn log_player_stats(&self, transform: &MIsometry<f32>) {
        /*let position = transform * MVector::origin();
        let cylinder = self.cylinder;
        let cosh_horizontal_distance =
            -position.mip(&cylinder.center) / sqrtf(1.0 + math::sqr(position.mip(&cylinder.axis)));
        tracing::info!("cosh_horizontal_distance: {}", cosh_horizontal_distance);*/
        self.cylinder
            .contains_point(&(transform * MPoint::origin()), true);
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
                    if self.cylinder.contains_point(&pos, false) {
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
    axis_point: MPoint<f32>,
    axis_direction: MDirection<f32>,
    tangent_plane: MDirection<f32>,
    axis_to_tangent: MDirection<f32>, // Vector starting at axis_point and pointing in the direction of the tangent point
}

impl StraightWallCylinder {
    pub fn contains_point(&self, point: &MPoint<f32>, debugging: bool) -> bool {
        /*
        Let axis_point = c (for center)
        Let axis_direction = a (for axis)
        Let tangent_plane = n (for normal)
        Let point = p (for point)
        Constraints: <c,c> = -1, <a,a> = 1, <n,n> = 1, <p,p> = -1, <c,a> = 0
        First question: What reflection (Householder reflection?) can keep c and a fixed but move p to be spanned by c, a, and n?
            Call this reflection vector r. Contraints: <r,c> = 0, <r,a> = 0, <r,r> = 1
            Applying the reflection to p yields p-2r<p,r>

        New idea: Store a point on the surface and the normal at that point. To find the distance from `p` to the surface,
            - Slide the point of tangency along the vertical direction to be as close to `p` as possible.
            - Consider the plane at this new point of tangency perpendicular to this vertical direction. The intersection between this plane and the cylinder is a conic.

        New idea:
            - Project the point so that it's on the plane containing the axis and point of tangency.
            - Then, put the point back on the hyperboloid by moving it away from the axis.

        New idea (current):
            - Find the distance from the point to the axis
            - Find the plane the point lies on that is orthogonal to the axis
            - Find the point on that plane in a suitable direction relative to tangent_plane that is the same distance from the axis as the point
            - See whether that point is in front or behind the tangent_plane
        */

        let closest_axis_point = (self.axis_point.as_ref() * -point.mip(&self.axis_point)
            + self.axis_direction.as_ref() * point.mip(&self.axis_direction))
        .normalized_point();

        let cosh_closest_axis_point_distance = -point.mip(&closest_axis_point);

        let test_point = closest_axis_point.as_ref() * cosh_closest_axis_point_distance
            + self.axis_to_tangent.as_ref()
                * libm::sqrtf(sqr(cosh_closest_axis_point_distance) - 1.0);

        if debugging {
            // Nothing to debug yet
        }

        test_point.mip(&self.tangent_plane) > 0.0
    }

    pub fn renormalize(&mut self) {}

    pub fn set_axis(&mut self, axis: MDirection<f32>) {
        // TODO
    }
}

impl std::ops::Mul<StraightWallCylinder> for &MIsometry<f32> {
    type Output = StraightWallCylinder;

    fn mul(self, rhs: StraightWallCylinder) -> Self::Output {
        StraightWallCylinder {
            axis_point: self * rhs.axis_point,
            axis_direction: self * rhs.axis_direction,
            tangent_plane: self * rhs.tangent_plane,
            axis_to_tangent: self * rhs.axis_to_tangent,
        }
    }
}

#[derive(Copy, Clone)]
pub struct StraightWallCylinderOld {
    center: MPoint<f32>,
    axis: MDirection<f32>,
    center_radius: f32,
}

impl StraightWallCylinderOld {
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

        let cosh_horizontal_distance =
            -point.mip(&self.center) / sqrtf(1.0 + math::sqr(point.mip(&self.axis)));
        cosh_horizontal_distance < coshf(self.center_radius)
            && (self.center_radius + 0.01 - cosh_horizontal_distance.acosh()).fract() < 0.1
    }

    pub fn renormalize(&mut self) {
        let mut center = *self.center.as_ref();

        // Try to make sure center is orthogonal to axis, but add a bit of extra logic to avoid
        // catastrophic cancellation
        let estimated_mip = center.mip(&self.axis);
        let w_product = center.w * self.axis.w;
        let safe_mip =
            (estimated_mip.abs() - w_product.abs() * 1.0e-5).max(0.0) * estimated_mip.signum();
        center -= self.axis.as_ref() * safe_mip;

        // I want to normalize "center" to ensure that x^2+y^2+z^2-w^2 = -1. The traditional way
        // to do this is to scale the vector "center" by "sqrt(-mip(center, center))", but this risks
        // catastrophic cancellation. Instead, we only modify the w-coordinate.
        center.w = libm::sqrtf(1.0 + center.xyz().norm_squared());
        self.center = MPoint::new_unchecked(center.x, center.y, center.z, center.w);
    }

    pub fn set_axis(&mut self, axis: MDirection<f32>) {
        self.axis = axis;
    }
}

impl std::ops::Mul<StraightWallCylinderOld> for &MIsometry<f32> {
    type Output = StraightWallCylinderOld;

    fn mul(self, rhs: StraightWallCylinderOld) -> Self::Output {
        StraightWallCylinderOld {
            center: self * rhs.center,
            axis: self * rhs.axis,
            center_radius: rhs.center_radius,
        }
    }
}
