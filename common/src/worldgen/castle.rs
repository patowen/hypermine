use libm::{coshf, sqrtf};
use na::ComplexField;

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
            cylinder: StraightWallCylinder::new(),
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
    horizontal_tangent_vector: MDirection<f32>,
    center_plane: MDirection<f32>,
    tanh_radius: f32,
}

impl StraightWallCylinder {
    pub fn new() -> Self {
        let radius = 30.0;
        StraightWallCylinder {
            axis_point: Vertex::A.dual_to_node()
                * MIsometry::translation_along(&na::Vector3::new(0.0, radius, 0.0))
                * MPoint::origin(),
            axis_direction: Vertex::A.dual_to_node() * MDirection::x(),
            tangent_plane: Vertex::A.dual_to_node() * MDirection::y(),
            horizontal_tangent_vector: Vertex::A.dual_to_node() * MDirection::z(),
            center_plane: Vertex::A.dual_to_node() * MDirection::x(),
            tanh_radius: libm::tanhf(radius),
        }
    }

    pub fn contains_point(&self, point: &MPoint<f32>, debugging: bool) -> bool {
        /*
            - It's a cylinder in the Beltrami Klein projection, so we have (x/w)^2 + (y/w)^2 = k^2 for some k
            - Expand this out, and it becomes x^2 + y^2 - k^2*w^2 = 0. This equation is suitable when close to the center
            - Note that an alternative equation is 0 = x^2 + y^2 - k^2*w^2 = x^2 + y^2 - k^2*(x^2+y^2+z^2+1) = (1-k^2)x^2 + (1-k^2)y^2 - k^2z^2 - k^2 = 0
            - Lorentz-boost to edge of surface: x -> ax+bw, w = bx+aw with a^2 = b^2+1
                (ax+bw)^2 + y^2 - k^2*(bx+aw)^2 = 0
                a^2x^2 + 2abxw + b^2w^2 + y^2 - k^2b^2x^2 - 2k^2abxw - k^2a^2w^2 = 0
                    Want: Origin = (0,0,0,1) to be on shape. Plug in origin and get 0 = b^2 - k^2a^2 = a^2 - 1 - k^2a^2 = -1 + (1-k^2)a^2. So, a^2 = 1/(1-k^2) and b^2 = k^2/(1-k^2). So, ab = k/(1-k^2)
                x^2 + 2kxw + k^2w^2 + y^2(1-k^2) - k^4x^2 - 2k^3xw - k^2w^2 = 0
                (1-k^4)x^2 + 2k(1-k^2)xw + (1-k^2)y^2 = 0
                (1-k^2)(1+k^2)x^2 + 2k(1-k^2)xw + (1-k^2)y^2 = 0
                (1+k^2)x^2 + 2kxw + y^2 = 0
            - Lorentz-boost away from surface: z -> az+bw, w = bz+aw
                (1+k^2)x^2 + 2kx(bz+aw) + y^2 = 0
            - There doesn't seem to be any signs of numerical instability.
        */
        let x = -point.mip(&self.tangent_plane);
        let y = point.mip(&self.horizontal_tangent_vector);
        let z = point.mip(&self.center_plane);
        let w = libm::sqrtf(x * x + y * y + z * z + 1.0);
        let k = self.tanh_radius;

        // This logic belongs elsewhere, but I need to have some kind of demo available, so I'll put it here for now.
        let value = (1.0 + k * k) * x * x + 2.0 * k * x * w + y * y;
        let value_cutoff = -0.2 * libm::sqrtf(z * z + 1.0);

        if debugging {
            // Put println debugging here
        }

        let distance_to_axis = libm::acoshf(libm::sqrtf(
            sqr(point.mip(&self.axis_point)) - sqr(point.mip(&self.axis_direction)),
        ));
        return distance_to_axis.fract() < 0.1;

        if value > 0.0 {
            return false; // Outside of cylinder
        }
        // TODO: value_cutoff is not yet confirmed to be invariant to the particular perspective chosen. Do some more math.
        if value > value_cutoff {
            // Outer wall
            return true;
        }
        // Add a floor pattern for easier debugging
        if z > -0.05 && z < 0.0 {
            let distance_to_axis = libm::acoshf(libm::sqrtf(
                sqr(point.mip(&self.axis_point)) - sqr(point.mip(&self.axis_direction)),
            ));
            if distance_to_axis.fract() < 0.1 {
                return true;
            }
        }

        false
    }

    pub fn renormalize(&mut self) {
        if self.center_plane.w.abs() < 100.0 {
            /*// Project origin onto center plane
            let center_plane_point = (MVector::origin()
                + self.center_plane.as_ref() * self.center_plane.w)
                .normalized_point();
            /*
               Move c (axis_point) along a (axis_direction) to make it as close as possible to p (center_plane_point).
               This math is already worked out below. t = -<a,p> / <c,p>

               Note: <c+ta,c+ta> = <c,c> + 2t<a,c> + t^2<a,a> = -1 + t^2
            */
            let required_distance = -self.axis_direction.mip(&center_plane_point)
                / self.axis_point.mip(&center_plane_point);
            let mut grounded_axis_point =
                self.axis_point.as_ref() + self.axis_direction.as_ref() * required_distance;
            grounded_axis_point /= libm::sqrtf(1.0 - sqr(required_distance));
            grounded_axis_point.w = libm::sqrtf(1.0 + grounded_axis_point.xyz().norm_squared());
            /*projected_axis_point /= libm::sqrtf(
                1.0 + 2.0 * sqr(self.axis_point.mip(&self.center_plane))
                    - self.axis_point.mip(&self.center_plane),
            );*/
            //projected_axis_point.w = libm::sqrtf(1.0 + projected_axis_point.xyz().norm_squared());
            self.axis_point = grounded_axis_point.to_point_unchecked();
            self.axis_direction = self.center_plane;
            tracing::info!("{}, {}", self.axis_point.mip(&self.axis_direction), center_plane_point.mip(&self.center_plane));*/
            let center_plane_point = (MVector::origin()
                + self.center_plane.as_ref() * self.center_plane.w)
                .normalized_point();
            let expected_cosh_dist = libm::sqrtf(
                sqr(center_plane_point.mip(&self.axis_point))
                    - sqr(center_plane_point.mip(&self.axis_direction)),
            );
            let mut projected_axis_point = self.axis_point.as_ref()
                - self.center_plane.as_ref() * self.axis_point.mip(&self.center_plane);
            let actual_cosh_dist = -center_plane_point.mip(&projected_axis_point);
            tracing::info!("{}, {}", expected_cosh_dist, actual_cosh_dist);
            projected_axis_point *= expected_cosh_dist / actual_cosh_dist;
            let actual_cosh_dist = -center_plane_point.mip(&projected_axis_point);
            tracing::info!("{}, {}", expected_cosh_dist, actual_cosh_dist);
            self.axis_point = projected_axis_point.to_point_unchecked();
            self.axis_direction = self.center_plane;
        }

        /*
           Minimize w in terms of t for `(c+ta)/sqrt(-<c+ta, c+ta>)`
           = (c+ta)/sqrt(-[<c,c> + t<c,a> + t<a,c> + t^2<a,a>])
           = (c+ta)/sqrt(1 - t^2)
           if w(t) = (c.w + t*a.w) / sqrt(1 - t^2) + (c.w + t*a.w) * (1 - t^2)^(-1/2)
               Derivative of `(1 - t^2)^(-1/2)` is `t * (1 - t^2)^(-3/2)`
           then w'(t) = (c.w + t*a.w) * t * (1 - t^2)^(-3/2) + a.w * (1 - t^2)^(-1/2).
           w'(t)*(1 - t^2)^(3/2) = (c.w + t*a.w) * t + a.w * (1 - t^2)
               = t*c.w + t^2*a.w + a.w - t^2*a.w
               = a.w + t * c.w

           Set to 0: t = -a.w / c.w
        */
        let required_distance = -self.axis_direction.w / self.axis_point.w;
        let mut closest_axis_point = (self.axis_point.as_ref()
            + self.axis_direction.as_ref() * required_distance)
            / libm::sqrtf(1.0 - sqr(required_distance));
        closest_axis_point.w = libm::sqrtf(1.0 + closest_axis_point.xyz().norm_squared());
        self.axis_point = closest_axis_point.to_point_unchecked();

        // normalize axis_direction
        let mut origin_to_closest_axis_point = closest_axis_point;
        origin_to_closest_axis_point.w = 0.0;
        let origin_to_closest_axis_point = origin_to_closest_axis_point.normalized_direction(); // TODO: What if closest_axis_point is at the origin?
        let mut closest_axis_direction = *self.axis_direction.as_ref();
        closest_axis_direction -= origin_to_closest_axis_point.as_ref()
            * closest_axis_direction.mip(&origin_to_closest_axis_point);
        closest_axis_direction.w = 0.0; // Also project away origin
        self.axis_direction = closest_axis_direction.normalized_direction();

        let new_horizontal_tangent_vector = self
            .axis_direction
            .as_ref()
            .xyz()
            .cross(&self.axis_point.as_ref().xyz())
            .normalize();
        let new_horizontal_tangent_vector =
            MVector::from(new_horizontal_tangent_vector.to_homogeneous());
        // TODO: Better handle things near the origin
        let rotation = MIsometry::rotation(
            &self.horizontal_tangent_vector,
            &new_horizontal_tangent_vector.to_direction_unchecked(),
        );
        /*self.tangent_plane = (self.tangent_plane.as_ref()
        + (new_axis_to_tangent - self.axis_to_tangent.as_ref()) * axis_to_tangent_projection)
        .to_direction_unchecked();*/
        self.tangent_plane = rotation * self.tangent_plane;
        self.horizontal_tangent_vector = new_horizontal_tangent_vector.to_direction_unchecked();
    }

    pub fn set_axis(&mut self, axis: MDirection<f32>) {
        self.center_plane = axis;
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
            horizontal_tangent_vector: -(self * rhs.horizontal_tangent_vector), // TODO: Handle reflections better than this ad-hoc negation
            center_plane: self * rhs.center_plane,
            tanh_radius: rhs.tanh_radius,
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
