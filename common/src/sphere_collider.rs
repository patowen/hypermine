use crate::{
    dodeca::Vertex,
    math,
    ray_tracing::{ChunkRayTracer, Ray, RayStatus, RtChunkContext},
    world::Material,
};

pub struct SphereCollider {
    pub radius: f32,
}

impl ChunkRayTracer for SphereCollider {
    fn trace_ray(&self, ctx: &RtChunkContext, status: &mut RayStatus) {
        for axis in 0..3 {
            self.find_face_collision(ctx, axis, status);
        }

        for axis in 0..3 {
            self.find_edge_collision(ctx, axis, status);
        }

        self.find_vertex_collision(ctx, status);
    }

    fn max_radius(&self) -> f32 {
        self.radius
    }
}

impl SphereCollider {
    /// Detect collisions where a sphere contacts the front side of a voxel face
    fn find_face_collision(&self, ctx: &RtChunkContext, t_axis: usize, status: &mut RayStatus) {
        let u_axis = (t_axis + 1) % 3;
        let v_axis = (t_axis + 2) % 3;

        // Loop through all grid planes overlapping the bounding box
        for t in ctx.bounding_box.grid_plane_iterator(t_axis) {
            // Find a normal to the grid plane. Note that (t, 0, 0, x) is a normal of the plane whose closest point
            // to the origin is (x, 0, 0, t), and we use that fact here.
            let normal = math::lorentz_normalize(&tuv_to_xyz(
                t_axis,
                na::Vector4::new(1.0, 0.0, 0.0, grid_to_dual(ctx, t)),
            ));

            let tanh_distance =
                solve_sphere_plane_intersection(ctx.ray, &normal, self.radius.sinh());

            // If tanh_distance is out of range or NaN, no collision occurred.
            if !(tanh_distance >= 0.0 && tanh_distance < status.tanh_distance) {
                continue;
            }

            // Whether we are approaching the front or back of the face. An approach from the positive t direction
            // is 1, and an approach from the negative t direction is -1.
            let collision_side = -math::mip(&ctx.ray.direction, &normal).signum();

            // Which side we approach the plane from affects which voxel we want to use for collision checking
            let voxel_t = if collision_side > 0.0 { t } else { t + 1 };

            let ray_endpoint = ctx.ray.ray_point(tanh_distance);
            let contact_point = ray_endpoint - normal * math::mip(&ray_endpoint, &normal);

            // Compute the u and v-coordinates of the voxels at the contact point
            let Some(voxel_u) = dual_to_voxel(ctx, contact_point[u_axis] / contact_point.w) else {
                continue;
            };
            let Some(voxel_v) = dual_to_voxel(ctx, contact_point[v_axis] / contact_point.w) else {
                continue;
            };

            // Ensure that the relevant voxel is solid
            if ctx.get_voxel(tuv_to_xyz(t_axis, [voxel_t, voxel_u, voxel_v])) == Material::Void {
                continue;
            }

            // A collision was found. Update the status.
            status.update(ctx, tanh_distance, normal * collision_side);
        }
    }

    /// Detect collisions where a sphere contacts a voxel edge
    fn find_edge_collision(&self, ctx: &RtChunkContext, t_axis: usize, status: &mut RayStatus) {
        let u_axis = (t_axis + 1) % 3;
        let v_axis = (t_axis + 2) % 3;

        // Loop through all grid lines overlapping the bounding box
        for (u, v) in ctx.bounding_box.grid_line_iterator(u_axis, v_axis) {
            let edge_pos = math::lorentz_normalize(&tuv_to_xyz(
                t_axis,
                na::Vector4::new(0.0, grid_to_dual(ctx, u), grid_to_dual(ctx, v), 1.0),
            ));
            let edge_dir = tuv_to_xyz(t_axis, na::Vector4::new(1.0, 0.0, 0.0, 0.0));

            let tanh_distance =
                solve_sphere_line_intersection(ctx.ray, &edge_pos, &edge_dir, self.radius.cosh());

            // If tanh_distance is out of range or NaN, no collision occurred.
            if !(tanh_distance >= 0.0 && tanh_distance < status.tanh_distance) {
                continue;
            }

            let ray_endpoint = ctx.ray.ray_point(tanh_distance);
            let contact_point = -edge_pos * math::mip(&ray_endpoint, &edge_pos)
                + edge_dir * math::mip(&ray_endpoint, &edge_dir);

            // Compute the t-coordinate of the voxels at the contact point
            let Some(voxel_t) = dual_to_voxel(ctx, contact_point[t_axis] / contact_point.w) else {
                continue;
            };

            // Ensure that the edge has a solid voxel adjacent to it
            if (0..2).all(|du| {
                (0..2).all(|dv| {
                    ctx.get_voxel(tuv_to_xyz(t_axis, [voxel_t, u + du, v + dv])) == Material::Void
                })
            }) {
                continue;
            }

            // A collision was found. Update the status.
            status.update(ctx, tanh_distance, ray_endpoint - contact_point);
        }
    }

    /// Detect collisions where a sphere contacts a voxel vertex
    fn find_vertex_collision(&self, ctx: &RtChunkContext, status: &mut RayStatus) {
        let float_dimension = ctx.dimension as f32;

        // Loop through all grid points contained in the bounding box
        for (x, y, z) in ctx.bounding_box.grid_point_iterator(0, 1, 2) {
            // Skip vertices that have no solid voxels adjacent to them
            if (0..2).all(|dx| {
                (0..2).all(|dy| {
                    (0..2).all(|dz| ctx.get_voxel([x + dx, y + dy, z + dz]) == Material::Void)
                })
            }) {
                continue;
            }

            // Determine the cube-centric coordinates of the vertex
            let vertex_position = math::lorentz_normalize(
                &na::Vector3::new(x as f32, y as f32, z as f32)
                    .scale(Vertex::chunk_to_dual_factor() as f32 / float_dimension)
                    .insert_row(3, 1.0),
            );

            let tanh_distance =
                solve_sphere_point_intersection(ctx.ray, &vertex_position, self.radius.cosh());

            // If tanh_distance is out of range or NaN, no collision occurred.
            if !(tanh_distance >= 0.0 && tanh_distance < status.tanh_distance) {
                continue;
            }

            // A collision was found. Update the status.
            let ray_endpoint = ctx.ray.ray_point(tanh_distance);
            status.update(ctx, tanh_distance, ray_endpoint - vertex_position);
        }
    }
}

/// Finds the tanh of the distance a sphere will have to travel along a ray before it
/// intersects the given point.
///
/// Returns NaN if there's no such intersection
fn solve_sphere_point_intersection(
    ray: &Ray,
    point_position: &na::Vector4<f32>,
    cosh_radius: f32,
) -> f32 {
    // This could be made more numerically stable by using a formula that depends on sinh_radius,
    // but the precision requirements of collision should be pretty lax.
    let mip_pos_a = math::mip(&ray.position, point_position);
    let mip_dir_a = math::mip(&ray.direction, point_position);

    solve_quadratic(
        mip_pos_a.powi(2) - cosh_radius.powi(2),
        mip_pos_a * mip_dir_a,
        mip_dir_a.powi(2) + cosh_radius.powi(2),
    )
}

/// Finds the tanh of the distance a sphere will have to travel along a ray before it
/// intersects the given line.
///
/// Returns NaN if there's no such intersection
fn solve_sphere_line_intersection(
    ray: &Ray,
    line_position: &na::Vector4<f32>,
    line_direction: &na::Vector4<f32>,
    cosh_radius: f32,
) -> f32 {
    // This could be made more numerically stable by using a formula that depends on sinh_radius,
    // but the precision requirements of collision should be pretty lax.
    let mip_pos_a = math::mip(&ray.position, line_position);
    let mip_dir_a = math::mip(&ray.direction, line_position);
    let mip_pos_b = math::mip(&ray.position, line_direction);
    let mip_dir_b = math::mip(&ray.direction, line_direction);

    solve_quadratic(
        mip_pos_a.powi(2) - mip_pos_b.powi(2) - cosh_radius.powi(2),
        mip_pos_a * mip_dir_a - mip_pos_b * mip_dir_b,
        mip_dir_a.powi(2) - mip_dir_b.powi(2) + cosh_radius.powi(2),
    )
}

/// Finds the tanh of the distance a sphere will have to travel along a ray before it
/// intersects the given plane.
///
/// Returns NaN if there's no such intersection
fn solve_sphere_plane_intersection(
    ray: &Ray,
    plane_normal: &na::Vector4<f32>,
    sinh_radius: f32,
) -> f32 {
    let mip_pos_a = math::mip(&ray.position, plane_normal);
    let mip_dir_a = math::mip(&ray.direction, plane_normal);

    solve_quadratic(
        mip_pos_a.powi(2) - sinh_radius.powi(2),
        mip_pos_a * mip_dir_a,
        mip_dir_a.powi(2) + sinh_radius.powi(2),
    )
}

/// Finds the lower solution `x` of `constant_term + 2 * half_linear_term * x + quadratic_term * x * x == 0`
///
/// Returns NaN if no such solution exists.
///
/// If a small perturbation to these terms would result in a solution of `x == 0.0`, this function has logic to
/// to return 0.0 if three conditions hold in the context of collision checking:
/// 1. The collider must be intersecting the object. This manifests as `constant_term <= 0.0`.
/// 2. The collider must not be too far inside the object. This manifests as `constant_term >= -EPSILON`.
/// 3. The direction of motion must be towards the collider. This manifests as `double_linear_term < 0.0`.
fn solve_quadratic(constant_term: f32, half_linear_term: f32, quadratic_term: f32) -> f32 {
    const EPSILON: f32 = 1e-4;

    // Extra logic to ensure precision issues don't allow a collider to clip through a surface
    if (-EPSILON..=0.0).contains(&constant_term) && half_linear_term < 0.0 {
        return 0.0;
    }

    let discriminant = half_linear_term * half_linear_term - quadratic_term * constant_term;

    // We use an alternative quadratic formula to ensure that we return a positive number if `constant_term > 0.0`.
    // Otherwise, the edge case of a small positive `constant_term` could be mishandled.
    // Note that discriminant can be negative, which allows this function to return NaN when there is no solution.
    constant_term / (-half_linear_term + discriminant.sqrt())
}

/// Converts from t-u-v coordinates to x-y-z coordinates. t-u-v coordinates are a permuted version of x-y-z coordinates.
/// `t_axis` determines which of the three x-y-z coordinates corresponds to the t-coordinate. This function works with
/// any indexable entity with at least three entries. Any entry after the third entry is ignored.
fn tuv_to_xyz<T: std::ops::IndexMut<usize, Output = N>, N: Copy>(t_axis: usize, tuv: T) -> T {
    let mut result = tuv;
    (
        result[t_axis],
        result[(t_axis + 1) % 3],
        result[(t_axis + 2) % 3],
    ) = (result[0], result[1], result[2]);
    result
}

/// Converts a single coordinate from dual coordinates in the Klein-Beltrami model to an integer coordinate
/// suitable for voxel lookup. Margins are included. Returns `None` if the coordinate is outside the chunk.
#[inline]
fn dual_to_voxel(ctx: &RtChunkContext, dual_coord: f32) -> Option<usize> {
    let voxel_coord =
        (dual_coord * Vertex::dual_to_chunk_factor() as f32 * ctx.dimension_f32).floor();

    if !(voxel_coord >= 0.0 && voxel_coord < ctx.dimension_f32) {
        None
    } else {
        Some(voxel_coord as usize + 1)
    }
}

#[inline]
fn grid_to_dual(ctx: &RtChunkContext, grid_coord: usize) -> f32 {
    grid_coord as f32 / ctx.dimension_f32 * Vertex::chunk_to_dual_factor() as f32
}