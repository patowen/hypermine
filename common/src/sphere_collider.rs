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
            self.find_side_collision(ctx, axis, status);
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
    fn find_side_collision(&self, ctx: &RtChunkContext, axis: usize, status: &mut RayStatus) {
        let float_dimension = ctx.dimension as f32;
        let plane0 = (axis + 1) % 3;
        let plane1 = (axis + 2) % 3;

        for i in ctx.bounding_box.voxel_plane_iterator(axis) {
            let normal = math::lorentz_normalize(&permuted_vector4(
                axis,
                1.0,
                0.0,
                0.0,
                i as f32 / float_dimension * Vertex::chunk_to_dual_factor() as f32,
            ));

            let tanh_distance =
                solve_sphere_plane_intersection(ctx.ray, &normal, self.radius.sinh());

            // If tanh_distance is out of range or NaN, no collision occurred.
            if !(tanh_distance >= 0.0 && tanh_distance < status.tanh_distance) {
                continue;
            }

            let mip_dir_norm = math::mip(&ctx.ray.direction, &normal);
            let i_with_offset = if mip_dir_norm < 0.0 { i } else { i + 1 };

            let collision_point = ctx.ray.point(tanh_distance);
            let projected_pos = collision_point - normal * math::mip(&collision_point, &normal);
            let projected_pos = projected_pos / projected_pos.w;
            let j0 =
                (projected_pos[plane0] * Vertex::dual_to_chunk_factor() as f32 * float_dimension)
                    .floor();
            let j1 =
                (projected_pos[plane1] * Vertex::dual_to_chunk_factor() as f32 * float_dimension)
                    .floor();

            if !(j0 >= 0.0 && j0 < float_dimension && j1 >= 0.0 && j1 < float_dimension) {
                continue;
            }
            let j0 = j0 as usize;
            let j1 = j1 as usize;

            if ctx.get_voxel(permuted_array3(axis, i_with_offset, j0 + 1, j1 + 1)) == Material::Void
            {
                continue;
            }

            status.update(ctx, tanh_distance, normal * -mip_dir_norm.signum());
        }
    }

    fn find_edge_collision(&self, ctx: &RtChunkContext, axis: usize, status: &mut RayStatus) {
        let float_dimension = ctx.dimension as f32;
        // TODO: Consider using alternate form of coordinates: t, u, and v, instead of axis, plane0, and plane1.
        let plane0 = (axis + 1) % 3;
        let plane1 = (axis + 2) % 3;

        for (i, j) in ctx.bounding_box.voxel_line_iterator(plane0, plane1) {
            let edge_pos = math::lorentz_normalize(&permuted_vector4(
                axis,
                0.0,
                i as f32 / float_dimension * Vertex::chunk_to_dual_factor() as f32,
                j as f32 / float_dimension * Vertex::chunk_to_dual_factor() as f32,
                1.0,
            ));
            let edge_dir = permuted_vector4(axis, 1.0, 0.0, 0.0, 0.0);

            let tanh_distance =
                solve_sphere_line_intersection(ctx.ray, &edge_pos, &edge_dir, self.radius.cosh());

            // If tanh_distance is out of range or NaN, no collision occurred.
            if !(tanh_distance >= 0.0 && tanh_distance < status.tanh_distance) {
                continue;
            }

            let collision_point = ctx.ray.point(tanh_distance);
            let projected_pos = -edge_pos * math::mip(&collision_point, &edge_pos)
                + edge_dir * math::mip(&collision_point, &edge_dir);
            let projected_pos = projected_pos / projected_pos.w;
            let k = (projected_pos[axis] * Vertex::dual_to_chunk_factor() as f32 * float_dimension)
                .floor();

            if !(k >= 0.0 && k < float_dimension) {
                continue;
            }
            let k = k as usize;

            if (0..2).all(|di| {
                (0..2).all(|dj| {
                    ctx.get_voxel(permuted_array3(axis, k + 1, i + di, j + dj)) == Material::Void
                })
            }) {
                continue;
            }

            status.update(ctx, tanh_distance, collision_point - projected_pos);
        }
    }

    fn find_vertex_collision(&self, ctx: &RtChunkContext, status: &mut RayStatus) {
        let float_dimension = ctx.dimension as f32;

        for (x, y, z) in ctx.bounding_box.voxel_iterator(0, 1, 2) {
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

            let collision_point = ctx.ray.point(tanh_distance);
            status.update(ctx, tanh_distance, collision_point - vertex_position);
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

/// Finds the lower solution `t` of `constant_term + 2 * double_linear_term * t + quadratic_term * t * t == 0`
///
/// Returns NaN if no such solution exists.
///
/// If a small perturbation to these terms would result in a solution of `t == 0.0`, this function has logic to
/// to return 0.0 if three conditions hold in the context of collision checking:
/// 1. The collider must be intersecting the object. This manifests as `constant_term <= 0.0`.
/// 2. The collider must not be too far inside the object. This manifests as `constant_term >= -EPSILON`.
/// 3. The direction of motion must be towards the collider. This manifests as `double_linear_term < 0.0`.
fn solve_quadratic(constant_term: f32, double_linear_term: f32, quadratic_term: f32) -> f32 {
    const EPSILON: f32 = 1e-4;

    // Extra logic to ensure precision issues don't allow a collider to clip through a surface
    if (-EPSILON..=0.0).contains(&constant_term) && double_linear_term < 0.0 {
        return 0.0;
    }

    let discriminant = double_linear_term * double_linear_term - quadratic_term * constant_term;

    // We use an alternative quadratic formula to ensure that we return a positive number if `constant_term > 0.0`.
    // Otherwise, the edge case of a small positive `constant_term` could be mishandled.
    // Note that discriminant can be negative, which allows this function to return NaN when there is no solution.
    constant_term / (-double_linear_term + discriminant.sqrt())
}

fn permuted_array3<N: Default>(axis: usize, axis_val: N, other_val0: N, other_val1: N) -> [N; 3] {
    let mut result: [N; 3] = Default::default();
    result[axis] = axis_val;
    result[(axis + 1) % 3] = other_val0;
    result[(axis + 2) % 3] = other_val1;
    result
}

fn permuted_vector4<N: na::RealField + Copy>(
    axis: usize,
    axis_coord: N,
    other_coord0: N,
    other_coord1: N,
    w: N,
) -> na::Vector4<N> {
    let mut result = na::Vector4::zeros();
    result[axis] = axis_coord;
    result[(axis + 1) % 3] = other_coord0;
    result[(axis + 2) % 3] = other_coord1;
    result[3] = w;
    result
}
