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
        let float_size = ctx.dimension as f32;
        let voxel_start = na::Point3::from_homogeneous(ctx.ray.position).unwrap()
            * Vertex::dual_to_chunk_factor() as f32
            * float_size;
        let voxel_end = na::Point3::from_homogeneous(ctx.ray.point(status.tanh_length)).unwrap()
            * Vertex::dual_to_chunk_factor() as f32
            * float_size;
        let max_voxel_radius = self.radius * Vertex::dual_to_chunk_factor() as f32 * float_size;
        let bbox = [0, 1, 2].map(|coord| {
            get_usize_range(
                ctx.dimension,
                voxel_start[coord],
                voxel_end[coord],
                max_voxel_radius,
            )
        });

        for axis in 0..3 {
            self.find_side_collision(ctx, &bbox, axis, status);
        }

        for axis in 0..3 {
            self.find_edge_collision(ctx, &bbox, axis, status);
        }

        self.find_vertex_collision(ctx, &bbox, status);
    }

    fn max_radius(&self) -> f32 {
        self.radius
    }
}

impl SphereCollider {
    fn find_side_collision(
        &self,
        ctx: &RtChunkContext,
        bbox: &[[usize; 2]; 3],
        axis: usize,
        status: &mut RayStatus,
    ) {
        let float_dimension = ctx.dimension as f32;
        let plane0 = (axis + 1) % 3;
        let plane1 = (axis + 2) % 3;

        for i in bbox[axis][0]..bbox[axis][1] {
            let normal = math::lorentz_normalize(&permuted_vector4(
                axis,
                1.0,
                0.0,
                0.0,
                i as f32 / float_dimension * Vertex::chunk_to_dual_factor() as f32,
            ));

            let tanh_length_candidate =
                solve_sphere_plane_intersection(ctx.ray, &normal, self.radius.sinh());

            // If t_candidate is out of range or NaN, don't continue collision checking
            if !(tanh_length_candidate >= 0.0 && tanh_length_candidate < status.tanh_length) {
                continue;
            }

            let mip_dir_norm = math::mip(&ctx.ray.direction, &normal);
            let i_with_offset = if mip_dir_norm < 0.0 { i } else { i + 1 };

            let translated_square_pos = ctx.ray.point(tanh_length_candidate);
            let projected_pos =
                translated_square_pos - normal * math::mip(&translated_square_pos, &normal);
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

            status.update(ctx, tanh_length_candidate, normal * -mip_dir_norm.signum());
        }
    }

    fn find_edge_collision(
        &self,
        ctx: &RtChunkContext,
        bbox: &[[usize; 2]; 3],
        axis: usize,
        status: &mut RayStatus,
    ) {
        let float_dimension = ctx.dimension as f32;
        let plane0 = (axis + 1) % 3;
        let plane1 = (axis + 2) % 3;

        for (i, j) in (bbox[plane0][0]..bbox[plane0][1])
            .flat_map(|i| (bbox[plane1][0]..bbox[plane1][1]).map(move |j| (i, j)))
        {
            let edge_pos = math::lorentz_normalize(&permuted_vector4(
                axis,
                0.0,
                i as f32 / float_dimension * Vertex::chunk_to_dual_factor() as f32,
                j as f32 / float_dimension * Vertex::chunk_to_dual_factor() as f32,
                1.0,
            ));
            let edge_dir = permuted_vector4(axis, 1.0, 0.0, 0.0, 0.0);

            let tanh_length_candidate =
                solve_sphere_line_intersection(ctx.ray, &edge_pos, &edge_dir, self.radius.cosh());

            // If t_candidate is out of range or NaN, don't continue collision checking
            if !(tanh_length_candidate >= 0.0 && tanh_length_candidate < status.tanh_length) {
                continue;
            }

            let translated_square_pos = ctx.ray.point(tanh_length_candidate);
            let projected_pos = -edge_pos * math::mip(&translated_square_pos, &edge_pos)
                + edge_dir * math::mip(&translated_square_pos, &edge_dir);
            let projected_pos = projected_pos / projected_pos.w;
            let k = (projected_pos[axis] * Vertex::dual_to_chunk_factor() as f32 * float_dimension)
                .floor();

            if !(k >= 0.0 && k < float_dimension) {
                continue;
            }
            let k = k as usize;

            if (0..4).all(|idx| {
                ctx.get_voxel(permuted_array3(
                    axis,
                    k + 1,
                    i + idx % 2,
                    j + (idx >> 1) % 2,
                )) == Material::Void
            }) {
                continue;
            }

            status.update(
                ctx,
                tanh_length_candidate,
                translated_square_pos - projected_pos,
            );
        }
    }

    fn find_vertex_collision(
        &self,
        ctx: &RtChunkContext,
        bbox: &[[usize; 2]; 3],
        status: &mut RayStatus,
    ) {
        let float_dimension = ctx.dimension as f32;

        for (i, j, k) in (bbox[0][0]..bbox[0][1]).flat_map(|i| {
            (bbox[1][0]..bbox[1][1])
                .flat_map(move |j| (bbox[2][0]..bbox[2][1]).map(move |k| (i, j, k)))
        }) {
            if (0..8).all(|idx| {
                ctx.get_voxel([i + idx % 2, j + (idx >> 1) % 2, k + (idx >> 2) % 2])
                    == Material::Void
            }) {
                continue;
            }

            let vert = math::lorentz_normalize(
                &na::Vector3::new(i as f32, j as f32, k as f32)
                    .scale(Vertex::chunk_to_dual_factor() as f32 / float_dimension)
                    .insert_row(3, 1.0),
            );

            let tanh_length_candidate =
                solve_sphere_point_intersection(ctx.ray, &vert, self.radius.cosh());

            // If t_candidate is out of range or NaN, don't continue collision checking
            if !(tanh_length_candidate >= 0.0 && tanh_length_candidate < status.tanh_length) {
                continue;
            }

            let translated_square_pos = ctx.ray.point(tanh_length_candidate);
            status.update(ctx, tanh_length_candidate, translated_square_pos - vert);
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

fn get_usize_range(dimension: usize, point0: f32, point1: f32, width: f32) -> [usize; 2] {
    if !point0.is_finite() || !point1.is_finite() {
        return [0, dimension + 1];
    }
    let result_min = (point0.min(point1) - width).max(-0.5) + 1.0;
    let result_max = (point0.max(point1) + width).min(dimension as f32 + 0.5) + 1.0;

    if result_min > result_max {
        // Empty range
        return [1, 0];
    }

    [result_min.floor() as usize, result_max.floor() as usize]
}
