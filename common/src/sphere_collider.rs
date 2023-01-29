use crate::{
    dodeca::Vertex,
    math,
    ray_tracing::{ChunkRayTracer, RayStatus, RtChunkContext},
    world::Material,
};

pub struct SphereCollider {
    pub radius: f32,
}

impl ChunkRayTracer for SphereCollider {
    fn trace_ray(&self, ctx: &RtChunkContext, status: &mut RayStatus) {
        let float_size = ctx.dimension as f32;
        let voxel_start = (ctx.ray.column(0) / ctx.ray[(3, 0)]).xyz()
            * Vertex::dual_to_chunk_factor() as f32
            * float_size;
        let end_pos = ctx.ray * na::vector![1.0, status.tanh_length];
        let voxel_end =
            (end_pos / end_pos[3]).xyz() * Vertex::dual_to_chunk_factor() as f32 * float_size;
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
            self.find_edge_collision(ctx, &bbox, axis, status);
        }

        self.find_vertex_collision(ctx, &bbox, status);
    }

    fn max_radius(&self) -> f32 {
        self.radius
    }
}

impl SphereCollider {
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
                find_intersection_two_vectors(&ctx.ray, &edge_pos, &edge_dir, self.radius.cosh());

            // If t_candidate is out of range or NaN, don't continue collision checking
            if !(tanh_length_candidate >= 0.0 && tanh_length_candidate < status.tanh_length) {
                continue;
            }

            let translated_square_pos = ctx.ray * na::vector![1.0, tanh_length_candidate];
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

            status.update(ctx, tanh_length_candidate, translated_square_pos - projected_pos);
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
                find_intersection_one_vector(&ctx.ray, &vert, self.radius.cosh());

            // If t_candidate is out of range or NaN, don't continue collision checking
            if !(tanh_length_candidate >= 0.0 && tanh_length_candidate < status.tanh_length) {
                continue;
            }

            let translated_square_pos = ctx.ray * na::vector![1.0, tanh_length_candidate];
            status.update(ctx, tanh_length_candidate, translated_square_pos - vert);
        }
    }
}

/// Find the smallest value of `t` where the point in the pos-dir line (v=pos+dir*t) satisfies
/// `<v,a>^2 / <v,v> == c^2`
///
/// If `a` is direction-like, this finds intersections with a surface that is `sinh(c)` units
/// away from the plane whose normal is `a`.  If `a` is point-like, this finds intersections
/// with a sphere centered at `a` with radius `cosh(c)`.
///
/// Returns NaN if there's no such intersection
fn find_intersection_one_vector(ray: &na::Matrix4x2<f32>, a: &na::Vector4<f32>, c: f32) -> f32 {
    let mip_pos_a = math::mip(&ray.column(0), a);
    let mip_dir_a = math::mip(&ray.column(1), a);

    // The following 3 variables are terms of the quadratic formula. We use double the linear
    // term because it removes the annoying constants from that formula.
    let quadratic_term = mip_dir_a.powi(2) + c.powi(2);
    let double_linear_term = mip_pos_a * mip_dir_a;
    let constant_term = mip_pos_a.powi(2) - c.powi(2);

    if (-1e-4..=0.0).contains(&constant_term) && double_linear_term < 0.0 {
        return 0.0;
    }

    // If the player is already close to the wall, this function can produce incorrect results, so
    // ensure that we record a collision as long as the player is moving towards the wall.
    let discriminant = double_linear_term * double_linear_term - quadratic_term * constant_term;

    // While discriminant can be negative, NaNs propagate the way we want to, so we don't have
    // to check for this.
    constant_term / (-double_linear_term + discriminant.sqrt())
}

/// Find the smallest value of `t` where the point in the pos-dir line (v=pos+dir*t) satisfies
/// `(<v,a>^2 - <v,b>^2) / <v,v> == c^2`
///
/// This finds intersections with a surface that is `cosh(c)` units away from the line
/// with a point at `a` with direction `b`, where `<a,b>==0`.
///
/// Returns NaN if there's no such intersection
fn find_intersection_two_vectors(
    ray: &na::Matrix4x2<f32>,
    a: &na::Vector4<f32>,
    b: &na::Vector4<f32>,
    c: f32,
) -> f32 {
    // This could be made more numerically stable, but the precision requirements of collision should be pretty lax.
    let mip_pos_a = math::mip(&ray.column(0), a);
    let mip_dir_a = math::mip(&ray.column(1), a);
    let mip_pos_b = math::mip(&ray.column(0), b);
    let mip_dir_b = math::mip(&ray.column(1), b);

    // The following 3 variables are terms of the quadratic formula. We use double the linear
    // term because it removes the annoying constants from that formula.
    let quadratic_term = mip_dir_a.powi(2) - mip_dir_b.powi(2) + c.powi(2);
    let double_linear_term = mip_pos_a * mip_dir_a - mip_pos_b * mip_dir_b;
    let constant_term = mip_pos_a.powi(2) - mip_pos_b.powi(2) - c.powi(2);

    // If the player is already close to the wall, this function can produce incorrect results, so
    // ensure that we record a collision as long as the player is moving towards the wall.
    if (-1e-4..=0.0).contains(&constant_term) && double_linear_term < 0.0 {
        return 0.0;
    }

    let discriminant = double_linear_term * double_linear_term - quadratic_term * constant_term;

    // While discriminant can be negative, NaNs propagate the way we want to, so we don't have
    // to check for this.
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
