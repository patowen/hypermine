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

        self.find_vertex_collision(ctx, &bbox, status);
    }

    fn max_radius(&self) -> f32 {
        self.radius
    }
}

impl SphereCollider {
    fn find_vertex_collision(
        &self,
        ctx: &RtChunkContext,
        bbox: &[[usize; 2]; 3],
        status: &mut RayStatus,
    ) {
        let float_dimension = ctx.dimension as f32;

        for i in bbox[0][0]..bbox[0][1] {
            for j in bbox[1][0]..bbox[1][1] {
                for k in bbox[2][0]..bbox[2][1] {
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
                    if !(tanh_length_candidate >= 0.0 && tanh_length_candidate < status.tanh_length)
                    {
                        continue;
                    }

                    let translated_square_pos = ctx.ray * na::vector![1.0, tanh_length_candidate];
                    status.update(ctx, tanh_length_candidate, translated_square_pos - vert);
                }
            }
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
