use crate::{
    collision::{ChunkShapeCastingContext, Ray, RayEndpoint},
    math,
    world::Material,
};

/// Handles collisions for a single chunk.
pub fn chunk_shape_cast(ctx: &ChunkShapeCastingContext, endpoint: &mut RayEndpoint) {
    for axis in 0..3 {
        find_face_collision(ctx, axis, endpoint);
    }

    for axis in 0..3 {
        find_edge_collision(ctx, axis, endpoint);
    }

    find_vertex_collision(ctx, endpoint);
}

/// Detect collisions where a sphere contacts the front side of a voxel face
fn find_face_collision(ctx: &ChunkShapeCastingContext, t_axis: usize, endpoint: &mut RayEndpoint) {
    let u_axis = (t_axis + 1) % 3;
    let v_axis = (t_axis + 2) % 3;

    // Loop through all grid planes overlapping the bounding box
    for t in ctx.bounding_box.grid_planes(t_axis) {
        // Find a normal to the grid plane. Note that (t, 0, 0, x) is a normal of the plane whose closest point
        // to the origin is (x, 0, 0, t), and we use that fact here.
        let normal = math::lorentz_normalize(&tuv_to_xyz(
            t_axis,
            na::Vector4::new(1.0, 0.0, 0.0, grid_to_dual(ctx, t)),
        ));

        let Some(tanh_distance) =
            solve_sphere_plane_intersection(ctx.ray, &normal, ctx.collider_radius.sinh()) else {
                continue;
            };

        // If tanh_distance is out of range or NaN, no collision occurred.
        if tanh_distance < 0.0 || tanh_distance >= endpoint.tanh_distance {
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

        // A collision was found. Update the endpoint.
        endpoint.update(ctx, tanh_distance, normal * collision_side);
    }
}

/// Detect collisions where a sphere contacts a voxel edge
fn find_edge_collision(ctx: &ChunkShapeCastingContext, t_axis: usize, endpoint: &mut RayEndpoint) {
    let u_axis = (t_axis + 1) % 3;
    let v_axis = (t_axis + 2) % 3;

    // Loop through all grid lines overlapping the bounding box
    for (u, v) in ctx.bounding_box.grid_lines(u_axis, v_axis) {
        let edge_pos = math::lorentz_normalize(&tuv_to_xyz(
            t_axis,
            na::Vector4::new(0.0, grid_to_dual(ctx, u), grid_to_dual(ctx, v), 1.0),
        ));
        let edge_dir = tuv_to_xyz(t_axis, na::Vector4::new(1.0, 0.0, 0.0, 0.0));

        let Some(tanh_distance) = solve_sphere_line_intersection(
            ctx.ray,
            &edge_pos,
            &edge_dir,
            ctx.collider_radius.cosh(),
        ) else {
            continue;
        };

        // If tanh_distance is out of range, no collision occurred.
        if tanh_distance < 0.0 || tanh_distance >= endpoint.tanh_distance {
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

        // A collision was found. Update the endpoint.
        endpoint.update(ctx, tanh_distance, ray_endpoint - contact_point);
    }
}

/// Detect collisions where a sphere contacts a voxel vertex
fn find_vertex_collision(ctx: &ChunkShapeCastingContext, endpoint: &mut RayEndpoint) {
    // Loop through all grid points contained in the bounding box
    for (x, y, z) in ctx.bounding_box.grid_points(0, 1, 2) {
        // Skip vertices that have no solid voxels adjacent to them
        if (0..2).all(|dx| {
            (0..2).all(|dy| {
                (0..2).all(|dz| ctx.get_voxel([x + dx, y + dy, z + dz]) == Material::Void)
            })
        }) {
            continue;
        }

        // Determine the cube-centric coordinates of the vertex
        let vertex_position = math::lorentz_normalize(&na::Vector4::new(
            grid_to_dual(ctx, x),
            grid_to_dual(ctx, y),
            grid_to_dual(ctx, z),
            1.0,
        ));

        let Some(tanh_distance) =
            solve_sphere_point_intersection(ctx.ray, &vertex_position, ctx.collider_radius.cosh()) else {
                continue;
            };

        // If tanh_distance is out of range or NaN, no collision occurred.
        if tanh_distance < 0.0 || tanh_distance >= endpoint.tanh_distance {
            continue;
        }

        // A collision was found. Update the endpoint.
        let ray_endpoint = ctx.ray.ray_point(tanh_distance);
        endpoint.update(ctx, tanh_distance, ray_endpoint - vertex_position);
    }
}

/// Finds the tanh of the distance a sphere will have to travel along a ray before it
/// intersects the given plane.
fn solve_sphere_plane_intersection(
    ray: &Ray,
    plane_normal: &na::Vector4<f32>,
    sinh_radius: f32,
) -> Option<f32> {
    let mip_pos_a = math::mip(&ray.position, plane_normal);
    let mip_dir_a = math::mip(&ray.direction, plane_normal);

    solve_quadratic(
        mip_pos_a.powi(2) - sinh_radius.powi(2),
        mip_pos_a * mip_dir_a,
        mip_dir_a.powi(2) + sinh_radius.powi(2),
    )
}

/// Finds the tanh of the distance a sphere will have to travel along a ray before it
/// intersects the given line.
fn solve_sphere_line_intersection(
    ray: &Ray,
    line_position: &na::Vector4<f32>,
    line_direction: &na::Vector4<f32>,
    cosh_radius: f32,
) -> Option<f32> {
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
/// intersects the given point.
fn solve_sphere_point_intersection(
    ray: &Ray,
    point_position: &na::Vector4<f32>,
    cosh_radius: f32,
) -> Option<f32> {
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

/// Finds the lower solution `x` of `constant_term + 2 * half_linear_term * x + quadratic_term * x * x == 0`
/// if such a solution exists. Assumes that `quadratic_term` is positive.
///
/// If a small perturbation to these terms would result in a solution of `x == 0.0`, this function has logic to
/// to return Some(0.0) if three conditions hold in the context of collision checking:
/// 1. The collider must be intersecting the object. This manifests as `constant_term <= 0.0`.
/// 2. The collider must not be too far inside the object. This manifests as `constant_term >= -EPSILON`.
/// 3. The direction of motion must be towards the collider. This manifests as `double_linear_term < 0.0`.
fn solve_quadratic(constant_term: f32, half_linear_term: f32, quadratic_term: f32) -> Option<f32> {
    const EPSILON: f32 = 1e-4;

    // Extra logic to ensure precision issues don't allow a collider to clip through a surface
    if (-EPSILON..=0.0).contains(&constant_term) && half_linear_term < 0.0 {
        return Some(0.0);
    }

    let discriminant = half_linear_term * half_linear_term - quadratic_term * constant_term;
    if discriminant < 0.0 {
        return None;
    }

    // We use an alternative quadratic formula to ensure that we return a positive number if `constant_term > 0.0`.
    // Otherwise, the edge case of a small positive `constant_term` could be mishandled.
    // Note that discriminant can be negative, which allows this function to return NaN when there is no solution.
    Some(constant_term / (-half_linear_term + discriminant.sqrt()))
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
fn dual_to_voxel(ctx: &ChunkShapeCastingContext, dual_coord: f32) -> Option<usize> {
    let floor_grid_coord = (dual_coord * ctx.dual_to_grid_factor).floor();

    if !(floor_grid_coord >= 0.0 && floor_grid_coord < ctx.dimension as f32) {
        None
    } else {
        Some(floor_grid_coord as usize + 1)
    }
}

/// Converts a single coordinate from grid coordinates to dual coordiantes in the Klein-Beltrami model. This
/// can be used to find the positions of voxel gridlines.
#[inline]
fn grid_to_dual(ctx: &ChunkShapeCastingContext, grid_coord: usize) -> f32 {
    grid_coord as f32 / ctx.dual_to_grid_factor
}

#[cfg(test)]
mod tests {
    use crate::{
        collision::VoxelAABB,
        dodeca::Vertex,
        graph::NodeId,
        node::{ChunkId, VoxelData},
    };

    use super::*;
    use approx::*;

    struct TestShapeCastContext {
        collider_radius: f32,
        dimension: usize,
        dual_to_grid_factor: f32,
        voxel_data: VoxelData,
        transform: na::Matrix4<f32>,
    }

    impl TestShapeCastContext {
        fn new(collider_radius: f32) -> Self {
            // Use an arbitrary transformation
            let transform = na::Rotation3::from_euler_angles(0.1, 0.2, 0.3).to_homogeneous()
                * math::translate_along(&na::Vector3::new(0.5, 0.3, 0.2));
            let dimension: usize = 12;

            let mut test_ctx = TestShapeCastContext {
                collider_radius,
                dimension,
                dual_to_grid_factor: Vertex::dual_to_chunk_factor() as f32 * dimension as f32,
                voxel_data: VoxelData::Solid(Material::Void),
                transform,
            };

            // Populate voxels. Consists of a single cube with grid coordinates from (1, 1, 1) to (2, 2, 2)
            test_ctx.set_voxel([2, 2, 2], Material::Dirt);

            test_ctx
        }

        fn set_voxel(&mut self, coords: [usize; 3], material: Material) {
            let dimension_with_margin = self.dimension + 2;
            debug_assert!(coords[0] < dimension_with_margin);
            debug_assert!(coords[1] < dimension_with_margin);
            debug_assert!(coords[2] < dimension_with_margin);
            self.voxel_data.data_mut(self.dimension as u8)[coords[0]
                + coords[1] * dimension_with_margin
                + coords[2] * dimension_with_margin.pow(2)] = material;
        }
    }

    fn cast_with_test_ray(
        test_ctx: &TestShapeCastContext,
        ray_start_grid_coords: [f32; 3],
        ray_end_grid_coords: [f32; 3],
        wrapped_fn: impl FnOnce(&ChunkShapeCastingContext, f32),
    ) {
        let ray_start = math::lorentz_normalize(&na::Vector4::new(
            ray_start_grid_coords[0] / test_ctx.dual_to_grid_factor,
            ray_start_grid_coords[1] / test_ctx.dual_to_grid_factor,
            ray_start_grid_coords[2] / test_ctx.dual_to_grid_factor,
            1.0,
        ));

        let ray_end = math::lorentz_normalize(&na::Vector4::new(
            ray_end_grid_coords[0] / test_ctx.dual_to_grid_factor,
            ray_end_grid_coords[1] / test_ctx.dual_to_grid_factor,
            ray_end_grid_coords[2] / test_ctx.dual_to_grid_factor,
            1.0,
        ));

        let ray = Ray::new(
            ray_start,
            math::lorentz_normalize(
                &((ray_end - ray_start)
                    + ray_start * math::mip(&ray_start, &(ray_end - ray_start))),
            ),
        );

        let tanh_distance = (-math::mip(&ray_start, &ray_end)).acosh();

        let bounding_box = VoxelAABB::from_ray_segment_and_radius(
            test_ctx.dimension,
            test_ctx.dual_to_grid_factor,
            &ray,
            tanh_distance,
            test_ctx.collider_radius,
        )
        .unwrap();

        let ctx = ChunkShapeCastingContext {
            dimension: test_ctx.dimension,
            dual_to_grid_factor: test_ctx.dual_to_grid_factor,
            chunk: ChunkId::new(NodeId::ROOT, Vertex::A),
            transform: test_ctx.transform,
            voxel_data: &test_ctx.voxel_data,
            collider_radius: test_ctx.collider_radius,
            ray: &ray,
            bounding_box,
        };

        wrapped_fn(&ctx, tanh_distance)
    }

    fn chunk_shape_cast_wrapper(ctx: &ChunkShapeCastingContext, tanh_distance: f32) -> RayEndpoint {
        let mut endpoint = RayEndpoint {
            tanh_distance,
            hit: None,
        };
        chunk_shape_cast(ctx, &mut endpoint);
        endpoint
    }

    fn find_face_collision_wrapper(
        ctx: &ChunkShapeCastingContext,
        t_axis: usize,
        tanh_distance: f32,
    ) -> RayEndpoint {
        let mut endpoint = RayEndpoint {
            tanh_distance,
            hit: None,
        };
        find_face_collision(ctx, t_axis, &mut endpoint);
        endpoint
    }

    fn find_edge_collision_wrapper(
        ctx: &ChunkShapeCastingContext,
        t_axis: usize,
        tanh_distance: f32,
    ) -> RayEndpoint {
        let mut endpoint = RayEndpoint {
            tanh_distance,
            hit: None,
        };
        find_edge_collision(ctx, t_axis, &mut endpoint);
        endpoint
    }

    fn find_vertex_collision_wrapper(
        ctx: &ChunkShapeCastingContext,
        tanh_distance: f32,
    ) -> RayEndpoint {
        let mut endpoint = RayEndpoint {
            tanh_distance,
            hit: None,
        };
        find_vertex_collision(ctx, &mut endpoint);
        endpoint
    }

    fn test_face_collision(ctx: &ChunkShapeCastingContext, t_axis: usize, tanh_distance: f32) {
        let endpoint = chunk_shape_cast_wrapper(ctx, tanh_distance);
        assert_endpoints_hit_and_eq(
            &endpoint,
            &find_face_collision_wrapper(ctx, t_axis, tanh_distance),
        );
        sanity_check_normal(ctx, &endpoint);
    }

    fn test_edge_collision(ctx: &ChunkShapeCastingContext, t_axis: usize, tanh_distance: f32) {
        let endpoint = chunk_shape_cast_wrapper(ctx, tanh_distance);
        assert_endpoints_hit_and_eq(
            &endpoint,
            &find_edge_collision_wrapper(ctx, t_axis, tanh_distance),
        );
        sanity_check_normal(ctx, &endpoint);
    }

    fn test_vertex_collision(ctx: &ChunkShapeCastingContext, tanh_distance: f32) {
        let endpoint = chunk_shape_cast_wrapper(ctx, tanh_distance);
        assert_endpoints_hit_and_eq(
            &endpoint,
            &find_vertex_collision_wrapper(ctx, tanh_distance),
        );
        sanity_check_normal(ctx, &endpoint);
    }

    fn assert_endpoints_hit_and_eq(endpoint0: &RayEndpoint, endpoint1: &RayEndpoint) {
        assert_eq!(endpoint0.tanh_distance, endpoint1.tanh_distance);
        assert!(endpoint0.hit.is_some());
        assert!(endpoint1.hit.is_some());
        assert_eq!(
            endpoint0.hit.as_ref().unwrap().normal,
            endpoint1.hit.as_ref().unwrap().normal
        );
    }

    /// Ensures that the normal is pointing outward, opposite the ray direction.
    fn sanity_check_normal(ctx: &ChunkShapeCastingContext, endpoint: &RayEndpoint) {
        // The ray we care about is after its start point has moved to the contact point.
        let ray = math::translate(
            &ctx.ray.position,
            &math::lorentz_normalize(&ctx.ray.ray_point(endpoint.tanh_distance)),
        ) * ctx.ray;

        // The ray should be in the normal's coordinate system.
        let ray = ctx.transform * &ray;

        let normal = endpoint.hit.as_ref().unwrap().normal;

        // Project normal to be perpendicular to the ray's position
        let corrected_normal =
            math::lorentz_normalize(&(normal + ray.position * math::mip(&normal, &ray.position)));

        // Check that the normal and ray are pointing opposite directions
        assert!(math::mip(&corrected_normal, &ray.direction) < 0.0);
    }

    /// Tests that a suitable collision is found when approaching a single voxel from various angles.
    #[test]
    fn chunk_shape_cast_examples() {
        let collider_radius = 0.02;
        let test_ctx = TestShapeCastContext::new(collider_radius);

        // Approach a single voxel from various angles. Ensure that a suitable collision is found each time.

        // Face collisions
        cast_with_test_ray(
            &test_ctx,
            [0.0, 1.5, 1.5],
            [1.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                test_face_collision(ctx, 0, tanh_distance);
            },
        );

        cast_with_test_ray(
            &test_ctx,
            [1.5, 1.5, 3.0],
            [1.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                test_face_collision(ctx, 2, tanh_distance);
            },
        );

        // Edge collisions
        cast_with_test_ray(
            &test_ctx,
            [1.5, 3.0, 0.0],
            [1.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                test_edge_collision(ctx, 0, tanh_distance);
            },
        );

        cast_with_test_ray(
            &test_ctx,
            [3.0, 1.5, 3.0],
            [1.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                test_edge_collision(ctx, 1, tanh_distance);
            },
        );

        // Vertex collisions
        cast_with_test_ray(
            &test_ctx,
            [0.0, 0.0, 0.0],
            [1.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                test_vertex_collision(ctx, tanh_distance);
            },
        );

        cast_with_test_ray(
            &test_ctx,
            [3.0, 3.0, 0.0],
            [1.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                test_vertex_collision(ctx, tanh_distance);
            },
        );
    }

    /// Tests that colliding with a face from the back side is impossible
    #[test]
    fn face_collisions_one_sided() {
        let collider_radius = 0.01;
        let test_ctx = TestShapeCastContext::new(collider_radius);

        cast_with_test_ray(
            &test_ctx,
            [1.5, 1.5, 1.5],
            [4.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                assert!(chunk_shape_cast_wrapper(ctx, tanh_distance).hit.is_none());
            },
        )
    }

    /// Tests that collisions aren't detected past the ray's endpoint
    #[test]
    fn no_collisions_past_ray_endpoint() {
        let collider_radius = 0.01;
        let test_ctx = TestShapeCastContext::new(collider_radius);

        cast_with_test_ray(
            &test_ctx,
            [8.0, 1.5, 1.5],
            [2.5, 1.5, 1.5],
            |ctx, tanh_distance| {
                assert!(chunk_shape_cast_wrapper(ctx, tanh_distance).hit.is_none());
            },
        )
    }

    #[test]
    fn solve_sphere_plane_intersection_example() {
        // Hit the z=0 plane with a radius of 0.2
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(math::origin(), na::Vector4::new(0.8, 0.0, 0.6, 0.0));
        let normal = -na::Vector4::z();
        let hit_point = math::lorentz_normalize(
            &ray.ray_point(solve_sphere_plane_intersection(&ray, &normal, 0.2_f32.sinh()).unwrap()),
        );
        assert_abs_diff_eq!(
            math::mip(&hit_point, &normal),
            0.2_f32.sinh(),
            epsilon = 1e-4
        );
    }

    #[test]
    fn solve_sphere_plane_intersection_direct_hit() {
        // Directly hit the z=0 plane with a ray 0.5 units away and a radius of 0.2.
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(math::origin(), na::Vector4::z());
        let normal = -na::Vector4::z();
        assert_abs_diff_eq!(
            solve_sphere_plane_intersection(&ray, &normal, 0.2_f32.sinh()).unwrap(),
            0.3_f32.tanh(),
            epsilon = 1e-4
        );
    }

    #[test]
    fn solve_sphere_plane_intersection_miss() {
        // No collision with the plane anywhere along the ray's line
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(math::origin(), na::Vector4::x());
        let normal = -na::Vector4::z();
        assert!(solve_sphere_plane_intersection(&ray, &normal, 0.2_f32.sinh()).is_none());
    }

    #[test]
    fn solve_sphere_plane_intersection_margin() {
        // Sphere is already contacting the plane, with some error
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.2))
            * &Ray::new(math::origin(), na::Vector4::z());
        let normal = -na::Vector4::z();
        assert_eq!(
            solve_sphere_plane_intersection(&ray, &normal, 0.2001_f32.sinh()).unwrap(),
            0.0
        );
    }

    #[test]
    fn solve_sphere_line_intersection_example() {
        // Hit the x=z=0 line with a radius of 0.2
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(
                math::origin(),
                na::Vector4::new(1.0, 2.0, 3.0, 0.0).normalize(),
            );
        let line_position = na::Vector4::w();
        let line_direction = na::Vector4::y();
        let hit_point = math::lorentz_normalize(
            &ray.ray_point(
                solve_sphere_line_intersection(
                    &ray,
                    &line_position,
                    &line_direction,
                    0.2_f32.cosh(),
                )
                .unwrap(),
            ),
        );
        // Measue the distance from hit_point to the line and ensure it's equal to the radius
        assert_abs_diff_eq!(
            (math::mip(&hit_point, &line_position).powi(2)
                - math::mip(&hit_point, &line_direction).powi(2))
            .sqrt(),
            0.2_f32.cosh(),
            epsilon = 1e-4
        );
    }

    #[test]
    fn solve_sphere_line_intersection_direct_hit() {
        // Directly hit the x=z=0 line with a ray 0.5 units away and a radius of 0.2.

        // Ensure the ray is slightly off-center so that the distance math is shown to be correct
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.7, 0.0))
            * math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(math::origin(), na::Vector4::z());
        let line_position = na::Vector4::w();
        let line_direction = na::Vector4::y();
        assert_abs_diff_eq!(
            solve_sphere_line_intersection(&ray, &line_position, &line_direction, 0.2_f32.cosh())
                .unwrap(),
            0.3_f32.tanh(),
            epsilon = 1e-4
        );
    }

    #[test]
    fn solve_sphere_line_intersection_miss() {
        // No collision with the line anywhere along the ray's line
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(math::origin(), na::Vector4::x());
        let line_position = na::Vector4::w();
        let line_direction = na::Vector4::y();
        assert!(solve_sphere_line_intersection(
            &ray,
            &line_position,
            &line_direction,
            0.2_f32.cosh()
        )
        .is_none());
    }

    #[test]
    fn solve_sphere_line_intersection_margin() {
        // Sphere is already contacting the line, with some error
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.2))
            * &Ray::new(math::origin(), na::Vector4::z());
        let line_position = na::Vector4::w();
        let line_direction = na::Vector4::y();
        assert_eq!(
            solve_sphere_line_intersection(
                &ray,
                &line_position,
                &line_direction,
                0.2001_f32.cosh()
            )
            .unwrap(),
            0.0
        );
    }

    #[test]
    fn solve_sphere_point_intersection_example() {
        // Hit the origin with a radius of 0.2
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(
                math::origin(),
                na::Vector4::new(1.0, 2.0, 6.0, 0.0).normalize(),
            );
        let point_position = math::origin();
        let hit_point = math::lorentz_normalize(&ray.ray_point(
            solve_sphere_point_intersection(&ray, &point_position, 0.2_f32.cosh()).unwrap(),
        ));
        assert_abs_diff_eq!(
            -math::mip(&hit_point, &point_position),
            0.2_f32.cosh(),
            epsilon = 1e-4
        );
    }

    #[test]
    fn solve_sphere_point_intersection_direct_hit() {
        // Directly hit the origin with a ray 0.5 units away and a radius of 0.2.
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(math::origin(), na::Vector4::z());
        let point_position = math::origin();
        assert_abs_diff_eq!(
            solve_sphere_point_intersection(&ray, &point_position, 0.2_f32.cosh()).unwrap(),
            0.3_f32.tanh(),
            epsilon = 1e-4
        );
    }

    #[test]
    fn solve_sphere_point_intersection_miss() {
        // No collision with the point anywhere along the ray's line
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.5))
            * &Ray::new(math::origin(), na::Vector4::x());
        let point_position = math::origin();
        assert!(solve_sphere_point_intersection(&ray, &point_position, 0.2_f32.cosh()).is_none());
    }

    #[test]
    fn solve_sphere_point_intersection_margin() {
        // Sphere is already contacting the point, with some error
        let ray = math::translate_along(&na::Vector3::new(0.0, 0.0, -0.2))
            * &Ray::new(math::origin(), na::Vector4::z());
        let point_position = math::origin();
        assert_eq!(
            solve_sphere_point_intersection(&ray, &point_position, 0.2001_f32.cosh()).unwrap(),
            0.0
        );
    }

    #[test]
    fn solve_quadratic_example() {
        let a = 1.0;
        let b = 2.0;
        let c = -5.0;
        let x = solve_quadratic(c, b / 2.0, a).unwrap();

        // x should be a solution
        assert_abs_diff_eq!(a * x * x + b * x + c, 0.0, epsilon = 1e-4);

        // x should be the smallest solution, less than the parabola's vertex.
        assert!(x < -b / (2.0 * a));
    }

    #[test]
    fn tuv_to_xyz_example() {
        assert_eq!(tuv_to_xyz(0, [2, 4, 6]), [2, 4, 6]);
        assert_eq!(tuv_to_xyz(1, [2, 4, 6]), [6, 2, 4]);
        assert_eq!(tuv_to_xyz(2, [2, 4, 6]), [4, 6, 2]);

        assert_eq!(tuv_to_xyz(1, [2, 4, 6, 8]), [6, 2, 4, 8]);
    }
}
