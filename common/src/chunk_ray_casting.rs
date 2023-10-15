use crate::{
    graph_collision::Ray,
    math,
    node::{ChunkLayout, Coords, VoxelData},
    world::Material,
};

pub struct ChunkCastHit {
    /// The tanh of the distance traveled along the ray to result in this hit.
    pub tanh_distance: f32,

    /// The coordinates of the block that was hit, including margins.
    pub voxel_coords: Coords,

    /// Which of the three axes is orthogonal to the face of the block that was hit.
    pub face_axis: u32,

    // Either +1 or -1, depending on whether the outside of the face that was hit was in the positive or
    // negative direction in `face_axis`.
    pub face_direction: i8,
}

/// Performs sphere casting (swept collision query) against the voxels in the chunk with the given `voxel_data`
///
/// The `ray` parameter is given and any resulting hit normals are given in the chunk's dual coordinate system.
///
/// The `tanh_distance` is the hyperbolic tangent of the distance along the ray to check for hits.
pub fn chunk_ray_cast(
    voxel_data: &VoxelData,
    layout: &ChunkLayout,
    ray: &Ray,
    tanh_distance: f32,
) -> Option<ChunkCastHit> {
    let mut hit: Option<ChunkCastHit> = None;

    let Some(bounding_box) = VoxelAABB::from_ray_segment(layout, ray, tanh_distance) else {
        return None;
    };

    for t_axis in 0..3 {
        hit = find_face_collision(
            voxel_data,
            layout,
            &bounding_box,
            t_axis,
            ray,
            hit.as_ref().map_or(tanh_distance, |hit| hit.tanh_distance),
        )
        .or(hit);
    }

    hit
}

/// Detect collisions where a sphere contacts the front side of a voxel face
fn find_face_collision(
    voxel_data: &VoxelData,
    layout: &ChunkLayout,
    bounding_box: &VoxelAABB,
    t_axis: usize,
    ray: &Ray,
    tanh_distance: f32,
) -> Option<ChunkCastHit> {
    let mut hit: Option<ChunkCastHit> = None;

    let u_axis = (t_axis + 1) % 3;
    let v_axis = (t_axis + 2) % 3;

    // Loop through all grid planes overlapping the bounding box
    for t in bounding_box.grid_planes(t_axis) {
        // Find a normal to the grid plane. Note that (t, 0, 0, x) is a normal of the plane whose closest point
        // to the origin is (x, 0, 0, t), and we use that fact here.
        let normal = math::lorentz_normalize(&tuv_to_xyz(
            t_axis,
            na::Vector4::new(1.0, 0.0, 0.0, layout.grid_to_dual(t)),
        ));

        let Some(new_tanh_distance) = solve_point_plane_intersection(ray, &normal) else {
            continue;
        };

        // If new_tanh_distance is out of range, no collision occurred.
        if new_tanh_distance >= hit.as_ref().map_or(tanh_distance, |hit| hit.tanh_distance) {
            continue;
        }

        // Whether we are approaching the front or back of the face. An approach from the positive t direction
        // is 1, and an approach from the negative t direction is -1.
        let collision_side = -math::mip(&ray.direction, &normal).signum();

        // Which side we approach the plane from affects which voxel we want to use for collision checking
        let voxel_t = if collision_side > 0.0 { t } else { t + 1 };

        let ray_endpoint = ray.ray_point(new_tanh_distance);
        let contact_point = ray_endpoint - normal * math::mip(&ray_endpoint, &normal);

        // Compute the u and v-coordinates of the voxels at the contact point
        let Some(voxel_u) = layout.dual_to_voxel(contact_point[u_axis] / contact_point.w) else {
            continue;
        };
        let Some(voxel_v) = layout.dual_to_voxel(contact_point[v_axis] / contact_point.w) else {
            continue;
        };

        // Ensure that the relevant voxel is solid
        if !voxel_is_solid(
            voxel_data,
            layout,
            tuv_to_xyz(t_axis, [voxel_t, voxel_u, voxel_v]),
        ) {
            continue;
        }

        // A collision was found. Update the hit.
        hit = Some(ChunkCastHit {
            tanh_distance: new_tanh_distance,
            voxel_coords: Coords(tuv_to_xyz(
                t_axis,
                [voxel_t as u8, voxel_u as u8, voxel_v as u8],
            )),
            face_axis: t_axis as u32,
            face_direction: collision_side as i8,
        });
    }

    hit
}

/// Finds the tanh of the distance a point will have to travel along a ray before it
/// intersects the given plane.
fn solve_point_plane_intersection(ray: &Ray, plane_normal: &na::Vector4<f32>) -> Option<f32> {
    let mip_pos_a = math::mip(&ray.position, plane_normal);
    let mip_dir_a = math::mip(&ray.direction, plane_normal);

    let result = -mip_pos_a / mip_dir_a;
    if result.is_finite() && result > 0.0 {
        Some(result)
    } else {
        None
    }
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

/// Checks whether a voxel can be collided with. Any non-void voxel falls under this category.
fn voxel_is_solid(voxel_data: &VoxelData, layout: &ChunkLayout, coords: [usize; 3]) -> bool {
    let dimension_with_margin = layout.dimension() + 2;
    debug_assert!(coords[0] < dimension_with_margin);
    debug_assert!(coords[1] < dimension_with_margin);
    debug_assert!(coords[2] < dimension_with_margin);
    voxel_data.get(
        coords[0] + coords[1] * dimension_with_margin + coords[2] * dimension_with_margin.pow(2),
    ) != Material::Void
}

/// Represents a discretized region in the voxel grid contained by an axis-aligned bounding box.
struct VoxelAABB {
    // The bounds are of the form [[x_min, x_max], [y_min, y_max], [z_min, z_max]], using voxel coordinates with margins.
    // Any voxel that intersects the cube of interest is included in these bounds. By adding or subtracting 1 in the right
    // places, these bounds can be used to find other useful info related to the cube of interset, such as what grid points
    // it contains.
    bounds: [[usize; 2]; 3],
}

impl VoxelAABB {
    /// Returns a bounding box that is guaranteed to cover a given radius around a ray segment. Returns None if the
    /// bounding box lies entirely outside the chunk, including its margins.
    pub fn from_ray_segment(
        layout: &ChunkLayout,
        ray: &Ray,
        tanh_distance: f32,
    ) -> Option<VoxelAABB> {
        // Convert the ray to grid coordinates
        let grid_start =
            na::Point3::from_homogeneous(ray.position).unwrap() * layout.dual_to_grid_factor();
        let grid_end = na::Point3::from_homogeneous(ray.ray_point(tanh_distance)).unwrap()
            * layout.dual_to_grid_factor();
        let mut bounds = [[0; 2]; 3];
        for axis in 0..3 {
            let grid_min = grid_start[axis].min(grid_end[axis]);
            let grid_max = grid_start[axis].max(grid_end[axis]);
            let voxel_min = (grid_min + 1.0).floor().max(0.0);
            let voxel_max = (grid_max + 1.0)
                .floor()
                .min(layout.dimension() as f32 + 1.0);

            // This will happen when voxel_min is greater than dimension+1 or voxel_max is less than 0, which
            // occurs when the cube is out of range.
            if voxel_min > voxel_max {
                return None;
            }

            // We convert to usize here instead of earlier because out-of-range voxel coordinates can be negative.
            bounds[axis] = [voxel_min.floor() as usize, voxel_max.floor() as usize];
        }

        Some(VoxelAABB { bounds })
    }

    /// Iterator over grid planes intersecting the region, represented as integers determining the plane's fixed coordinate
    pub fn grid_planes(&self, axis: usize) -> impl Iterator<Item = usize> {
        self.bounds[axis][0]..self.bounds[axis][1]
    }
}
