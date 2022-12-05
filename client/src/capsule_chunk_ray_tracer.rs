use common::{dodeca::Vertex, math, world::Material};

use crate::chunk_ray_tracer::{ChunkRayTracer, RayTracingResultHandle, VoxelDataWrapper};

pub struct CapsuleChunkRayTracer {
    pub radius: f64,
    pub height: f64,
}

impl ChunkRayTracer for CapsuleChunkRayTracer {
    fn trace_ray_in_chunk(
        &self,
        voxel_data: VoxelDataWrapper,
        pos: &na::Vector4<f64>,
        dir: &na::Vector4<f64>,
        handle: &mut RayTracingResultHandle,
    ) {
        CapsuleChunkRayTracingPass::new(self.radius, self.height, voxel_data, pos, dir, handle)
            .trace_ray_in_chunk();
    }

    fn max_radius(&self) -> f64 {
        self.radius + self.height
    }
}

struct CapsuleChunkRayTracingPass<'a, 'b> {
    radius: f64,
    height: f64,
    voxel_data: VoxelDataWrapper<'a>,
    pos: &'a na::Vector4<f64>,
    dir: &'a na::Vector4<f64>,
    handle: &'a mut RayTracingResultHandle<'b>,

    // Start and end of region to check in voxel coordinates
    // TODO: These can be used for more fine-tuned pre-collision-check filtering
    #[allow(dead_code)]
    voxel_start: na::Vector3<f64>,
    #[allow(dead_code)]
    voxel_end: na::Vector3<f64>,

    // Bounding box of all voxels that can be collided with
    // [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
    bbox: [[usize; 2]; 3],
}

impl CapsuleChunkRayTracingPass<'_, '_> {
    fn new<'a, 'b>(
        radius: f64,
        height: f64,
        voxel_data: VoxelDataWrapper<'a>,
        pos: &'a na::Vector4<f64>,
        dir: &'a na::Vector4<f64>,
        handle: &'a mut RayTracingResultHandle<'b>,
    ) -> CapsuleChunkRayTracingPass<'a, 'b> {
        let float_size = voxel_data.dimension() as f64;
        let voxel_start = (pos / pos[3]).xyz() * Vertex::dual_to_chunk_factor() * float_size;
        let end_pos = pos + dir * handle.t();
        let voxel_end = (end_pos / end_pos[3]).xyz() * Vertex::dual_to_chunk_factor() * float_size;
        let max_voxel_radius = (radius + height) * Vertex::dual_to_chunk_factor() * float_size;
        let bbox = [0, 1, 2].map(|coord| {
            get_usize_range(
                voxel_data.dimension(),
                voxel_start[coord],
                voxel_end[coord],
                max_voxel_radius,
            )
        });

        CapsuleChunkRayTracingPass {
            radius,
            height,
            voxel_data,
            pos,
            dir,
            handle,
            voxel_start,
            voxel_end,
            bbox,
        }
    }

    fn trace_ray_in_chunk(&mut self) {
        for coord_axis in 0..3 {
            self.trace_ray_for_sides(coord_axis);
        }

        for coord_axis in 0..3 {
            self.trace_ray_for_edges(coord_axis);
        }

        self.trace_ray_for_vertices();
    }

    fn trace_ray_for_quad(
        &mut self,
        vertex0: &na::Vector4<f64>,
        vertex1: &na::Vector4<f64>,
        vertex2: &na::Vector4<f64>,
        vertex3: &na::Vector4<f64>,
    ) {
        // Compute triangle normal
        let normal = math::lorentz_normalize(&math::triangle_normal(vertex0, vertex1, vertex2));
        self.trace_ray_for_sphere_polygon(&normal, [vertex0, vertex1, vertex2, vertex3]);
    }

    fn trace_ray_for_segment(
        &mut self,
        endpoint0: &na::Vector4<f64>,
        endpoint1: &na::Vector4<f64>,
    ) {
        self.trace_ray_for_sphere_segment(endpoint0, endpoint1);

        let lower_endpoint0 = self.get_lower_point(endpoint0);
        let lower_endpoint1 = self.get_lower_point(endpoint1);

        self.trace_ray_for_sphere_segment(&lower_endpoint0, &lower_endpoint1);

        let normal = math::lorentz_normalize(&math::triangle_normal(endpoint0, endpoint1, &lower_endpoint1));
        self.trace_ray_for_sphere_polygon(&normal, [endpoint0, endpoint1, &lower_endpoint1, &lower_endpoint0]);
        self.trace_ray_for_sphere_polygon(&(-normal), [endpoint1, endpoint0, &lower_endpoint0, &lower_endpoint1]);
    }

    fn trace_ray_for_point(&mut self, point: &na::Vector4<f64>) {
        self.trace_ray_for_sphere_point(point);
    }

    fn get_lower_point(&mut self, point: &na::Vector4<f64>) -> na::Vector4<f64> {
        let up = math::lorentz_normalize(
            &(self.voxel_data.up() + point * math::mip(self.voxel_data.up(), point)),
        );

        point * self.height.cosh() - up * self.height.sinh()
    }

    fn trace_ray_for_sphere_polygon<const N: usize>(
        &mut self,
        normal: &na::Vector4<f64>,
        vertices: [&na::Vector4<f64>; N],
    ) {
        if math::mip(self.dir, normal) >= 0.0 {
            return;
        }

        let t_candidate =
            find_intersection_one_vector(self.pos, self.dir, normal, self.radius.sinh());

        if !(t_candidate >= 0.0 && t_candidate < self.handle.t()) {
            return;
        }

        let new_pos = self.pos + self.dir * t_candidate;
        let projected_pos = math::project_ortho(&new_pos, normal);

        // Check if inside the triangle. Project to 2D by removing w coordinate and largest of xyz in the normal.
        let largest_xyz_coord = (0..3)
            .max_by(|&x, &y| {
                normal[x]
                    .abs()
                    .partial_cmp(&normal[y].abs())
                    .expect("not NaN")
            })
            .unwrap();

        let i0 = (largest_xyz_coord + 1) % 3;
        let i1 = (largest_xyz_coord + 2) % 3;

        let is_on_correct_side = |v0: &na::Vector4<f64>, v1: &na::Vector4<f64>| {
            ((v1[i0] - v0[i0]) * (projected_pos[i1] - v0[i1])
                - (v1[i1] - v0[i1]) * (projected_pos[i0] - v0[i0]))
                .signum()
                == normal[largest_xyz_coord].signum()
        };

        if (0..N).all(|n| is_on_correct_side(vertices[n], vertices[(n + 1) % N])) {
            self.handle.update(t_candidate, [0, 0, 0], 0, 0, *normal);
        }
    }

    fn trace_ray_for_sphere_segment(
        &mut self,
        endpoint0: &na::Vector4<f64>,
        endpoint1: &na::Vector4<f64>,
    ) {
        let segment_dir =
            math::lorentz_normalize(&(endpoint1 + endpoint0 * math::mip(endpoint1, endpoint0)));

        let t_candidate = find_intersection_two_vectors(
            self.pos,
            self.dir,
            endpoint0,
            &segment_dir,
            self.radius.cosh(),
        );

        if !(t_candidate >= 0.0 && t_candidate < self.handle.t()) {
            return;
        }

        let new_pos = self.pos + self.dir * t_candidate;
        let projected_pos = math::lorentz_normalize(
            &(-endpoint0 * math::mip(&new_pos, endpoint0)
                + segment_dir * math::mip(&new_pos, &segment_dir)),
        );

        let location_on_segment = math::mip(&projected_pos, &segment_dir);
        if location_on_segment >= 0.0 && location_on_segment <= math::mip(endpoint1, &segment_dir) {
            self.handle
                .update(t_candidate, [0, 0, 0], 0, 0, new_pos - projected_pos);
        }
    }

    fn trace_ray_for_sphere_point(&mut self, point: &na::Vector4<f64>) {
        let t_candidate =
            find_intersection_one_vector(self.pos, self.dir, point, self.radius.cosh());

        if !(t_candidate >= 0.0 && t_candidate < self.handle.t()) {
            return;
        }

        let translated_square_pos = self.pos + self.dir * t_candidate;
        self.handle
            .update(t_candidate, [0, 0, 0], 0, 0, translated_square_pos - point);
    }

    fn trace_ray_for_sides(&mut self, coord_axis: usize) {
        let float_size = self.voxel_data.dimension() as f64;
        let coord_plane0 = (coord_axis + 1) % 3;
        let coord_plane1 = (coord_axis + 2) % 3;

        let coords_to_vertex = |i, j, k| {
            let mut vertex = na::Vector4::zeros();
            vertex[coord_axis] = i as f64 / float_size * Vertex::chunk_to_dual_factor();
            vertex[coord_plane0] = j as f64 / float_size * Vertex::chunk_to_dual_factor();
            vertex[coord_plane1] = k as f64 / float_size * Vertex::chunk_to_dual_factor();
            vertex.w = 1.0;
            math::lorentz_normalize(&vertex)
        };

        let get_voxel_data = |data: &VoxelDataWrapper, i, j, k| {
            let mut coords = [0, 0, 0];
            coords[coord_axis] = i;
            coords[coord_plane0] = j;
            coords[coord_plane1] = k;
            data.get(coords)
        };

        for i in self.bbox[coord_axis][0]..self.bbox[coord_axis][1] {
            for j in self.bbox[coord_plane0][0].max(1) - 1
                ..self.bbox[coord_plane0][1].min(self.voxel_data.dimension())
            {
                for k in self.bbox[coord_plane1][0].max(1) - 1
                    ..self.bbox[coord_plane1][1].min(self.voxel_data.dimension())
                {
                    let mut vertices = [[na::Vector4::zeros(); 2]; 2];
                    for (square0, vertices) in vertices.iter_mut().enumerate() {
                        for (square1, vertex) in vertices.iter_mut().enumerate() {
                            *vertex = coords_to_vertex(i, j + square0, k + square1);
                        }
                    }

                    if i > 0 && get_voxel_data(&self.voxel_data, i - 1, j, k) != Material::Void {
                        self.trace_ray_for_quad(
                            &vertices[0][0],
                            &vertices[1][0],
                            &vertices[1][1],
                            &vertices[0][1],
                        );
                    }

                    if i < self.voxel_data.dimension()
                        && get_voxel_data(&self.voxel_data, i, j, k) != Material::Void
                    {
                        self.trace_ray_for_quad(
                            &vertices[0][0],
                            &vertices[0][1],
                            &vertices[1][1],
                            &vertices[1][0],
                        );
                    }
                }
            }
        }
    }

    fn trace_ray_for_edges(&mut self, coord_axis: usize) {
        let float_size = self.voxel_data.dimension() as f64;
        let coord_plane0 = (coord_axis + 1) % 3;
        let coord_plane1 = (coord_axis + 2) % 3;

        let coords_to_vertex = |i, j, k| {
            let mut vertex = na::Vector4::zeros();
            vertex[coord_axis] = i as f64 / float_size * Vertex::chunk_to_dual_factor();
            vertex[coord_plane0] = j as f64 / float_size * Vertex::chunk_to_dual_factor();
            vertex[coord_plane1] = k as f64 / float_size * Vertex::chunk_to_dual_factor();
            vertex.w = 1.0;
            math::lorentz_normalize(&vertex)
        };

        let get_voxel_data = |data: &VoxelDataWrapper, i, j, k| {
            let mut coords = [0, 0, 0];
            coords[coord_axis] = i;
            coords[coord_plane0] = j;
            coords[coord_plane1] = k;
            data.get(coords)
        };

        for i in self.bbox[coord_axis][0].max(1) - 1
            ..self.bbox[coord_axis][1].min(self.voxel_data.dimension())
        {
            for j in self.bbox[coord_plane0][0]..self.bbox[coord_plane0][1] {
                for k in self.bbox[coord_plane1][0]..self.bbox[coord_plane1][1] {
                    let mut vertices = [na::Vector4::zeros(); 2];
                    for (square0, vertex) in vertices.iter_mut().enumerate() {
                        *vertex = coords_to_vertex(i + square0, j, k);
                    }

                    if (j > 0
                        && k > 0
                        && get_voxel_data(&self.voxel_data, i, j - 1, k - 1) != Material::Void)
                        || (j < self.voxel_data.dimension()
                            && k > 0
                            && get_voxel_data(&self.voxel_data, i, j, k - 1) != Material::Void)
                        || (j > 0
                            && k < self.voxel_data.dimension()
                            && get_voxel_data(&self.voxel_data, i, j - 1, k) != Material::Void)
                        || (j < self.voxel_data.dimension()
                            && k < self.voxel_data.dimension()
                            && get_voxel_data(&self.voxel_data, i, j, k) != Material::Void)
                    {
                        self.trace_ray_for_segment(&vertices[0], &vertices[1]);
                    }
                }
            }
        }
    }

    fn trace_ray_for_vertices(&mut self) {
        let size = self.voxel_data.dimension();
        let float_size = size as f64;

        for i in self.bbox[0][0]..self.bbox[0][1] {
            for j in self.bbox[1][0]..self.bbox[1][1] {
                for k in self.bbox[2][0]..self.bbox[2][1] {
                    if (i == 0
                        || j == 0
                        || k == 0
                        || self.voxel_data.get([i - 1, j - 1, k - 1]) == Material::Void)
                        && (i == size
                            || j == 0
                            || k == 0
                            || self.voxel_data.get([i, j - 1, k - 1]) == Material::Void)
                        && (i == 0
                            || j == size
                            || k == 0
                            || self.voxel_data.get([i - 1, j, k - 1]) == Material::Void)
                        && (i == size
                            || j == size
                            || k == 0
                            || self.voxel_data.get([i, j, k - 1]) == Material::Void)
                        && (i == 0
                            || j == 0
                            || k == size
                            || self.voxel_data.get([i - 1, j - 1, k]) == Material::Void)
                        && (i == size
                            || j == 0
                            || k == size
                            || self.voxel_data.get([i, j - 1, k]) == Material::Void)
                        && (i == 0
                            || j == size
                            || k == size
                            || self.voxel_data.get([i - 1, j, k]) == Material::Void)
                        && (i == size
                            || j == size
                            || k == size
                            || self.voxel_data.get([i, j, k]) == Material::Void)
                    {
                        continue;
                    }

                    let vert = math::lorentz_normalize(&na::Vector4::new(
                        i as f64 / float_size * Vertex::chunk_to_dual_factor(),
                        j as f64 / float_size * Vertex::chunk_to_dual_factor(),
                        k as f64 / float_size * Vertex::chunk_to_dual_factor(),
                        1.0,
                    ));

                    self.trace_ray_for_point(&vert);
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
fn find_intersection_one_vector(
    pos: &na::Vector4<f64>,
    dir: &na::Vector4<f64>,
    a: &na::Vector4<f64>,
    c: f64,
) -> f64 {
    let mip_pos_a = math::mip(pos, a);
    let mip_dir_a = math::mip(dir, a);

    // The following 3 variables are terms of the quadratic formula. We use double the linear
    // term because it removes the annoying constants from that formula.
    let quadratic_term = mip_dir_a.powi(2) + c.powi(2);
    let double_linear_term = mip_pos_a * mip_dir_a;
    let constant_term = mip_pos_a.powi(2) - c.powi(2);

    let discriminant = double_linear_term * double_linear_term - quadratic_term * constant_term;

    // While discriminant can be negative, NaNs propagate the way we want to, so we don't have
    // to check for this.
    (-double_linear_term - discriminant.sqrt()) / quadratic_term
}

/// Find the smallest value of `t` where the point in the pos-dir line (v=pos+dir*t) satisfies
/// `(<v,a>^2 - <v,b>^2) / <v,v> == c^2`
///
/// This finds intersections with a surface that is `cosh(c)` units away from the line
/// with a point at `a` with direction `b`, where `<a,b>==0`.
///
/// Returns NaN if there's no such intersection
fn find_intersection_two_vectors(
    pos: &na::Vector4<f64>,
    dir: &na::Vector4<f64>,
    a: &na::Vector4<f64>,
    b: &na::Vector4<f64>,
    c: f64,
) -> f64 {
    // TODO: Improve numerical stability
    let mip_pos_a = math::mip(pos, a);
    let mip_dir_a = math::mip(dir, a);
    let mip_pos_b = math::mip(pos, b);
    let mip_dir_b = math::mip(dir, b);

    // The following 3 variables are terms of the quadratic formula. We use double the linear
    // term because it removes the annoying constants from that formula.
    let quadratic_term = mip_dir_a.powi(2) - mip_dir_b.powi(2) + c.powi(2);
    let double_linear_term = mip_pos_a * mip_dir_a - mip_pos_b * mip_dir_b;
    let constant_term = mip_pos_a.powi(2) - mip_pos_b.powi(2) - c.powi(2);

    let discriminant = double_linear_term * double_linear_term - quadratic_term * constant_term;

    // While discriminant can be negative, NaNs propagate the way we want to, so we don't have
    // to check for this.
    (-double_linear_term - discriminant.sqrt()) / quadratic_term
}

fn get_usize_range(dimension: usize, point0: f64, point1: f64, width: f64) -> [usize; 2] {
    if !point0.is_finite() || !point1.is_finite() {
        return [0, dimension + 1];
    }
    let result_min = (point0.min(point1) - width).max(-0.5) + 1.0;
    let result_max = (point0.max(point1) + width).min(dimension as f64 + 0.5) + 1.0;

    if result_min > result_max {
        // Empty range
        return [1, 0];
    }

    [result_min.floor() as usize, result_max.floor() as usize]
}
