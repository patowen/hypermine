use common::math::{MIsometry, MPoint, MVector, PermuteXYZ, sqr};
use libm::{coshf, logf, powf, sinhf, sqrtf, tanhf};

use crate::graphics::{
    Mesh,
    asset_loader::AssetLoadContext,
    meshes::{MeshGeometryDefinition, Vertex},
};

fn voxel_to_mvector_simple(voxel: na::Vector3<f32>) -> MVector<f32> {
    let factor = 1.0 / sqrtf(1.0 - sqr(voxel.x) - sqr(voxel.y));
    //MVector::new(voxel.x, voxel.y, factor * tanhf(voxel.z), 1.0)

    // This is already pre-scaled
    MVector::new(
        voxel.x * coshf(voxel.z) * factor,
        voxel.y * coshf(voxel.z) * factor,
        sinhf(voxel.z),
        coshf(voxel.z) * factor,
    )
}

fn pseudo_chunk_to_isometry(pseudo_chunk: na::Vector3<f32>, boost: f32) -> MIsometry<f32> {
    let w = voxel_to_mvector_boosted(pseudo_chunk, boost).to_point_unchecked();
    let x = voxel_to_mvector_boosted_partial_x(pseudo_chunk, boost).normalized_direction();
    let z = voxel_to_mvector_boosted_partial_z(pseudo_chunk, boost).to_direction_unchecked();
    // TODO: There might be a better way to get `y`, especially since voxel_to_mvector_boosted_partial_y will be
    // in the wrong direction.
    let mut y = voxel_to_mvector_boosted_partial_y(pseudo_chunk, boost);
    y -= MVector::from(x) * y.mip(&x);
    let y = y.normalized_direction();
    MIsometry::from_columns_unchecked(&[x, y, z], w)
}

fn chunk_to_mvector_simple(
    klein_coords: na::Vector2<f32>,
    chunk: na::Vector3<f32>,
    boost: f32,
) -> MVector<f32> {
    pseudo_chunk_to_isometry(
        na::Vector3::new(
            klein_coords.x * coshf(boost),
            klein_coords.y * coshf(boost),
            0.0,
        ),
        boost,
    )
    .inverse()
        * voxel_to_mvector_boosted(
            na::Vector3::new(
                klein_coords.x * coshf(boost) + chunk.x,
                klein_coords.y * coshf(boost) + chunk.y,
                chunk.z,
            ),
            boost,
        )
}

fn chunk_to_pseudo_chunk(
    klein_coords: na::Vector2<f32>,
    chunk: na::Vector3<f32>,
    boost: f32,
) -> na::Vector3<f32> {
    let origin = na::Vector3::new(klein_coords[0], klein_coords[1], 1.0);
    let origin_norm = sqrtf(-(sqr(origin[0]) + sqr(origin[1]) - sqr(origin[2])));
    let normalized_origin = origin / origin_norm;
    let z = normalized_origin;
    let mut y = z.cross(&na::Vector3::x());
    y.z *= -1.0;
    y /= sqrtf(sqr(y.x) + sqr(y.y) - sqr(y.z));
    let mut x = y.cross(&z);
    x.z *= -1.0;
    let mut conversion = na::Matrix3::from_columns(&[x, y, z]).try_inverse().unwrap();
    conversion *= na::Matrix3::new_translation(&klein_coords); // Applying skew

    let horizontal_coords = conversion * chunk.xy().push(1.0);
    na::Vector3::new(
        horizontal_coords[0] / horizontal_coords[2],
        horizontal_coords[1] / horizontal_coords[2],
        chunk.z,
    )

    // TODO: Apply boost
}

/// Computes
/// `translation_along([0, 0, -boost]) * voxel_to_mvector_simple([x / cosh(boost), y / cosh(boost), boost + z])`
fn voxel_to_mvector_boosted(pseudo_chunk: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        pseudo_chunk.x * factor * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        pseudo_chunk.y * factor * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        sinhf(pseudo_chunk.z) * (1.0 - scaled_factor_delta * sqr(tanhf(boost)))
            - scaled_factor_delta * tanhf(boost) * coshf(pseudo_chunk.z),
        coshf(pseudo_chunk.z) * (1.0 + scaled_factor_delta)
            + scaled_factor_delta * tanhf(boost) * sinhf(pseudo_chunk.z),
    )
}

fn voxel_to_mvector_boosted_partial_x(pseudo_chunk: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let dist_squared_partial_x = pseudo_chunk.x * 2.0;
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor_radicand_partial_x = -dist_squared_partial_x / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let factor_partial_x = -0.5 * factor_radicand_partial_x * factor / factor_radicand;
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    let scaled_factor_delta_partial_x = (dist_squared_partial_x
        * (factor_radicand + sqrtf(factor_radicand))
        - dist_squared
            * (factor_radicand_partial_x
                + 0.5 / sqrtf(factor_radicand) * factor_radicand_partial_x))
        / sqr(factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        (factor + pseudo_chunk.x * factor_partial_x)
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        pseudo_chunk.y
            * factor_partial_x
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        sinhf(pseudo_chunk.z) * (-scaled_factor_delta_partial_x * sqr(tanhf(boost)))
            - scaled_factor_delta_partial_x * tanhf(boost) * coshf(pseudo_chunk.z),
        coshf(pseudo_chunk.z) * (scaled_factor_delta_partial_x)
            + scaled_factor_delta_partial_x * tanhf(boost) * sinhf(pseudo_chunk.z),
    )
}

fn voxel_to_mvector_boosted_partial_y(pseudo_chunk: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let dist_squared_partial_y = pseudo_chunk.y * 2.0;
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor_radicand_partial_y = -dist_squared_partial_y / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let factor_partial_y = -0.5 * factor_radicand_partial_y * factor / factor_radicand;
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    let scaled_factor_delta_partial_y = (dist_squared_partial_y
        * (factor_radicand + sqrtf(factor_radicand))
        - dist_squared
            * (factor_radicand_partial_y
                + 0.5 / sqrtf(factor_radicand) * factor_radicand_partial_y))
        / sqr(factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        pseudo_chunk.x
            * factor_partial_y
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        (factor + pseudo_chunk.y * factor_partial_y)
            * (coshf(pseudo_chunk.z) + sinhf(pseudo_chunk.z) * tanhf(boost)),
        sinhf(pseudo_chunk.z) * (-scaled_factor_delta_partial_y * sqr(tanhf(boost)))
            - scaled_factor_delta_partial_y * tanhf(boost) * coshf(pseudo_chunk.z),
        coshf(pseudo_chunk.z) * (scaled_factor_delta_partial_y)
            + scaled_factor_delta_partial_y * tanhf(boost) * sinhf(pseudo_chunk.z),
    )
}

fn voxel_to_mvector_boosted_partial_z(pseudo_chunk: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(pseudo_chunk.x) + sqr(pseudo_chunk.y);
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        pseudo_chunk.x * factor * (sinhf(pseudo_chunk.z) + coshf(pseudo_chunk.z) * tanhf(boost)),
        pseudo_chunk.y * factor * (sinhf(pseudo_chunk.z) + coshf(pseudo_chunk.z) * tanhf(boost)),
        coshf(pseudo_chunk.z) * (1.0 - scaled_factor_delta * sqr(tanhf(boost)))
            - scaled_factor_delta * tanhf(boost) * sinhf(pseudo_chunk.z),
        sinhf(pseudo_chunk.z) * (1.0 + scaled_factor_delta)
            + scaled_factor_delta * tanhf(boost) * coshf(pseudo_chunk.z),
    )
}

fn add_quad(
    chunk: &BltChunk,
    layout: &BltLayout,
    transform: &MIsometry<f32>,
    geometry: &mut MeshGeometryDefinition,
    points: [na::Vector3<i32>; 4],
    texture: usize,
) {
    let vertices: Vec<_> = points
        .into_iter()
        .enumerate()
        .map(|(i, point)| {
            let len = geometry.vertices.len();
            geometry.vertices.push(Vertex {
                position: common::dodeca::Side::A.reflection()
                    * common::dodeca::Vertex::A.dual_to_node()
                    * (transform * chunk.point_from_voxel(layout, point.cast())).tuv_to_xyz(1),
                texcoords: na::Vector3::new((i & 1) as f32, ((i >> 1) & 1) as f32, texture as f32),
                normal: common::math::MDirection::x(),
            });
            len as u32
        })
        .collect();
    geometry.indices.extend(&[
        vertices[0],
        vertices[1],
        vertices[2],
        vertices[1],
        vertices[3],
        vertices[2],
    ]);
}

fn add_voxel(
    chunk: &BltChunk,
    layout: &BltLayout,
    transform: &MIsometry<f32>,
    geometry: &mut MeshGeometryDefinition,
    coords: na::Vector3<i32>,
) {
    for x_axis in 0..3 {
        let t = na::Vector3::x().tuv_to_xyz(x_axis);
        let u = na::Vector3::y().tuv_to_xyz(x_axis);
        let v = na::Vector3::z().tuv_to_xyz(x_axis);
        add_quad(
            chunk,
            layout,
            transform,
            geometry,
            [coords, coords + t, coords + u, coords + t + u],
            x_axis,
        );
        add_quad(
            chunk,
            layout,
            transform,
            geometry,
            [
                coords + v,
                coords + u + v,
                coords + t + v,
                coords + t + u + v,
            ],
            x_axis,
        );
    }
}

pub struct SampleSurface {
    pub geometry: MeshGeometryDefinition,
}

impl SampleSurface {
    pub fn new() -> Self {
        let mut graph = BltGraph::new();
        let mut current_chunk = graph.root_chunk;
        let mut current_transform = MIsometry::<f32>::identity();

        let mut geometry = MeshGeometryDefinition {
            vertices: Vec::new(),
            indices: Vec::new(),
        };
        for k in 0..2 {
            for x in (0..(graph.layout.horizontal_size as i32)).step_by(2) {
                for y in (0..(graph.layout.horizontal_size as i32)).step_by(2) {
                    for z in (0..(graph.layout.outer_vertical_size as i32)).step_by(2) {
                        add_voxel(
                            graph.chunk(current_chunk),
                            &graph.layout,
                            &current_transform,
                            &mut geometry,
                            na::Vector3::new(x, y, z),
                        );
                    }
                }
            }
            let index = 3;
            current_transform *= graph
                .chunk(current_chunk)
                .outer_isometry(&graph.layout, index);
            current_chunk = graph.add_outer(current_chunk, index);
        }
        SampleSurface { geometry }
    }
}

impl skid_steer::Source for SampleSurface {
    type Output = Mesh;

    async fn load<'a>(self, context: &'a skid_steer::Context<'a>) -> Option<Mesh> {
        let ctx: &AssetLoadContext = context.get().unwrap();

        let colors = ctx.load_cached(crate::graphics::PngArray {
            path: "materials".into(),
            size: common::world::Material::COUNT - 1,
        });

        Mesh::from_definition(ctx, self.geometry, colors).await
    }

    fn free(mut output: Self::Output, context: &skid_steer::Context) {
        let ctx: &AssetLoadContext = context.get().unwrap();
        unsafe { output.destroy(ctx.device()) };
    }
}

struct BltGraph {
    chunks: Vec<BltChunk>,
    root_chunk: u32,
    layout: BltLayout,
}

impl BltGraph {
    fn new() -> Self {
        BltGraph {
            chunks: vec![BltChunk::new_central()],
            root_chunk: 0,
            layout: BltLayout::default(),
        }
    }

    fn add_outer(&mut self, inner_chunk: u32, index: u8) -> u32 {
        let outer_chunk = self.chunks.len() as u32;
        let inner = &mut self.chunks[inner_chunk as usize];
        let mut outer = inner.new_outer(&self.layout, index);
        inner.outer_neighbors[index as usize] = Some(outer_chunk);
        outer.inner_neighbor = Some(inner_chunk);
        self.chunks.push(outer);
        outer_chunk
    }

    fn chunk(&self, chunk: u32) -> &BltChunk {
        &self.chunks[chunk as usize]
    }
}

struct BltLayout {
    horizontal_size: u8,
    central_vertical_size: u8,
    outer_vertical_size: u8,
    central_voxel_width: f32,
    voxel_height: f32,
}

impl Default for BltLayout {
    fn default() -> Self {
        Self {
            horizontal_size: 11,
            central_vertical_size: 11,
            outer_vertical_size: 11,
            central_voxel_width: 0.7 / 11.0,
            voxel_height: logf(2.0) / 11.0,
        }
    }
}

struct BltChunk {
    inner_neighbor: Option<u32>,
    inner_neighbor_index: u8,
    outer_neighbors: [Option<u32>; 4],
    klein_coords: na::Vector2<f32>,
    voxel_coords_conversion: na::Matrix3<f32>,
    boost: f32,
}

impl BltChunk {
    fn new_central() -> Self {
        BltChunk {
            inner_neighbor: None,
            inner_neighbor_index: 0,
            outer_neighbors: [None; 4],
            klein_coords: na::Vector2::zeros(),
            voxel_coords_conversion: na::Matrix3::identity(),
            boost: 0.0,
        }
    }

    fn point_from_voxel(&self, layout: &BltLayout, voxel_coords: na::Vector3<f32>) -> MPoint<f32> {
        self.point_from_chunk(na::Vector3::new(
            voxel_coords[0] * layout.central_voxel_width,
            voxel_coords[1] * layout.central_voxel_width,
            voxel_coords[2] * layout.voxel_height,
        ))
    }

    fn point_from_chunk(&self, chunk_coords: na::Vector3<f32>) -> MPoint<f32> {
        let horizontal_coords = self.voxel_coords_conversion * chunk_coords.xy().push(1.0);
        voxel_to_mvector_boosted(
            na::Vector3::new(
                horizontal_coords[0] / horizontal_coords[2],
                horizontal_coords[1] / horizontal_coords[2],
                chunk_coords.z,
            ),
            self.boost,
        )
        .to_point_unchecked()
    }

    fn outer_isometry(&self, layout: &BltLayout, index: u8) -> MIsometry<f32> {
        let scale = layout.central_voxel_width * layout.horizontal_size as f32 * 0.5;
        let chunk_pos = na::Vector3::new(
            scale * (index & 1) as f32,
            scale * (index >> 1) as f32,
            layout.voxel_height * layout.outer_vertical_size as f32,
        );
        self.isometry_from_chunk(chunk_pos)
    }

    fn isometry_from_chunk(&self, chunk: na::Vector3<f32>) -> MIsometry<f32> {
        let horizontal_coords = self.voxel_coords_conversion * chunk.xy().push(1.0);
        let pseudo_chunk = na::Vector3::new(
            horizontal_coords[0] / horizontal_coords[2],
            horizontal_coords[1] / horizontal_coords[2],
            chunk.z,
        );
        let w = voxel_to_mvector_boosted(pseudo_chunk, self.boost).to_point_unchecked();
        let x = voxel_to_mvector_boosted_partial_x(pseudo_chunk, self.boost).normalized_direction();
        let z =
            voxel_to_mvector_boosted_partial_z(pseudo_chunk, self.boost).to_direction_unchecked();
        // TODO: There might be a better way to get `y`, especially since voxel_to_mvector_boosted_partial_y will be
        // in the wrong direction.
        let mut y = voxel_to_mvector_boosted_partial_y(pseudo_chunk, self.boost);
        y -= MVector::from(x) * y.mip(&x);
        let y = y.normalized_direction();
        let result = MIsometry::from_columns_unchecked(&[x, y, z], w);
        println!("sanity check: {:?}", result.inverse() * result);
        println!("result: {:?}", result);
        result
    }

    fn new_outer(&self, layout: &BltLayout, index: u8) -> Self {
        // TODO: Support non-zero `index`
        let new_boost = self.boost + layout.voxel_height * layout.outer_vertical_size as f32;
        let scale_factor = coshf(new_boost) / coshf(self.boost);
        let scale = layout.central_voxel_width * layout.horizontal_size as f32 * 0.5;
        let chunk_pos = na::Vector3::new(
            scale * (index & 1) as f32,
            scale * (index >> 1) as f32,
            layout.voxel_height * layout.outer_vertical_size as f32,
        );
        // Find a set of m-orthogonal coordinates that moves (0,0,1) to (chunk_pos.x,chunk_pos.y,1), with the x-coordinate aligned
        let origin = na::Vector3::new(chunk_pos[0], chunk_pos[1], 1.0);
        let origin_norm = sqrtf(-(sqr(origin[0]) + sqr(origin[1]) - sqr(origin[2])));
        let normalized_origin = origin / origin_norm;
        let z = normalized_origin;
        let mut y = z.cross(&na::Vector3::x());
        y.z *= -1.0;
        y /= sqrtf(sqr(y.x) + sqr(y.y) - sqr(y.z));
        let mut x = y.cross(&z);
        x.z *= -1.0;
        let mut conversion = na::Matrix3::from_columns(&[x, y, z]).try_inverse().unwrap();
        conversion.set_column(2, &na::Vector3::new(0.0, 0.0, 0.9395)); // Applying skew
        println!("conversion: {:?}", conversion);
        BltChunk {
            inner_neighbor: None,
            inner_neighbor_index: index,
            outer_neighbors: [None; 4],
            klein_coords: na::Vector2::zeros(),
            voxel_coords_conversion: self.voxel_coords_conversion
                * conversion
                * na::Matrix3::new_scaling(scale_factor * 0.5),
            boost: new_boost,
        }
    }
}

#[cfg(test)]
mod tests {
    use common::math::MIsometry;

    use super::*;

    #[test]
    fn example() {
        let example = na::vector![0.1, 0.2, 2.0];
        let boost = 1.0;
        println!(
            "{:?}",
            MIsometry::translation_along(&(na::Vector3::z() * -boost))
                * voxel_to_mvector_simple(na::Vector3::new(
                    example.x / coshf(boost),
                    example.y / coshf(boost),
                    boost + example.z
                ))
        );
        println!("{:?}", voxel_to_mvector_boosted(example, boost));

        println!(
            "{:?}",
            (voxel_to_mvector_boosted(example + na::Vector3::z() * 0.001, boost)
                - voxel_to_mvector_boosted(example + na::Vector3::z() * -0.001, boost))
                / 0.002
        );
        println!(
            "{:?}",
            (voxel_to_mvector_boosted(example + na::Vector3::z() * 0.0001, boost)
                - voxel_to_mvector_boosted(example + na::Vector3::z() * -0.0001, boost))
                / 0.0002
        );
        println!("{:?}", voxel_to_mvector_boosted_partial_z(example, boost));
        BltChunk::new_central().isometry_from_chunk(example);
    }

    #[test]
    fn test_chunk_to_mvector_simple() {
        let klein_coords = na::Vector2::new(0.2, 0.3);
        let chunk = na::Vector3::new(0.25, 0.35, 0.0);
        let boost = 0.0;
        println!("{:?}", chunk_to_mvector_simple(klein_coords, chunk, boost));
        println!(
            "{:?}",
            voxel_to_mvector_boosted(chunk_to_pseudo_chunk(klein_coords, chunk, boost), boost)
        );
    }
}
