use common::math::{MPoint, MVector, PermuteXYZ, sqr};
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

/// Computes
/// `translation_along([0, 0, -boost]) * voxel_to_mvector_simple([x / cosh(boost), y / cosh(boost), boost + z])`
fn voxel_to_mvector_boosted(voxel: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(voxel.x) + sqr(voxel.y);
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        voxel.x * factor * (coshf(voxel.z) + sinhf(voxel.z) * tanhf(boost)),
        voxel.y * factor * (coshf(voxel.z) + sinhf(voxel.z) * tanhf(boost)),
        sinhf(voxel.z) * (1.0 - scaled_factor_delta * sqr(tanhf(boost)))
            - scaled_factor_delta * tanhf(boost) * coshf(voxel.z),
        coshf(voxel.z) * (1.0 + scaled_factor_delta)
            + scaled_factor_delta * tanhf(boost) * sinhf(voxel.z),
    )
}

fn voxel_to_mvector_boosted_partial_x(voxel: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(voxel.x) + sqr(voxel.y);
    let dist_squared_partial_x = voxel.x * 2.0;
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
        (factor + voxel.x * factor_partial_x) * (coshf(voxel.z) + sinhf(voxel.z) * tanhf(boost)),
        voxel.y * factor_partial_x * (coshf(voxel.z) + sinhf(voxel.z) * tanhf(boost)),
        sinhf(voxel.z) * (-scaled_factor_delta_partial_x * sqr(tanhf(boost)))
            - scaled_factor_delta_partial_x * tanhf(boost) * coshf(voxel.z),
        coshf(voxel.z) * (scaled_factor_delta_partial_x)
            + scaled_factor_delta_partial_x * tanhf(boost) * sinhf(voxel.z),
    )
}

fn voxel_to_mvector_boosted_partial_y(voxel: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(voxel.x) + sqr(voxel.y);
    let dist_squared_partial_y = voxel.y * 2.0;
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
        voxel.x * factor_partial_y * (coshf(voxel.z) + sinhf(voxel.z) * tanhf(boost)),
        (factor + voxel.y * factor_partial_y) * (coshf(voxel.z) + sinhf(voxel.z) * tanhf(boost)),
        sinhf(voxel.z) * (-scaled_factor_delta_partial_y * sqr(tanhf(boost)))
            - scaled_factor_delta_partial_y * tanhf(boost) * coshf(voxel.z),
        coshf(voxel.z) * (scaled_factor_delta_partial_y)
            + scaled_factor_delta_partial_y * tanhf(boost) * sinhf(voxel.z),
    )
}

fn voxel_to_mvector_boosted_partial_z(voxel: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    // `factor = 1 + scaled_factor_delta/cosh(boost)^2`
    let dist_squared = sqr(voxel.x) + sqr(voxel.y);
    let factor_radicand = 1.0 - dist_squared / sqr(coshf(boost));
    let factor = 1.0 / sqrtf(factor_radicand);
    let scaled_factor_delta = dist_squared / (factor_radicand + sqrtf(factor_radicand));
    MVector::new(
        voxel.x * factor * (sinhf(voxel.z) + coshf(voxel.z) * tanhf(boost)),
        voxel.y * factor * (sinhf(voxel.z) + coshf(voxel.z) * tanhf(boost)),
        coshf(voxel.z) * (1.0 - scaled_factor_delta * sqr(tanhf(boost)))
            - scaled_factor_delta * tanhf(boost) * sinhf(voxel.z),
        sinhf(voxel.z) * (1.0 + scaled_factor_delta)
            + scaled_factor_delta * tanhf(boost) * coshf(voxel.z),
    )
}

fn coords_to_mvector(coords: na::Vector3<i32>, width_factor: f32) -> MVector<f32> {
    voxel_to_mvector_boosted(
        na::Vector3::new(
            coords[0] as f32 * 0.04 * width_factor,
            coords[1] as f32 * 0.04 * width_factor,
            coords[2] as f32 * logf(2.0) / 20.0,
        ),
        0.0,
    )
}

fn add_quad(
    geometry: &mut MeshGeometryDefinition,
    points: [na::Vector3<i32>; 4],
    width_factor: f32,
    texture: usize,
) {
    let vertices: Vec<_> = points
        .into_iter()
        .enumerate()
        .map(|(i, point)| {
            let len = geometry.vertices.len();
            geometry.vertices.push(Vertex {
                position: common::dodeca::Vertex::A.dual_to_node()
                    * coords_to_mvector(point, width_factor)
                        .to_point_unchecked()
                        .tuv_to_xyz(1),
                texcoords: na::Vector3::new((i & 1) as f32, ((i >> 1) & 1) as f32, texture as f32),
                normal: common::math::MDirection::x(),
            });
            len as u32
        })
        .collect();
    geometry.indices.extend(&[
        vertices[0],
        vertices[2],
        vertices[1],
        vertices[1],
        vertices[2],
        vertices[3],
    ]);
}

fn add_voxel(geometry: &mut MeshGeometryDefinition, coords: na::Vector3<i32>, width_factor: f32) {
    for x_axis in 0..3 {
        let t = na::Vector3::x().tuv_to_xyz(x_axis);
        let u = na::Vector3::y().tuv_to_xyz(x_axis);
        let v = na::Vector3::z().tuv_to_xyz(x_axis);
        add_quad(
            geometry,
            [coords, coords + t, coords + u, coords + t + u],
            width_factor,
            x_axis,
        );
        add_quad(
            geometry,
            [
                coords + v,
                coords + u + v,
                coords + t + v,
                coords + t + u + v,
            ],
            width_factor,
            x_axis,
        );
    }
}

pub struct SampleSurface {
    pub geometry: MeshGeometryDefinition,
}

impl SampleSurface {
    pub fn new() -> Self {
        let mut geometry = MeshGeometryDefinition {
            vertices: Vec::new(),
            indices: Vec::new(),
        };
        for x in -6..=7 {
            for y in -6..=7 {
                for z in 0..10 {
                    for k in 0..5 {
                        add_voxel(
                            &mut geometry,
                            na::Vector3::new(
                                x * 2 + (2i32.pow(k as u32) - 1) * 15,
                                y * 2 + (2i32.pow(k as u32) - 1) * 15,
                                -1 - z * 2 - 20 * k,
                            ),
                            powf(0.5, k as f32),
                        );
                    }
                }
            }
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

    fn add_outer(&mut self, inner_chunk: u32, index: u8) {
        let outer_chunk = self.chunks.len() as u32;
        let inner = &mut self.chunks[inner_chunk as usize];
        let mut outer = inner.new_outer(&self.layout, index);
        inner.outer_neighbors[index as usize] = Some(outer_chunk);
        outer.inner_neighbor = Some(inner_chunk);
        self.chunks.push(outer);
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
            horizontal_size: 12,
            central_vertical_size: 12,
            outer_vertical_size: 12,
            central_voxel_width: 0.7 / 12.0,
            voxel_height: logf(2.0) / 12.0,
        }
    }
}

struct BltChunk {
    inner_neighbor: Option<u32>,
    inner_neighbor_index: u8,
    outer_neighbors: [Option<u32>; 4],
    voxel_coords_conversion: na::Matrix3<f32>,
    boost: f32,
}

impl BltChunk {
    fn new_central() -> Self {
        BltChunk {
            inner_neighbor: None,
            inner_neighbor_index: 0,
            outer_neighbors: [None; 4],
            voxel_coords_conversion: na::Matrix3::identity(),
            boost: 0.0,
        }
    }

    fn point_from_voxel(&self, voxel: na::Vector3<f32>) -> MPoint<f32> {
        let horizontal_coords = self.voxel_coords_conversion * voxel.xy().push(1.0);
        voxel_to_mvector_boosted(
            na::Vector3::new(
                horizontal_coords[0] / horizontal_coords[2],
                horizontal_coords[1] / horizontal_coords[2],
                voxel.z,
            ),
            self.boost,
        )
        .to_point_unchecked()
    }

    fn new_outer(&self, layout: &BltLayout, index: u8) -> Self {
        if index != 0 {
            unimplemented!();
        }
        BltChunk {
            inner_neighbor: None,
            inner_neighbor_index: 0,
            outer_neighbors: [None; 4],
            voxel_coords_conversion: na::Matrix3::identity(),
            boost: self.boost + layout.voxel_height * layout.outer_vertical_size as f32,
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
    }
}
