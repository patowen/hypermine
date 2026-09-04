use common::math::{MVector, PermuteXYZ, sqr};
use libm::{coshf, logf, powf, sinhf, sqrtf, tanhf};

use crate::graphics::{
    Mesh,
    asset_loader::AssetLoadContext,
    meshes::{MeshGeometryDefinition, Vertex},
};

fn voxel_to_mvector_simple(voxel: na::Vector3<f32>) -> MVector<f32> {
    let factor = sqrtf(1.0 - sqr(voxel.x) - sqr(voxel.y));
    MVector::new(voxel.x, voxel.y, factor * tanhf(voxel.z), 1.0)
}

fn voxel_to_mvector_boosted(voxel: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    let f = 1.0 - sqrtf(1.0 - sqr(voxel.x) - sqr(voxel.y));
    let g = f * sqr(coshf(boost));
    let u = voxel.z + boost;
    MVector::new(
        voxel.x * coshf(boost),
        voxel.y * coshf(boost),
        (coshf(u) * tanhf(boost) * g + sinhf(u) * (1.0 - g)) / (coshf(u) - tanhf(boost) * sinhf(u)),
        0.0, // TODO
    )
}

fn coords_to_mvector(coords: na::Vector3<i32>, width_factor: f32) -> MVector<f32> {
    voxel_to_mvector_simple(na::Vector3::new(
        coords[0] as f32 * 0.04 * width_factor,
        coords[1] as f32 * 0.04 * width_factor,
        coords[2] as f32 * logf(2.0) / 20.0,
    ))
}

fn add_quad(
    geometry: &mut MeshGeometryDefinition,
    points: [na::Vector3<i32>; 4],
    width_factor: f32,
) {
    let vertices: Vec<_> = points
        .into_iter()
        .enumerate()
        .map(|(i, point)| {
            let len = geometry.vertices.len();
            geometry.vertices.push(Vertex {
                position: common::dodeca::Vertex::A.dual_to_node()
                    * coords_to_mvector(point, width_factor)
                        .normalized_point()
                        .tuv_to_xyz(1),
                texcoords: na::Vector3::new((i & 1) as f32, ((i >> 1) & 1) as f32, 0.0),
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
        for x in -3..=3 {
            for y in -3..=3 {
                for z in 0..=10 {
                    for k in 0..6 {
                        add_voxel(
                            &mut geometry,
                            na::Vector3::new(x * 2, y * 2, -z * 2 - 20 * k),
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

#[cfg(test)]
mod tests {
    use common::math::MIsometry;

    use super::*;

    #[test]
    fn example() {
        let example = na::vector![0.1, 0.2, 2.0];
        let boost = 1.0;
        // I'll have to think about this later.
        println!(
            "{:?}",
            MIsometry::translation_along(&(na::Vector3::z() * boost))
                * voxel_to_mvector_simple(example)
                * coshf(boost)
        );
        println!("{:?}", voxel_to_mvector_boosted(example, boost));
    }
}
