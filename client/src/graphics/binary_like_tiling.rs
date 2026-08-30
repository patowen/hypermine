use common::math::{MVector, sqr};
use libm::{coshf, sinhf, sqrtf, tanhf};

use crate::graphics::{Mesh, asset_loader::AssetLoadContext, meshes::MeshGeometryDefinition};

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

pub struct SampleSurface {
    pub geometry: MeshGeometryDefinition,
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
