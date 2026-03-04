use common::math::{MVector, sqr};
use libm::{coshf, sinhf, sqrtf, tanhf};

fn voxel_to_mvector_simple(voxel: na::Vector3<f32>) -> MVector<f32> {
    let factor = sqrtf(1.0 - sqr(voxel.x) - sqr(voxel.y));
    MVector::new(voxel.x, voxel.y, factor * tanhf(voxel.z), 1.0)
}

fn voxel_to_mvector_boosted(voxel: na::Vector3<f32>, boost: f32) -> MVector<f32> {
    let f = 1.0 - sqrtf(1.0 - sqr(voxel.x) - sqr(voxel.y));
    let g = f * sqr(coshf(boost));
    let u = voxel.z + boost;
    println!("{}", g);
    MVector::new(
        voxel.x * coshf(boost),
        voxel.y * coshf(boost),
        sqr(coshf(boost)) * tanhf(-boost + u) + coshf(boost) * sinhf(boost),
        //(coshf(u) * tanhf(boost) * g + sinhf(u) * (1.0 - g)) / (coshf(u) - tanhf(boost) * sinhf(u)),
        0.0, // TODO
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn example() {
        let example = na::vector![0.0, 0.0, 2.0];
        let boost = 1.0;
        // I'll have to think about this later.
        println!("{:?}", voxel_to_mvector_simple(example) * coshf(boost));
        println!("{:?}", voxel_to_mvector_boosted(example, boost));
    }
}
