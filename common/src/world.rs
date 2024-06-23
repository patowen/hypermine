use serde::{Deserialize, Serialize};

#[derive(
    Debug, Copy, Clone, Default, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize,
)]
#[repr(u16)]
pub enum Material {
    #[default]
    Void,
    Dirt,
    Sand,
    Silt,
    Clay,
    Mud,
    SandyLoam,
    SiltyLoam,
    ClayLoam,
    RedSand,
    Limestone,
    Shale,
    Dolomite,
    Sandstone,
    RedSandstone,
    Marble,
    Slate,
    Granite,
    Diorite,
    Andesite,
    Gabbro,
    Basalt,
    Olivine,
    Water,
    Lava,
    Wood,
    Leaves,
    WoodPlanks,
    GreyBrick,
    WhiteBrick,
    Ice,
    IceSlush,
    Gravel,
    Snow,
    CoarseGrass,
    TanGrass,
    LushGrass,
    MudGrass,
    Grass,
    CaveGrass,
}

impl Material {
    pub const COUNT: usize = 40;

    pub const ALL: [Self; Self::COUNT] = [
        Material::Void,
        Material::Dirt,
        Material::Sand,
        Material::Silt,
        Material::Clay,
        Material::Mud,
        Material::SandyLoam,
        Material::SiltyLoam,
        Material::ClayLoam,
        Material::RedSand,
        Material::Limestone,
        Material::Shale,
        Material::Dolomite,
        Material::Sandstone,
        Material::RedSandstone,
        Material::Marble,
        Material::Slate,
        Material::Granite,
        Material::Diorite,
        Material::Andesite,
        Material::Gabbro,
        Material::Basalt,
        Material::Olivine,
        Material::Water,
        Material::Lava,
        Material::Wood,
        Material::Leaves,
        Material::WoodPlanks,
        Material::GreyBrick,
        Material::WhiteBrick,
        Material::Ice,
        Material::IceSlush,
        Material::Gravel,
        Material::Snow,
        Material::CoarseGrass,
        Material::TanGrass,
        Material::LushGrass,
        Material::MudGrass,
        Material::Grass,
        Material::CaveGrass,
    ];
}

impl TryFrom<u16> for Material {
    type Error = ();

    fn try_from(value: u16) -> Result<Self, Self::Error> {
        Material::ALL.get(value as usize).ok_or(()).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::Material;

    #[test]
    fn u16_to_material_consistency_check() {
        for i in 0..Material::COUNT {
            let index = u16::try_from(i).unwrap();
            let material =
                Material::try_from(index).expect("no missing entries in try_from match statement");
            assert_eq!(index, material as u16);
        }
    }
}
