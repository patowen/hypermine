use serde::{Deserialize, Serialize};

#[derive(
    Debug, Copy, Clone, Default, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize,
)]
#[repr(transparent)]
pub struct Material(u16);

pub mod materials {
    use super::Material;

    pub const VOID: Material = Material(0);
    pub const DIRT: Material = Material(1);
    pub const SAND: Material = Material(2);
    pub const SILT: Material = Material(3);
    pub const CLAY: Material = Material(4);
    pub const MUD: Material = Material(5);
    pub const SANDY_LOAM: Material = Material(6);
    pub const SILTY_LOAM: Material = Material(7);
    pub const CLAY_LOAM: Material = Material(8);
    pub const RED_SAND: Material = Material(9);
    pub const LIMESTONE: Material = Material(10);
    pub const SHALE: Material = Material(11);
    pub const DOLOMITE: Material = Material(12);
    pub const SANDSTONE: Material = Material(13);
    pub const RED_SANDSTONE: Material = Material(14);
    pub const MARBLE: Material = Material(15);
    pub const SLATE: Material = Material(16);
    pub const GRANITE: Material = Material(17);
    pub const DIORITE: Material = Material(18);
    pub const ANDESITE: Material = Material(19);
    pub const GABBRO: Material = Material(20);
    pub const BASALT: Material = Material(21);
    pub const OLIVINE: Material = Material(22);
    pub const WATER: Material = Material(23);
    pub const LAVA: Material = Material(24);
    pub const WOOD: Material = Material(25);
    pub const LEAVES: Material = Material(26);
    pub const WOOD_PLANKS: Material = Material(27);
    pub const GREY_BRICK: Material = Material(28);
    pub const WHITE_BRICK: Material = Material(29);
    pub const ICE: Material = Material(30);
    pub const ICE_SLUSH: Material = Material(31);
    pub const GRAVEL: Material = Material(32);
    pub const SNOW: Material = Material(33);
    pub const COARSE_GRASS: Material = Material(34);
    pub const TAN_GRASS: Material = Material(35);
    pub const LUSH_GRASS: Material = Material(36);
    pub const MUD_GRASS: Material = Material(37);
    pub const GRASS: Material = Material(38);
    pub const CAVE_GRASS: Material = Material(39);
}

impl Material {
    pub const COUNT: usize = 40;
}
