use std::time::Duration;

use serde::{Deserialize, Serialize};

use crate::{dodeca, math::MVector};

/// Manually specified simulation config parameters
#[derive(Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct SimConfigRaw {
    /// Number of steps per second
    pub rate: Option<u16>,
    /// Maximum distance at which nodes will be rendered in meters
    pub view_distance: Option<f32>,
    /// Maximum distance at which new chunks will be generated in meters
    pub chunk_generation_distance: Option<f32>,
    /// Distance at which fog becomes completely opaque in meters
    pub fog_distance: Option<f32>,
    pub input_queue_size_ms: Option<u16>,
    /// Whether gameplay-like restrictions exist, such as limited inventory
    pub gameplay_enabled: Option<bool>,
    /// Number of voxels along the edge of a chunk
    pub chunk_size: Option<u8>,
    /// Approximate length of the edge of a voxel in meters
    ///
    /// Curved spaces have a notion of absolute distance, defined with respect to the curvature. We
    /// mostly work in those units, but the conversion from meters, used to scale the player and assets
    /// and so forth, is configurable. This effectively allows for configurable curvature.
    ///
    /// Note that exact voxel size varies within each chunk. We reference the mean width of the voxels
    /// along the X axis through the center of a chunk.
    pub voxel_size: Option<f32>,
    /// Static configuration information relevant to character physics
    #[serde(default)]
    pub character: CharacterConfigRaw,
}

/// Complete simulation config parameters
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SimConfig {
    /// Amount of time between each step. Inverse of the rate
    pub step_interval: Duration,
    /// Maximum distance at which nodes will be rendered in absolute units
    pub view_distance: f32,
    /// Maximum distance at which new chunks will be generate in absolute units
    pub chunk_generation_distance: f32,
    /// Distance at which fog becomes completely opaque in absolute units
    pub fog_distance: f32,
    pub input_queue_size: Duration,
    /// Whether gameplay-like restrictions exist, such as limited inventory
    pub gameplay_enabled: bool,
    /// Number of voxels along the edge of a chunk
    pub chunk_size: u8,
    /// Static configuration information relevant to character physics
    pub character: CharacterConfig,
    /// Scaling factor converting meters to absolute units
    pub meters_to_absolute: f32,
}

impl SimConfig {
    pub fn from_raw(x: &SimConfigRaw) -> Self {
        let chunk_size = x.chunk_size.unwrap_or(12);
        let voxel_size = x.voxel_size.unwrap_or(1.0);
        let meters_to_absolute = meters_to_absolute(chunk_size, voxel_size);
        SimConfig {
            step_interval: Duration::from_secs(1) / x.rate.unwrap_or(30) as u32,
            view_distance: x.view_distance.unwrap_or(75.0) * meters_to_absolute,
            chunk_generation_distance: x.chunk_generation_distance.unwrap_or(60.0)
                * meters_to_absolute,
            fog_distance: x.fog_distance.unwrap_or(90.0) * meters_to_absolute,
            input_queue_size: Duration::from_millis(x.input_queue_size_ms.unwrap_or(50).into()),
            gameplay_enabled: x.gameplay_enabled.unwrap_or(false),
            chunk_size,
            character: CharacterConfig::from_raw(&x.character, meters_to_absolute),
            meters_to_absolute,
        }
    }
}

/// Compute the scaling factor from meters to absolute units, given the number of voxels in a chunk
/// and the approximate size of a voxel in meters.
fn meters_to_absolute(chunk_size: u8, voxel_size: f32) -> f32 {
    let a = MVector::from(dodeca::Vertex::A.chunk_to_node() * na::Vector4::new(1.0, 0.5, 0.5, 1.0))
        .normalized_point();
    let b = MVector::from(dodeca::Vertex::A.chunk_to_node() * na::Vector4::new(0.0, 0.5, 0.5, 1.0))
        .normalized_point();
    let minimum_chunk_face_separation = a.distance(&b);
    let absolute_voxel_size = minimum_chunk_face_separation / f32::from(chunk_size);
    absolute_voxel_size / voxel_size
}

/// Static configuration information relevant to character physics as provided in configuration files
#[derive(Default, Clone, Debug, Serialize, Deserialize)]
pub struct CharacterConfigRaw {
    /// Character movement speed in m/s during no-clip
    pub no_clip_movement_speed: Option<f32>,
    /// Character maximumum movement speed while on the ground in m/s
    pub max_ground_speed: Option<f32>,
    /// Character artificial speed cap to avoid overloading the server in m/s
    pub speed_cap: Option<f32>,
    /// Maximum ground slope (0=horizontal, 1=45 degrees)
    pub max_ground_slope: Option<f32>,
    /// Character acceleration while on the ground in m/s^2
    pub ground_acceleration: Option<f32>,
    /// Character acceleration while in the air in m/s^2
    pub air_acceleration: Option<f32>,
    /// Acceleration of gravity in m/s^2
    pub gravity_acceleration: Option<f32>,
    /// Air resistance in (m/s^2) per (m/s); scales linearly with respect to speed
    pub air_resistance: Option<f32>,
    /// How fast the player jumps off the ground in m/s
    pub jump_speed: Option<f32>,
    /// How far away the player needs to be from the ground in meters to be considered in the air in meters
    pub ground_distance_tolerance: Option<f32>,
    /// Radius of the character in meters
    pub character_radius: Option<f32>,
    /// How far a character can reach when placing blocks in meters
    pub block_reach: Option<f32>,
}

/// Static configuration information relevant to character physics. Most fields are based on
/// absolute units and seconds.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CharacterConfig {
    /// Character movement speed in units/s during no-clip
    pub no_clip_movement_speed: f32,
    /// Character maximumum movement speed while on the ground in units/s
    pub max_ground_speed: f32,
    /// Character artificial speed cap to avoid overloading the server in units/s
    pub speed_cap: f32,
    /// Maximum ground slope (0=horizontal, 1=45 degrees)
    pub max_ground_slope: f32,
    /// Character acceleration while on the ground in units/s^2
    pub ground_acceleration: f32,
    /// Character acceleration while in the air in units/s^2
    pub air_acceleration: f32,
    /// Acceleration of gravity in units/s^2
    pub gravity_acceleration: f32,
    /// Air resistance in (units/s^2) per (units/s); scales linearly with respect to speed
    pub air_resistance: f32,
    /// How fast the player jumps off the ground in units/s
    pub jump_speed: f32,
    /// How far away the player needs to be from the ground in meters to be considered in the air in absolute units
    pub ground_distance_tolerance: f32,
    /// Radius of the character in absolute units
    pub character_radius: f32,
    /// How far a character can reach when placing blocks in absolute units
    pub block_reach: f32,
}

impl CharacterConfig {
    pub fn from_raw(x: &CharacterConfigRaw, meters_to_absolute: f32) -> Self {
        CharacterConfig {
            no_clip_movement_speed: x.no_clip_movement_speed.unwrap_or(12.0) * meters_to_absolute,
            max_ground_speed: x.max_ground_speed.unwrap_or(4.0) * meters_to_absolute,
            speed_cap: x.speed_cap.unwrap_or(30.0) * meters_to_absolute,
            max_ground_slope: x.max_ground_slope.unwrap_or(1.73), // 60 degrees
            ground_acceleration: x.ground_acceleration.unwrap_or(20.0) * meters_to_absolute,
            air_acceleration: x.air_acceleration.unwrap_or(2.0) * meters_to_absolute,
            gravity_acceleration: x.gravity_acceleration.unwrap_or(20.0) * meters_to_absolute,
            air_resistance: x.air_resistance.unwrap_or(0.2),
            jump_speed: x.jump_speed.unwrap_or(8.0) * meters_to_absolute,
            ground_distance_tolerance: x.ground_distance_tolerance.unwrap_or(0.2)
                * meters_to_absolute,
            character_radius: x.character_radius.unwrap_or(0.4) * meters_to_absolute,
            block_reach: x.block_reach.unwrap_or(10.0) * meters_to_absolute,
        }
    }
}
