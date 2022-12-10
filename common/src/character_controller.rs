use std::ops::Mul;

use crate::{
    math,
    proto::{Character, Position},
    sanitize_motion_input, SimConfig,
};

struct CharacterInput {
    movement: na::Vector3<f32>,
    orientation: na::UnitQuaternion<f32>,
    attempt_jump: bool,
    no_clip: bool,
}

struct CharacterControllerPass<'a> {
    position: &'a mut Position,
    character: &'a mut Character,
    input: &'a CharacterInput,
    config: &'a SimConfig,
    dt_seconds: f32,
}

impl CharacterControllerPass<'_> {
    fn step(&mut self) {
        let (direction, speed) = sanitize_motion_input(self.input.movement);

        if self.input.no_clip {
            self.character.velocity = na::Vector3::zeros();
            self.position.local *= math::translate_along(
                &direction,
                speed * self.config.no_clip_movement_speed * self.dt_seconds,
            );
        } else {
            // Update velocity
            let current_to_target_velocity = direction.mul(speed) - self.character.velocity;
            let max_delta_velocity = self.config.ground_acceleration * self.dt_seconds;
            if current_to_target_velocity.norm_squared() > max_delta_velocity.powi(2) {
                self.character.velocity +=
                    current_to_target_velocity.normalize() * max_delta_velocity;
            } else {
                self.character.velocity += current_to_target_velocity;
            }

            // Update position
            let (mut velocity_direction, velocity_speed) =
                na::Unit::new_and_get(self.character.velocity);
            if speed == 0.0 {
                // Use an arbitrary direction rather than NaN
                velocity_direction = -na::Vector3::z_axis();
            }

            self.position.local *=
                math::translate_along(&velocity_direction, velocity_speed * self.dt_seconds);
        }
    }
}
