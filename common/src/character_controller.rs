use crate::{
    graph::Graph,
    math,
    proto::{CharacterInput, Position},
    sanitize_motion_input, SimConfig,
};

pub fn run_character_step<T>(
    config: &SimConfig,
    graph: &Graph<T>,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    input: &CharacterInput,
    dt_seconds: f32,
) {
    CharacterControllerPass {
        config,
        graph,
        position,
        velocity,
        input,
        dt_seconds,
    }
    .step();
}

struct CharacterControllerPass<'a, T> {
    config: &'a SimConfig,
    graph: &'a Graph<T>,
    position: &'a mut Position,
    velocity: &'a mut na::Vector3<f32>,
    input: &'a CharacterInput,
    dt_seconds: f32,
}

impl<T> CharacterControllerPass<'_, T> {
    fn step(&mut self) {
        let movement = sanitize_motion_input(self.input.movement);

        if self.input.no_clip {
            *self.velocity = na::Vector3::zeros();
            self.position.local *= math::translate_along(
                &(movement * self.config.no_clip_movement_speed * self.dt_seconds),
            );
        } else {
            let old_velocity = *self.velocity;

            // Update velocity
            let current_to_target_velocity =
                movement * self.config.max_ground_speed - *self.velocity;
            let max_delta_velocity = self.config.ground_acceleration * self.dt_seconds;
            if current_to_target_velocity.norm_squared() > max_delta_velocity.powi(2) {
                *self.velocity += current_to_target_velocity.normalize() * max_delta_velocity;
            } else {
                *self.velocity += current_to_target_velocity;
            }

            // Update position using the average between the old and new velocity to avoid the following two issues:
            // 1. Input lag, which would occur if only the old velocity was used
            // 2. Discontinuities, which would occur if only the new velocity was used
            self.position.local *=
                math::translate_along(&((*self.velocity + old_velocity) * 0.5 * self.dt_seconds));
        }

        // Renormalize
        self.position.local = math::renormalize_isometry(&self.position.local);
        let (next_node, transition_xf) = self
            .graph
            .normalize_transform(self.position.node, &self.position.local);
        if next_node != self.position.node {
            self.position.node = next_node;
            self.position.local = transition_xf * self.position.local;
        }
    }
}
