use crate::{
    dodeca::Vertex,
    math,
    node::DualGraph,
    proto::{CharacterInput, Position},
    ray_tracing, sanitize_motion_input,
    sphere_collider::SphereCollider,
    SimConfig,
};

pub fn run_character_step(
    cfg: &SimConfig,
    graph: &DualGraph,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    input: &CharacterInput,
    dt_seconds: f32,
) {
    CharacterControllerPass {
        cfg,
        graph,
        position,
        velocity,
        input,
        dt_seconds,
    }
    .step();
}

struct CharacterControllerPass<'a> {
    cfg: &'a SimConfig,
    graph: &'a DualGraph,
    position: &'a mut Position,
    velocity: &'a mut na::Vector3<f32>,
    input: &'a CharacterInput,
    dt_seconds: f32,
}

impl CharacterControllerPass<'_> {
    fn step(&mut self) {
        let movement = sanitize_motion_input(self.input.movement);

        if self.input.no_clip {
            // If no-clip is on, the velocity field is useless, and we don't want to accidentally
            // save velocity from when no-clip was off.
            *self.velocity = na::Vector3::zeros();
            self.position.local *= math::translate_along(
                &(movement * self.cfg.no_clip_movement_speed * self.dt_seconds),
            );
        } else {
            let old_velocity = *self.velocity;

            // Update velocity
            let current_to_target_velocity = movement * self.cfg.max_ground_speed - *self.velocity;
            let max_delta_velocity = self.cfg.ground_acceleration * self.dt_seconds;
            if current_to_target_velocity.norm_squared() > math::sqr(max_delta_velocity) {
                *self.velocity += current_to_target_velocity.normalize() * max_delta_velocity;
            } else {
                *self.velocity += current_to_target_velocity;
            }

            // Update position by using the average of the old velocity and new velocity, which has
            // the effect of modeling a velocity that changes linearly over the timestep. This is
            // necessary to avoid the following two issues:
            // 1. Input lag, which would occur if only the old velocity was used
            // 2. Movement artifacts, which would occur if only the new velocity was used. One
            //    example of such an artifact is the player moving backwards slightly when they
            //    stop moving after releasing a direction key.
            self.position.local *=
                self.trace_ray(&((*self.velocity + old_velocity) * 0.5 * self.dt_seconds));
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

    fn trace_ray(&self, relative_displacement: &na::Vector3<f32>) -> na::Matrix4<f32> {
        let relative_displacement = relative_displacement.to_homogeneous();
        let displacement_sqr = relative_displacement.norm_squared();
        if displacement_sqr < 1e-16 {
            return na::Matrix4::identity();
        }

        let displacement_norm = displacement_sqr.sqrt();
        let displacement_normalized = relative_displacement / displacement_norm;

        let ray_status = ray_tracing::trace_ray(
            self.graph,
            self.cfg.chunk_size as usize,
            &SphereCollider { radius: 0.02 },
            (self.position.node, Vertex::A).into(),
            self.position.local,
            ray_tracing::Ray::new(
                self.position.local * math::origin(),
                self.position.local * displacement_normalized,
            ),
            displacement_norm.tanh(),
        );

        if let ray_tracing::RayTracingResult::Inconclusive = ray_status.result {
            return na::Matrix4::identity();
        }

        math::translate(
            &math::origin(),
            &math::lorentz_normalize(
                &(math::origin() + displacement_normalized * ray_status.tanh_length),
            ),
        )
    }
}
