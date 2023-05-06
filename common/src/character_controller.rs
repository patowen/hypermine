use rand_distr::num_traits::Zero;
use tracing::warn;

use crate::{
    character_controller::{
        bound_vector::{BoundVector, VectorBound},
        collision::{check_collision, CollisionCheckingResult, CollisionContext},
    },
    math,
    node::{ChunkLayout, DualGraph},
    proto::{CharacterInput, Position},
    sanitize_motion_input, SimConfig,
};

pub fn run_character_step(
    cfg: &SimConfig,
    graph: &DualGraph,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    on_ground: &mut bool,
    input: &CharacterInput,
    dt_seconds: f32,
) {
    let movement = sanitize_motion_input(input.movement);

    if input.no_clip {
        // If no-clip is on, the velocity field is useless, and we don't want to accidentally
        // save velocity from when no-clip was off.
        *velocity = na::Vector3::zeros();
        position.local *=
            math::translate_along(&(movement * cfg.no_clip_movement_speed * dt_seconds));
    } else {
        let collision_context = CollisionContext {
            graph,
            chunk_layout: ChunkLayout::new(cfg.chunk_size as usize),
            radius: cfg.character_radius,
        };

        let up = get_relative_up(graph, position);
        let max_slope = cfg.max_floor_slope;

        // Initialize ground_normal
        let mut ground_normal = None;
        if *on_ground {
            let ground_result = get_ground_normal(
                &collision_context,
                &up,
                max_slope,
                cfg.ground_distance_tolerance,
                position,
            );
            if let Some(ground_collision) = ground_result.collision {
                ground_normal = Some(ground_collision.normal);
            }
        }

        // Jump if appropriate
        if input.jump && ground_normal.is_some() {
            let horizontal_velocity = *velocity - *up * up.dot(velocity);
            *velocity = horizontal_velocity + *up * cfg.jump_speed;
            ground_normal = None;
        }

        let old_velocity = *velocity;

        // Update velocity
        if let Some(ground_normal) = ground_normal {
            apply_ground_controls(
                cfg.ground_acceleration,
                cfg.max_ground_speed,
                dt_seconds,
                &movement,
                &up,
                &ground_normal,
                velocity,
            );
        } else {
            apply_air_controls(cfg.air_acceleration, dt_seconds, &movement, velocity);

            // Apply air resistance
            *velocity *= (-cfg.air_resistance * dt_seconds).exp();
        }

        // Apply gravity
        *velocity -= *up * cfg.gravity_acceleration * dt_seconds;

        // Apply speed cap
        *velocity = velocity.cap_magnitude(cfg.speed_cap);

        // Estimate the average velocity by using the average of the old velocity and new velocity,
        // which has the effect of modeling a velocity that changes linearly over the timestep.
        // This is necessary to avoid the following two issues:
        // 1. Input lag, which would occur if only the old velocity was used
        // 2. Movement artifacts, which would occur if only the new velocity was used. One
        //    example of such an artifact is the character moving backwards slightly when they
        //    stop moving after releasing a direction key.
        let estimated_average_velocity = (*velocity + old_velocity) * 0.5;

        apply_velocity(
            &collision_context,
            &up,
            max_slope,
            estimated_average_velocity * dt_seconds,
            position,
            velocity,
            &mut ground_normal,
        );

        *on_ground = ground_normal.is_some();
    }

    // Renormalize
    position.local = math::renormalize_isometry(&position.local);
    let (next_node, transition_xf) = graph.normalize_transform(position.node, &position.local);
    if next_node != position.node {
        position.node = next_node;
        position.local = transition_xf * position.local;
    }
}

fn apply_ground_controls(
    ground_acceleration: f32,
    max_ground_speed: f32,
    dt_seconds: f32,
    movement: &na::Vector3<f32>,
    up: &na::UnitVector3<f32>,
    ground_normal: &na::UnitVector3<f32>,
    velocity: &mut na::Vector3<f32>,
) {
    let movement_norm = movement.norm();
    let target_ground_velocity = if movement_norm < 1e-16 {
        na::Vector3::zeros()
    } else {
        let mut unit_movement = movement / movement_norm;
        math::project_to_plane(&mut unit_movement, ground_normal, up, 0.0);
        unit_movement.try_normalize_mut(1e-16);
        unit_movement * movement_norm * max_ground_speed
    };
    let mut ground_velocity = *velocity;
    math::project_to_plane(&mut ground_velocity, ground_normal, up, 0.0);
    let current_to_target_velocity = target_ground_velocity - ground_velocity;
    let max_delta_velocity = ground_acceleration * dt_seconds;
    if current_to_target_velocity.norm_squared() > max_delta_velocity.powi(2) {
        *velocity += current_to_target_velocity.normalize() * max_delta_velocity;
    } else {
        *velocity += current_to_target_velocity;
    }
}

fn apply_air_controls(
    air_acceleration: f32,
    dt_seconds: f32,
    movement: &na::Vector3<f32>,
    velocity: &mut na::Vector3<f32>,
) {
    *velocity += movement * air_acceleration * dt_seconds;
}

/// Updates the position based on the given average velocity while handling collisions. Also updates the velocity
/// based on collisions that occur.
fn apply_velocity(
    collision_context: &CollisionContext,
    up: &na::UnitVector3<f32>,
    max_slope: f32,
    expected_displacement: na::Vector3<f32>,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    ground_normal: &mut Option<na::UnitVector3<f32>>,
) {
    // To prevent an unbounded runtime, we only allow a limited number of collisions to be processed in
    // a single step. If the player encounters excessively complex geometry, it is possible to hit this limit,
    // in which case further movement processing is delayed until the next time step.
    const MAX_COLLISION_ITERATIONS: u32 = 6;
    let min_slope_up_component = 1.0 / (max_slope.powi(2) + 1.0).sqrt();

    let mut remaining_displacement = BoundVector::new(expected_displacement, Some(*velocity));
    let mut vertical_correction_direction = BoundVector::new(-up.into_inner(), None);

    let mut all_collisions_resolved = false;
    for _ in 0..MAX_COLLISION_ITERATIONS {
        let collision_result = check_collision(
            collision_context,
            position,
            &remaining_displacement.inner.vector,
        );
        position.local *= collision_result.displacement_transform;

        if let Some(collision) = collision_result.collision {
            // Update the expected displacement to whatever is remaining.
            remaining_displacement.inner.vector -= collision_result.displacement_vector;

            let mut push_direction = &collision.normal;

            if collision.normal.dot(up) > min_slope_up_component {
                if let Some(ground_normal) = ground_normal {
                    if vertical_correction_direction.is_facing(ground_normal) {
                        vertical_correction_direction.inner.vector.normalize_mut();
                        let vertical_correction_bound = VectorBound::new_push(
                            *ground_normal,
                            na::UnitVector3::new_normalize(
                                vertical_correction_direction.inner.vector,
                            ),
                        );
                        remaining_displacement.apply(|v| v.apply_bound(&vertical_correction_bound));
                        vertical_correction_direction.inner.vector.set_zero();
                    }

                    remaining_displacement.apply(|v| {
                        apply_ground_normal_change(up, &collision.normal, &mut v.vector);
                    });
                }
                *ground_normal = Some(collision.normal);
                push_direction = up;
            }

            if let Some(ground_normal) = ground_normal {
                remaining_displacement
                    .add_temporary_bound(VectorBound::new_pull(*ground_normal, *up));
            }
            remaining_displacement
                .add_and_apply_bound(VectorBound::new_push(collision.normal, *push_direction));
            remaining_displacement.remove_temporary_bounds();

            if vertical_correction_direction.is_facing(&collision.normal) {
                vertical_correction_direction
                    .add_and_apply_bound(VectorBound::new_push(collision.normal, collision.normal));
            }
        } else {
            all_collisions_resolved = true;
            break;
        }
    }

    if !all_collisions_resolved {
        warn!("A character entity processed too many collisions and collision resolution was cut short.");
    }

    *velocity = remaining_displacement.tagalong.unwrap().vector;
}

fn get_ground_normal(
    collision_context: &CollisionContext,
    up: &na::UnitVector3<f32>,
    max_slope: f32,
    allowed_distance: f32,
    position: &Position,
) -> CollisionCheckingResult {
    const MAX_COLLISION_ITERATIONS: u32 = 6;
    let min_slope_up_component = 1.0 / (max_slope.powi(2) + 1.0).sqrt();

    let mut allowed_displacement = BoundVector::new(-up.into_inner() * allowed_distance, None);

    for _ in 0..MAX_COLLISION_ITERATIONS {
        let collision_result = check_collision(
            collision_context,
            position,
            &allowed_displacement.inner.vector,
        );
        if let Some(collision) = collision_result.collision.as_ref() {
            if collision.normal.dot(up) > min_slope_up_component {
                return collision_result;
            }
            allowed_displacement
                .add_and_apply_bound(VectorBound::new_push(collision.normal, collision.normal));
        } else {
            return CollisionCheckingResult::stationary();
        }
    }
    CollisionCheckingResult::stationary()
}

/// Returns the up-direction relative to the given position
fn get_relative_up(graph: &DualGraph, position: &Position) -> na::UnitVector3<f32> {
    na::UnitVector3::new_normalize(
        (math::mtranspose(&position.local)
            * graph
                .get(position.node)
                .as_ref()
                .unwrap()
                .state
                .up_direction())
        .xyz(),
    )
}

fn apply_ground_normal_change(
    up: &na::UnitVector3<f32>,
    new_ground_normal: &na::UnitVector3<f32>,
    subject: &mut na::Vector3<f32>,
) {
    // Try to keep the horizontal component constant. Update the vertical component to ensure
    // the subject is parallel to the ground, and scale the vector to undo any change in magnitude.
    let subject_norm = subject.norm();
    if subject_norm > 1e-16 {
        let mut unit_subject = *subject / subject_norm;
        math::project_to_plane(&mut unit_subject, new_ground_normal, up, 0.0);
        unit_subject.try_normalize_mut(1e-16);
        *subject = unit_subject * subject_norm;
    }
}

mod bound_vector {
    use rand_distr::num_traits::Zero;
    use tracing::warn;

    use crate::math;

    pub struct BoundVector {
        pub inner: VectorWithErrorMargin,
        pub tagalong: Option<VectorWithErrorMargin>,
        bounds: VectorBounds,
    }

    impl BoundVector {
        pub fn new(inner: na::Vector3<f32>, tagalong: Option<na::Vector3<f32>>) -> Self {
            BoundVector {
                inner: VectorWithErrorMargin::new(inner),
                tagalong: tagalong.map(VectorWithErrorMargin::new),
                bounds: VectorBounds {
                    permanent_bounds: vec![],
                    temporary_bounds: vec![],
                },
            }
        }

        pub fn apply(&mut self, mut f: impl FnMut(&mut VectorWithErrorMargin)) {
            f(&mut self.inner);
            if let Some(ref mut tagalong) = self.tagalong {
                f(tagalong);
            }
        }

        pub fn is_facing(&self, vector: &na::Vector3<f32>) -> bool {
            self.inner.vector.dot(vector) < 0.0
        }

        pub fn add_and_apply_bound(&mut self, new_bound: VectorBound) {
            self.apply_bound(&new_bound);
            self.bounds.permanent_bounds.push(new_bound);
        }

        fn apply_bound(&mut self, new_bound: &VectorBound) {
            if self.inner.vector.is_zero() {
                return;
            }

            self.apply(|v| v.apply_bound(new_bound));

            // Check if all constraints are satisfied
            if self.bounds.iter().all(|b| self.inner.check_bound(b)) {
                return;
            }

            // If not all constraints are satisfied, find the first constraint that if applied will satisfy
            // the remaining constriants
            for bound in self.bounds.iter().filter(|b| !self.inner.check_bound(b)) {
                let Some(ortho_bound) = bound.get_constrained_with_bound(new_bound)
                else {
                    warn!("Unsatisfied existing bound is parallel to new bound. Is the character squeezed between two walls?");
                    continue;
                };

                let mut candidate = self.inner.clone();
                candidate.apply_bound(&ortho_bound);

                if self.bounds.iter().all(|b| candidate.check_bound(b)) {
                    self.inner = candidate;
                    if let Some(ref mut tagalong) = self.tagalong {
                        tagalong.apply_bound(&ortho_bound);
                    }
                    return;
                }
            }

            // If no choice satisfies all constraints, keep all bounds and set the vector to 0
            self.apply(|v| v.vector.set_zero());
        }

        pub fn add_temporary_bound(&mut self, new_bound: VectorBound) {
            self.bounds.temporary_bounds.push(new_bound);
        }

        pub fn remove_temporary_bounds(&mut self) {
            self.bounds.temporary_bounds.clear();
        }
    }

    #[derive(Clone)]
    pub struct VectorWithErrorMargin {
        pub vector: na::Vector3<f32>,
        pub error_margin: f32,
    }

    impl VectorWithErrorMargin {
        fn new(vector: na::Vector3<f32>) -> Self {
            // Corrective term to ensure that normals face away from any potential collision surfaces
            const RELATIVE_EPSILON: f32 = 1e-4;

            let error_margin = vector.magnitude() * RELATIVE_EPSILON;
            VectorWithErrorMargin {
                vector,
                error_margin,
            }
        }

        pub fn apply_bound(&mut self, bound: &VectorBound) {
            math::project_to_plane(
                &mut self.vector,
                &bound.normal,
                &bound.push_direction,
                self.error_margin * bound.target_distance_factor,
            );
        }

        fn check_bound(&self, bound: &VectorBound) -> bool {
            self.vector.dot(&bound.normal) >= self.error_margin * bound.checked_distance_factor()
        }
    }

    struct VectorBounds {
        permanent_bounds: Vec<VectorBound>,
        temporary_bounds: Vec<VectorBound>,
    }

    impl VectorBounds {
        fn iter(&self) -> impl Iterator<Item = &VectorBound> {
            self.permanent_bounds.iter().chain(&self.temporary_bounds)
        }
    }

    pub struct VectorBound {
        normal: na::UnitVector3<f32>,
        push_direction: na::UnitVector3<f32>,
        target_distance_factor: f32, // Margin of error when the bound is applied
    }

    impl VectorBound {
        pub fn new_push(
            normal: na::UnitVector3<f32>,
            push_direction: na::UnitVector3<f32>,
        ) -> Self {
            VectorBound {
                normal,
                push_direction,
                target_distance_factor: 1.0,
            }
        }

        pub fn new_pull(
            normal: na::UnitVector3<f32>,
            push_direction: na::UnitVector3<f32>,
        ) -> Self {
            VectorBound {
                normal: na::UnitVector3::new_unchecked(-normal.as_ref()),
                push_direction,
                target_distance_factor: -1.0,
            }
        }

        // An additional margin of error is needed when the bound is checked to ensure that an
        // applied bound always passes the check.
        fn checked_distance_factor(&self) -> f32 {
            self.target_distance_factor - 1.0
        }

        fn get_constrained_with_bound(&self, bound: &VectorBound) -> Option<VectorBound> {
            const MIN_ORTHO_NORM: f32 = 1e-5;

            let mut ortho_bound_push_direction = self.push_direction.into_inner();
            math::project_to_plane(
                &mut ortho_bound_push_direction,
                &bound.normal,
                &bound.push_direction,
                0.0,
            );

            na::UnitVector3::try_new(ortho_bound_push_direction, MIN_ORTHO_NORM).map(|d| {
                VectorBound {
                    normal: self.normal,
                    push_direction: d,
                    target_distance_factor: self.target_distance_factor,
                }
            })
        }
    }
}

mod collision {
    use tracing::error;

    use crate::{
        graph_collision, math,
        node::{ChunkLayout, DualGraph},
        proto::Position,
    };

    /// Checks for collisions when a character moves with a character-relative displacement vector of `relative_displacement`.
    pub fn check_collision(
        collision_context: &CollisionContext,
        position: &Position,
        relative_displacement: &na::Vector3<f32>,
    ) -> CollisionCheckingResult {
        // Split relative_displacement into its norm and a unit vector
        let relative_displacement = relative_displacement.to_homogeneous();
        let displacement_sqr = relative_displacement.norm_squared();
        if displacement_sqr < 1e-16 {
            // Fallback for if the displacement vector isn't large enough to reliably be normalized.
            // Any value that is sufficiently large compared to f32::MIN_POSITIVE should work as the cutoff.
            return CollisionCheckingResult::stationary();
        }

        let displacement_norm = displacement_sqr.sqrt();
        let displacement_normalized = relative_displacement / displacement_norm;

        let ray = graph_collision::Ray::new(math::origin(), displacement_normalized);
        let tanh_distance = displacement_norm.tanh();

        let cast_hit = graph_collision::sphere_cast(
            collision_context.radius,
            collision_context.graph,
            &collision_context.chunk_layout,
            position,
            &ray,
            tanh_distance,
        );

        let cast_hit = match cast_hit {
            Ok(r) => r,
            Err(e) => {
                error!("Collision checking returned {:?}", e);
                return CollisionCheckingResult::stationary();
            }
        };

        let distance = cast_hit
            .as_ref()
            .map_or(tanh_distance, |hit| hit.tanh_distance)
            .atanh();

        let displacement_vector = displacement_normalized.xyz() * distance;
        let displacement_transform = math::translate_along(&displacement_vector);

        CollisionCheckingResult {
            displacement_vector,
            displacement_transform,
            collision: cast_hit.map(|hit| Collision {
                // `CastEndpoint` has its `normal` given relative to the character's original position,
                // but we want the normal relative to the character after the character moves to meet the wall.
                // This normal now represents a contact point at the origin, so we omit the w-coordinate
                // to ensure that it's orthogonal to the origin.
                normal: na::UnitVector3::new_normalize(
                    (math::mtranspose(&displacement_transform) * hit.normal).xyz(),
                ),
            }),
        }
    }

    pub struct CollisionContext<'a> {
        pub graph: &'a DualGraph,
        pub chunk_layout: ChunkLayout,
        pub radius: f32,
    }

    pub struct CollisionCheckingResult {
        /// The displacement allowed for the character before hitting a wall. The result of
        /// `math::translate_along(&displacement_vector)` is `displacement_transform`.
        pub displacement_vector: na::Vector3<f32>,

        /// Multiplying the character's position by this matrix will move the character as far as it can up to its intended
        /// displacement until it hits the wall.
        pub displacement_transform: na::Matrix4<f32>,

        pub collision: Option<Collision>,
    }

    impl CollisionCheckingResult {
        /// Return a CollisionCheckingResult with no movement and no collision; useful if the character is not moving
        /// and has nothing to check collision against. Also useful as a last resort fallback if an unexpected error occurs.
        pub fn stationary() -> CollisionCheckingResult {
            CollisionCheckingResult {
                displacement_vector: na::Vector3::zeros(),
                displacement_transform: na::Matrix4::identity(),
                collision: None,
            }
        }
    }

    pub struct Collision {
        /// This collision normal faces away from the collision surface and is given in the perspective of the character
        /// _after_ it is transformed by `allowed_displacement`. The 4th coordinate of this normal vector is assumed to be
        /// 0.0 and is therefore omitted.
        pub normal: na::UnitVector3<f32>,
    }
}
