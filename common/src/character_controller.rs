use std::mem::replace;

use tracing::warn;

use crate::{
    character_controller::{
        collision::{check_collision, Collision, CollisionContext},
        vector_bounds::{VectorBound, VectorBoundGroup},
    },
    math,
    node::{ChunkLayout, DualGraph},
    proto::{CharacterInput, Position},
    sanitize_motion_input, SimConfig,
};

/// Runs a single step of character movement
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
        *velocity = movement * cfg.no_clip_movement_speed;
        position.local *= math::translate_along(&(*velocity * dt_seconds));
    } else {
        // Initialize current state
        let collision_context = CollisionContext {
            graph,
            chunk_layout: ChunkLayout::new(cfg.chunk_size as usize),
            radius: cfg.character_radius,
        };

        let up = graph.get_relative_up(position).unwrap();

        let mut ground_normal = None;
        if *on_ground {
            ground_normal = get_ground_normal(
                &collision_context,
                &up,
                cfg.max_ground_slope,
                cfg.ground_distance_tolerance,
                position,
            );
        }

        // Handle jumping
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
        let average_velocity = (*velocity + old_velocity) * 0.5;

        // Handle actual movement
        apply_velocity(
            &collision_context,
            &up,
            cfg.max_ground_slope,
            average_velocity,
            dt_seconds,
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

/// Returns the normal corresponding to the ground below the character, up to the `allowed_distance`. If
/// no such ground exists, returns `None`.
fn get_ground_normal(
    collision_context: &CollisionContext,
    up: &na::UnitVector3<f32>,
    max_slope: f32,
    allowed_distance: f32,
    position: &Position,
) -> Option<na::UnitVector3<f32>> {
    // Since the character can be at a corner between a slanted wall and the ground, the first collision
    // directly below the character is not guaranteed to be part of the ground regardless of whether the
    // character is on the ground. To handle this, we repeatedly redirect the direction we search to be
    // parallel to walls we collide with to ensure that we find the ground if is indeed below the character.
    const MAX_COLLISION_ITERATIONS: u32 = 6;
    let mut allowed_displacement = -up.into_inner() * allowed_distance;
    let mut bounds = VectorBoundGroup::new(&allowed_displacement);

    for _ in 0..MAX_COLLISION_ITERATIONS {
        let collision_result = check_collision(collision_context, position, &allowed_displacement);
        if let Some(collision) = collision_result.collision.as_ref() {
            if is_ground(up, max_slope, &collision.normal) {
                // We found the ground, so return its normal.
                return Some(collision.normal);
            }
            bounds.apply_and_add_bound(
                VectorBound::new_push(collision.normal, collision.normal),
                &[],
                &mut allowed_displacement,
                None,
            );
        } else {
            // Return `None` if we travel the whole `allowed_displacement` and don't find the ground.
            return None;
        }
    }
    // Return `None` if we fail to find the ground after the maximum number of attempts
    None
}

/// Checks whether the given normal is flat enough to be considered part of the ground
fn is_ground(up: &na::UnitVector3<f32>, max_slope: f32, normal: &na::UnitVector3<f32>) -> bool {
    let min_slope_up_component = 1.0 / (max_slope.powi(2) + 1.0).sqrt();
    normal.dot(up) > min_slope_up_component
}

/// Updates the velocity based on user input assuming the character is on the ground
fn apply_ground_controls(
    ground_acceleration: f32,
    max_ground_speed: f32,
    dt_seconds: f32,
    movement: &na::Vector3<f32>,
    up: &na::UnitVector3<f32>,
    ground_normal: &na::UnitVector3<f32>,
    velocity: &mut na::Vector3<f32>,
) {
    // Set `target_ground_velocity` to have a consistent magnitude regardless
    // of the movement direction, but ensure that the horizontal direction matches
    // the horizontal direction of the intended movement direction.
    let movement_norm = movement.norm();
    let target_ground_velocity = if movement_norm < 1e-16 {
        na::Vector3::zeros()
    } else {
        let mut unit_movement = movement / movement_norm;
        math::project_to_plane(&mut unit_movement, ground_normal, up, 0.0);
        unit_movement.try_normalize_mut(1e-16);
        unit_movement * movement_norm * max_ground_speed
    };

    // Set `ground_velocity` to be the current velocity's ground-parallel component,
    // using a basis that contains the up vector to ensure that the result is unaffected
    // by gravity.
    let mut ground_velocity = *velocity;
    math::project_to_plane(&mut ground_velocity, ground_normal, up, 0.0);

    // Adjust the ground-parallel component of the velocity vector to be closer to the
    // target velocity.
    let current_to_target_velocity = target_ground_velocity - ground_velocity;
    let max_delta_velocity = ground_acceleration * dt_seconds;
    if current_to_target_velocity.norm_squared() > max_delta_velocity.powi(2) {
        *velocity += current_to_target_velocity.normalize() * max_delta_velocity;
    } else {
        *velocity += current_to_target_velocity;
    }
}

/// Updates the velocity based on user input assuming the character is in the air
fn apply_air_controls(
    air_acceleration: f32,
    dt_seconds: f32,
    movement: &na::Vector3<f32>,
    velocity: &mut na::Vector3<f32>,
) {
    *velocity += movement * air_acceleration * dt_seconds;
}

/// Updates the character's position based on the given average velocity while handling collisions.
/// Also updates the velocity and ground normal based on collisions that occur.
#[allow(clippy::too_many_arguments)] // TODO: Reduce argument count (Peer review feedback needed on how best to organize this)
fn apply_velocity(
    collision_context: &CollisionContext,
    up: &na::UnitVector3<f32>,
    max_slope: f32,
    average_velocity: na::Vector3<f32>,
    dt_seconds: f32,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    ground_normal: &mut Option<na::UnitVector3<f32>>,
) {
    // To prevent an unbounded runtime, we only allow a limited number of collisions to be processed in
    // a single step. If the character encounters excessively complex geometry, it is possible to hit this limit,
    // in which case further movement processing is delayed until the next time step.
    const MAX_COLLISION_ITERATIONS: u32 = 6;

    let mut remaining_dt_seconds = dt_seconds;

    let initial_velocity_info = VelocityInfo {
        bounds: VectorBoundGroup::new(&average_velocity),
        average_velocity,
        final_velocity: *velocity,
    };
    let mut velocity_info = initial_velocity_info.clone();

    let mut ground_collision_handled = false;

    let mut all_collisions_resolved = false;
    for _ in 0..MAX_COLLISION_ITERATIONS {
        let expected_displacement = velocity_info.average_velocity * remaining_dt_seconds;

        let collision_result = check_collision(collision_context, position, &expected_displacement);
        position.local *= collision_result.displacement_transform;

        if let Some(collision) = collision_result.collision {
            // Update the expected dt to whatever is remaining.
            remaining_dt_seconds *= 1.0
                - collision_result.displacement_vector.magnitude()
                    / expected_displacement.magnitude();

            handle_collision(
                collision,
                up,
                max_slope,
                &initial_velocity_info,
                &mut velocity_info,
                ground_normal,
                &mut ground_collision_handled,
            );
        } else {
            all_collisions_resolved = true;
            break;
        }
    }

    if !all_collisions_resolved {
        warn!("A character entity processed too many collisions and collision resolution was cut short.");
    }

    *velocity = velocity_info.final_velocity;
}

/// Updates character information based on the results of a single collision
fn handle_collision(
    collision: Collision,
    up: &na::UnitVector3<f32>,
    max_slope: f32,
    initial_velocity_info: &VelocityInfo,
    velocity_info: &mut VelocityInfo,
    ground_normal: &mut Option<na::UnitVector3<f32>>,
    ground_collision_handled: &mut bool,
) {
    // Collisions are divided into two categories: Ground collisions and wall collisions.
    // Ground collisions will only affect vertical movement of the character, while wall collisions will
    // push the character away from the wall in a perpendicular direction. If the character is on the ground,
    // we have extra logic to ensure that slanted wall collisions do not lift the character off the ground.
    if is_ground(up, max_slope, &collision.normal) {
        let stay_on_ground_bounds = [VectorBound::new_pull(collision.normal, *up)];
        if !*ground_collision_handled {
            // Wall collisions can turn vertical momentum into unwanted horizontal momentum. This can
            // occur if the character jumps at the corner between the ground and a slanted wall. If the wall
            // collision is handled first, this horizontal momentum will push the character away from the wall.
            // This can also occur if the character is on the ground and walks into a slanted wall. A single frame
            // of downward momentum caused by gravity can turn into unwanted horizontal momentum that pushes
            // the character away from the wall. Neither of these issues can occur if the ground collision is
            // handled first, so when computing how the velocity vectors change, we rewrite history as if
            // the ground collision was first. This is only necessary for the first ground collision, since
            // afterwards, there is no more unexpected vertical momentum.
            let old_velocity_info = replace(velocity_info, initial_velocity_info.clone());
            velocity_info.bounds.apply_and_add_bound(
                VectorBound::new_push(collision.normal, *up),
                &stay_on_ground_bounds,
                &mut velocity_info.average_velocity,
                Some(&mut velocity_info.final_velocity),
            );
            for bound in old_velocity_info.bounds.bounds() {
                velocity_info.bounds.apply_and_add_bound(
                    bound.clone(),
                    &stay_on_ground_bounds,
                    &mut velocity_info.average_velocity,
                    Some(&mut velocity_info.final_velocity),
                );
            }

            *ground_collision_handled = true;
        } else {
            velocity_info.bounds.apply_and_add_bound(
                VectorBound::new_push(collision.normal, *up),
                &stay_on_ground_bounds,
                &mut velocity_info.average_velocity,
                Some(&mut velocity_info.final_velocity),
            );
        }

        *ground_normal = Some(collision.normal);
    } else {
        let mut stay_on_ground_bounds = Vec::new();
        if let Some(ground_normal) = ground_normal {
            stay_on_ground_bounds.push(VectorBound::new_pull(*ground_normal, *up));
        }
        velocity_info.bounds.apply_and_add_bound(
            VectorBound::new_push(collision.normal, collision.normal),
            &stay_on_ground_bounds,
            &mut velocity_info.average_velocity,
            Some(&mut velocity_info.final_velocity),
        );
    }
}

/// Contains info related to the average velocity over the timestep and the current velocity at
/// the end of the timestep.
#[derive(Clone)]
struct VelocityInfo {
    bounds: VectorBoundGroup,
    average_velocity: na::Vector3<f32>,
    final_velocity: na::Vector3<f32>,
}

/// This module is used to transform vectors to ensure that they fit constraints discovered during collision checking.
mod vector_bounds {
    use rand_distr::num_traits::Zero;
    use tracing::warn;

    use crate::math;

    /// Encapsulates all the information needed to constrain a vector based on a set of `VectorBound`s.
    #[derive(Clone)]
    pub struct VectorBoundGroup {
        bounds: Vec<VectorBound>,
        error_margin: f32,
    }

    impl VectorBoundGroup {
        /// Initializes a `VectorBoundGroup` with an empty list of bounds. The `initial_vector` is the first vector
        /// we expect these bounds to be applied to, a hint to determine what kind of error margin is needed
        /// to prevent floating point approximation limits from causing phantom collisions. Note that this
        /// error margin is not needed if the resulting vector is zero, since no phantom collision can occur
        /// if the character is stopped.
        pub fn new(initial_vector: &na::Vector3<f32>) -> Self {
            let error_margin = initial_vector.magnitude() * 1e-4;

            VectorBoundGroup {
                bounds: vec![],
                error_margin,
            }
        }

        /// Returns the internal list of `VectorBound`s contained in the `VectorBoundGroup` struct.
        pub fn bounds(&self) -> &[VectorBound] {
            &self.bounds
        }

        /// Constrains `vector` with `new_bound` while keeping the existing constraints and any constraints in
        /// `temporary_bounds` satisfied. All projection transformations applied to `vector` are also applied
        /// to `tagalong` to allow two vectors to be transformed consistently with each other.
        pub fn apply_and_add_bound(
            &mut self,
            new_bound: VectorBound,
            temporary_bounds: &[VectorBound],
            vector: &mut na::Vector3<f32>,
            tagalong: Option<&mut na::Vector3<f32>>,
        ) {
            self.apply_bound(&new_bound, temporary_bounds, vector, tagalong);
            self.bounds.push(new_bound);
        }

        /// Helper function to logically separate the "add" and the "apply" in `apply_and_add_bound` function.
        fn apply_bound(
            &self,
            new_bound: &VectorBound,
            temporary_bounds: &[VectorBound],
            vector: &mut na::Vector3<f32>,
            mut tagalong: Option<&mut na::Vector3<f32>>,
        ) {
            // There likely isn't a perfect way to get a vector properly constrained with a list of bounds. The main
            // difficulty is finding which set of linearly independent bounds need to be applied so that all bounds are
            // satisfied. Since bounds are one-sided and not guaranteed to be linearly independent from each other, this
            // requires some ad-hoc choices. The algorithm we choose here is to (1) assume that `new_bound` is one of these
            // linearly independent bounds, (2) if necessary, pair it up with each existing bound to find the first such
            // bound that allows all bounds to be satisfied, and (3) zero out the vector if no such pairing works, as we
            // assume that we need to apply three linearly independent bounds.

            // Combine existing bounds with temporary bounds into an iterator
            let bounds_iter = self.bounds.iter().chain(temporary_bounds.iter());

            // Apply new_bound if necessary.
            if !new_bound.check_vector(vector, self.error_margin) {
                new_bound.constrain_vector(vector, self.error_margin);
                if let Some(ref mut tagalong) = tagalong {
                    // Note: The tagalong vector does not need an error margin.
                    new_bound.constrain_vector(tagalong, 0.0);
                }
            }

            // Check if all constraints are satisfied
            if (bounds_iter.clone()).all(|b| b.check_vector(vector, self.error_margin)) {
                return;
            }

            // If not all constraints are satisfied, find the first constraint that if applied will satisfy
            // the remaining constriants
            for bound in
                (bounds_iter.clone()).filter(|b| !b.check_vector(vector, self.error_margin))
            {
                let Some(ortho_bound) = bound.get_self_constrained_with_bound(new_bound)
                else {
                    warn!("Unsatisfied existing bound is parallel to new bound. Is the character squeezed between two walls?");
                    continue;
                };

                let mut candidate = *vector;
                ortho_bound.constrain_vector(&mut candidate, self.error_margin);

                if (bounds_iter.clone()).all(|b| b.check_vector(&candidate, self.error_margin)) {
                    *vector = candidate;
                    if let Some(ref mut tagalong) = tagalong {
                        ortho_bound.constrain_vector(tagalong, 0.0);
                    }
                    return;
                }
            }

            // If no choice satisfies all constraints, keep all bounds and set the vector to 0
            vector.set_zero();
            if let Some(ref mut tagalong) = tagalong {
                tagalong.set_zero();
            }
        }
    }

    /// Represents a single constraint for a vector. `VectorBound`s alone conceptually contain
    /// enough information to apply to a vector, but practically, one other piece of information
    /// is needed: `error_margin`, which exists in `VectorBoundGroup`.
    #[derive(Clone)]
    pub struct VectorBound {
        normal: na::UnitVector3<f32>,
        projection_direction: na::UnitVector3<f32>,
        error_margin_factor: f32, // Margin of error when the bound is applied
    }

    impl VectorBound {
        /// Creates a `VectorBound` that pushes vectors away from the plane given
        /// by the normal in `projection_direction`. After applying such a bound to
        /// a vector, its dot product with `normal` should be positive even counting
        /// floating point approximation limitations.
        pub fn new_push(
            normal: na::UnitVector3<f32>,
            projection_direction: na::UnitVector3<f32>,
        ) -> Self {
            VectorBound {
                normal,
                projection_direction,
                error_margin_factor: 1.0,
            }
        }

        /// Creates a `VectorBound` that pulls vectors towards the plane given
        /// by the normal in `projection_direction`. Even after applying such a bound to
        /// a vector, its dot product with `normal` should still be positive even counting
        /// floating point approximation limitations. This ensures that `new_push` and
        /// `new_pull` don't conflict with each other even with equal parameters.
        pub fn new_pull(
            normal: na::UnitVector3<f32>,
            projection_direction: na::UnitVector3<f32>,
        ) -> Self {
            VectorBound {
                normal: na::UnitVector3::new_unchecked(-normal.as_ref()),
                projection_direction,
                error_margin_factor: -1.0,
            }
        }

        /// Updates `subject` with a projection transformation based on the constraint given by `self`.
        /// This function does not check whether such a constraint is needed.
        fn constrain_vector(&self, subject: &mut na::Vector3<f32>, error_margin: f32) {
            math::project_to_plane(
                subject,
                &self.normal,
                &self.projection_direction,
                error_margin * self.error_margin_factor,
            );
        }

        /// Checks whether `subject` satisfies the constraint given by `self`.
        fn check_vector(&self, subject: &na::Vector3<f32>, error_margin: f32) -> bool {
            // An additional margin of error is needed when the bound is checked to ensure that an
            // applied bound always passes the check.
            let error_margin_factor_for_check = self.error_margin_factor - 0.5;
            subject.is_zero()
                || subject.dot(&self.normal) >= error_margin * error_margin_factor_for_check
        }

        /// Returns a `VectorBound` that is an altered version of `self` so that it no longer interferes
        /// with `bound`. This is achieved by altering the projection direction by a factor of
        /// `bound`'s projection direction to be orthogonal to `bound`'s normal. If this is not
        /// possible, returns `None`.
        fn get_self_constrained_with_bound(&self, bound: &VectorBound) -> Option<VectorBound> {
            let mut ortho_bound_projection_direction = self.projection_direction.into_inner();
            math::project_to_plane(
                &mut ortho_bound_projection_direction,
                &bound.normal,
                &bound.projection_direction,
                0.0,
            );

            na::UnitVector3::try_new(ortho_bound_projection_direction, 1e-5).map(|d| VectorBound {
                normal: self.normal,
                projection_direction: d,
                error_margin_factor: self.error_margin_factor,
            })
        }
    }
}

/// This module is used to encapsulate character collision checking
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

    /// Contains information about the character and the world that is only relevant for collision checking
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
