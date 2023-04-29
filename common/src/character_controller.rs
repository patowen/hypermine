use tracing::{error, warn};

use crate::{
    graph_collision, math,
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
        let max_slope_angle = cfg.max_floor_slope_angle;

        // Initialize ground_normal
        let mut ground_normal = None;
        if *on_ground {
            let ground_result = get_ground_normal(
                &collision_context,
                &up,
                max_slope_angle,
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
            max_slope_angle,
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
    ground_normal: &na::Vector3<f32>,
    velocity: &mut na::Vector3<f32>,
) {
    let movement_norm = movement.norm();
    let target_velocity = if movement_norm < 1e-16 {
        na::Vector3::zeros()
    } else {
        let mut unit_movement = movement / movement_norm;
        let upward_correction = -unit_movement.dot(ground_normal) / up.dot(ground_normal);
        unit_movement += **up * upward_correction;
        unit_movement.try_normalize_mut(1e-16);
        unit_movement * movement_norm
    };
    let vertical_component = velocity.dot(ground_normal) / up.dot(ground_normal);
    let current_to_target_velocity =
        target_velocity * max_ground_speed - (*velocity - **up * vertical_component);
    let max_delta_velocity = ground_acceleration * dt_seconds;
    if current_to_target_velocity.norm_squared() > math::sqr(max_delta_velocity) {
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
    max_slope_angle: f32,
    mut expected_displacement: na::Vector3<f32>,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    ground_normal: &mut Option<na::UnitVector3<f32>>,
) {
    // To prevent an unbounded runtime, we only allow a limited number of collisions to be processed in
    // a single step. If the player encounters excessively complex geometry, it is possible to hit this limit,
    // in which case further movement processing is delayed until the next time step.
    const MAX_COLLISION_ITERATIONS: u32 = 5;
    let cos_max_slope = max_slope_angle.cos();

    let mut active_normals = Vec::<na::UnitVector3<f32>>::new();
    let mut active_normals_horizontal = Vec::<na::UnitVector3<f32>>::new();

    let mut expected_displacement_horizontal = expected_displacement;
    let mut velocity_horizontal = *velocity;

    if let Some(ground_normal) = ground_normal {
        expected_displacement_horizontal -=
            **up * (expected_displacement.dot(ground_normal) / up.dot(ground_normal));
        velocity_horizontal -= **up * (velocity.dot(ground_normal) / up.dot(ground_normal));
    }

    let mut all_collisions_resolved = false;
    for _ in 0..MAX_COLLISION_ITERATIONS {
        let collision_result = check_collision(collision_context, position, &expected_displacement);
        position.local *= collision_result.displacement_transform;

        if let Some(collision) = collision_result.collision {
            // Update the expected displacement to whatever is remaining.
            let displacement_factor = 1.0
                - collision_result.displacement_vector.magnitude()
                    / expected_displacement.magnitude();
            expected_displacement *= displacement_factor;
            expected_displacement_horizontal *= displacement_factor;

            if collision.normal.dot(up) > cos_max_slope {
                if ground_normal.is_some() {
                    expected_displacement = expected_displacement_horizontal;
                    *velocity = velocity_horizontal;
                }
                apply_ground_normal_change(
                    up,
                    ground_normal.is_some(),
                    &collision.normal,
                    &mut expected_displacement,
                );
                apply_ground_normal_change(
                    up,
                    ground_normal.is_some(),
                    &collision.normal,
                    velocity,
                );
                *ground_normal = Some(collision.normal);
                expected_displacement_horizontal = expected_displacement;
                velocity_horizontal = *velocity;
                active_normals.retain(|n| n.dot(velocity) < 0.0);
                active_normals_horizontal.retain(|n| n.dot(velocity) < 0.0);
            } else {
                active_normals.retain(|n| n.dot(&collision.normal) < 0.0);
                active_normals_horizontal.retain(|n| n.dot(&collision.normal) < 0.0);
            }

            active_normals.push(collision.normal);
            if collision.normal.dot(&expected_displacement_horizontal) < 0.0 {
                active_normals_horizontal.push(collision.normal);
            }

            apply_normals(&active_normals, &mut expected_displacement);
            apply_normals(&active_normals, velocity);

            if let Some(ground_normal) = ground_normal {
                let mut active_normals_horizontal_with_ground = active_normals_horizontal.clone();
                active_normals_horizontal_with_ground.push(*ground_normal);
                apply_normals(
                    &active_normals_horizontal_with_ground,
                    &mut expected_displacement_horizontal,
                );
                apply_normals(
                    &active_normals_horizontal_with_ground,
                    &mut velocity_horizontal,
                );

                if velocity.dot(ground_normal) > 0.0 {
                    expected_displacement = expected_displacement_horizontal;
                    *velocity = velocity_horizontal;
                }
            }
        } else {
            all_collisions_resolved = true;
            break;
        }
    }

    if !all_collisions_resolved {
        warn!("A character entity processed too many collisions and collision resolution was cut short.");
    }
}

fn get_ground_normal(
    collision_context: &CollisionContext,
    up: &na::UnitVector3<f32>,
    max_slope_angle: f32,
    allowed_distance: f32,
    position: &Position,
) -> CollisionCheckingResult {
    const MAX_COLLISION_ITERATIONS: u32 = 5;
    let cos_max_slope = max_slope_angle.cos();

    let mut allowed_displacement = -up.into_inner() * allowed_distance;
    let mut active_normals = Vec::<na::UnitVector3<f32>>::with_capacity(3);

    for _ in 0..MAX_COLLISION_ITERATIONS {
        let collision_result = check_collision(collision_context, position, &allowed_displacement);
        if let Some(collision) = collision_result.collision.as_ref() {
            if collision.normal.dot(up) > cos_max_slope {
                return collision_result;
            }
            active_normals.retain(|n| n.dot(&collision.normal) < 0.0);
            active_normals.push(collision.normal);
            apply_normals(&active_normals, &mut allowed_displacement);
        } else {
            return CollisionCheckingResult::stationary();
        }
    }
    CollisionCheckingResult::stationary()
}

/// Checks for collisions when a character moves with a character-relative displacement vector of `relative_displacement`.
fn check_collision(
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

/// Modifies the `subject` by a linear combination of the `normals` to ensure that it is approximately
/// orthogonal to all the normals. The normals are assumed to be linearly independent, and, assuming the final
/// result is nonzero, a small correction is applied to ensure that the subject is moving away from the surfaces
/// the normals represent even when floating point approximation is involved.
fn apply_normals(normals: &[na::UnitVector3<f32>], subject: &mut na::Vector3<f32>) {
    if normals.len() >= 3 {
        // The normals are assumed to be linearly independent, so applying all of them will zero out the subject.
        // There is no need to do any extra logic to handle precision limitations in this case.
        *subject = na::Vector3::zeros();
    }

    // Corrective term to ensure that normals face away from any potential collision surfaces
    const RELATIVE_EPSILON: f32 = 1e-4;
    apply_normals_internal(normals, subject, subject.magnitude() * RELATIVE_EPSILON);
}

/// Modifies the `subject` by a linear combination of the `normals` so that the dot product with each normal is
/// `distance`. The `normals` must be linearly independent for this function to work as expected.
fn apply_normals_internal(
    normals: &[na::UnitVector3<f32>],
    subject: &mut na::Vector3<f32>,
    distance: f32,
) {
    let mut ortho_normals: Vec<na::Vector3<f32>> = normals.iter().map(|n| n.into_inner()).collect();
    for i in 0..normals.len() {
        // Perform the Gram-Schmidt process on `normals` to produce `ortho_normals`.
        for j in i + 1..normals.len() {
            ortho_normals[j] = (ortho_normals[j]
                - ortho_normals[i] * ortho_normals[j].dot(&ortho_normals[i]))
            .normalize();
        }

        // The following formula ensures that the dot product of `subject` and `normals[i]` is `distance`.
        // Because we only move the subject along `ortho_normals[i]`, this adjustment does not affect the
        // subject's dot product with any earlier normals.
        *subject += ortho_normals[i]
            * ((distance - subject.dot(&normals[i])) / ortho_normals[i].dot(&normals[i]));
    }
}

fn apply_ground_normal_change(
    up: &na::UnitVector3<f32>,
    was_on_ground: bool,
    new_ground_normal: &na::UnitVector3<f32>,
    subject: &mut na::Vector3<f32>,
) {
    if was_on_ground {
        let subject_norm = subject.norm();
        if subject_norm > 1e-16 {
            let mut unit_subject = *subject / subject_norm;
            let upward_correction =
                -unit_subject.dot(new_ground_normal) / up.dot(new_ground_normal);
            unit_subject += **up * upward_correction;
            unit_subject.try_normalize_mut(1e-16);
            *subject = unit_subject * subject_norm;
        }
    } else {
        // TODO: Consider using fancier formula for max_upward_correction, one that makes
        // new_ground_normal and subject as collinear as possible.
        let mut upward_correction = -subject.dot(new_ground_normal) / up.dot(new_ground_normal);
        let max_upward_correction = -subject.dot(up);
        if upward_correction > max_upward_correction {
            upward_correction = max_upward_correction;
        }
        if upward_correction >= 0.0 {
            *subject += **up * upward_correction;
        }
    }
}

struct CollisionContext<'a> {
    graph: &'a DualGraph,
    chunk_layout: ChunkLayout,
    radius: f32,
}

struct CollisionCheckingResult {
    /// The displacement allowed for the character before hitting a wall. The result of
    /// `math::translate_along(&displacement_vector)` is `displacement_transform`.
    displacement_vector: na::Vector3<f32>,

    /// Multiplying the character's position by this matrix will move the character as far as it can up to its intended
    /// displacement until it hits the wall.
    displacement_transform: na::Matrix4<f32>,

    collision: Option<Collision>,
}

struct Collision {
    /// This collision normal faces away from the collision surface and is given in the perspective of the character
    /// _after_ it is transformed by `allowed_displacement`. The 4th coordinate of this normal vector is assumed to be
    /// 0.0 and is therefore omitted.
    normal: na::UnitVector3<f32>,
}

impl CollisionCheckingResult {
    /// Return a CollisionCheckingResult with no movement and no collision; useful if the character is not moving
    /// and has nothing to check collision against. Also useful as a last resort fallback if an unexpected error occurs.
    fn stationary() -> CollisionCheckingResult {
        CollisionCheckingResult {
            displacement_vector: na::Vector3::zeros(),
            displacement_transform: na::Matrix4::identity(),
            collision: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn apply_normals_internal_examples() {
        // Zero vectors (No-op but should not panic)
        test_apply_normals_internal(&[], [0.60, -0.85, 0.90], 0.2);

        // One vector
        test_apply_normals_internal(&[[-0.48, -0.10, -0.67]], [0.85, -0.53, -0.61], 0.2);

        // Two vectors
        test_apply_normals_internal(
            &[[-0.17, 0.07, -0.38], [-0.85, 0.19, -0.84]],
            [0.19, -0.84, -0.62],
            0.2,
        );

        // Three vectors (Not in use as of the creation of this test but should work anyways)
        test_apply_normals_internal(
            &[
                [-0.24, 0.90, -0.06],
                [-0.91, 0.01, 0.44],
                [0.02, -0.65, -0.12],
            ],
            [0.91, -0.01, -0.61],
            0.2,
        );
    }

    fn test_apply_normals_internal(normals: &[[f32; 3]], subject: [f32; 3], distance: f32) {
        let normals: Vec<na::UnitVector3<f32>> = normals
            .iter()
            .map(|n| na::UnitVector3::new_normalize(na::Vector3::new(n[0], n[1], n[2])))
            .collect();
        let mut subject = na::Vector3::new(subject[0], subject[1], subject[2]);

        apply_normals_internal(&normals, &mut subject, distance);
        for normal in normals {
            assert_abs_diff_eq!(subject.dot(&normal), distance, epsilon = 1.0e-5);
        }
    }
}
