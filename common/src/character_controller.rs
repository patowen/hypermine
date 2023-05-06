use rand_distr::num_traits::Zero;
use tracing::{error, warn};

use crate::{
    character_controller::bound_vector::BoundVector,
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
        let collision_result =
            check_collision(collision_context, position, &remaining_displacement.inner);
        position.local *= collision_result.displacement_transform;

        if let Some(collision) = collision_result.collision {
            // Update the expected displacement to whatever is remaining.
            remaining_displacement.inner -= collision_result.displacement_vector;

            let mut push_direction = &collision.normal;

            if collision.normal.dot(up) > min_slope_up_component {
                if let Some(ground_normal) = ground_normal {
                    if vertical_correction_direction.inner.dot(ground_normal) < 0.0 {
                        vertical_correction_direction.inner.normalize_mut();
                        remaining_displacement.apply(|v| {
                            *v -= vertical_correction_direction.inner * (v.dot(ground_normal))
                                / vertical_correction_direction.inner.dot(ground_normal);
                        });
                        vertical_correction_direction.inner.set_zero();
                    }
                }
                remaining_displacement.apply(|v| {
                    apply_ground_normal_change(up, ground_normal.is_some(), &collision.normal, v);
                });
                *ground_normal = Some(collision.normal);
                push_direction = up;
            }

            if let Some(ground_normal) = ground_normal {
                remaining_displacement.add_temporary_bound(
                    na::UnitVector3::new_unchecked(-ground_normal.as_ref()),
                    *up,
                );
            }
            remaining_displacement.add_and_apply_bound(collision.normal, *push_direction);
            remaining_displacement.remove_temporary_bounds();

            if vertical_correction_direction.inner.dot(&collision.normal) < 0.0 {
                vertical_correction_direction
                    .add_and_apply_bound(collision.normal, collision.normal);
            }
        } else {
            all_collisions_resolved = true;
            break;
        }
    }

    if !all_collisions_resolved {
        warn!("A character entity processed too many collisions and collision resolution was cut short.");
    }

    *velocity = remaining_displacement.tagalong.unwrap();
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
        let collision_result =
            check_collision(collision_context, position, &allowed_displacement.inner);
        if let Some(collision) = collision_result.collision.as_ref() {
            if collision.normal.dot(up) > min_slope_up_component {
                return collision_result;
            }
            allowed_displacement.add_and_apply_bound(collision.normal, collision.normal);
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

mod bound_vector {
    use rand_distr::num_traits::Zero;
    use tracing::warn;

    pub struct BoundVector {
        pub inner: na::Vector3<f32>,
        pub tagalong: Option<na::Vector3<f32>>,
        bounds: VectorBounds,
        inner_distance_factor: f32,
        tagalong_distance_factor: f32,
    }

    impl BoundVector {
        pub fn new(inner: na::Vector3<f32>, tagalong: Option<na::Vector3<f32>>) -> Self {
            // Corrective term to ensure that normals face away from any potential collision surfaces
            const RELATIVE_EPSILON: f32 = 1e-4;

            let inner_distance_factor = inner.magnitude() * RELATIVE_EPSILON;
            let tagalong_distance_factor = tagalong
                .as_ref()
                .map_or(0.0, |t| t.magnitude() * RELATIVE_EPSILON);

            BoundVector {
                inner,
                tagalong,
                bounds: VectorBounds {
                    permanent_bounds: vec![],
                    temporary_bounds: vec![],
                },
                inner_distance_factor,
                tagalong_distance_factor,
            }
        }

        pub fn apply(&mut self, mut f: impl FnMut(&mut na::Vector3<f32>)) {
            f(&mut self.inner);
            if let Some(ref mut tagalong) = self.tagalong {
                f(tagalong);
            }
        }

        pub fn add_and_apply_bound(
            &mut self,
            new_bound_normal: na::UnitVector3<f32>,
            new_bound_push_direction: na::UnitVector3<f32>,
        ) {
            let new_bound = VectorBound {
                normal: new_bound_normal,
                push_direction: new_bound_push_direction,
                target_distance_factor: 1.0,
            };
            self.apply_bound(&new_bound);
            self.bounds.permanent_bounds.push(new_bound);
        }

        fn apply_bound(&mut self, new_bound: &VectorBound) {
            ensure_dot_product(
                self.inner_distance_factor * new_bound.target_distance_factor,
                &new_bound.push_direction,
                &new_bound.normal,
                &mut self.inner,
            );
            if let Some(ref mut tagalong) = self.tagalong {
                ensure_dot_product(
                    self.tagalong_distance_factor * new_bound.target_distance_factor,
                    &new_bound.push_direction,
                    &new_bound.normal,
                    tagalong,
                );
            }

            // Check if all constraints are satisfied
            if self.bounds.iter().all(|b| {
                self.inner.dot(&b.normal)
                    >= self.inner_distance_factor * b.checked_distance_factor()
            }) {
                return;
            }

            // If not all constraints are satisfied, find the first constraint that if applied will satisfy
            // the remaining constriants
            for bound in self.bounds.iter().filter(|b| {
                self.inner.dot(&b.normal) < self.inner_distance_factor * b.checked_distance_factor()
            }) {
                const MIN_ORTHO_NORM: f32 = 1e-5;

                let mut candidate = self.inner;
                let mut ortho_bound_push_direction = bound.push_direction.into_inner();
                ensure_dot_product(
                    0.0,
                    &new_bound.push_direction,
                    &new_bound.normal,
                    &mut ortho_bound_push_direction,
                );

                let Some(ortho_bound_push_direction) =
                    na::UnitVector3::try_new(ortho_bound_push_direction, MIN_ORTHO_NORM)
                else {
                    warn!("Unsatisfied existing bound is parallel to new bound. Is the character squeezed between two walls?");
                    continue;
                };

                ensure_dot_product(
                    self.inner_distance_factor * bound.target_distance_factor,
                    &ortho_bound_push_direction,
                    &bound.normal,
                    &mut candidate,
                );

                if self.bounds.iter().all(|b| {
                    candidate.dot(&b.normal)
                        > self.inner_distance_factor * b.checked_distance_factor()
                }) {
                    self.inner = candidate;
                    if let Some(ref mut tagalong) = self.tagalong {
                        ensure_dot_product(
                            self.tagalong_distance_factor * bound.target_distance_factor,
                            &ortho_bound_push_direction,
                            &bound.normal,
                            tagalong,
                        );
                    }
                    return;
                }
            }

            // If no choice satisfies all constraints, keep all bounds and set the vector to 0
            self.inner.set_zero();
            if let Some(ref mut tagalong) = self.tagalong {
                tagalong.set_zero();
            }
        }

        pub fn add_temporary_bound(
            &mut self,
            normal: na::UnitVector3<f32>,
            push_direction: na::UnitVector3<f32>,
        ) {
            self.bounds.temporary_bounds.push(VectorBound {
                normal,
                push_direction,
                target_distance_factor: -1.0,
            });
        }

        pub fn remove_temporary_bounds(&mut self) {
            self.bounds.temporary_bounds.clear();
        }
    }

    // Updates `moving_vector` by moving it along the line determined by `adjustment_direction` so that
    // its dot product with `fixed_vector` is `dot_product`. This effectively projects vectors onto a plane,
    // where the plane does not need to go through the origin, and the projection can be non-orthogonal.
    // Precondition: For this to be possible, the adjustment direction cannot be orthogonal to the fixed vector.
    fn ensure_dot_product(
        dot_product: f32,
        adjustment_direction: &na::UnitVector3<f32>,
        fixed_vector: &na::UnitVector3<f32>,
        moving_vector: &mut na::Vector3<f32>,
    ) {
        *moving_vector += adjustment_direction.as_ref()
            * ((dot_product - moving_vector.dot(fixed_vector))
                / adjustment_direction.dot(fixed_vector));
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

    struct VectorBound {
        normal: na::UnitVector3<f32>,
        push_direction: na::UnitVector3<f32>,
        target_distance_factor: f32, // Margin of error when the bound is applied
    }

    impl VectorBound {
        // An additional margin of error is needed when the bound is checked to ensure that an
        // applied bound always passes the check.
        fn checked_distance_factor(&self) -> f32 {
            self.target_distance_factor - 1.0
        }
    }

    #[cfg(test)]
    mod tests {
        use approx::assert_abs_diff_eq;

        use super::*;

        #[test]
        fn ensure_dot_product_example() {
            let dot_product = 4.0;
            let adjustment_direction: na::UnitVector3<f32> =
                na::UnitVector3::new_normalize(na::Vector3::new(3.0, -2.0, 7.0));
            let fixed_vector: na::UnitVector3<f32> =
                na::UnitVector3::new_normalize(na::Vector3::new(3.0, -2.0, 7.0));
            let mut moving_vector = na::Vector3::new(-6.0, -3.0, 4.0);
            ensure_dot_product(
                dot_product,
                &adjustment_direction,
                &fixed_vector,
                &mut moving_vector,
            );
            assert_abs_diff_eq!(
                fixed_vector.dot(&moving_vector),
                dot_product,
                epsilon = 1.0e-5
            );
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
