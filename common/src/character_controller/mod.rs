mod collision;
mod vector_bounds;

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
    sim_config: &SimConfig,
    graph: &DualGraph,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    on_ground: &mut bool,
    input: &CharacterInput,
    dt_seconds: f32,
) {
    let movement = sanitize_motion_input(input.movement);

    let ctx = CharacterControllerContext {
        cfg: CharacterConfig::new(sim_config),
        collision_context: CollisionContext {
            graph,
            chunk_layout: ChunkLayout::new(sim_config.chunk_size as usize),
            radius: sim_config.character_radius,
        },
        up: graph.get_relative_up(position).unwrap(),
        dt_seconds,
        movement_input: sanitize_motion_input(input.movement),
        jump_input: input.jump,
    };

    if input.no_clip {
        *velocity = movement * ctx.cfg.no_clip_movement_speed;
        *on_ground = false;
        position.local *= math::translate_along(&(*velocity * dt_seconds));
    } else {
        let mut ground_normal = None;
        if *on_ground {
            ground_normal = get_ground_normal(&ctx, position);
        }

        // Handle jumping
        if input.jump && ground_normal.is_some() {
            let horizontal_velocity = *velocity - *ctx.up * ctx.up.dot(velocity);
            *velocity = horizontal_velocity + *ctx.up * ctx.cfg.jump_speed;
            ground_normal = None;
        }

        let old_velocity = *velocity;

        // Update velocity
        if let Some(ground_normal) = ground_normal {
            apply_ground_controls(&ctx, &movement, &ground_normal, velocity);
        } else {
            apply_air_controls(&ctx, &movement, velocity);

            // Apply air resistance
            *velocity *= (-ctx.cfg.air_resistance * dt_seconds).exp();
        }

        // Apply gravity
        *velocity -= *ctx.up * ctx.cfg.gravity_acceleration * dt_seconds;

        // Apply speed cap
        *velocity = velocity.cap_magnitude(ctx.cfg.speed_cap);

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
            &ctx,
            average_velocity,
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
    ctx: &CharacterControllerContext,
    position: &Position,
) -> Option<na::UnitVector3<f32>> {
    // Since the character can be at a corner between a slanted wall and the ground, the first collision
    // directly below the character is not guaranteed to be part of the ground regardless of whether the
    // character is on the ground. To handle this, we repeatedly redirect the direction we search to be
    // parallel to walls we collide with to ensure that we find the ground if is indeed below the character.
    const MAX_COLLISION_ITERATIONS: u32 = 6;
    let mut allowed_displacement = -ctx.up.into_inner() * ctx.cfg.ground_distance_tolerance;
    let mut bounds = VectorBoundGroup::new(&allowed_displacement);

    for _ in 0..MAX_COLLISION_ITERATIONS {
        let collision_result =
            check_collision(&ctx.collision_context, position, &allowed_displacement);
        if let Some(collision) = collision_result.collision.as_ref() {
            if is_ground(ctx, &collision.normal) {
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
fn is_ground(ctx: &CharacterControllerContext, normal: &na::UnitVector3<f32>) -> bool {
    let min_slope_up_component = 1.0 / (ctx.cfg.max_ground_slope.powi(2) + 1.0).sqrt();
    normal.dot(&ctx.up) > min_slope_up_component
}

/// Updates the velocity based on user input assuming the character is on the ground
fn apply_ground_controls(
    ctx: &CharacterControllerContext,
    movement: &na::Vector3<f32>,
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
        math::project_to_plane(&mut unit_movement, ground_normal, &ctx.up, 0.0);
        unit_movement.try_normalize_mut(1e-16);
        unit_movement * movement_norm * ctx.cfg.max_ground_speed
    };

    // Set `ground_velocity` to be the current velocity's ground-parallel component,
    // using a basis that contains the up vector to ensure that the result is unaffected
    // by gravity.
    let mut ground_velocity = *velocity;
    math::project_to_plane(&mut ground_velocity, ground_normal, &ctx.up, 0.0);

    // Adjust the ground-parallel component of the velocity vector to be closer to the
    // target velocity.
    let current_to_target_velocity = target_ground_velocity - ground_velocity;
    let max_delta_velocity = ctx.cfg.ground_acceleration * ctx.dt_seconds;
    if current_to_target_velocity.norm_squared() > max_delta_velocity.powi(2) {
        *velocity += current_to_target_velocity.normalize() * max_delta_velocity;
    } else {
        *velocity += current_to_target_velocity;
    }
}

/// Updates the velocity based on user input assuming the character is in the air
fn apply_air_controls(
    ctx: &CharacterControllerContext,
    movement: &na::Vector3<f32>,
    velocity: &mut na::Vector3<f32>,
) {
    *velocity += movement * ctx.cfg.air_acceleration * ctx.dt_seconds;
}

/// Updates the character's position based on the given average velocity while handling collisions.
/// Also updates the velocity and ground normal based on collisions that occur.
fn apply_velocity(
    ctx: &CharacterControllerContext,
    average_velocity: na::Vector3<f32>,
    position: &mut Position,
    velocity: &mut na::Vector3<f32>,
    ground_normal: &mut Option<na::UnitVector3<f32>>,
) {
    // To prevent an unbounded runtime, we only allow a limited number of collisions to be processed in
    // a single step. If the character encounters excessively complex geometry, it is possible to hit this limit,
    // in which case further movement processing is delayed until the next time step.
    const MAX_COLLISION_ITERATIONS: u32 = 6;

    let mut remaining_dt_seconds = ctx.dt_seconds;

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

        let collision_result =
            check_collision(&ctx.collision_context, position, &expected_displacement);
        position.local *= collision_result.displacement_transform;

        if let Some(collision) = collision_result.collision {
            // Update the expected dt to whatever is remaining.
            remaining_dt_seconds *= 1.0
                - collision_result.displacement_vector.magnitude()
                    / expected_displacement.magnitude();

            handle_collision(
                ctx,
                collision,
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
    ctx: &CharacterControllerContext,
    collision: Collision,
    initial_velocity_info: &VelocityInfo,
    velocity_info: &mut VelocityInfo,
    ground_normal: &mut Option<na::UnitVector3<f32>>,
    ground_collision_handled: &mut bool,
) {
    // Collisions are divided into two categories: Ground collisions and wall collisions.
    // Ground collisions will only affect vertical movement of the character, while wall collisions will
    // push the character away from the wall in a perpendicular direction. If the character is on the ground,
    // we have extra logic to ensure that slanted wall collisions do not lift the character off the ground.
    if is_ground(ctx, &collision.normal) {
        let stay_on_ground_bounds = [VectorBound::new_pull(collision.normal, ctx.up)];
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
                VectorBound::new_push(collision.normal, ctx.up),
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
                VectorBound::new_push(collision.normal, ctx.up),
                &stay_on_ground_bounds,
                &mut velocity_info.average_velocity,
                Some(&mut velocity_info.final_velocity),
            );
        }

        *ground_normal = Some(collision.normal);
    } else {
        let mut stay_on_ground_bounds = Vec::new();
        if let Some(ground_normal) = ground_normal {
            stay_on_ground_bounds.push(VectorBound::new_pull(*ground_normal, ctx.up));
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

/// Contains static configuration information relevant to character physics that generally doesn't change at all
/// for a character, even across timesteps
struct CharacterConfig {
    no_clip_movement_speed: f32,
    max_ground_speed: f32,
    speed_cap: f32,
    max_ground_slope: f32,
    ground_acceleration: f32,
    air_acceleration: f32,
    gravity_acceleration: f32,
    air_resistance: f32,
    jump_speed: f32,
    ground_distance_tolerance: f32,
}

impl CharacterConfig {
    fn new(cfg: &SimConfig) -> Self {
        CharacterConfig {
            no_clip_movement_speed: cfg.no_clip_movement_speed,
            max_ground_speed: cfg.max_ground_speed,
            speed_cap: cfg.speed_cap,
            max_ground_slope: cfg.max_ground_slope,
            ground_acceleration: cfg.ground_acceleration,
            air_acceleration: cfg.air_acceleration,
            gravity_acceleration: cfg.gravity_acceleration,
            air_resistance: cfg.air_resistance,
            jump_speed: cfg.jump_speed,
            ground_distance_tolerance: cfg.ground_distance_tolerance,
        }
    }
}

/// Contains all information about a character that the character controller doesn't change during
/// one of its simulation steps
struct CharacterControllerContext<'a> {
    collision_context: CollisionContext<'a>,
    up: na::UnitVector3<f32>,
    cfg: CharacterConfig,
    dt_seconds: f32,
    movement_input: na::Vector3<f32>,
    jump_input: bool,
}
