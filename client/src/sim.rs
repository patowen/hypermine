use std::time::Duration;

use fxhash::FxHashMap;
use hecs::Entity;
use tracing::{debug, error, trace};

use crate::{net, prediction::PredictedMotion, Net};
use common::{
    character_controller,
    graph::{Graph, NodeId},
    node::{populate_fresh_nodes, DualGraph},
    proto::{self, Character, CharacterInput, CharacterState, Command, Component, Position},
    sanitize_motion_input, EntityId, GraphEntities, SimConfig, Step,
};

/// Game state
pub struct Sim {
    net: Net,

    // World state
    pub graph: DualGraph,
    pub graph_entities: GraphEntities,
    entity_ids: FxHashMap<EntityId, Entity>,
    pub world: hecs::World,
    pub params: Option<Parameters>,
    pub local_character: Option<Entity>,
    orientation: na::UnitQuaternion<f32>,
    step: Option<Step>,

    // Input state
    since_input_sent: Duration,
    /// Most recent input
    ///
    /// Units are relative to movement speed.
    movement_input: na::Vector3<f32>,
    /// Average input over the current time step. The portion of the timestep which has not yet
    /// elapsed is considered to have zero input.
    ///
    /// Units are relative to movement speed.
    average_movement_input: na::Vector3<f32>,
    no_clip: bool,
    /// Whether no_clip will be toggled next step
    toggle_no_clip: bool,
    jump: bool,
    jump_next_step: bool,
    jump_next_step_sticky: bool,
    prediction: PredictedMotion,
    /// The last extrapolated inter-frame view position, used for rendering and gravity-specific
    /// orientation computations
    view_position: Position,
}

impl Sim {
    pub fn new(net: Net) -> Self {
        Self {
            net,

            graph: Graph::new(),
            graph_entities: GraphEntities::new(),
            entity_ids: FxHashMap::default(),
            world: hecs::World::new(),
            params: None,
            local_character: None,
            orientation: na::one(),
            step: None,

            since_input_sent: Duration::new(0, 0),
            movement_input: na::zero(),
            average_movement_input: na::zero(),
            no_clip: true,
            toggle_no_clip: false,
            jump: false,
            jump_next_step: false,
            jump_next_step_sticky: false,
            prediction: PredictedMotion::new(proto::Position {
                node: NodeId::ROOT,
                local: na::one(),
            }),
            view_position: Position::origin(),
        }
    }

    pub fn rotate(&mut self, delta: &na::UnitQuaternion<f32>) {
        self.orientation *= delta;
    }

    pub fn look(&mut self, delta_yaw: f32, delta_pitch: f32) {
        if self.no_clip {
            self.look_free(delta_yaw, delta_pitch);
        } else {
            self.look_with_gravity(delta_yaw, delta_pitch);
        }
    }

    fn look_free(&mut self, delta_yaw: f32, delta_pitch: f32) {
        self.orientation *= na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), delta_yaw)
            * na::UnitQuaternion::from_axis_angle(&na::Vector3::x_axis(), delta_pitch);
    }

    fn look_with_gravity(&mut self, delta_yaw: f32, delta_pitch: f32) {
        let Some(up) = self
            .get_relative_up(&self.view_position)
            .map(|up| self.orientation.conjugate() * up)
        else {
            return;
        };

        self.orientation *= na::UnitQuaternion::from_axis_angle(&up, delta_yaw);

        if up.x.abs() < 0.9 {
            // Full pitch implementation with logic to prevent turning upside-down
            let current_pitch = -up.z.atan2(up.y);
            let mut target_pitch = current_pitch + delta_pitch;
            if delta_pitch > 0.0 {
                target_pitch = target_pitch
                    .min(std::f32::consts::FRAC_PI_2) // Don't allow pitching up far enough to be upside-down
                    .max(current_pitch); // But if already upside-down, don't make any corrections.
            } else {
                target_pitch = target_pitch
                    .max(-std::f32::consts::FRAC_PI_2) // Don't allow pitching down far enough to be upside-down
                    .min(current_pitch); // But if already upside-down, don't make any corrections.
            }

            self.orientation *= na::UnitQuaternion::from_axis_angle(
                &na::Vector3::x_axis(),
                target_pitch - current_pitch,
            );
        } else {
            // Player is rolled about 90 degrees. Since player view is sideways, we just
            // allow them to pitch as far as they want.
            self.orientation *=
                na::UnitQuaternion::from_axis_angle(&na::Vector3::x_axis(), delta_pitch);
        }
    }

    fn get_horizontal_orientation(&self) -> na::UnitQuaternion<f32> {
        let Some(up) = self
            .get_relative_up(&self.view_position)
            .map(|up| self.orientation.conjugate() * up)
        else {
            return self.orientation;
        };

        if up.x.abs() < 0.9 {
            let current_pitch = -up.z.atan2(up.y);
            self.orientation
                * na::UnitQuaternion::from_axis_angle(&na::Vector3::x_axis(), -current_pitch)
        } else {
            self.orientation
        }
    }

    /// Returns the up-direction relative to the given position
    fn get_relative_up(&self, position: &Position) -> Option<na::UnitVector3<f32>> {
        self.graph.get(position.node).as_ref().map(|node| {
            na::UnitVector3::new_normalize(
                (common::math::mtranspose(&position.local) * node.state.up_direction()).xyz(),
            )
        })
    }

    pub fn set_movement_input(&mut self, mut raw_movement_input: na::Vector3<f32>) {
        if !self.no_clip {
            // Vertical movement keys shouldn't do anything unless no-clip is on.
            raw_movement_input.y = 0.0;
        }
        if raw_movement_input.norm_squared() >= 1.0 {
            // Cap movement input at 1
            raw_movement_input.normalize_mut();
        }
        self.movement_input = raw_movement_input;
    }

    pub fn toggle_no_clip(&mut self) {
        // We prepare to toggle no_clip after the next step instead of immediately, as otherwise,
        // there would be a discontinuity when predicting the player's position within a given step,
        // causing an undesirable jolt.
        self.toggle_no_clip = true;
    }

    pub fn set_jump(&mut self, jump: bool) {
        self.jump_next_step = jump;
        self.jump_next_step_sticky = jump || self.jump_next_step_sticky;
    }

    pub fn params(&self) -> Option<&Parameters> {
        self.params.as_ref()
    }

    pub fn step(&mut self, dt: Duration) {
        self.orientation.renormalize_fast();

        while let Ok(msg) = self.net.incoming.try_recv() {
            self.handle_net(msg);
        }

        if let Some(step_interval) = self.params.as_ref().map(|x| x.cfg.step_interval) {
            self.since_input_sent += dt;
            if let Some(overflow) = self.since_input_sent.checked_sub(step_interval) {
                // At least one step interval has passed since we last sent input, so it's time to
                // send again.

                // Update average movement input for the time between the last input sample and the end of
                // the previous step. dt > overflow because we check whether a step has elapsed
                // after each increment.
                self.average_movement_input += self.movement_input * (dt - overflow).as_secs_f32()
                    / step_interval.as_secs_f32();

                // Send fresh input
                self.send_input();

                // Toggle no clip at the start of a new step
                if self.toggle_no_clip {
                    self.no_clip = !self.no_clip;
                    self.toggle_no_clip = false;
                }

                self.jump = self.jump_next_step || self.jump_next_step_sticky;
                self.jump_next_step_sticky = false;

                // Reset state for the next step
                if overflow > step_interval {
                    // If it's been more than two timesteps since we last sent input, skip ahead
                    // rather than spamming the server.
                    self.average_movement_input = na::zero();
                    self.since_input_sent = Duration::new(0, 0);
                } else {
                    self.average_movement_input =
                        self.movement_input * overflow.as_secs_f32() / step_interval.as_secs_f32();
                    // Send the next input a little sooner if necessary to stay in sync
                    self.since_input_sent = overflow;
                }
            } else {
                // Update average movement input for the time within the current step
                self.average_movement_input +=
                    self.movement_input * dt.as_secs_f32() / step_interval.as_secs_f32();
            }
            self.update_view_position();
            if !self.no_clip {
                self.align_to_gravity();
            }
        }
    }

    fn align_to_gravity(&mut self) {
        let Some(up) = self
            .get_relative_up(&self.view_position)
            .map(|up| self.orientation.conjugate() * up)
        else {
            return;
        };

        if up.z.abs() < 0.9 {
            // If facing not too vertically, roll the camera to make it vertical.
            let delta_roll = -up.x.atan2(up.y);
            self.orientation *=
                na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), delta_roll);
        } else if up.y > 0.0 {
            // Otherwise, if not upside-down, pan the camera to make it vertical.
            let delta_yaw = (up.x / up.z).atan();
            self.orientation *=
                na::UnitQuaternion::from_axis_angle(&na::Vector3::y_axis(), delta_yaw);
        } else {
            // Otherwise, rotate the camera to look straight up or down.
            self.orientation *=
                na::UnitQuaternion::rotation_between(&(na::Vector3::z() * up.z.signum()), &up)
                    .unwrap();
        }
    }

    fn handle_net(&mut self, msg: net::Message) {
        use net::Message::*;
        match msg {
            ConnectionLost(e) => {
                error!("connection lost: {}", e);
            }
            Hello(msg) => {
                self.params = Some(Parameters {
                    character_id: msg.character,
                    cfg: msg.sim_config,
                });
                // Populate the root node
                populate_fresh_nodes(&mut self.graph);
            }
            Spawns(msg) => self.handle_spawns(msg),
            StateDelta(msg) => {
                // Discard out-of-order messages, taking care to account for step counter wrapping.
                if self.step.map_or(false, |x| x.wrapping_sub(msg.step) >= 0) {
                    return;
                }
                self.step = Some(msg.step);
                for &(id, ref new_pos) in &msg.positions {
                    self.update_position(id, new_pos);
                }
                for &(id, ref new_state) in &msg.character_states {
                    self.update_character_state(id, new_state);
                }
                self.reconcile_prediction(msg.latest_input);
            }
        }
    }

    fn update_position(&mut self, id: EntityId, new_pos: &Position) {
        match self.entity_ids.get(&id) {
            None => debug!(%id, "position update for unknown entity"),
            Some(&entity) => match self.world.get::<&mut Position>(entity) {
                Ok(mut pos) => {
                    if pos.node != new_pos.node {
                        self.graph_entities.remove(pos.node, entity);
                        self.graph_entities.insert(new_pos.node, entity);
                    }
                    *pos = *new_pos;
                }
                Err(e) => error!(%id, "position update error: {}", e),
            },
        }
    }

    fn update_character_state(&mut self, id: EntityId, new_character_state: &CharacterState) {
        match self.entity_ids.get(&id) {
            None => debug!(%id, "character state update for unknown entity"),
            Some(&entity) => match self.world.get::<&mut Character>(entity) {
                Ok(mut ch) => {
                    ch.state = new_character_state.clone();
                }
                Err(e) => {
                    error!(%id, "character state update error: {}", e)
                }
            },
        }
    }

    fn reconcile_prediction(&mut self, latest_input: u16) {
        let Some(params) = self.params.as_ref() else {
            return;
        };
        let id = params.character_id;
        let Some(&entity) = self.entity_ids.get(&id) else {
            debug!(%id, "reconciliation attempted for unknown entity");
            return;
        };
        let pos = match self.world.get::<&Position>(entity) {
            Ok(pos) => pos,
            Err(e) => {
                error!(%id, "reconciliation error: {}", e);
                return;
            }
        };
        let ch = match self.world.get::<&Character>(entity) {
            Ok(ch) => ch,
            Err(e) => {
                error!(%id, "reconciliation error: {}", e);
                return;
            }
        };
        self.prediction.reconcile(
            &params.cfg,
            &self.graph,
            latest_input,
            *pos,
            ch.state.velocity,
            ch.state.on_ground,
        );
    }

    fn handle_spawns(&mut self, msg: proto::Spawns) {
        self.step = self.step.max(Some(msg.step));
        let mut builder = hecs::EntityBuilder::new();
        for (id, components) in msg.spawns {
            self.spawn(&mut builder, id, components);
        }
        for &id in &msg.despawns {
            match self.entity_ids.get(&id) {
                Some(&entity) => self.destroy(entity),
                None => error!(%id, "despawned unknown entity"),
            }
        }
        if !msg.nodes.is_empty() {
            trace!(count = msg.nodes.len(), "adding nodes");
        }
        for node in &msg.nodes {
            self.graph.insert_child(node.parent, node.side);
        }
        populate_fresh_nodes(&mut self.graph);
    }

    fn spawn(
        &mut self,
        builder: &mut hecs::EntityBuilder,
        id: EntityId,
        components: Vec<Component>,
    ) {
        trace!(%id, "spawning entity");
        builder.add(id);
        let mut node = None;
        for component in components {
            use common::proto::Component::*;
            match component {
                Character(x) => {
                    builder.add(x);
                }
                Position(x) => {
                    node = Some(x.node);
                    builder.add(x);
                }
            };
        }
        let entity = self.world.spawn(builder.build());
        if let Some(node) = node {
            self.graph_entities.insert(node, entity);
        }
        if id == self.params.as_ref().unwrap().character_id {
            self.local_character = Some(entity);
        }
        if let Some(x) = self.entity_ids.insert(id, entity) {
            self.destroy_idless(x);
            error!(%id, "id collision");
        }
    }

    fn send_input(&mut self) {
        let params = self.params.as_ref().unwrap();
        let orientation = if self.no_clip {
            self.orientation
        } else {
            self.get_horizontal_orientation()
        };
        let character_input = CharacterInput {
            movement: sanitize_motion_input(orientation * self.average_movement_input),
            jump: self.jump,
            no_clip: self.no_clip,
        };
        let generation = self
            .prediction
            .push(&params.cfg, &self.graph, &character_input);

        // Any failure here will be better handled in handle_net's ConnectionLost case
        let _ = self.net.outgoing.send(Command {
            generation,
            character_input,
            orientation: self.orientation,
        });
    }

    fn update_view_position(&mut self) {
        let mut view_position = *self.prediction.predicted_position();
        let mut view_velocity = *self.prediction.predicted_velocity();
        let mut view_on_ground = *self.prediction.predicted_on_ground();
        let orientation = if self.no_clip {
            self.orientation
        } else {
            self.get_horizontal_orientation()
        };
        if let Some(ref params) = self.params {
            // Apply input that hasn't been sent yet
            let predicted_input = CharacterInput {
                // We divide by how far we are through the timestep because self.average_movement_input
                // is always over the entire timestep, filling in zeroes for the future, and we
                // want to use the average over what we have so far. Dividing by zero is handled
                // by the character_controller sanitizing this input.
                movement: orientation * self.average_movement_input
                    / (self.since_input_sent.as_secs_f32()
                        / params.cfg.step_interval.as_secs_f32()),
                jump: self.jump,
                no_clip: self.no_clip,
            };
            character_controller::run_character_step(
                &params.cfg,
                &self.graph,
                &mut view_position,
                &mut view_velocity,
                &mut view_on_ground,
                &predicted_input,
                self.since_input_sent.as_secs_f32(),
            );
        }

        // Rotate the player orientation to stay consistent with changes in gravity
        if !self.no_clip {
            if let (Some(old_up), Some(new_up)) = (
                self.get_relative_up(&self.view_position),
                self.get_relative_up(&view_position),
            ) {
                self.orientation = na::UnitQuaternion::rotation_between_axis(&old_up, &new_up)
                    .unwrap_or(na::UnitQuaternion::identity())
                    * self.orientation;
            }
        }

        self.view_position = view_position;
    }

    pub fn view(&self) -> Position {
        Position {
            node: self.view_position.node,
            local: self.view_position.local * self.orientation.to_homogeneous(),
        }
    }

    /// Destroy all aspects of an entity
    fn destroy(&mut self, entity: Entity) {
        let id = *self
            .world
            .get::<&EntityId>(entity)
            .expect("destroyed nonexistent entity");
        self.entity_ids.remove(&id);
        self.destroy_idless(entity);
    }

    /// Destroy an entity without an EntityId mapped
    fn destroy_idless(&mut self, entity: Entity) {
        if let Ok(position) = self.world.get::<&Position>(entity) {
            self.graph_entities.remove(position.node, entity);
        }
        self.world
            .despawn(entity)
            .expect("destroyed nonexistent entity");
    }
}

/// Simulation details received on connect
pub struct Parameters {
    pub cfg: SimConfig,
    pub character_id: EntityId,
}
