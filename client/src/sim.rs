use std::time::Duration;

use fxhash::FxHashMap;
use hecs::Entity;
use tracing::{debug, error, trace};

use crate::{
    local_character_controller::LocalCharacterController, net, prediction::PredictedMotion, Net,
};
use common::{
    character_controller,
    graph::{Graph, NodeId},
    graph_collision::Ray,
    graph_ray_casting,
    node::{populate_fresh_nodes, Chunk, ChunkId, ChunkLayout, DualGraph},
    proto::{
        self, BlockUpdate, Character, CharacterInput, CharacterState, Command, Component, Position,
    },
    sanitize_motion_input,
    world::Material,
    EntityId, GraphEntities, SimConfig, Step,
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
    /// Whether the current step starts with a jump
    is_jumping: bool,
    /// Whether the jump button has been pressed since the last step
    jump_pressed: bool,
    /// Whether the jump button is currently held down
    jump_held: bool,
    /// Whether the place-block button has been pressed since the last step
    place_block_pressed: bool,
    /// Whether the break-block button has been pressed since the last step
    break_block_pressed: bool,
    prediction: PredictedMotion,
    local_character_controller: LocalCharacterController,
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
            step: None,

            since_input_sent: Duration::new(0, 0),
            movement_input: na::zero(),
            average_movement_input: na::zero(),
            no_clip: true,
            toggle_no_clip: false,
            is_jumping: false,
            jump_pressed: false,
            jump_held: false,
            place_block_pressed: false,
            break_block_pressed: false,
            prediction: PredictedMotion::new(proto::Position {
                node: NodeId::ROOT,
                local: na::one(),
            }),
            local_character_controller: LocalCharacterController::new(),
        }
    }

    /// Rotates the camera's view in a context-dependent manner based on the desired yaw and pitch angles.
    pub fn look(&mut self, delta_yaw: f32, delta_pitch: f32, delta_roll: f32) {
        if self.no_clip {
            self.local_character_controller
                .look_free(delta_yaw, delta_pitch, delta_roll);
        } else {
            self.local_character_controller
                .look_level(delta_yaw, delta_pitch);
        }
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

    pub fn set_jump_held(&mut self, jump_held: bool) {
        self.jump_held = jump_held;
        self.jump_pressed = jump_held || self.jump_pressed;
    }

    pub fn set_jump_pressed_true(&mut self) {
        self.jump_pressed = true;
    }

    pub fn set_place_block_pressed_true(&mut self) {
        self.place_block_pressed = true;
    }

    pub fn set_break_block_pressed_true(&mut self) {
        self.break_block_pressed = true;
    }

    pub fn params(&self) -> Option<&Parameters> {
        self.params.as_ref()
    }

    pub fn step(&mut self, dt: Duration) {
        self.local_character_controller.renormalize_orientation();

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

                self.is_jumping = self.jump_held || self.jump_pressed;
                self.jump_pressed = false;

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
                self.local_character_controller.align_to_gravity();
            }
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
                for block_update in msg.block_updates.into_iter() {
                    let Some(node_id) = self.graph.from_hash(block_update.node_hash) else {
                        tracing::warn!("Block update received from unknown node hash");
                        continue;
                    };
                    let Some(Chunk::Populated { voxels, .. }) = self
                        .graph
                        .get_chunk_mut(ChunkId::new(node_id, block_update.vertex))
                    else {
                        tracing::warn!("Block update received from ungenerated chunk");
                        continue;
                    };
                    let Some(voxel) = voxels
                        .data_mut(self.params.as_ref().unwrap().cfg.chunk_size)
                        .get_mut(block_update.coords as usize)
                    else {
                        tracing::warn!("Block update received for out-of-bounds block");
                        continue;
                    };
                    *voxel = block_update.new_material;
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
            self.local_character_controller.orientation()
        } else {
            self.local_character_controller.horizontal_orientation()
        };
        let character_input = CharacterInput {
            movement: sanitize_motion_input(orientation * self.average_movement_input),
            jump: self.is_jumping,
            no_clip: self.no_clip,
            block_update: self.get_block_update(),
        };
        let generation = self
            .prediction
            .push(&params.cfg, &self.graph, &character_input);

        // Any failure here will be better handled in handle_net's ConnectionLost case
        let _ = self.net.outgoing.send(Command {
            generation,
            character_input,
            orientation: self.local_character_controller.orientation(),
        });
    }

    fn update_view_position(&mut self) {
        let mut view_position = *self.prediction.predicted_position();
        let mut view_velocity = *self.prediction.predicted_velocity();
        let mut view_on_ground = *self.prediction.predicted_on_ground();
        let orientation = if self.no_clip {
            self.local_character_controller.orientation()
        } else {
            self.local_character_controller.horizontal_orientation()
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
                jump: self.is_jumping,
                no_clip: self.no_clip,
                block_update: None,
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

        self.local_character_controller.update_position(
            view_position,
            self.graph.get_relative_up(&view_position).unwrap(),
            !self.no_clip,
        )
    }

    pub fn view(&self) -> Position {
        self.local_character_controller.oriented_position()
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

    fn get_block_update(&self) -> Option<BlockUpdate> {
        let placing = if self.place_block_pressed {
            true
        } else if self.break_block_pressed {
            false
        } else {
            return None;
        };
        let dimension = self.params.as_ref().unwrap().cfg.chunk_size;

        let view_position = self.local_character_controller.oriented_position();
        let ray_casing_result = graph_ray_casting::ray_cast(
            &self.graph,
            &ChunkLayout::new(dimension as usize),
            &view_position,
            &Ray::new(na::Vector4::w(), -na::Vector4::z()),
            0.5,
        );

        let Ok(ray_casting_result) = ray_casing_result else {
            tracing::warn!("Tried to run a raycast beyond generated terrain.");
            return None;
        };

        let hit = ray_casting_result?;

        let block_pos = if placing {
            self.get_block_neighbor(
                hit.chunk,
                hit.voxel_coords,
                hit.face_axis as usize,
                hit.face_direction as isize,
            )?
        } else {
            (hit.chunk, hit.voxel_coords)
        };

        let lwm = dimension as u32 + 2;
        let voxel_coords = (block_pos.1[0] as u32 + 1)
            + (block_pos.1[1] as u32 + 1) * lwm
            + (block_pos.1[2] as u32 + 1) * lwm * lwm;

        let material = if placing {
            Material::Wood
        } else {
            Material::Void
        };

        Some(BlockUpdate {
            node_hash: self.graph.hash_of(block_pos.0.node),
            vertex: block_pos.0.vertex,
            coords: voxel_coords,
            new_material: material,
        })
    }

    fn get_block_neighbor(
        &self,
        mut chunk: ChunkId,
        mut coords: [usize; 3],
        coord_axis: usize,
        coord_direction: isize,
    ) -> Option<(ChunkId, [usize; 3])> {
        let dimension = self.params.as_ref().unwrap().cfg.chunk_size as usize;
        if coords[coord_axis] == dimension - 1 && coord_direction == 1 {
            let new_vertex = chunk.vertex.adjacent_vertices()[coord_axis];
            let coord_plane0 = (coord_axis + 1) % 3;
            let coord_plane1 = (coord_axis + 2) % 3;
            let mut new_coords: [usize; 3] = [0; 3];
            for (i, new_coord) in new_coords.iter_mut().enumerate() {
                if new_vertex.canonical_sides()[i] == chunk.vertex.canonical_sides()[coord_plane0] {
                    *new_coord = coords[coord_plane0];
                } else if new_vertex.canonical_sides()[i]
                    == chunk.vertex.canonical_sides()[coord_plane1]
                {
                    *new_coord = coords[coord_plane1];
                } else {
                    *new_coord = coords[coord_axis];
                }
            }
            coords = new_coords;
            chunk.vertex = new_vertex;
        } else if coords[coord_axis] == 0 && coord_direction == -1 {
            chunk.node = self
                .graph
                .neighbor(chunk.node, chunk.vertex.canonical_sides()[coord_axis])?;
        } else {
            coords[coord_axis] = (coords[coord_axis] as isize + coord_direction) as usize;
        }

        Some((chunk, coords))
    }
}

/// Simulation details received on connect
pub struct Parameters {
    pub cfg: SimConfig,
    pub character_id: EntityId,
}
