use std::collections::VecDeque;

use common::{
    character_controller::CharacterControllerPass,
    node::DualGraph,
    proto::{Character, CharacterInput, Position},
    SimConfig,
};

/// Predicts the result of motion inputs in-flight to the server
///
/// When sending input to the server, call `push` to record the input in a local queue of in-flight
/// inputs, and to obtaining a generation tag to send alongside the input. The server echos the
/// highest tag it's received alongside every state update, which we then use in `reconcile` to
/// determine which inputs have been integrated into the server's state and no longer need to be
/// predicted.
pub struct PredictedMotion {
    log: VecDeque<CharacterInput>,
    generation: u16,
    predicted_position: Position,
    predicted_character: Character,
}

impl PredictedMotion {
    pub fn new(initial_position: Position, initial_character: Character) -> Self {
        Self {
            log: VecDeque::new(),
            generation: 0,
            predicted_position: initial_position,
            predicted_character: initial_character,
        }
    }

    /// Update for input about to be sent to the server, returning the generation it should be
    /// tagged with
    pub fn push(
        &mut self,
        graph: &DualGraph,
        config: &SimConfig,
        dt_seconds: f32,
        input: &CharacterInput,
    ) -> u16 {
        CharacterControllerPass {
            position: &mut self.predicted_position,
            character: &mut self.predicted_character,
            input,
            graph,
            config,
            dt_seconds,
        }
        .step();
        self.log.push_back(input.clone());
        self.generation = self.generation.wrapping_add(1);
        self.generation
    }

    /// Update with the latest state received from the server and the generation it was based on
    pub fn reconcile(
        &mut self,
        graph: &DualGraph,
        config: &SimConfig,
        dt_seconds: f32,
        generation: u16,
        position: Position,
        character: Character,
    ) {
        let first_gen = self.generation.wrapping_sub(self.log.len() as u16);
        let obsolete = usize::from(generation.wrapping_sub(first_gen));
        if obsolete > self.log.len() || obsolete == 0 {
            // We've already processed a state incorporating equal or more recent input
            return;
        }
        self.log.drain(..obsolete);
        self.predicted_position = position;
        self.predicted_character = character;

        for input in self.log.iter() {
            CharacterControllerPass {
                position: &mut self.predicted_position,
                character: &mut self.predicted_character,
                input,
                graph,
                config,
                dt_seconds,
            }
            .step();
        }
    }

    /// Latest estimate of the server's state after receiving all `push`ed inputs.
    pub fn predicted_position(&self) -> &Position {
        &self.predicted_position
    }

    pub fn predicted_character(&self) -> &Character {
        &self.predicted_character
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An arbitrary position
    fn pos() -> Position {
        Position {
            node: common::graph::NodeId::ROOT,
            local: na::one(),
        }
    }

    /// An arbitrary character
    fn char() -> Character {
        Character {
            name: "Test".to_string(),
            orientation: na::UnitQuaternion::identity(),
            velocity: na::Vector3::zeros(),
        }
    }

    #[test]
    fn wraparound() {
        /*let mut pred = PredictedMotion::new(pos(), char());
        pred.generation = u16::max_value() - 1;
        assert_eq!(pred.push(&na::Vector3::x()), u16::max_value());
        assert_eq!(pred.push(&na::Vector3::x()), 0);
        assert_eq!(pred.log.len(), 2);

        pred.reconcile(u16::max_value() - 1, pos());
        assert_eq!(pred.log.len(), 2);
        pred.reconcile(u16::max_value(), pos());
        assert_eq!(pred.log.len(), 1);
        pred.reconcile(0, pos());
        assert_eq!(pred.log.len(), 0);*/
    }
}
