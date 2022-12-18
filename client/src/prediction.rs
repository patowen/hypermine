use std::collections::VecDeque;

use common::proto::Position;

/// Predicts the result of motion inputs in-flight to the server
///
/// When sending input to the server, call `push` to record the input in a local queue of in-flight
/// inputs, and to obtaining a generation tag to send alongside the input. The server echos the
/// highest tag it's received alongside every state update, which we then use in `reconcile` to
/// determine which inputs have been integrated into the server's state and no longer need to be
/// predicted.
pub struct PredictedMotion<T: Clone> {
    log: VecDeque<T>,
    generation: u16,
    predicted_position: Position,
    predicted_velocity: na::Vector3<f32>,
}

impl<T: Clone> PredictedMotion<T> {
    pub fn new(initial_position: Position) -> Self {
        Self {
            log: VecDeque::new(),
            generation: 0,
            predicted_position: initial_position,
            predicted_velocity: na::Vector3::zeros(),
        }
    }

    /// Update for input about to be sent to the server, returning the generation it should be
    /// tagged with
    pub fn push(
        &mut self,
        input: &T,
        mut step_function: impl FnMut(&mut Position, &mut na::Vector3<f32>, &T),
    ) -> u16 {
        step_function(
            &mut self.predicted_position,
            &mut self.predicted_velocity,
            input,
        );
        self.log.push_back(input.clone());
        self.generation = self.generation.wrapping_add(1);
        self.generation
    }

    /// Update with the latest state received from the server and the generation it was based on
    pub fn reconcile(
        &mut self,
        generation: u16,
        position: Position,
        velocity: na::Vector3<f32>,
        mut step_function: impl FnMut(&mut Position, &mut na::Vector3<f32>, &T),
    ) {
        let first_gen = self.generation.wrapping_sub(self.log.len() as u16);
        let obsolete = usize::from(generation.wrapping_sub(first_gen));
        if obsolete > self.log.len() || obsolete == 0 {
            // We've already processed a state incorporating equal or more recent input
            return;
        }
        self.log.drain(..obsolete);
        self.predicted_position = position;
        self.predicted_velocity = velocity;

        for input in self.log.iter() {
            step_function(
                &mut self.predicted_position,
                &mut self.predicted_velocity,
                input,
            );
        }
    }

    /// Latest estimate of the server's state after receiving all `push`ed inputs.
    pub fn predicted_position(&self) -> &Position {
        &self.predicted_position
    }

    pub fn predicted_velocity(&self) -> &na::Vector3<f32> {
        &self.predicted_velocity
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

    /// An arbitrary velocity
    fn vel() -> na::Vector3<f32> {
        na::Vector3::zeros()
    }

    #[test]
    fn wraparound() {
        let mut pred: PredictedMotion<nalgebra::Vector3<f32>> = PredictedMotion::new(pos());
        pred.generation = u16::max_value() - 1;
        assert_eq!(pred.push(&na::Vector3::x(), |_, _, _| {}), u16::max_value());
        assert_eq!(pred.push(&na::Vector3::x(), |_, _, _| {}), 0);
        assert_eq!(pred.log.len(), 2);

        pred.reconcile(u16::max_value() - 1, pos(), vel(), |_, _, _| {});
        assert_eq!(pred.log.len(), 2);
        pred.reconcile(u16::max_value(), pos(), vel(), |_, _, _| {});
        assert_eq!(pred.log.len(), 1);
        pred.reconcile(0, pos(), vel(), |_, _, _| {});
        assert_eq!(pred.log.len(), 0);
    }
}
