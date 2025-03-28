#![allow(clippy::needless_borrowed_reference)]

use rand::{
    Rng,
    distr::{Distribution, StandardUniform},
};

#[macro_use]
mod id;

extern crate nalgebra as na;
pub mod character_controller;
pub mod chunk_collision;
mod chunk_ray_casting;
mod chunks;
pub mod codec;
pub mod collision_math;
pub mod cursor;
pub mod dodeca;
pub mod graph;
pub mod graph_collision;
mod graph_entities;
pub mod graph_ray_casting;
pub mod lru_slab;
mod margins;
pub mod math;
pub mod node;
mod plane;
pub mod proto;
mod sim_config;
pub mod terraingen;
pub mod traversal;
pub mod voxel_math;
pub mod world;
pub mod worldgen;

pub use chunks::Chunks;
pub use graph_entities::GraphEntities;
pub use lru_slab::LruSlab;
pub use sim_config::{SimConfig, SimConfigRaw};

// Stable IDs made of 8 random bytes for easy persistent references
mkid!(EntityId: u64);

impl std::fmt::Display for EntityId {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:016x}", self.0)
    }
}

impl Distribution<EntityId> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> EntityId {
        EntityId(rng.random())
    }
}

pub type Step = i32;

pub fn defer<F: FnOnce()>(f: F) -> Defer<F> {
    Defer::new(f)
}

pub struct Defer<F: FnOnce()>(Option<F>);

impl<F: FnOnce()> Defer<F> {
    pub fn new(f: F) -> Self {
        Self(Some(f))
    }

    pub fn invoke(self) {}

    pub fn cancel(mut self) {
        self.0 = None;
    }
}

impl<F: FnOnce()> Drop for Defer<F> {
    fn drop(&mut self) {
        if let Some(f) = self.0.take() {
            f()
        }
    }
}

/// Clamp speed to to 1.0 and graceful NaN handling
pub fn sanitize_motion_input(v: na::Vector3<f32>) -> na::Vector3<f32> {
    if !v.iter().all(|x| x.is_finite()) {
        return na::Vector3::zeros();
    }
    v / v.norm().max(1.0)
}

pub fn tracing_guard() -> tracing::dispatcher::DefaultGuard {
    use tracing_subscriber::util::SubscriberInitExt;
    tracing_subscriber().set_default()
}

pub fn init_tracing() {
    use tracing_subscriber::util::SubscriberInitExt;
    tracing_subscriber().init();
}

fn tracing_subscriber() -> impl tracing::Subscriber {
    use tracing_subscriber::{filter, fmt, layer::SubscriberExt, registry};

    registry()
        .with(
            fmt::layer()
                .with_target(false)
                .with_ansi(cfg!(not(windows))),
        )
        .with(
            filter::EnvFilter::from_default_env()
                .add_directive(tracing_subscriber::filter::LevelFilter::INFO.into()),
        )
}
