use std::{sync::Arc, time::Instant};

use common::{
    dodeca::{self, Vertex},
    graph::NodeId,
    math::{MIsometry, MPoint},
    node::{Chunk, ChunkId, VoxelData},
};
use metrics::histogram;

use crate::{
    Config, Sim,
    graphics::{Base, Frustum},
    loader::{Cleanup, LoadCtx, LoadFuture, Loadable, Loader, WorkQueue},
};

pub struct WorldgenDriver {
    work_queue: WorkQueue<ChunkDesc>,
}

impl WorldgenDriver {
    pub fn new(config: Arc<Config>, loader: &mut Loader) -> Self {
        Self {
            work_queue: loader.make_queue(config.chunk_load_parallelism as usize),
        }
    }

    pub fn drive(
        &mut self,
        sim: &mut Sim,
        nearby_nodes: &[(NodeId, MIsometry<f32>)],
        frustum: &Frustum,
    ) {
        let drive_worldgen_started = Instant::now();

        // Check for chunks that have finished generating
        while let Some(chunk) = self.work_queue.poll() {
            sim.add_chunk_to_graph(ChunkId::new(chunk.node, chunk.chunk), chunk.voxels);
        }

        // Determine what to load/render
        let view = sim.view();
        if !sim.graph.contains(view.node) {
            // Graph is temporarily out of sync with the server; we don't know where we are, so
            // there's no point trying to draw.
            return;
        }
        let frustum_planes = frustum.planes();
        let local_to_view = view.local.inverse();

        'nearby_nodes: for &(node, ref node_transform) in nearby_nodes {
            let node_to_view = local_to_view * node_transform;
            let origin = node_to_view * MPoint::origin();
            if !frustum_planes.contain(&origin, dodeca::BOUNDING_SPHERE_RADIUS) {
                // Don't bother generating or drawing chunks from nodes that are wholly outside the frustum.
                continue;
            }

            for vertex in Vertex::iter() {
                let chunk_id = ChunkId::new(node, vertex);

                if !matches!(sim.graph[chunk_id], Chunk::Fresh) {
                    continue;
                }

                // Skip chunks beyond the chunk generation distance
                if (node_to_view * vertex.chunk_bounding_sphere_center())
                    .distance(&MPoint::origin())
                    > sim.cfg.chunk_generation_distance + dodeca::CHUNK_BOUNDING_SPHERE_RADIUS
                {
                    continue;
                }

                // Generate voxel data
                let params = common::worldgen::ChunkParams::new(
                    sim.graph.layout().dimension(),
                    &mut sim.graph,
                    chunk_id,
                );
                if let Some(voxel_data) = sim.preloaded_voxel_data.remove(&chunk_id) {
                    sim.add_chunk_to_graph(chunk_id, voxel_data);
                } else if self.work_queue.load(ChunkDesc { node, params }).is_ok() {
                    sim.graph[chunk_id] = Chunk::Generating;
                } else {
                    // No capacity is available in the work queue. Stop trying to prepare chunks to generate.
                    break 'nearby_nodes;
                }
            }
        }
        histogram!("frame.cpu.drive_worldgen").record(drive_worldgen_started.elapsed());
    }
}

struct ChunkDesc {
    node: NodeId,
    params: common::worldgen::ChunkParams,
}

struct LoadedChunk {
    node: NodeId,
    chunk: Vertex,
    voxels: VoxelData,
}

impl Cleanup for LoadedChunk {
    // TODO: `Base` is unused. Try to uncouple Loader from graphics.
    unsafe fn cleanup(self, _gfx: &Base) {}
}

impl Loadable for ChunkDesc {
    type Output = LoadedChunk;
    fn load(self, _ctx: &LoadCtx) -> LoadFuture<'_, Self::Output> {
        Box::pin(async move {
            Ok(LoadedChunk {
                node: self.node,
                chunk: self.params.chunk(),
                voxels: self.params.generate_voxels(),
            })
        })
    }
}
