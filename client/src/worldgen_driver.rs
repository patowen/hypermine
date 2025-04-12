use std::time::Instant;

use common::{
    dodeca::{self, Vertex},
    graph::{Graph, NodeId},
    math::MPoint,
    node::{Chunk, ChunkId, VoxelData},
    proto::{BlockUpdate, Position},
    traversal,
};
use fxhash::FxHashMap;
use metrics::histogram;

use crate::{
    graphics::Base,
    loader::{Cleanup, LoadCtx, LoadFuture, Loadable, Loader, WorkQueue},
};

pub struct WorldgenDriver {
    work_queue: Option<WorkQueue<ChunkDesc>>, // TODO: This WorkQueue likely suffers from use-after-free due to graphics coupling
    /// Voxel data that have been downloaded from the server for chunks not yet introduced to the graph
    preloaded_block_updates: FxHashMap<ChunkId, Vec<BlockUpdate>>,
    /// Voxel data that has been fetched from the server but not yet introduced to the graph
    preloaded_voxel_data: FxHashMap<ChunkId, VoxelData>,
}

impl WorldgenDriver {
    pub fn new() -> Self {
        Self {
            work_queue: None,
            preloaded_block_updates: FxHashMap::default(),
            preloaded_voxel_data: FxHashMap::default(),
        }
    }

    pub fn init_work_queue(&mut self, loader: &mut Loader, chunk_load_parallelism: usize) {
        self.work_queue = Some(loader.make_queue(chunk_load_parallelism));
    }

    pub fn drive(&mut self, view: Position, chunk_generation_distance: f32, graph: &mut Graph) {
        if self.work_queue.is_none() {
            // If there's no work queue to use, it is impossible to generate chunks.
            return;
        };

        let drive_worldgen_started = Instant::now();

        // Check for chunks that have finished generating
        while let Some(chunk) = self.work_queue.as_mut().unwrap().poll() {
            self.add_chunk_to_graph(graph, ChunkId::new(chunk.node, chunk.chunk), chunk.voxels);
        }

        if !graph.contains(view.node) {
            // Graph is temporarily out of sync with the server; we don't know where we are, so
            // there's no point trying to generate chunks.
            return;
        }
        let local_to_view = view.local.inverse();

        let nearby_nodes = traversal::nearby_nodes(graph, &view, chunk_generation_distance);

        'nearby_nodes: for &(node, ref node_transform) in &nearby_nodes {
            let node_to_view = local_to_view * node_transform;
            let origin = node_to_view * MPoint::origin();

            // Skip nodes beyond the chunk generation distance
            if origin.distance(&MPoint::origin())
                > chunk_generation_distance + dodeca::BOUNDING_SPHERE_RADIUS
            {
                continue;
            }

            for vertex in Vertex::iter() {
                let chunk_id = ChunkId::new(node, vertex);

                if !matches!(graph[chunk_id], Chunk::Fresh) {
                    continue;
                }

                // Skip chunks beyond the chunk generation distance
                if (node_to_view * vertex.chunk_bounding_sphere_center())
                    .distance(&MPoint::origin())
                    > chunk_generation_distance + dodeca::CHUNK_BOUNDING_SPHERE_RADIUS
                {
                    continue;
                }

                // Generate voxel data
                let params =
                    common::worldgen::ChunkParams::new(graph.layout().dimension(), graph, chunk_id);
                if let Some(voxel_data) = self.preloaded_voxel_data.remove(&chunk_id) {
                    self.add_chunk_to_graph(graph, chunk_id, voxel_data);
                } else if (self.work_queue.as_mut().unwrap())
                    .load(ChunkDesc { node, params })
                    .is_ok()
                {
                    graph[chunk_id] = Chunk::Generating;
                } else {
                    // No capacity is available in the work queue. Stop trying to prepare chunks to generate.
                    break 'nearby_nodes;
                }
            }
        }
        histogram!("frame.cpu.drive_worldgen").record(drive_worldgen_started.elapsed());
    }

    /// Adds established voxel data to the graph. This could come from world generation or sent from the server,
    /// depending on whether the chunk has been modified.
    pub fn add_chunk_to_graph(
        &mut self,
        graph: &mut Graph,
        chunk_id: ChunkId,
        voxel_data: VoxelData,
    ) {
        graph.populate_chunk(chunk_id, voxel_data);

        if let Some(block_updates) = self.preloaded_block_updates.remove(&chunk_id) {
            for block_update in block_updates {
                // The chunk was just populated, so a block update should always succeed.
                assert!(graph.update_block(&block_update));
            }
        }
    }

    pub fn apply_block_update(&mut self, graph: &mut Graph, block_update: BlockUpdate) {
        if !graph.update_block(&block_update) {
            self.preloaded_block_updates
                .entry(block_update.chunk_id)
                .or_default()
                .push(block_update);
        }
    }

    pub fn apply_voxel_data(
        &mut self,
        graph: &mut Graph,
        chunk_id: ChunkId,
        voxel_data: VoxelData,
    ) {
        if graph.contains(chunk_id.node) {
            self.add_chunk_to_graph(graph, chunk_id, voxel_data);
        } else {
            self.preloaded_voxel_data.insert(chunk_id, voxel_data);
        }
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
