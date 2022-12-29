use common::{
    dodeca::Vertex,
    graph::NodeId,
    node::{Chunk, DualGraph, VoxelData},
    worldgen::ChunkParams,
};
use tokio::{sync::mpsc, runtime::Runtime};

pub struct ChunkLoader {
    send: mpsc::Sender<ChunkDesc>,
    recv: mpsc::Receiver<LoadedChunk>,
    capacity: usize,
    fill: usize,
}

impl ChunkLoader {
    pub fn new(runtime: &mut Runtime, capacity: usize) -> Self {
        let (input_send, mut input_recv) = mpsc::channel::<ChunkDesc>(capacity);
        let (output_send, output_recv) = mpsc::channel::<LoadedChunk>(capacity);
        runtime.spawn(async move {
            while let Some(chunk_desc) = input_recv.recv().await {
                let out = output_send.clone();
                tokio::spawn(async move {
                    let _ = out
                        .send(LoadedChunk {
                            node: chunk_desc.node,
                            chunk: chunk_desc.params.chunk(),
                            voxels: chunk_desc.params.generate_voxels(),
                        })
                        .await;
                });
            }
        });
        ChunkLoader {
            send: input_send,
            recv: output_recv,
            capacity,
            fill: 0,
        }
    }

    /// Begin loading a single chunk, if capacity is available
    pub fn load(&mut self, node: NodeId, params: ChunkParams) -> bool {
        if self.fill == self.capacity {
            return false;
        }
        self.fill += 1;
        if self.send.try_send(ChunkDesc { node, params }).is_err() {
            self.fill -= 1;
            return false;
        }

        true
    }

    /// Move load results into graph data structure, freeing capacity
    pub fn drive(&mut self, graph: &mut DualGraph) {
        while let Ok(chunk) = self.recv.try_recv() {
            self.fill -= 1;
            graph.get_mut(chunk.node).as_mut().unwrap().chunks[chunk.chunk] = Chunk::Populated {
                surface: None,
                voxels: chunk.voxels,
            };
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
