use crate::worldgen::ChunkParams;
use crate::{
    dodeca::Vertex,
    graph::NodeId,
    math,
    node::{Chunk, DualGraph, VoxelData},
    proto::Position,
    traversal::nearby_nodes,
};
use tokio::{runtime::Handle, sync::mpsc};

pub struct ChunkLoader {
    send: mpsc::Sender<ChunkDesc>,
    recv: mpsc::Receiver<LoadedChunk>,
    capacity: usize,
    fill: usize,
}

impl ChunkLoader {
    pub fn new(runtime: Handle, capacity: usize) -> Self {
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

    pub fn load_chunks(
        &mut self,
        graph: &mut DualGraph,
        dimension: u8,
        position: &Position,
        distance: f64,
    ) {
        let mut nodes = nearby_nodes(graph, position, distance);
        // Sort nodes by distance to the view to prioritize loading closer data and improve early Z
        // performance
        let view_pos = position.local * math::origin();
        nodes.sort_unstable_by(|&(_, ref xf_a), &(_, ref xf_b)| {
            math::mip(&view_pos, &(xf_a * math::origin()))
                .partial_cmp(&math::mip(&view_pos, &(xf_b * math::origin())))
                .unwrap_or(std::cmp::Ordering::Less)
        });

        for &(node, _) in &nodes {
            for chunk in Vertex::iter() {
                if let Chunk::Fresh = graph
                    .get(node)
                    .as_ref()
                    .expect("all nodes must be populated before rendering")
                    .chunks[chunk]
                {
                    if let Some(params) = ChunkParams::new(dimension, graph, node, chunk) {
                        if self.load(node, params) {
                            graph.get_mut(node).as_mut().unwrap().chunks[chunk] = Chunk::Generating;
                        }
                    }
                }
            }
        }
    }

    /// Begin loading a single chunk, if capacity is available
    fn load(&mut self, node: NodeId, params: ChunkParams) -> bool {
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
    params: ChunkParams,
}

struct LoadedChunk {
    node: NodeId,
    chunk: Vertex,
    voxels: VoxelData,
}
