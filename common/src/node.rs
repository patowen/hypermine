/*the name of this module is pretty arbitrary at the moment*/

use std::ops::{Index, IndexMut};

use serde::{Deserialize, Serialize};

use crate::dodeca::Vertex;
use crate::graph::{Graph, NodeId};
use crate::lru_slab::SlotId;
use crate::proto::Position;
use crate::world::Material;
use crate::worldgen::NodeState;
use crate::{math, Chunks};

pub type DualGraph = Graph<Node>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ChunkId {
    pub node: NodeId,
    pub vertex: Vertex,
}

impl ChunkId {
    pub fn new(node: NodeId, vertex: Vertex) -> Self {
        ChunkId { node, vertex }
    }
}

impl DualGraph {
    pub fn get_chunk_mut(&mut self, chunk: ChunkId) -> Option<&mut Chunk> {
        Some(&mut self.get_mut(chunk.node).as_mut()?.chunks[chunk.vertex])
    }

    pub fn get_chunk(&self, chunk: ChunkId) -> Option<&Chunk> {
        Some(&self.get(chunk.node).as_ref()?.chunks[chunk.vertex])
    }

    pub fn get_chunk_neighbor(
        &self,
        chunk: ChunkId,
        coord_axis: usize,
        coord_direction: i8,
    ) -> Option<ChunkId> {
        if coord_direction == 1 {
            Some(ChunkId::new(
                chunk.node,
                chunk.vertex.adjacent_vertices()[coord_axis],
            ))
        } else {
            Some(ChunkId::new(
                self.neighbor(chunk.node, chunk.vertex.canonical_sides()[coord_axis])?,
                chunk.vertex,
            ))
        }
    }

    pub fn get_block_neighbor(
        &self,
        chunk_size: u8,
        mut chunk: ChunkId,
        mut coords: Coords,
        coord_axis: usize,
        coord_direction: i8,
    ) -> Option<(ChunkId, Coords)> {
        if coords[coord_axis] == chunk_size && coord_direction == 1 {
            let new_vertex = chunk.vertex.adjacent_vertices()[coord_axis];
            let coord_plane0 = (coord_axis + 1) % 3;
            let coord_plane1 = (coord_axis + 2) % 3;
            let mut new_coords = Coords([0; 3]);
            for i in 0..3 {
                if new_vertex.canonical_sides()[i] == chunk.vertex.canonical_sides()[coord_plane0] {
                    new_coords[i] = coords[coord_plane0];
                } else if new_vertex.canonical_sides()[i]
                    == chunk.vertex.canonical_sides()[coord_plane1]
                {
                    new_coords[i] = coords[coord_plane1];
                } else {
                    new_coords[i] = coords[coord_axis];
                }
            }
            coords = new_coords;
            chunk.vertex = new_vertex;
        } else if coords[coord_axis] == 1 && coord_direction == -1 {
            chunk.node = self.neighbor(chunk.node, chunk.vertex.canonical_sides()[coord_axis])?;
        } else {
            coords[coord_axis] = coords[coord_axis].wrapping_add_signed(coord_direction);
        }

        Some((chunk, coords))
    }

    pub fn update_block(
        &mut self,
        chunk_size: u8,
        chunk: ChunkId,
        coords: Coords,
        new_material: Material,
    ) {
        let Some(Chunk::Populated {
            voxels,
            surface,
            old_surface,
        }) = self.get_chunk_mut(chunk)
        else {
            panic!("Tried to update block in nonexistent chunk");
        };
        let voxel = voxels
            .data_mut(chunk_size)
            .get_mut(coords.to_index(chunk_size))
            .expect("coords are in-bounds");

        *voxel = new_material;
        *old_surface = surface.take().or(*old_surface);

        // Remove margins from any adjacent chunks
        for coord_axis in 0..3 {
            if coords[coord_axis] == chunk_size - 1 {
                self.remove_margins(
                    chunk_size,
                    self.get_chunk_neighbor(chunk, coord_axis, 1)
                        .expect("neighboring chunk exists"),
                );
            }
            if coords[coord_axis] == 0 {
                self.remove_margins(
                    chunk_size,
                    self.get_chunk_neighbor(chunk, coord_axis, -1)
                        .expect("neighboring chunk exists"),
                );
            }
        }
    }

    fn remove_margins(&mut self, chunk_size: u8, chunk: ChunkId) {
        let Some(Chunk::Populated {
            voxels,
            surface,
            old_surface,
        }) = self.get_chunk_mut(chunk)
        else {
            panic!("Tried to remove margins of unpopulated chunk");
        };

        if voxels.is_solid() {
            voxels.data_mut(chunk_size);
            *old_surface = surface.take().or(*old_surface);
        }
    }

    /// Returns the up-direction relative to the given position, or `None` if the
    /// position is in an unpopulated node.
    pub fn get_relative_up(&self, position: &Position) -> Option<na::UnitVector3<f32>> {
        self.get(position.node).as_ref().map(|n| {
            na::UnitVector3::new_normalize(
                (math::mtranspose(&position.local) * n.state.up_direction()).xyz(),
            )
        })
    }
}

impl Index<ChunkId> for DualGraph {
    type Output = Chunk;

    fn index(&self, chunk: ChunkId) -> &Chunk {
        self.get_chunk(chunk).unwrap()
    }
}

impl IndexMut<ChunkId> for DualGraph {
    fn index_mut(&mut self, chunk: ChunkId) -> &mut Chunk {
        self.get_chunk_mut(chunk).unwrap()
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Coords(pub [u8; 3]);

impl Coords {
    /// Returns the array index corresponding to these coordinates, including margins
    pub fn to_index(&self, chunk_size: u8) -> usize {
        let chunk_size_with_margin = chunk_size as usize + 2;
        (self.0[0] as usize + 1)
            + (self.0[1] as usize + 1) * chunk_size_with_margin
            + (self.0[2] as usize + 1) * chunk_size_with_margin.pow(2)
    }
}

impl Index<usize> for Coords {
    type Output = u8;

    fn index(&self, coord_axis: usize) -> &u8 {
        self.0.index(coord_axis)
    }
}

impl IndexMut<usize> for Coords {
    fn index_mut(&mut self, coord_axis: usize) -> &mut u8 {
        self.0.index_mut(coord_axis)
    }
}

pub struct Node {
    pub state: NodeState,
    /// We can only populate chunks which lie within a cube of populated nodes, so nodes on the edge
    /// of the graph always have some `None` chunks.
    pub chunks: Chunks<Chunk>,
}

#[derive(Default)]
pub enum Chunk {
    #[default]
    Fresh,
    Generating,
    Populated {
        voxels: VoxelData,
        surface: Option<SlotId>,
        old_surface: Option<SlotId>,
    },
}

/// Like `VoxelData` but designed for sending or receiving over a network where the deserialized version
/// may fail to obey certain constriants. Its purpose is to avoid client-side crashes due to a misbehaving
/// server
#[derive(Serialize, Deserialize)]
pub struct UncheckedVoxelData {
    inner: VoxelData,
}

impl std::fmt::Debug for UncheckedVoxelData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UncheckedVoxelData").finish()
    }
}

impl UncheckedVoxelData {
    pub fn new(inner: VoxelData) -> Self {
        UncheckedVoxelData { inner }
    }

    pub fn validate(self, dimension: u8) -> Option<VoxelData> {
        match self.inner {
            VoxelData::Solid(_) => Some(self.inner),
            VoxelData::Dense(ref d) if d.len() == (usize::from(dimension) + 2).pow(3) => {
                Some(self.inner)
            }
            _ => None,
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub enum VoxelData {
    Solid(Material),
    Dense(Box<[Material]>),
}

impl VoxelData {
    pub fn data_mut(&mut self, dimension: u8) -> &mut [Material] {
        match *self {
            VoxelData::Dense(ref mut d) => d,
            VoxelData::Solid(mat) => {
                let lwm = usize::from(dimension) + 2;
                let mut data = vec![Material::Void; lwm.pow(3)];

                // Populate all blocks except the margins, as margins are not fully implemented yet.
                for i in 1..(lwm - 1) {
                    for j in 1..(lwm - 1) {
                        for k in 1..(lwm - 1) {
                            data[i + j * lwm + k * lwm.pow(2)] = mat;
                        }
                    }
                }
                *self = VoxelData::Dense(data.into_boxed_slice());
                self.data_mut(dimension)
            }
        }
    }

    pub fn get(&self, index: usize) -> Material {
        match *self {
            VoxelData::Dense(ref d) => d[index],
            VoxelData::Solid(mat) => mat,
        }
    }

    pub fn is_solid(&self) -> bool {
        match *self {
            VoxelData::Dense(_) => false,
            VoxelData::Solid(_) => true,
        }
    }
}

/// Contains the context needed to know the locations of individual cubes within a chunk in the chunk's coordinate
/// system. A given `ChunkLayout` is uniquely determined by its dimension.
pub struct ChunkLayout {
    dimension: usize,
    dual_to_grid_factor: f32,
}

impl ChunkLayout {
    pub fn new(dimension: usize) -> Self {
        ChunkLayout {
            dimension,
            dual_to_grid_factor: Vertex::dual_to_chunk_factor() as f32 * dimension as f32,
        }
    }

    /// Number of cubes on one axis of the chunk. Margins are not included.
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Scale by this to convert dual coordinates to homogeneous grid coordinates.
    #[inline]
    pub fn dual_to_grid_factor(&self) -> f32 {
        self.dual_to_grid_factor
    }

    /// Converts a single coordinate from dual coordinates in the Klein-Beltrami model to an integer coordinate
    /// suitable for voxel lookup. Margins are included. Returns `None` if the coordinate is outside the chunk.
    #[inline]
    pub fn dual_to_voxel(&self, dual_coord: f32) -> Option<usize> {
        let floor_grid_coord = (dual_coord * self.dual_to_grid_factor).floor();

        if !(floor_grid_coord >= 0.0 && floor_grid_coord < self.dimension as f32) {
            None
        } else {
            Some(floor_grid_coord as usize + 1)
        }
    }

    /// Converts a single coordinate from grid coordinates to dual coordiantes in the Klein-Beltrami model. This
    /// can be used to find the positions of voxel gridlines.
    #[inline]
    pub fn grid_to_dual(&self, grid_coord: usize) -> f32 {
        grid_coord as f32 / self.dual_to_grid_factor
    }
}

/// Ensures that every new node of the given DualGraph is populated with a [Node] and is
/// ready for world generation.
pub fn populate_fresh_nodes(graph: &mut DualGraph) {
    let fresh = graph.fresh().to_vec();
    graph.clear_fresh();
    for &node in &fresh {
        populate_node(graph, node);
    }
}

fn populate_node(graph: &mut DualGraph, node: NodeId) {
    *graph.get_mut(node) = Some(Node {
        state: graph
            .parent(node)
            .and_then(|i| {
                let parent_state = &graph.get(graph.neighbor(node, i)?).as_ref()?.state;
                Some(parent_state.child(graph, node, i))
            })
            .unwrap_or_else(NodeState::root),
        chunks: Chunks::default(),
    });
}
