use std::collections::VecDeque;

use crate::{
    cursor::{ChunkDirection, CoordAxis, CoordDirection, Coords, SimpleChunkOrientation},
    dodeca::Vertex,
    graph::Graph,
    math,
    node::{Chunk, ChunkId, VoxelData},
};

pub fn fix_margins(graph: &mut Graph, chunk: ChunkId, direction: ChunkDirection) {
    let dimension = graph.layout().dimension();
    let (neighbor_chunk, neighbor_orientation) = match direction.direction {
        CoordDirection::Plus => (
            ChunkId::new(
                chunk.node,
                chunk.vertex.adjacent_vertices()[direction.axis as usize],
            ),
            chunk.vertex.adjacent_chunk_orientations()[direction.axis as usize],
        ),
        CoordDirection::Minus => (
            ChunkId::new(
                graph
                    .neighbor(
                        chunk.node,
                        chunk.vertex.canonical_sides()[direction.axis as usize],
                    )
                    .unwrap(),
                chunk.vertex,
            ),
            SimpleChunkOrientation::identity(),
        ),
    };

    let Some(Chunk::Populated {
        voxels: neighbor_chunk_data,
        ..
    }) = graph.get_chunk(neighbor_chunk)
    else {
        // TODO: Decide best way to fix margins when neighboring chunk isn't generated yet
        return;
    };

    let neighbor_direction = neighbor_orientation * direction; // TODO: Double-check that change of coordinates is in correct direction
                                                               // TODO: This is wrong. we don't want the margins of the neighbor. We want the edges of the neighbor.
    let neighbor_margin_coord = match neighbor_direction.direction {
        CoordDirection::Plus => dimension + 1,
        CoordDirection::Minus => 0,
    };
    let mut margin_contents = VecDeque::with_capacity((dimension as usize).pow(2));
    for j in 0..dimension {
        for i in 0..dimension {
            let neighbor_coords = CoordsWithMargins(math::tuv_to_xyz(
                neighbor_direction.axis as usize,
                [i + 1, j + 1, neighbor_margin_coord],
            ));
            margin_contents.push_back(neighbor_chunk_data.get(neighbor_coords.to_index(dimension)));
        }
    }

    let Some(Chunk::Populated {
        voxels: chunk_data, ..
    }) = graph.get_chunk_mut(chunk)
    else {
        panic!();
    };
    let chunk_data = chunk_data.data_mut(dimension);
    let margin_coord = match direction.direction {
        CoordDirection::Plus => dimension + 1,
        CoordDirection::Minus => 0,
    };
    for j in 0..dimension {
        for i in 0..dimension {
            let coords = CoordsWithMargins(math::tuv_to_xyz(
                direction.axis as usize,
                [i + 1, j + 1, margin_coord],
            ));
            chunk_data[coords.to_index(dimension)] = margin_contents.pop_front().unwrap();
        }
    }
}

pub fn fix_margins2(
    dimension: u8,
    destination_vertex: Vertex,
    destination: &mut VoxelData,
    direction: ChunkDirection,
    source: &mut VoxelData,
) {
    let neighbor_orientation = match direction.direction {
        CoordDirection::Plus => {
            destination_vertex.adjacent_chunk_orientations()[direction.axis as usize]
        }
        CoordDirection::Minus => SimpleChunkOrientation::identity(),
    };

    let margin_coord = match direction.direction {
        CoordDirection::Plus => dimension + 1,
        CoordDirection::Minus => 0,
    };
    let edge_coord = match direction.direction {
        CoordDirection::Plus => dimension,
        CoordDirection::Minus => 1,
    };
    let chunk_data = destination.data_mut(dimension);
    let neighbor_chunk_data = source.data_mut(dimension);
    for j in 0..dimension {
        for i in 0..dimension {
            chunk_data[CoordsWithMargins(math::tuv_to_xyz(
                direction.axis as usize,
                [margin_coord, i + 1, j + 1],
            ))
            .to_index(dimension)] = neighbor_chunk_data[(neighbor_orientation
                * CoordsWithMargins(math::tuv_to_xyz(
                    direction.axis as usize,
                    [edge_coord, i + 1, j + 1],
                )))
            .to_index(dimension)];

            neighbor_chunk_data[(neighbor_orientation
                * CoordsWithMargins(math::tuv_to_xyz(
                    direction.axis as usize,
                    [margin_coord, i + 1, j + 1],
                )))
            .to_index(dimension)] = chunk_data[CoordsWithMargins(math::tuv_to_xyz(
                direction.axis as usize,
                [edge_coord, i + 1, j + 1],
            ))
            .to_index(dimension)];
        }
    }
}

pub fn initialize_margins(dimension: u8, voxels: &mut VoxelData, direction: ChunkDirection) {
    // If voxels is solid, the margins are already set up the way they should be.
    if voxels.is_solid() {
        return;
    }

    let margin_coord = match direction.direction {
        CoordDirection::Plus => dimension + 1,
        CoordDirection::Minus => 0,
    };
    let edge_coord = match direction.direction {
        CoordDirection::Plus => dimension,
        CoordDirection::Minus => 1,
    };
    let chunk_data = voxels.data_mut(dimension);
    for j in 0..dimension {
        for i in 0..dimension {
            chunk_data[CoordsWithMargins(math::tuv_to_xyz(
                direction.axis as usize,
                [margin_coord, i + 1, j + 1],
            ))
            .to_index(dimension)] = chunk_data[CoordsWithMargins(math::tuv_to_xyz(
                direction.axis as usize,
                [edge_coord, i + 1, j + 1],
            ))
            .to_index(dimension)];
        }
    }
}

/// Coordinates for a discrete voxel within a chunk, including margins
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoordsWithMargins(pub [u8; 3]);

impl CoordsWithMargins {
    /// Returns the array index in `VoxelData` corresponding to these coordinates
    pub fn to_index(&self, chunk_size: u8) -> usize {
        let chunk_size_with_margin = chunk_size as usize + 2;
        (self.0[0] as usize)
            + (self.0[1] as usize) * chunk_size_with_margin
            + (self.0[2] as usize) * chunk_size_with_margin.pow(2)
    }
}

impl std::ops::Index<CoordAxis> for CoordsWithMargins {
    type Output = u8;

    fn index(&self, coord_axis: CoordAxis) -> &u8 {
        self.0.index(coord_axis as usize)
    }
}

impl std::ops::IndexMut<CoordAxis> for CoordsWithMargins {
    fn index_mut(&mut self, coord_axis: CoordAxis) -> &mut u8 {
        self.0.index_mut(coord_axis as usize)
    }
}

impl std::ops::Mul<CoordsWithMargins> for SimpleChunkOrientation {
    type Output = CoordsWithMargins;

    fn mul(self, rhs: CoordsWithMargins) -> Self::Output {
        let mut result = CoordsWithMargins([0; 3]);
        for axis in CoordAxis::iter() {
            result[self[axis]] = rhs[axis];
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        cursor::{CoordAxis, CoordDirection, Coords},
        dodeca::Vertex,
        graph::{Graph, NodeId},
        node::ChunkId,
        world::Material,
    };

    use super::*;

    #[derive(Clone)]
    struct Cell {
        chunk: ChunkId,
        coords: Coords,
        axes: [CoordAxis; 3],
        directions: [CoordDirection; 3],
    }

    fn cell_neighbor(
        graph: &Graph,
        cell: &Cell,
        cell_coord_axis: CoordAxis,
        cell_coord_direction: CoordDirection,
    ) -> Option<Cell> {
        let coord_axis = cell.axes[cell_coord_axis as usize];
        let coord_direction: CoordDirection =
            cell.directions[cell_coord_axis as usize] * cell_coord_direction;
        if cell.coords[coord_axis] == graph.layout().dimension() - 1
            && coord_direction == CoordDirection::Plus
        {
            let new_vertex = cell.chunk.vertex.adjacent_vertices()[coord_axis as usize];
            // Permute coordinates based on differences in the canonical orders between the old
            // and new vertex
            let [coord_plane0, coord_plane1] = coord_axis.other_axes();
            let mut new_coords = Coords([0; 3]);
            for current_axis in CoordAxis::iter() {
                if new_vertex.canonical_sides()[current_axis as usize]
                    == cell.chunk.vertex.canonical_sides()[coord_plane0 as usize]
                {
                    new_coords[current_axis] = cell.coords[coord_plane0];
                } else if new_vertex.canonical_sides()[current_axis as usize]
                    == cell.chunk.vertex.canonical_sides()[coord_plane1 as usize]
                {
                    new_coords[current_axis] = cell.coords[coord_plane1];
                } else {
                    new_coords[current_axis] = cell.coords[coord_axis];
                }
            }
            Some(Cell {
                chunk: ChunkId::new(cell.chunk.node, new_vertex),
                coords: new_coords,
                axes: cell.axes, // TODO: Correct axes and directions
                directions: cell.directions,
            })
        } else if cell.coords[coord_axis] == 0 && coord_direction == CoordDirection::Minus {
            let mut new_cell = cell.clone();
            new_cell.chunk.node = graph.neighbor(
                cell.chunk.node,
                cell.chunk.vertex.canonical_sides()[coord_axis as usize],
            )?;
            new_cell.directions[cell_coord_axis as usize] *= CoordDirection::Minus;
            Some(new_cell)
        } else {
            let mut new_cell = cell.clone();
            new_cell.coords[coord_axis] =
                new_cell.coords[coord_axis].wrapping_add_signed(coord_direction as i8);
            Some(new_cell)
        }
    }

    fn bar() {
        let graph = Graph::new(12);
        let chunk = graph.get_chunk(ChunkId::new(NodeId::ROOT, Vertex::A));

        //graph.get_block_neighbor(chunk, coords, coord_axis, coord_direction);
    }

    #[test]
    fn foo() {}

    #[test]
    fn test_fix_margins2() {
        let mut destination = VoxelData::Solid(Material::Void);
        destination.data_mut(12)[Coords([11, 2, 10]).to_index(12)] = Material::WoodPlanks;

        let mut source = VoxelData::Solid(Material::Void);
        source.data_mut(12)[Coords([2, 10, 11]).to_index(12)] = Material::WoodPlanks;

        println!("{:?}", Vertex::F.adjacent_vertices()[0]);
        println!("{:?}", Vertex::J.adjacent_vertices()[2]);

        fix_margins2(
            12,
            Vertex::F,
            &mut destination,
            ChunkDirection {
                axis: CoordAxis::X,
                direction: CoordDirection::Plus,
            },
            &mut source,
        );

        println!(
            "{:?}",
            destination.get(CoordsWithMargins([12, 3, 11]).to_index(12))
        );
        println!(
            "{:?}",
            source.get(CoordsWithMargins([3, 11, 12]).to_index(12))
        );
        println!(
            "{:?}",
            destination.get(CoordsWithMargins([13, 3, 11]).to_index(12))
        );
        println!(
            "{:?}",
            source.get(CoordsWithMargins([3, 11, 13]).to_index(12))
        );
    }
}
