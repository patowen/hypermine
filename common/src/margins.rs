use crate::{
    dodeca::Vertex,
    graph::Graph,
    math,
    node::{Chunk, ChunkId, VoxelData},
    voxel_math::{ChunkDirection, CoordAxis, CoordSign, Coords, SimpleChunkOrientation},
    world::Material,
};

pub fn fix_margins(
    dimension: u8,
    destination_vertex: Vertex,
    destination: &mut VoxelData,
    direction: ChunkDirection,
    source: &mut VoxelData,
) {
    let neighbor_orientation = match direction.sign {
        CoordSign::Plus => {
            destination_vertex.adjacent_chunk_orientations()[direction.axis as usize]
        }
        CoordSign::Minus => SimpleChunkOrientation::identity(),
    };

    let margin_coord = match direction.sign {
        CoordSign::Plus => dimension + 1,
        CoordSign::Minus => 0,
    };
    let edge_coord = match direction.sign {
        CoordSign::Plus => dimension,
        CoordSign::Minus => 1,
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

/// Updates the margins of a given VoxelData to match the voxels they're next to. This is a good assumption to start
/// with before taking into account neighboring chunks because it results in the least rendering and is generally accurate when
/// the neighboring chunks are solid.
pub fn initialize_margins(dimension: u8, voxels: &mut VoxelData) {
    // If voxels is solid, the margins are already set up the way they should be.
    if voxels.is_solid() {
        return;
    }

    for direction in ChunkDirection::iter() {
        let margin_coord = match direction.sign {
            CoordSign::Plus => dimension + 1,
            CoordSign::Minus => 0,
        };
        let edge_coord = match direction.sign {
            CoordSign::Plus => dimension,
            CoordSign::Minus => 1,
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
}

pub fn update_margin_voxel(
    graph: &mut Graph,
    chunk: ChunkId,
    coords: Coords,
    direction: ChunkDirection,
    material: Material,
) {
    let dimension = graph.layout().dimension();
    let edge_coord = match direction.sign {
        CoordSign::Plus => dimension - 1,
        CoordSign::Minus => 0,
    };
    if coords[direction.axis] != edge_coord {
        // There is nothing to do if we're not on an edge voxel.
        return;
    }
    let Some(Chunk::Populated {
        modified: _neighbor_modified,
        voxels: neighbor_voxels,
        surface: neighbor_surface,
        old_surface: neighbor_old_surface,
    }) = graph
        .get_chunk_neighbor(chunk, direction.axis, direction.sign)
        .map(|chunk_id| &mut graph[chunk_id])
    else {
        // If the neighboring chunk to check is not populated, there is nothing to do.
        return;
    };

    let margin_coord = match direction.sign {
        CoordSign::Plus => dimension + 1,
        CoordSign::Minus => 0,
    };
    let neighbor_orientation = match direction.sign {
        CoordSign::Plus => chunk.vertex.adjacent_chunk_orientations()[direction.axis as usize],
        CoordSign::Minus => SimpleChunkOrientation::identity(),
    };
    let mut neighbor_coords = CoordsWithMargins::from(coords);
    neighbor_coords[direction.axis] = margin_coord;
    neighbor_coords = neighbor_orientation * neighbor_coords;

    neighbor_voxels.data_mut(dimension)[neighbor_coords.to_index(dimension)] = material;
    *neighbor_old_surface = neighbor_surface.take().or(*neighbor_old_surface);
}

/// Coordinates for a discrete voxel within a chunk, including margins
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoordsWithMargins(pub [u8; 3]);

impl CoordsWithMargins {
    /// Returns the array index in `VoxelData` corresponding to these coordinates
    pub fn to_index(self, chunk_size: u8) -> usize {
        let chunk_size_with_margin = chunk_size as usize + 2;
        (self.0[0] as usize)
            + (self.0[1] as usize) * chunk_size_with_margin
            + (self.0[2] as usize) * chunk_size_with_margin.pow(2)
    }
}

impl From<Coords> for CoordsWithMargins {
    fn from(value: Coords) -> Self {
        CoordsWithMargins([value.0[0] + 1, value.0[1] + 1, value.0[2] + 1])
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
        dodeca::Vertex,
        voxel_math::{CoordAxis, CoordSign, Coords},
        world::Material,
    };

    use super::*;

    #[test]
    fn test_fix_margins2() {
        let mut destination = VoxelData::Solid(Material::Void);
        destination.data_mut(12)[Coords([11, 2, 10]).to_index(12)] = Material::WoodPlanks;

        let mut source = VoxelData::Solid(Material::Void);
        source.data_mut(12)[Coords([2, 10, 11]).to_index(12)] = Material::WoodPlanks;

        println!("{:?}", Vertex::F.adjacent_vertices()[0]);
        println!("{:?}", Vertex::J.adjacent_vertices()[2]);

        fix_margins(
            12,
            Vertex::F,
            &mut destination,
            ChunkDirection {
                axis: CoordAxis::X,
                sign: CoordSign::Plus,
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
