#[cfg(test)]
mod tests {
    use crate::{
        dodeca::Vertex,
        graph::{Graph, NodeId},
        node::{ChunkId, CoordAxis, CoordDirection, Coords},
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
}
