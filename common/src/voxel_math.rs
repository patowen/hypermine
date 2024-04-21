use std::ops::{Index, IndexMut};

use serde::{Deserialize, Serialize};

use crate::dodeca::Side;

/// Represents a particular axis in a voxel grid
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordAxis {
    X = 0,
    Y = 1,
    Z = 2,
}

/// Trying to convert a `usize` to a `CoordAxis` returns this struct if the provided
/// `usize` is out-of-bounds
#[derive(Debug, Clone, Copy)]
pub struct CoordAxisOutOfBounds;

impl CoordAxis {
    /// Iterates through the the axes in ascending order
    pub fn iter() -> impl ExactSizeIterator<Item = Self> {
        [Self::X, Self::Y, Self::Z].iter().copied()
    }

    /// Returns the pair axes orthogonal to the current axis
    pub fn other_axes(self) -> [Self; 2] {
        match self {
            Self::X => [Self::Y, Self::Z],
            Self::Y => [Self::Z, Self::X],
            Self::Z => [Self::X, Self::Y],
        }
    }
}

impl TryFrom<usize> for CoordAxis {
    type Error = CoordAxisOutOfBounds;

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::X),
            1 => Ok(Self::Y),
            2 => Ok(Self::Z),
            _ => Err(CoordAxisOutOfBounds),
        }
    }
}

/// Represents a direction in a particular axis. This struct is meant to be used with a coordinate axis,
/// so when paired with the X-axis, it represents the postitive X-direction when set to Plus and the
/// negative X-direction when set to Minus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordSign {
    Plus = 1,
    Minus = -1,
}

impl CoordSign {
    /// Iterates through the two possible coordinate directions
    pub fn iter() -> impl ExactSizeIterator<Item = Self> {
        [CoordSign::Plus, CoordSign::Minus].iter().copied()
    }
}

impl std::ops::Mul for CoordSign {
    type Output = CoordSign;

    fn mul(self, rhs: Self) -> Self::Output {
        if self == rhs {
            CoordSign::Plus
        } else {
            CoordSign::Minus
        }
    }
}

impl std::ops::MulAssign for CoordSign {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

/// Coordinates for a discrete voxel within a chunk, not including margins
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Coords(pub [u8; 3]);

impl Coords {
    /// Returns the array index in `VoxelData` corresponding to these coordinates
    pub fn to_index(self, chunk_size: u8) -> usize {
        let chunk_size_with_margin = chunk_size as usize + 2;
        (self.0[0] as usize + 1)
            + (self.0[1] as usize + 1) * chunk_size_with_margin
            + (self.0[2] as usize + 1) * chunk_size_with_margin.pow(2)
    }
}

impl Index<CoordAxis> for Coords {
    type Output = u8;

    fn index(&self, coord_axis: CoordAxis) -> &u8 {
        self.0.index(coord_axis as usize)
    }
}

impl IndexMut<CoordAxis> for Coords {
    fn index_mut(&mut self, coord_axis: CoordAxis) -> &mut u8 {
        self.0.index_mut(coord_axis as usize)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChunkDirection {
    pub axis: CoordAxis,
    pub sign: CoordSign,
}

impl ChunkDirection {
    pub fn iter() -> impl ExactSizeIterator<Item = ChunkDirection> {
        [
            ChunkDirection {
                axis: CoordAxis::X,
                sign: CoordSign::Plus,
            },
            ChunkDirection {
                axis: CoordAxis::Y,
                sign: CoordSign::Plus,
            },
            ChunkDirection {
                axis: CoordAxis::Z,
                sign: CoordSign::Plus,
            },
            ChunkDirection {
                axis: CoordAxis::X,
                sign: CoordSign::Minus,
            },
            ChunkDirection {
                axis: CoordAxis::Y,
                sign: CoordSign::Minus,
            },
            ChunkDirection {
                axis: CoordAxis::Z,
                sign: CoordSign::Minus,
            },
        ]
        .iter()
        .copied()
    }
}

/// Represents one of the 48 possible orientations a chunk can be viewed from, including reflections.
/// This is analogous to a 3x3 rotation matrix with a restricted domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChunkOrientation {
    directions: [ChunkDirection; 3],
}

impl Index<CoordAxis> for ChunkOrientation {
    type Output = ChunkDirection;

    fn index(&self, index: CoordAxis) -> &Self::Output {
        &self.directions[index as usize]
    }
}

impl std::ops::Mul<ChunkDirection> for ChunkOrientation {
    type Output = ChunkDirection;

    fn mul(self, rhs: ChunkDirection) -> Self::Output {
        ChunkDirection {
            axis: self[rhs.axis].axis,
            sign: self[rhs.axis].sign * rhs.sign,
        }
    }
}

impl std::ops::Mul<ChunkOrientation> for ChunkOrientation {
    type Output = ChunkOrientation;

    fn mul(self, rhs: ChunkOrientation) -> Self::Output {
        ChunkOrientation {
            directions: rhs.directions.map(|d| self * d),
        }
    }
}

impl std::ops::MulAssign for ChunkOrientation {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

/// Represents one of the 6 possible orientations a chunk can be viewed from, including reflections, that preserve the origin.
/// This is analogous to a 3x3 rotation matrix with a restricted domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SimpleChunkOrientation {
    axes: [CoordAxis; 3],
}

impl SimpleChunkOrientation {
    pub fn identity() -> Self {
        SimpleChunkOrientation {
            axes: [CoordAxis::X, CoordAxis::Y, CoordAxis::Z],
        }
    }

    pub fn from_permutation(from: [Side; 3], to: [Side; 3]) -> Self {
        assert!(from[0] != from[1] && from[0] != from[2] && from[1] != from[2]);
        assert!(to[0] != to[1] && to[0] != to[2] && to[1] != to[2]);
        SimpleChunkOrientation {
            axes: from.map(|f| {
                CoordAxis::try_from(
                    to.iter()
                        .position(|&t| f == t)
                        .expect("from and to must have same set of sides"),
                )
                .unwrap()
            }),
        }
    }
}

impl Index<CoordAxis> for SimpleChunkOrientation {
    type Output = CoordAxis;

    fn index(&self, index: CoordAxis) -> &Self::Output {
        &self.axes[index as usize]
    }
}

impl std::ops::Mul<Coords> for SimpleChunkOrientation {
    type Output = Coords;

    fn mul(self, rhs: Coords) -> Self::Output {
        let mut result = Coords([0; 3]);
        for axis in CoordAxis::iter() {
            result[self[axis]] = rhs[axis];
        }
        result
    }
}

impl std::ops::Mul<ChunkDirection> for SimpleChunkOrientation {
    type Output = ChunkDirection;

    fn mul(self, rhs: ChunkDirection) -> Self::Output {
        ChunkDirection {
            axis: self[rhs.axis],
            sign: rhs.sign,
        }
    }
}
