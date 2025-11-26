World generation in Hypermine has the following principles:
* The generated contents of each chunk must depend only on the parameters of the world generation algorithm, not on gameplay.
* Everything must be generated in its settled state to ensure that unmodified terrain does not need to be stored on disk, and to help with immersion.

Hypermine generates worlds in three main stages:
* Generate a `NodeState` for each node
* Generate a `ChunkParams` for each chunk, based on the `NodeState` of all 8 nodes adjacent to the chunk's origin
* Using the information in `ChunkParams`, generate the voxel data for each chunk

Every procedural world generation algorithm needs to start with some way of generating noise. While Perlin noise is often used as a standard approach, it is not clear how it would be adapted to a large hyperbolic world. Instead, Hypermine uses a different approach specifically designed for a hyperbolic tiling, inspired by Hyperrogue.

To understand this algorithm, we use a 2D analogy.
TOOD: Describe pentagonal tessellation, interpreting it as a bunch of lines, and making it so that crossing each line adds a particular number to the value.