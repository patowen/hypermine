World generation in Hypermine has the following principles:
* The generated contents of each chunk must depend only on the parameters of the world generation algorithm, not on gameplay.
* Everything must be generated in its settled state to ensure that unmodified terrain does not need to be stored on disk, and to help with immersion.

Hypermine generates worlds in three main stages:
* Generate a `NodeState` for each node
* Generate a `ChunkParams` for each chunk, based on the `NodeState` of all 8 nodes adjacent to the chunk's origin
* Using the information in `ChunkParams`, generate the voxel data for each chunk

Every procedural world generation algorithm needs to start with some way of generating noise. This noise is needed for determining properties like temperature and rainfall (which affects what material the ground is made of), and it is also used in generating hills and valleys. While Perlin noise is often used as a standard approach, it is not clear how it would be adapted to a large hyperbolic world. Instead, Hypermine uses a different approach specifically designed for a hyperbolic tiling, inspired by Hyperrogue.

To understand this algorithm, we'll use a 2D analogy, as the 3D version is much the same. We will also focus on temperature as a concrete example of something controlled by noise.

The first step is to determine the temperature associated with the center of each node. The process behind this is described later.

TODO: Picture of pentagonal tiling with a color at the center of each node to represent the temperature.

Then, bilinearly interpolate these temperatures based on voxel coordinates to create a continuous temperature function across the world. Note that in 3D, this would be trilinear interpolation.

TODO: Picture of the same pentagonal tiling with gradients added. Decorate the center of each node with a dot to highlight the control points of the interpolation.

Finally, for each voxel, add a random offset to its temperature, drawn independently from a normal distribution.

TODO: Picture of same pentagonal tiling with the final temperature of each voxel.

## Determining the temperature of each node
To decide on a temperature for each node, we break the dodecahedral tiling up into the planes that divide one group of nodes from another. We associate each plane with a temperature offset, such that crossing a plane increases or decreases the temperature by a specific amount. Once we decide on a temperature of the root node, every other node has its temperature fully determined by the planes crossed by going from the root node to that node.

TODO: Add 2D analogy and suitable pictures