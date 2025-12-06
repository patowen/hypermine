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
The pentagonal tiling of the hyperbolic plane can be thought of as a set of lines dividing the hyperbolic plane instead of individual pentagons.

TODO: Picture of pentagonal tiling with lines colored to distinguish them from each other

Similarly, the dodecahedral tiling can be thought of as a set of planes dividing hyperbolic space. This interpretation of the dodecahedral tiling is important for understanding how the noise function works between nodes.

To decide on a temperature for each node, we break the dodecahedral tiling up into this planes. We associate each plane with a randomly chosen temperature offset, such that crossing a specific plane in one direction increases or decreases the temperature by a specific amount, and crossing the same plane from the other side has the opposite effect. Once we decide on a temperature of the root node, this definition fully determines the temperature of every other node.

The following diagram shows an example the 2D equivalent of this algorithm.

TODO: Picture of pentagonal tiling with each line labeled with the temperature offset, using arrows or something similar to show how this offset applies. The center of each pentagon is also labeled with a number with its current temperature. Integers are used everywhere to allow the reader to verify the math easily in their head.

In this diagram, the randomly-chosen temperature offset of each line, along with the derived temperature of each node, is shown. Note how the difference in temperature between any two adjacent nodes matches the temperature offset of the line dividing them.

The resulting noise generation algorithm allows for random variation while keeping nearby nodes similar to each other, which proves useful for multiple aspects of world generation.

## Hills and Valleys
TODO
