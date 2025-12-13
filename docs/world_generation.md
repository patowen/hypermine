World generation in Hypermine has the following principles:
* The generated contents of each chunk must depend only on the parameters of the world generation algorithm, not on gameplay.
* Everything must be generated in its settled state to ensure that unmodified terrain does not need to be stored on disk, and to help with immersion.
* World generation should be as isotropic as reasonably possible, so that it is not obvious which direction leads to the origin.

Hypermine generates worlds in three main stages:
* Generate a `NodeState` for each node
* Generate a `ChunkParams` for each chunk, based on the `NodeState` of all 8 nodes adjacent to the chunk's origin
* Using the information in `ChunkParams`, generate the voxel data for each chunk

Every procedural world generation algorithm needs to start with some way of generating noise. This noise is needed for determining properties like temperature and rainfall (which affects what material the ground is made of), and it is also used in generating hills and valleys. While Perlin noise is often used as a standard approach, it is not clear how it would be adapted to a large hyperbolic world. Instead, Hypermine uses a different approach specifically designed for a hyperbolic tiling, inspired by Hyperrogue.

To understand this algorithm, we'll use a 2D analogy, as the 3D version is much the same.

## Coarse noise function output
The first step is to form a coarse approximation of the noise function, deciding on one value for each node. To do this, we take advantage of the tiling itself.

For a 2D analogy, the pentagonal tiling of the hyperbolic plane can be thought of as a set of lines dividing the hyperbolic plane instead of individual pentagons.

TODO: Picture of pentagonal tiling with lines colored to distinguish them from each other

Similarly, in 3D, the dodecahedral tiling can be thought of as a set of planes dividing hyperbolic space. This interpretation of the dodecahedral tiling is important for understanding how the noise function works between nodes.

To decide on a noise value for each node, we break the dodecahedral tiling up into this planes. We associate each plane with a randomly chosen noise value offset, such that crossing a specific plane in one direction increases or decreases the noise value by a specific amount, and crossing the same plane from the other side has the opposite effect. Once we decide on a noise value for the root node, this definition fully determines the noise value of every other node.

The following diagram shows an example the 2D equivalent of this algorithm.

TODO: Picture of pentagonal tiling with each line labeled with the noise value offset, using arrows or something similar to show how this offset applies. The center of each pentagon is also labeled with a number with its current noise value. Integers are used everywhere to allow the reader to verify the math easily in their head.

In this diagram, the randomly-chosen noise value offset of each line, along with the derived noise value of each node, is shown. Note how the difference in noise values between any two adjacent nodes matches the noise value offset of the line dividing them.

The algorithm allows for random variation while keeping nearby nodes similar to each other, which proves useful for multiple aspects of world generation.

## Fine noise function output
The next step is to use this coarse approximation of the noise function to produce the actual noise function. To begin, set the noise value at the center of each node to the coarse output we computed earlier.

TODO: Picture of pentagonal tiling with a color at the center of each node to represent the noise value.

Then, bilinearly interpolate these values based on voxel coordinates to create a continuous function across the world. Note that in 3D, this would be trilinear interpolation.

TODO: Picture of the same pentagonal tiling with gradients added. Decorate the center of each node with a dot to highlight the control points of the interpolation.

Finally, for each voxel, add a random offset to its noise value, drawn independently from some distribution.

TODO: Picture of same pentagonal tiling with the final noise value of each voxel.

## Terrain shape
Hypermine uses the 3D noise function to determine the shape of the terrain. This may seem surprising, as it is arguably simpler to use a 2D noise function instead to form a heightmap of the world. However, using 3D noise instead is a useful way of generating more interesting terrain with overhangs, and more importantly, a 2D heightmap works less well in hyperbolic space because of the way that space expands as you move away from the ground plane. The naive approach would cause hills and valleys to have significantly less detail than terrain near the ground plane.

Instead, the basic algorithm is as follows: Using a 3D noise function, we determine the hypothetical elevation of the terrain at each voxel. We then subtract this elevation from the actual height of the voxel above the ground plane to determine a value that can be roughly translated to voxel's height relative to the terrain's surface. If this value is below zero, we are inside the terrain, and the voxel should be solid, and otherwise, it should be void. If the value is below zero but close to zero, one of the surface materials (like dirt) will be used, while if the value is far below zero, a material like stone will be used instead.

Note that the above is a simplification of the actual algorithm. It is recommended to read the implementation of `worldgen::ChunkParams::generate_terrain` to understand all the details. For instance, terracing is used to add flatter terrain layers with steeper hills between them, and the strength of this terracing effect is controlled by another parameter affected by noise called `blockiness`. In addition, some measures are taken to make the terrain surface smoother than the dirt/stone interface.

## Terrain material
TODO
