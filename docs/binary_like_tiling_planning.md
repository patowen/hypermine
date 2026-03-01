# Summary of the tiling
Divide the ground plane into the order-5 square tiling (likely with the same chunk system as before, so these squares start in a subdivided state). Use the Beltrami Klein model to divide the squares evenly (for dividing chunks into voxels).

Then, create a layer of chunks above and below the ground plane. The upper layer of the upper chunk and the lower layer of the lower chunk will be equidistant surfaces (curved). To divide each chunk into voxels, use the Beltrami Klein model for the horizontal coordinates and the actual height above the ground plane for the vertical coordinates.

Note that in general, voxels' top and bottom faces will not be able to have all four vertices be coplanar, so the voxels will effectively be curved. We should pick the right triangulation out of the two possible triangulations to avoid issues with this.

To expand further upward, each chunk will need to split into four chunks to avoid their widths from growing out of control.

# Finding a function between voxel coordinates and hyperboloid coordinates

First, find the formula for the global coordinate system, as if a single chunk was infinitely big. This formula can then be modified with change-of-coordinates formulas to hopefully be adaptable for every chunk.

Let `(x, y, z)` be the voxel coordinates, where `z` is the height. Then, if `z = 0`, we have `v = [x, y, 0, 1] / sqrt(1 - x^2 - y^2)`. Using the sphere analogy, we can derive that modifying the height is as simple as modifying `v.z` and scaling the result to fit on the hyperboloid again. The amount we need to move `v.z` will end up being `tanh(z)`.

`v = [x, y, 0, 1] / sqrt(1 - x^2 - y^2) + [0, 0, tanh(z), 0]` but scaled.

This becomes
```
v = [x, y, sqrt(1 - x^2 - y^2) * tanh(z), 1] / sqrt(1 - x^2 - y^2 - (1 - x^2 - y^2)*tanh(z)^2)
  = [x, y, sqrt(1 - x^2 - y^2) * tanh(z), 1] / sqrt((1 - x^2 - y^2)(1 - tanh(z)^2))
  = [x, y, sqrt(1 - x^2 - y^2) * tanh(z), 1] / sqrt((1 - x^2 - y^2)(sech(z)^2))
  = [x, y, sqrt(1 - x^2 - y^2) * tanh(z), 1] * cosh(z) / sqrt(1 - x^2 - y^2)
```

We can likely ignore this denominator in most circumstances, only normalizing when necessary, which should simplify the formula to `v = [x, y, sqrt(1 - x^2 - y^2) * tanh(z), 1]`.

As a sanity check, we can `mip` this with `[0, 0, 1, 0]` to confirm that the result is what we expect: `sinh(z)`.

```
mip([x, y, sqrt(1 - x^2 - y^2) * tanh(z), 1] * cosh(z) / sqrt(1 - x^2 - y^2), [0, 0, 1, 0])
= sqrt(1 - x^2 - y^2) * tanh(z) * cosh(z) / sqrt(1 - x^2 - y^2)
= tanh(z) * cosh(z)
= sinh(z)
```

As another sanity check, we can `mip` this with the scaled `[x, y, 0, 1]` to confirm that we are going directly upwards without slanting, making the distance to that point equal to the distance to the plane. The resulting `mip` should be `-cosh(z)`.
```
mip([x, y, sqrt(1 - x^2 - y^2) * tanh(z), 1] * cosh(z) / sqrt(1 - x^2 - y^2), [x, y, 0, 1] / sqrt(1 - x^2 - y^2))
= (x * x + y * y - 1) * cosh(z) / (1 - x^2 - y^2)
= -cosh(z)
```

All sanity checks have passed, so the formula is accurate. The next step will be to generalize it.
