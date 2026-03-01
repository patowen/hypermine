# Summary of the tiling
Divide the ground plane into the order-5 square tiling (likely with the same chunk system as before, so these squares start in a subdivided state). Use the Beltrami Klein model to divide the squares evenly (for dividing chunks into voxels).

Then, create a layer of chunks above and below the ground plane. The upper layer of the upper chunk and the lower layer of the lower chunk will be equidistant surfaces (curved). To divide each chunk into voxels, use the Beltrami Klein model for the horizontal coordinates and the actual height above the ground plane for the vertical coordinates.

Note that in general, voxels' top and bottom faces will not be able to have all four vertices be coplanar, so the voxels will effectively be curved. We should pick the right triangulation out of the two possible triangulations to avoid issues with this.

To expand further upward, each chunk will need to split into four chunks to avoid their widths from growing out of control.

# Finding a function between voxel coordinates and hyperboloid coordinates

First, find the formula for the global coordinate system, as if a single chunk was infinitely big. This formula can then be modified with change-of-coordinates formulas to hopefully be adaptable for every chunk.

Let `(X, Y, Z)` be the voxel coordinates, where `Z` is the height. Then, if `Z = 0`, we have `v = [X, Y, 0, 1] / sqrt(1 - X^2 - Y^2)`. Using the sphere analogy, we can derive that modifying the height is as simple as modifying `z` (lowercase to note that this is hyperboloid coordinates and not voxel coordinates) and scaling the result to fit on the hyperboloid again. The amount we need to move `z` will end up being `tanh(Z)`.

`v = [X, Y, 0, 1] / sqrt(1 - X^2 - Y^2) + [0, 0, tanh(Z), 0]` but scaled.

This becomes
```
v = [X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] / sqrt(1 - X^2 - Y^2 - (1 - X^2 - Y^2)*tanh(Z)^2)
  = [X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] / sqrt((1 - X^2 - Y^2)(1 - tanh(Z)^2))
  = [X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] / sqrt((1 - X^2 - Y^2)(sech(Z)^2))
  = [X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] * cosh(Z) / sqrt(1 - X^2 - Y^2)
```

We can likely ignore this denominator in most circumstances, only normalizing when necessary, which should simplify the formula to `v = [X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1]`.

As a sanity check, we can `mip` this with `[0, 0, 1, 0]` to confirm that the result is what we expect: `sinh(Z)`.

```
mip([X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] * cosh(Z) / sqrt(1 - X^2 - Y^2), [0, 0, 1, 0])
= sqrt(1 - X^2 - Y^2) * tanh(Z) * cosh(Z) / sqrt(1 - X^2 - Y^2)
= tanh(Z) * cosh(Z)
= sinh(Z)
```

As another sanity check, we can `mip` this with the scaled `[X, Y, 0, 1]` to confirm that we are going directly upwards without slanting, making the distance to that point equal to the distance to the plane. The resulting `mip` should be `-cosh(Z)`.
```
mip([X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] * cosh(Z) / sqrt(1 - X^2 - Y^2), [X, Y, 0, 1] / sqrt(1 - X^2 - Y^2))
= (X * X + Y * Y - 1) * cosh(Z) / (1 - X^2 - Y^2)
= -cosh(Z)
```

All sanity checks have passed, so the formula is accurate. The next step will be to generalize it. We can use the simpler un-normalized formula `[X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1]` for now.

# Increasing height
To perform a Lorentz-boost in the z-direction, we leave x and y fixed, while updating `z -> cz + sw` and `w -> sz + cw` where `c = cosh(boost)` and `s = sinh(boost)`.
This results in `[X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] -> [X, Y, c * sqrt(1 - X^2 - Y^2) * tanh(Z) + s, c + s * sqrt(1 - X^2 - Y^2) * tanh(Z)]`

There isn't an obvious way to simplify that formula, so we will need to be systematic.

Observation: `X` and `Y` will never get that far from the origin due to the tiling, so `sqrt(1 - X^2 - Y^2)` will remain numerically stable. We can keep that separate in our calculations.

Note: High up chunks will have `X` and `Y` values that barely change, so we will need a different representation to avoid precision issues. However, I believe this can be solved entirely separately.

As for height, one approach might be to make the w-coordinate equal to 1 when X and Y are 0, and Z is equal to the negative of the boost. In such a situation, `tanh(Z) = -s/c`.
```
w = c + s * sqrt(1 - X^2 - Y^2) * tanh(Z)
  = c + s * sqrt(1 - 0^2 - 0^2) * -s/c
  = c - s^2 / c
  = c - (c^2 - 1) / c
  = c^2/c - (c^2 - 1) / c
  = 1 / c
```

In other words, normalizing `w` to 1 in this case can be done by dividing by `1 / cosh(boost)`, or, in other words, multiplying by `cosh(boost)`. This resembles the simplified normalized formula `[X, Y, sqrt(1 - X^2 - Y^2) * tanh(Z), 1] * cosh(Z) / sqrt(1 - X^2 - Y^2)`, since that formula also has a `cosh(Z)` factor in it, so it checks out.

This leaves us with
`[X * c, Y * c, c^2 * sqrt(1 - X^2 - Y^2) * tanh(Z) + c * s, c^2 + c * s * sqrt(1 - X^2 - Y^2) * tanh(Z)]`

We're not done. c and s can all get unreasonably large, and `tanh(Z)` can get incredibly close to 1, so catastrophic cancellation is guaranteed without more simplifying. The factor of `sqrt(1 - X^2 - Y^2)` is also a potential issue, but that is likely solvable later on by ensuring that we transform things to allow the factor to stay very close to 1. Because of that, we will ignore that factor for now and assume it's equal to 1 to see what kinds of simplifications we would be able to do.

If `sqrt(1 - X^2 - Y^2) = 1`, and `Z = -b + U` then
```
z = c^2 * sqrt(1 - X^2 - Y^2) * tanh(Z) + c * s
  = c^2 * tanh(Z) + c * s
  = cosh(b)^2 * tanh(-b + U) + cosh(b) * sinh(b)
  = cosh(b)^2 * sinh(-b + U) / cosh(-b + U) + cosh(b) * sinh(b)
  = cosh(b)^2 * (sinh(-b)cosh(U) + cosh(-b)sinh(U)) / (cosh(-b)cosh(U) + sinh(-b)sinh(U)) + cosh(b) * sinh(b)
  = cosh(b)^2 * (-sinh(b)cosh(U) + cosh(b)sinh(U)) / (cosh(b)cosh(U) - sinh(b)sinh(U)) + cosh(b) * sinh(b)

Try multiplying this by (cosh(b)cosh(U) - sinh(b)sinh(U)) to see if we can cancel anything in the numerator

cosh(b)^2 * (-sinh(b)cosh(U) + cosh(b)sinh(U)) + cosh(b)sinh(b) * (cosh(b)cosh(U) - sinh(b)sinh(U))
= -sinh(b)cosh(b)^2cosh(U) + cosh(b)^3sinh(U) + cosh(b)^2sinh(b)cosh(U) - cosh(b)sinh(b)^2sinh(U)
= cosh(U)[ -sinh(b)cosh(b)^2 + cosh(b)^2sinh(b) ] + sinh(U)[ cosh(b)^3 - cosh(b)sinh(b)^2 ]
= sinh(U)[ cosh(b)^3 - cosh(b)sinh(b) ]
= sinh(U)cosh(b)[ cosh(b)^2 - sinh(b)^2 ]
= sinh(U)cosh(b)

Divide by the denominator again

z = sinh(U)cosh(b) / (cosh(b)cosh(U) - sinh(b)sinh(U))
  = sinh(U) / (cosh(U) - tanh(b)sinh(U))
```

This final formula is numerically stable, approaching `sinh(U) / (cosh(U) - sinh(U))` as b approaches infinity. Since `U` remains small, this is fine.
