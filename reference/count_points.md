# Point density

This function calculates the point density of the input spherical points
at a series of counter stations. Creates `gridsize` regular grid of
counter stations, calculates the distance to all input points at each
counter station, and then calculates the density using `FUN`. Each input
point is weighted by the corresponding item of `weights`. The weights
are normalized to 1 before calculation.

## Usage

``` r
count_points(azi, inc, FUN, sigma, ngrid, weights, r)
```

## Arguments

- azi, inc:

  degrees.

- FUN:

  The method of density estimation to use. Defaults to
  [`exponential_kamb()`](https://tobiste.github.io/structr/reference/density-funs.md).

- sigma:

  (optional) numeric. The number of standard deviations defining the
  expected number of standard deviations by which a random sample from a
  uniform distribution of points would be expected to vary from being
  evenly distributed across the hemisphere. This controls the size of
  the counting circle, and therefore the degree of smoothing. Higher
  `sigma`s will lead to more smoothing of the resulting density
  distribution. This parameter only applies to Kamb-based methods.
  Defaults to `3`.

- ngrid:

  numeric. The size of the grid that the density is estimated on.

- weights:

  (optional) numeric vector of length of `azi`. The relative weight to
  be applied to each input measurement. The array will be normalized to
  sum to 1, so absolute value of the `weights` do not affect the result.
  Defaults to `NULL`

- r:

  numeric. radius of stereonet circle

## Value

list
