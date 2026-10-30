# Statistics of a Bingham Distribution

Statistics of a Bingham Distribution

## Usage

``` r
bingham_statistics(x, w = NULL, na.rm = TRUE)
```

## Source

Borradaile, G. (2003). Spherical-Orientation Data. In: Statistics of
Earth Science Data. Springer, Berlin, Heidelberg.
https://doi.org/10.1007/978-3-662-05223-5_10

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- w:

  numeric. Optional weights for each observation.

- na.rm:

  logical. Whether `NA` values should be removed before the computation
  proceeds.

## Value

list

- `"mu"`:

  Maximum concentration vector

- `"axes"`:

  Axes of Bingham ellipse

- `"k"`:

  two-column vector containing the estimates for the maximum
  (\\k\_\text{min}\\) and minimum concentration (\\k\_\text{max}\\).

- `"a95"`:

  two-column vector containing the estimates for the maximum and minimum
  95% confidence cone. In radians when `x` is an object of class
  `"Vec3"`,otherwise in degrees.

- `"beta"`:

  The shape factor of the distribution given by the ratio
  \\\frac{k\_\text{min}}{k\_\text{max}}\\.

## See also

[`inertia_tensor.spherical()`](https://tobiste.github.io/structr/reference/inertia.md)

## Examples

``` r
set.seed(1234)
x <- rfb(100, mu = Line(120, 50), k = 15, A = diag(c(-5, 0, 5)))

stereoplot()
stereo_point(x)


bingham_statistics(x)
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 118.97488  48.74017 
#> 
#> $axes
#> Line object (n = 2):
#>        azimuth   plunge
#> [1,]  13.12232 13.47624
#> [2,] 272.29509 38.09204
#> 
#> $k
#> [1] 3.487357 2.466239
#> 
#> $a95
#> [1] 7.496868 8.914776
#> 
#> $beta
#> [1] 0.7071944
#> 
```
