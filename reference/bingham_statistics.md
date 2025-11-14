# Elliptical concentration and confidence cone estimation

Elliptical concentration and confidence cone estimation

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

- `k`:

  two-column vector containing the estimates for the minimum
  (\\k\_\text{min}\\) and maximum concentration (\\k\_\text{max}\\).

- `a95`:

  two-column vector containing the estimates for the minimum and maximum
  95% confidence cone.

- `beta`:

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
#> $k
#> [1] 2.559359 3.930289
#> 
#> $a95
#> [1] 8.751095 7.061806
#> 
#> $beta
#> [1] 0.6511886
#> 
```
