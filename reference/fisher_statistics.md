# Fisher's Statistics

Estimates concentration parameter, angular standard deviation, and
confidence limit.

## Usage

``` r
fisher_statistics(x, w = NULL, conf.level = 0.95, na.rm = TRUE)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- w:

  numeric. Optional weights for each observation.

- conf.level:

  numeric. Level of confidence.

- na.rm:

  logical. Whether `NA` values should be removed before the computation
  proceeds.

## Value

list, with

- `"mu"`:

  Mean vector

- `"k"`:

  estimated concentration parameter \\\kappa\\ for the von Mises-Fisher
  distribution

- `"csd"`:

  estimated angular standard deviation enclosing 63% of the orientation
  data. Angle is in degrees if `x` is a spherical object, and radians if
  otherwise.

- `"alpha"`:

  Confidence limit for given `conf.level`. Angle is in degrees if `x` is
  a spherical object, and radians if otherwise.

## Examples

``` r
set.seed(20250411)
x <- rvmf(100, mu = Line(120, 50), k = 5)
fisher_statistics(x)
#> $mu
#> Ray object (n = 1):
#>   azimuth    plunge 
#> 118.32579  50.88713 
#> 
#> $k
#> [1] 4.952381
#> 
#> $csd
#> [1] 36.39804
#> 
#> $csd_2s
#> [1] 62.91019
#> 
#> $alpha
#> [1] 7.103673
#> 
```
