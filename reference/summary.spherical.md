# Summary statistics

Calculates the arithmetic mean, variance, 68% cone, and the confidence
cone around the mean.

## Usage

``` r
# S3 method for class 'spherical'
summary(object, ...)
```

## Arguments

- object:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- ...:

  parameters passed to
  [`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md),
  [`sph_var()`](https://tobiste.github.io/structr/reference/stats.md),
  [`delta()`](https://tobiste.github.io/structr/reference/stats.md), and
  [`sph_confidence_angle()`](https://tobiste.github.io/structr/reference/stats.md)

## Value

named vector

## Examples

``` r
set.seed(20250411)
summary(rvmf(100, mu = Line(90, 20), k = 20))
#>         azimuth          plunge        variance        68% cone confidence cone 
#>     90.34420929     21.23577351      0.04613671     17.47209466      3.07781518 
```
