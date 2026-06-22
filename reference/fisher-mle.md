# Maximum likelihood estimation of the Fisher parameters.

MLE parameters describing a Fisher distribution for isotropic
directional vectors. Based on Mardia and Jupp (2000, p. 198).

## Usage

``` r
fisher_MLE(x)

# S3 method for class 'Vec3'
fisher_MLE(x)

# S3 method for class 'Ray'
fisher_MLE(x)

# S3 method for class 'Plane'
fisher_MLE(x)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Ray"`, or `"Plane"`, where the rows are
  the observations and the columns are the coordinates.

## Value

A list with members

- `muHat`:

  the mean ray, identical to
  [`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md)

- `rBar`:

  non-negative real number. Mean resultant length.

- `kappaHat`:

  a positive real number. Concentration parameter \\\kappa\\ of the
  Fisher distribution

## See also

Other distribution-MLE:
[`bingham-mle`](https://tobiste.github.io/structr/reference/bingham-mle.md),
[`dist.mle`](https://tobiste.github.io/structr/reference/dist.mle.md),
[`watson-mle`](https://tobiste.github.io/structr/reference/watson-mle.md)

## Examples

``` r
r <- fisher_MLE(Ray(example_lines))
print(r)
#> $muHat
#> Ray object (n = 1):
#>  azimuth   plunge 
#> 68.51277 20.49587 
#> 
#> $rBar
#> [1] 0.7831482
#> 
#> $kappaHat
#> [1] 4.609632
#> 

plot(example_lines)
points(r$muHat, col = 'red')
```
