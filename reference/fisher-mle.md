# Maximum likelihood estimation of the Fisher parameters.

MLE parameters describing a von Mises-Fisher distribution for isotropic
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

[`fisher_inference()`](https://tobiste.github.io/structr/reference/fisher-inference.md)
for confidence regions, and
[`rvmf()`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)
to simulate a distribution.
[`vmf_MLE()`](https://tobiste.github.io/structr/reference/dist.mle.md)
is an alternative MLE function.

Other distribution-MLE:
[`bingham-mle`](https://tobiste.github.io/structr/reference/bingham-mle.md),
[`dist.mle`](https://tobiste.github.io/structr/reference/dist.mle.md),
[`watson-mle`](https://tobiste.github.io/structr/reference/watson-mle.md)

## Examples

``` r
set.seed(20250411)
x <- rvmf(100, mu = Ray(120, 50), k = 5)

fisher_MLE(x)
#> $muHat
#> Ray object (n = 1):
#>   azimuth    plunge 
#> 118.32579  50.88713 
#> 
#> $rBar
#> [1] 0.8000961
#> 
#> $kappaHat
#> [1] 5.000538
#> 
```
