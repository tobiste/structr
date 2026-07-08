# Maximum likelihood estimation of Spherical Rotational Symmetric Distributions

MLE parameters of a von Mises-Fisher or Kent distribution.

## Usage

``` r
kent_MLE(x)

vmf_MLE(x)
```

## Source

Adapted from
[`Directional::kent.mle()`](https://rdrr.io/pkg/Directional/man/kent.mle.html)
and `Directional::vmf.mle()`

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

## See also

[`fisher_inference()`](https://tobiste.github.io/structr/reference/fisher-inference.md)
for confidence regions, and
[`rvmf()`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)
to simulate a distribution.
[`fisher_MLE()`](https://tobiste.github.io/structr/reference/fisher-mle.md)
is an alternative MLE function for the Fisher distribution.

Other distribution-MLE:
[`bingham-mle`](https://tobiste.github.io/structr/reference/bingham-mle.md),
[`fisher-mle`](https://tobiste.github.io/structr/reference/fisher-mle.md),
[`watson-mle`](https://tobiste.github.io/structr/reference/watson-mle.md)

## Examples

``` r
set.seed(20250411)
x <- rvmf(100, mu = Ray(120, 50), k = 5)
vmf_MLE(x)
#> $loglik
#> [1] -122.7913
#> 
#> $mu
#> Ray object (n = 1):
#>   azimuth    plunge 
#> 118.32579  50.88713 
#> 
#> $kappa
#> [1] 5.000134
#> 
 
x2 <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
kent_MLE(x2)
#> $G
#> Ray object (n = 3):
#>       azimuth    plunge
#> [1,] 122.5000  54.81324
#> [2,] 185.9610 -17.48617
#> [3,] 265.7139  29.45251
#> 
#> $param
#>       kappa        beta         psi 
#>  4.63302477 -0.71444711 -0.06849438 
#> 
#> $logcon
#> [1] 4.964508
#> 
#> $loglik
#> [1] -131.3965
#> 
```
