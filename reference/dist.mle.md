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
#> [1,] 123.2236  56.17187
#> [2,] 195.9225 -11.27133
#> [3,] 278.9303  31.41707
#> 
#> $param
#>      kappa       beta        psi 
#>  4.8028979 -1.0240553  0.1395559 
#> 
#> $logcon
#> [1] 5.121457
#> 
#> $loglik
#> [1] -128.4912
#> 
```
