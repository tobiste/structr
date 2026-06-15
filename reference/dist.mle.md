# MLE of spherical rotational symmetric distributions

Estimates the parameters of a von Mises-Fisher or Kent distribution.

## Usage

``` r
kent_mle(x)

vmf_mle(x)
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

## Examples

``` r
x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
kent_mle(x)
#> $G
#> Line object (n = 3):
#>        azimuth   plunge
#> mean  113.2616 55.11472
#> major 346.9808 22.42006
#> minor 245.7722 25.22680
#> 
#> $param
#>      kappa       beta        psi 
#>  3.6438422  1.0098947 -0.3909601 
#> 
#> $logcon
#> [1] 4.254582
#> 
#> $loglik
#> [1] -154.3166
#> 
vmf_mle(x)
#> $loglik
#> [1] -160.4942
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 113.26159  55.11472 
#> 
#> $kappa
#> [1] 3.401661
#> 
```
