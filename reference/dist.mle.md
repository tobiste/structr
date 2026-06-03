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
#> mean  125.9532 57.20234
#> major 353.8098 23.38288
#> minor 253.9362 21.63231
#> 
#> $param
#>      kappa       beta        psi 
#>  4.0158159  0.6416725 -0.2747487 
#> 
#> $logcon
#> [1] 4.489065
#> 
#> $loglik
#> [1] -145.0527
#> 
vmf_mle(x)
#> $loglik
#> [1] -146.973
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 125.95319  57.20234 
#> 
#> $kappa
#> [1] 3.914275
#> 
```
