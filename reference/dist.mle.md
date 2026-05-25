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
#> mean  127.0481 55.22457
#> major 225.5824  5.88340
#> minor 319.5876 34.13016
#> 
#> $param
#>     kappa      beta       psi 
#> 4.9185437 1.0855948 0.7357289 
#> 
#> $logcon
#> [1] 5.21757
#> 
#> $loglik
#> [1] -126.3626
#> 
vmf_mle(x)
#> $loglik
#> [1] -131.1206
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 127.04807  55.22457 
#> 
#> $kappa
#> [1] 4.59807
#> 
```
