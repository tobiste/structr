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
#> mean  124.1607 54.98262
#> major 333.0631 31.52448
#> minor 234.4801 13.67472
#> 
#> $param
#>      kappa       beta        psi 
#>  4.7701882  0.3511715 -0.6389088 
#> 
#> $logcon
#> [1] 5.054357
#> 
#> $loglik
#> [1] -128.0912
#> 
vmf_mle(x)
#> $loglik
#> [1] -128.3083
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 124.16069  54.98262 
#> 
#> $kappa
#> [1] 4.730259
#> 
```
