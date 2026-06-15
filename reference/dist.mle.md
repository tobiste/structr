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
#> mean  121.8427 56.85973
#> major 213.6296  1.16627
#> minor 304.3904 33.11434
#> 
#> $param
#>     kappa      beta       psi 
#> 4.3553407 0.8952561 0.5167863 
#> 
#> $logcon
#> [1] 4.765223
#> 
#> $loglik
#> [1] -137.6312
#> 
vmf_mle(x)
#> $loglik
#> [1] -141.3212
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 121.84275  56.85973 
#> 
#> $kappa
#> [1] 4.146815
#> 
```
