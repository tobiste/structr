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
#>        azimuth    plunge
#> mean  122.9352 55.622103
#> major  27.0212  4.032176
#> minor 294.2884 34.073311
#> 
#> $param
#>     kappa      beta       psi 
#> 3.9093067 0.8316299 0.3661257 
#> 
#> $logcon
#> [1] 4.426098
#> 
#> $loglik
#> [1] -147.7472
#> 
vmf_mle(x)
#> $loglik
#> [1] -151.3769
#> 
#> $mu
#> Line object (n = 1):
#>  azimuth   plunge 
#> 122.9352  55.6221 
#> 
#> $kappa
#> [1] 3.74091
#> 
```
