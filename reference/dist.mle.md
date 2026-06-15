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
#>         azimuth    plunge
#> mean  122.14196 53.916880
#> major  24.74456  5.360276
#> minor 290.89906 35.556277
#> 
#> $param
#>     kappa      beta       psi 
#> 4.6060597 0.7474375 0.3105677 
#> 
#> $logcon
#> [1] 4.94582
#> 
#> $loglik
#> [1] -132.0177
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.2794
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 122.14196  53.91688 
#> 
#> $kappa
#> [1] 4.453727
#> 
```
