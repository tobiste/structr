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
#>          azimuth   plunge
#> mean  124.208057 54.13234
#> major   1.359726 21.41403
#> minor 259.694309 27.27438
#> 
#> $param
#>      kappa       beta        psi 
#>  4.4788690  0.6431361 -0.1692155 
#> 
#> $logcon
#> [1] 4.840511
#> 
#> $loglik
#> [1] -134.5738
#> 
vmf_mle(x)
#> $loglik
#> [1] -136.2243
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 124.20806  54.13234 
#> 
#> $kappa
#> [1] 4.366956
#> 
```
