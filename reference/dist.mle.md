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
#> mean  122.4002 49.42138
#> major 347.0914 31.33563
#> minor 242.1126 23.00121
#> 
#> $param
#>      kappa       beta        psi 
#>  4.5105035 -0.5049675 -0.4772698 
#> 
#> $logcon
#> [1] 4.857341
#> 
#> $loglik
#> [1] -133.7262
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.6002
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 122.40019  49.42138 
#> 
#> $kappa
#> [1] 4.439309
#> 
```
