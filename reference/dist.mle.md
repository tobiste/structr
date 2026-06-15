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
#> mean  109.6736 57.83627
#> major 336.6182 23.23475
#> minor 237.1613 20.94249
#> 
#> $param
#>      kappa       beta        psi 
#>  4.0830047  0.9072054 -0.5407599 
#> 
#> $logcon
#> [1] 4.561777
#> 
#> $loglik
#> [1] -143.7888
#> 
vmf_mle(x)
#> $loglik
#> [1] -147.9456
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 109.67356  57.83627 
#> 
#> $kappa
#> [1] 3.875419
#> 
```
