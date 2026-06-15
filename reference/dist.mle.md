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
#> mean  112.0672 54.57495
#> major 223.1084 14.32664
#> minor 322.1524 31.61201
#> 
#> $param
#>      kappa       beta        psi 
#>  4.6654842 -0.9296244  0.7601233 
#> 
#> $logcon
#> [1] 5.006212
#> 
#> $loglik
#> [1] -131.1158
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.753
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 112.06722  54.57495 
#> 
#> $kappa
#> [1] 4.432456
#> 
```
