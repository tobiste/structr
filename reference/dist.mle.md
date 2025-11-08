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
#> mean  120.0617 53.16399
#> major 352.3547 24.61499
#> minor 249.7065 25.54456
#> 
#> $param
#>      kappa       beta        psi 
#>  4.6533355  0.6255517 -0.3342611 
#> 
#> $logcon
#> [1] 4.974807
#> 
#> $loglik
#> [1] -130.8535
#> 
vmf_mle(x)
#> $loglik
#> [1] -132.3082
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 120.06169  53.16399 
#> 
#> $kappa
#> [1] 4.543299
#> 
```
