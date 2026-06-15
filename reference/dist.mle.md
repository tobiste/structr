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
#> mean  111.0004 57.79568
#> major 344.7370 20.43284
#> minor 245.3041 23.74557
#> 
#> $param
#>      kappa       beta        psi 
#>  3.7967715  1.1077790 -0.4002009 
#> 
#> $logcon
#> [1] 4.376522
#> 
#> $loglik
#> [1] -150.7201
#> 
vmf_mle(x)
#> $loglik
#> [1] -157.9062
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 111.00045  57.79568 
#> 
#> $kappa
#> [1] 3.495424
#> 
```
