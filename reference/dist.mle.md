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
#> mean  115.1643 56.87300
#> major 341.3469 24.31396
#> minor 241.3277 21.06023
#> 
#> $param
#>      kappa       beta        psi 
#>  3.9814468  0.7455562 -0.4783929 
#> 
#> $logcon
#> [1] 4.471553
#> 
#> $loglik
#> [1] -145.9664
#> 
vmf_mle(x)
#> $loglik
#> [1] -148.7276
#> 
#> $mu
#> Line object (n = 1):
#>  azimuth   plunge 
#> 115.1643  56.8730 
#> 
#> $kappa
#> [1] 3.844413
#> 
```
