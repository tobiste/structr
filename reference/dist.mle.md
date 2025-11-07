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
#> mean  113.0885 55.594756
#> major 205.0742  1.359349
#> minor 296.0041 34.370653
#> 
#> $param
#>     kappa      beta       psi 
#> 4.1644817 0.9926607 0.3802058 
#> 
#> $logcon
#> [1] 4.6313
#> 
#> $loglik
#> [1] -142.0666
#> 
vmf_mle(x)
#> $loglik
#> [1] -147.0126
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 113.08852  55.59476 
#> 
#> $kappa
#> [1] 3.912686
#> 
```
