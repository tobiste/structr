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
#>         azimuth     plunge
#> mean  110.74320 52.5605852
#> major  20.24138  0.3842039
#> minor 289.94724 37.4367461
#> 
#> $param
#>     kappa      beta       psi 
#> 5.1762594 0.7622751 0.2810781 
#> 
#> $logcon
#> [1] 5.396428
#> 
#> $loglik
#> [1] -120.6208
#> 
vmf_mle(x)
#> $loglik
#> [1] -122.6146
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 110.74320  52.56059 
#> 
#> $kappa
#> [1] 5.00902
#> 
```
