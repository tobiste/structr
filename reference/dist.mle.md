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
#> mean  113.7037 52.77888
#> major 355.6945 19.63281
#> minor 253.7122 30.19843
#> 
#> $param
#>      kappa       beta        psi 
#>  4.6018292  0.9363895 -0.2525779 
#> 
#> $logcon
#> [1] 4.957601
#> 
#> $loglik
#> [1] -132.4459
#> 
vmf_mle(x)
#> $loglik
#> [1] -136.2161
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 113.70372  52.77888 
#> 
#> $kappa
#> [1] 4.367318
#> 
```
