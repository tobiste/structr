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
#>         azimuth    plunge
#> mean  126.80323 57.156559
#> major  26.05902  6.862124
#> minor 291.75658 31.939310
#> 
#> $param
#>     kappa      beta       psi 
#> 4.2283761 0.8677502 0.3390515 
#> 
#> $logcon
#> [1] 4.666714
#> 
#> $loglik
#> [1] -140.4084
#> 
vmf_mle(x)
#> $loglik
#> [1] -143.9973
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 126.80323  57.15656 
#> 
#> $kappa
#> [1] 4.035245
#> 
```
