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
#> mean  122.8168 54.71873
#> major  16.5715 11.19626
#> minor 279.1990 32.95466
#> 
#> $param
#>     kappa      beta       psi 
#> 4.7782467 0.5033282 0.1417142 
#> 
#> $logcon
#> [1] 5.066478
#> 
#> $loglik
#> [1] -128.0865
#> 
vmf_mle(x)
#> $loglik
#> [1] -128.8756
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 122.81684  54.71873 
#> 
#> $kappa
#> [1] 4.703308
#> 
```
