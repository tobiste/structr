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
#> mean  124.32382 54.269728
#> major  28.55937  4.132617
#> minor 295.61435 35.415840
#> 
#> $param
#>     kappa      beta       psi 
#> 4.5803989 0.5055935 0.3823792 
#> 
#> $logcon
#> [1] 4.91166
#> 
#> $loglik
#> [1] -132.2352
#> 
vmf_mle(x)
#> $loglik
#> [1] -133.0894
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 124.32382  54.26973 
#> 
#> $kappa
#> [1] 4.507605
#> 
```
