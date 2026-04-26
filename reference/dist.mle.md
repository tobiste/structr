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
#> mean  127.00917 54.202676
#> major  33.08843  2.822952
#> minor 301.06155 35.650733
#> 
#> $param
#>     kappa      beta       psi 
#> 4.7439679 0.7291103 0.4644638 
#> 
#> $logcon
#> [1] 5.052034
#> 
#> $loglik
#> [1] -129.1097
#> 
vmf_mle(x)
#> $loglik
#> [1] -131.1539
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 127.00917  54.20268 
#> 
#> $kappa
#> [1] 4.596529
#> 
```
