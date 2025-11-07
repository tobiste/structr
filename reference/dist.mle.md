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
#> mean  122.7104 56.01447
#> major 326.6274 31.64367
#> minor 229.6679 11.12364
#> 
#> $param
#>      kappa       beta        psi 
#>  4.6647142  0.2695082 -0.7290731 
#> 
#> $logcon
#> [1] 4.969102
#> 
#> $loglik
#> [1] -130.2265
#> 
vmf_mle(x)
#> $loglik
#> [1] -130.2242
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 122.71042  56.01447 
#> 
#> $kappa
#> [1] 4.639821
#> 
```
