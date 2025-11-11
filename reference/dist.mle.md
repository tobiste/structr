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
#>         azimuth   plunge
#> mean  122.83131 57.51928
#> major  10.31148 13.70244
#> minor 272.60426 28.81272
#> 
#> $param
#>      kappa       beta        psi 
#> 4.44105598 0.24023893 0.04162708 
#> 
#> $logcon
#> [1] 4.794063
#> 
#> $loglik
#> [1] -135.0334
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.9761
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 122.83131  57.51928 
#> 
#> $kappa
#> [1] 4.422466
#> 
```
