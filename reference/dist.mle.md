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
#> mean  112.2454 55.36948
#> major 344.5699 22.88518
#> minor 243.4891 24.48041
#> 
#> $param
#>      kappa       beta        psi 
#>  3.8309294  0.8876933 -0.4290218 
#> 
#> $logcon
#> [1] 4.374518
#> 
#> $loglik
#> [1] -149.6879
#> 
vmf_mle(x)
#> $loglik
#> [1] -154.017
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 112.24543  55.36948 
#> 
#> $kappa
#> [1] 3.640049
#> 
```
