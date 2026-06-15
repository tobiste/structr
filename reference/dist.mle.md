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
#> mean  113.0823 60.63655
#> major 344.2094 19.44834
#> minor 246.3790 21.09856
#> 
#> $param
#>      kappa       beta        psi 
#>  4.0027910  0.7675148 -0.3908010 
#> 
#> $logcon
#> [1] 4.489236
#> 
#> $loglik
#> [1] -145.4621
#> 
vmf_mle(x)
#> $loglik
#> [1] -148.3971
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 113.08231  60.63655 
#> 
#> $kappa
#> [1] 3.857489
#> 
```
