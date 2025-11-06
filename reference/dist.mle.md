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
#> mean  120.9548 54.96717
#> major 357.5319 21.11430
#> minor 256.4078 26.54808
#> 
#> $param
#>      kappa       beta        psi 
#>  4.6840493  0.8985300 -0.2218555 
#> 
#> $logcon
#> [1] 5.01796
#> 
#> $loglik
#> [1] -130.6707
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.0249
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 120.95481  54.96717 
#> 
#> $kappa
#> [1] 4.465198
#> 
```
