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
#> mean  130.85627 52.952878
#> major  32.49197  6.266239
#> minor 297.85960 36.334762
#> 
#> $param
#>     kappa      beta       psi 
#> 4.4346232 0.9681485 0.4220158 
#> 
#> $logcon
#> [1] 4.832381
#> 
#> $loglik
#> [1] -136.0532
#> 
vmf_mle(x)
#> $loglik
#> [1] -140.3348
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 130.85627  52.95288 
#> 
#> $kappa
#> [1] 4.188621
#> 
```
