# One-sample inference about the mean of the Watson distribution.

Assumes large concentration — either kappa \>\> 0 or kappa \<\< 0. From
Mardia and Jupp (2000, Section 10.7.3).

## Usage

``` r
watson_inference(x, alpha, shape)

# S3 method for class 'Vec3'
watson_inference(x, alpha = 0.05, shape = NULL)

# S3 method for class 'Line'
watson_inference(x, alpha = 0.05, shape = NULL)

# S3 method for class 'Plane'
watson_inference(x, alpha = 0.05, shape = NULL)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, or `"Plane"`, where the rows are
  the observations and the columns are the coordinates.

- alpha:

  A real number, between 0 and 1. The significance level for the
  confidence region.

- shape:

  `NULL` or character, either 'bipolar' or 'girdle'. If `NULL`, then
  this function chooses automatically.

## Value

A list with members `$shape`, `$tBar`, `$rhs`, `$pvalue`.

- `shape`:

  is either 'bipolar' or 'girdle'. If 'bipolar', then the confidence
  region consists of all lines u such that `u^T %*% $tBar %*% u > $rhs`.
  If 'girdle', then the confidence region consists of all lines u such
  that `u^T %*% $tBar %*% u < $rhs`.

- `tBar`:

  orientation tesor

- `rhs`:

- `pvalue`:

  is an R function that takes as input a line u0 and produces as output
  a real number in `[0, 1]` — the p-value for the null hypothesis that
  the Watson mean is u0.

## See also

[`rwatson()`](https://tobiste.github.io/structr/reference/rwatson.md)
for simulating a Watson distribution, and
[`watson_MLE()`](https://tobiste.github.io/structr/reference/watson-mle.md)
to estimate distribution parameters.

Other distribution-inference:
[`bingham-inference`](https://tobiste.github.io/structr/reference/bingham-inference.md),
[`fisher-inference`](https://tobiste.github.io/structr/reference/fisher-inference.md)

## Examples

``` r
r <- watson_inference(example_lines)
print(r)
#> $shape
#> [1] "bipolar"
#> 
#> $tBar
#> Orientation tensor
#>            [,1]      [,2]       [,3]
#> [1,] 0.16204564 0.2552464 0.08096582
#> [2,] 0.25524643 0.7382838 0.19306540
#> [3,] 0.08096582 0.1930654 0.09967054
#> 
#> $rhs
#> [1] 0.8864225
#> 
#> $pvalue
#> function (u0) 
#> {
#>     u0 <- as.vector(Vec3(u0))
#>     f <- as.numeric((t1 - u0 %*% tBar %*% u0) * (n - 1)/(1 - 
#>         t1))
#>     1 - stats::pf(f, 2, 2 * n - 2)
#> }
#> <bytecode: 0x55ac3a8eabd8>
#> <environment: 0x55ac3a8e9818>
#> 

r$pvalue(Line(60, 10))
#> [1] 2.499618e-08
```
