# Permutation Tests of Regressions

`perm_rsq` retrieves the \\R^2\\ values from successful `n_perm`
permutations of a regression function. `perm_rsq_pvalue` post-processes
this vector to count how many of its entries exceed the original \\R^2\\
value. That fraction is a p-value for the null hypothesis that the
regression result is meaningless.

## Usage

``` r
perm_rsq(n_perm, FUN, x, ...)

perm_rsq_pvalue(r_perm, r)
```

## Arguments

- n_perm:

  A real number (positive integer). The number of permutations to try.

- FUN:

  An R function. Must return a result list with fields `convergence`,
  `min_eigenvalue`, `r_squared` Typically
  [`regression_greatcircle()`](https://tobiste.github.io/structr/reference/best_fit.md)
  or a similar regression function.

- x:

  A real vector. The values of the independent variable. Assumed to be
  the first argument passed to `FUN`.

- ...:

  Other arguments passed to `FUN`, after `x`.

- r_perm:

  vector of permutated \\R^2\\ values

- r:

  original \\R^2\\ value

## Value

`perm_rsq` returns a a real vector. The maximum length is `n_perms`.
Often the length is less than `n_perms`, because the regression failed
(as signaled by error != 0 or `min_eigenvalue` \<= 0).

`perm_rsq_pvalue` returns the p-value.

## See also

Other geodesic-regression:
[`best_fit`](https://tobiste.github.io/structr/reference/best_fit.md)

## Examples

``` r
set.seed(20250411)
data("gray_example")

# original regression
bestgc_clea <- regression_greatcircle(gray_example[1:8, ])

# permutation test
pr <- perm_rsq(100, FUN = regression_greatcircle, x = gray_example[1:8, ])

# p-value
perm_rsq_pvalue(pr, bestgc_clea)
#> [1] 0.03125
```
