# Tests of rotational symmetry of a Bingham distribution

Tests of rotational symmetry of a Bingham distribution

## Usage

``` r
watson_symmetry_test(lambda, k, n, level = 0.95)
```

## Arguments

- lambda:

  numeric.

- k:

  Bingham concentrations \\(0, \kappa\_\text{major},
  \kappa\_\text{minor})\\

- n:

  integer. Number of vectors

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

## Value

object of class `"Test"`

## Details

Whether a Bingham distribution is a Watson distribution: \$\$w_b = n
\sum((k_i - \text{mean}(k)) (\lambda_i - \text{mean}(\lambda)))\$\$ over
the two smallest eigenvalues tests the bipolar form, \\w_g\\ over the
two largest the girdle form, each asymptotically chi-square with 2
degrees of freedom (MJ eq. 10.7.8-10.7.18). \\k_i\\ are the
maximum-likelihood Bingham concentrations \\(0, \kappa\_\text{major},
\kappa\_\text{minor})\\ and \\\lambda_i\\ the eigenvalues of the
orientation tensor, both in descending order. Mardia & Jupp's calcite
example gives \\w_g = 1.33\\ (girdle accepted) and \\w_b = 9.97\\.

## Examples

``` r
r <- rbingham_eig(100, c(100, 1, 0), "Line")
e <- ot_eigen(r)$values
k <- c(0, bingham_statistics(r)$k)
watson_symmetry_test(e, k, n = nrow(r))
#> Error in watson_symmetry_test(e, k, n = nrow(r)): could not find function "watson_symmetry_test"
```
