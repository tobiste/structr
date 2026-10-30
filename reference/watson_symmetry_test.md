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

## See also

Other tests:
[`bingham_test()`](https://tobiste.github.io/structr/reference/bingham_test.md),
[`fisher_ftest()`](https://tobiste.github.io/structr/reference/fisher_ftest.md),
[`fisher_goodness_of_fit()`](https://tobiste.github.io/structr/reference/fisher_goodness_of_fit.md),
[`fisher_kent_test()`](https://tobiste.github.io/structr/reference/fisher_kent_test.md),
[`rayleigh_test3()`](https://tobiste.github.io/structr/reference/rayleigh_test3.md),
[`wellner()`](https://tobiste.github.io/structr/reference/wellner.md)

## Examples

``` r
r <- rbingham_eig(100, c(100, 1, 0), "Line")
e <- ot_eigen(r)$values
k <- c(0, bingham_statistics(r)$k)
watson_symmetry_test(e, k, n = nrow(r))
#> $bipolar
#> 
#>  Watson bipolar symmetry
#> 
#> statistic = 3290.3, df = 2, p-value < 2.2e-16
#> critical value (95% confidence level) = 5.9915
#> null hypothesis rejected at the 95% confidence level: yes
#> source: Mardia & Jupp (2000), eq. 10.7.8-10.7.18 
#> 
#> 
#> $girdle
#> 
#>  Watson girdle symmetry
#> 
#> statistic = -1046.9, df = 2, p-value = 1
#> critical value (95% confidence level) = 5.9915
#> null hypothesis rejected at the 95% confidence level: no
#> source: Mardia & Jupp (2000), eq. 10.7.8-10.7.18 
#> 
#> 
```
