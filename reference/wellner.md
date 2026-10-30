# Wellner's Two-Sample Test

Wellner's (1979, Example 1a) Rayleigh-style T-statistic, which
quantifies the dissimilarity between two sets of vectors. The statistic
increases with the degree of difference between the datasets.

## Usage

``` r
wellner(x, y, n_perm, level)

# S3 method for class 'Ray'
wellner(x, y, n_perm = 1000L, level = 0.95)

# S3 method for class 'Line'
wellner(x, y, n_perm = 1000L, level = 0.95)

# S3 method for class 'Vec3'
wellner(x, y, n_perm = 1000L, level = 0.95)

# S3 method for class 'Plane'
wellner(x, y, n_perm = 1000L, level = 0.95)
```

## Source

modified after `geologyGeometry` (J.R. Davis):
http://www.joshuadavis.us/software/

## Arguments

- x, y:

  objects of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- n_perm:

  integer. Number of permutations.

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

## Value

`wellner()` computes Wellner's T-statistic, a non-negative measure of
dissimilarity between two datasets. The value is zero when the datasets
are identical.

Inference estimates the fraction of permutation tests in which the
computed T-statistic exceeds the observed T for the original data (a
number between 0 and 1). This value can be interpreted as a p-value for
the null hypothesis that the two populations are identical (not merely
that their means coincide). Thus, smaller p-values indicate stronger
evidence that the two populations differ in a meaningful way.

## Details

P-values are estimated using permutation-based inference to assess
whether the two datasets are drawn from the same population.

## References

Jon A. Wellner. "Permutation Tests for Directional Data." Ann. Statist.
7(5) 929-943, September, 1979.
[doi:10.1214/aos/1176344779](https://doi.org/10.1214/aos/1176344779)

## See also

Other tests:
[`bingham_test()`](https://tobiste.github.io/structr/reference/bingham_test.md),
[`fisher_ftest()`](https://tobiste.github.io/structr/reference/fisher_ftest.md),
[`rayleigh_test3()`](https://tobiste.github.io/structr/reference/rayleigh_test3.md)

## Examples

``` r
set.seed(20250411)
test <- rvmf(100)
wellner(test, Line(120, 50))
#> 
#>  Wellner permutation test
#> 
#> statistic = 8.2327, p-value = 0.03896
#> critical value (95% confidence level) = 8.2246
#> null hypothesis rejected at the 95% confidence level: yes
#> further statistics:
#>       n_perm   df_uniform    p_uniform 
#> 1000.0000000    5.0000000    0.1438671 
#> source: Wellner (1979), Ann. Statist. 7, 929-943 
#> 
```
