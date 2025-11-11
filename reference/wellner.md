# Wellner's Two-Sample Test

Wellner's (1979, Example 1a) Rayleigh-style T-statistic, which
quantifies the dissimilarity between two sets of vectors. The statistic
increases with the degree of difference between the datasets.

## Usage

``` r
wellner(x, y)

# S3 method for class 'Ray'
wellner(x, y)

# S3 method for class 'Line'
wellner(x, y)

# S3 method for class 'Vec3'
wellner(x, y)

# S3 method for class 'Plane'
wellner(x, y)

wellner_inference(x, y, n_perm)

# S3 method for class 'Vec3'
wellner_inference(x, y, n_perm = 1000)

# S3 method for class 'Line'
wellner_inference(x, y, n_perm = 1000)

# S3 method for class 'Ray'
wellner_inference(x, y, n_perm = 1000)

# S3 method for class 'Plane'
wellner_inference(x, y, n_perm = 1000)
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

## Value

`wellner()` computes Wellner's T-statistic, a non-negative measure of
dissimilarity between two datasets. The value is zero when the datasets
are identical.

`wellner_inference()` estimates the fraction of permutation tests in
which the computed T-statistic exceeds the observed T for the original
data (a number between 0 and 1). This value can be interpreted as a
p-value for the null hypothesis that the two populations are identical
(not merely that their means coincide). Thus, smaller p-values indicate
stronger evidence that the two populations differ in a meaningful way.

## Details

`wellner_inference()` performs a permutation-based inference using
Wellner's T-statistic to assess whether the two datasets are drawn from
the same population.

## References

Jon A. Wellner. "Permutation Tests for Directional Data." Ann. Statist.
7(5) 929-943, September, 1979.
[doi:10.1214/aos/1176344779](https://doi.org/10.1214/aos/1176344779)

## Examples

``` r
test <- rvmf(100)
wellner(test, Line(120, 50))
#> [1] 7.620339
wellner_inference(test, Line(120, 50))
#> [1] 0.022
```
