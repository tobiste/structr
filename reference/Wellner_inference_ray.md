# Two-sample test, based on permutations and Wellner's Rayleigh-style T-statistic (Wellner, 1979, Example 1a).

Assumes large sample sizes, specifically that choose(nx + ny, nx) \>\>
numPerms. For small sample sizes, see rayWellnerExactInference.

## Usage

``` r
Wellner_inference_ray(xs, ys, numPerms)
```

## Arguments

- xs:

  A list of rays.

- ys:

  A list of rays.

- numPerms:

  A real number (positive integer). The number of permutations, say
  1,000 or 10,000.

## Value

A real number, between 0 and 1 inclusive. The fraction of tests in which
T exceeds the original T for the data. You can interpret this as a
p-value for the null hypothesis that the two populations are identical
(not just that their means are identical). In other words, small values
of p indicate that the distinction between the two populations is
meaningful.
