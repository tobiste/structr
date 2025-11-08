# Two-sample test, based on permutations and Wellner's Rayleigh-style T-statistic (Wellner, 1979, Example 1a).

Assumes small sample sizes. Deterministically generates all choose(nx +
ny, nx) reassignments of the data to the two groups. For large sample
sizes, see rayWellnerInference.

## Usage

``` r
wellner_exact_inference_ray(xs, ys)
```

## Arguments

- xs:

  A list of rays.

- ys:

  A list of rays.

## Value

A real number, between 0 and 1 inclusive. The fraction of tests in which
T exceeds the original T for the data. You can interpret this as a
p-value for the null hypothesis that the two populations are identical
(not just that their means are identical). In other words, small values
of p indicate that the distinction between the two populations is
meaningful.
