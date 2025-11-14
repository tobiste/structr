# Sample covariance matrix, approximated in the tangent space at a given rotation.

Appropriate only if the sample is tightly concentrated near the center.

## Usage

``` r
rotLeftCovariance(rs, center)
```

## Arguments

- rs:

  A list of rotation matrices.

- center:

  A rotation matrix. Typically the geodesic mean of the `rs`.

## Value

A 3x3 real matrix (symmetric, non-negative eigenvalues).
