# Bootstrap uncertainty for a wissi() result.

Yields B stress tensors from resampled datasets. The dispersion
Theta-bar on S^5 approximates the data noise level (Eq. 37: Theta ~
d-bar).

## Usage

``` r
wissi_bootstrap(normals, slips, weights = NULL, B = 500L, seed = NULL, ...)
```

## Arguments

- normals, slips, weights:

  As for wissi().

- B:

  Number of bootstrap replicates. Default 500.

- seed:

  Optional RNG seed.

- ...:

  Additional arguments passed to wissi().

## Value

A named list with: optimal : wissi() result for the full dataset thetas
: length-B vector of angular stress distances from optimal dispersion :
mean Theta (approximates noise level p of data) sd : standard deviation
of Theta values D_bar : mean Orife-Lisle distance from optimal DM_bar :
mean Michael distance from optimal
