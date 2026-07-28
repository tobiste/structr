# Bootstrap uncertainty for a WISSI result.

Yields `n_iter` stress tensors from resampled datasets. The dispersion
Theta-bar on S^5 approximates the data noise level (Eq. 37: Theta ~
d-bar).

## Usage

``` r
slip_inversion_wissi_boot(x, n_iter = 500L, seed = NULL, ...)
```

## Arguments

- x:

  object of class `"Pair"` or `"Fault"` with at least 4 rows.

- n_iter:

  Number of bootstrap replicates. Default `500`.

- seed:

  Optional RNG seed.

- ...:

  Additional arguments passed to
  [`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md).

## Value

A named list with:

- `optimal`:

  slip_inversion_wissi() result for the full dataset

- `thetas`:

  length-B vector of angular stress distances from optimal

- `dispersion`:

  mean Theta (approximates noise level p of data)

- `sd`:

  standard deviation of Theta values

- `D_bar`:

  mean Orife-Lisle distance from optimal

- `DM_bar`:

  mean Michael distance from optimal

## See also

Other wissi:
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md),
[`slip_inversion_wissi_polyphase()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi_polyphase.md)

## Examples

``` r
res <- slip_inversion_wissi_boot(angelier1990$AVB, n_iter = 4)
#> Error in slip_inversion_wissi_boot(angelier1990$AVB, n_iter = 4): could not find function "slip_inversion_wissi_boot"

stereoplot()
angelier(angelier1990$KAM, col = 'grey')

lines(res$optimal$principal_axes, res$sd, col = 2:4)
#> Error: object 'res' not found
points(res$optimal$principal_axes, pch = 16:18, cex = 2, col= 2:4)
#> Error: object 'res' not found
text(res$optimal$principal_axes, 
 label = rownames(res$optimal$principal_axes), col= 2:4, adj = -.25)
#> Error: object 'res' not found
```
