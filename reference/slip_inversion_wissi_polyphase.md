# Polyphase stress inversion via spectral clustering on \\S^5\\ (Stage 5).

Identifies `k` stress phases automatically using the *eigengap
heuristic*, then runs
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md)
on each phase subset.

## Usage

``` r
slip_inversion_wissi_polyphase(
  x,
  weights = NULL,
  k_max = 4L,
  sigma_K_deg = 30,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates. Object must be complete, i.e. no `NA` values. For
  Michael's, Angelier's, and Yamaji-Sato's methods, at least 4 rows of
  fault measurements are required, while Hansen's method requires at
  least 7.

- weights:

  Optional length-`N` numeric vector of non-negative data quality
  weights. `NA` values are replaced with the observed mean before
  scaling (neutral imputation). The vector is normalised internally so
  that `mean(weights) = 1`, meaning only the ratios between weights
  matter. Use
  [`scale_weights()`](https://tobiste.github.io/structr/reference/scale_weights.md)
  to construct weights from field quality ranks, measurement errors, and
  prior RUP values. Default: uniform weights.

- k_max:

  Maximum number of phases to consider. Default `4`.

- sigma_K_deg:

  Affinity bandwidth in degrees of angular stress distance. Faults
  within this distance are considered similar. Default `30.` Increase to
  merge nearby phases, decrease to split them.

- seed:

  Optional RNG seed for k-means reproducibility.

- ...:

  Additional arguments passed to
  [`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md)
  for each phase.

## Value

A named list with:

- `assignment`:

  integer vector of length N (phase label 1..k per fault)

- `k_opt`:

  number of phases identified

- `gaps`:

  Laplacian eigenvalue gaps (*eigengap* criterion)

- `phase_results`:

  list of `k`
  [`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md)
  results, one per phase

- `D_mat`:

  N x N pairwise ASD matrix between fault poles (degrees)

## See also

Other wissi:
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md),
[`slip_inversion_wissi_boot()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi_boot.md)

## Examples

``` r
res <- slip_inversion_wissi_polyphase(angelier1990$KAM, sigma_K_deg = 15)

# check amount of clusters detected
res$k_opt
#> [1] 3


# plot the results in lambert projection:
cols <- assign_col_d(seq_len(res$k_opt))#' 
par(mfrow = c(1, 2))

stereoplot()
angelier(angelier1990$KAM, col = assign_col_d(res$assignment))
legend('topleft', title = 'Phases', fill = cols, legend = 1:3)


stereoplot()
angelier(angelier1990$KAM, col = 'grey')
for(k in seq_len(res$k_opt)){
   res_k <- res$phase_results[[k]]
   points(res_k$principal_axes, col = cols[k], pch = 16:18, cex = 2) 
 }
legend('topright', title = "Principal axes", pch = 16:18, col = 'black', 
   legend = c('S1', 'S2', 'S3'))

 
 dev.off()
#> null device 
#>           1 
```
