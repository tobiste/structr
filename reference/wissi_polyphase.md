# Polyphase stress inversion via spectral clustering on \\S^5\\ (Stage 5).

Identifies k stress phases automatically using the eigengap heuristic,
then runs wissi() on each phase subset.

## Usage

``` r
wissi_polyphase(
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

  object of class `"Pair"` or `"Fault"` with at least 4 rows.

- weights:

  numeric. Weightings for the faults. Must have the same length as `x`

- k_max:

  Maximum number of phases to consider. Default `4`.

- sigma_K_deg:

  Affinity bandwidth in degrees of angular stress distance. Faults
  within this distance are considered similar. Default `30.` Increase to
  merge nearby phases, decrease to split them.

- seed:

  Optional RNG seed for k-means reproducibility.

- ...:

  Additional arguments passed to wissi() for each phase.

## Value

A named list with: assignment : integer vector of length N (phase label
1..k per fault) k_opt : number of phases identified gaps : Laplacian
eigenvalue gaps (eigengap criterion) phase_results : list of k wissi()
results, one per phase D_mat : N x N pairwise ASD matrix between fault
poles (degrees)

## Examples

``` r
res <- wissi_polyphase(angelier1990$TYM)
res$k_opt
#> [1] 1
```
