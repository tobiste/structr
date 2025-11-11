# Direction of maximum horizontal stress from the stress tensor

Calculates the direction of maximum horizontal stress using only the
directions of the principal stress and \\R = \frac{S1 - S2}{S1 - S3}\\.
This function implements equations 11 and 10 from Lund and Townend
(2007).

## Usage

``` r
SH(S1, S2, S3, R, tol = .Machine$double.eps^0.5, ortho.tol = 0.005)
```

## Arguments

- S1, S2, S3:

  The principal stress orientations. The variables hold the coordinates
  in the North, East and Down geographical coordinate system, e.g.
  `S1 = c(s1N,s1E,s1D)`. Given as object of class `"Vec3"` or `"Line"`

- R:

  numeric. Relative magnitude of `S2` with respect to `S1` and `S3`: \\R
  = \frac{S1 - S2}{S1 - S3}\\. Values ranging from 0 to 1, with 0 being
  `S1==S2` and 1 being `S2==S3`. Equivalent to the stress shape ratio of
  Gephart & Forsyth (1984).

- tol:

  Tolerance of comparison.

- ortho.tol:

  tolerance angle (in degree) for orthogonality check of the three
  principal stress vectors.

## Value

The direction of SH from North as numeric angle in degrees (radians if
all principal stress axes were given as `"Vec3"` objects).

## References

Lund and Townend, (2007). Calculating horizontal stress orientations
with full or partial knowledge of the tectonic stress tensor, Geophys.
J. Int.,
doi:[doi:10.1111/j.1365-246X.2007.03468.x](https://doi.org/10.1111/j.1365-246X.2007.03468.x)
.

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)
for stress inversion of fault slip data.

Other SH-from-tensor:
[`SH_from_tensor()`](https://tobiste.github.io/structr/reference/SH_from_tensor.md)

## Examples

``` r
# first example from https://www.snsn.se/SH/SHcode/benchmark.out
S1 <- Line(250.89, 70.07)
S3 <- Line(103.01, 17.07)
S2 <- crossprod(S3, S1)
SH(S1, S2, S3, R = 1) #  70.89
#> [1] 70.89

R <- seq(0, 1, .05)
cbind(R, SH = SH(S1, S2, S3, R = R))
#>          R       SH
#>  [1,] 0.00 13.01021
#>  [2,] 0.05 13.18337
#>  [3,] 0.10 13.37695
#>  [4,] 0.15 13.59476
#>  [5,] 0.20 13.84162
#>  [6,] 0.25 14.12371
#>  [7,] 0.30 14.44908
#>  [8,] 0.35 14.82843
#>  [9,] 0.40 15.27621
#> [10,] 0.45 15.81249
#> [11,] 0.50 16.46586
#> [12,] 0.55 17.27842
#> [13,] 0.60 18.31445
#> [14,] 0.65 19.67656
#> [15,] 0.70 21.53704
#> [16,] 0.75 24.20154
#> [17,] 0.80 28.23884
#> [18,] 0.85 34.69668
#> [19,] 0.90 45.01043
#> [20,] 0.95 58.66746
#> [21,] 1.00 70.89000
```
