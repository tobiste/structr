# Angelier's Ratio Upsilon (RUP)

The per-fault "ratio upsilon" (RUP) parameter after Angelier (1990) as
an estimator for the quality-of-fit of fault-slip inversions.

## Usage

``` r
tau2rup(tau, fault, lambda = sqrt(3)/2)
```

## Arguments

- tau:

  symmetric 3x3 matrix. The (reduced) stress tensor.

- fault:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- lambda:

  numeric. The maximum shear stress. \\\sqrt(3)/2\\ by default.

## Value

numeric. The per-fault "ratio upsilon" (RUP) parameter in percent.

## Details

The RUP estimator varies from 0%, indicating maximum shear stress
parallel to fault slip and with same sense (i.e. perfect fit), to 200%,
corresponding to maximum shear stress parallel to fault slip but with
opposite sense (i.e. largest misfit). The quality of the fit is good if
RUP \\\le\\ 50%, (potentially) acceptable if 50%\<RUP\\\le\\ 75%, and
poor otherwise.

Angelier (1990) introduced the *ratio upsilon*, which is the ratio of
\\\upsilon_i\\ to the largest shear stress on a fault \\i\\ expressed in
percentage: \$\$\text{RUP} = \upsilon_i / \lambda \times 100\$\$ where
\\\lambda = \frac{\sqrt{3}}{2}\\ is the maximum shear stress.

## References

Angelier, J. (1990). Inversion of field data in fault tectonics to
obtain the regional stress—III. A new rapid direct inversion method by
analytical means. Geophys. J. Int, 103, 363–376.
<https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)

Other stress-tensor:
[`fault_instability_criterion()`](https://tobiste.github.io/structr/reference/fault_instability_criterion.md),
[`reduced_stress()`](https://tobiste.github.io/structr/reference/reduced_stress.md),
[`stress_shape()`](https://tobiste.github.io/structr/reference/stress_shape.md),
[`tau-comp`](https://tobiste.github.io/structr/reference/tau-comp.md),
[`tau2stress()`](https://tobiste.github.io/structr/reference/tau2stress.md)

## Examples

``` r
f <- angelier1990$TYM
tau <- reduced_stress(f)
tau2rup(tau, f)
#>  [1] 39.266859 22.364271 42.794333 14.800803 20.461318 64.718606 38.120249
#>  [8] 18.789513 33.810029 26.242500 23.850645 35.834642 34.723973 22.583185
#> [15] 17.643757 28.502268 21.267580 16.202854 51.018881 66.347714 36.709113
#> [22] 50.411469 68.209016 34.135993 25.437568 50.886127 42.192872 53.871795
#> [29]  9.421844 47.329481 30.113931 13.797249 27.532282 46.592528  9.067210
#> [36] 20.237287 13.242937 19.749635
```
