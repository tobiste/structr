# Reduced stress tensor

Calculates the reduced stress tensor (without optimization).

## Usage

``` r
reduced_stress(fault, method = c("michael", "angelier"))
```

## Arguments

- fault:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- method:

  character. The inversion algorithm. One of `"michael"` (the default)
  for a bootstrapped linear inversion, or `"angelier"` for an iterative
  direct inversion.

## Value

`"ellipsoid"` object. 3x3 matrix (i.e. 2nd rank tensor).

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)

Other stress-tensor:
[`fault_instability_criterion()`](https://tobiste.github.io/structr/reference/fault_instability_criterion.md),
[`stress_shape()`](https://tobiste.github.io/structr/reference/stress_shape.md),
[`tau-comp`](https://tobiste.github.io/structr/reference/tau-comp.md),
[`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md),
[`tau2stress()`](https://tobiste.github.io/structr/reference/tau2stress.md)

## Examples

``` r
f <- angelier1990$TYM
reduced_stress(f)
#> Ellipsoid tensor
#>            [,1]       [,2]       [,3]
#> [1,] -0.8168187  0.1113011 -0.1570288
#> [2,]  0.1113011 -0.6797825 -0.1514454
#> [3,] -0.1570288 -0.1514454  1.4966011
```
